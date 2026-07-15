"""Model sources — one per git checkout of the QSP model.

The model is *fully described by git-tracked files*: every model-change commit
ships the regenerated ``cpp/qsp/ode/ODE_system.cpp``, ``PDAC_model.sbml`` and
``model_structure.json`` alongside the ``.m`` sources. So pointing the explorer
at a different branch needs **no MATLAB and no codegen** — only a checkout and a
C++ compile (~10-20 s with ccache warm; the SUNDIALS deps are seeded from an
existing checkout so there is no network fetch).

Two kinds of source:

- ``worktree`` — an existing ``git worktree`` (primary + any feature worktrees).
  Reflects the *working tree*, so it shows uncommitted edits. Flagged ``dirty``
  when it has them, because then it is NOT the commit it claims to be.
- ``ref`` — an arbitrary branch/tag/sha the user asked for. Materialised as a
  detached worktree under ``SCRATCH_ROOT``. Detached on purpose: git refuses to
  check out a branch that is already checked out in another worktree, and a
  detached HEAD at the same ref sidesteps that entirely.

Build rules follow ``notes/workflows/iterative_model_debugging.md``:
``cpp/sim/build/`` is **always per-worktree, never symlinked** (pool hashes fold
in binary bytes), and only ``_deps/*-src`` is path-portable when seeding — the
``*-build``/``*-subbuild`` dirs bake the donor's absolute path into their
CMakeCache and cmake bails if you copy them.
"""
from __future__ import annotations

import re
import shutil
import subprocess
import threading
from dataclasses import dataclass, field
from pathlib import Path

from qsp_model_explorer.config import ExplorerConfig

# cmake progress lines look like "[ 42%] Building CXX object ..."
_PCT_RE = re.compile(r"^\[\s*(\d+)%\]")

MAX_LOG = 400


def git(args: list[str], cwd: Path) -> str:
    res = subprocess.run(["git", *args], cwd=str(cwd), capture_output=True,
                         text=True, check=True)
    return res.stdout


def git_ok(args: list[str], cwd: Path) -> bool:
    return subprocess.run(["git", *args], cwd=str(cwd),
                          capture_output=True, text=True).returncode == 0


@dataclass
class ModelSource:
    """A checkout of the model + the loaded state derived from it."""

    id: str
    repo: Path
    kind: str            # "worktree" | "ref"
    branch: str          # branch name, or "detached"
    sha: str             # short sha of HEAD
    subject: str         # HEAD commit subject
    cfg: ExplorerConfig  # the model contract (paths, build recipe, scenarios)
    dirty: bool = False  # working tree has uncommitted changes
    ref: str = ""        # the ref the user asked for (kind == "ref")

    status: str = "cold"   # cold | building | loading | ready | error
    message: str = ""
    progress: int = 0
    log: list[str] = field(default_factory=list)
    state: dict = field(default_factory=dict)
    fingerprint: str = ""
    lock: threading.Lock = field(default_factory=threading.Lock, repr=False)

    # --- paths (everything the model needs hangs off `repo` + the config) --
    def _p(self, rel: str | None) -> Path | None:
        return (self.repo / rel) if rel else None

    @property
    def binary(self) -> Path:
        return self.repo / self.cfg.binary

    @property
    def build_dir(self) -> Path:
        return self.repo / self.cfg.build_dir

    @property
    def template(self) -> Path:
        return self.repo / self.cfg.template

    @property
    def drug_meta(self) -> Path | None:
        return self._p(self.cfg.drug_metadata)

    @property
    def healthy(self) -> Path | None:
        return self._p(self.cfg.healthy_state)

    @property
    def submodel(self) -> Path | None:
        return self._p(self.cfg.submodel_priors)

    @property
    def priors_csv(self) -> Path:
        return self.repo / self.cfg.priors_csv

    @property
    def model_structure(self) -> Path:
        return self.repo / self.cfg.model_structure

    @property
    def ode_system(self) -> Path | None:
        return self._p(self.cfg.ode_system)

    def theta_flavors(self) -> list[dict]:
        """Which parameter vectors this branch can be run at.

        θ is a SEPARATE axis from the equations: a commit can move the model by
        changing the ODEs, or by moving where θ sits. Each flavor is a layer stacked
        on the template defaults. Declared in the model's ``[views] theta_flavors``.
        """
        return self.cfg.theta_flavors

    def scenarios(self) -> dict:
        out = {}
        for s in self.cfg.scenarios:
            out[s.id] = {
                "yaml": self.repo / s.yaml,
                "t_end": float(s.t_end),
                "target_dirs": [self.repo / d for d in s.target_dirs],
            }
        return out

    # --- build ------------------------------------------------------------
    @property
    def built(self) -> bool:
        return self.binary.exists()

    @property
    def stale(self) -> bool:
        """Binary older than the generated ODE source it should have been built
        from. This is exactly the failure that had the explorer serving a
        9-day-old model: a binary can exist and still be the wrong model. With no
        ``ode_system`` declared (prebuilt models), only presence is checked."""
        if not self.binary.exists():
            return True
        ode = self.ode_system
        if ode is None or not ode.exists():
            return False
        return self.binary.stat().st_mtime < ode.stat().st_mtime

    def say(self, line: str) -> None:
        self.log.append(line.rstrip())
        if len(self.log) > MAX_LOG:
            del self.log[: len(self.log) - MAX_LOG]

    def label(self) -> str:
        name = self.branch if self.branch != "detached" else self.sha
        tag = " *" if self.dirty else ""
        return f"{name}{tag}"

    def public(self) -> dict:
        return {"id": self.id, "kind": self.kind, "branch": self.branch,
                "sha": self.sha, "subject": self.subject, "dirty": self.dirty,
                "path": str(self.repo), "status": self.status,
                "message": self.message, "progress": self.progress,
                "built": self.built, "stale": self.stale,
                "label": self.label(), "ref": self.ref}


def _parse_worktrees(porcelain: str) -> list[dict]:
    blocks, cur = [], {}
    for line in porcelain.splitlines():
        if not line.strip():
            if cur:
                blocks.append(cur)
                cur = {}
            continue
        key, _, val = line.partition(" ")
        cur[key] = val
    if cur:
        blocks.append(cur)
    return blocks


def discover_worktrees(home: Path, cfg: ExplorerConfig,
                       scratch_root: Path | None = None) -> list[ModelSource]:
    """Every existing git worktree is a ready-made source.

    Scratch worktrees we created for a ref are git worktrees too, so they come back
    through this scan after a restart. Classify them by location, or they'd be
    indistinguishable from the user's real worktrees and refuse to be removed.
    """
    out = []
    for b in _parse_worktrees(git(["worktree", "list", "--porcelain"], home)):
        path = Path(b.get("worktree", ""))
        if not path.exists():
            continue
        is_scratch = scratch_root is not None and path.parent == scratch_root
        branch = "detached"
        if "branch" in b:
            branch = b["branch"].replace("refs/heads/", "")
        sha = (b.get("HEAD") or "")[:7]
        subject = ""
        try:
            subject = git(["log", "-1", "--format=%s", "HEAD"], path).strip()
        except subprocess.CalledProcessError:
            pass
        dirty = bool(git(["status", "--porcelain"], path).strip())
        out.append(ModelSource(id=path.name, repo=path,
                               kind="ref" if is_scratch else "worktree",
                               branch=branch, sha=sha, subject=subject,
                               cfg=cfg, dirty=dirty))
    return out


def list_refs(home: Path) -> list[dict]:
    """Branches (local + remote) the user can point the explorer at."""
    fmt = "%(refname:short)%09%(objectname:short)%09%(committerdate:relative)%09%(contents:subject)"
    out = []
    raw = git(["for-each-ref", f"--format={fmt}", "--sort=-committerdate",
               "refs/heads", "refs/remotes"], home)
    for line in raw.splitlines():
        parts = line.split("\t")
        if len(parts) < 4 or parts[0].endswith("/HEAD"):
            continue
        out.append({"ref": parts[0], "sha": parts[1], "when": parts[2],
                    "subject": parts[3]})
    return out


def resolve_ref(home: Path, ref: str) -> tuple[str, str]:
    """(full sha, subject) — raises CalledProcessError on an unknown ref."""
    sha = git(["rev-parse", ref], home).strip()
    subject = git(["log", "-1", "--format=%s", sha], home).strip()
    return sha, subject


def find_dep_donor(sources: list[ModelSource]) -> Path | None:
    """A checkout with populated SUNDIALS/yaml-cpp sources, to seed from."""
    for s in sources:
        if (s.build_dir / "_deps" / "sundials-src").is_dir():
            return s.build_dir / "_deps"
    return None


def add_ref_worktree(home: Path, scratch_root: Path, ref: str,
                     cfg: ExplorerConfig) -> ModelSource:
    """Materialise `ref` as a detached worktree. Detached, not a branch checkout:
    git refuses a branch already checked out elsewhere, and we only need the tree."""
    sha, subject = resolve_ref(home, ref)
    short = sha[:7]
    path = scratch_root / f"{re.sub(r'[^A-Za-z0-9._-]', '-', ref)}-{short}"
    if not path.exists():
        scratch_root.mkdir(parents=True, exist_ok=True)
        git(["worktree", "add", "--detach", str(path), sha], home)
    return ModelSource(id=path.name, repo=path, kind="ref", branch="detached",
                       sha=short, subject=subject, cfg=cfg, dirty=False, ref=ref)


def remove_ref_worktree(home: Path, src: ModelSource) -> None:
    if src.kind != "ref":
        raise ValueError("refusing to remove a non-scratch worktree")
    git(["worktree", "remove", "--force", str(src.repo)], home)


def seed_deps(src: ModelSource, donor: Path | None) -> None:
    """Copy only `_deps/*-src` from a donor checkout. The `*-build`/`*-subbuild`
    dirs hard-code the donor's absolute path in their CMakeCache and make cmake
    bail; the src copy is what actually skips the slow part (the network clone).
    rsync, not cp -r: sundials-src/doc has relative symlinks pointing outside the
    copied subtree, which cp follows and dies on."""
    deps = src.build_dir / "_deps"
    if donor is None or (deps / "sundials-src").is_dir():
        return
    deps.mkdir(parents=True, exist_ok=True)
    for srcdir in sorted(donor.glob("*-src")):
        if not srcdir.is_dir():
            continue
        src.say(f"seeding {srcdir.name} ...")
        subprocess.run(["rsync", "-a", str(srcdir), str(deps) + "/"],
                       capture_output=True, text=True)


def build(src: ModelSource, python_exe: Path, donor: Path | None,
          jobs: int = 4) -> None:
    """Configure + compile the simulator into the source's OWN build dir, per the
    model's ``[build]`` config.

    Never runs codegen here — the generated ODE source / param template /
    model_structure.json are git-tracked and already current for whatever ref is
    checked out. A model with ``cmake_source`` unset ships a prebuilt binary and is
    never compiled."""
    mc = src.cfg
    if not mc.cmake_source:
        raise RuntimeError(
            f"binary missing at {src.binary} and this model declares no "
            "[build] cmake_source (prebuilt-only) — build it out of band")
    src.status = "building"
    src.progress = 0
    if mc.seed_deps:
        seed_deps(src, donor)

    src_dir = str(src.repo / mc.cmake_source)
    build_dir = str(src.build_dir)

    # A cache from another checkout (copied binary, seeded deps) poisons configure.
    cache = src.build_dir / "CMakeCache.txt"
    if cache.exists():
        txt = cache.read_text(errors="ignore")
        if str(src.build_dir) not in txt:
            src.say("stale CMakeCache from another checkout — discarding")
            cache.unlink()
            shutil.rmtree(src.build_dir / "CMakeFiles", ignore_errors=True)

    cfgcmd = ["cmake", "-S", src_dir, "-B", build_dir,
              "-DCMAKE_BUILD_TYPE=Release",
              "-DCMAKE_CXX_COMPILER_LAUNCHER=ccache",
              "-DCMAKE_C_COMPILER_LAUNCHER=ccache",
              f"-DPython3_EXECUTABLE={python_exe}"]
    _run_streaming(src, cfgcmd, "configure")
    bld = ["cmake", "--build", build_dir, "-j", str(jobs)]
    if mc.build_target:
        bld += ["--target", mc.build_target]
    _run_streaming(src, bld, "build")
    if not src.binary.exists():
        raise RuntimeError("build finished but the simulator binary is missing")
    src.progress = 100


def _run_streaming(src: ModelSource, cmd: list[str], phase: str) -> None:
    src.say(f"$ {' '.join(cmd)}")
    proc = subprocess.Popen(cmd, cwd=str(src.repo), stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, text=True, bufsize=1)
    tail: list[str] = []
    for line in proc.stdout:  # type: ignore[union-attr]
        line = line.rstrip()
        tail.append(line)
        if len(tail) > 40:
            tail.pop(0)
        m = _PCT_RE.match(line)
        if m:
            src.progress = int(m.group(1))
        if m or "error" in line.lower() or line.startswith("--"):
            src.say(line)
    proc.wait()
    if proc.returncode != 0:
        for line in tail:
            src.say(line)
        raise RuntimeError(f"{phase} failed (rc={proc.returncode}) — see log")
