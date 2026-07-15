"""The model contract — everything that was hardcoded for one model, in a file.

A model repo ships an ``explorer.toml`` (at its root or under ``.model_explorer/``)
that tells the explorer where the simulator binary lives, how to build it, which
data files describe the model, and which scenarios to offer. Optional ``[views]``
curate the default slider set / module lenses / trajectory panel; when omitted the
server falls back to auto-generated defaults so a bare model still loads.

Paths in the config are interpreted relative to the *model checkout* (so they hold
across branches/worktrees), EXCEPT ``glossary``/``interventions``, which are
UI-content files interpreted relative to the config file itself.

The config is the Level-1 decoupling boundary: it captures the QSPIO/``qsp_sim``
framework's layout without baking one model's names into the code. A future
non-``qsp_sim`` backend would slot in below this.
"""
from __future__ import annotations

import tomllib
from dataclasses import dataclass, field
from pathlib import Path

CONFIG_NAMES = ("explorer.toml", ".model_explorer/explorer.toml")


@dataclass
class Scenario:
    id: str
    yaml: str
    t_end: float
    target_dirs: list[str] = field(default_factory=list)


@dataclass
class ExplorerConfig:
    """Parsed ``explorer.toml`` — all paths are repo-relative strings."""

    name: str
    config_path: Path

    # --- build / simulator -------------------------------------------------
    binary: str
    build_dir: str
    ode_system: str | None
    cmake_source: str | None          # None => prebuilt binary, never compile
    build_target: str
    seed_deps: bool

    # --- model data files --------------------------------------------------
    template: str
    priors_csv: str
    model_structure: str
    drug_metadata: str | None
    healthy_state: str | None
    submodel_priors: str | None

    # --- scenarios + sim invocation ---------------------------------------
    scenarios: list[Scenario]
    sim_extra_args: list[str]

    # --- optional curated UI ----------------------------------------------
    theta_flavors: list[dict]
    curated: list[str]
    modules: dict[str, list[str]]
    traj_species: list[str]
    glossary_path: Path | None
    interventions_path: Path | None


def find_config(repo: Path, explicit: Path | None = None) -> Path:
    """Locate the explorer.toml for a model checkout."""
    if explicit is not None:
        p = Path(explicit).expanduser().resolve()
        if not p.is_file():
            raise FileNotFoundError(f"config not found: {p}")
        return p
    for name in CONFIG_NAMES:
        cand = repo / name
        if cand.is_file():
            return cand.resolve()
    raise FileNotFoundError(
        f"no explorer.toml under {repo} (looked for {', '.join(CONFIG_NAMES)}). "
        "Pass --config PATH or add one; see the example in the qsp-model-explorer README."
    )


# Fallback θ flavors when the config omits [[theta_flavors]]. "template" always
# works (raw defaults); the CSV/submodel overlays are only meaningful if those
# files are declared, so they're pruned in load_config accordingly.
_DEFAULT_FLAVORS = [
    {"id": "template", "label": "template defaults", "desc": "raw param template — no overlay"},
    {"id": "csv", "label": "CSV prior medians", "desc": "priors CSV median column"},
    {"id": "submodel", "label": "submodel medians",
     "desc": "submodel_priors.yaml — the SBI-faithful prior center (default)"},
]


def load_config(config_path: Path) -> ExplorerConfig:
    doc = tomllib.loads(config_path.read_text())
    cfg_dir = config_path.parent

    model = doc.get("model", {})
    build = doc.get("build", {})
    paths = doc.get("paths", {})
    sim = doc.get("sim", {})
    views = doc.get("views", {})

    def req(section: dict, key: str, where: str) -> str:
        if key not in section:
            raise KeyError(f"explorer.toml [{where}] missing required key '{key}'")
        return str(section[key])

    scenarios = []
    for s in doc.get("scenarios", []):
        scenarios.append(Scenario(
            id=req(s, "id", "[[scenarios]]"),
            yaml=req(s, "yaml", "[[scenarios]]"),
            t_end=float(s.get("t_end", 1.0)),
            target_dirs=[str(d) for d in s.get("target_dirs", [])],
        ))
    if not scenarios:
        raise ValueError("explorer.toml declares no [[scenarios]]")

    submodel_priors = paths.get("submodel_priors")
    priors_csv = req(paths, "priors_csv", "paths")

    # Prune θ flavors whose backing file is absent.
    flavors = views.get("theta_flavors") or _DEFAULT_FLAVORS
    have = {"template"}
    if priors_csv:
        have.add("csv")
    if submodel_priors:
        have.add("submodel")
    flavors = [f for f in flavors if f.get("id") in have]

    def content_path(key: str) -> Path | None:
        v = views.get(key)
        return (cfg_dir / v).resolve() if v else None

    modules = {k: list(v) for k, v in (views.get("modules") or {}).items()}

    return ExplorerConfig(
        name=str(model.get("name", config_path.parent.name)),
        config_path=config_path,
        binary=req(build, "binary", "build"),
        build_dir=str(build.get("build_dir", "build")),
        ode_system=build.get("ode_system"),
        cmake_source=build.get("cmake_source"),
        build_target=str(build.get("target", "")),
        seed_deps=bool(build.get("seed_deps", False)),
        template=req(paths, "template", "paths"),
        priors_csv=priors_csv,
        model_structure=req(paths, "model_structure", "paths"),
        drug_metadata=paths.get("drug_metadata"),
        healthy_state=paths.get("healthy_state"),
        submodel_priors=submodel_priors,
        scenarios=scenarios,
        sim_extra_args=[str(a) for a in sim.get("extra_args", [])],
        theta_flavors=flavors,
        curated=list(views.get("curated", [])),
        modules=modules,
        traj_species=list(views.get("traj_species", [])),
        glossary_path=content_path("glossary"),
        interventions_path=content_path("interventions"),
    )
