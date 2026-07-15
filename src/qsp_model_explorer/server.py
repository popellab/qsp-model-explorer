#!/usr/bin/env python
"""Interactive QSP model explorer — local backend.

Wraps the compiled ``qsp_sim`` binary in a tiny stdlib HTTP server so a
browser frontend can slide parameters and watch calibration-target
readouts + Hill-function saturation update in near-real-time (a single
sim is ~0.4 s). Results are cached by (source, scenario, overrides) so
revisiting a parameter state is instant.

The model contract lives in a config
------------------------------------
The model is served from a checkout that ships an ``explorer.toml`` (see
``config.py``) declaring the simulator binary, its build recipe, the model data
files, and the scenarios. Nothing about one model's file layout or species names
is baked into this code — point ``--repo`` at any QSPIO/``qsp_sim`` checkout that
carries an ``explorer.toml``.

Multiple model sources
----------------------
When the model is fully described by git-tracked files (the generated
``ODE_system.cpp`` / param template / ``model_structure.json`` committed alongside
the sources), the explorer can point at ANY branch: check it out, compile
(~10-20 s, ccache warm), load. See ``sources.py``.

Everything model-derived hangs off a ``ModelSource`` rather than a module-level
``REPO`` constant. Two sources can be live at once, which is the point: set the
delta-baseline to branch A and the current run to branch B and the diagram shows
exactly what a commit did.

Design notes / provenance:
- Baseline operating point overlays the model's ``submodel_priors.yaml`` medians
  onto its param template via ``submodel_medians.load_submodel_medians`` (handles
  lognormal / gamma / invgamma / normal correctly — a naive exp(mu) blows up the
  gamma-distributed rates).
- A target's observable maps CSV columns -> a value: the CSV columns (species,
  compartment volumes, aggregates, plus the H_* Hill rules) ARE the
  ``species_dict`` keys. We reuse
  ``qsp_hpc.calibration.yaml_loader.load_calibration_targets`` to compile each
  target's ``compute_test_statistic(time, species_dict)`` wrapper.

Run:  qsp-model-explorer --repo /path/to/model [--port 8765]
  or  python -m qsp_model_explorer.server --repo /path/to/model
Then open http://localhost:8765/
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import subprocess
import sys
import threading
import time
import traceback
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs, urlparse

import pandas as pd
import yaml

from qsp_model_explorer import sources as SRC
from qsp_model_explorer.config import ExplorerConfig, find_config, load_config
from qsp_model_explorer.hill_deps import HillDeps
from qsp_model_explorer.submodel_medians import load_submodel_medians
# qsp-hpc-tools provides the calibration-target loader + param-XML renderer for the
# QSPIO/qsp_sim framework. It is a required dependency (see pyproject).
from qsp_hpc.calibration.yaml_loader import load_calibration_targets
from qsp_hpc.cpp.param_xml import ParamXMLRenderer

STATIC = Path(__file__).resolve().parent / "static"

# --- Runtime config (set in main() from --repo / --config) ---------------
# HOME is the MODEL checkout being served (not the tool's own location); CFG is its
# parsed explorer.toml. CACHE_ROOT / SCRATCH_ROOT hang off HOME.
HOME: Path = None            # type: ignore[assignment]
CFG: ExplorerConfig = None   # type: ignore[assignment]
CACHE_ROOT: Path = None      # type: ignore[assignment]
SCRATCH_ROOT: Path = None    # type: ignore[assignment]

# Whitelisted math functions for rate-law evaluation (flux). Rate laws are
# model-generated (model_structure.json), not user input.
_FLUX_FUNCS = {"__builtins__": {}, "log": math.log, "log10": math.log10,
               "exp": math.exp, "sqrt": math.sqrt, "pow": pow,
               "max": max, "min": min, "abs": abs}

# Hill regulators in a rate law: H_<name> columns whose driver species are hidden
# behind the column (recovered via hill_deps from the generated C++).
HILL_RE = re.compile(r"\bH_[A-Za-z0-9_]+\b")


def _safe_expr(rate_law: str) -> str:
    """Rewrite a rate law so it is eval-safe: `^`->`**` and dotted species
    names (V_T.C1) -> underscore (V_T__C1) so they parse as plain identifiers."""
    expr = rate_law.replace("^", "**")
    return re.sub(r"[A-Za-z_][A-Za-z0-9_.]*",
                  lambda m: m.group(0).replace(".", "__"), expr)


# The curated slider set, module "lenses", and trajectory small-multiples are all
# MODEL content — a model supplies them via its explorer.toml [views]. When a model
# omits them, load_model() auto-generates sensible fallbacks from the live species,
# so a bare model still loads. Species/param names in [views] are filtered against
# the live model at load, so stale entries drop silently.
N_AUTO_CURATED = 12   # fallback slider count when [views].curated is empty
N_AUTO_TRAJ = 12      # fallback trajectory panel count

# --- Source registry ----------------------------------------------------
SOURCES: dict[str, SRC.ModelSource] = {}
REGISTRY_LOCK = threading.Lock()
DEFAULT_SRC = ""


def refresh_registry() -> None:
    """Re-scan git worktrees. Keeps already-loaded sources (and their state)."""
    with REGISTRY_LOCK:
        for s in SRC.discover_sources(HOME, CFG, SCRATCH_ROOT):
            old = SOURCES.get(s.id)
            if old is None:
                SOURCES[s.id] = s
            else:
                # refresh the git facts; keep loaded model state
                old.branch, old.sha = s.branch, s.sha
                old.subject, old.dirty = s.subject, s.dirty


def get_source(sid: str | None) -> SRC.ModelSource:
    sid = sid or DEFAULT_SRC
    src = SOURCES.get(sid)
    if src is None:
        raise KeyError(f"unknown source {sid}")
    return src


# --- Per-source load ----------------------------------------------------
# --- θ: which parameter VALUES to run the equations at ---------------------
# θ is a separate axis from the equations. A commit moves the model either by changing
# the ODEs or by moving where θ sits (a submodel-prior refresh moves medians), and the
# diagram only ever showed the first. Making θ selectable — and selectable from a
# DIFFERENT branch than the equations — is what lets you set that split instead of
# inferring it.
THETA_CACHE: dict[tuple, dict] = {}
THETA_LOCK = threading.Lock()


def _csv_medians(src: SRC.ModelSource) -> tuple[dict, dict]:
    df = pd.read_csv(src.priors_csv)
    vals, marg = {}, {}
    for _, r in df.iterrows():
        vals[r["name"]] = float(r["median"])
        marg[r["name"]] = {"distribution": str(r["distribution"]),
                           "sigma": float(r["dist_param2"])}
    return vals, marg


def _yaml_marginals(path: Path) -> dict:
    """Per-param sigma/distribution from a submodel_priors-shaped yaml."""
    doc = yaml.safe_load(path.read_text())
    out = {}
    for p in doc.get("parameters", []):
        m = p.get("marginal", {})
        out[p["name"]] = {"distribution": m.get("distribution", "lognormal"),
                          "sigma": float(m.get("sigma", 0.0))}
    return out


def load_theta(src: SRC.ModelSource, tsrc: SRC.ModelSource, flavor: str) -> dict:
    """The effective parameter vector for running `src`'s equations.

    Values come from `tsrc` (possibly a DIFFERENT branch), matched BY NAME onto
    `src`'s template. A param present in tsrc but absent from src's template is
    dropped and reported — that's the honest behaviour when the param set differs
    across commits.
    """
    key = (src.id, tsrc.id, flavor)
    with THETA_LOCK:
        if key in THETA_CACHE:
            return THETA_CACHE[key]

    # The template must be the EQUATIONS' own — it has to match the binary. Only the
    # VALUES come from tsrc.
    renderer = ParamXMLRenderer(src.template)
    defaults = dict(renderer.template_defaults)
    xml = src.template.read_text()

    vals: dict[str, float] = {}
    marg: dict[str, dict] = {}
    if flavor == "csv":
        vals, marg = _csv_medians(tsrc)
    elif flavor == "submodel":
        if tsrc.submodel is None:
            raise RuntimeError("submodel θ flavor requested but the model declares "
                               "no [paths] submodel_priors")
        sm, _ = load_submodel_medians(tsrc.submodel)
        vals.update(sm)
        marg.update(_yaml_marginals(tsrc.submodel))
    # flavor == "template": no overlay

    applied, dropped = {}, []
    for name, v in vals.items():
        if name not in defaults:
            dropped.append(name)          # param doesn't exist in these equations
            continue
        xml, n = re.subn(rf"<{re.escape(name)}>.*?</{re.escape(name)}>",
                         f"<{name}>{v}</{name}>", xml)
        if n:
            applied[name] = float(v)

    params = {**defaults, **applied}
    theta = {"key": key, "flavor": flavor, "theta_src": tsrc.id,
             "xml": xml, "params": params, "marginals": marg,
             "n_applied": len(applied), "dropped": sorted(dropped),
             "label": f"{tsrc.label()} · {flavor}"}
    src.say(f"θ [{flavor} from {tsrc.id}]: applied {len(applied)}/{len(vals)}"
            + (f", dropped {len(dropped)} not in these equations" if dropped else ""))
    with THETA_LOCK:
        THETA_CACHE[key] = theta
    return theta


def load_slider_meta(src: SRC.ModelSource, theta: dict) -> dict:
    """Per-parameter slider center + [lo, hi] range + units, AT THE CHOSEN θ.

    Center = the value θ actually runs at (so the slider starts at the operating
    point, whichever flavor you picked). Range for lognormals = exp(mu ± 2σ), σ taken
    from θ's own marginal when it has one (an SBI posterior σ is far tighter than the
    prior's — which is the point) and the CSV prior σ otherwise.
    """
    df = pd.read_csv(src.priors_csv)
    marg = theta["marginals"]
    params = theta["params"]

    out = {}
    for _, row in df.iterrows():
        name = row["name"]
        units = str(row.get("units", "") or "")
        center = float(params.get(name, row["median"]))
        m = marg.get(name, {})
        dist = str(m.get("distribution") or row["distribution"])
        log = dist == "lognormal"
        if log and center > 0:
            sigma = float(m.get("sigma") or row["dist_param2"])
            mu = math.log(center)
            lo, hi = math.exp(mu - 2 * sigma), math.exp(mu + 2 * sigma)
        else:
            lo, hi = center / 5.0, center * 5.0
        out[name] = {"name": name, "center": center, "lo": lo, "hi": hi,
                     "units": units, "log": log, "submodel": name in marg}
    return out


def load_targets(src: SRC.ModelSource) -> dict:
    """Compile compute_test_statistic per target, keyed by scenario."""
    by_scenario = {}
    for scen, cfg in src.state["scenarios"].items():
        dirs = [d for d in cfg["target_dirs"] if d.exists()]
        if not dirs:
            by_scenario[scen] = []   # a scenario may declare no calibration targets
            src.say(f"scenario {scen}: 0 targets (no target dirs)")
            continue
        tdf = load_calibration_targets(dirs)
        compiled = []
        for _, r in tdf.iterrows():
            ns: dict = {}
            try:
                exec(r["model_output_code"], ns)
                fn = ns["compute_test_statistic"]
            except Exception as e:  # noqa: BLE001
                src.say(f"target {r['test_statistic_id']} failed to compile: {e}")
                continue
            compiled.append({
                "id": r["test_statistic_id"],
                "fn": fn,
                "median": float(r["median"]),
                "ci_lo": float(r["ci95_lower"]),
                "ci_hi": float(r["ci95_upper"]),
                "units": str(r["units"] or ""),
            })
        by_scenario[scen] = compiled
        src.say(f"scenario {scen}: {len(compiled)} targets")
    return by_scenario


def load_graph(src: SRC.ModelSource) -> dict:
    ms = json.loads(src.model_structure.read_text())
    compartments = [{"name": c["name"], "description": c.get("description", "")}
                    for c in ms.get("compartments", [])]
    species = [{"name": s["name"], "compartment": s.get("compartment", ""),
                "base_name": s.get("base_name", s["name"]),
                "units": s.get("units", ""), "description": s.get("description", "")}
               for s in ms.get("species", [])]
    reactions = []
    for rx in ms.get("reactions", []):
        reactions.append({"name": rx.get("name", ""),
                          "reactants": rx.get("reactants", []),
                          "products": rx.get("products", []),
                          "rate_law": rx.get("rate_law", ""),
                          "parameters": rx.get("parameters", [])})
    # Hill parameters (static descriptions; live values come from sim rules).
    hills = [{"name": p["name"], "description": p.get("description", "")}
             for p in ms.get("parameters", []) if p["name"].startswith("H_")]
    return {"compartments": compartments, "species": species,
            "reactions": reactions, "hills": hills}


def compute_fluxes(src: SRC.ModelSource, final_row: dict, eff_params: dict) -> dict:
    """Evaluate every reaction's rate law at the endpoint. Returns
    {reaction_name: flux or None}. Namespace = params overlaid with the
    (time-varying) sim columns so V_T, H_* etc. take their live values."""
    ns = {k.replace(".", "__"): v for k, v in eff_params.items()}
    ns.update({k.replace(".", "__"): v for k, v in final_row.items()})
    out = {}
    for name, expr in src.state["rx_safe"]:
        try:
            out[name] = float(eval(expr, _FLUX_FUNCS, ns))  # noqa: S307
        except Exception:  # noqa: BLE001
            out[name] = None
    return out


def compute_signs(src: SRC.ModelSource, row: dict, eff_params: dict,
                  hill_driver_signs: dict | None = None):
    """Regulation signs + influencer species per reaction at the baseline
    operating point. Returns (signs, influencers).

    Two paths, both keyed by the real reaction name:
    - Direct: local sign of d(flux)/d(species) for each species that appears in
      the rate law (perturb ~1%; small additive bump when 0). Products rise with
      flux, so this IS the regulator→product sign: +1 activates, -1 inhibits.
    - Hill-mediated: for each H_* the rate law references, compose d(flux)/d(Hill)
      (perturb the Hill column here) with d(Hill)/d(species) (from hill_deps, the
      generated C++ Hill formula). This recovers cytokine→target edges (IL-10 ⊣
      maturation, TGF-β → Treg, …) whose driver is hidden inside the Hill column.

    Local at the diagnosis operating point — a Hill that swaps numerator/denominator
    regime elsewhere could differ, but this reliably surfaces a backwards arrow."""
    hill_driver_signs = hill_driver_signs or {}
    species_set = {s["name"] for s in src.state["graph"]["species"]}
    ns0 = {k.replace(".", "__"): v for k, v in eff_params.items()}
    ns0.update({k.replace(".", "__"): v for k, v in row.items()})
    signs: dict = {}
    influencers: dict = {}

    def _sign(expr, key, base, rel):
        old = ns0.get(key)
        if old is None:
            return None
        ns0[key] = old + (abs(old) * rel if old != 0 else (1e-3 if rel > 0.01 else 1e-6))
        try:
            up = float(eval(expr, _FLUX_FUNCS, ns0))  # noqa: S307
        except Exception:  # noqa: BLE001
            up = None
        ns0[key] = old
        if up is None:
            return None
        d = up - base
        thresh = max(abs(base) * 1e-6, 1e-30)
        return 1 if d > thresh else (-1 if d < -thresh else 0)

    for rx in src.state["graph"]["reactions"]:
        expr = _safe_expr(rx["rate_law"])
        try:
            base = float(eval(expr, _FLUX_FUNCS, ns0))  # noqa: S307
        except Exception:  # noqa: BLE001
            continue
        rxout = {}
        # direct species in the rate law
        for sp in {t for t in re.findall(r"[A-Za-z_][A-Za-z0-9_.]*", rx["rate_law"])
                   if t in species_set}:
            s = _sign(expr, sp.replace(".", "__"), base, 0.01)
            if s is not None:
                rxout[sp] = s
        inf = set(rxout.keys())
        # Hill-mediated drivers: d(flux)/d(Hill) × d(Hill)/d(species)
        for h in set(HILL_RE.findall(rx["rate_law"])):
            drivers = hill_driver_signs.get(h)
            if not drivers:
                continue
            s_dfdh = _sign(expr, h.replace(".", "__"), base, 0.02)
            if not s_dfdh:  # None or 0
                continue
            for sp, sds in drivers.items():
                if sp not in species_set:
                    continue
                inf.add(sp)
                rxout.setdefault(sp, s_dfdh * sds)  # don't override a direct sign
        if rxout:
            signs[rx["name"]] = rxout
        if inf:
            influencers[rx["name"]] = sorted(inf)
    return signs, influencers


def build_signs(src: SRC.ModelSource) -> dict:
    """Run one baseline sim; compute regulation signs (direct + Hill-mediated) and
    attach per-reaction influencer species lists for the frontend edge builder.

    Computed at the source's DEFAULT θ. Signs are a structural property of the rate
    laws (a positive unit factor can't flip a monotonic sign), so they don't need
    recomputing when you switch θ flavor."""
    try:
        theta = src.state["theta"]
        # Signs are structural (rate-law based), so any scenario's operating point
        # works — use the first declared scenario rather than a fixed name.
        scen = next(iter(src.state["scenarios"]))
        base = run_sim(src, scen, {}, theta)
        if base.get("error") or not base.get("species"):
            src.say(f"sign baseline sim failed: {base.get('error')}")
            return {}
        hd = src.state.get("hill_deps")
        driver_signs = hd.driver_signs(base["species"], theta["params"]) if hd else {}
        signs, influencers = compute_signs(src, base["species"],
                                           theta["params"], driver_signs)
        for rx in src.state["graph"]["reactions"]:
            rx["influencers"] = influencers.get(rx["name"], [])
        n_hill = sum(len(v) for v in driver_signs.values())
        src.say(f"hill drivers: {len(driver_signs)} hills, {n_hill} (hill,species) "
                f"pairs {'(C++ parse ok)' if hd and hd.ok else '(disabled)'}")
        return signs
    except Exception as e:  # noqa: BLE001
        src.say(f"sign computation failed: {e}")
        return {}


def load_model(src: SRC.ModelSource) -> None:
    """Everything the old module-level startup() did, but for one source."""
    src.status = "loading"
    st = src.state
    st.clear()
    st["scenarios"] = src.scenarios()
    st["theta_flavors"] = src.theta_flavors()
    # Default θ = this branch's own richest overlay: submodel medians if the model
    # ships them, else CSV medians, else the raw template.
    flavor_ids = {f["id"] for f in st["theta_flavors"]}
    default_flavor = next((f for f in ("submodel", "csv", "template")
                           if f in flavor_ids), "template")
    theta = load_theta(src, src, default_flavor)
    st["theta"] = theta
    st["base_params"] = theta["params"]

    # Fingerprint the *equations*: binary bytes + scenario definitions. θ is NOT
    # folded in here — it varies per request, so it's folded into the cache key
    # instead (see run_sim). This is what stops branch A's cached run from being
    # served for branch B, or θ_A's run for θ_B.
    h = hashlib.sha256()
    h.update(src.binary.read_bytes())
    for cfg in st["scenarios"].values():
        if cfg["yaml"].exists():
            h.update(cfg["yaml"].read_bytes())
    src.fingerprint = h.hexdigest()[:12]

    st["sliders"] = load_slider_meta(src, theta)
    st["targets"] = load_targets(src)
    st["graph"] = load_graph(src)
    sp = {s["name"] for s in st["graph"]["species"]}
    species_order = [s["name"] for s in st["graph"]["species"]]

    # Curated slider set: from [views].curated, filtered to live params; else the
    # first N params in the table.
    if CFG.curated:
        st["curated"] = [p for p in CFG.curated if p in st["sliders"]]
    else:
        st["curated"] = list(st["sliders"].keys())[:N_AUTO_CURATED]

    # Module lenses: from [views].modules, filtered to live species; else a single
    # "all" lens over every species.
    if CFG.modules:
        st["modules"] = {m: [s for s in lst if s in sp]
                         for m, lst in CFG.modules.items()}
    else:
        st["modules"] = {"all": species_order}

    # Trajectory panel: from [views].traj_species, filtered; else the first N species.
    if CFG.traj_species:
        st["traj_species"] = [s for s in CFG.traj_species if s in sp]
    else:
        st["traj_species"] = species_order[:N_AUTO_TRAJ]

    st["rx_safe"] = [(rx["name"], _safe_expr(rx["rate_law"]))
                     for rx in st["graph"]["reactions"]]
    st["glossary"] = (json.loads(CFG.glossary_path.read_text())
                      if CFG.glossary_path and CFG.glossary_path.exists() else {})
    inter = (json.loads(CFG.interventions_path.read_text()).get("interventions", [])
             if CFG.interventions_path and CFG.interventions_path.exists() else [])
    # keep only interventions whose params are all in the template (overrides that
    # miss a param would raise at sim time) — the param set differs across branches
    tmpl = set(st["base_params"])
    st["interventions"] = [iv for iv in inter
                           if all(p in tmpl for p in iv.get("overrides", {}))]
    dropped = len(inter) - len(st["interventions"])
    if dropped:
        src.say(f"interventions: dropped {dropped} with params not in template")
    st["hill_deps"] = HillDeps(src.ode_system,
                               [s["name"] for s in st["graph"]["species"]],
                               list(st["base_params"].keys()))
    st["signs"] = build_signs(src)  # needs rx_safe + targets + hill_deps above
    n_sign_edges = sum(len(v) for v in st["signs"].values())
    src.say(f"regulation signs: {len(st['signs'])} reactions, "
            f"{n_sign_edges} (reaction,species) pairs")
    src.say(f"ready. {len(st['sliders'])} sliders, {len(st['graph']['species'])} species, "
            f"{len(st['graph']['reactions'])} reactions, {len(st['modules'])} modules")
    src.status = "ready"
    src.message = ""


def activate(src: SRC.ModelSource, force_build: bool = False) -> None:
    """Build (if needed) then load. Safe to call concurrently; serialised per source."""
    with src.lock:
        if src.status == "ready" and not force_build:
            return
        try:
            src.log.clear()
            if force_build or src.stale:
                why = "forced" if force_build else (
                    "no binary" if not src.built else "binary older than ODE_system.cpp")
                src.say(f"build needed ({why})")
                refresh_registry()
                donor = SRC.find_dep_donor(list(SOURCES.values()))
                SRC.build(src, Path(sys.executable), donor)
                src.say("build ok")
            load_model(src)
        except Exception as e:  # noqa: BLE001
            src.status = "error"
            src.message = str(e)[:400]
            src.say(f"ERROR: {e}")
            src.say(traceback.format_exc()[-800:])


def activate_async(src: SRC.ModelSource, force_build: bool = False) -> None:
    if src.status in ("building", "loading"):
        return
    src.status = "building" if (force_build or src.stale) else "loading"
    threading.Thread(target=activate, args=(src, force_build), daemon=True).start()


# --- Simulation ---------------------------------------------------------
def _cache_key(scenario: str, overrides: dict) -> str:
    payload = json.dumps({"s": scenario, "o": {k: round(float(v), 12)
                          for k, v in sorted(overrides.items())}}, sort_keys=True)
    return hashlib.sha256(payload.encode()).hexdigest()[:16]


def _theta_tag(theta: dict) -> str:
    return hashlib.sha256(
        f"{theta['theta_src']}|{theta['flavor']}".encode()
    ).hexdigest()[:8]


def run_sim(src: SRC.ModelSource, scenario: str, overrides: dict,
            theta: dict) -> dict:
    cfg = src.state["scenarios"][scenario]
    # Cache namespace folds in BOTH the equations (fingerprint = binary + scenarios)
    # and θ, so neither a different branch nor a different parameter vector can ever
    # be served from a stale entry.
    cdir = (CACHE_ROOT / f"{src.id}-{src.fingerprint}-{_theta_tag(theta)}"
            / _cache_key(scenario, overrides))
    csv_path = cdir / "out.csv"
    rules_path = cdir / "rules.txt"
    cache_hit = csv_path.exists() and rules_path.exists()
    t0 = time.time()
    if not cache_hit:
        cdir.mkdir(parents=True, exist_ok=True)
        xml = theta["xml"]
        for name, val in overrides.items():
            xml, n = re.subn(rf"<{re.escape(name)}>.*?</{re.escape(name)}>",
                             f"<{name}>{float(val)!r}</{name}>", xml)
            if not n:
                raise KeyError(f"param not in template: {name}")
        (cdir / "param.xml").write_text(xml)
        cmd = [str(src.binary), "--param", str(cdir / "param.xml"),
               "--scenario", str(cfg["yaml"]),
               "--csv-out", str(csv_path), "--rules-out", str(rules_path),
               "--t-end-days", str(cfg["t_end"])]
        if src.drug_meta is not None:
            cmd += ["--drug-metadata", str(src.drug_meta)]
        if src.healthy is not None:
            cmd += ["--evolve-to-diagnosis", str(src.healthy)]
        cmd += CFG.sim_extra_args
        res = subprocess.run(cmd, capture_output=True, text=True, timeout=180)
        if res.returncode != 0:
            for p in (csv_path, rules_path):
                p.unlink(missing_ok=True)
            return {"error": res.stderr.strip()[-800:] or "qsp_sim failed",
                    "wall_s": round(time.time() - t0, 2), "cache_hit": False}
    wall = time.time() - t0

    df = pd.read_csv(csv_path)
    time_arr = df["Time"].to_numpy(float)
    species_dict = {c: df[c].to_numpy(float) for c in df.columns if c != "Time"}

    # Calibration-target readouts.
    targets_out = []
    for t in src.state["targets"].get(scenario, []):
        try:
            val = float(t["fn"](time_arr, species_dict))
        except Exception as e:  # noqa: BLE001
            targets_out.append({"id": t["id"], "error": str(e)[:120],
                                "median": t["median"], "units": t["units"]})
            continue
        lo, hi, med = t["ci_lo"], t["ci_hi"], t["median"]
        z = None
        if val > 0 and med > 0 and hi > lo > 0:
            sigma_log = (math.log(hi) - math.log(lo)) / (2 * 1.96)
            if sigma_log > 0:
                z = (math.log(val) - math.log(med)) / sigma_log
        in_ci = (not math.isnan(lo)) and lo <= val <= hi
        targets_out.append({"id": t["id"], "value": val, "median": med,
                            "ci_lo": lo, "ci_hi": hi, "z": z, "in_ci": in_ci,
                            "units": t["units"]})
    # Stable alphabetical order so a slider drag doesn't reshuffle the list;
    # the client can optionally re-sort by |z| via a toggle.
    targets_out.sort(key=lambda r: r["id"])

    # Hill-function saturation at the final (diagnosis / endpoint) row.
    hills_out = []
    for c in df.columns:
        if c.startswith("H_"):
            v = float(df[c].iloc[-1])
            hills_out.append({"name": c, "value": v, "saturated": v >= 0.95})
    hills_out.sort(key=lambda r: -r["value"])

    # Full trajectory for the time slider: every column over time, plus each
    # reaction's rate law evaluated at every timestep. Endpoint values (species,
    # fluxes) are just the last index — kept as top-level keys for compatibility.
    eff_params = {**theta["params"], **overrides}

    def _clean(x):
        return None if (x is None or math.isnan(x) or math.isinf(x)) else float(x)

    series = {c: [_clean(x) for x in arr] for c, arr in species_dict.items()}
    n = len(time_arr)
    fluxes_t = {name: [None] * n for name, _ in src.state["rx_safe"]}
    for ti in range(n):
        row = {c: arr[ti] for c, arr in species_dict.items()}
        for name, v in compute_fluxes(src, row, eff_params).items():
            fluxes_t[name][ti] = None if v is None else _clean(v)
    species_vals = {c: series[c][-1] for c in series}
    fluxes = {name: fluxes_t[name][-1] for name in fluxes_t}

    return {"cache_hit": cache_hit, "wall_s": round(wall, 3), "n_rows": len(df),
            "t_end": float(time_arr[-1]) if len(time_arr) else 0.0,
            "source": src.id, "fingerprint": src.fingerprint,
            "theta": {"flavor": theta["flavor"], "theta_src": theta["theta_src"],
                      "label": theta["label"], "n_applied": theta["n_applied"],
                      "dropped": theta["dropped"][:20],
                      "n_dropped": len(theta["dropped"])},
            "times": [float(t) for t in time_arr],
            "series": series, "fluxes_t": fluxes_t,
            "targets": targets_out, "hills": hills_out, "species": species_vals,
            "fluxes": fluxes}


# --- HTTP ---------------------------------------------------------------
def resolve_theta(src: SRC.ModelSource, req: dict) -> dict:
    """θ for a request. Defaults to the source's own submodel medians — i.e. the
    behaviour before θ was selectable. `theta_src` may name ANOTHER branch, which is
    the whole point: branch B's equations at branch A's θ."""
    # Default to the source's own default flavor (submodel if it ships submodel
    # priors, else csv, else template) rather than assuming "submodel".
    flavor = req.get("theta") or src.state["theta"]["flavor"]
    tsid = req.get("theta_src") or src.id
    tsrc = get_source(tsid)
    if tsrc.status != "ready":
        raise RuntimeError(f"θ source {tsid} not loaded")
    return load_theta(src, tsrc, flavor)


class Handler(BaseHTTPRequestHandler):
    def log_message(self, *a):  # quieter
        pass

    def _send(self, code, body, ctype="application/json"):
        data = body if isinstance(body, bytes) else json.dumps(body).encode()
        self.send_response(code)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(data)))
        self.end_headers()
        self.wfile.write(data)

    def _body(self) -> dict:
        n = int(self.headers.get("Content-Length", 0))
        return json.loads(self.rfile.read(n) or b"{}")

    def do_GET(self):
        u = urlparse(self.path)
        path, q = u.path, parse_qs(u.query)
        if path in ("/", "/index.html"):
            return self._send(200, (STATIC / "index.html").read_bytes(),
                              "text/html; charset=utf-8")

        if path == "/sources":
            refresh_registry()
            return self._send(200, {"sources": [s.public() for s in SOURCES.values()],
                                    "default": DEFAULT_SRC})

        if path == "/source":
            try:
                src = get_source((q.get("id") or [None])[0])
            except KeyError as e:
                return self._send(404, {"error": str(e)})
            return self._send(200, {**src.public(), "log": src.log[-60:]})

        if path == "/refs":
            return self._send(200, {"refs": SRC.list_refs(HOME)})

        if path == "/sliders":
            # Just the param table for a source — lets the client show the OTHER
            # branch's centers alongside the current ones without refetching the
            # whole graph (a prior refresh moves medians, and that is invisible in
            # the diagram: it changes where θ sits, not what the equations are).
            try:
                src = get_source((q.get("src") or [None])[0])
            except KeyError as e:
                return self._send(404, {"error": str(e)})
            if src.status != "ready":
                return self._send(409, {"status": src.status, "id": src.id})
            # θ-aware: comparing branch A's params to branch B's is only meaningful at
            # a stated θ, otherwise you're diffing two different flavors by accident.
            try:
                theta = resolve_theta(src, {k: v[0] for k, v in q.items()})
            except (KeyError, RuntimeError) as e:
                return self._send(400, {"error": str(e)})
            sliders = (src.state["sliders"] if theta["key"] == src.state["theta"]["key"]
                       else load_slider_meta(src, theta))
            return self._send(200, {"source": src.id, "sliders": sliders,
                                    "theta": theta["label"]})

        if path == "/meta":
            try:
                src = get_source((q.get("src") or [None])[0])
            except KeyError as e:
                return self._send(404, {"error": str(e)})
            if src.status != "ready":
                # not an error — the client polls /source and shows build progress
                return self._send(409, {"status": src.status, "id": src.id,
                                        "message": src.message})
            st = src.state
            g = dict(st["graph"])
            g["signs"] = st["signs"]
            # Sliders are θ-dependent (the center IS the operating point), so /meta
            # takes the same θ selector as /simulate.
            try:
                theta = resolve_theta(src, {k: v[0] for k, v in q.items()})
            except (KeyError, RuntimeError) as e:
                return self._send(400, {"error": str(e)})
            sliders = (st["sliders"] if theta["key"] == st["theta"]["key"]
                       else load_slider_meta(src, theta))
            return self._send(200, {
                "source": src.public(),
                "scenarios": list(st["scenarios"].keys()),
                "theta_flavors": st["theta_flavors"],
                "theta": {"flavor": theta["flavor"], "theta_src": theta["theta_src"],
                          "label": theta["label"], "n_applied": theta["n_applied"],
                          "n_dropped": len(theta["dropped"]),
                          "dropped": theta["dropped"][:20]},
                "sliders": sliders,
                "curated": st["curated"],
                "modules": st["modules"],
                "default_view": ("Overview" if "Overview" in st["modules"]
                                 else next(iter(st["modules"]), "all")),
                "traj_species": st["traj_species"],
                "glossary": st["glossary"],
                "interventions": st["interventions"],
                "targets_by_scenario": {s: [t["id"] for t in ts]
                                        for s, ts in st["targets"].items()},
                "graph": g,
            })
        return self._send(404, {"error": "not found"})

    def do_POST(self):
        path = urlparse(self.path).path

        if path == "/source/activate":
            req = self._body()
            try:
                src = get_source(req.get("id"))
            except KeyError as e:
                return self._send(404, {"error": str(e)})
            activate_async(src, force_build=bool(req.get("rebuild")))
            return self._send(200, src.public())

        if path == "/source/add":
            req = self._body()
            ref = (req.get("ref") or "").strip()
            if not ref:
                return self._send(400, {"error": "ref required"})
            try:
                src = SRC.add_ref_worktree(HOME, SCRATCH_ROOT, ref, CFG)
            except subprocess.CalledProcessError:
                return self._send(400, {"error": f"unknown ref: {ref}"})
            except Exception as e:  # noqa: BLE001
                return self._send(500, {"error": str(e)[:300]})
            with REGISTRY_LOCK:
                SOURCES.setdefault(src.id, src)
            src = SOURCES[src.id]
            activate_async(src)
            return self._send(200, src.public())

        if path == "/source/remove":
            req = self._body()
            try:
                src = get_source(req.get("id"))
            except KeyError as e:
                return self._send(404, {"error": str(e)})
            if src.kind != "ref":
                return self._send(400, {"error": "only added refs can be removed"})
            try:
                SRC.remove_ref_worktree(HOME, src)
            except Exception as e:  # noqa: BLE001
                return self._send(500, {"error": str(e)[:300]})
            with REGISTRY_LOCK:
                SOURCES.pop(src.id, None)
            return self._send(200, {"removed": src.id})

        if path == "/simulate":
            req = self._body()
            try:
                src = get_source(req.get("src"))
            except KeyError as e:
                return self._send(404, {"error": str(e)})
            if src.status != "ready":
                return self._send(409, {"status": src.status, "id": src.id,
                                        "message": src.message})
            scenario = req.get("scenario", "baseline_no_treatment")
            overrides = {k: float(v) for k, v in (req.get("overrides") or {}).items()}
            if scenario not in src.state["scenarios"]:
                return self._send(400, {"error": f"unknown scenario {scenario}"})
            try:
                theta = resolve_theta(src, req)
            except (KeyError, RuntimeError) as e:
                return self._send(400, {"error": str(e)})
            # A param may not exist on this branch (priors differ across commits);
            # drop it rather than 500, and tell the client what was dropped.
            known = set(theta["params"])
            dropped = sorted(set(overrides) - known)
            overrides = {k: v for k, v in overrides.items() if k in known}
            try:
                out = run_sim(src, scenario, overrides, theta)
            except Exception as e:  # noqa: BLE001
                return self._send(500, {"error": repr(e)})
            if dropped:
                out["dropped_params"] = dropped
            return self._send(200, out)

        return self._send(404, {"error": "not found"})


def run_check(repo: Path, cfg: ExplorerConfig) -> int:
    """Dry-run the config: verify every referenced file exists and actually loads,
    without booting the server or building the binary. Returns a process exit code
    (0 = usable, 1 = at least one hard failure). Optional-but-malformed files are
    warnings; missing required files / unparseable schemas are failures."""
    fails: list[str] = []
    warns: list[str] = []

    def ok(msg: str) -> None:
        print(f"  \033[32m[ok]\033[0m   {msg}")

    def warn(msg: str) -> None:
        warns.append(msg)
        print(f"  \033[33m[warn]\033[0m {msg}")

    def fail(msg: str) -> None:
        fails.append(msg)
        print(f"  \033[31m[FAIL]\033[0m {msg}")

    def P(rel: str) -> Path:
        return repo / rel

    print(f"\nchecking model '{cfg.name}'")
    print(f"config: {cfg.config_path}")
    print(f"repo:   {repo}\n")

    # --- simulator binary + build ----------------------------------------
    binary = P(cfg.binary)
    if binary.exists():
        ok(f"binary present: {cfg.binary}")
    elif cfg.cmake_source:
        if P(cfg.cmake_source).exists():
            warn(f"binary missing ({cfg.binary}) — will build from [build] "
                 f"cmake_source '{cfg.cmake_source}' on launch")
        else:
            fail(f"binary missing AND [build] cmake_source '{cfg.cmake_source}' "
                 "does not exist — cannot build")
    else:
        fail(f"binary missing ({cfg.binary}) and no [build] cmake_source to build it")

    if cfg.ode_system:
        ok("ode_system declared: Hill-driver signs enabled") if P(cfg.ode_system).exists() \
            else warn(f"ode_system '{cfg.ode_system}' missing — Hill signs degrade to "
                      "direct-rate-law only")

    # --- param template ---------------------------------------------------
    tmpl = P(cfg.template)
    if not tmpl.exists():
        fail(f"template missing: {cfg.template}")
    else:
        try:
            n = len(ParamXMLRenderer(tmpl).template_defaults)
            ok(f"template parses: {n} parameters")
        except Exception as e:  # noqa: BLE001
            fail(f"template '{cfg.template}' failed to parse: {e}")

    # --- priors CSV -------------------------------------------------------
    pcsv = P(cfg.priors_csv)
    if not pcsv.exists():
        fail(f"priors_csv missing: {cfg.priors_csv}")
    else:
        try:
            df = pd.read_csv(pcsv)
            need = {"name", "median", "distribution", "dist_param2"}
            missing = need - set(df.columns)
            if missing:
                fail(f"priors_csv missing columns: {', '.join(sorted(missing))}")
            else:
                ok(f"priors_csv parses: {len(df)} params"
                   + ("" if "units" in df.columns else " (no 'units' column — sliders unit-less)"))
        except Exception as e:  # noqa: BLE001
            fail(f"priors_csv '{cfg.priors_csv}' failed to read: {e}")

    # --- model_structure.json --------------------------------------------
    msp = P(cfg.model_structure)
    if not msp.exists():
        fail(f"model_structure missing: {cfg.model_structure}")
    else:
        try:
            ms = json.loads(msp.read_text())
            nsp = len(ms.get("species", []))
            nrx = len(ms.get("reactions", []))
            if not nsp:
                fail("model_structure has no 'species'")
            elif not nrx:
                warn("model_structure has no 'reactions' — diagram will be node-only")
            else:
                ok(f"model_structure parses: {nsp} species, "
                   f"{len(ms.get('compartments', []))} compartments, {nrx} reactions")
        except Exception as e:  # noqa: BLE001
            fail(f"model_structure '{cfg.model_structure}' failed to parse: {e}")

    # --- optional simulator inputs ---------------------------------------
    for label, rel in (("drug_metadata", cfg.drug_metadata),
                       ("healthy_state", cfg.healthy_state)):
        if rel:
            ok(f"{label}: {rel}") if P(rel).exists() \
                else fail(f"{label} declared but missing: {rel}")

    # --- submodel priors --------------------------------------------------
    if cfg.submodel_priors:
        smp = P(cfg.submodel_priors)
        if not smp.exists():
            fail(f"submodel_priors declared but missing: {cfg.submodel_priors} "
                 "(the 'submodel' θ flavor will error)")
        else:
            try:
                med, _ = load_submodel_medians(smp)
                ok(f"submodel_priors parses: {len(med)} medians")
            except Exception as e:  # noqa: BLE001
                fail(f"submodel_priors '{cfg.submodel_priors}' failed to load: {e}")

    # --- scenarios + calibration targets ---------------------------------
    print()
    for s in cfg.scenarios:
        yp = P(s.yaml)
        if not yp.exists():
            fail(f"scenario '{s.id}': yaml missing ({s.yaml})")
            continue
        dirs, missing_dirs = [], []
        for d in s.target_dirs:
            (dirs if P(d).exists() else missing_dirs).append(P(d))
        ntgt = 0
        try:
            if dirs:
                tdf = load_calibration_targets(dirs)
                ntgt = len(tdf)
        except Exception as e:  # noqa: BLE001
            warn(f"scenario '{s.id}': target load raised {type(e).__name__}: "
                 f"{str(e)[:100]}")
        ok(f"scenario '{s.id}': yaml ok, {ntgt} targets"
           + (f" (t_end={s.t_end})" if s.t_end else ""))
        for md in missing_dirs:
            warn(f"scenario '{s.id}': target_dir missing "
                 f"({md.relative_to(repo)}) — skipped")

    # --- UI content -------------------------------------------------------
    for label, path in (("glossary", cfg.glossary_path),
                        ("interventions", cfg.interventions_path)):
        if path is None:
            continue
        if not path.exists():
            warn(f"{label} declared but missing: {path}")
            continue
        try:
            json.loads(path.read_text())
            ok(f"{label} parses: {path.name}")
        except Exception as e:  # noqa: BLE001
            fail(f"{label} '{path.name}' is not valid JSON: {e}")

    # --- verdict ----------------------------------------------------------
    print()
    if fails:
        print(f"\033[31m✗ {len(fails)} failure(s)\033[0m"
              + (f", {len(warns)} warning(s)" if warns else "")
              + " — fix the failures before serving.")
        return 1
    print(f"\033[32m✓ config is usable\033[0m"
          + (f" ({len(warns)} warning(s))" if warns else "")
          + f" — run: qsp-model-explorer --repo {repo}")
    return 0


def main():
    global DEFAULT_SRC, HOME, CFG, CACHE_ROOT, SCRATCH_ROOT
    ap = argparse.ArgumentParser(
        description="Interactive QSP model explorer. Point it at a model checkout "
                    "that ships an explorer.toml.")
    ap.add_argument("--repo", type=Path, default=Path.cwd(),
                    help="model checkout to serve (default: current directory)")
    ap.add_argument("--config", type=Path, default=None,
                    help="explorer.toml path (default: <repo>/explorer.toml or "
                         "<repo>/.model_explorer/explorer.toml)")
    ap.add_argument("--check", action="store_true",
                    help="validate the config + referenced files and exit "
                         "(no server, no build)")
    ap.add_argument("--port", type=int, default=8765)
    args = ap.parse_args()

    HOME = args.repo.expanduser().resolve()
    if not HOME.is_dir():
        raise SystemExit(f"--repo is not a directory: {HOME}")
    CFG = load_config(find_config(HOME, args.config))
    CACHE_ROOT = HOME / "cache/model_explorer"
    SCRATCH_ROOT = HOME.parent / f"{HOME.name}-explorer-refs"

    if args.check:
        raise SystemExit(run_check(HOME, CFG))

    print(f"[startup] model '{CFG.name}' at {HOME}")
    print(f"[startup] config: {CFG.config_path}")

    refresh_registry()
    # Default = the checkout being served (its own worktree).
    DEFAULT_SRC = HOME.name if HOME.name in SOURCES else next(iter(SOURCES))
    print(f"[startup] {len(SOURCES)} model sources: "
          + ", ".join(f"{s.id}[{s.label()}]" for s in SOURCES.values()))
    print(f"[startup] loading default source '{DEFAULT_SRC}' ...")
    src = SOURCES[DEFAULT_SRC]
    activate(src)
    for line in src.log:
        print(f"[{src.id}] {line}")
    if src.status != "ready":
        print(f"[startup] WARNING: default source not ready ({src.message})")

    srv = ThreadingHTTPServer(("127.0.0.1", args.port), Handler)
    print(f"\n  Model explorer running -> http://localhost:{args.port}/\n")
    try:
        srv.serve_forever()
    except KeyboardInterrupt:
        print("\nbye")


if __name__ == "__main__":
    main()
