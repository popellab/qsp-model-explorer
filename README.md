# qsp-model-explorer

An interactive, in-browser explorer for QSPIO / `qsp_sim` QSP models. Slide
parameters and watch calibration-target readouts, Hill-function saturation, the
reaction diagram, and full trajectories update in near-real-time (a single sim is
~0.4 s; results are cached). Point it at **any** model checkout that ships an
`explorer.toml` — including two branches at once, to see exactly what a commit did.

## What it needs from a model

The explorer drives a compiled `qsp_sim`-style simulator. A model is a git checkout
that provides:

- the compiled simulator **binary** (or a CMake project to build it);
- a **param template** (`param_all.xml`-style) — the parameter defaults;
- **`model_structure.json`** — species / compartments / reactions / rate laws;
- a **priors CSV** — parameter medians + distributions (drives slider ranges);
- optionally the generated **`ODE_system.cpp`** (for Hill-driver sign recovery and a
  staleness check), a **`submodel_priors.yaml`** (the SBI-faithful operating point),
  **scenario YAMLs**, and **calibration-target** directories.

All of that is declared in an `explorer.toml`, so nothing about one model's file
layout or species names is baked into the tool.

## Install

Requires Python ≥ 3.11 and [`qsp-hpc-tools`](https://github.com/) (the calibration
loader + param-XML renderer — not on PyPI, install it from its own repo first):

```bash
pip install -e /path/to/qsp-hpc-tools     # provides the qsp_hpc package
pip install -e /path/to/qsp-model-explorer
```

## Run

```bash
qsp-model-explorer --repo /path/to/your/model      # reads <repo>/explorer.toml
# or
python -m qsp_model_explorer.server --repo /path/to/your/model --port 8765
```

Then open http://localhost:8765/. Add other branches from the UI (they get a
detached worktree, compile, and load) and set one as the Δ-baseline to diff two
commits on the diagram.

## Validate a config before serving

When you first write an `explorer.toml` (or after moving files), dry-run it. `--check`
verifies that every referenced file exists and actually loads — it parses the param
template, the priors CSV, `model_structure.json`, each scenario's targets, the submodel
priors, and the UI content — then exits **without booting the server or building the
binary**:

```bash
qsp-model-explorer --repo /path/to/your/model --check
```

Missing required files or unparseable schemas are `[FAIL]` (exit code 1); optional
gaps (a missing target dir, no `ode_system`, an absent binary that will be compiled on
launch) are `[warn]` (exit 0). Sample output:

```
checking model 'QSPIO PDAC'

  [ok]   binary present: cpp/sim/build/qsp_sim
  [ok]   template parses: 719 parameters
  [ok]   priors_csv parses: 264 params
  [ok]   model_structure parses: 177 species, 11 compartments, 330 reactions
  [ok]   submodel_priors parses: 172 medians
  [ok]   scenario 'baseline_no_treatment': yaml ok, 27 targets (t_end=1.0)
  [warn] scenario '...': target_dir missing (...) — skipped

✓ config is usable (1 warning(s)) — run: qsp-model-explorer --repo ...
```

Because it exits non-zero on failure, it also works as a CI / pre-commit guard that a
model repo's `explorer.toml` still points at real files.

## `explorer.toml`

Put it at the model repo root, or under `.model_explorer/explorer.toml`, or pass
`--config PATH`. Paths are **repo-relative** (so they hold across branches), except
`[views].glossary` / `interventions`, which are relative to the config file.

```toml
[model]
name = "My QSP model"

[build]
binary = "cpp/sim/build/qsp_sim"     # required: the simulator binary
build_dir = "cpp/sim/build"
cmake_source = "cpp/sim"             # omit for a prebuilt binary (never compiled)
target = "qsp_sim"                   # cmake --build target (optional)
ode_system = "cpp/qsp/ode/ODE_system.cpp"   # optional: Hill signs + staleness check
seed_deps = true                     # seed SUNDIALS/yaml-cpp _deps/*-src from a donor checkout

[paths]
template = "resources/cpp/param_all.xml"     # required
priors_csv = "parameters/pdac_priors.csv"    # required
model_structure = "model_structure.json"     # required
drug_metadata = "resources/cpp/drug_metadata.yaml"     # optional → --drug-metadata
healthy_state = "resources/cpp/healthy_state.yaml"     # optional → --evolve-to-diagnosis
submodel_priors = "notes/calibration/submodel_priors.yaml"   # optional (enables the "submodel" θ flavor)

[sim]
extra_args = ["--min-cadence-hours", "12"]   # appended to every qsp_sim call

[[scenarios]]                        # one or more
id = "baseline_no_treatment"
yaml = "scenarios/baseline_no_treatment.yaml"
t_end = 1.0
target_dirs = ["calibration_targets/baseline_no_treatment"]

[views]                              # all optional — omit for auto-generated defaults
curated = ["k_C1_growth", "k_C1_death"]   # default slider set (else first 12 params)
traj_species = ["V_T.C1", "V_T.CD8"]      # trajectory panel (else first 12 species)
glossary = "glossary.json"                # UI content, relative to this file
interventions = "interventions.json"

[[views.theta_flavors]]              # parameter-vector overlays (else template/csv/submodel)
id = "submodel"
label = "submodel medians"
desc = "the SBI-faithful prior center"

[views.modules]                      # named "lenses" over subsets of species
"Overview" = ["V_T.C1", "V_T.CD8", "V_T.Treg"]
```

Unknown species/param names in `[views]` are filtered against the live model at load,
so stale entries drop silently. A model that omits every `[views]` key still loads
with auto-generated defaults.

## The simulator contract

The tool invokes the binary as:

```
<binary> --param P.xml --scenario S.yaml --csv-out out.csv --rules-out rules.txt \
         --t-end-days T [--drug-metadata M.yaml] [--evolve-to-diagnosis H.yaml] <extra_args>
```

`out.csv` must have a `Time` column plus one column per observable (species amounts,
compartment volumes, aggregates, and `H_*` Hill rules); `rules.txt` carries the
evaluated rules. This is the QSPIO `qsp_sim` CLI — a genuinely different simulator
would need a backend adapter (not yet abstracted).
