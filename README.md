# qsp-model-explorer

A browser UI for poking at QSPIO / `qsp_sim` models. Move a parameter slider and the
reaction diagram, calibration-target scores, Hill-gate saturation, and time courses
re-render in about half a second, so you can see what a parameter does without a
rebuild-and-replot loop. It can also load two branches at once and show the difference
between them.

It works with any model checkout that includes an `explorer.toml` (below), so it isn't
tied to one model.

## What's in the UI

- Sliders for the model parameters, with the reaction diagram redrawn on each change and
  edges scaled by flux.
- Calibration targets shown as a z-score against the measured value, colored by whether
  the value falls inside the CI.
- Hill functions drawn as valves on the edges they gate. A saturated valve on a low-flux
  edge marks a chokepoint.
- Sparklines for the key species over the run, with a time slider.
- A second "baseline" branch to diff against: run branch B's equations at branch A's
  parameters and the diagram and target scores show the delta.
- Preset interventions (block a cytokine, deplete a cell type) as one-click overrides.

## Running it

You need Python 3.11+ and `qsp-hpc-tools` installed.

```bash
pip install -e /path/to/qsp-hpc-tools
pip install -e /path/to/qsp-model-explorer

qsp-model-explorer --repo /path/to/your/model     # serves at http://localhost:8765/
```

If the model already has an `explorer.toml`, there's nothing else to set up. You can add
more branches from the UI and pick one as the Δ-baseline to compare two commits.

## Pointing it at your own model

The `explorer.toml` tells the explorer where the model's files are. A minimal one needs a
binary, three data files, and a scenario:

```toml
[build]
binary = "build/qsp_sim"                     # or set cmake_source to build it

[paths]
template = "param_all.xml"                   # parameter defaults (XML)
priors_csv = "priors.csv"                    # medians + distributions, used for slider ranges
model_structure = "model_structure.json"     # species / reactions / rate laws

[[scenarios]]
id = "baseline"
yaml = "scenarios/baseline.yaml"
```

Put it at the repo root or in `.model_explorer/explorer.toml`. Everything past those four
keys is optional, and nothing about the model's layout or species names is baked into the
tool.

Before starting the server, run `--check` to confirm the config points at files that
exist and parse. It validates the config without booting the server or building anything:

```bash
qsp-model-explorer --repo /path/to/your/model --check
```
```
  [ok]   template parses: 719 parameters
  [ok]   priors_csv parses: 264 params
  [ok]   model_structure parses: 177 species, 11 compartments, 330 reactions
  [warn] scenario '...': target_dir missing (...) — skipped
✓ config is usable (1 warning(s))
```

Missing or unparseable required files are `[FAIL]` and exit 1; optional gaps are `[warn]`
and exit 0, so `--check` also works as a CI or pre-commit check.

<details>
<summary>Full <code>explorer.toml</code> reference</summary>

<br>

Put it at the model repo root, in `.model_explorer/explorer.toml`, or pass `--config
PATH`. Paths are repo-relative (so they hold across branches), except `[views].glossary`
and `interventions`, which are relative to the config file.

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
drug_metadata = "resources/cpp/drug_metadata.yaml"     # optional -> --drug-metadata
healthy_state = "resources/cpp/healthy_state.yaml"     # optional -> --evolve-to-diagnosis
submodel_priors = "notes/calibration/submodel_priors.yaml"   # optional -> "submodel" theta flavor

[sim]
extra_args = ["--min-cadence-hours", "12"]   # appended to every qsp_sim call

[[scenarios]]                        # one or more
id = "baseline_no_treatment"
yaml = "scenarios/baseline_no_treatment.yaml"
t_end = 1.0
target_dirs = ["calibration_targets/baseline_no_treatment"]

[views]                              # all optional; omitted keys fall back to defaults
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

Unknown species or param names in `[views]` are filtered against the live model at load,
so stale entries drop silently, and a model that omits every `[views]` key still loads
with generated defaults.

</details>

<details>
<summary>What the explorer expects from a model, and the simulator CLI</summary>

<br>

The explorer drives a compiled `qsp_sim`-style simulator. A model checkout provides:

- the compiled simulator binary (or a CMake project to build it);
- a param template (`param_all.xml`-style) with the parameter defaults;
- `model_structure.json` with species, compartments, reactions, and rate laws;
- a priors CSV of parameter medians and distributions (drives the slider ranges);
- optionally the generated `ODE_system.cpp` (Hill-driver sign recovery and a staleness
  check), a `submodel_priors.yaml` (the SBI-faithful operating point), scenario YAMLs, and
  calibration-target directories.

The tool invokes the binary as:

```
<binary> --param P.xml --scenario S.yaml --csv-out out.csv --rules-out rules.txt \
         --t-end-days T [--drug-metadata M.yaml] [--evolve-to-diagnosis H.yaml] <extra_args>
```

`out.csv` needs a `Time` column plus one column per observable (species amounts,
compartment volumes, aggregates, and `H_*` Hill rules); `rules.txt` carries the evaluated
rules. This is the QSPIO `qsp_sim` CLI. A different simulator would need a backend adapter,
which isn't abstracted yet.

</details>
