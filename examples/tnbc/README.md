# TNBC example

A self-contained model to try the explorer on, built from the published, MIT-licensed
[QSPIO-TNBC](https://github.com/popellab/QSPIO-TNBC) model (Wang et al., Popel lab) — a
triple-negative breast cancer immuno-oncology QSP model with 154 species and 234
reactions.

From a clone of this repo:

```bash
pip install -e /path/to/qsp-hpc-tools     # dependency
pip install -e .                          # the explorer
qsp-model-explorer --repo examples/tnbc   # → http://localhost:8765/
```

That's the whole setup — the compiled simulator is bundled, so there's no build step.

## What's here

| File | Origin |
|------|--------|
| `cpp/sim/build/qsp_sim` | compiled simulator, **macOS arm64**, from QSPIO-TNBC |
| `resources/cpp/param_all.xml` | parameter template, from QSPIO-TNBC |
| `cpp/qsp/ode/ODE_system.cpp` | generated ODE source (drives the signed diagram edges), from QSPIO-TNBC |
| `model_structure.json` | species/reactions, generated with `qsp-export-model` |
| `priors.csv` | 88 rate constants pulled from the template (slider ranges) |
| `scenarios/baseline.yaml` | a no-treatment, 400-day natural-history run |
| `resources/cpp/drug_metadata.yaml` | empty (no dosing in this scenario) |
| `explorer.toml` | ties it together |

See `NOTICE` for attribution and license.

## Not on macOS?

The bundled binary is macOS arm64. On Linux/Windows, clone
[QSPIO-TNBC](https://github.com/popellab/QSPIO-TNBC), build `qsp_sim` there, and copy it
to `cpp/sim/build/qsp_sim` — everything else in this directory is platform-independent.
