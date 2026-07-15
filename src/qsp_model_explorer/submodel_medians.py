#!/usr/bin/env python3
"""Compose submodel posterior medians — the single source of truth for the
"submodel shadow" layer that sits between the base CSV priors and any joint
NPE posterior.

`submodel_priors.yaml` (produced by `scripts/regen_submodel_priors.py`) holds
a posterior marginal per QSP parameter that has submodel-target data. Several
consumers need the *point* value (median) of those marginals overlaid onto the
base priors:

  - the C++ trace path (`scripts/trace_at_submodel_medians.py`, which overlays
    onto `param_all.xml`), and
  - the MATLAB median simulation (`scripts/run_median_simulation.m`, used by the
    cascade figures), which overlays onto the SimBiology model.

Both call `load_submodel_medians()` here so the value-selection rules (median,
with an `exp(mu)` fallback when a tiny-mu lognormal underflows its stored
median, and an invgamma approximation) live in exactly one place.

As a CLI, this writes `{"medians": {name: value, ...}}` JSON — the same schema
`run_median_simulation.m` already consumes for joint posterior medians:

    python scripts/submodel_medians.py --out-json /tmp/submodel_medians.json
"""
from __future__ import annotations

import argparse
import json
import sys
from math import exp
from pathlib import Path

import yaml

REPO = Path(__file__).resolve().parent.parent
DEFAULT_SUBMODEL = REPO / "notes" / "calibration" / "submodel_priors.yaml"


def load_submodel_medians(submodel_yaml: Path) -> tuple[dict[str, float], dict]:
    """Return ({param_name: median}, meta) from a submodel_priors.yaml.

    Value selection per parameter marginal:
      - use ``marginal.median`` when non-zero;
      - for a lognormal whose stored median underflowed to 0.0, fall back to
        ``exp(mu)`` (see feedback_submodel_priors_median_underflow);
      - for an invgamma, approximate the median as ``scale / (shape + 1)``;
      - otherwise skip the parameter (leave the base prior in place).

    ``meta`` carries ``underflow_fallbacks`` and ``skipped`` name lists for
    diagnostics by callers.
    """
    sm = yaml.safe_load(Path(submodel_yaml).read_text())
    medians: dict[str, float] = {}
    underflow_fallbacks: list[str] = []
    skipped: list[str] = []
    for entry in sm.get("parameters", []):
        name = entry["name"]
        marg = entry.get("marginal", {}) or {}
        # Marginal fields are sometimes serialized as strings (e.g. '1e-06');
        # coerce to float so downstream consumers get real numbers.
        med = _as_float(marg.get("median", 0.0))
        if med == 0.0 and marg.get("distribution") == "lognormal" and "mu" in marg:
            med = exp(_as_float(marg["mu"]))
            underflow_fallbacks.append(name)
        if med == 0.0:
            if marg.get("distribution") == "invgamma" and "shape" in marg and "scale" in marg:
                med = _as_float(marg["scale"]) / (_as_float(marg["shape"]) + 1.0)
            else:
                skipped.append(name)
                continue
        medians[name] = med
    return medians, {"underflow_fallbacks": underflow_fallbacks, "skipped": skipped}


def _as_float(x) -> float:
    """Coerce a marginal field to float; non-numeric / missing -> 0.0."""
    try:
        return float(x)
    except (TypeError, ValueError):
        return 0.0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--submodel-yaml", type=Path, default=DEFAULT_SUBMODEL,
                    help="submodel_priors.yaml to read (default: notes/calibration/...)")
    ap.add_argument("--out-json", type=Path, required=True,
                    help="output JSON path, written as {\"medians\": {name: value}}")
    args = ap.parse_args()

    if not args.submodel_yaml.is_file():
        print(f"submodel yaml not found: {args.submodel_yaml}", file=sys.stderr)
        return 1

    medians, meta = load_submodel_medians(args.submodel_yaml)
    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.write_text(json.dumps({"medians": medians}, indent=2))
    print(f"wrote {len(medians)} submodel medians to {args.out_json}"
          f" ({len(meta['underflow_fallbacks'])} exp(mu) fallbacks,"
          f" {len(meta['skipped'])} skipped)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
