"""Hill-mediated regulation: recover which species drive each Hill function, and
the sign of that drive, from the generated C++ ODE right-hand side.

`model_structure.json` carries reaction rate laws but *not* the Hill (`H_*`) rule
formulas — a reaction's rate law references `H_IL10_APC` as an opaque column, so
the regulator (IL-10) that actually drives it is invisible to a bare-token scan.
The generated `cpp/qsp/ode/ODE_system.cpp` `f()` body, however, spells every Hill
out as `realtype AUX_VAR_H_X = <expr over SPVAR(...)/PARAM(...)>;`. We parse those,
translate them to evaluable Python, and expose the local sign of d(Hill)/d(species)
at a given operating point. The explorer composes that with d(flux)/d(Hill) — read
straight from the rate law, where the Hill *is* a column — to get the full
species → reaction regulatory edge and its sign, with no C++ flux→reaction mapping.

Signs, not magnitudes: unit-conversion factors baked into the C++ params are
positive constants, so they cannot flip the sign of a monotonic derivative — which
is all a Hill's driver relationship is. We only report the sign.
"""
from __future__ import annotations

import math
import re
from pathlib import Path

_FUNCS = {"pow": pow, "exp": math.exp, "log": math.log, "log10": math.log10,
          "sqrt": math.sqrt, "abs": abs, "max": max, "min": min, "__builtins__": {}}
_STD = (("std::pow", "pow"), ("std::exp", "exp"), ("std::log10", "log10"),
        ("std::log", "log"), ("std::sqrt", "sqrt"), ("std::fabs", "abs"),
        ("std::abs", "abs"), ("std::max", "max"), ("std::min", "min"))
_STMT = re.compile(r"^\s*realtype (AUX_VAR_\w+)\s*=\s*(.*);\s*$")


class HillDeps:
    """Parses ODE_system.cpp `f()` for the `AUX_VAR_*` assignments and exposes the
    Hill formulas + their direct driver species. `ok` is False if parsing found no
    Hills (missing/renamed file, changed codegen) — callers should degrade to
    direct-token signs only."""

    def __init__(self, ode_path: Path, species_names, param_names):
        self.ok = False
        self.aux: list[tuple[str, object]] = []      # (name, compiled) in eval order
        self.hill_code: dict[str, object] = {}       # H_X -> compiled expr
        self.hill_drivers: dict[str, set] = {}       # H_X -> {species names}
        self._par = set(param_names)
        self._sp_enum = {"SP_" + n.replace(".", "_"): n for n in species_names}
        try:
            self._parse(Path(ode_path))
        except Exception as e:  # noqa: BLE001
            print(f"[hill_deps] parse failed ({e}); Hill-mediated edges disabled")

    # -- parsing -----------------------------------------------------------
    def _translate(self, expr: str) -> str:
        def sp(m):
            e = m.group(1)
            return f"SP[{self._sp_enum[e]!r}]" if e in self._sp_enum else "0.0"
        expr = re.sub(r"SPVAR\((SP_\w+)\)", sp, expr)
        expr = re.sub(r"_species_var\[(SP_\w+)\]", sp, expr)
        expr = re.sub(r"PARAM\((P_\w+)\)", lambda m: f"PR.get({m.group(1)[2:]!r},1.0)", expr)
        expr = re.sub(r"_class_parameter\[(P_\w+)\]", lambda m: f"PR.get({m.group(1)[2:]!r},1.0)", expr)
        for a, b in _STD:
            expr = expr.replace(a, b)
        return expr

    def _parse(self, ode_path: Path):
        lines = ode_path.read_text().splitlines()
        start = next(i for i, l in enumerate(lines) if "ODE_system::f(" in l)
        end = next(i for i, l in enumerate(lines) if i > start and "ODE_system::jac(" in l)
        for l in lines[start:end]:
            m = _STMT.match(l)
            if not m:
                continue
            name, cpp = m.group(1), m.group(2)
            py = self._translate(cpp)
            self.aux.append((name, compile(py, "<aux>", "eval")))
            if name.startswith("AUX_VAR_H_"):
                hname = name[len("AUX_VAR_"):]
                self.hill_code[hname] = compile(py, "<hill>", "eval")
                self.hill_drivers[hname] = set(re.findall(r"SP\['([^']+)'\]", py))
        self.ok = bool(self.hill_code)

    # -- evaluation --------------------------------------------------------
    def _base_aux(self, spval: dict, params: dict) -> dict:
        ns = dict(_FUNCS); ns["SP"] = spval; ns["PR"] = params
        for name, code in self.aux:
            try:
                ns[name] = eval(code, ns)  # noqa: S307
            except Exception:  # noqa: BLE001
                ns[name] = float("nan")
        return ns

    def driver_signs(self, spval: dict, params: dict) -> dict:
        """{H_X: {species: +1/-1}} — local sign of d(Hill)/d(species) at this
        operating point. Only the Hill's own formula is recomputed under the
        perturbation (V_T / totals held at base), so this is the direct
        regulator effect, not volume dilution."""
        if not self.ok:
            return {}
        base = self._base_aux(spval, params)
        out: dict = {}
        for hname, code in self.hill_code.items():
            bv = base.get("AUX_VAR_" + hname)
            if bv is None or bv != bv:
                continue
            dd = {}
            for s in self.hill_drivers.get(hname, ()):
                v0 = spval.get(s)
                if v0 is None:
                    continue
                sp2 = dict(spval); sp2[s] = v0 + (abs(v0) * 0.02 if v0 != 0 else 1e-9)
                ns = dict(base); ns["SP"] = sp2; ns["PR"] = params
                try:
                    nv = eval(code, ns)  # noqa: S307
                except Exception:  # noqa: BLE001
                    continue
                if nv != nv:
                    continue
                d = nv - bv
                thr = max(abs(bv) * 1e-6, 1e-30)
                if d > thr:
                    dd[s] = 1
                elif d < -thr:
                    dd[s] = -1
            if dd:
                out[hname] = dd
        return out
