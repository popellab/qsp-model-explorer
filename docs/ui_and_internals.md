# QSP Model Explorer

An interactive, browser-based tool for poking at the PDAC QSP model: drag a
parameter and watch the perturbation cascade through the species diagram, the
calibration-target readouts, and the Hill-function saturation strip — all backed
by live `qsp_sim` runs (~0.4 s each, cached so revisits are instant).

Built for two things: **debugging** (spot chokepoints and oversaturated Hill
gates, see how far each target is off) and **explaining the model** to
collaborators/reviewers (a living version of the cascade cartoon; share a link to
any counterfactual state).

## Run it

```bash
.venv/bin/python tools/model_explorer/server.py          # defaults to :8765
```

Then open <http://localhost:8765/>. It's a local single-user tool — no external
services, no CDN, works offline.

## What you see

- **Left — parameters.** A curated dozen sliders to start; type in the box to add
  any of the 262. Cyan names are submodel-constrained (their range is
  `exp(mu ± 2σ)` from the submodel posterior); grey names use the CSV prior.
- **Center — the model.** Every species as a node, auto-generated from
  `model_structure.json` (nothing hand-drawn, so it never drifts). Nodes color by
  **Δ from baseline** (blue = down, red = up, log₂) so a perturbation visibly
  spreads; switch to **abs (log)** for raw magnitudes. Scroll to zoom, drag to pan,
  **fit** to reset. Hover a node for its value, units, and Δ.
  - **Module lens** (the view selector, top-left of the diagram). The default is
    **Overview** — a clean whiteboard-level cancer-immunity cascade (cancer → antigen →
    DC → CD8 → killing, with Treg / TAM / MDSC / CAF / collagen / cytokine regulators) —
    so you start on the cartoon and drill down into a specific module or the full map
    only when you want to (progressive disclosure). "All compartments"
    is the full map: a **force-directed organic layout** where species cluster by
    compartment into soft blob outlines and reaction edges thread within and across
    them — connected species pull together, so the trafficking compartments end up
    adjacent where they exchange cells (positions are seeded deterministically, so the
    layout is stable across redraws). The **edge-basis** dropdown (all-view only)
    chooses what an edge means: **regulation** (default) draws influence edges — every
    rate-law regulator + reactant → product, the real mechanistic wiring; **conversions**
    draws only stoichiometric reactant → product mass transfers, a much sparser subset
    (most QSP reactions are synthesis/secretion/degradation with one side empty, or
    drive their product through the rate law, so they carry no reactant→product edge).
    Regulation includes **Hill-mediated** drivers (cytokines acting through `H_*`
    columns), recovered from the generated C++ Hill formulas — so the cytokine→target
    wiring (IL-10, TGF-β, IL-12, IFN-γ, IL-6, CXCL9 …) is drawn, not hidden. The six presets
    are the
    mechanistic cascades from `scripts/cascade_figures.m` (Input & APC, T-cell
    activation, Tumor killing, Checkpoint / suppression, CAF subtypes, Chemokine
    trafficking): each shows only that module's species, laid out **left→right in
    signal-flow order** (topological layering over the module's *influence* edges —
    a regulator/reactant → product, since the real signal flow lives in the rate
    laws, not the sparse stoichiometry; feedback edges are drawn but don't drive the
    layering). "＋ new custom module…" opens a builder — name it, add any species,
    save (persists in the browser via localStorage). Synapse-complex species
    (the `syn_CD8_C1` / `syn_CD8_APC` / `syn_M_C` compartments — dense
    receptor-occupancy books, 30–80 states each) are **hidden by default** from both
    the preset modules and the all-compartments map; the **syn** toggle next to the
    view selector brings them (and their compartment boxes) back.
  - **Plain labels** (toggle, on by default — *review mode*). Module-view nodes and
    every tooltip show human-readable biology names (`CD8 T cell`, `M2 macrophage
    (TAM)`, `TGF-β`) from `glossary.json` instead of model IDs (`V_T.CD8`) — for
    reviewing the mechanism with biologists/clinicians who don't speak modeler. Flip
    it off for model IDs; the tooltip always shows both. Applies to the
    all-compartments view too — its blobs are pushed apart so no two compartments
    overlap, and each compartment label sits above its blob with a halo so species
    nodes never cover it. Anything unlisted falls back to its model ID.
  - **Activation vs inhibition** (module views). Each influence edge carries an
    SBGN-style head — a pointed **arrow** where the source *promotes* the product,
    a flat **tee bar** where it *inhibits* — derived server-side from the local sign
    of ∂flux/∂source at the baseline operating point (`graph.signs`). The edge
    tooltip reads it as a sentence ("M2 macrophage inhibits M1 macrophage"). This
    covers **both** direct rate-law regulators **and Hill-mediated ones** — a cytokine
    that acts only through a Hill column (`H_IL10`, `H_TGFb`) is recovered by reading
    the Hill's formula out of the generated C++ (`hill_deps.py`) and composing
    d(flux)/d(Hill) with d(Hill)/d(species), so e.g. IL-10 ⊣ DC maturation, TGF-β →
    Treg, IFN-γ → CXCL9 all draw as signed edges.
  - **Flux edges** (toggle; on by default in a module view). Directed reactant/
    regulator → product edges weighted by reaction throughput: **thickness ∝ |flux|**,
    color by **Δ flux from baseline** (or **mag** = absolute |flux|). This is how you
    spot chokepoints — an edge whose flux swings hard when you nudge an upstream
    parameter. Hover an edge for the reaction, its rate law, and the flux.
    - **⟿ flow** (toggle, on): animates each edge as a **traveling sine wave** rolling
      in the flux direction. All waves travel at the **same speed**; **|flux| sets the
      wave amplitude** — a high-flux reaction ripples with tall waves, a throttled one
      stays nearly flat (amplitude tracks throughput, not velocity). The animation
      **freezes at the endpoint** (the default resting view) so a static snapshot stays
      calm — scrub back into the trajectory to see it flow. Turn the toggle off to
      freeze it everywhere.
    - **⊘ gates** (toggle, on in module views): Hill-saturation "valves" on the
      reaction edges. Each reaction that carries Hill terms gets a **row of valves at
      its midpoint — one per `H_*`** — each a ring with an inner disc sized + colored
      (green → amber → **red** ≥ 0.95) by that Hill's saturation at the current
      timepoint. A red valve on a low-flux edge is a chokepoint: that gate is closed
      and throttling the flow. Each valve has a generous hover target — mouse over it
      for the gate name and value.
- **Perturb strip** (above the diagram). One-click **knockout / blockade** presets
  for mechanistic review — **Block TGF-β / IL-10 / IL-6 / VEGF / collagen / IFN-γ /
  CXCL9 / CXCL12 / CCL2**, **Deplete Treg / TAM / MDSC**, **Block CD8 expansion**.
  Each is honest: just a set of secretion/recruitment/proliferation rate constants
  set to 0, applied **on top of** the sliders and re-simulated. They reshape the
  evolve-to-diagnosis run too, so the effect shows even in the diagnosis snapshot.
  This is the review verb a wet-lab/clinical expert reasons in ("block TGF-β — does
  exhaustion fall? does collagen drop?") — flip to **baseline (your Δ)** coloring to
  read the answer straight off the diagram. Toggle several to stack them; **clear**
  resets. Each slider also carries a **⊘** to knock that one parameter to 0.
- **Bottom — time slider.** A single `qsp_sim` run carries the whole trajectory,
  so scrubbing needs **no new simulation** — it just re-indexes the frame you're
  looking at. Drag to move nodes/edges/Hills to any timepoint; **▶** plays the
  trajectory, **⇥** jumps back to the endpoint. **Δ vs** picks the node/edge
  coloring reference: **t₀ (dynamics)** — the default — colors vs the diagnosis
  start, so scrubbing shows the system *evolving* (the T-cell cascade lighting up
  red as priming ramps under treatment); **baseline (your Δ)** colors vs the
  no-override baseline at the same time, isolating what *your* parameter change did
  (only meaningful once you've moved a slider). Calibration targets don't move with
  the slider — each is evaluated at its own fixed measurement time server-side.
- **Right — readouts.** Calibration targets in stable alphabetical order (so a
  slider drag doesn't reshuffle the rows) — flip the **sort by |z|** toggle to
  rank by worst-fit instead. Each has a CI bar (band + median tick +
  current-value marker) colored green (inside CI) → amber → red. Below that,
  **Key trajectories** — small-multiple sparklines of the readouts an expert
  critiques by *shape* (tumor burden, CD8 / exhausted CD8 / Treg / TAM / MDSC, IFN-γ /
  TGF-β / IL-10, CAF / collagen). Each shows the current run (cyan) over the faint
  no-override baseline (grey) with an amber marker at the scrubbed time and the
  current value + ×(vs baseline); drawn from the trajectory already in hand, so it's
  free. This is where "block TGF-β and the Treg curve drops below baseline" reads at
  a glance. Below, the
  Hill-saturation strip: every `H_*` at the scrubbed timepoint,
  red when ≥ 0.95 (a dead/oversaturated gate). The header shows how many targets
  sit outside their CI.

## Scenarios

Switch scenario (top) between the diagnosis snapshot and the neoadjuvant GVAX /
GVAX+nivo treatment runs. Each scenario carries its own target set and endpoint
time; the baseline reference for Δ-coloring is re-established per scenario.

## Shareable counterfactual links

Four URL params compose into a shareable link:

- `?scenario=<name>` — pick the scenario (diagnosis snapshot / GVAX / GVAX+nivo).
- `?view=<module name>` — open a specific module lens (`all`, a preset, or `custom:<name>`).
- `?ovr=name:value,name:value` — preset slider overrides.
- `?int=name,name` — preset active interventions (e.g. `block_tgfb,deplete_treg`).
- `?t=<0..1>` — park the time slider at a fraction of the trajectory (0 = t₀, 1 = endpoint).

e.g. `http://localhost:8765/?scenario=gvax_nivo_neoadjuvant_zheng2022&view=Tumor%20killing&ovr=k_C_CD8_exh:0.18&t=0.2`
points a collaborator at a specific "what if," in a specific cascade, at a specific
day — without walking them through the controls.

## How it works

`server.py` (stdlib HTTP, no framework) builds an SBI-faithful baseline by
overlaying `notes/calibration/submodel_priors.yaml` medians onto
`resources/cpp/param_all.xml` (via the maintained `scripts/submodel_medians`
selector — distribution-correct, not a naive `exp(mu)`). Each `/simulate`
renders the param XML, runs `cpp/sim/build/qsp_sim` for the chosen scenario,
parses the trajectory, and evaluates every calibration target's
`compute_test_statistic` (reused from `qsp_hpc.calibration.yaml_loader`). Results
cache under `cache/model_explorer/` keyed by (scenario, overrides).

Design notes and the scoping rationale: `notes/tools/interactive_model_explorer_scoping.md`.

## Caveats / next steps

- Module layout uses the curated cascade order to break feedback cycles; a species
  that legitimately sits both up- and downstream lands by its curated position.
- Params that change the diagnosis burn-in shift the whole baseline snapshot —
  that's the intended counterfactual, not a bug.
- Very large growth-rate bumps can trip the evolve-to-diagnosis sanity floor;
  those return a sim error in the status bar rather than a result.
