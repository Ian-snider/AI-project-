# Refinery Pipe Network: Flow Simulator

Fluid mechanics AI project, group of 4: Ian Snider (lead), Drew Corbitt,
Perry Meredith, Zac Neikamp. Working build, iteration 3 (November 24, 2026),
the finalized handoff baseline for the next AI improvement pass.
Dev labels are visible in this build (SHOW_DEV_NOTES = true in index.html).

## Run it
Open index.html in a browser, or serve the folder:
    python3 -m http.server 8901
then visit http://localhost:8901

## Verify the numbers (run before showing anyone)
    python3 tests/run_tests.py     # Python reference core: known-answer cases
    node tests/crosscheck_js.js    # the page solver must reproduce the core

Both must print ALL PASS. The laminar case is checked against the
Hagen-Poiseuille result from class; the preset network must balance flow
at every junction; a 10-pipe, 4-pump stress network must converge; a
deliberately reversed branch must be flagged.

## Layout
- index.html — the site (GitHub Pages ready, computes in the browser).
  Four tabs: Application (builder and demo), Theory and Math,
  Code and Decisions, Video and Slides.
- solver.js — the page solver, a mirror of the Python core
- examples.js — the use-case library: every network is plain data
  (nodes and pipes tables). Six examples: refinery transfer, parallel
  pipes, gravity feed, the starved branch, laminar viscous oil, and
  Example 6, the final demo case: two source tanks fill one open
  receiving tank over the top, to its max fill level.
- core/network.py — Python reference implementation (verification standard)
- tests/ — the checker above, plus test_vectors.json it generates
- checks/ — screenshots from the visual check (desktop and phone)

## The builder
Edit any row of the node and pipe tables, or add and remove rows.
Limits enforced by guardrails: 12 nodes, 10 pipes, 4 pumps; every
connected part of the network needs at least one fixed-head node;
a pipe cannot run from a node to itself; diameters and lengths must
be positive. Reversed flow is flagged, never hidden.

## The tank-fill panel
Under the results on the Application tab. It works for any sink node:
pick the receiving tank, give its inside diameter, bottom elevation,
start level and max fill level (levels are heights above the tank
bottom). Three solve-for modes: fill time from the solved network
flow; the flow a target fill time requires, checked honestly against
what the network and pump can deliver; and the diameter that carries
the required flow at a target velocity. The discharge is modelled
over the top of an open tank: the sink node's head is the discharge
elevation at the top of the tank, so the fill flow is constant and
the liquid level is volume bookkeeping only.

## Next iterations (assigned to the next AI improvement pass)
3D view of the same model with elevation, economics module (total
cost summed over the network), guided teaching mode, P&ID page, the
presentation video and slides on Tab 4, Drew's second-model iteration
for the comparison, then the final build with dev labels stripped.
