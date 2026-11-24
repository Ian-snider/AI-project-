// Use-case library for the generic pipe-network builder. Iteration 2.
// Every network is data: a list of nodes and a list of pipes. The solver
// (solver.js) does not know which example is loaded; that is the point.
// Symbol key (SI units): H total head at a node (m), fixed for sources and
// sinks, solved for junctions; Q pipe flow (m3/s); D diameter (m);
// L length (m); eps roughness (m); K fitting loss sum (dimensionless);
// s pump speed ratio (dimensionless); x, y diagram position (per cent of
// the diagram width and height, display only, not physics).
"use strict";
const STEEL_EPS = 4.5e-5; // commercial steel [R5]

const FLUIDS = {
  "Crude oil (preset)": { rho: 850, mu: 0.008 },
  "Water, 20 C (check case)": { rho: 998, mu: 0.001002 },
  "Viscous oil (preset)": { rho: 900, mu: 0.35 }
};

// Stated illustrative pump presets [R7], not manufacturer data.
const PUMPS = {
  "Pump 1": { H0: 40, kp: 900, eta: 0.75 },
  "Pump 2": { H0: 55, kp: 1400, eta: 0.75 },
  "Pump 3": { H0: 30, kp: 600, eta: 0.75 },
  "Pump W (weak)": { H0: 28, kp: 700, eta: 0.75 },
  "Pump L (loading)": { H0: 15, kp: 400, eta: 0.75 }
};

function N(name, kind, head, x, y) { return { name, kind, head, x, y }; }
function P(name, a, b, D, L, K, pumpName, s) {
  return { name, a, b, D, L, eps: STEEL_EPS, K, pumpName: pumpName || null, s: s || 1 };
}

const EXAMPLES = {
  "refinery": {
    title: "Example 1 — Refinery transfer network",
    fluid: "Crude oil (preset)",
    blurb: "Two crude tanks at different levels feed one header, a pump pushes the combined flow, and the stream splits to two process units at different heights. The original project preset.",
    lesson: "Watch Pipe 2: the low tank contributes far less than the high tank, and the pump, not gravity, sets the split between the two units.",
    nodes: [
      N("Tank A", "source", 25, 8, 16), N("Tank B", "source", 15, 8, 48),
      N("Junction J", "junction", null, 38, 32), N("Split S", "junction", null, 64, 32),
      N("Unit 1", "sink", 18, 92, 16), N("Unit 2", "sink", 8, 92, 48)
    ],
    pipes: [
      P("Pipe 1", "Tank A", "Junction J", 0.20, 120, 1.25),
      P("Pipe 2", "Tank B", "Junction J", 0.20, 90, 0.95),
      P("Pipe 3", "Junction J", "Split S", 0.25, 200, 2.95, "Pump 1"),
      P("Pipe 4", "Split S", "Unit 1", 0.15, 150, 2.55),
      P("Pipe 5", "Split S", "Unit 2", 0.15, 110, 2.25)
    ]
  },
  "parallel": {
    title: "Example 2 — Two pipes in parallel",
    fluid: "Water, 20 C (check case)",
    blurb: "One high tank, one low tank, and two pipes side by side connecting them. No pump, no junction: the simplest network there is, and the split still surprises people.",
    lesson: "The wide pipe is 1.5 times the diameter, so it has 2.25 times the area, but it carries about 2.85 times the flow: a bigger pipe also loses less head per unit of flow. Area alone does not set the split.",
    nodes: [
      N("Tank High", "source", 30, 10, 32), N("Tank Low", "sink", 10, 90, 32)
    ],
    pipes: [
      P("Pipe 1", "Tank High", "Tank Low", 0.10, 100, 1.0),
      P("Pipe 2", "Tank High", "Tank Low", 0.15, 100, 1.0)
    ]
  },
  "gravity": {
    title: "Example 3 — Gravity feed, no pump",
    fluid: "Water, 20 C (check case)",
    blurb: "A roof tank feeds two outlets at different heights through one junction. Nothing but elevation drives this network.",
    lesson: "The low outlet draws more than the high outlet even though its pipe run is not the shortest path on the page: what matters is the head difference each branch actually gets.",
    nodes: [
      N("Roof Tank", "source", 40, 10, 20), N("Junction J", "junction", null, 50, 32),
      N("High Outlet", "sink", 15, 90, 16), N("Low Outlet", "sink", 5, 90, 48)
    ],
    pipes: [
      P("Pipe 1", "Roof Tank", "Junction J", 0.12, 60, 1.5),
      P("Pipe 2", "Junction J", "High Outlet", 0.08, 40, 2.0),
      P("Pipe 3", "Junction J", "Low Outlet", 0.08, 30, 2.0)
    ]
  },
  "starved": {
    title: "Example 4 — The starved branch",
    fluid: "Crude oil (preset)",
    blurb: "A weak pump lifts crude from a low supply tank to a junction feeding a high unit and a low unit. At full speed the junction head sits just below what the high unit needs.",
    lesson: "At speed ratio 1.00 the high unit is starved and its pipe runs slightly backward. Drag the pump speed to 1.15, then 1.30, and watch the branch come alive. This is the operating-point lesson: a pump and a network settle together.",
    nodes: [
      N("Supply Tank", "source", 10, 10, 32), N("Junction J", "junction", null, 55, 32),
      N("High Unit", "sink", 34, 90, 16), N("Low Unit", "sink", 12, 90, 48)
    ],
    pipes: [
      P("Pipe 1", "Supply Tank", "Junction J", 0.15, 80, 1.5, "Pump W (weak)"),
      P("Pipe 2", "Junction J", "High Unit", 0.10, 60, 2.0),
      P("Pipe 3", "Junction J", "Low Unit", 0.10, 60, 2.0)
    ]
  },
  "laminar": {
    title: "Example 5 — Laminar flow, viscous oil",
    fluid: "Viscous oil (preset)",
    blurb: "A thick oil creeps through two pipes in series. The Reynolds number stays near 18, deep in the laminar regime, so this whole network can also be solved exactly by hand.",
    lesson: "In laminar flow the head loss grows in direct proportion to the velocity, not its square, and the solver reproduces the hand result from the Hagen-Poiseuille equation. This is the teaching and checking case.",
    nodes: [
      N("Tank A", "source", 14, 10, 32), N("Middle M", "junction", null, 50, 32),
      N("Tank B", "sink", 10, 90, 32)
    ],
    pipes: [
      P("Pipe 1", "Tank A", "Middle M", 0.06, 40, 0.5),
      P("Pipe 2", "Middle M", "Tank B", 0.06, 60, 0.5)
    ]
  },
  "topfill": {
    title: "Example 6 — Top fill: two source tanks to one receiving tank",
    fluid: "Crude oil (preset)",
    blurb: "The team's final demo case. Two source tanks at different levels feed one junction; a loading pump pushes the combined flow through a fill line that discharges over the top of an open receiving tank, above its max fill level. The receiving tank node's head (21 m) is the discharge elevation at the top of the tank, not the liquid level: the lift is measured from the network up to the top of the tank. Use the tank-fill panel below the results to solve for fill time, for the flow a target time requires, or for a diameter.",
    lesson: "Because the discharge is over the top at a fixed elevation, the fill flow does not decay as the tank fills: the level never pushes back. The higher source tank contributes about three times what the lower one does, and the fill time is an exact volume-over-flow calculation, not an estimate.",
    fill: { tank: "Receiving Tank", tankD: 3, bottom: 15, start: 1.0, max: 5.2 },
    nodes: [
      N("Source A", "source", 14, 8, 16), N("Source B", "source", 11, 8, 50),
      N("Junction J", "junction", null, 42, 32),
      N("Receiving Tank", "sink", 21, 90, 32)
    ],
    pipes: [
      P("A Feed", "Source A", "Junction J", 0.12, 90, 1.2),
      P("B Feed", "Source B", "Junction J", 0.12, 110, 1.2),
      P("Fill Line", "Junction J", "Receiving Tank", 0.15, 160, 2.45, "Pump L (loading)")
    ]
  }
};

if (typeof module !== "undefined") module.exports = { EXAMPLES, FLUIDS, PUMPS };
