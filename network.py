# Refinery pipe network - reference calculation core (Python, verification standard).
# Iteration 1, November 19, 2026. Minimal code: numpy + scipy library routines only.
#
# Symbol key (SI units throughout):
#   H    total head at a node (metres above the shared datum)
#   Q    volumetric flow rate in a pipe (cubic metres per second); signed,
#        positive from the pipe's start node to its end node
#   D    pipe inside diameter (metres)
#   L    pipe length (metres)
#   eps  pipe wall roughness (metres)
#   K    sum of minor-loss coefficients of a pipe's fittings (dimensionless)
#        [CLASS-W13] values from a standard textbook table [R4]
#   rho  fluid density (kilograms per cubic metre)
#   mu   fluid dynamic viscosity (pascal seconds)
#   V    mean velocity in a pipe (metres per second)
#   Re   Reynolds number (dimensionless) = rho * V * D / mu
#   fF   Fanning friction factor (dimensionless); laminar fF = 16 / Re
#        [CLASS-HW4]; turbulent from the Churchill correlation [R3],
#        stated in Darcy form and divided by 4 (Darcy factor = 4 * Fanning)
#        [CLASS-W11]
#   hL   head loss in a pipe (metres) = (4 * fF * L / D + K) * V * |V| / (2 * g)
#        [CLASS-W10/W11 engineering Bernoulli with frictional heating;
#         CLASS-W13 fitting term]
#   g    gravitational acceleration (metres per second squared) = 9.81
#   H0   pump shutoff head at reference speed (metres) - stated preset [R7]
#   kp   pump curve coefficient (seconds squared per metre to the fifth):
#        pump head gain = s^2 * H0 - kp * Q * |Q|
#   s    pump speed ratio (dimensionless, 1 = reference speed) [NEW-EXT affinity]
#   eta  pump efficiency (dimensionless), stated constant preset = 0.75 [R7]
#   Pp   pump hydraulic power delivered to the fluid divided by eta, i.e.
#        shaft power (watts) = rho * g * Q * (pump head gain) / eta [CLASS-W11]
#
# Method [NEW-EXT]: node-head method. Unknowns are the heads at free
# junctions. Each pipe's flow comes from its end heads by a scalar root
# find (scipy brentq); junction flow balance is solved by scipy fsolve.
# Flow balance at a junction is continuity [CLASS pipe-flow weeks].

import numpy as np
from scipy.optimize import brentq, fsolve

g = 9.81


def fanning_f(Re, D, eps):
    # [CLASS-HW4] laminar limit; [R3] Churchill (1977) for all other regimes.
    if Re <= 0:
        return 0.0
    if Re < 2300:
        return 16.0 / Re
    a = (2.457 * np.log(1.0 / ((7.0 / Re) ** 0.9 + 0.27 * eps / D))) ** 16
    b = (37530.0 / Re) ** 16
    f_darcy = 8.0 * ((8.0 / Re) ** 12 + 1.0 / (a + b) ** 1.5) ** (1.0 / 12.0)
    return f_darcy / 4.0  # [CLASS-W11] Darcy factor = 4 * Fanning factor


def pipe_flow(dH, pipe, rho, mu):
    # dH is start-node head minus end-node head (metres), pump excluded.
    D, L, eps, K = pipe["D"], pipe["L"], pipe["eps"], pipe["K"]
    A = np.pi * D ** 2 / 4.0

    def resid(Q):
        V = Q / A
        Re = rho * abs(V) * D / mu
        fF = fanning_f(max(Re, 1e-9), D, eps)
        loss = (4.0 * fF * L / D + K) * V * abs(V) / (2 * g)  # V*|V| carries the sign
        gain = 0.0
        if pipe.get("pump"):
            p = pipe["pump"]
            gain = p["s"] ** 2 * p["H0"] - p["kp"] * Q * abs(Q)
        return dH + gain - loss

    lo, hi = -3.0, 3.0
    if resid(lo) * resid(hi) > 0:
        return 0.0, False  # no root in bracket: caller must flag, never hide
    return brentq(resid, lo, hi, xtol=1e-12), True


def solve(nodes, pipes, rho, mu):
    # nodes: {name: fixed head in metres, or None for a free junction}
    # pipes: [{name, a, b, D, L, eps, K, pump: None or {H0, kp, s, eta}}]
    free = [n for n, h in nodes.items() if h is None]
    idx = {n: i for i, n in enumerate(free)}

    def heads(x):
        return {n: (nodes[n] if nodes[n] is not None else x[idx[n]]) for n in nodes}

    def balance(x):
        H = heads(x)
        r = np.zeros(len(free))
        for p in pipes:
            Q, ok = pipe_flow(H[p["a"]] - H[p["b"]], p, rho, mu)
            if p["a"] in idx:
                r[idx[p["a"]]] -= Q
            if p["b"] in idx:
                r[idx[p["b"]]] += Q
        return r

    fixed = [h for h in nodes.values() if h is not None]
    guess = np.full(len(free), float(np.mean(fixed))) if free else np.zeros(0)
    x, info, ier, _ = fsolve(balance, guess, full_output=True)
    H = heads(x)
    out = {"ier": int(ier),
           "max_residual": float(np.max(np.abs(balance(x)))) if free else 0.0,
           "heads": {k: float(v) for k, v in H.items()}, "pipes": []}
    for p in pipes:
        Q, ok = pipe_flow(H[p["a"]] - H[p["b"]], p, rho, mu)
        A = np.pi * p["D"] ** 2 / 4.0
        V = Q / A
        Re = rho * abs(V) * p["D"] / mu
        fF = fanning_f(max(Re, 1e-9), p["D"], p["eps"])
        hL = (4.0 * fF * p["L"] / p["D"] + p["K"]) * V * abs(V) / (2 * g)
        row = {"name": p["name"], "a": p["a"], "b": p["b"], "Q": float(Q),
               "V": float(V), "Re": float(Re), "fF": float(fF), "hL": float(hL),
               "reversed": bool(Q < 0), "solved": bool(ok),
               "Leq": float(p["K"] * p["D"] / (4.0 * fF)) if fF > 0 else 0.0}
        if p.get("pump"):
            pu = p["pump"]
            gain = pu["s"] ** 2 * pu["H0"] - pu["kp"] * Q * abs(Q)
            row["pump_head"] = float(gain)
            row["pump_power_W"] = float(rho * g * Q * gain / pu["eta"]) if Q > 0 else 0.0
        out["pipes"].append(row)
    return out


# ---- Shared preset (identical numbers drive the website) ----
# Heads are liquid-surface or required-inlet total heads (metres).
PRESET_NODES = {"Tank A": 25.0, "Tank B": 15.0, "Junction J": None,
                "Split S": None, "Unit 1": 18.0, "Unit 2": 8.0}

PUMPS = {  # stated illustrative presets [R7], not manufacturer data
    "Pump 1": {"H0": 40.0, "kp": 900.0, "s": 1.0, "eta": 0.75},
    "Pump 2": {"H0": 55.0, "kp": 1400.0, "s": 1.0, "eta": 0.75},
    "Pump 3": {"H0": 30.0, "kp": 600.0, "s": 1.0, "eta": 0.75},
    "Pump W (weak)": {"H0": 28.0, "kp": 700.0, "s": 1.0, "eta": 0.75},
    "Pump L (loading)": {"H0": 15.0, "kp": 400.0, "s": 1.0, "eta": 0.75},
}

STEEL_EPS = 4.5e-5  # commercial steel roughness (metres) [R5]


def preset_pipes(pump_name="Pump 1", speed=1.0):
    pump = dict(PUMPS[pump_name])
    pump["s"] = speed
    # K compositions [R4]: sharp entrance 0.5; gate valve, fully open, 0.15;
    # 90-degree flanged elbow 0.3; swing check valve 2.5; tee, branch flow, 1.8.
    return [
        {"name": "Pipe 1", "a": "Tank A", "b": "Junction J", "D": 0.20, "L": 120.0,
         "eps": STEEL_EPS, "K": 0.5 + 0.15 + 2 * 0.3, "pump": None},
        {"name": "Pipe 2", "a": "Tank B", "b": "Junction J", "D": 0.20, "L": 90.0,
         "eps": STEEL_EPS, "K": 0.5 + 0.15 + 0.3, "pump": None},
        {"name": "Pipe 3", "a": "Junction J", "b": "Split S", "D": 0.25, "L": 200.0,
         "eps": STEEL_EPS, "K": 0.3 + 2.5 + 0.15, "pump": pump},
        {"name": "Pipe 4", "a": "Split S", "b": "Unit 1", "D": 0.15, "L": 150.0,
         "eps": STEEL_EPS, "K": 1.8 + 2 * 0.3 + 0.15, "pump": None},
        {"name": "Pipe 5", "a": "Split S", "b": "Unit 2", "D": 0.15, "L": 110.0,
         "eps": STEEL_EPS, "K": 1.8 + 0.3 + 0.15, "pump": None},
    ]


FLUIDS = {  # crude is a stated demo preset [R7]; water is a standard table [R5]
    "Crude oil (preset)": {"rho": 850.0, "mu": 0.008},
    "Water, 20 C (check case)": {"rho": 998.0, "mu": 0.001002},
    "Viscous oil (preset)": {"rho": 900.0, "mu": 0.35},
}
