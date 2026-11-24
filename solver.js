// Refinery pipe network - browser solver. Iteration 1, November 19, 2026.
// Mirror of core/network.py; cross-checked against it in tests/crosscheck_js.js.
// Symbol key (SI units): H total head at a node (m); Q pipe flow (m3/s), signed
// start node to end node; D diameter (m); L length (m); eps roughness (m);
// K sum of fitting loss coefficients (dimensionless); rho density (kg/m3);
// mu viscosity (Pa s); V mean velocity (m/s); Re Reynolds number
// (dimensionless); fF Fanning friction factor (dimensionless);
// hL head loss (m); g = 9.81 m/s2; H0 pump shutoff head (m);
// kp pump curve coefficient (s2/m5); s pump speed ratio (dimensionless);
// eta pump efficiency (dimensionless).
"use strict";
const G = 9.81;

function fanningF(Re, D, eps) {
  if (Re <= 0) return 0;
  if (Re < 2300) return 16 / Re;
  const a = Math.pow(2.457 * Math.log(1 / (Math.pow(7 / Re, 0.9) + 0.27 * eps / D)), 16);
  const b = Math.pow(37530 / Re, 16);
  const fDarcy = 8 * Math.pow(Math.pow(8 / Re, 12) + 1 / Math.pow(a + b, 1.5), 1 / 12);
  return fDarcy / 4;
}

function pipeFlow(dH, p, rho, mu) {
  const A = Math.PI * p.D * p.D / 4;
  function resid(Q) {
    const V = Q / A;
    const Re = rho * Math.abs(V) * p.D / mu;
    const fF = fanningF(Math.max(Re, 1e-9), p.D, p.eps);
    const loss = (4 * fF * p.L / p.D + p.K) * V * Math.abs(V) / (2 * G);
    let gain = 0;
    if (p.pump) gain = p.pump.s * p.pump.s * p.pump.H0 - p.pump.kp * Q * Math.abs(Q);
    return dH + gain - loss;
  }
  let lo = -3, hi = 3, flo = resid(lo), fhi = resid(hi);
  if (flo === 0) return { Q: lo, ok: true };
  if (flo * fhi > 0) return { Q: 0, ok: false };
  for (let i = 0; i < 200; i++) {
    const mid = (lo + hi) / 2, fm = resid(mid);
    if (flo * fm <= 0) { hi = mid; } else { lo = mid; flo = fm; }
    if (hi - lo < 1e-13) break;
  }
  return { Q: (lo + hi) / 2, ok: true };
}

function solveNetwork(nodes, pipes, rho, mu) {
  // nodes: {name: fixed head (m) or null}; pipes as in core/network.py preset.
  const names = Object.keys(nodes);
  const free = names.filter(n => nodes[n] === null);
  const idx = {}; free.forEach((n, i) => idx[n] = i);
  const fixedVals = names.filter(n => nodes[n] !== null).map(n => nodes[n]);
  const meanFixed = fixedVals.reduce((a, b) => a + b, 0) / fixedVals.length;

  function heads(x) {
    const H = {};
    names.forEach(n => H[n] = nodes[n] !== null ? nodes[n] : x[idx[n]]);
    return H;
  }
  function balance(x) {
    const H = heads(x), r = new Array(free.length).fill(0);
    for (const p of pipes) {
      const q = pipeFlow(H[p.a] - H[p.b], p, rho, mu).Q;
      if (p.a in idx) r[idx[p.a]] -= q;
      if (p.b in idx) r[idx[p.b]] += q;
    }
    return r;
  }
  function solveLinear(A, b) { // Gaussian elimination, small systems only
    const n = b.length, M = A.map((row, i) => row.concat(b[i]));
    for (let c = 0; c < n; c++) {
      let piv = c;
      for (let r = c + 1; r < n; r++) if (Math.abs(M[r][c]) > Math.abs(M[piv][c])) piv = r;
      if (Math.abs(M[piv][c]) < 1e-14) return null;
      [M[c], M[piv]] = [M[piv], M[c]];
      for (let r = 0; r < n; r++) {
        if (r === c) continue;
        const f = M[r][c] / M[c][c];
        for (let k = c; k <= n; k++) M[r][k] -= f * M[c][k];
      }
    }
    return M.map((row, i) => row[n] / row[i]);
  }

  let x = free.map(() => meanFixed), converged = free.length === 0;
  let r = balance(x);
  for (let iter = 0; iter < 60 && free.length; iter++) {
    const maxr = Math.max(...r.map(Math.abs));
    if (maxr < 1e-10) { converged = true; break; }
    const J = [];
    for (let j = 0; j < free.length; j++) {
      const xj = x.slice(); const h = 1e-6 * Math.max(1, Math.abs(x[j]));
      xj[j] += h;
      const rj = balance(xj);
      J.push(r.map((v, i) => (rj[i] - v) / h));
    }
    const Jt = free.map((_, i) => free.map((_, j) => J[j][i])); // transpose to rows
    const step = solveLinear(Jt, r.map(v => -v));
    if (!step) break;
    let t = 1, improved = false;
    for (let ls = 0; ls < 20; ls++) {
      const xn = x.map((v, i) => v + t * step[i]);
      const rn = balance(xn);
      if (Math.max(...rn.map(Math.abs)) < maxr) { x = xn; r = rn; improved = true; break; }
      t /= 2;
    }
    if (!improved) break;
  }
  const maxRes = free.length ? Math.max(...r.map(Math.abs)) : 0;
  if (free.length && maxRes < 1e-8) converged = true;

  const H = heads(x);
  const rows = pipes.map(p => {
    const q = pipeFlow(H[p.a] - H[p.b], p, rho, mu);
    const A = Math.PI * p.D * p.D / 4, V = q.Q / A;
    const Re = rho * Math.abs(V) * p.D / mu;
    const fF = fanningF(Math.max(Re, 1e-9), p.D, p.eps);
    const hL = (4 * fF * p.L / p.D + p.K) * V * Math.abs(V) / (2 * G);
    const row = { name: p.name, a: p.a, b: p.b, Q: q.Q, V, Re, fF, hL,
      reversed: q.Q < 0, solved: q.ok, Leq: fF > 0 ? p.K * p.D / (4 * fF) : 0 };
    if (p.pump) {
      const gain = p.pump.s * p.pump.s * p.pump.H0 - p.pump.kp * q.Q * Math.abs(q.Q);
      row.pump_head = gain;
      row.pump_power_W = q.Q > 0 ? rho * G * q.Q * gain / p.pump.eta : 0;
    }
    return row;
  });
  return { converged, max_residual: maxRes, heads: H, pipes: rows };
}

// ---- Adapter: builder state (examples.js) -> solver inputs ----
// Nodes arrive as a list with kind and head; junctions become null (free).
// Pipes arrive with a pump preset name and speed ratio; the solver wants
// the pump object itself. PUMPS and the examples live in examples.js.
function buildInputs(nodeList, pipeList) {
  const nodes = {};
  for (const n of nodeList) nodes[n.name] = n.kind === "junction" ? null : n.head;
  const pipes = pipeList.map(p => {
    let pump = null;
    if (p.pumpName) {
      const base = PUMPS[p.pumpName];
      pump = { H0: base.H0, kp: base.kp, s: p.s || 1, eta: base.eta };
    }
    return { name: p.name, a: p.a, b: p.b, D: p.D, L: p.L, eps: p.eps, K: p.K, pump };
  });
  return { nodes, pipes };
}
if (typeof module !== "undefined") module.exports = { solveNetwork, buildInputs, fanningF };
