// The graph solver of graphsolver.js on WebGPU: one compute kernel per pass, the same
// 32-bit arithmetic, neighbours from the same compressed rows. As with gpu.js, results
// agree with the CPU reference to rounding, not bitwise (see docs/PLAN.md).
//
// Ping-pong passes (fill, area, sediment, erosion sweeps) alternate between two buffers;
// their sweep counts must be even, so each result lands back in the first buffer.
//
// The interface is that of cpuGraphSolver (regionmesh.js): pass
// `(win) => gpuGraphSolver(device, win)` as the solver factory.

import { FLAG_OUTLET, FLAG_FIXED, BIG, CLASSES, GRAPH_DEFAULTS } from './graphsolver.js';

const WG = 64;

// WebGPU guarantees only 8 storage buffers per kernel, so the read-only data is packed:
//   graph  off (n + 1) | nbr (E) | dist (E) | lapw (E)      (floats as bit patterns)
//   node   U | P | area | inflow | flags | K                 (n each; K is rewritten by setK)
//   sed    sedIn (3n) | frac (3n)
//   R      rec, recDist interleaved
//   Qs     3n of flux, then n of deposition rate (the sediment pass writes both)
const SHADER = /* wgsl */ `
struct Params { n: u32, E: u32, dt: f32, m: f32, eps: f32, D: f32, big: f32, g0: f32, g1: f32, g2: f32 };
@group(0) @binding(0) var<uniform> p: Params;
@group(0) @binding(1) var<storage, read_write> z: array<f32>;
@group(0) @binding(2) var<storage, read_write> w: array<f32>;
@group(0) @binding(3) var<storage, read_write> w2: array<f32>;
@group(0) @binding(4) var<storage, read_write> R: array<i32>;
@group(0) @binding(5) var<storage, read_write> A: array<f32>;
@group(0) @binding(6) var<storage, read_write> A2: array<f32>;
@group(0) @binding(7) var<storage, read_write> Qs: array<f32>;
@group(0) @binding(8) var<storage, read_write> Qs2: array<f32>;
@group(0) @binding(9) var<storage, read_write> prod: array<f32>;
@group(0) @binding(10) var<storage, read_write> zr: array<f32>;
@group(0) @binding(11) var<storage, read_write> zr2: array<f32>;
@group(0) @binding(12) var<storage, read> graph: array<u32>;
@group(0) @binding(13) var<storage, read> node: array<u32>;
@group(0) @binding(14) var<storage, read> sed: array<f32>;

const OUTLET: u32 = ${FLAG_OUTLET}u;
const FIXED: u32 = ${FLAG_FIXED}u;
const CL: u32 = ${CLASSES}u;

fn off(i: u32) -> u32 { return graph[i]; }
fn nbr(e: u32) -> u32 { return graph[p.n + 1u + e]; }
fn dist(e: u32) -> f32 { return bitcast<f32>(graph[p.n + 1u + p.E + e]); }
fn lapw(e: u32) -> f32 { return bitcast<f32>(graph[p.n + 1u + 2u * p.E + e]); }
fn U(i: u32) -> f32 { return bitcast<f32>(node[i]); }
fn P(i: u32) -> f32 { return bitcast<f32>(node[p.n + i]); }
fn area(i: u32) -> f32 { return bitcast<f32>(node[2u * p.n + i]); }
fn inflow(i: u32) -> f32 { return bitcast<f32>(node[3u * p.n + i]); }
fn flags(i: u32) -> u32 { return node[4u * p.n + i]; }
fn K(i: u32) -> f32 { return bitcast<f32>(node[5u * p.n + i]); }
fn rec(i: u32) -> i32 { return R[2u * i]; }
fn recDist(i: u32) -> f32 { return bitcast<f32>(R[2u * i + 1u]); }
fn held(i: u32) -> bool { return (flags(i) & (OUTLET | FIXED)) != 0u; }
fn blocked(q: u32) -> bool { return (flags(q) & FIXED) != 0u && (flags(q) & OUTLET) == 0u; }
fn G(k: u32) -> f32 { if (k == 0u) { return p.g0; } if (k == 1u) { return p.g1; } return p.g2; }

@compute @workgroup_size(${WG})
fn uplift(@builtin(global_invocation_id) g: vec3<u32>) {
  let i = g.x; if (i >= p.n || held(i)) { return; }
  let up = U(i) * p.dt;
  z[i] = z[i] + up;
  if (w[i] < p.big) { w[i] = w[i] + up; }
}

@compute @workgroup_size(${WG})
fn fill(@builtin(global_invocation_id) g: vec3<u32>) {
  let i = g.x; if (i >= p.n) { return; }
  if (held(i)) { w2[i] = z[i]; return; }
  var lo = w[i];
  for (var e = off(i); e < off(i + 1u); e++) {
    let j = nbr(e);
    if (blocked(j)) { continue; }
    lo = min(lo, w[j] + p.eps * dist(e));
  }
  w2[i] = max(z[i], lo);
}

@compute @workgroup_size(${WG})
fn receivers(@builtin(global_invocation_id) g: vec3<u32>) {
  let i = g.x; if (i >= p.n) { return; }
  var r = -1; var rd = 0.0;
  if (!held(i) && w[i] < p.big) {
    var best = 0.0;
    for (var e = off(i); e < off(i + 1u); e++) {
      let j = nbr(e);
      if (blocked(j)) { continue; }
      let s = (w[i] - w[j]) / dist(e);
      if (s > best) { best = s; r = i32(j); rd = dist(e); }
    }
  }
  R[2u * i] = r; R[2u * i + 1u] = bitcast<i32>(rd);
}

@compute @workgroup_size(${WG})
fn gather(@builtin(global_invocation_id) g: vec3<u32>) {
  let i = g.x; if (i >= p.n) { return; }
  var a = area(i) * P(i) + inflow(i);
  for (var e = off(i); e < off(i + 1u); e++) { let j = nbr(e); if (rec(j) == i32(i)) { a = a + A[j]; } }
  A2[i] = a;
}

@compute @workgroup_size(${WG})
fn sediment(@builtin(global_invocation_id) g: vec3<u32>) {
  let i = g.x; if (i >= p.n) { return; }
  var d = 0.0;
  for (var k = 0u; k < CL; k++) {
    var q = sed[k * p.n + i];
    for (var e = off(i); e < off(i + 1u); e++) { let j = nbr(e); if (rec(j) == i32(i)) { q = q + Qs[k * p.n + j]; } }
    var laid = 0.0;
    if (!held(i)) { laid = q * (1.0 - exp(-G(k) * area(i) / max(A[i], area(i)))); }
    d = d + laid;
    Qs2[k * p.n + i] = (q - laid) + prod[k * p.n + i];
  }
  Qs2[CL * p.n + i] = d / area(i);
}

@compute @workgroup_size(${WG})
fn erodeInit(@builtin(global_invocation_id) g: vec3<u32>) {
  let i = g.x; if (i >= p.n) { return; }
  if (held(i) || rec(i) < 0) { zr[i] = w[i]; } else { zr[i] = w[i] - U(i) * p.dt; }
}

@compute @workgroup_size(${WG})
fn erode(@builtin(global_invocation_id) g: vec3<u32>) {
  let i = g.x; if (i >= p.n) { return; }
  let r = rec(i);
  if (held(i) || r < 0) { zr2[i] = z[i]; return; }
  let F = K(i) * p.dt * pow(A[i], p.m) / recDist(i);
  zr2[i] = ((w[i] + p.dt * Qs[CL * p.n + i]) + F * zr[u32(r)]) / (1.0 + F);
}

@compute @workgroup_size(${WG})
fn produce(@builtin(global_invocation_id) g: vec3<u32>) {
  let i = g.x; if (i >= p.n) { return; }
  let r = rec(i);
  var E = 0.0;
  if (!held(i) && r >= 0) { E = max(0.0, K(i) * pow(A[i], p.m) * (zr[i] - zr[u32(r)]) / recDist(i)); }
  for (var k = 0u; k < CL; k++) { prod[k * p.n + i] = E * area(i) * sed[(CL + k) * p.n + i]; }
}

@compute @workgroup_size(${WG})
fn creep(@builtin(global_invocation_id) g: vec3<u32>) {
  let i = g.x; if (i >= p.n) { return; }
  if (held(i) || p.D == 0.0) { z[i] = zr[i]; return; }
  var s = 0.0;
  for (var e = off(i); e < off(i + 1u); e++) { s = s + lapw(e) * (zr[nbr(e)] - zr[i]); }
  z[i] = zr[i] + (p.D * p.dt / area(i)) * s;
}

@compute @workgroup_size(${WG})
fn restart(@builtin(global_invocation_id) g: vec3<u32>) {
  let i = g.x; if (i >= p.n) { return; }
  if (held(i) || rec(i) >= 0) { w[i] = z[i]; } else { w[i] = p.big; }
}
`;

// bindings each kernel uses (layout 'auto' includes only these)
const USES = {
  uplift: [0, 1, 2, 13],
  fill: [0, 1, 2, 3, 12, 13],
  receivers: [0, 2, 4, 12, 13],
  gather: [0, 4, 5, 6, 12, 13],
  sediment: [0, 4, 5, 7, 8, 9, 12, 13, 14],
  erodeInit: [0, 2, 4, 10, 13],
  erode: [0, 1, 2, 4, 5, 7, 10, 11, 13],
  produce: [0, 4, 5, 9, 10, 13, 14],
  creep: [0, 1, 10, 12, 13],
  restart: [0, 1, 2, 4, 13],
};

/**
 * A WebGPU graph solver for a mesh window (regionmesh.js buildWindow, or any object
 * with g {n, off, nbr, dist, lapw, area} and the fields of GraphSolver).
 * @param {GPUDevice} device
 */
export function gpuGraphSolver(device, win) {
  const g = win.g, n = g.n, E = g.nbr.length;
  if (win.relax && win.relax.some((r) => r > 0)) throw new Error('gpuGraphSolver: relaxation zones (regionmesh relax option) run on the CPU only');
  const usage = GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_SRC | GPUBufferUsage.COPY_DST;
  const buf = (len, data, Type = Float32Array) => {
    const b = device.createBuffer({ size: Math.max(4, len * 4), usage });
    if (data) device.queue.writeBuffer(b, 0, Type.from(data));
    return b;
  };
  const wInit = new Float32Array(n);
  for (let i = 0; i < n; i++) wInit[i] = win.flags[i] & (FLAG_OUTLET | FLAG_FIXED) ? win.z[i] : BIG;
  const C = CLASSES * n;
  const graph = new Uint32Array(n + 1 + 3 * E);
  graph.set(g.off, 0); graph.set(g.nbr, n + 1);
  new Float32Array(graph.buffer, 4 * (n + 1 + E), E).set(g.dist);
  new Float32Array(graph.buffer, 4 * (n + 1 + 2 * E), E).set(g.lapw);
  const nodeData = new Uint32Array(6 * n), nf = new Float32Array(nodeData.buffer);
  nf.set(win.U, 0); nf.set(win.P, n); nf.set(g.area, 2 * n); if (win.inflow) nf.set(win.inflow, 3 * n);
  nodeData.set(win.flags, 4 * n); nf.set(win.K, 5 * n);
  const sedData = new Float32Array(2 * C);
  if (win.sedIn) sedData.set(win.sedIn, 0);
  sedData.set(win.frac, C);
  const B = {
    z: buf(n, win.z), wA: buf(n, wInit), wB: buf(n), R: buf(2 * n, new Int32Array(2 * n).fill(-1), Int32Array),
    aA: buf(n, new Float32Array(n)), aB: buf(n), qA: buf(C + n, new Float32Array(C + n)), qB: buf(C + n), prod: buf(C, new Float32Array(C)),
    zA: buf(n), zB: buf(n), graph: buf(graph.length, graph, Uint32Array), node: buf(6 * n, nodeData, Uint32Array), sed: buf(2 * C, sedData),
  };
  const params = device.createBuffer({ size: 48, usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST });
  const module = device.createShaderModule({ code: SHADER });
  const pipes = {};
  for (const name of Object.keys(USES)) pipes[name] = device.createComputePipeline({ layout: 'auto', compute: { module, entryPoint: name } });
  const slots = (o) => ({ 0: params, 1: B.z, 2: o.w, 3: o.w2, 4: B.R, 5: o.A, 6: o.A2, 7: o.Qs, 8: o.Qs2, 9: B.prod, 10: o.zr, 11: o.zr2, 12: B.graph, 13: B.node, 14: B.sed });
  const base = { w: B.wA, w2: B.wB, A: B.aA, A2: B.aB, Qs: B.qA, Qs2: B.qB, zr: B.zA, zr2: B.zB };
  const FLIP = { fill: { w: B.wB, w2: B.wA }, gather: { A: B.aB, A2: B.aA }, sediment: { Qs: B.qB, Qs2: B.qA }, erode: { zr: B.zB, zr2: B.zA } };
  const group = (name, map) => device.createBindGroup({
    layout: pipes[name].getBindGroupLayout(0),
    entries: USES[name].map((b) => ({ binding: b, resource: { buffer: map[b] } })),
  });
  const groups = {};
  for (const name of Object.keys(USES)) groups[name] = [group(name, slots(base)), group(name, slots({ ...base, ...(FLIP[name] || {}) }))];
  const groupsX = Math.ceil(n / WG);
  let iterations = 0, lastParams = '';
  const setParams = (q) => {
    const key = JSON.stringify(q);
    if (key === lastParams) return;
    lastParams = key;
    const data = new ArrayBuffer(48), u = new Uint32Array(data), fl = new Float32Array(data);
    u[0] = n; u[1] = E; fl[2] = q.dt; fl[3] = q.m; fl[4] = q.eps; fl[5] = q.D; fl[6] = BIG; fl[7] = q.G[0]; fl[8] = q.G[1]; fl[9] = q.G[2];
    device.queue.writeBuffer(params, 0, data);
  };
  const read = async (buffer, len, Type = Float32Array) => {
    const staging = device.createBuffer({ size: len * 4, usage: GPUBufferUsage.MAP_READ | GPUBufferUsage.COPY_DST });
    const enc = device.createCommandEncoder();
    enc.copyBufferToBuffer(buffer, 0, staging, 0, len * 4);
    device.queue.submit([enc.finish()]);
    await staging.mapAsync(GPUMapMode.READ);
    const out = new Type(staging.getMappedRange().slice(0));
    staging.unmap();
    staging.destroy();
    return out;
  };
  return {
    async run(count, p0) {
      const q = { ...GRAPH_DEFAULTS, ...p0 };
      for (const k of ['fillSweeps', 'areaSweeps', 'sedSweeps', 'sweeps']) if (q[k] % 2) throw new Error(`gpuGraphSolver: ${k} must be even`);
      setParams(q);
      const enc = device.createCommandEncoder();
      const pass = enc.beginComputePass();
      const go = (name, flip = 0) => { pass.setPipeline(pipes[name]); pass.setBindGroup(0, groups[name][flip]); pass.dispatchWorkgroups(groupsX); };
      for (let it = 0; it < count; it++) {
        go('uplift');
        for (let s = 0; s < q.fillSweeps; s++) go('fill', s % 2);
        go('receivers');
        for (let s = 0; s < q.areaSweeps; s++) go('gather', s % 2);
        for (let s = 0; s < q.sedSweeps; s++) go('sediment', s % 2);
        go('erodeInit');
        for (let s = 0; s < q.sweeps; s++) go('erode', s % 2);
        go('produce');
        go('creep');
        go('restart');
      }
      pass.end();
      device.queue.submit([enc.finish()]);
      iterations += count;
    },
    readZ: () => read(B.z, n),
    setK: (K) => device.queue.writeBuffer(B.node, 20 * n, Float32Array.from(K)),
    async readResult() {
      const Rr = await read(B.R, 2 * n, Int32Array), q = await read(B.qA, C + n);
      const rec = new Int32Array(n), recDist = new Float32Array(n), rd = new Float32Array(Rr.buffer);
      for (let i = 0; i < n; i++) { rec[i] = Rr[2 * i]; recDist[i] = rd[2 * i + 1]; }
      return { z: await read(B.z, n), rec, recDist, A: await read(B.aA, n), Qs: q.slice(0, C), dep: q.slice(C), iterations };
    },
    destroy() { for (const b of Object.values(B)) b.destroy(); params.destroy(); },
  };
}
