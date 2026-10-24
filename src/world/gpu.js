// The parallel stream-power solver of jacobi.js on WebGPU.
//
// Each pass of jacobi.js is one compute kernel here, with the same arithmetic in
// 32-bit floats. The results agree with the CPU reference to rounding (WGSL's pow and
// division are not required to round exactly like JavaScript's), not bitwise; see
// docs/PLAN.md on what that means for determinism.
//
// Ping-pong passes (fill, area, erosion sweeps) alternate between two buffers; their
// sweep counts must be even, so every result lands back in the first buffer of its pair.

import { OFF8, OPPOSITE, FLAG_OUTLET, FLAG_FIXED, BIG, JACOBI_DEFAULTS } from './jacobi.js';

const WG = 8;

function shaderSource() {
  const off = OFF8.map(([a, b]) => `vec2<i32>(${a}, ${b})`).join(', ');
  const opp = OPPOSITE.join(', ');
  return /* wgsl */ `
struct Params { nx: u32, ny: u32, dx: f32, dt: f32, m: f32, eps: f32, kd: f32, big: f32 };
@group(0) @binding(0) var<uniform> p: Params;
@group(0) @binding(1) var<storage, read_write> z: array<f32>;
@group(0) @binding(2) var<storage, read_write> w: array<f32>;
@group(0) @binding(3) var<storage, read_write> w2: array<f32>;
@group(0) @binding(4) var<storage, read_write> rec: array<i32>;
@group(0) @binding(5) var<storage, read_write> A: array<f32>;
@group(0) @binding(6) var<storage, read_write> A2: array<f32>;
@group(0) @binding(7) var<storage, read_write> zr: array<f32>;
@group(0) @binding(8) var<storage, read_write> zr2: array<f32>;
@group(0) @binding(9) var<storage, read> U: array<f32>;
@group(0) @binding(10) var<storage, read> K: array<f32>;
@group(0) @binding(11) var<storage, read> R: array<f32>;
@group(0) @binding(12) var<storage, read> flags: array<u32>;

const OFF = array<vec2<i32>, 8>(${off});
const OPP = array<i32, 8>(${opp});
const OUTLET: u32 = ${FLAG_OUTLET}u;
const FIXED: u32 = ${FLAG_FIXED}u;

fn inside(i: i32, j: i32) -> bool { return i >= 0 && j >= 0 && i < i32(p.nx) && j < i32(p.ny); }
fn idx(i: i32, j: i32) -> u32 { return u32(j) * p.nx + u32(i); }
fn held(c: u32) -> bool { return (flags[c] & (OUTLET | FIXED)) != 0u; }
fn blocked(q: u32) -> bool { return (flags[q] & FIXED) != 0u && (flags[q] & OUTLET) == 0u; }
fn dist(k: i32) -> f32 { if (k < 4) { return p.dx; } return p.dx * 1.41421356; }

@compute @workgroup_size(${WG}, ${WG})
fn uplift(@builtin(global_invocation_id) g: vec3<u32>) {
  if (g.x >= p.nx || g.y >= p.ny) { return; }
  let c = g.y * p.nx + g.x;
  if (held(c)) { return; }
  let up = U[c] * p.dt;
  z[c] = z[c] + up;
  if (w[c] < p.big) { w[c] = w[c] + up; }
}

@compute @workgroup_size(${WG}, ${WG})
fn fill(@builtin(global_invocation_id) g: vec3<u32>) {
  if (g.x >= p.nx || g.y >= p.ny) { return; }
  let c = g.y * p.nx + g.x;
  if (held(c)) { w2[c] = z[c]; return; }
  var lo = w[c];
  for (var k = 0; k < 8; k++) {
    let a = i32(g.x) + OFF[k].x; let b = i32(g.y) + OFF[k].y;
    if (!inside(a, b)) { continue; }
    let q = idx(a, b);
    if (blocked(q)) { continue; }
    lo = min(lo, w[q] + p.eps * dist(k));
  }
  w2[c] = max(z[c], lo);
}

@compute @workgroup_size(${WG}, ${WG})
fn receivers(@builtin(global_invocation_id) g: vec3<u32>) {
  if (g.x >= p.nx || g.y >= p.ny) { return; }
  let c = g.y * p.nx + g.x;
  var r = -1;
  if (!held(c) && w[c] < p.big) {
    var best = 0.0;
    for (var k = 0; k < 8; k++) {
      let a = i32(g.x) + OFF[k].x; let b = i32(g.y) + OFF[k].y;
      if (!inside(a, b)) { continue; }
      let q = idx(a, b);
      if (blocked(q)) { continue; }
      let s = (w[c] - w[q]) / dist(k);
      if (s > best) { best = s; r = k; }
    }
  }
  rec[c] = r;
}

@compute @workgroup_size(${WG}, ${WG})
fn area(@builtin(global_invocation_id) g: vec3<u32>) {
  if (g.x >= p.nx || g.y >= p.ny) { return; }
  let c = g.y * p.nx + g.x;
  var s = p.dx * p.dx * R[c];
  for (var k = 0; k < 8; k++) {
    let a = i32(g.x) + OFF[k].x; let b = i32(g.y) + OFF[k].y;
    if (!inside(a, b)) { continue; }
    let q = idx(a, b);
    if (rec[q] == OPP[k]) { s = s + A[q]; }
  }
  A2[c] = s;
}

@compute @workgroup_size(${WG}, ${WG})
fn erodeInit(@builtin(global_invocation_id) g: vec3<u32>) {
  if (g.x >= p.nx || g.y >= p.ny) { return; }
  let c = g.y * p.nx + g.x;
  if (held(c) || rec[c] < 0) { zr[c] = w[c]; } else { zr[c] = w[c] - U[c] * p.dt; }
}

@compute @workgroup_size(${WG}, ${WG})
fn erode(@builtin(global_invocation_id) g: vec3<u32>) {
  if (g.x >= p.nx || g.y >= p.ny) { return; }
  let c = g.y * p.nx + g.x;
  let k = rec[c];
  if (held(c) || k < 0) { zr2[c] = z[c]; return; }
  let r = idx(i32(g.x) + OFF[k].x, i32(g.y) + OFF[k].y);
  let F = K[c] * p.dt * pow(A[c], p.m) / dist(k);
  zr2[c] = (w[c] + F * zr[r]) / (1.0 + F);
}

@compute @workgroup_size(${WG}, ${WG})
fn creep(@builtin(global_invocation_id) g: vec3<u32>) {
  if (g.x >= p.nx || g.y >= p.ny) { return; }
  let c = g.y * p.nx + g.x;
  if (held(c) || g.x == 0u || g.y == 0u || g.x == p.nx - 1u || g.y == p.ny - 1u || p.kd == 0.0) { z[c] = zr[c]; return; }
  z[c] = zr[c] + p.kd * ((zr[c - 1u] + zr[c + 1u]) + (zr[c - p.nx] + zr[c + p.nx]) - 4.0 * zr[c]);
}

@compute @workgroup_size(${WG}, ${WG})
fn restart(@builtin(global_invocation_id) g: vec3<u32>) {
  if (g.x >= p.nx || g.y >= p.ny) { return; }
  let c = g.y * p.nx + g.x;
  if (held(c) || rec[c] >= 0) { w[c] = z[c]; } else { w[c] = p.big; }
}
`;
}

// bindings each kernel uses (layout 'auto' includes only these)
const USES = {
  uplift: [0, 1, 2, 9, 12],
  fill: [0, 1, 2, 3, 12],
  receivers: [0, 2, 4, 12],
  area: [0, 4, 5, 6, 11],
  erodeInit: [0, 2, 4, 7, 9, 12],
  erode: [0, 1, 2, 4, 5, 7, 8, 10, 12],
  creep: [0, 1, 7, 12],
  restart: [0, 1, 2, 4, 12],
};

/** WebGPU device, or null where WebGPU is unavailable. */
export async function getDevice() {
  if (typeof navigator === 'undefined' || !navigator.gpu) return null;
  const adapter = await navigator.gpu.requestAdapter();
  if (!adapter) return null;
  return adapter.requestDevice();
}

/**
 * A GPU solver with the interface runRegion expects (see region.js).
 * @param {GPUDevice} device
 * @param {object} input regionInput(...)
 */
export function gpuSolver(device, input) {
  const { nx, ny, dx } = input;
  const n = nx * ny;
  const usage = GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_SRC | GPUBufferUsage.COPY_DST;
  const buf = (data, type = Float32Array) => {
    const b = device.createBuffer({ size: n * 4, usage });
    if (data) device.queue.writeBuffer(b, 0, type.from(data));
    return b;
  };
  const wInit = new Float32Array(n);
  for (let c = 0; c < n; c++) wInit[c] = input.flags[c] & (FLAG_OUTLET | FLAG_FIXED) ? input.z[c] : BIG;
  const B = {
    z: buf(input.z), wA: buf(wInit), wB: buf(), rec: buf(new Int32Array(n).fill(-1), Int32Array),
    aA: buf(new Float32Array(n)), aB: buf(), zA: buf(), zB: buf(),
    U: buf(input.U), K: buf(input.K), R: buf(input.runoff), flags: buf(input.flags, Uint32Array),
  };
  const params = device.createBuffer({ size: 32, usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST });
  const module = device.createShaderModule({ code: shaderSource() });
  const pipes = {};
  for (const name of Object.keys(USES)) pipes[name] = device.createComputePipeline({ layout: 'auto', compute: { module, entryPoint: name } });
  const slots = (o) => ({ 0: params, 1: B.z, 2: o.w, 3: o.w2, 4: B.rec, 5: o.A, 6: o.A2, 7: o.zr, 8: o.zr2, 9: B.U, 10: B.K, 11: B.R, 12: B.flags });
  const base = { w: B.wA, w2: B.wB, A: B.aA, A2: B.aB, zr: B.zA, zr2: B.zB };
  const AB = slots(base);
  // the flipped group swaps only the pair a ping-pong kernel writes
  const FLIP = { fill: { w: B.wB, w2: B.wA }, area: { A: B.aB, A2: B.aA }, erode: { zr: B.zB, zr2: B.zA } };
  const group = (name, map) => device.createBindGroup({
    layout: pipes[name].getBindGroupLayout(0),
    entries: USES[name].map((b) => ({ binding: b, resource: { buffer: map[b] } })),
  });
  const groups = {};
  for (const name of Object.keys(USES)) groups[name] = [group(name, AB), group(name, slots({ ...base, ...(FLIP[name] || {}) }))];
  const gx = Math.ceil(nx / WG), gy = Math.ceil(ny / WG);
  let iterations = 0;
  let lastParams = '';

  const setParams = (p) => {
    const key = JSON.stringify(p);
    if (key === lastParams) return;
    lastParams = key;
    const data = new ArrayBuffer(32);
    const u = new Uint32Array(data), fl = new Float32Array(data);
    u[0] = nx; u[1] = ny;
    fl[2] = dx; fl[3] = p.dt; fl[4] = p.m; fl[5] = p.eps; fl[6] = (p.D * p.dt) / (dx * dx); fl[7] = BIG;
    device.queue.writeBuffer(params, 0, data);
  };

  const read = async (buffer, Type = Float32Array) => {
    const staging = device.createBuffer({ size: n * 4, usage: GPUBufferUsage.MAP_READ | GPUBufferUsage.COPY_DST });
    const enc = device.createCommandEncoder();
    enc.copyBufferToBuffer(buffer, 0, staging, 0, n * 4);
    device.queue.submit([enc.finish()]);
    await staging.mapAsync(GPUMapMode.READ);
    const out = new Type(staging.getMappedRange().slice(0));
    staging.unmap();
    staging.destroy();
    return out;
  };

  return {
    async run(count, p0) {
      const p = { ...JACOBI_DEFAULTS, ...p0 };
      for (const k of ['fillSweeps', 'areaSweeps', 'sweeps']) if (p[k] % 2) throw new Error(`gpuSolver: ${k} must be even`);
      setParams(p);
      const enc = device.createCommandEncoder();
      const pass = enc.beginComputePass();
      const go = (name, flip = 0) => { pass.setPipeline(pipes[name]); pass.setBindGroup(0, groups[name][flip]); pass.dispatchWorkgroups(gx, gy); };
      for (let it = 0; it < count; it++) {
        go('uplift');
        for (let s = 0; s < p.fillSweeps; s++) go('fill', s % 2);
        go('receivers');
        for (let s = 0; s < p.areaSweeps; s++) go('area', s % 2);
        go('erodeInit');
        for (let s = 0; s < p.sweeps; s++) go('erode', s % 2);
        go('creep');
        go('restart');
      }
      pass.end();
      device.queue.submit([enc.finish()]);
      iterations += count;
    },
    readZ: () => read(B.z),
    setK: (K) => device.queue.writeBuffer(B.K, 0, Float32Array.from(K)),
    async readResult() {
      const [z, rec, A] = [await read(B.z), await read(B.rec, Int32Array), await read(B.aA)];
      return { z, w: z, rec, A, iterations };
    },
    destroy() { for (const b of Object.values(B)) b.destroy(); params.destroy(); },
  };
}
