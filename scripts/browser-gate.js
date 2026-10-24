#!/usr/bin/env node
// The browser gates: run the generator in headless Chromium (Playwright) and check it
// against the design's budgets and the golden hashes stored by the Node tests.
//   node scripts/browser-gate.js [seed] [--only=planet|region|zoom] [--tiles=3]
// Needs the `playwright` package and a Chromium it can launch.
import { createServer } from 'node:http';
import { readFile } from 'node:fs/promises';
import { readFileSync } from 'node:fs';
import { extname, join, normalize, resolve } from 'node:path';
import { fileURLToPath } from 'node:url';
import { createRequire } from 'node:module';

const root = resolve(fileURLToPath(new URL('..', import.meta.url)));
const types = { '.html': 'text/html', '.js': 'text/javascript', '.json': 'application/json' };
const server = createServer(async (req, res) => {
  let path = normalize(decodeURIComponent(new URL(req.url, 'http://x').pathname));
  if (path.endsWith('/')) path += 'index.html';
  const file = join(root, path);
  if (!file.startsWith(root)) return res.writeHead(403).end();
  try {
    const body = await readFile(file);
    res.writeHead(200, { 'content-type': types[extname(file)] || 'application/octet-stream' });
    res.end(body);
  } catch { res.writeHead(404).end(); }
});
await new Promise((ok) => server.listen(0, ok));
const port = server.address().port;

const require = createRequire(import.meta.url);
let chromium;
try { ({ chromium } = require('playwright')); } catch { ({ chromium } = await import('playwright')); }
const browser = await chromium.launch({ args: ['--enable-unsafe-webgpu'] });
const args = process.argv.slice(2);
const seed = args.find((a) => !a.startsWith('--')) || 'gate';
const only = args.find((a) => a.startsWith('--only='))?.slice(7);
const tiles = Number(args.find((a) => a.startsWith('--tiles='))?.slice(8) || 3);

async function runPage(path) {
  const page = await browser.newPage();
  page.on('pageerror', (e) => console.error('page error:', e.message));
  await page.goto(`http://localhost:${port}/web/gate/${path}`);
  await page.waitForFunction(() => window.gateResult, null, { timeout: 1800000 });
  const r = await page.evaluate(() => window.gateResult);
  await page.close();
  if (r.error) { console.error(r.error); process.exitCode = 1; }
  return r;
}

const golden = JSON.parse(readFileSync(join(root, 'test/golden.json'), 'utf8'));
let ok = true;
const check = (pass, msg) => { console.log(`${pass ? 'PASS' : 'FAIL'}  ${msg}`); ok &&= pass; };

if (!only || only === 'planet') {
  console.log('-- milestone 1: planet');
  const r = await runPage(`?seed=${encodeURIComponent(seed)}`);
  if (!r.error) {
    check(r.seconds < 10, `level 0 in a browser worker: ${r.seconds.toFixed(1)} s for ${r.cells} cells, ${r.events} events (budget 10 s)`);
    const same = JSON.stringify(r.golden) === JSON.stringify(golden.l0.hashes);
    check(same, `golden level-0 hashes in Chromium ${same ? 'match' : 'differ from'} Node's${same ? '' : `: ${JSON.stringify(r.golden)}`}`);
  } else ok = false;
}

if (!only || only === 'region') {
  console.log('-- milestone 2: region');
  const r = await runPage(`region.html?seed=${encodeURIComponent(seed)}&tiles=${tiles}`);
  if (!r.error) {
    const soft = r.adapter.architecture === 'swiftshader';
    console.log(`      WebGPU adapter: ${r.adapter.vendor} ${r.adapter.architecture}${soft ? ' (software emulation on the CPU: times say nothing about a real GPU)' : ''}`);
    for (const t of r.results) {
      console.log(`      tile ${t.tile} (planet ${t.planetElevation.toFixed(0)} m): ${t.seconds.toFixed(2)} s, Hack ${t.hack.toFixed(2)} (r2 ${t.hackR2.toFixed(2)}), ` +
        `concavity ${t.concavity.toFixed(2)} (r2 ${t.concavityR2.toFixed(2)}), parent block means within ${t.consistency.toFixed(1)} m, ${t.riverNodes} river nodes, Strahler order ${t.maxOrder}`);
      if (t.cpu) console.log(`      same tile on the CPU reference: ${t.cpuSeconds.toFixed(1)} s, Hack ${t.cpu.hack.toFixed(2)}, concavity ${t.cpu.concavity.toFixed(2)}; GPU-CPU elevation difference ${t.gpuCpuRms.toFixed(2)} m rms over ${t.relief.toFixed(0)} m relief`);
    }
    const worst = Math.max(...r.results.map((t) => t.seconds));
    check(worst < 1, `each tile under 1 s: slowest ${worst.toFixed(2)} s${soft ? ' [not measurable here: software GPU]' : ''}`);
    check(r.results.every((t) => t.hack > 0.45 && t.hack < 0.7 && t.hackR2 > 0.9), `Hack's law: exponents ${r.results.map((t) => t.hack.toFixed(2)).join(', ')} within 0.45-0.7 (r2 > 0.9)`);
    check(r.results.every((t) => t.concavity > 0.3 && t.concavity < 0.7 && t.concavityR2 > 0.6), `slope-area: concavities ${r.results.map((t) => t.concavity.toFixed(2)).join(', ')} within 0.3-0.7 (r2 > 0.6)`);
    const cmp = r.results.find((t) => t.cpu);
    if (cmp) check(cmp.gpuCpuRms < 0.02 * cmp.relief, `GPU matches the CPU reference: ${cmp.gpuCpuRms.toFixed(2)} m rms (< 2% of relief)`);
  } else ok = false;
}

if (!only || only === 'zoom') {
  console.log('-- milestone 5: zoom (planet to 1 m in the viewer)');
  const r = await runPage(`../world/?auto=1&store=0&seed=${encodeURIComponent(seed)}`);
  if (!r.error) {
    for (const st of r.steps) console.log(`      ${st.label.padEnd(28)} ${st.key.padEnd(18)} ${st.seconds.toFixed(2)} s  ${st.source}`);
    const by = (label) => r.steps.find((st) => st.label === label);
    check(by('site') && by('site').info.includes('trees'), 'zoomed from the planet to a 1.2-m site tile');
    check(by('site again (cache)').seconds < 0.1 && by('landscape again (cache)').seconds < 0.1, `revisits come from the cache (${by('site again (cache)').seconds.toFixed(3)} s, ${by('landscape again (cache)').seconds.toFixed(3)} s)`);
    const nb = by('site neighbour (prefetched)');
    check(nb.source.startsWith('from cache'), `a neighbouring site was prefetched (${nb.seconds.toFixed(3)} s, ${nb.source})`);
  } else ok = false;
}

await browser.close();
server.close();
process.exit(ok && !process.exitCode ? 0 : 1);
