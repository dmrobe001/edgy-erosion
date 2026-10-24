#!/usr/bin/env node
// The browser gates: run the generator in headless Chromium (Playwright) and check it
// against the design's budgets and the golden hashes stored by the Node tests.
//   node scripts/browser-gate.js [seed]
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
const page = await browser.newPage();
page.on('pageerror', (e) => console.error('page error:', e.message));
const seed = process.argv[2] || 'gate';
await page.goto(`http://localhost:${port}/web/gate/?seed=${encodeURIComponent(seed)}`);
await page.waitForFunction(() => window.gateResult, null, { timeout: 120000 });
const r = await page.evaluate(() => window.gateResult);
await browser.close();
server.close();

const golden = JSON.parse(readFileSync(join(root, 'test/golden.json'), 'utf8'));
let ok = true;
const check = (pass, msg) => { console.log(`${pass ? 'PASS' : 'FAIL'}  ${msg}`); ok &&= pass; };
if (r.error) { console.error(r.error); process.exit(1); }
check(r.seconds < 10, `level 0 in a browser worker: ${r.seconds.toFixed(1)} s for ${r.cells} cells, ${r.events} events (budget 10 s)`);
const same = JSON.stringify(r.golden) === JSON.stringify(golden.l0.hashes);
check(same, `golden level-0 hashes in Chromium ${same ? 'match' : 'differ from'} Node's${same ? '' : `: ${JSON.stringify(r.golden)}`}`);
process.exit(ok ? 0 : 1);
