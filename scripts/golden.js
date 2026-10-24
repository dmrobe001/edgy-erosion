#!/usr/bin/env node
// Regenerate test/golden.json after a deliberate model change (and bump GENERATOR_VERSION).
//   node scripts/golden.js
import { readFileSync, writeFileSync } from 'node:fs';
import { generatePlanet, GENERATOR_VERSION } from '../src/world/l0.js';
import { fieldHash } from '../src/world/hash.js';

const file = new URL('../test/golden.json', import.meta.url);
const g = JSON.parse(readFileSync(file, 'utf8'));
const hashes = (t) => Object.fromEntries(Object.entries(t.fields).map(([k, v]) => [k, fieldHash(v)]));
g.version = GENERATOR_VERSION;
g.l0.hashes = hashes(generatePlanet(g.l0.seed, g.l0.options));
writeFileSync(file, JSON.stringify(g, null, 2) + '\n');
console.log(`wrote ${file.pathname} for ${GENERATOR_VERSION}`);
