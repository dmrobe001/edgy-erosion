#!/usr/bin/env node
// Bundle the browser app into ONE self-contained HTML file: CSS and every JS module
// inlined, no server needed (it opens from file:// too).
//
//   node scripts/build-single.js                 -> dist/edgy-erosion.html
//   node scripts/build-single.js out.html --fragment --no-downloads
//
// --fragment      omit <!doctype>/<html>/<head>/<body> (for hosts that add their own)
// --no-downloads  hide the PNG/snapshot download buttons (for sandboxes that block downloads)
//
// The bundler is deliberately dumb: modules are concatenated in dependency order with
// their `import` lines removed and `export` keywords stripped. That works because every
// top-level name in src/ is unique; the check below fails loudly if that stops being true.

import { readFileSync, writeFileSync, mkdirSync } from 'node:fs';
import { dirname, join } from 'node:path';
import { fileURLToPath } from 'node:url';

const root = join(dirname(fileURLToPath(import.meta.url)), '..');
const args = process.argv.slice(2);
const out = args.find((a) => !a.startsWith('--')) || join(root, 'dist', 'edgy-erosion.html');
const fragment = args.includes('--fragment');
const noDownloads = args.includes('--no-downloads');

const MODULES = [
  'src/rng.js', 'src/hexgrid.js', 'src/params.js', 'src/lithology.js', 'src/flow.js',
  'src/processes/fluvial.js', 'src/processes/hillslope.js', 'src/processes/rockfall.js', 'src/processes/coast.js',
  'src/mesh.js', 'src/simulation.js', 'src/nest.js', 'src/render.js', 'web/app.js',
];

// Every relative import must point at a bundled module, or the page would break at runtime.
for (const file of MODULES) {
  const code = readFileSync(join(root, file), 'utf8');
  for (const m of code.matchAll(/^import\s[^;]*?from\s+['"](\.[^'"]+)['"]/gms)) {
    const target = join(dirname(file), m[1]).replace(/\\/g, '/');
    if (target.endsWith('index.js')) continue; // the public API file only re-exports
    if (!MODULES.includes(target)) throw new Error(`${file} imports ${target}, which is not in MODULES`);
  }
}

const seen = new Map();
const chunks = MODULES.map((file) => {
  let code = readFileSync(join(root, file), 'utf8');
  code = code.replace(/^import\s[^;]*?from\s+['"][^'"]+['"];?\s*$/gms, '');
  if (/^\s*import\s/m.test(code)) throw new Error(`${file}: unhandled import statement`);
  if (/^export\s*\{/m.test(code)) throw new Error(`${file}: re-exports are not supported`);
  code = code.replace(/^export\s+(?=(?:async\s+)?(?:function|class|const|let|var)\b)/gm, '');
  for (const m of code.matchAll(/^(?:async\s+)?(?:function\*?|class|const|let|var)\s+([A-Za-z_$][\w$]*)/gm)) {
    if (seen.has(m[1])) throw new Error(`top-level name "${m[1]}" defined in both ${seen.get(m[1])} and ${file}`);
    seen.set(m[1], file);
  }
  if (code.includes('</script')) throw new Error(`${file}: contains "</script"`);
  return `// ---- ${file} ----\n${code}`;
});

const css = readFileSync(join(root, 'web/style.css'), 'utf8');
const html = readFileSync(join(root, 'index.html'), 'utf8');
const bodyMarkup = html.slice(html.indexOf('<body>') + 6, html.indexOf('</body>'))
  .replace(/\s*<script[^>]*src=[^>]*><\/script>/g, '');
const title = (html.match(/<title>(.*?)<\/title>/) || [, 'Edgy Erosion'])[1];

const script = `<script type="module">
${noDownloads ? 'window.EDGY_NO_DOWNLOADS = true;\n' : ''}${chunks.join('\n')}
</script>`;
const page = fragment
  ? `<title>${title}</title>\n<style>\n${css}</style>\n${bodyMarkup}\n${script}\n`
  : `<!doctype html>\n<html lang="en">\n<head>\n<meta charset="utf-8">\n<meta name="viewport" content="width=device-width, initial-scale=1">\n<title>${title}</title>\n<style>\n${css}</style>\n</head>\n<body>${bodyMarkup}\n${script}\n</body>\n</html>\n`;

mkdirSync(dirname(out), { recursive: true });
writeFileSync(out, page);
console.log(`wrote ${out} (${(page.length / 1024).toFixed(0)} KB)`);
