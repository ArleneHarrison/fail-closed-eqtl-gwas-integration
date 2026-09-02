import fs from "node:fs";
import path from "node:path";
import crypto from "node:crypto";
import { createRequire } from "node:module";

const require = createRequire(import.meta.url);
const sharp = require("sharp");

const BLUE = "#0072B2";
const TEAL = "#009E73";
const VERMILION = "#D55E00";
const INK = "#18324B";
const MUTED = "#4D6172";
const PANEL = "#F7F9FA";
const BORDER = "#B9C5CE";
const GRID = "#E3E9ED";

function args() {
  const result = {};
  for (let i = 2; i < process.argv.length; i += 2) result[process.argv[i].replace(/^--/, "")] = process.argv[i + 1];
  return result;
}
function esc(value) {
  return String(value).replaceAll("&", "&amp;").replaceAll("<", "&lt;").replaceAll(">", "&gt;");
}
function sha(file) {
  return crypto.createHash("sha256").update(fs.readFileSync(file)).digest("hex");
}
function parseTsv(file) {
  const lines = fs.readFileSync(file, "utf8").trim().split(/\r?\n/);
  const header = lines[0].split("\t");
  return lines.slice(1).map(line => Object.fromEntries(line.split("\t").map((value, index) => [header[index], value])));
}
function text(x, y, value, size = 18, weight = 400, fill = INK, anchor = "start") {
  return `<text x="${x}" y="${y}" font-size="${size}" font-weight="${weight}" fill="${fill}" text-anchor="${anchor}">${esc(value)}</text>`;
}
function tile(x, y, width, height, title, detail, accent) {
  return `<g><rect x="${x}" y="${y}" width="${width}" height="${height}" rx="3" fill="${PANEL}" stroke="${BORDER}" stroke-width="1.4"/>` +
    `<rect x="${x}" y="${y}" width="${width}" height="6" rx="2" fill="${accent}"/>` +
    text(x + 16, y + 31, title, 16, 700) + text(x + 16, y + 57, detail, 11.5, 600, accent) + `</g>`;
}
function logMap(value, min, max, left, right) {
  const fraction = (Math.log10(value) - Math.log10(min)) / (Math.log10(max) - Math.log10(min));
  return left + fraction * (right - left);
}

const a = args();
for (const key of ["mutation-summary", "boundary", "scaling-summary", "out-png", "out-svg", "manifest"]) {
  if (!a[key]) throw new Error(`missing --${key}`);
}
const mutation = parseTsv(a["mutation-summary"]);
const boundary = parseTsv(a.boundary);
const scaling = parseTsv(a["scaling-summary"]);
const labels = {
  valid_complete_contract: "Complete contract",
  provenance_incomplete: "Provenance incomplete",
  variant_id_missing: "Variant ID missing",
  variant_duplicate: "Variant duplicate",
  eqtl_se_invalid: "eQTL SE invalid",
  gwas_se_invalid: "GWAS SE invalid",
  maf_invalid: "MAF invalid",
  eqtl_n_invalid: "eQTL N invalid",
  gwas_n_invalid: "GWAS N invalid",
  allele_unresolved: "Allele unresolved",
  vcf_absent: "VCF absent",
  case_fraction_heterogeneous: "Case fraction heterogeneous",
  ld_variant_ids_missing: "LD IDs missing",
  ld_variant_set_mismatch: "LD set mismatch",
  eqtl_trait_invalid: "eQTL trait invalid",
  gwas_trait_invalid: "GWAS trait invalid",
  case_fraction_invalid: "Case fraction invalid",
  ld_asymmetric: "LD asymmetric",
  ld_diagonal_invalid: "LD diagonal invalid",
  ld_nonfinite: "LD non-finite",
  ld_ill_conditioned: "LD ill-conditioned",
  ld_threshold_unspecified: "LD threshold unspecified",
};
let body = `<svg xmlns="http://www.w3.org/2000/svg" width="1400" height="1350" viewBox="0 0 1400 1350">` +
  `<rect width="1400" height="1350" fill="white"/>` +
  `<style>text{font-family:Arial,Helvetica,'DejaVu Sans',sans-serif;letter-spacing:.05px}</style>` +
  text(70, 66, "Executable validation of fail-closed preflight contracts", 30, 700) +
  text(70, 100, "Synthetic workflow validation; no association, posterior, biological, or comparative-performance result.", 16, 400, MUTED) +
  `<line x1="70" y1="124" x2="1330" y2="124" stroke="${BORDER}" stroke-width="1.4"/>` +
  text(70, 164, "a", 25, 700) + text(106, 164, "Contract-mutation coverage", 22, 700) +
  text(106, 190, "100 replicates per case; each tile reports the expected READY decision or failure code.", 15, 400, MUTED);

mutation.forEach((row, index) => {
  const column = index % 4;
  const r = Math.floor(index / 4);
  const x = 70 + column * 320;
  const y = 215 + r * 91;
  const accent = row.case === "valid_complete_contract" ? BLUE : TEAL;
  const outcome = row.case === "valid_complete_contract" ? "100/100 READY" : `100/100 ${row.expected_code}`;
  body += tile(x, y, 292, 73, labels[row.case] ?? row.case, outcome, accent);
});

const lowerTop = 215 + Math.ceil(mutation.length / 4) * 91 + 28;
body += `<line x1="70" y1="${lowerTop}" x2="1330" y2="${lowerTop}" stroke="${BORDER}" stroke-width="1.2"/>`;

// Panel b: normalized boundary behavior.
const bx0 = 70, bx1 = 665, by0 = lowerTop + 40;
body += text(bx0, by0, "b", 25, 700) + text(bx0 + 36, by0, "Boundary behavior", 22, 700);
body += text(bx0 + 36, by0 + 27, "Diagnostic values are normalized to the prospectively declared threshold.", 14, 400, MUTED);
const plotLeft = bx0 + 180, plotRight = bx1 - 20, axisTop = by0 + 80, axisBottom = by0 + 305;
const metrics = ["Asymmetry", "Diagonal error", "Condition number"];
metrics.forEach((metric, index) => {
  const y = axisTop + index * 92;
  body += text(plotLeft - 20, y + 6, metric, 14, 600, INK, "end");
  body += `<line x1="${plotLeft}" y1="${y}" x2="${plotRight}" y2="${y}" stroke="${GRID}" stroke-width="2"/>`;
});
const stopX = logMap(1, 0.35, 4.8, plotLeft, plotRight);
body += `<line x1="${stopX}" y1="${axisTop - 32}" x2="${stopX}" y2="${axisBottom}" stroke="${INK}" stroke-width="1.5" stroke-dasharray="6 5"/>`;
body += text(stopX + 8, axisTop - 40, "stop boundary", 13, 400, MUTED);
for (const row of boundary) {
  let metric, value, threshold;
  if (row.case.startsWith("asymmetry")) { metric = "Asymmetry"; value = Number(row.symmetry_max_error); threshold = 1e-8; }
  else if (row.case.startsWith("diagonal")) { metric = "Diagonal error"; value = Number(row.diagonal_max_error); threshold = 1e-8; }
  else { metric = "Condition number"; value = Number(row.condition_number); threshold = 1e12; }
  const x = logMap(value / threshold, 0.35, 4.8, plotLeft, plotRight);
  const y = axisTop + metrics.indexOf(metric) * 92;
  if (row.observed_status === "READY") body += `<circle cx="${x}" cy="${y}" r="10" fill="${BLUE}"/>`;
  else body += `<g stroke="${VERMILION}" stroke-width="7" stroke-linecap="round"><line x1="${x - 9}" y1="${y - 9}" x2="${x + 9}" y2="${y + 9}"/><line x1="${x - 9}" y1="${y + 9}" x2="${x + 9}" y2="${y - 9}"/></g>`;
}
[0.5, 1, 2, 4].forEach(value => {
  const x = logMap(value, 0.35, 4.8, plotLeft, plotRight);
  body += text(x, axisBottom + 26, `${value}×`, 12, 400, MUTED, "middle");
});
body += text((plotLeft + plotRight) / 2, axisBottom + 58, "Observed diagnostic / declared threshold", 14, 600, INK, "middle");
body += `<circle cx="${plotLeft}" cy="${axisBottom + 93}" r="8" fill="${BLUE}"/>${text(plotLeft + 18, axisBottom + 99, "READY", 13, 600, BLUE)}`;
body += `<g stroke="${VERMILION}" stroke-width="5" stroke-linecap="round"><line x1="${plotLeft + 100}" y1="${axisBottom + 85}" x2="${plotLeft + 116}" y2="${axisBottom + 101}"/><line x1="${plotLeft + 100}" y1="${axisBottom + 101}" x2="${plotLeft + 116}" y2="${axisBottom + 85}"/></g>${text(plotLeft + 126, axisBottom + 99, "INELIGIBLE", 13, 600, VERMILION)}`;

// Panel c: log-log runtime scaling.
const cx0 = 735, cx1 = 1330;
body += text(cx0, by0, "c", 25, 700) + text(cx0 + 36, by0, "Computational scaling", 22, 700);
body += text(cx0 + 36, by0 + 27, "Five complete-contract runs at each variant count; all emitted READY.", 14, 400, MUTED);
const cLeft = cx0 + 95, cRight = cx1 - 25, cTop = axisTop - 10, cBottom = axisBottom;
const xMin = 10, xMax = 1000, yMin = 0.05, yMax = 250;
for (const value of [0.1, 1, 10, 100]) {
  const y = cBottom - logMap(value, yMin, yMax, 0, cBottom - cTop);
  body += `<line x1="${cLeft}" y1="${y}" x2="${cRight}" y2="${y}" stroke="${GRID}" stroke-width="1.5"/>`;
  body += text(cLeft - 14, y + 5, String(value), 12, 400, MUTED, "end");
}
for (const value of [10, 50, 100, 250, 500, 1000]) {
  const x = logMap(value, xMin, xMax, cLeft, cRight);
  body += text(x, cBottom + 26, String(value), 12, 400, MUTED, "middle");
}
body += `<line x1="${cLeft}" y1="${cTop}" x2="${cLeft}" y2="${cBottom}" stroke="${BORDER}" stroke-width="1.4"/><line x1="${cLeft}" y1="${cBottom}" x2="${cRight}" y2="${cBottom}" stroke="${BORDER}" stroke-width="1.4"/>`;
const points = scaling.map(row => ({
  x: logMap(Number(row.n_variants), xMin, xMax, cLeft, cRight),
  y: cBottom - logMap(Number(row.median_runtime_ms), yMin, yMax, 0, cBottom - cTop),
  yMin: cBottom - logMap(Number(row.min_runtime_ms), yMin, yMax, 0, cBottom - cTop),
  yMax: cBottom - logMap(Number(row.max_runtime_ms), yMin, yMax, 0, cBottom - cTop),
}));
body += `<polyline points="${points.map(p => `${p.x},${p.y}`).join(" ")}" fill="none" stroke="${BLUE}" stroke-width="4"/>`;
for (const p of points) {
  body += `<line x1="${p.x}" y1="${p.yMin}" x2="${p.x}" y2="${p.yMax}" stroke="${BLUE}" stroke-opacity=".35" stroke-width="8"/><circle cx="${p.x}" cy="${p.y}" r="7" fill="${BLUE}"/>`;
}
body += text((cLeft + cRight) / 2, cBottom + 58, "Variants in complete contract", 14, 600, INK, "middle");
body += `<text x="${cx0 + 24}" y="${(cTop + cBottom) / 2}" font-size="14" font-weight="600" fill="${INK}" text-anchor="middle" transform="rotate(-90 ${cx0 + 24} ${(cTop + cBottom) / 2})">Gate runtime (ms)</text>`;

body += `<line x1="70" y1="1265" x2="1330" y2="1265" stroke="${BORDER}" stroke-width="1.3"/>`;
body += text(70, 1300, "Reading rule: colour supplements text and marker shape; 95% exact binomial intervals for 100/100 equal 96.4%–100%.", 13, 400, MUTED);
body += text(70, 1327, "READY is a technical eligibility decision only. The benchmark does not evaluate statistical power, biology, or superiority.", 13, 400, MUTED);
body += `</svg>`;

for (const file of [a["out-png"], a["out-svg"], a.manifest]) fs.mkdirSync(path.dirname(file), {recursive: true});
fs.writeFileSync(a["out-svg"], body, "utf8");
await sharp(Buffer.from(body)).png().toFile(a["out-png"]);
const sources = [a["mutation-summary"], a.boundary, a["scaling-summary"]].map(file => ({path: file, sha256: sha(file), byte_count: fs.statSync(file).size}));
const manifest = {
  figure_number: "5",
  title: "Executable validation of fail-closed preflight contracts",
  workflow_only: true,
  figure_contract: "workflow-only; no association, posterior, biological, clinical, causal, or comparative-performance inference",
  visual_style: "journal-neutral high-impact bioinformatics: white canvas, square panels, thin rules, Arial/Helvetica-compatible typography, and colour-blind-safe semantic accents",
  sources,
  outputs: {
    png: {path: a["out-png"], sha256: sha(a["out-png"]), byte_count: fs.statSync(a["out-png"]).size},
    svg: {path: a["out-svg"], sha256: sha(a["out-svg"]), byte_count: fs.statSync(a["out-svg"]).size},
  },
  caption: "Figure 5. Executable contract validation and scaling. Synthetic workflow benchmark only; no biological or model result.",
};
fs.writeFileSync(a.manifest, `${JSON.stringify(manifest, null, 2)}\n`, "utf8");
process.stdout.write(`${JSON.stringify(manifest.outputs)}\n`);
