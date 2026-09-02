import fs from "node:fs";
import path from "node:path";
import crypto from "node:crypto";
import { createRequire } from "node:module";

const require = createRequire(import.meta.url);
const sharp = require("sharp");

const BLUE = "#0072B2";
const SKY = "#56B4E9";
const TEAL = "#009E73";
const ORANGE = "#E69F00";
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
function logY(value, min, max, top, bottom) {
  const f = (Math.log10(value) - Math.log10(min)) / (Math.log10(max) - Math.log10(min));
  return bottom - f * (bottom - top);
}

const a = args();
for (const key of ["summary-tsv", "summary-json", "input-verification", "out-png", "out-svg", "manifest"]) {
  if (!a[key]) throw new Error(`missing --${key}`);
}
const rows = parseTsv(a["summary-tsv"]);
const aggregate = JSON.parse(fs.readFileSync(a["summary-json"], "utf8"));
const codes = [
  "E_NO_ELIGIBLE_GENE",
  "E_VARIANT_DUPLICATE",
  "E_MAF_INVALID",
  "E_CASE_FRACTION_HETEROGENEOUS",
  "E_LD_VARIANT_SET_MISMATCH",
  "E_LD_ILL_CONDITIONED",
];
const shortCodes = {
  E_NO_ELIGIBLE_GENE: "No eligible identifier",
  E_VARIANT_DUPLICATE: "Duplicate variant",
  E_MAF_INVALID: "MAF invalid",
  E_CASE_FRACTION_HETEROGENEOUS: "Case fraction heterogeneous",
  E_LD_VARIANT_SET_MISMATCH: "LD variant-set mismatch",
  E_LD_ILL_CONDITIONED: "LD ill-conditioned",
};

let body = `<svg xmlns="http://www.w3.org/2000/svg" width="1400" height="1070" viewBox="0 0 1400 1070">` +
  `<rect width="1400" height="1070" fill="white"/>` +
  `<style>text{font-family:Arial,Helvetica,'DejaVu Sans',sans-serif;letter-spacing:.05px}</style>` +
  text(70, 62, "Predeclared real-data transfer benchmark", 30, 700) +
  text(70, 96, "Twelve deterministic chr11 windows; workflow eligibility only, with no model or biological interpretation.", 16, 400, MUTED) +
  `<line x1="70" y1="120" x2="1330" y2="120" stroke="${BORDER}" stroke-width="1.4"/>`;

// Panel a: predeclared sampling and retained region outcomes.
body += text(70, 160, "a", 25, 700) + text(106, 160, "Deterministic region transfer", 22, 700);
body += text(106, 186, "500-kb windows were fixed every 10 Mb before regional statistics or gate outcomes were inspected.", 14, 400, MUTED);
const x0 = 110, x1 = 1310, axisY = 245;
body += `<line x1="${x0}" y1="${axisY}" x2="${x1}" y2="${axisY}" stroke="${INK}" stroke-width="3"/>`;
rows.forEach((row, index) => {
  const x = x0 + index * ((x1 - x0) / (rows.length - 1));
  const noGene = row.status_codes.includes("E_NO_ELIGIBLE_GENE");
  const color = noGene ? ORANGE : VERMILION;
  body += `<line x1="${x}" y1="${axisY - 12}" x2="${x}" y2="${axisY + 12}" stroke="${INK}" stroke-width="1.2"/>`;
  body += `<circle cx="${x}" cy="${axisY}" r="11" fill="${color}" stroke="white" stroke-width="2"/>`;
  body += text(x, axisY + 36, row.region_id, 13, 700, INK, "middle");
  body += text(x, axisY + 56, `${Math.round(Number(row.source_region.split(":")[1].split("-")[0]) / 1e6)} Mb`, 11, 400, MUTED, "middle");
});
body += `<rect x="110" y="322" width="1200" height="56" rx="3" fill="${PANEL}" stroke="${BORDER}"/>`;
body += text(135, 347, "12/12 retained", 15, 700) + text(282, 347, "11 reached LD preflight", 15, 700) + text(514, 347, "1 retained as E_NO_ELIGIBLE_GENE", 15, 700, ORANGE);
body += text(887, 347, "0 READY; 12 INELIGIBLE", 15, 700, VERMILION);
body += text(135, 368, "Every declared window has a machine-readable terminal status; no replacement window was selected.", 13, 400, MUTED);

// Panel b: region-by-code heatmap.
body += `<line x1="70" y1="410" x2="1330" y2="410" stroke="${BORDER}" stroke-width="1.2"/>`;
body += text(70, 450, "b", 25, 700) + text(106, 450, "Observed status-code matrix", 22, 700);
body += text(106, 476, "Filled cells indicate that the preflight record retained the named technical status.", 14, 400, MUTED);
const heatLeft = 400, heatTop = 505, cellW = 66, cellH = 46;
rows.forEach((row, i) => body += text(heatLeft + i * cellW + cellW / 2, heatTop - 14, row.region_id, 12, 700, INK, "middle"));
codes.forEach((code, r) => {
  body += text(heatLeft - 18, heatTop + r * cellH + 28, shortCodes[code], 13, 600, INK, "end");
  rows.forEach((row, c) => {
    const hit = row.status_codes.split(";").includes(code);
    body += `<rect x="${heatLeft + c * cellW}" y="${heatTop + r * cellH}" width="${cellW - 5}" height="${cellH - 5}" rx="2" fill="${hit ? (code === "E_NO_ELIGIBLE_GENE" ? ORANGE : TEAL) : PANEL}" stroke="${hit ? "white" : GRID}"/>`;
    if (hit) body += text(heatLeft + c * cellW + (cellW - 5) / 2, heatTop + r * cellH + 28, "•", 24, 700, "white", "middle");
  });
  const observed = aggregate.status_code_distribution[code]?.regions ?? 0;
  body += text(heatLeft + rows.length * cellW + 24, heatTop + r * cellH + 28, `${observed}/12`, 13, 700, code === "E_NO_ELIGIBLE_GENE" ? ORANGE : TEAL);
});

// Panel c: numerical stop values for regions reaching LD.
body += `<line x1="70" y1="810" x2="1330" y2="810" stroke="${BORDER}" stroke-width="1.2"/>`;
body += text(70, 850, "c", 25, 700) + text(106, 850, "LD conditioning stop", 22, 700);
body += text(106, 876, "Condition number for the 11 regions reaching LD preflight; the supplied procedural threshold was 10¹².", 14, 400, MUTED);
const numRows = rows.filter(row => row.ld_condition_number);
const plotLeft = 180, plotRight = 1305, plotTop = 902, plotBottom = 988;
for (const tick of [1e12, 1e13, 1e14, 1e15]) {
  const y = logY(tick, 5e11, 4e15, plotTop, plotBottom);
  body += `<line x1="${plotLeft}" y1="${y}" x2="${plotRight}" y2="${y}" stroke="${tick === 1e12 ? VERMILION : GRID}" stroke-width="${tick === 1e12 ? 2 : 1.2}" stroke-dasharray="${tick === 1e12 ? "6 5" : "none"}"/>`;
  body += text(plotLeft - 14, y + 4, `10${tick === 1e12 ? "¹²" : tick === 1e13 ? "¹³" : tick === 1e14 ? "¹⁴" : "¹⁵"}`, 11, 400, MUTED, "end");
}
numRows.forEach((row, i) => {
  const x = plotLeft + 70 + i * ((plotRight - plotLeft - 140) / (numRows.length - 1));
  const y = logY(Number(row.ld_condition_number), 5e11, 4e15, plotTop, plotBottom);
  body += `<circle cx="${x}" cy="${y}" r="7" fill="${BLUE}"/>`;
  body += text(x, plotBottom + 28, row.region_id, 11, 700, INK, "middle");
});
body += text(plotRight, plotTop - 9, "No repair, pruning, ridge adjustment, or model call", 12.5, 600, VERMILION, "end");

body += `<line x1="70" y1="1025" x2="1330" y2="1025" stroke="${BORDER}" stroke-width="1.2"/>`;
body += text(70, 1050, "Observed frequencies are restricted to this locked QTD000216–GCST90132314–1000G-EUR chr11 benchmark.", 12.5, 400, MUTED);
body += `</svg>`;

for (const file of [a["out-png"], a["out-svg"], a.manifest]) fs.mkdirSync(path.dirname(file), {recursive: true});
fs.writeFileSync(a["out-svg"], body, "utf8");
await sharp(Buffer.from(body)).png().toFile(a["out-png"]);
const sources = [a["summary-tsv"], a["summary-json"], a["input-verification"]].map(file => ({path: file, sha256: sha(file), byte_count: fs.statSync(file).size}));
const manifest = {
  figure_number: "6",
  title: "Predeclared real-data transfer benchmark",
  figure_contract: "workflow-only; no association, posterior, biological, clinical, causal, or comparative-performance inference",
  visual_style: "journal-neutral high-impact bioinformatics: white canvas, compact aligned panels, thin rules, Arial/Helvetica-compatible typography, and colour-blind-safe semantic accents",
  sources,
  outputs: {
    png: {path: a["out-png"], sha256: sha(a["out-png"]), byte_count: fs.statSync(a["out-png"]).size},
    svg: {path: a["out-svg"], sha256: sha(a["out-svg"]), byte_count: fs.statSync(a["out-svg"]).size},
  },
  caption: "Figure 6. Predeclared 12-region real-data transfer benchmark. All windows were retained with a machine-readable terminal status. Eleven windows reached LD preflight and were ineligible; one contained no identifier satisfying the locked completeness rule. Status-code frequencies and condition numbers describe workflow behavior only within this fixed chr11 benchmark. No model was run and no association, posterior, biological, causal, clinical, or method-superiority inference is supported.",
};
fs.writeFileSync(a.manifest, `${JSON.stringify(manifest, null, 2)}\n`, "utf8");
process.stdout.write(`${JSON.stringify(manifest.outputs)}\n`);
