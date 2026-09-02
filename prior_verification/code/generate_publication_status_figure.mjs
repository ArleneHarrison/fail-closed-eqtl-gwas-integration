/** Render the workflow-only status figure from the locked TSV, with no analysis inputs. */
import crypto from "node:crypto";
import fs from "node:fs/promises";
import path from "node:path";
import { createRequire } from "node:module";
import { fileURLToPath } from "node:url";

const require = createRequire(import.meta.url);
const sharp = require("sharp");

const EXPECTED = [
  "synthetic_ready_positive_control",
  "synthetic_ill_conditioned_control",
  "fixed_ifitm2_cad_case",
  "independent_coordinate_control",
];
// Journal-neutral high-impact bioinformatics palette; colour is never the
// sole carrier of the READY/INELIGIBLE distinction because markers and labels remain explicit.
const READY = "#0072B2";
const STOP = "#D55E00";
const INK = "#18324B";
const MUTED = "#4D6172";

function arg(name) {
  const index = process.argv.indexOf(name);
  if (index < 0 || !process.argv[index + 1]) throw new Error(`missing ${name}`);
  return path.resolve(process.argv[index + 1]);
}

function esc(value) {
  return String(value).replaceAll("&", "&amp;").replaceAll("<", "&lt;").replaceAll(">", "&gt;");
}

async function sha256(file) {
  return crypto.createHash("sha256").update(await fs.readFile(file)).digest("hex");
}

function parseTsv(text) {
  const lines = text.trim().split(/\r?\n/);
  const header = lines.shift().split("\t");
  return lines.map(line => Object.fromEntries(header.map((key, i) => [key, line.split("\t")[i] ?? ""])));
}

function validate(rows) {
  const required = ["case_id", "case_class", "summary_rows", "unique_summary_variants", "ld_variant_ids", "status", "status_codes", "condition_number", "interpretation_boundary", "source_artifact"];
  if (rows.length !== 4 || rows.some(row => Object.keys(row).length !== required.length || required.some(key => !(key in row)))) {
    throw new Error("publication table schema differs from the locked figure contract");
  }
  if (JSON.stringify(rows.map(row => row.case_id)) !== JSON.stringify(EXPECTED)) throw new Error("case order differs from locked figure contract");
  if (JSON.stringify(rows.map(row => row.status)) !== JSON.stringify(["READY", "INELIGIBLE", "INELIGIBLE", "INELIGIBLE"])) throw new Error("unexpected gate status");
}

function marker(x, y, ready) {
  if (ready) return `<circle cx="${x}" cy="${y}" r="14" fill="${READY}" stroke="#FFFFFF" stroke-width="2"/>`;
  return `<path d="M ${x - 11} ${y - 11} L ${x + 11} ${y + 11} M ${x + 11} ${y - 11} L ${x - 11} ${y + 11}" fill="none" stroke="${STOP}" stroke-width="7" stroke-linecap="round"/>`;
}

function rowSvg(row, y, label) {
  const ready = row.status === "READY";
  const color = ready ? READY : STOP;
  return `
    ${marker(70, y, ready)}
    <text x="105" y="${y - 8}" class="label">${esc(label)}</text>
    <text x="105" y="${y + 18}" class="status" fill="${color}">${esc(row.status)}</text>
    <text x="620" y="${y - 5}" class="detail">summary rows: ${esc(row.summary_rows)} | LD IDs: ${esc(row.ld_variant_ids)}</text>
    <text x="620" y="${y + 20}" class="detail">condition number: ${esc(row.condition_number)}</text>`;
}

function svg(rows) {
  return `<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg" width="1400" height="900" viewBox="0 0 1400 900" role="img" aria-labelledby="title description">
  <title id="title">Executable gate outcomes</title>
  <desc id="description">Workflow-only status display. It shows two synthetic contract controls and two real workflow cases. It contains no association, posterior, or biological interpretation.</desc>
  <style>
    .title { font: 700 32px Arial, Helvetica, 'DejaVu Sans', sans-serif; fill: ${INK}; letter-spacing: 0.2px; }
    .subtitle { font: 18px Arial, Helvetica, 'DejaVu Sans', sans-serif; fill: ${MUTED}; }
    .section { font: 700 23px Arial, Helvetica, 'DejaVu Sans', sans-serif; fill: ${INK}; }
    .label { font: 700 20px Arial, Helvetica, 'DejaVu Sans', sans-serif; fill: ${INK}; }
    .status { font: 700 18px Arial, Helvetica, 'DejaVu Sans', sans-serif; letter-spacing: 0.1px; }
    .detail { font: 17px Arial, Helvetica, 'DejaVu Sans', sans-serif; fill: ${MUTED}; }
    .note { font: 16px Arial, Helvetica, 'DejaVu Sans', sans-serif; fill: ${MUTED}; }
    .legend { font: 16px Arial, Helvetica, 'DejaVu Sans', sans-serif; fill: ${INK}; }
  </style>
  <rect width="1400" height="900" fill="#FFFFFF"/>
  <text x="48" y="62" class="title">Executable gate outcomes</text>
  <text x="48" y="94" class="subtitle">Workflow-only status display from the locked publication table; no association, posterior, or biological interpretation.</text>
  <line x1="48" y1="124" x2="1352" y2="124" stroke="#D6DEE5" stroke-width="2"/>
  <rect x="48" y="156" width="1304" height="228" rx="3" fill="#F7F9FA" stroke="#B9C5CE" stroke-width="1.5"/>
  <rect x="48" y="156" width="1304" height="7" rx="2" fill="${READY}"/>
  <text x="72" y="192" class="section">A. Synthetic contract controls</text>
  ${rowSvg(rows[0], 255, "Complete-contract positive control")}
  ${rowSvg(rows[1], 330, "Numerical-stop negative control")}
  <rect x="48" y="414" width="1304" height="228" rx="3" fill="#F7F9FA" stroke="#B9C5CE" stroke-width="1.5"/>
  <rect x="48" y="414" width="1304" height="7" rx="2" fill="${STOP}"/>
  <text x="72" y="450" class="section">B. Real workflow cases</text>
  ${rowSvg(rows[2], 513, "Fixed regional processing case")}
  ${rowSvg(rows[3], 588, "Independent coordinate-control benchmark")}
  <line x1="48" y1="676" x2="1352" y2="676" stroke="#B9C5CE" stroke-width="1.5"/>
  <text x="72" y="715" class="note">Reading rule: READY means the supplied technical contracts passed. INELIGIBLE preserves a stop record.</text>
  <text x="72" y="744" class="note">Neither status is a biological result, a model result, a posterior, or a comparative-performance result.</text>
  ${marker(72, 793, true)}<text x="98" y="799" class="legend">READY</text>
  ${marker(245, 793, false)}<text x="271" y="799" class="legend">INELIGIBLE</text>
</svg>`;
}

async function main() {
  const source = arg("--source-tsv");
  const outSvg = arg("--out-svg");
  const outPng = arg("--out-png");
  const manifest = arg("--manifest");
  const rows = parseTsv(await fs.readFile(source, "utf8"));
  validate(rows);
  const content = svg(rows);
  await fs.mkdir(path.dirname(outSvg), { recursive: true });
  await fs.mkdir(path.dirname(outPng), { recursive: true });
  await fs.mkdir(path.dirname(manifest), { recursive: true });
  await fs.writeFile(outSvg, content, "utf8");
  await sharp(Buffer.from(content)).png().toFile(outPng);
  const script = fileURLToPath(import.meta.url);
  const payload = {
    schema_version: "1.0", source_tsv: source, source_tsv_sha256: await sha256(source),
    script, script_sha256: await sha256(script),
    outputs: { png: { path: outPng, sha256: await sha256(outPng) }, svg: { path: outSvg, sha256: await sha256(outSvg) } },
    figure_contract: "workflow-only; no association, posterior, biological, or comparative-performance claim",
    visual_style: "journal-neutral high-impact bioinformatics: white canvas, square panels, thin rules, Helvetica/Arial, and colour-blind-safe semantic accents",
    case_ids: rows.map(row => row.case_id),
    status_codes: Object.fromEntries(rows.map(row => [row.case_id, row.status_codes])),
  };
  await fs.writeFile(manifest, `${JSON.stringify(payload, null, 2)}\n`, "utf8");
}

await main();
