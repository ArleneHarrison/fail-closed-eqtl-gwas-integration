/**
 * Render Figures 1--4 from locked audit records only.
 * The figures deliberately visualise provenance, row-accounting and gate state;
 * they never calculate or display an association, posterior, biological result,
 * model comparison or repaired LD matrix.
 */
import crypto from "node:crypto";
import fs from "node:fs/promises";
import path from "node:path";
import { createRequire } from "node:module";
import { fileURLToPath } from "node:url";

const require = createRequire(import.meta.url);
const sharp = require("sharp");

// Journal-neutral high-impact bioinformatics figure system: white canvas,
// square information panels, thin rules, and a colour-blind-safe semantic palette.
const BLUE = "#0072B2";       // Technical provenance / retained process.
const VERMILION = "#D55E00";  // Fail-closed stop record.
const TEAL = "#009E73";       // Contract boundary / permitted state.
const INK = "#18324B";
const MUTED = "#4D6172";
const PANEL = "#F7F9FA";
const BORDER = "#B9C5CE";
const WIDTH = 1600;
const HEIGHT = 1120;

function requiredArg(name) {
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

function assertContains(text, expected, source) {
  if (!text.includes(expected)) throw new Error(`${source} lacks locked evidence: ${expected}`);
}

function frame(title, subtitle, body) {
  return `<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg" width="${WIDTH}" height="${HEIGHT}" viewBox="0 0 ${WIDTH} ${HEIGHT}" role="img" aria-labelledby="title description">
  <title id="title">${esc(title)}</title>
  <desc id="description">${esc(subtitle)} Workflow-only audit visualisation; no association, posterior, biological, or comparative-performance result.</desc>
  <style>
    .title { font: 700 34px Arial, Helvetica, 'DejaVu Sans', sans-serif; fill: ${INK}; letter-spacing: 0.2px; }
    .subtitle { font: 18px Arial, Helvetica, 'DejaVu Sans', sans-serif; fill: ${MUTED}; }
    .head { font: 700 22px Arial, Helvetica, 'DejaVu Sans', sans-serif; fill: ${INK}; }
    .body { font: 18px Arial, Helvetica, 'DejaVu Sans', sans-serif; fill: ${INK}; }
    .small { font: 16px Arial, Helvetica, 'DejaVu Sans', sans-serif; fill: ${MUTED}; }
    .tiny { font: 14px Arial, Helvetica, 'DejaVu Sans', sans-serif; fill: ${MUTED}; }
    .status { font: 700 19px Arial, Helvetica, 'DejaVu Sans', sans-serif; letter-spacing: 0.1px; }
  </style>
  <rect width="${WIDTH}" height="${HEIGHT}" fill="#FFFFFF"/>
  <text x="58" y="62" class="title">${esc(title)}</text>
  <text x="58" y="94" class="subtitle">${esc(subtitle)}</text>
  <line x1="58" y1="122" x2="1542" y2="122" stroke="${BORDER}" stroke-width="2"/>
  ${body}
</svg>`;
}

function box(x, y, width, height, heading, lines, accent = BLUE) {
  const lineSvg = lines.map((line, index) => `<text x="${x + 24}" y="${y + 68 + index * 25}" class="small">${esc(line)}</text>`).join("\n");
  return `<rect x="${x}" y="${y}" width="${width}" height="${height}" rx="3" fill="${PANEL}" stroke="${BORDER}" stroke-width="1.5"/>
  <rect x="${x}" y="${y}" width="${width}" height="7" rx="2" fill="${accent}"/>
  <text x="${x + 24}" y="${y + 36}" class="head">${esc(heading)}</text>${lineSvg}`;
}

function arrow(x1, y1, x2, y2, label = "") {
  const labelSvg = label ? `<text x="${(x1 + x2) / 2}" y="${(y1 + y2) / 2 - 10}" text-anchor="middle" class="tiny">${esc(label)}</text>` : "";
  return `<defs><marker id="arrow" markerWidth="10" markerHeight="10" refX="8" refY="3" orient="auto"><path d="M0,0 L0,6 L9,3 z" fill="${MUTED}"/></marker></defs>
  <line x1="${x1}" y1="${y1}" x2="${x2}" y2="${y2}" stroke="${MUTED}" stroke-width="2.5" marker-end="url(#arrow)"/>${labelSvg}`;
}

function stopMarker(x, y) {
  return `<path d="M ${x - 12} ${y - 12} L ${x + 12} ${y + 12} M ${x + 12} ${y - 12} L ${x - 12} ${y + 12}" fill="none" stroke="${VERMILION}" stroke-width="7" stroke-linecap="round"/>`;
}

function footerNote() {
  return `<line x1="58" y1="1014" x2="1542" y2="1014" stroke="${BORDER}" stroke-width="1.5"/>
  <text x="58" y="1045" class="small">Reading rule: blue/teal panels document supplied technical records; vermilion marks an INELIGIBLE stop.</text>
  <text x="58" y="1071" class="small">No panel reports an association, posterior, biological result, clinical implication, or method-performance comparison.</text>`;
}

function figure1() {
  const stages = [
    [70, "Source provenance", ["URL, version/build, access term", "retrieval date, bytes, SHA-256"]],
    [410, "Coordinate / allele audit", ["retain mapped, unmapped, ambiguous", "no automatic strand or indel repair"]],
    [750, "LD / numerical audit", ["variant order, finiteness, symmetry", "diagonal and declared stop rule"]],
    [1090, "Gate decision", ["READY: supplied contracts pass", "INELIGIBLE: preserve stop record"]],
  ];
  let body = "";
  for (const [x, title, lines] of stages) body += box(x, 220, 300, 184, title, lines, x === 1090 ? VERMILION : BLUE);
  body += arrow(370, 312, 404, 312) + arrow(710, 312, 744, 312) + arrow(1050, 312, 1084, 312);
  body += box(242, 530, 520, 210, "Gate contract", ["No imputation, LD regularisation, pruning,", "region reduction, or unrecorded allele flip", "inside the gate."], TEAL);
  body += box(838, 530, 520, 210, "Allowed downstream boundary", ["A separately specified model may be considered", "only after READY. This package records no", "model execution for either real workflow case."], VERMILION);
  body += arrow(800, 450, 800, 520, "technical status only");
  body += `<text x="58" y="160" class="head">Fail-closed architecture: provenance before model eligibility</text>${footerNote()}`;
  return frame("Figure 1. Fail-closed workflow architecture", "Source manifest to technical gate; workflow-only architecture, not an analysis result.", body);
}

function figure2() {
  const stages = [
    [70, 205, "Locked source region", ["295,093 source rows", "5,182 fixed-case rows retained"]],
    [430, 205, "Coordinate audit", ["5,047 mapped", "135 unmapped (retained in audit)"]],
    [790, 205, "Outcome scan", ["20,073,070 rows scanned", "10,587 regional rows retained"]],
    [1150, 205, "Harmonisation audit", ["4,657 reconciled", "2,863 aligned | 1,794 explicit flips", "525 failed or ambiguous retained"]],
    [250, 590, "EUR LD audit", ["503 samples", "42,414 VCF records inspected", "42,192 eligible | 222 excluded"]],
    [700, 590, "Locked LD archive", ["4,365 unique harmonized variants", "mean-impute and complete-genotype", "sensitivity archives byte-identical"]],
    [1150, 590, "Numerical stop", ["condition number 7.63 x 10^14", "predeclared stop 1 x 10^12", "no model was run"]],
  ];
  let body = `<text x="58" y="160" class="head">Fixed regional processing case: every failed or ambiguous row remains auditable</text>`;
  for (const [x, y, title, lines] of stages) body += box(x, y, 320, 190, title, lines, title === "Numerical stop" ? VERMILION : BLUE);
  body += arrow(390, 300, 424, 300) + arrow(750, 300, 784, 300) + arrow(1110, 300, 1144, 300);
  body += arrow(1310, 400, 1310, 530, "harmonized records") + arrow(1310, 530, 410, 585);
  body += arrow(570, 685, 694, 685) + arrow(1020, 685, 1144, 685);
  body += `${stopMarker(1190, 852)}<text x="1222" y="858" class="status" fill="${VERMILION}">INELIGIBLE</text>
  <text x="1222" y="886" class="status" fill="${VERMILION}">E_LD_ILL_CONDITIONED</text>
  <text x="250" y="870" class="small">Workflow behaviour only: no association, shared signal, causal, treatment, or biological inference.</text>${footerNote()}`;
  return frame("Figure 2. Fixed-case full-row audit flow", "Audited row accounting through the numerical stop; no model output is shown.", body);
}

function figure3() {
  let body = `<text x="58" y="160" class="head">Unmodified 4,365-variant EUR LD archive: technical diagnostics trigger a no-repair stop</text>`;
  body += box(70, 205, 430, 260, "Archive inputs", ["4,365 unique harmonized variants", "503 EUR samples", "no retained variant required realised", "mean imputation in the sensitivity audit"], BLUE);
  body += box(585, 205, 430, 260, "Structural diagnostics", ["finite matrix", "symmetry maximum error: 2.22 x 10^-16", "diagonal maximum error: 2.22 x 10^-16", "minimum eigenvalue: -3.69 x 10^-13"], BLUE);
  body += box(1100, 205, 430, 260, "Prospective numerical gate", ["condition number: 7.63 x 10^14", "declared stop threshold: 1 x 10^12", "threshold exceeded; no matrix repair"], VERMILION);
  body += `<rect x="250" y="560" width="1100" height="250" rx="3" fill="${PANEL}" stroke="${BORDER}" stroke-width="1.5"/>
  <rect x="250" y="560" width="1100" height="7" rx="2" fill="${TEAL}"/>
  <text x="280" y="610" class="head">Sensitivity record</text>
  <text x="280" y="654" class="body">Mean-imputed and complete-genotype-only matrices had the same SHA-256.</text>
  <text x="280" y="688" class="body">This records zero realised imputation among retained variants; it does not repair the ill-conditioned matrix.</text>
  ${stopMarker(294, 750)}<text x="325" y="757" class="status" fill="${VERMILION}">INELIGIBLE — E_LD_ILL_CONDITIONED; no SuSiE or colocalisation run</text>${footerNote()}`;
  return frame("Figure 3. Fixed-case numerical stop diagnostic", "LD integrity and conditioning record; no association or posterior plot.", body);
}

function figure4() {
  const stages = [
    [70, 200, "Coordinate control", ["GRCh38 chr11:2,000,000–2,100,000", "deterministic non-statistic selection", "251 selected rows"]],
    [440, 200, "Mapping and reconciliation", ["251 mapped", "133 aligned | 110 explicit reversals", "5 coordinate-absent | 3 unresolved"]],
    [810, 200, "Summary audit", ["243 rows", "234 unique target variants", "all source rows retained in audit"]],
    [1180, 200, "EUR VCF / LD audit", ["503 samples; 2,885 regional records", "2,871 biallelic | 14 excluded", "225 present | 9 unavailable"]],
  ];
  let body = `<text x="58" y="160" class="head">Independent coordinate-control benchmark: selected without inspecting association statistics</text>`;
  for (const [x, y, title, lines] of stages) body += box(x, y, 330, 210, title, lines, BLUE);
  body += arrow(400, 305, 434, 305) + arrow(770, 305, 804, 305) + arrow(1140, 305, 1174, 305);
  body += box(210, 575, 1180, 250, "Machine-readable stop record", [
    "E_VARIANT_DUPLICATE (9) | E_MAF_INVALID (9) | E_CASE_FRACTION_HETEROGENEOUS (243)",
    "E_LD_VARIANT_SET_MISMATCH (9) | E_LD_ILL_CONDITIONED",
    "condition number: 810317569661433.5; prospective stop threshold: 1 x 10^12; no repair or model run",
  ], VERMILION);
  body += `${stopMarker(246, 905)}<text x="278" y="912" class="status" fill="${VERMILION}">INELIGIBLE — non-biological workflow benchmark only</text>
  <text x="210" y="955" class="small">The benchmark demonstrates auditable stop behaviour; it does not estimate failure frequency.</text>
  <text x="210" y="978" class="small">It does not support an association, causal, cell-type, or clinical conclusion.</text>${footerNote()}`;
  return frame("Figure 4. Independent coordinate-control workflow benchmark", "Independently selected coordinate input through an auditable INELIGIBLE stop; no biological interpretation.", body);
}

async function writeFigure({ number, title, caption, content, sourceFiles, outDir }) {
  const slug = `figure_${number.toLowerCase().replaceAll(" ", "_").replaceAll(".", "")}`;
  const svgPath = path.join(outDir, `${slug}.svg`);
  const pngPath = path.join(outDir, `${slug}.png`);
  const manifestPath = path.join(outDir, `${slug}_manifest.json`);
  await fs.writeFile(svgPath, content, "utf8");
  await sharp(Buffer.from(content)).png().toFile(pngPath);
  const script = fileURLToPath(import.meta.url);
  const sourceHashes = Object.fromEntries(await Promise.all(sourceFiles.map(async file => [file, await sha256(file)])));
  const payload = {
    schema_version: "1.0", figure_number: number, title, caption,
    visual_style: "journal-neutral high-impact bioinformatics: white canvas, square panels, thin rules, Helvetica/Arial, and colour-blind-safe semantic accents",
    figure_contract: "workflow-only; no association, posterior, biological, clinical, causal, or comparative-performance claim",
    source_files: sourceHashes, script, script_sha256: await sha256(script),
    outputs: { png: { path: pngPath, sha256: await sha256(pngPath) }, svg: { path: svgPath, sha256: await sha256(svgPath) } },
  };
  await fs.writeFile(manifestPath, `${JSON.stringify(payload, null, 2)}\n`, "utf8");
}

async function main() {
  const fixed = requiredArg("--fixed-report");
  const coordinate = requiredArg("--coordinate-report");
  const contract = requiredArg("--contract-spec");
  const outDir = requiredArg("--out-dir");
  const [fixedText, coordinateText, contractText] = await Promise.all([fs.readFile(fixed, "utf8"), fs.readFile(coordinate, "utf8"), fs.readFile(contract, "utf8")]);
  for (const value of ["295,093 indexed source rows", "4,657 eQTL--CAD rows", "condition number `7.63e14`", "No ridge, pruning", "No coloc-SuSiE model was run."]) assertContains(fixedText, value, "fixed report");
  for (const value of ["251 retained selected-identifier rows", "243 rows entered the summary audit", "condition number 810317569661433.5", "E_VARIANT_DUPLICATE", "no coloc or SuSiE command was run"]) assertContains(coordinateText, value, "coordinate report");
  for (const value of ["E_LD_ILL_CONDITIONED", "does not invoke coloc, SuSiE, imputation, LD regularisation, LD pruning, or region reduction"]) assertContains(contractText, value, "contract specification");
  await fs.mkdir(outDir, { recursive: true });
  await writeFigure({ number: "1", title: "Fail-closed workflow architecture", caption: "Figure 1. Fail-closed workflow from immutable source record to technical gate. Workflow architecture only; a READY status is not an association or model result, and an INELIGIBLE status preserves a stop record.", content: figure1(), sourceFiles: [contract], outDir });
  await writeFigure({ number: "2", title: "Fixed-case full-row audit flow", caption: "Figure 2. Fixed-case full-row audit flow. Counts are rendered from the locked audit record and terminate at the numerical INELIGIBLE stop. No model, association, posterior, biological, or treatment result is shown.", content: figure2(), sourceFiles: [fixed], outDir });
  await writeFigure({ number: "3", title: "Fixed-case numerical stop diagnostic", caption: "Figure 3. Numerical stop diagnostic; no model execution.", content: figure3(), sourceFiles: [fixed, contract], outDir });
  await writeFigure({ number: "4", title: "Independent coordinate-control workflow benchmark", caption: "Figure 4. Independent coordinate-control workflow benchmark. The predeclared coordinate input produces an auditable INELIGIBLE stop record only and has no biological, association, causal, cell-type, or clinical interpretation.", content: figure4(), sourceFiles: [coordinate, contract], outDir });
}

await main();
