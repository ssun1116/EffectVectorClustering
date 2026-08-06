import base64
import html
import json
import textwrap

import numpy as np


def load_plotly_heatmap_payload(path):
    text = path.read_text()
    call_start = text.rfind("Plotly.newPlot(")
    if call_start < 0:
        raise ValueError(f"Could not find Plotly.newPlot in {path}")
    data_start = text.find("[", call_start)
    data, _ = json.JSONDecoder().raw_decode(text[data_start:])
    heatmap = next(trace for trace in data if trace.get("type") == "heatmap")
    z = heatmap["z"]
    if not (isinstance(z, dict) and "bdata" in z):
        matrix = np.asarray(z, dtype=np.float64)
        z = {
            "dtype": "f8",
            "bdata": base64.b64encode(matrix.tobytes()).decode("ascii"),
            "shape": ",".join(map(str, matrix.shape)),
        }
    return {
        "x": heatmap["x"],
        "y": heatmap["y"],
        "z": z,
        "zmin": heatmap.get("zmin", -1),
        "zmax": heatmap.get("zmax", 1),
    }


def as_float(row, key):
    val = row.get(key, "")
    return None if val == "" else float(val)


def build_modules(summaries, members):
    by_module = {}
    for row in members:
        by_module.setdefault(str(row["module_id"]), []).append(row)

    modules = []
    for row in summaries:
        module_key = str(row["module_id"])
        module_id = int(row["module_id"])
        module_members = by_module.get(module_key, [])
        modules.append(
            {
                "module_id": module_id,
                "start": int(row["start"]),
                "end": int(row["end"]),
                "size": int(row["size"]),
                "within_mean": as_float(row, "within_mean"),
                "neighbor_mean": as_float(row, "neighbor_mean"),
                "contrast": as_float(row, "contrast"),
                "n_unique_genes": int(row["n_unique_genes"]),
                "top_gene_fraction": as_float(row, "top_gene_fraction"),
                "gene_entropy": as_float(row, "gene_entropy"),
                "n_unique_cell_types": int(row["n_unique_cell_types"]),
                "top_genes": row["top_genes"],
                "top_cells": row["top_cells"],
                "auto_theme": row.get("auto_theme", ""),
                "theme_hits": row.get("theme_hits", ""),
                "members": [
                    {
                        "position": int(m["position"]),
                        "label": m["label"],
                        "cell_type": m["cell_type"],
                        "gene": m["gene"],
                    }
                    for m in module_members
                ],
            }
        )
    return modules


def write_browser(clustered_html, summaries, members, output_html):
    payload = load_plotly_heatmap_payload(clustered_html)
    modules = build_modules(summaries, members)
    data = {
        "labels": payload["x"],
        "yLabels": payload["y"],
        "z": payload["z"],
        "zmin": payload["zmin"],
        "zmax": payload["zmax"],
        "modules": modules,
    }
    data_json = json.dumps(data, separators=(",", ":")).replace("</", "<\\/")
    html_text = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Effect Vector Heatmap Module Browser</title>
<style>
  :root {{
    color-scheme: light;
    --ink: #17202a;
    --muted: #667085;
    --line: #d0d5dd;
    --panel: #f8fafc;
    --accent: #b42318;
    --blue: #275f9f;
  }}
  * {{ box-sizing: border-box; }}
  body {{
    margin: 0;
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
    color: var(--ink);
    background: white;
  }}
  header {{
    padding: 16px 20px 10px;
    border-bottom: 1px solid var(--line);
  }}
  h1 {{
    margin: 0 0 6px;
    font-size: 20px;
    font-weight: 650;
  }}
  .sub {{
    margin: 0;
    color: var(--muted);
    font-size: 13px;
    line-height: 1.35;
  }}
  main {{
    display: grid;
    grid-template-columns: minmax(0, 1fr) minmax(0, 1fr);
    gap: 14px;
    padding: 14px;
    min-height: calc(100vh - 74px);
  }}
  .heatmapPanel, .sidePanel {{
    min-width: 0;
  }}
  .toolbar {{
    display: flex;
    gap: 8px;
    align-items: center;
    flex-wrap: wrap;
    margin-bottom: 10px;
    width: min(100%, calc(100vh - 126px));
  }}
  button, select, input {{
    border: 1px solid var(--line);
    background: white;
    color: var(--ink);
    border-radius: 6px;
    padding: 6px 8px;
    font-size: 13px;
  }}
  button {{
    cursor: pointer;
  }}
  #searchBox {{
    flex: 1 1 240px;
    min-width: 220px;
  }}
  button.active {{
    border-color: var(--blue);
    color: var(--blue);
    font-weight: 600;
  }}
  .canvasWrap {{
    position: relative;
    width: min(100%, calc(100vh - 126px));
    aspect-ratio: 1 / 1;
    border: 1px solid var(--line);
    background: #fff;
  }}
  .heatmapArea {{
    display: flex;
    align-items: flex-start;
    gap: 10px;
  }}
  .heatmapArea.zoomed {{
    display: grid;
    grid-template-columns: 230px minmax(0, 1fr) 40px;
    grid-template-rows: 150px auto;
    gap: 8px 10px;
    align-items: stretch;
  }}
  .heatmapArea.zoomed .canvasWrap {{
    grid-column: 2;
    grid-row: 2;
    width: min(100%, calc(100vh - 270px));
  }}
  .heatmapArea.zoomed .legend {{
    grid-column: 3;
    grid-row: 2;
  }}
  .xAxisLabels, .yAxisLabels {{
    display: none;
    color: var(--ink);
    font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace;
    font-size: 9px;
    line-height: 1.1;
  }}
  .heatmapArea.zoomed .xAxisLabels {{
    display: block;
    grid-column: 2;
    grid-row: 1;
    position: relative;
    overflow: hidden;
  }}
  .heatmapArea.zoomed .yAxisLabels {{
    display: block;
    grid-column: 1;
    grid-row: 2;
    position: relative;
    overflow: hidden;
  }}
  .xTick {{
    position: absolute;
    bottom: 0;
    width: 150px;
    transform: translateX(-2px) rotate(-58deg);
    transform-origin: bottom left;
    white-space: nowrap;
  }}
  .yTick {{
    position: absolute;
    right: 0;
    width: 225px;
    transform: translateY(-50%);
    white-space: nowrap;
    text-align: right;
    padding-right: 4px;
  }}
  canvas {{
    position: absolute;
    inset: 0;
    width: 100%;
    height: 100%;
    image-rendering: pixelated;
  }}
  .legend {{
    display: flex;
    flex-direction: column;
    align-items: center;
    gap: 6px;
    min-width: 40px;
    color: var(--muted);
    font-size: 12px;
  }}
  .ramp {{
    width: 12px;
    height: 160px;
    border: 1px solid var(--line);
    background: linear-gradient(180deg, rgb(185,91,90), rgb(239,237,240), rgb(87,129,188));
  }}
  .status {{
    writing-mode: vertical-rl;
    text-orientation: mixed;
  }}
  .sidePanel {{
    display: flex;
    flex-direction: column;
    gap: 12px;
    max-height: calc(100vh - 102px);
  }}
  .detail, .tableBox {{
    border: 1px solid var(--line);
    background: var(--panel);
    border-radius: 8px;
    overflow: hidden;
  }}
  .detail {{
    padding: 12px;
  }}
  .detail h2 {{
    margin: 0 0 8px;
    font-size: 16px;
  }}
  .metrics {{
    display: grid;
    grid-template-columns: repeat(5, minmax(0, 1fr));
    gap: 8px;
    margin: 8px 0 10px;
  }}
  .metric {{
    background: white;
    border: 1px solid var(--line);
    border-radius: 6px;
    padding: 7px;
  }}
  .metric b {{
    display: block;
    font-size: 13px;
  }}
  .metric span {{
    color: var(--muted);
    font-size: 11px;
  }}
  .chips {{
    display: flex;
    gap: 6px;
    flex-wrap: wrap;
    margin: 6px 0;
  }}
  .chip {{
    border: 1px solid var(--line);
    background: white;
    border-radius: 999px;
    padding: 3px 8px;
    font-size: 12px;
  }}
  .tableBox {{
    flex: 1;
    overflow: auto;
    background: white;
  }}
  table {{
    width: 100%;
    border-collapse: collapse;
    font-size: 12px;
  }}
  th, td {{
    border-bottom: 1px solid #eef2f6;
    padding: 7px 8px;
    text-align: left;
    vertical-align: top;
  }}
  th {{
    position: sticky;
    top: 0;
    z-index: 1;
    background: #f8fafc;
    font-weight: 650;
  }}
  tr {{
    cursor: pointer;
  }}
  tr.selected {{
    background: #fff1f0;
  }}
  .small {{
    color: var(--muted);
    font-size: 12px;
  }}
  .members {{
    max-height: 220px;
    overflow: auto;
    background: white;
    border: 1px solid var(--line);
    border-radius: 6px;
    padding: 6px 8px;
    font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace;
    font-size: 11px;
    line-height: 1.45;
  }}
  @media (max-width: 1100px) {{
    main {{
      grid-template-columns: 1fr;
    }}
    .canvasWrap {{
      width: 100%;
    }}
    .heatmapArea {{
      flex-direction: column;
    }}
    .heatmapArea.zoomed {{
      display: flex;
    }}
    .heatmapArea.zoomed .canvasWrap {{
      width: 100%;
    }}
    .heatmapArea.zoomed .xAxisLabels, .heatmapArea.zoomed .yAxisLabels {{
      display: none;
    }}
    .legend {{
      flex-direction: row;
      min-width: 0;
    }}
    .ramp {{
      width: 160px;
      height: 12px;
      background: linear-gradient(90deg, rgb(87,129,188), rgb(239,237,240), rgb(185,91,90));
    }}
    .status {{
      writing-mode: horizontal-tb;
    }}
    .sidePanel {{
      max-height: none;
    }}
  }}
</style>
</head>
<body>
<header>
  <h1>Effect Vector Heatmap Module Browser</h1>
  <p class="sub">Data-driven modules are drawn as diagonal boundaries. The table shows unbiased module composition; the Theme column is a separate interpretation layer based on module gene composition.</p>
</header>
<main>
  <section class="heatmapPanel">
    <div class="toolbar">
      <button id="toggleBoundaries" class="active">Module Boundaries</button>
      <button id="toggleLabels" class="active">Axis Labels</button>
      <button id="zoomButton">Zoom Module</button>
      <label class="small">Sort table
        <select id="sortSelect">
          <option value="module_id" selected>module number</option>
          <option value="size">module size</option>
          <option value="within_mean">within similarity</option>
          <option value="neighbor_mean">neighbor similarity</option>
          <option value="contrast">contrast (within - neighbor)</option>
        </select>
      </label>
      <input id="searchBox" type="search" placeholder="Search perturbations/celltypes">
    </div>
    <div class="heatmapArea">
      <div id="xAxisLabels" class="xAxisLabels"></div>
      <div id="yAxisLabels" class="yAxisLabels"></div>
      <div class="canvasWrap" id="canvasWrap">
        <canvas id="heatmap"></canvas>
        <canvas id="overlay"></canvas>
      </div>
      <div class="legend"><span>+1</span><div class="ramp"></div><span>-1</span><span id="status" class="small status"></span></div>
    </div>
  </section>
  <aside class="sidePanel">
    <section class="detail" id="detail"></section>
    <section class="tableBox">
      <table id="moduleTable"></table>
    </section>
  </aside>
</main>
<script id="heatmap-data" type="application/json">{data_json}</script>
<script>
const DATA = JSON.parse(document.getElementById('heatmap-data').textContent);
const labels = DATA.labels;
const modules = DATA.modules;
const N = labels.length;
const heatmap = document.getElementById('heatmap');
const overlay = document.getElementById('overlay');
const hctx = heatmap.getContext('2d');
const octx = overlay.getContext('2d');
const statusEl = document.getElementById('status');
const heatmapArea = document.querySelector('.heatmapArea');
const xAxisLabelsEl = document.getElementById('xAxisLabels');
const yAxisLabelsEl = document.getElementById('yAxisLabels');
let selectedId = modules[0].module_id;
let showBoundaries = true;
let showLabels = true;
let zoomMode = false;
let sortKey = 'module_id';
let filterText = '';
let matrixValues = null;
let yPos = null;

function b64ToBytes(b64) {{
  const bin = atob(b64);
  const bytes = new Uint8Array(bin.length);
  for (let i = 0; i < bin.length; i++) bytes[i] = bin.charCodeAt(i);
  return bytes;
}}

function colorRamp(v) {{
  const t = Math.max(0, Math.min(1, (v + 1) / 2));
  const stops = [
    [0.0, [87,129,188]],
    [0.5, [239,237,240]],
    [1.0, [185,91,90]]
  ];
  let lo = stops[0], hi = stops[2];
  if (t <= 0.5) {{ lo = stops[0]; hi = stops[1]; }}
  else {{ lo = stops[1]; hi = stops[2]; }}
  const u = (t - lo[0]) / (hi[0] - lo[0]);
  return [
    Math.round(lo[1][0] + (hi[1][0] - lo[1][0]) * u),
    Math.round(lo[1][1] + (hi[1][1] - lo[1][1]) * u),
    Math.round(lo[1][2] + (hi[1][2] - lo[1][2]) * u)
  ];
}}

function decodeMatrix() {{
  const z = DATA.z;
  const bytes = b64ToBytes(z.bdata);
  if (z.dtype === 'f8') matrixValues = new Float64Array(bytes.buffer);
  else if (z.dtype === 'f4') matrixValues = new Float32Array(bytes.buffer);
  else throw new Error('Unsupported dtype: ' + z.dtype);
  yPos = new Map(DATA.yLabels.map((label, i) => [label, i]));
  renderHeatmap();
}}

function matrixValue(row, col) {{
  const sourceRow = yPos.get(labels[row]);
  return matrixValues[sourceRow * N + col];
}}

function currentExtent() {{
  if (!zoomMode) return {{start: 0, end: N, size: N}};
  const m = selectedModule();
  return {{start: m.start, end: m.end, size: m.size}};
}}

function renderHeatmap() {{
  const extent = currentExtent();
  const image = hctx.createImageData(extent.size, extent.size);
  for (let row = 0; row < extent.size; row++) {{
    for (let col = 0; col < extent.size; col++) {{
      const v = matrixValue(extent.start + row, extent.start + col);
      const [r, g, b] = colorRamp(v);
      const p = (row * extent.size + col) * 4;
      image.data[p] = r;
      image.data[p + 1] = g;
      image.data[p + 2] = b;
      image.data[p + 3] = 255;
    }}
  }}
  heatmap.width = extent.size;
  heatmap.height = extent.size;
  overlay.width = extent.size;
  overlay.height = extent.size;
  hctx.putImageData(image, 0, 0);
  statusEl.textContent = '';
  heatmapArea.classList.toggle('zoomed', zoomMode);
  drawOverlay();
  renderAxisLabels();
}}

function fmt(x, digits = 3) {{
  return x === null || Number.isNaN(x) ? '' : Number(x).toFixed(digits);
}}

function esc(s) {{
  return String(s ?? '').replace(/[&<>"']/g, c => ({{'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}}[c]));
}}

function topCount(summaryText) {{
  const match = String(summaryText ?? '').match(/:(\\d+)(?:;|$)/);
  return match ? Number(match[1]) : null;
}}

function formatCounts(summaryText) {{
  return String(summaryText ?? '').split(';').slice(0, 5).map(part => {{
    const trimmed = part.trim();
    const match = trimmed.match(/^(.*):(\\d+)$/);
    if (!match) return trimmed;
    return match[1].trim() + ' (n=' + match[2] + ')';
  }}).join('; ');
}}

function stripCellCode(cellType) {{
  return String(cellType ?? '').replace(/^\\d+\\s+/, '');
}}

function formatCellCounts(summaryText) {{
  return String(summaryText ?? '').split(';').slice(0, 5).map(part => {{
    const trimmed = part.trim();
    const match = trimmed.match(/^(.*):(\\d+)$/);
    if (!match) return stripCellCode(trimmed);
    return stripCellCode(match[1].trim()) + ' (n=' + match[2] + ')';
  }}).join('; ');
}}

function formatMemberLabel(member) {{
  return stripCellCode(member.cell_type) + ' - ' + member.gene;
}}

function topCellTypeFraction(m) {{
  const count = topCount(m.top_cells);
  return count === null ? null : count / Math.max(1, m.size);
}}

function selectedModule() {{
  return modules.find(m => m.module_id === selectedId) || modules[0];
}}

function drawOverlay() {{
  const extent = currentExtent();
  octx.clearRect(0, 0, extent.size, extent.size);
  if (zoomMode) {{
    octx.save();
    octx.lineWidth = 3;
    octx.strokeStyle = 'rgba(180,35,24,0.95)';
    octx.strokeRect(0.5, 0.5, extent.size - 1, extent.size - 1);
    octx.restore();
  }} else if (showBoundaries) {{
    octx.save();
    octx.lineWidth = 2;
    octx.strokeStyle = 'rgba(23,32,42,0.72)';
    for (const m of modules) {{
      octx.strokeRect(m.start + 0.5, m.start + 0.5, m.size, m.size);
    }}
    octx.restore();
  }}
  const m = selectedModule();
  octx.save();
  octx.lineWidth = 5;
  octx.strokeStyle = 'rgba(180,35,24,0.95)';
  if (zoomMode) {{
    octx.strokeRect(0.5, 0.5, extent.size - 1, extent.size - 1);
  }} else {{
    octx.strokeRect(m.start + 0.5, m.start + 0.5, m.size, m.size);
    octx.fillStyle = 'rgba(180,35,24,0.11)';
    octx.fillRect(m.start, m.start, m.size, m.size);
  }}
  octx.restore();

  if (showLabels) {{
    octx.save();
    octx.fillStyle = 'rgba(23,32,42,0.9)';
    octx.font = '18px sans-serif';
    octx.textBaseline = 'top';
    if (zoomMode) {{
      octx.fillText(String(selectedModule().module_id), 5, 5);
    }} else {{
      for (const m of modules) {{
        if (m.size < 16) continue;
        octx.fillText(String(m.module_id), m.start + 5, m.start + 5);
      }}
    }}
    octx.restore();
  }}
}}

function renderAxisLabels() {{
  if (!zoomMode) {{
    xAxisLabelsEl.innerHTML = '';
    yAxisLabelsEl.innerHTML = '';
    return;
  }}
  const m = selectedModule();
  const size = Math.max(1, m.members.length);
  xAxisLabelsEl.innerHTML = m.members.map((x, i) => {{
    const left = ((i + 0.5) / size * 100).toFixed(4);
    return `<span class="xTick" style="left:${{left}}%">${{esc(formatMemberLabel(x))}}</span>`;
  }}).join('');
  yAxisLabelsEl.innerHTML = m.members.map((x, i) => {{
    const top = ((i + 0.5) / size * 100).toFixed(4);
    return `<span class="yTick" style="top:${{top}}%">${{esc(formatMemberLabel(x))}}</span>`;
  }}).join('');
}}

function renderDetail() {{
  const m = selectedModule();
  const memberLines = m.members.slice(0, 140).map(x => `${{x.position}}  ${{formatMemberLabel(x)}}`).join('\\n');
  const more = m.members.length > 140 ? `\\n... ${{m.members.length - 140}} more` : '';
  document.getElementById('detail').innerHTML = `
    <h2>Module ${{m.module_id}} <span class="small">positions ${{m.start}}-${{m.end - 1}}</span></h2>
    <div class="chips">
      ${{m.auto_theme ? `<span class="chip">Theme: ${{esc(m.auto_theme)}}</span>` : ''}}
      <span class="chip">${{m.size}} perturbation-cell type combinations</span>
      <span class="chip">${{m.n_unique_genes}} perturbations</span>
      <span class="chip">${{m.n_unique_cell_types}} cell types</span>
    </div>
    <div class="metrics">
      <div class="metric"><b>${{fmt(m.within_mean)}}</b><span>within similarity</span></div>
      <div class="metric"><b>${{fmt(m.neighbor_mean)}}</b><span>neighbor similarity</span></div>
      <div class="metric"><b>${{fmt(m.contrast)}}</b><span>contrast</span></div>
      <div class="metric"><b>${{fmt(m.top_gene_fraction, 2)}}</b><span>top perturbation frac</span></div>
      <div class="metric"><b>${{fmt(topCellTypeFraction(m), 2)}}</b><span>top cell type frac</span></div>
    </div>
    <p class="small"><b>Top perturbations:</b> ${{esc(formatCounts(m.top_genes))}}</p>
    <p class="small"><b>Top cell types:</b> ${{esc(formatCellCounts(m.top_cells))}}</p>
    <pre class="members">${{esc(memberLines + more)}}</pre>
  `;
}}

function moduleMatches(m) {{
  if (!filterText) return true;
  const hay = [m.module_id, m.auto_theme, m.top_genes, m.top_cells, m.theme_hits].join(' ').toLowerCase();
  return hay.includes(filterText);
}}

function sortedModules() {{
  const rows = modules.filter(moduleMatches);
  return rows.sort((a, b) => {{
    if (sortKey === 'module_id') return a.module_id - b.module_id;
    const av = a[sortKey] ?? -Infinity;
    const bv = b[sortKey] ?? -Infinity;
    return bv - av;
  }});
}}

function renderTable() {{
  const rows = sortedModules();
  const body = rows.map(m => `
    <tr data-id="${{m.module_id}}" class="${{m.module_id === selectedId ? 'selected' : ''}}">
      <td>${{m.module_id}}</td>
      <td>${{m.start}}-${{m.end - 1}}</td>
      <td>${{m.size}}</td>
      <td>${{fmt(m.within_mean)}}</td>
      <td>${{fmt(m.neighbor_mean)}}</td>
      <td>${{fmt(m.contrast)}}</td>
      <td>${{esc(m.auto_theme || '')}}</td>
      <td>${{esc(formatCounts(m.top_genes))}}</td>
      <td>${{esc(formatCellCounts(m.top_cells))}}</td>
    </tr>`).join('');
  document.getElementById('moduleTable').innerHTML = `
    <thead>
      <tr><th>module number</th><th>range</th><th>module size</th><th>within similarity</th><th>neighbor similarity</th><th>contrast</th><th>Theme</th><th>top perturbations</th><th>top cell types</th></tr>
    </thead>
    <tbody>${{body}}</tbody>
  `;
  document.querySelectorAll('#moduleTable tbody tr').forEach(tr => {{
    tr.addEventListener('click', () => {{
      selectedId = Number(tr.dataset.id);
      renderAll();
    }});
  }});
}}

function renderAll() {{
  renderHeatmap();
  renderDetail();
  renderTable();
}}

document.getElementById('toggleBoundaries').addEventListener('click', e => {{
  showBoundaries = !showBoundaries;
  e.currentTarget.classList.toggle('active', showBoundaries);
  drawOverlay();
}});
document.getElementById('toggleLabels').addEventListener('click', e => {{
  showLabels = !showLabels;
  e.currentTarget.classList.toggle('active', showLabels);
  drawOverlay();
}});
document.getElementById('zoomButton').addEventListener('click', e => {{
  zoomMode = !zoomMode;
  e.currentTarget.classList.toggle('active', zoomMode);
  e.currentTarget.textContent = zoomMode ? 'Full Heatmap' : 'Zoom Module';
  renderHeatmap();
}});
document.getElementById('sortSelect').addEventListener('change', e => {{
  sortKey = e.target.value;
  renderTable();
}});
document.getElementById('searchBox').addEventListener('input', e => {{
  filterText = e.target.value.trim().toLowerCase();
  renderTable();
}});
overlay.addEventListener('mousemove', e => {{
  const rect = overlay.getBoundingClientRect();
  const extent = currentExtent();
  const x = extent.start + Math.floor((e.clientX - rect.left) / rect.width * extent.size);
  const y = extent.start + Math.floor((e.clientY - rect.top) / rect.height * extent.size);
  const m = modules.find(mm => x >= mm.start && x < mm.end && y >= mm.start && y < mm.end);
  overlay.title = zoomMode ? `${{labels[y]}} x ${{labels[x]}}` : (m ? `Module ${{m.module_id}}: ${{m.top_genes}}` : `${{x}}, ${{y}}`);
}});
overlay.addEventListener('click', e => {{
  if (zoomMode) return;
  const rect = overlay.getBoundingClientRect();
  const x = Math.floor((e.clientX - rect.left) / rect.width * N);
  const y = Math.floor((e.clientY - rect.top) / rect.height * N);
  const m = modules.find(mm => x >= mm.start && x < mm.end && y >= mm.start && y < mm.end);
  if (m) {{
    selectedId = m.module_id;
    renderAll();
  }}
}});

decodeMatrix();
renderAll();
</script>
</body>
</html>
"""
    output_html.write_text(html_text)
    print(output_html)
