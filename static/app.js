// Dengue Sequence Analyzer - frontend
// Fetches /api/analyze, renders an illustrated results page: summary tiles,
// serotype-confidence bars, an SVG genome map (protein blocks, coverage,
// variant lollipops, epitope bands), and sortable/filterable tables.

const IMPACT_ORDER = ['HIGH', 'MODERATE', 'LOW', 'MODIFIER'];
const IMPACT_COLOR = {
  HIGH: 'var(--status-critical)',
  MODERATE: 'var(--status-serious)',
  LOW: 'var(--status-warning)',
  MODIFIER: 'var(--status-good)',
};
const SEROTYPE_COLOR = { 1: 'var(--series-1)', 2: 'var(--series-2)', 3: 'var(--series-3)', 4: 'var(--series-4)' };
const EPITOPE_COLOR = 'var(--series-5)';

let currentResult = null;

const $ = (sel) => document.querySelector(sel);
const tooltip = $('#tooltip');

function showTooltip(html, x, y) {
  tooltip.innerHTML = html;
  tooltip.style.left = (x + 14) + 'px';
  tooltip.style.top = (y + 14) + 'px';
  tooltip.style.display = 'block';
}
function hideTooltip() { tooltip.style.display = 'none'; }

// ---------- input handling ----------

$('#file-btn').addEventListener('click', () => $('#file-input').click());
$('#file-input').addEventListener('change', (e) => {
  const file = e.target.files[0];
  if (!file) return;
  const reader = new FileReader();
  reader.onload = () => { $('#seq-input').value = reader.result; };
  reader.readAsText(file);
});

document.querySelectorAll('[data-example]').forEach((btn) => {
  btn.addEventListener('click', async () => {
    const serotype = btn.getAttribute('data-example');
    setStatus(`Loading a DENV-${serotype} example...`);
    try {
      const res = await fetch(`/api/example?serotype=${serotype}`);
      const data = await res.json();
      if (!res.ok) throw new Error(data.error || 'Failed to load example.');
      $('#seq-input').value = data.fasta;
      setStatus(`Loaded ${data.accession} (DENV-${serotype}, ${data.country || 'unknown origin'}).`);
    } catch (err) {
      setStatus(err.message, true);
    }
  });
});

$('#analyze-btn').addEventListener('click', analyze);

function setStatus(msg, isError) {
  const el = $('#status-line');
  el.textContent = msg;
  el.style.color = isError ? 'var(--status-critical)' : 'var(--text-secondary)';
}

async function analyze() {
  const sequence = $('#seq-input').value;
  if (!sequence.trim()) { setStatus('Paste or upload a sequence first.', true); return; }

  $('#analyze-btn').disabled = true;
  setStatus('Analyzing... (aligning against the reference genome and cross-referencing the database)');
  $('#results').classList.remove('show');

  try {
    const res = await fetch('/api/analyze', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ sequence }),
    });
    const data = await res.json();
    if (!res.ok) throw new Error(data.error || 'Analysis failed.');
    currentResult = data;
    renderAll(data);
    setStatus('Done.');
    $('#results').classList.add('show');
    $('#results').scrollIntoView({ behavior: 'smooth', block: 'start' });
  } catch (err) {
    setStatus(err.message, true);
  } finally {
    $('#analyze-btn').disabled = false;
  }
}

// ---------- render orchestration ----------

function renderAll(r) {
  renderWarnings(r.warnings);
  renderSummaryTiles(r);
  renderSerotypeBars(r.serotype_call);
  renderGenomeMap(r);
  renderVariantFilters(r.variants);
  renderVariantTable(r.variants);
  renderNeighborTable(r.nearest_neighbors);
}

function renderWarnings(warnings) {
  const el = $('#warnings');
  el.innerHTML = '';
  if (!warnings || !warnings.length) { el.classList.remove('show'); return; }
  el.classList.add('show');
  warnings.forEach((w) => {
    const div = document.createElement('div');
    div.className = 'warning-item';
    div.innerHTML = `<span class="icon">&#9888;</span><span>${escapeHtml(w)}</span>`;
    el.appendChild(div);
  });
}

function renderSummaryTiles(r) {
  const tiles = $('#summary-tiles');
  const s = r.summary_counts;
  const impactBreakdown = IMPACT_ORDER
    .filter((k) => s.by_impact[k])
    .map((k) => `<span><span class="dot" style="background:${IMPACT_COLOR[k]}"></span>${k} ${s.by_impact[k]}</span>`)
    .join('');

  tiles.innerHTML = `
    <div class="tile">
      <div class="label">Serotype call</div>
      <div class="value">DENV-${r.serotype_call.serotype}</div>
    </div>
    <div class="tile">
      <div class="label">Genome coverage</div>
      <div class="value">${r.coverage.pct}<span class="unit">%</span></div>
      <div class="breakdown">positions ${r.coverage.start.toLocaleString()}&ndash;${r.coverage.end.toLocaleString()} of ${r.coverage.reference_length.toLocaleString()}</div>
    </div>
    <div class="tile">
      <div class="label">Input length / GC%</div>
      <div class="value">${r.input_summary.length.toLocaleString()}<span class="unit"> bp</span></div>
      <div class="breakdown">GC ${r.input_summary.gc_pct}% &middot; N ${r.input_summary.n_pct}%</div>
    </div>
    <div class="tile">
      <div class="label">Total variants</div>
      <div class="value">${s.total.toLocaleString()}</div>
      <div class="breakdown">${impactBreakdown}</div>
    </div>
    <div class="tile">
      <div class="label">Novel variants</div>
      <div class="value">${s.novel.toLocaleString()}</div>
      <div class="breakdown">not in variant_catalog</div>
    </div>
    <div class="tile">
      <div class="label">Epitope-overlapping</div>
      <div class="value">${s.epitope_overlap.toLocaleString()}</div>
      <div class="breakdown">IEDB B-/T-cell epitopes</div>
    </div>
  `;
}

function renderSerotypeBars(call) {
  const wrap = $('#sero-bars');
  wrap.innerHTML = '';
  const maxScore = Math.max(...Object.values(call.scores_by_serotype), 0.01);

  [1, 2, 3, 4].forEach((st) => {
    const score = call.scores_by_serotype[st] || 0;
    const isCalled = st === call.serotype;
    const row = document.createElement('div');
    row.className = 'sero-bar-row' + (isCalled ? ' called' : '');
    row.innerHTML = `
      <div class="sero-label${isCalled ? ' called' : ''}">DENV-${st}${isCalled ? ' ✓' : ''}</div>
      <div class="sero-bar-track">
        <div class="sero-bar-fill" style="width:${(score / maxScore * 100).toFixed(1)}%; background:${SEROTYPE_COLOR[st]}"></div>
      </div>
      <div class="sero-value">${score}%</div>
    `;
    wrap.appendChild(row);
  });

  const note = document.createElement('div');
  note.style.cssText = 'font-size:11.5px;color:var(--text-muted);margin-top:6px;';
  note.textContent = call.margin_ratio
    ? `Winner leads the runner-up by ${call.margin_ratio}× (orientation: ${call.orientation.replace('_', ' ')}).`
    : `Orientation: ${call.orientation.replace('_', ' ')}.`;
  wrap.appendChild(note);
}

// ---------- genome map ----------

function renderGenomeMap(r) {
  const proteins = r.reference_info.proteins;
  const refLen = r.reference_info.length;
  const variants = r.variants;
  const covStart = r.coverage.start - 1; // 0-idx
  const covEnd = r.coverage.end;

  const svgWidth = 1000;
  const marginX = 16;
  const plotWidth = svgWidth - 2 * marginX;
  const scaleX = (pos0) => marginX + (pos0 / refLen) * plotWidth;

  const lollipopY = 34;
  const stemBottomY = 66;
  const proteinY = 66, proteinH = 28;
  const epitopeY = 100, epitopeH = 10;
  const axisY = 122;
  const svgHeight = 138;

  let svg = `<svg id="genome-map" viewBox="0 0 ${svgWidth} ${svgHeight}" width="${svgWidth}" height="${svgHeight}" xmlns="http://www.w3.org/2000/svg">`;
  svg += `<defs>
    <pattern id="hatch-uncovered" width="6" height="6" patternTransform="rotate(45)" patternUnits="userSpaceOnUse">
      <rect width="6" height="6" fill="var(--gridline)"></rect>
      <line x1="0" y1="0" x2="0" y2="6" stroke="var(--text-muted)" stroke-width="1.5"></line>
    </pattern>
  </defs>`;

  // baseline genome bar
  svg += `<rect x="${marginX}" y="${proteinY - 3}" width="${plotWidth}" height="${proteinH + 6}" rx="4" fill="var(--gridline)"></rect>`;

  // protein blocks, alternating neutral shades, direct text labels
  const proteinEntries = Object.entries(proteins).sort((a, b) => a[1].start - b[1].start);
  proteinEntries.forEach(([name, coords], i) => {
    const x1 = scaleX(coords.start), x2 = scaleX(coords.end);
    const w = Math.max(x2 - x1, 1);
    const shade = i % 2 === 0 ? 'var(--baseline)' : 'var(--text-muted)';
    svg += `<g class="protein-block" data-tip="${name}: ref ${coords.start + 1}-${coords.end} (${coords.length} nt)">
      <rect x="${x1}" y="${proteinY}" width="${w}" height="${proteinH}" fill="${shade}" opacity="0.55"></rect>`;
    if (w > 18) {
      svg += `<text x="${(x1 + x2) / 2}" y="${proteinY + proteinH / 2 + 3}" text-anchor="middle">${name}</text>`;
    }
    svg += `</g>`;
  });

  // epitope bands (single overlay series)
  const epitopePositions = collectEpitopeRanges(r.variants, proteins);
  epitopePositions.forEach((range) => {
    const x1 = scaleX(range.start), x2 = scaleX(range.end);
    svg += `<rect class="epitope-band" x="${x1}" y="${epitopeY}" width="${Math.max(x2 - x1, 2)}" height="${epitopeH}"
              fill="${EPITOPE_COLOR}" opacity="0.55" rx="2"
              data-tip="Epitope region: ${range.protein} ${range.startAA}-${range.endAA}"></rect>`;
  });

  // coverage: hatch the uncovered flanks over the full track height
  if (covStart > 0) {
    const w = scaleX(covStart) - marginX;
    svg += `<rect class="coverage-uncovered" x="${marginX}" y="${proteinY - 3}" width="${w}" height="${epitopeY + epitopeH - proteinY + 6}"
              data-tip="Not covered by query sequence"></rect>`;
  }
  if (covEnd < refLen) {
    const x = scaleX(covEnd);
    const w = scaleX(refLen) - x;
    svg += `<rect class="coverage-uncovered" x="${x}" y="${proteinY - 3}" width="${w}" height="${epitopeY + epitopeH - proteinY + 6}"
              data-tip="Not covered by query sequence"></rect>`;
  }

  // variant lollipops
  variants.forEach((v) => {
    const x = scaleX(v.position - 1);
    const color = IMPACT_COLOR[v.impact] || IMPACT_COLOR.MODIFIER;
    const label = v.aa_change ? `${v.protein || ''} ${v.aa_change}` : `pos ${v.position} ${v.ref_allele}>${v.alt_allele}`;
    const tip = `${label} &middot; ${v.impact}${v.known ? ` &middot; known (${v.frequency_pct}% of DENV-${r.serotype_call.serotype} genomes)` : ' &middot; novel'}${v.epitope_hits.length ? ' &middot; epitope overlap' : ''}`;
    svg += `<g class="variant-marker" data-tip="${tip}">
      <line x1="${x}" y1="${lollipopY + 5}" x2="${x}" y2="${stemBottomY}" stroke="${color}" stroke-width="1.5" opacity="0.7"></line>
      <circle cx="${x}" cy="${lollipopY}" r="4.5" fill="${color}"></circle>
    </g>`;
  });

  // axis ticks
  svg += `<line x1="${marginX}" y1="${axisY}" x2="${svgWidth - marginX}" y2="${axisY}" stroke="var(--baseline)" stroke-width="1"></line>`;
  [0, 0.25, 0.5, 0.75, 1].forEach((f) => {
    const pos = Math.round(f * refLen);
    const x = scaleX(pos);
    svg += `<line x1="${x}" y1="${axisY}" x2="${x}" y2="${axisY + 4}" stroke="var(--baseline)" stroke-width="1"></line>
            <text x="${x}" y="${axisY + 14}" text-anchor="middle" font-size="9.5" fill="var(--text-muted)">${pos.toLocaleString()}</text>`;
  });

  svg += `</svg>`;

  $('#genome-map-wrap').innerHTML = svg;

  // legend
  $('#map-legend').innerHTML = IMPACT_ORDER.map((k) =>
    `<span class="item"><span class="swatch" style="background:${IMPACT_COLOR[k]}"></span>${k}</span>`
  ).join('') + `<span class="item"><span class="band-swatch" style="background:${EPITOPE_COLOR};opacity:0.55"></span>Known epitope region</span>`
    + `<span class="item"><span class="band-swatch" style="background:repeating-linear-gradient(45deg,var(--gridline),var(--gridline) 3px,var(--text-muted) 3px,var(--text-muted) 4px)"></span>Not covered by query</span>`;

  // tooltips (event delegation)
  const mapSvg = $('#genome-map');
  mapSvg.querySelectorAll('[data-tip]').forEach((node) => {
    node.addEventListener('mousemove', (e) => showTooltip(node.getAttribute('data-tip'), e.clientX, e.clientY));
    node.addEventListener('mouseleave', hideTooltip);
  });
}

function collectEpitopeRanges(variants, proteins) {
  // Approximate epitope band positions from the variant list's own epitope
  // hits (protein-relative aa coordinates aren't directly on the variant, so
  // draw a small marker centered on the variant's own genome position for
  // every variant that has at least one epitope hit).
  const ranges = [];
  variants.forEach((v) => {
    if (v.epitope_hits && v.epitope_hits.length && v.protein && proteins[v.protein]) {
      const proteinStart = proteins[v.protein].start;
      const genomePos = proteinStart + (v.protein_position - 1) * 3;
      ranges.push({
        start: genomePos, end: genomePos + 3,
        protein: v.protein, startAA: v.protein_position, endAA: v.protein_position,
      });
    }
  });
  return ranges;
}

// ---------- variant table ----------

function renderVariantFilters(variants) {
  const proteinSelect = $('#filter-protein');
  const proteins = Array.from(new Set(variants.map((v) => v.protein).filter(Boolean))).sort();
  proteinSelect.innerHTML = '<option value="">All proteins</option>' +
    proteins.map((p) => `<option value="${p}">${p}</option>`).join('');
  ['filter-impact', 'filter-protein', 'filter-novel'].forEach((id) => {
    $(`#${id}`).onchange = () => renderVariantTable(currentResult.variants);
  });
}

let sortState = { col: 'position', dir: 1 };

function getFilteredVariants(variants) {
  const impact = $('#filter-impact').value;
  const protein = $('#filter-protein').value;
  const novel = $('#filter-novel').value;
  return variants.filter((v) => {
    if (impact && v.impact !== impact) return false;
    if (protein && v.protein !== protein) return false;
    if (novel === 'novel' && v.known) return false;
    if (novel === 'known' && !v.known) return false;
    return true;
  });
}

const VARIANT_COLUMNS = [
  { key: 'position', label: 'Position' },
  { key: 'protein', label: 'Protein' },
  { key: 'aa_change', label: 'AA change' },
  { key: 'variant_type', label: 'Type' },
  { key: 'impact', label: 'Impact' },
  { key: 'known', label: 'Known?' },
  { key: 'frequency_pct', label: 'Freq %' },
  { key: 'country_count', label: '# Countries' },
  { key: 'epitope_hits', label: 'Epitope' },
];

function renderVariantTable(allVariants) {
  const rows = getFilteredVariants(allVariants);
  rows.sort((a, b) => {
    const av = a[sortState.col], bv = b[sortState.col];
    if (av == null && bv == null) return 0;
    if (av == null) return 1;
    if (bv == null) return -1;
    if (av < bv) return -1 * sortState.dir;
    if (av > bv) return 1 * sortState.dir;
    return 0;
  });

  $('#variant-count-sub').textContent = `${rows.length} of ${allVariants.length} shown`;

  const thead = document.querySelector('#variant-table thead');
  thead.innerHTML = '<tr>' + VARIANT_COLUMNS.map((c) =>
    `<th data-col="${c.key}">${c.label}${sortState.col === c.key ? (sortState.dir === 1 ? ' ▲' : ' ▼') : ''}</th>`
  ).join('') + '</tr>';
  thead.querySelectorAll('th').forEach((th) => {
    th.onclick = () => {
      const col = th.getAttribute('data-col');
      sortState.dir = (sortState.col === col) ? -sortState.dir : 1;
      sortState.col = col;
      renderVariantTable(allVariants);
    };
  });

  const tbody = document.querySelector('#variant-table tbody');
  tbody.innerHTML = rows.map((v) => `
    <tr>
      <td>${v.position.toLocaleString()}</td>
      <td>${v.protein || '&mdash;'}${v.protein_position ? ' ' + v.protein_position : ''}</td>
      <td>${v.aa_change || '&mdash;'}</td>
      <td>${v.variant_type}${v.low_confidence ? ' <span title="Alt allele is an IUPAC ambiguity code, not a clean base call">&#9888;</span>' : ''}</td>
      <td><span class="impact-badge"><span class="dot" style="background:${IMPACT_COLOR[v.impact]}"></span>${v.impact}</span></td>
      <td>${v.known ? 'known' : '<span class="novel-flag">novel</span>'}</td>
      <td>${v.known ? v.frequency_pct + '%' : '&mdash;'}</td>
      <td>${v.known ? (v.country_count ?? '&mdash;') : '&mdash;'}</td>
      <td>${summarizeEpitopeHits(v.epitope_hits)}</td>
    </tr>
  `).join('');
}

function summarizeEpitopeHits(hits) {
  // IEDB has many independently-submitted, overlapping entries per real
  // epitope site, so a busy site can carry dozens of hits -- collapse to
  // "type ×count" per distinct epitope_type rather than listing every hit.
  if (!hits || !hits.length) return '&mdash;';
  const counts = {};
  hits.forEach((h) => { counts[h.epitope_type] = (counts[h.epitope_type] || 0) + 1; });
  return Object.entries(counts).map(([type, n]) => `${type} &times;${n}`).join(', ');
}

// ---------- nearest neighbor table ----------

function renderNeighborTable(neighbors) {
  const thead = document.querySelector('#neighbor-table thead');
  const tbody = document.querySelector('#neighbor-table tbody');
  thead.innerHTML = '<tr><th>Accession</th><th>Country</th><th>Collection date</th><th>Lineage</th><th>Shared variants</th><th>Similarity</th></tr>';
  if (!neighbors || !neighbors.length) {
    tbody.innerHTML = '<tr><td colspan="6" style="color:var(--text-muted)">No comparable genomes (query had no SNPs to match on, or fingerprint cache is unavailable).</td></tr>';
    return;
  }
  tbody.innerHTML = neighbors.map((n) => `
    <tr>
      <td>${n.accession}</td>
      <td>${n.country || '&mdash;'}</td>
      <td>${n.collection_date || '&mdash;'}</td>
      <td>${n.lineage || n.genotype || '&mdash;'}</td>
      <td>${n.shared_variants} / ${n.query_variant_count} (it has ${n.candidate_variant_count} total)</td>
      <td>${n.similarity_pct}%</td>
    </tr>
  `).join('');
}

// ---------- export ----------

$('#export-csv').addEventListener('click', () => {
  if (!currentResult) return;
  const rows = getFilteredVariants(currentResult.variants);
  const cols = ['position', 'protein', 'protein_position', 'aa_change', 'variant_type', 'impact',
    'is_synonymous', 'known', 'frequency_pct', 'country_count', 'first_observed_date', 'last_observed_date'];
  const csv = [cols.join(',')].concat(
    rows.map((v) => cols.map((c) => csvEscape(v[c])).join(','))
  ).join('\n');
  downloadFile(csv, 'variants.csv', 'text/csv');
});

$('#export-json').addEventListener('click', () => {
  if (!currentResult) return;
  downloadFile(JSON.stringify(currentResult, null, 2), 'dengue_analysis.json', 'application/json');
});

function csvEscape(val) {
  if (val === null || val === undefined) return '';
  const s = String(val);
  return /[",\n]/.test(s) ? '"' + s.replace(/"/g, '""') + '"' : s;
}

function downloadFile(content, filename, mime) {
  const blob = new Blob([content], { type: mime });
  const url = URL.createObjectURL(blob);
  const a = document.createElement('a');
  a.href = url; a.download = filename;
  document.body.appendChild(a); a.click(); document.body.removeChild(a);
  URL.revokeObjectURL(url);
}

function escapeHtml(s) {
  const div = document.createElement('div');
  div.textContent = s;
  return div.innerHTML;
}
