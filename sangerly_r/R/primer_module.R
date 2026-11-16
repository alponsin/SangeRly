# primer_module.R - "Primer mining and design" tab (step 8) of SangeRly
#
#   1. Find sequences and primers: BOLD, NCBI, FASTA, MSA tab, assembly consensus; published primers
#   2. Design: conventional PCR, Sanger sequencing (walking), mRNA exon-exon junction, degenerate (+ probe)
#   3. Quality: Tm, GC, hairpin, dimers, 3' end, binding in a sequence set
#   4. BLAST: NCBI web BLAST of the primers and predicted amplicons
# The calculations are in primer_core.R.

.pr_msa_js <- r"-----(
(function () {
  var S = {};
  var COL = { A: ['#d9efd4', '#2f9e57'], C: ['#d3e3f6', '#2f6fc4'], G: ['#fbe6cb', '#e58a2c'], T: ['#f7d6d6', '#d64545'] };
  var IU = { A: 'A', C: 'C', G: 'G', T: 'T', AG: 'R', CT: 'Y', CG: 'S', AT: 'W', GT: 'K', AC: 'M', CGT: 'B', AGT: 'D', ACT: 'H', ACG: 'V', ACGT: 'N' };
  var NW = 150, RULER = 20, TRK = 20, BOT = 60, CH = 16;
  function clamp(v, a, b) { return Math.min(Math.max(v, a), b); }
  function trunc(s, n) { return s.length > n ? s.substr(0, n - 1) + '…' : s; }
  function stats(st) {
    var seqs = st.d.seqs, n = seqs.length, L = st.L, mf = st.minfreq || 0.1, B = 'ACGT';
    st.maj = new Array(L); st.cons = new Array(L); st.freq = new Float32Array(L); st.pc = new Array(L);
    for (var j = 0; j < L; j++) {
      var c = [0, 0, 0, 0], gap = 0;
      for (var i = 0; i < n; i++) {
        var ch = seqs[i].charCodeAt(j);
        if (ch === 65) c[0]++; else if (ch === 67) c[1]++; else if (ch === 71) c[2]++; else if (ch === 84 || ch === 85) c[3]++; else if (ch === 45 || ch === 46) gap++;
      }
      var nb = c[0] + c[1] + c[2] + c[3], mx = 0, mi = -1, k;
      for (k = 0; k < 4; k++) if (c[k] > mx) { mx = c[k]; mi = k; }
      st.maj[j] = mi >= 0 ? B[mi] : '-'; st.freq[j] = mx / n;
      var set = ''; for (k = 0; k < 4; k++) if (nb && c[k] / nb >= mf) set += B[k];
      st.cons[j] = (gap >= n / 2 || !nb) ? '-' : (set ? IU[set] : st.maj[j]);
      var txt = []; for (k = 0; k < 4; k++) if (c[k]) txt.push(B[k] + ' ' + Math.round(100 * c[k] / n) + '%'); if (gap) txt.push('gap ' + Math.round(100 * gap / n) + '%');
      st.pc[j] = txt.join(', ');
    }
  }
  function layout(st) {
    var W = st.box.clientWidth, n = st.d ? st.d.seqs.length : 0;
    var top = RULER + st.tracks.length * TRK + 6;
    var H = clamp(top + n * CH + BOT + 8, 230, 560);
    return { W: W, H: H, top: top, botY: H - BOT, seqH: H - BOT - top };
  }
  function draw(st) {
    if (!st.d || !st.box.clientWidth) return;
    var o = layout(st), W = o.W, H = o.H, top = o.top, botY = o.botY, seqH = o.seqH, cv = st.cv;
    var dpr = window.devicePixelRatio || 1;
    if (cv.width !== Math.round(W * dpr) || cv.height !== Math.round(H * dpr)) { cv.width = Math.round(W * dpr); cv.height = Math.round(H * dpr); cv.style.width = W + 'px'; cv.style.height = H + 'px'; }
    var g = cv.getContext('2d'); g.setTransform(dpr, 0, 0, dpr, 0, 0); g.clearRect(0, 0, W, H);
    var n = st.d.seqs.length, L = st.L, cw = st.cw, ch = CH;
    st.x0 = clamp(st.x0, 0, Math.max(0, L * cw - (W - NW))); st.y0 = clamp(st.y0, 0, Math.max(0, n * ch - seqH));
    var c0 = Math.max(0, Math.floor(st.x0 / cw)), c1 = Math.min(L - 1, Math.ceil((st.x0 + W - NW) / cw));
    var r0 = Math.max(0, Math.floor(st.y0 / ch)), r1 = Math.min(n - 1, Math.ceil((st.y0 + seqH) / ch)), i, j, x, y, b;
    var showL = cw >= 8, wd = function (c) { return Math.max(cw - (cw > 4 ? 1 : 0), 1); };
    g.textAlign = 'center'; g.textBaseline = 'middle'; g.font = Math.min(13, Math.floor(cw * 0.85)) + 'px monospace';
    // sequences
    g.save(); g.beginPath(); g.rect(NW, top, W - NW, seqH); g.clip();
    for (i = r0; i <= r1; i++) {
      y = top + i * ch - st.y0; var s = st.d.seqs[i];
      for (j = c0; j <= c1; j++) {
        x = NW + j * cw - st.x0; b = s[j]; var fill, tc = '#333';
        if (b === '-' || b === '.') fill = '#ececec';
        else if (b === st.maj[j]) fill = COL[b] ? COL[b][0] : '#e5e5e5';
        else { fill = COL[b] ? COL[b][1] : '#8e8e8e'; tc = '#fff'; }
        g.fillStyle = fill; g.fillRect(x, y, wd(), ch - 1);
        if (showL && b !== '-' && b !== '.') { g.fillStyle = tc; g.fillText(b, x + cw / 2, y + ch / 2); }
      }
    }
    st.tracks.forEach(function (t) { g.globalAlpha = 0.13; g.fillStyle = t.color; g.fillRect(NW + (t.start - 1) * cw - st.x0, top, (t.end - t.start + 1) * cw, seqH); g.globalAlpha = 1; });
    g.restore();
    // names
    g.save(); g.beginPath(); g.rect(0, top, NW, seqH); g.clip(); g.fillStyle = '#fff'; g.fillRect(0, top, NW, seqH);
    g.fillStyle = '#333'; g.font = '11px sans-serif'; g.textAlign = 'left';
    for (i = r0; i <= r1; i++) g.fillText(trunc(st.d.names[i], 23), 6, top + i * ch - st.y0 + ch / 2);
    g.restore();
    // ruler + tracks (move with the columns only)
    g.save(); g.beginPath(); g.rect(NW, 0, W - NW, top); g.clip();
    g.fillStyle = '#555'; g.font = '10px sans-serif'; g.textAlign = 'center'; g.textBaseline = 'alphabetic';
    var steps = [1, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000], step = 5000; for (i = 0; i < steps.length; i++) if (steps[i] * cw >= 42) { step = steps[i]; break; }
    for (j = Math.ceil((c0 + 1) / step) * step; j <= c1 + 1; j += step) { x = NW + (j - 0.5) * cw - st.x0; g.fillRect(x, RULER - 6, 1, 5); g.fillText(j, x, RULER - 8); }
    g.textBaseline = 'middle';
    st.tracks.forEach(function (t, k) {
      var x1 = NW + (t.start - 1) * cw - st.x0, w = (t.end - t.start + 1) * cw, yy = RULER + k * TRK;
      g.fillStyle = t.color; g.fillRect(x1, yy, w, TRK - 4);
      g.font = 'bold 11px sans-serif'; var tw = g.measureText(t.label).width;
      if (tw + 6 < w) { g.fillStyle = '#fff'; g.textAlign = 'center'; g.fillText(t.label, x1 + w / 2, yy + (TRK - 4) / 2); }
      else { g.fillStyle = t.color; g.textAlign = 'left'; g.fillText(t.label, x1 + w + 4, yy + (TRK - 4) / 2); }
    });
    g.restore();
    g.fillStyle = '#555'; g.font = 'bold 11px sans-serif'; g.textAlign = 'left'; g.textBaseline = 'middle';
    st.tracks.forEach(function (t, k) { g.fillStyle = t.color; g.fillText(t.name, 6, RULER + k * TRK + (TRK - 4) / 2); });
    // consensus + conservation
    g.save(); g.beginPath(); g.rect(NW, botY, W - NW, BOT); g.clip();
    g.font = 'bold ' + Math.min(13, Math.floor(cw * 0.85)) + 'px monospace'; g.textAlign = 'center';
    for (j = c0; j <= c1; j++) {
      x = NW + j * cw - st.x0;
      g.fillStyle = '#fff6c8'; g.fillRect(x, botY + 2, wd(), 17);
      if (showL) { g.fillStyle = '#222'; g.fillText(st.cons[j], x + cw / 2, botY + 11); }
      var h = st.freq[j] * 34; g.fillStyle = st.freq[j] >= 0.9 ? '#7a9e7e' : (st.freq[j] >= 0.6 ? '#e0b34a' : '#e07b73'); g.fillRect(x, H - 4 - h, wd(), h);
    }
    g.restore();
    g.fillStyle = '#333'; g.font = '11px sans-serif'; g.textAlign = 'left'; g.fillText('Consensus', 6, botY + 11); g.fillText('Conservation', 6, botY + 38);
    // frames
    g.strokeStyle = '#d0d4d9'; g.strokeRect(0.5, 0.5, W - 1, H - 1);
    st.info.textContent = st.hover || (n + ' sequences x ' + L + ' columns | zoom ' + cw.toFixed(1) + ' px/column');
  }
  function fit(st) { var W = st.box.clientWidth; st.cw = clamp((W - NW) / st.L, 0.4, 40); st.x0 = 0; draw(st); }
  function focus(st, a, b) {
    var W = st.box.clientWidth; if (!W) { st.pending = [a, b]; return; }
    var span = b - a + 1; st.cw = clamp((W - NW) / (span + 20), 3, 30); st.x0 = Math.max(0, (a - 11) * st.cw); st.pending = null; draw(st);
  }
  function get(id) {
    if (S[id]) return S[id];
    var box = document.getElementById(id); if (!box) return null;
    var st = S[id] = { box: box, cw: 10, x0: 0, y0: 0, tracks: [], minfreq: 0.1, hover: '' };
    box.innerHTML = '';
    var bar = document.createElement('div'); bar.style.cssText = 'display:flex;gap:6px;align-items:center;margin-bottom:6px;flex-wrap:wrap';
    function btn(t, f, title) { var x = document.createElement('button'); x.type = 'button'; x.className = 'btn btn-outline-secondary btn-sm'; x.textContent = t; x.title = title || ''; x.onclick = f; bar.appendChild(x); }
    function zoom(f) { var W = st.box.clientWidth, mid = (st.x0 + (W - NW) / 2) / st.cw; st.cw = clamp(st.cw * f, 0.4, 40); st.x0 = mid * st.cw - (W - NW) / 2; draw(st); }
    btn('+', function () { zoom(1.4); }, 'Zoom in'); btn('−', function () { zoom(1 / 1.4); }, 'Zoom out'); btn('Fit all', function () { fit(st); }, 'Show the whole alignment');
    btn('Go to primers', function () { if (st.tracks.length) { var a = Math.min.apply(null, st.tracks.map(function (t) { return t.start; })), b = Math.max.apply(null, st.tracks.map(function (t) { return t.end; })); focus(st, a, b); } }, 'Centre on the primers and the probe');
    var info = document.createElement('span'); info.style.cssText = 'font-size:12px;color:#555;margin-left:8px'; bar.appendChild(info); st.info = info;
    box.appendChild(bar);
    var cv = document.createElement('canvas'); cv.style.cssText = 'display:block;width:100%;cursor:grab;border-radius:4px;touch-action:none'; box.appendChild(cv); st.cv = cv;
    var drag = null;
    cv.addEventListener('mousedown', function (e) { drag = { x: e.clientX, y: e.clientY, x0: st.x0, y0: st.y0 }; cv.style.cursor = 'grabbing'; });
    window.addEventListener('mouseup', function () { drag = null; cv.style.cursor = 'grab'; });
    window.addEventListener('mousemove', function (e) { if (drag) { st.x0 = drag.x0 - (e.clientX - drag.x); st.y0 = drag.y0 - (e.clientY - drag.y); draw(st); } });
    cv.addEventListener('mousemove', function (e) {
      if (drag || !st.d) return; var r = cv.getBoundingClientRect(), mx = e.clientX - r.left, my = e.clientY - r.top, o = layout(st);
      var h = '';
      if (mx > NW) { var j = Math.floor((mx - NW + st.x0) / st.cw); if (j >= 0 && j < st.L) {
        if (my >= o.top && my < o.botY) { var i = Math.floor((my - o.top + st.y0) / CH); if (i >= 0 && i < st.d.seqs.length) h = st.d.names[i] + ' | column ' + (j + 1) + ' | ' + st.d.seqs[i][j]; }
        else if (my >= o.botY) h = 'column ' + (j + 1) + ' | consensus ' + st.cons[j] + ' | ' + st.pc[j];
      } }
      if (h !== st.hover) { st.hover = h; st.info.textContent = h || (st.d.seqs.length + ' sequences x ' + st.L + ' columns | zoom ' + st.cw.toFixed(1) + ' px/column'); }
    });
    cv.addEventListener('mouseleave', function () { st.hover = ''; if (st.d) draw(st); });
    cv.addEventListener('wheel', function (e) {
      if (!st.d) return; e.preventDefault();
      if (e.shiftKey) { st.x0 += (e.deltaY || e.deltaX); draw(st); return; }
      if (e.altKey) { st.y0 += e.deltaY; draw(st); return; }
      var r = cv.getBoundingClientRect(), mx = clamp(e.clientX - r.left, NW, 1e6), f = e.deltaY < 0 ? 1.18 : 1 / 1.18, col = (mx - NW + st.x0) / st.cw;
      st.cw = clamp(st.cw * f, 0.4, 40); st.x0 = col * st.cw - (mx - NW); draw(st);
    }, { passive: false });
    if (window.ResizeObserver) new ResizeObserver(function () { if (st.pending) focus(st, st.pending[0], st.pending[1]); else draw(st); }).observe(box);
    return st;
  }
  function reg() {
    Shiny.addCustomMessageHandler('prMsaSet', function (m) {
      var st = get(m.id); if (!st) return;
      if (m.seqs) { st.d = { names: [].concat(m.names), seqs: [].concat(m.seqs) }; st.L = st.d.seqs[0].length; st.cw = 10; st.x0 = 0; st.y0 = 0; stats(st); }
      if (m.minfreq && m.minfreq !== st.minfreq) { st.minfreq = m.minfreq; if (st.d) stats(st); }
      if (m.tracks !== undefined) st.tracks = [].concat(m.tracks || []);
      if (m.focus) focus(st, m.focus[0], m.focus[1]); else if (m.seqs && st.box.clientWidth === 0) { /* drawn when shown */ } else draw(st);
    });
  }
  if (window.Shiny && Shiny.addCustomMessageHandler) reg(); else if (window.jQuery) jQuery(document).on('shiny:connected', reg);
})();

)-----"
.pr_card <- function(title, ...) div(class = "card mb-3", div(class = "card-header fw-semibold", title), div(class = "card-body", ...))
.pr_nz <- function(x, d) if (is.null(x) || length(x) == 0L || is.na(x[1])) d else x[1]
.pr_aln_vec <- function(x) {
  if (is.null(x)) return(NULL)
  y <- if (is.character(x)) x else try(as.character(x), silent = TRUE)
  if (inherits(y, "try-error") || !length(y)) return(NULL)
  if (is.null(names(y))) names(y) <- paste0("seq", seq_along(y))
  y
}
.pr_consensus <- function(aln) {
  ctx <- pr_aln_ctx(aln); keep <- ctx$gapfrac <= 0.5 & ctx$tot > 0
  paste(ctx$major[keep], collapse = "")
}
.pr_read_txt <- function(txt, file) {
  x <- character()
  if (!is.null(file) && nrow(file) > 0L) x <- c(x, readLines(file$datapath[1], warn = FALSE))
  c(x, if (!is.null(txt) && nzchar(txt)) txt)
}
# sequence of the template with the primers highlighted; the numbers are positions on the record
.pr_seq_html <- function(seq, win, a, b, marks, width = 50L, block = 10L, junc = NULL) {
  chars <- strsplit(seq, "")[[1]]; cls <- character(length(chars))
  for (k in seq_len(nrow(marks))) cls[marks$a[k]:marks$b[k]] <- marks$cls[k]
  if (length(junc)) { j <- junc[junc >= 1L & junc <= length(chars)]; cls[j] <- trimws(paste(cls[j], "pr-jb")) }     # last base of an exon
  starts <- seq(a, b, by = width); out <- character(length(starts))
  for (u in seq_along(starts)) {
    s0 <- starts[u]; e0 <- min(b, s0 + width - 1L)
    ix_all <- s0:e0
    blocks <- split(ix_all, (seq_along(ix_all) - 1L) %/% block)
    bl <- vapply(blocks, function(ix) {
      r <- rle(cls[ix]); pos <- 0L
      paste(vapply(seq_along(r$lengths), function(k) {
        txt <- paste(chars[ix[pos + seq_len(r$lengths[k])]], collapse = ""); pos <<- pos + r$lengths[k]
        if (nzchar(r$values[k])) sprintf('<span class="%s">%s</span>', r$values[k], txt) else txt
      }, ""), collapse = "")
    }, "")
    out[u] <- paste0('<span class="pr-pos">', formatC(as.numeric(pr_rec_pos(win, s0)), width = 9L, format = "d", big.mark = ""), '</span>  ', paste(bl, collapse = " "))
  }
  tags$pre(class = "pr-seq", HTML(paste(out, collapse = "\n")))
}
.pr_style_verdict <- function(dt, col = "Verdict")
  DT::formatStyle(dt, col, backgroundColor = DT::styleEqual(c("OK", "CHECK", "FAIL"), c("#d4edda", "#fff3cd", "#f8d7da")))

# ---- Primer3Plus-style report blocks ------------------------------------------------------------------------------------
.pr_metric <- function(label, value) tags$span(class = "p3-m", tags$u(label), ": ", value)
.pr_p3_head <- function(label, cls, seq) tags$div(class = paste("p3-head", cls), tags$b(paste0(label, ":")), tags$span(class = "p3-seq", seq))
.pr_tmv <- function(st) if (is.null(st)) 0 else max(0, st$tm)
.pr_struct_block <- function(title, st) {
  tags$div(class = "p3-st", tags$div(class = "p3-st-title", paste0(title, ":")),
    if (is.null(st)) tags$p(class = "text-muted", "No stable structure found.")
    else if (st$tm < 0) tags$p(class = "text-muted", sprintf("No stable structure at the primer concentration (best one: dG %.0f cal/mol, Tm %.0f C).", st$dg * 1000, st$tm))
    else tags$pre(class = "p3-pre", sprintf("Tm: %.1f C   dG: %.0f cal/mol   dH: %.0f cal/mol   dS: %.1f cal/mol*K\n\n%s",
                                            st$tm, st$dg * 1000, st$dh * 1000, st$ds, st$art)))
}
.pr_struct_set <- function(st) tagList(
  .pr_struct_block("Self Any Secondary Structure", st$sd), .pr_struct_block("Self End Secondary Structure", st$sd3),
  .pr_struct_block("Hairpin Secondary Structure", st$hp))
.pr_p3_metrics <- function(st, seq, start = NULL, len = nchar(seq), tm = NULL) tags$div(class = "p3-met",
  if (!is.null(start)) .pr_metric("Start", format(start, big.mark = ",")), .pr_metric("Length", paste(len, "bp")),
  if (!is.null(tm)) .pr_metric("Tm", sprintf("%.1f C", tm)), .pr_metric("GC", sprintf("%.1f %%", .pr_gc(seq))),
  .pr_metric("Any", sprintf("%.1f", .pr_tmv(st$sd))), .pr_metric("End", sprintf("%.1f", .pr_tmv(st$sd3))),
  .pr_metric("HP", sprintf("%.1f", .pr_tmv(st$hp))), .pr_metric("3' Stab", sprintf("%.1f", .pr_stab3(seq))))

primerUI <- function(id) {
  ns <- NS(id)
  cond <- function(expr, ...) conditionalPanel(sprintf(expr, ns("src")), ...)
  condm <- function(expr, ...) conditionalPanel(sprintf(expr, ns("mode")), ...)
  tagList(
    tags$script(HTML(.pr_msa_js)),
    tags$style(HTML(".pr-seq{font-family:Consolas,Menlo,monospace;font-size:13px;word-break:break-all}
                     .pr-box{background:#f6f8fa;border-radius:8px;padding:8px 12px;margin-bottom:6px}
                     .pr-small .form-group{margin-bottom:6px}
                     .pr-f{background:#ccccff;font-weight:700}.pr-r{background:#ffee55;font-weight:700}.pr-p{background:#b9f0c8;font-weight:700}
                     .pr-pos{color:#6b7785}pre.pr-seq{background:#fffff4;border:1px solid #d8d8b0;border-radius:4px;padding:10px 14px;font-size:14px;line-height:1.7;white-space:pre;overflow-x:auto;word-break:normal}.pr-pos{margin-right:2px}
                     .p3-head{padding:6px 10px;margin-top:8px;font-size:15px}.p3-head .p3-seq{font-family:Consolas,Menlo,monospace;background:#fff;border:1px solid #bbb;border-radius:3px;padding:2px 8px;margin-left:8px;letter-spacing:.5px}
                     .pr-jb{border-right:3px solid #8e44ad;padding-right:1px;margin-right:1px}
                     .p3-l{background:#ccccff}.p3-r{background:#ffee55}.p3-pair{background:#d0d0d0}
                     .p3-met{padding:4px 10px;font-size:13px}.p3-m{margin-right:16px;white-space:nowrap}
                     .p3-st{margin:6px 10px}.p3-st-title{color:#1a0dab;text-decoration:underline;font-size:13px;margin-top:8px}
                     pre.p3-pre{background:#fff;border:0;font-family:Consolas,Menlo,monospace;font-size:12.5px;line-height:1.25;margin:2px 0 6px 12px;padding:0}")),
    div(class = "container-fluid",
      h3(icon("dna"), " ", span(class = "i18n", `data-en` = "Primer mining and design", `data-id` = "Penambangan dan desain primer", "Primer mining and design")),
      p(class = "text-muted", "Find sequences (NCBI) and primers (BOLD, PR2, PrimerBank, qPrimerDB), design primers for PCR, sequencing, exon junctions or degenerate primers, check their quality and BLAST them."),
      tabsetPanel(id = ns("pt"), type = "tabs",
        # ------------------------------------------------------------------------------------------------ 1: sequences + primers
        tabPanel("Find sequences and primers", value = "mine", br(),
          h4(icon("magnifying-glass"), " Find sequences (NCBI)"),
          fluidRow(
            column(4,
              .pr_card("NCBI nucleotide (GenBank)",
                textInput(ns("ncbi_term"), "NCBI search", placeholder = "Danaus plexippus[Organism] AND COI[Gene]"),
                numericInput(ns("ncbi_max"), "Max sequences", 200, 10, 2000),
                helpText("Entrez syntax. The contact email of the BLAST tab is sent to NCBI if you filled it. Your own FASTA / GenBank / alignment files are loaded in the Design tab."),
                checkboxInput(ns("append"), "Add to the sequences already listed", FALSE),
                actionButton(ns("go_fetch"), "Search NCBI", icon = icon("magnifying-glass"), class = "btn-primary"))),
            column(8,
              .pr_card("Sequences",
                uiOutput(ns("pool_info")),
                DT::DTOutput(ns("pool_tbl")),
                div(class = "mt-2",
                  actionButton(ns("to_design"), "Send selected to Design", icon = icon("arrow-right"), class = "btn-success"), " ",
                  actionButton(ns("gb_get"), "Get GenBank annotation", icon = icon("download")), " ",
                  actionButton(ns("pool_rm"), "Remove selected"), " ",
                  actionButton(ns("pool_clear"), "Clear"), " ",
                  downloadButton(ns("pool_dl"), "FASTA"))),
              .pr_card("Annotation (GenBank): choose the target gene",
                uiOutput(ns("gb_info")),
                uiOutput(ns("gb_viz")),
                selectInput(ns("gb_type"), NULL, c("Genes and CDS" = "gene_cds", "mRNA / exons" = "mrna", "All features" = "all"), width = "260px"),
                DT::DTOutput(ns("gb_tbl")),
                div(class = "mt-2",
                  actionButton(ns("gb_use"), "Use selected gene in Design", icon = icon("bullseye"), class = "btn-success"), " ",
                  actionButton(ns("gb_splice"), "Design exon-exon junction primers for it", icon = icon("scissors")), " ",
                  downloadButton(ns("gb_dl"), "GenBank file"))))),
          hr(),
          h4(icon("flask"), " Find primers (BOLD and other databases)"),
          fluidRow(
            column(4,
              .pr_card("BOLD primer database",
                uiOutput(ns("bp_info")),
                helpText("Type a word in the filter above the table on the right (insect, plant, COI, 16S, rbcL, a submitter name ...). To use a newer export of ",
                         tags$a(href = "https://portal.boldsystems.org/primer", target = "_blank", "portal.boldsystems.org/primer"),
                         " (its FASTA button), load that file here; it replaces the built-in table."),
                fileInput(ns("bp_file"), NULL, buttonLabel = "BOLD FASTA ...", placeholder = "FASTA / CSV from the BOLD primer page")),
              .pr_card("BOLD sequence records (optional)",
                tags$details(tags$summary("Search records: the primer names written in them are listed, and the records can be used as sequences"),
                  textInput(ns("bold_taxon"), "Taxon", placeholder = "e.g. Danaus plexippus, Coleoptera, Zingiberaceae"),
                  selectizeInput(ns("bold_marker"), "Marker (empty = any)",
                                 choices = c("COI-5P", "COI-3P", "rbcLa", "matK", "ITS", "16S", "18S", "12S", "CYTB", "trnH-psbA"),
                                 selected = "COI-5P", multiple = TRUE, options = list(create = TRUE)),
                  fluidRow(column(6, numericInput(ns("bold_max"), "Max records", 300, 10, 5000)),
                           column(6, numericInput(ns("bold_per_sp"), "Max per species (0 = all)", 5, 0, 100))),
                  actionButton(ns("bold_go"), "Search BOLD records", icon = icon("magnifying-glass")), " ",
                  actionButton(ns("bold_add"), "Add records to the sequence list", icon = icon("plus")))),
              .pr_card("Load primers from another database",
                helpText("PR2-primers, PrimerBank and qPrimerDB have no open interface that SangeRly can query. Open the site, then bring the primers here: ",
                         "load the table you downloaded (CSV / TSV / XLSX: the sequence columns are recognised automatically, also primer sets with a forward and a reverse column) ",
                         "or paste primers copied from a result page (one per line: name  sequence)."),
                fileInput(ns("pdb_file"), NULL, buttonLabel = "Table ...", placeholder = "CSV / TSV / XLSX"),
                textAreaInput(ns("pdb_txt"), NULL, rows = 4, placeholder = "or paste:  name  sequence  (one per line)"),
                actionButton(ns("pdb_load"), "Load primers", icon = icon("table"), class = "btn-primary")),
              .pr_card("Primer databases on the web",
                tags$ul(class = "pr-small",
                  tags$li(tags$a(href = "https://app.pr2-primers.org/pr2-primers/", target = "_blank", "PR2-primers"), ": 18S / eukaryote rRNA primers and primer sets for metabarcoding (download the Primers or Primer sets table)."),
                  tags$li(tags$a(href = "https://pga.mgh.harvard.edu/primerbank/", target = "_blank", "PrimerBank"), ": human and mouse qPCR primers for gene expression (search by gene, copy the sequences)."),
                  tags$li(tags$a(href = "https://qprimerdb.biodb.org", target = "_blank", "qPrimerDB"), ": gene-specific qPCR primers for more than 1,000 organisms (search by gene, copy the sequences)."),
                  tags$li("In this app: the reference list on the right, and primers named inside BOLD records.")))),
            column(8,
              .pr_card("Primers (BOLD database and published primers)",
                fluidRow(column(5, textInput(ns("lib_q"), NULL, placeholder = "Filter: COI, fish, ITS, plants ...")),
                         column(7, checkboxGroupInput(ns("lib_src"), NULL, c("BOLD", "PR2", "Published"), selected = c("BOLD", "PR2", "Published"), inline = TRUE))),
                uiOutput(ns("lib_info")),
                DT::DTOutput(ns("lib_tbl")),
                div(class = "mt-2", actionButton(ns("lib_send"), "Send selected to Quality", icon = icon("arrow-right")))),
              .pr_card("Primers loaded from a database table",
                DT::DTOutput(ns("pdb_tbl")),
                div(class = "mt-2", actionButton(ns("pdb_send"), "Send selected to Quality", icon = icon("arrow-right"))))))),
        # ------------------------------------------------------------------------------------------------ 2
        tabPanel("Design", value = "design", br(),
          fluidRow(
            column(4,
              .pr_card("What to design",
                radioButtons(ns("mode"), NULL, c("Conventional PCR (primer pair)" = "pcr", "Sequencing primers (Sanger, with walking)" = "seq",
                                                  "mRNA exon-exon junction" = "junction", "Degenerate primers (from an alignment)" = "degen"))),
              .pr_card("Template",
                conditionalPanel(sprintf("(['pcr','seq'].indexOf(input['%s']) >= 0 && input['%s'] == 'gb') || (input['%s'] == 'junction' && input['%s'] == 'gb')", ns("mode"), ns("tpl_src"), ns("mode"), ns("tx_src")),
                  fileInput(ns("gb_file"), NULL, buttonLabel = "GenBank file ...", placeholder = "load your own .gb / .gbk (or use a record of step 1)")),
                conditionalPanel(sprintf("input['%s'] == 'pcr' || input['%s'] == 'seq'", ns("mode"), ns("mode")),
                  radioButtons(ns("tpl_src"), NULL, c("Paste or upload" = "paste", "A sequence of step 1" = "pool",
                                                       "Annotated GenBank record (step 1)" = "gb",
                                                       "Consensus of the alignment" = "cons", "Assembly consensus" = "asm")),
                  conditionalPanel(sprintf("input['%s'] == 'paste'", ns("tpl_src")),
                    textAreaInput(ns("tpl_txt"), NULL, rows = 5, placeholder = "FASTA or plain sequence (cDNA / mRNA for junctions)"),
                    fileInput(ns("tpl_file"), NULL, buttonLabel = "File ...", placeholder = "no file")),
                  conditionalPanel(sprintf("input['%s'] == 'pool'", ns("tpl_src")), selectInput(ns("tpl_pick"), "Sequence", choices = NULL)),
                  conditionalPanel(sprintf("input['%s'] == 'gb'", ns("tpl_src")),
                    selectizeInput(ns("gb_feat"), "Target gene / feature", choices = NULL, options = list(placeholder = "type to search, e.g. rbcL")),
                    fluidRow(column(6, numericInput(ns("gb_flank"), "Extra bases around it", 300, 0, 5000)),
                             column(6, checkboxInput(ns("gb_orient"), "Read it in the gene direction", TRUE))),
                    helpText("Or type a region below (positions of the GenBank record). Positions in the results are positions on that record.")),
                  uiOutput(ns("tpl_info")),
                  fluidRow(column(6, numericInput(ns("tgt_from"), "Region from", NA, 1)), column(6, numericInput(ns("tgt_to"), "to", NA, 1))),
                  conditionalPanel(sprintf("input['%s'] == 'pcr'", ns("mode")),
                    radioButtons(ns("tgt_role"), NULL, c("Primers flank the region (amplify all of it)" = "flank", "Primers inside the region (amplify part of it)" = "inside"))),
                  helpText("The region that must be inside the product (PCR) or read (sequencing). Empty = the whole sequence.")),
                condm("input['%s'] == 'junction'",
                  radioButtons(ns("tx_src"), NULL, c("Transcript of the GenBank record from step 1 (exons are read from the annotation)" = "gb",
                                                      "Paste or upload an mRNA / cDNA sequence" = "paste")),
                  conditionalPanel(sprintf("input['%s'] == 'gb'", ns("tx_src")), uiOutput(ns("tx_pick_ui"))),
                  conditionalPanel(sprintf("input['%s'] == 'paste'", ns("tx_src")),
                    textAreaInput(ns("tx_txt"), NULL, rows = 5, placeholder = "mRNA / cDNA sequence (FASTA or plain)"),
                    fileInput(ns("tx_file"), NULL, buttonLabel = "File ...", placeholder = "no file"),
                    radioButtons(ns("jn_type"), "Where are the exon boundaries?", c("Positions of the last base of each exon" = "pos", "Exon lengths" = "len",
                                                                                   "Each FASTA record is one exon" = "fasta")),
                    conditionalPanel(sprintf("input['%s'] != 'fasta'", ns("jn_type")), textInput(ns("jn_txt"), NULL, placeholder = "e.g. 120, 340, 512")),
                    helpText("Only for a sequence without annotation. With a GenBank record the exons are found automatically."))),
                condm("input['%s'] == 'degen'",
                  selectInput(ns("aln_src"), "Alignment", c("Sequences of step 1 (aligned here)" = "pool", "MSA tab: kept columns (trimmed)" = "trim",
                                                             "MSA tab: whole alignment" = "all", "Paste or upload FASTA (aligned here if needed)" = "paste")),
                  conditionalPanel(sprintf("input['%s'] == 'paste'", ns("aln_src")), textAreaInput(ns("aln_txt"), NULL, rows = 5),
                                   fileInput(ns("aln_file"), NULL, buttonLabel = "FASTA file ...", placeholder = "no file")),
                  uiOutput(ns("aln_info")),
                  fluidRow(column(6, numericInput(ns("dg_from"), "Columns from", NA, 1)), column(6, numericInput(ns("dg_to"), "to", NA, 1))))),
              .pr_card("Settings",
                div(class = "pr-small",
                  numericInput(ns("n_out"), "Max number of designs (primer pairs)", 5, 1, 50, step = 1),
                  tags$details(open = NA, tags$summary(strong("Primer")),
                    fluidRow(column(4, numericInput(ns("len_min"), "Length min", 18, 10, 40)), column(4, numericInput(ns("len_opt"), "opt", 20, 10, 40)),
                             column(4, numericInput(ns("len_max"), "max", 25, 10, 45))),
                    fluidRow(column(4, numericInput(ns("tm_min"), "Tm min", 57)), column(4, numericInput(ns("tm_opt"), "opt", 60)),
                             column(4, numericInput(ns("tm_max"), "max", 63))),
                    fluidRow(column(6, numericInput(ns("gc_min"), "GC % min", 40, 0, 100)), column(6, numericInput(ns("gc_max"), "GC % max", 60, 0, 100))),
                    fluidRow(column(6, numericInput(ns("prod_min"), "Product min (bp)", 100, 20)), column(6, numericInput(ns("prod_max"), "Product max (bp)", 700, 20))),
                    numericInput(ns("max_dtm"), "Max Tm difference F / R", 3, 0, 20)),
                  tags$details(tags$summary(strong("Reaction (for Tm)")),
                    fluidRow(column(6, numericInput(ns("na"), "Na+ / K+ (mM)", 50)), column(6, numericInput(ns("mg"), "Mg2+ (mM)", 1.5))),
                    fluidRow(column(6, numericInput(ns("dntp"), "dNTP (mM)", 0.2)), column(6, numericInput(ns("conc"), "Primer (nM)", 200)))),
                  condm("input['%s'] == 'seq'",
                    tags$details(open = NA, tags$summary(strong("Sequencing")),
                      numericInput(ns("sq_offset"), "Bases lost after the primer", 50, 0, 200),
                      numericInput(ns("sq_read"), "Read length (nt)", 700, 200, 1200),
                      numericInput(ns("sq_overlap"), "Overlap between reads (nt)", 100, 0, 400),
                      checkboxInput(ns("sq_both"), "Also reverse primers", TRUE))),
                  condm("input['%s'] == 'junction'",
                    tags$details(open = NA, tags$summary(strong("Exon junctions")),
                      selectInput(ns("jn_span"), "Where the junction is", c("Either primer spans a junction" = "either", "Forward primer spans" = "forward",
                                                                            "Reverse primer spans" = "reverse", "Primers in two different exons (no span)" = "flank")),
                      fluidRow(column(6, numericInput(ns("jn_ov5"), "Bases on 5' side", 7, 3, 15)), column(6, numericInput(ns("jn_ov3"), "Bases on 3' side", 4, 2, 15))))),
                  condm("input['%s'] == 'degen'",
                    tags$details(open = NA, tags$summary(strong("Degeneracy")),
                      fluidRow(column(6, numericInput(ns("dg_max"), "Max degeneracy", 64, 1, 4096)), column(6, numericInput(ns("dg_minfreq"), "Min base frequency", 0.1, 0.01, 0.5, 0.01))),
                      fluidRow(column(6, numericInput(ns("dg_gap"), "Max gap fraction", 0.2, 0, 0.9, 0.05)), column(6, numericInput(ns("dg_mm"), "Mismatches allowed", 0, 0, 3))),
                      numericInput(ns("dg_cov"), "Min coverage of the sequences", 0.8, 0.3, 1, 0.05))),
                  condm("input['%s'] != 'seq'",
                    checkboxInput(ns("probe"), "Also design a probe (hydrolysis / qPCR)", FALSE),
                    conditionalPanel(sprintf("input['%s']", ns("probe")), numericInput(ns("probe_delta"), "Probe Tm above the primers", 8, 3, 15)))),
                actionButton(ns("design_go"), "Design primers", icon = icon("wand-magic-sparkles"), class = "btn-primary"))),
            column(8,
              conditionalPanel(sprintf("input['%s'] == 'junction'", ns("mode")),
                .pr_card("Transcript and its exons", plotOutput(ns("tx_map"), height = "230px"), DT::DTOutput(ns("tx_tbl")))),
              conditionalPanel(sprintf("input['%s'] == 'degen'", ns("mode")),
                .pr_card("Alignment view: consensus, variation and conservation",
                  uiOutput(ns("msa_ctl")),
                  div(id = ns("msa_box"), style = "min-height:60px"))),
              .pr_card("Result",
                uiOutput(ns("design_msg")),
                DT::DTOutput(ns("design_tbl")),
                div(class = "mt-2",
                  actionButton(ns("send_ana"), "Send selected to Quality", icon = icon("arrow-right")), " ",
                  downloadButton(ns("design_dl"), "CSV"))),
              .pr_card("Primer positions", uiOutput(ns("design_detail")), uiOutput(ns("map_ui"))),
              .pr_card("Template sequence", uiOutput(ns("seq_view")))))),
        # ------------------------------------------------------------------------------------------------ 3
        tabPanel("Quality", value = "ana", br(),
          fluidRow(
            column(4,
              .pr_card("Primers",
                textAreaInput(ns("prim_txt"), NULL, rows = 9, placeholder = "One per line:  name  sequence   (or just the sequence, or FASTA)\nDegenerate letters (R Y S W K M B D H V N) are accepted."),
                helpText("Tm uses the reaction settings of the Design tab."),
                checkboxInput(ns("bind_on"), "Check binding in a set of sequences", FALSE),
                conditionalPanel(sprintf("input['%s']", ns("bind_on")),
                  selectInput(ns("bind_src"), NULL, c("Sequences of step 1" = "pool", "Alignment (step 1 / MSA tab)" = "aln", "Paste FASTA" = "paste")),
                  conditionalPanel(sprintf("input['%s'] == 'paste'", ns("bind_src")), textAreaInput(ns("bind_txt"), NULL, rows = 4)),
                  numericInput(ns("bind_max"), "Max sequences scanned", 300, 10, 2000)),
                actionButton(ns("ana_go"), "Analyse", icon = icon("flask"), class = "btn-primary"))),
            column(8,
              .pr_card("Primer quality", uiOutput(ns("ana_msg")), DT::DTOutput(ns("ana_tbl")),
                       helpText("dG values are estimates (kcal/mol, 37 degC, ungapped duplexes with up to 2 internal mismatches); more negative = stronger. Tm = melting temperature of that structure. ",
                                "Degeneracy = number of different sequences the primer stands for: 1 means one single, unambiguous sequence (see 'Ambiguous bases' = 0). ",
                                "CHECK is raised for: no G/C in the last 2 nt (no GC clamp), a run of 4 or more identical bases, 3 or more 3' bases pairing with the primer itself, Tm / GC / hairpin / dimer limits. Click a row to see the structures. ",
                                "A sequence whose name contains 'probe' (the probes sent from Design are named PairN_probe) is judged as a probe, not as a primer: length 18-30 nt, GC 30-80 %, no G at the 5' end (it quenches the dye), Tm about 6-10 C above the primers, no run of G; the GC clamp and 3' end rules are not applied to it.")),
              .pr_card("Pairs (hetero-dimers)", DT::DTOutput(ns("pair_tbl"))),
              .pr_card("Primer report (Primer3Plus style): primers, hairpins, self- and cross dimers",
                       selectInput(ns("rep_sel"), "Show", choices = character(), width = "460px"), uiOutput(ns("qrep"))),
              conditionalPanel(sprintf("input['%s']", ns("bind_on")),
                .pr_card("Binding in the sequence set", DT::DTOutput(ns("bind_tbl")), br(), DT::DTOutput(ns("bind_detail"))))))),
        # ------------------------------------------------------------------------------------------------ 4
        tabPanel("BLAST", value = "blast", br(),
          fluidRow(
            column(5,
              .pr_card("Primer parameters (as in NCBI Primer-BLAST)",
                helpText("Choose primers from the Quality list or type / paste them; edit any base (for a degenerate primer replace the ambiguous letter, or pick one variant) before sending them to NCBI."),
                selectizeInput(ns("bl_fpick"), "Forward primer from the list", choices = NULL, options = list(placeholder = "optional")),
                tags$label("Forward primer (5'->3' on the plus strand)"),
                div(style = "display:flex; gap:6px; align-items:flex-start;", div(style = "flex:1", textInput(ns("bl_fseq"), NULL, width = "100%")),
                    actionButton(ns("bl_fclear"), "Clear", class = "btn-outline-primary btn-sm")),
                uiOutput(ns("bl_fvar")),
                selectizeInput(ns("bl_rpick"), "Reverse primer from the list", choices = NULL, options = list(placeholder = "optional")),
                tags$label("Reverse primer (5'->3' on the minus strand)"),
                div(style = "display:flex; gap:6px; align-items:flex-start;", div(style = "flex:1", textInput(ns("bl_rseq"), NULL, width = "100%")),
                    actionButton(ns("bl_rclear"), "Clear", class = "btn-outline-primary btn-sm")),
                uiOutput(ns("bl_rvar")),
                div(class = "mt-2", uiOutput(ns("bl_pb"))))),
            column(7,
              .pr_card("Primers of the Quality tab (IUPAC)",
                DT::DTOutput(ns("bl_seq")),
                helpText("R = A/G, Y = C/T, S = G/C, W = A/T, K = G/T, M = A/C, B = C/G/T, D = A/G/T, H = A/C/T, V = A/C/G, N = any. ",
                         "The reverse complement is the primer as it appears on the template strand.")),
              .pr_card("How it works",
                tags$ul(class = "pr-small",
                  tags$li("The search runs on the NCBI website (your browser); the app only fills in the primers. Database, organism and product size are chosen there."),
                  tags$li(strong("Primer-BLAST"), ": checks where the pair amplifies and what else it amplifies (specificity)."),
                  tags$li(strong("BLAST"), ": searches each primer separately (short-sequence settings) to see all places it binds."),
                  tags$li("If NCBI refuses a degenerate primer, pick one variant (or replace the ambiguous letter) and send it again.")))))))))
}

primerServer <- function(id, consensus = reactive(NULL), aln_trim = reactive(NULL), aln_all = reactive(NULL), email = reactive("")) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    rv <- reactiveValues(bpdb = NULL, bold = NULL, pool = NULL, aln = NULL, design = NULL, prim = NULL, ana = NULL, pairs = NULL, bind = NULL,
                         hits = NULL, amp = NULL, lib_sel = NULL, gb = NULL, gb_err = NULL, spliced = NULL, tx_sel = NULL, bl = NULL, bl_err = NULL)
    notify <- function(msg, type = "error", dur = 12) showNotification(msg, type = type, duration = dur)
    prog <- function(m) tryCatch(incProgress(0, detail = m), error = function(e) NULL)

    prm <- reactive({
      pr_params(len_min = as.integer(.pr_nz(input$len_min, 18)), len_opt = as.integer(.pr_nz(input$len_opt, 20)), len_max = as.integer(.pr_nz(input$len_max, 25)),
                tm_min = .pr_nz(input$tm_min, 57), tm_opt = .pr_nz(input$tm_opt, 60), tm_max = .pr_nz(input$tm_max, 63),
                gc_min = .pr_nz(input$gc_min, 40), gc_max = .pr_nz(input$gc_max, 60), gc_opt = mean(c(.pr_nz(input$gc_min, 40), .pr_nz(input$gc_max, 60))),
                prod_min = as.integer(.pr_nz(input$prod_min, 100)), prod_max = as.integer(.pr_nz(input$prod_max, 700)),
                max_dtm = .pr_nz(input$max_dtm, 3), na = .pr_nz(input$na, 50), mg = .pr_nz(input$mg, 1.5),
                dntp = .pr_nz(input$dntp, 0.2), conc = .pr_nz(input$conc, 200),
                n_out = max(1L, min(50L, as.integer(.pr_nz(input$n_out, 5)))))
    })

    # ---------------------------------------------------------------------------------------------- 1. sequences
    mk_recs <- function(seqs, source) {
      nm <- unname(names(seqs)); if (is.null(nm)) nm <- paste0("seq", seq_along(seqs))
      seqs <- unname(seqs)
      data.frame(ID = sub("\\s.*$", "", nm), Species = ifelse(grepl("\\s", nm), sub("^\\S+\\s+", "", nm), ""), Marker = "", BIN = "",
                 Accession = "", Primers = "", Sequence = toupper(gsub("[^A-Za-z]", "", seqs)), Length = nchar(gsub("[^A-Za-z]", "", seqs)),
                 Source = source, stringsAsFactors = FALSE)
    }
    observeEvent(input$go_fetch, {
      res <- tryCatch(withProgress(message = "Looking for sequences ...", value = 0.3,
                        pr_ncbi_fetch(input$ncbi_term, as.integer(.pr_nz(input$ncbi_max, 200)), email(), prog)),
                      error = function(e) { notify(conditionMessage(e)); NULL })
      req(!is.null(res))
      rv$pool <- if (isTRUE(input$append) && !is.null(rv$pool)) rbind(rv$pool, res[, names(rv$pool)]) else res
      rv$aln <- NULL
      notify(sprintf("%d sequence(s) listed.", nrow(rv$pool)), "message", 6)
    })
    # BOLD: the records are searched for the primer names they carry; they only enter the sequence list when asked
    output$bp_info <- renderUI({
      b <- rv$bpdb; src <- "loaded from your file"
      if (is.null(b)) { b <- pr_bold_builtin(); src <- "built into the app" }
      if (is.null(b)) return(tags$div(class = "alert alert-warning", "The built-in BOLD table could not be read (the jsonlite package is needed). Load the FASTA of the BOLD primer page below instead."))
      helpText(strong(format(nrow(b), big.mark = ",")), sprintf(" BOLD primers (%s): %d forward, %d reverse, %d markers. Identical name + sequence entries were merged.", src, sum(b$Direction == "F"), sum(b$Direction == "R"), length(unique(b$Marker))))
    })
    observeEvent(input$bp_file, {
      res <- tryCatch(pr_read_primer_table(input$bp_file$datapath[1], input$bp_file$name[1]), error = function(e) { notify(conditionMessage(e), dur = 15); NULL })
      req(!is.null(res)); res$Source <- "BOLD primer database"; rv$bpdb <- res
      notify(sprintf("%d BOLD primer(s) loaded; they are in the table on the right.", nrow(res)), "message", 6)
    })
    observeEvent(input$bold_go, {
      res <- tryCatch(withProgress(message = "Searching BOLD ...", value = 0.3,
                        pr_bold_fetch(input$bold_taxon, input$bold_marker, as.integer(.pr_nz(input$bold_max, 300)), as.integer(.pr_nz(input$bold_per_sp, 0)), prog)),
                      error = function(e) { notify(conditionMessage(e)); NULL })
      req(!is.null(res)); rv$bold <- res
      nb <- pr_bold_primers(res)
      notify(sprintf("%d BOLD record(s); %s.", nrow(res), if (is.null(nb)) "no primer name found in them" else sprintf("%d primer(s) named in them", nrow(nb))), "message", 8)
    })
    observeEvent(input$bold_add, {
      b <- rv$bold; req(!is.null(b))
      rv$pool <- if (!is.null(rv$pool)) rbind(rv$pool, b[, names(rv$pool)]) else b
      rv$aln <- NULL; notify(sprintf("%d BOLD record(s) added to the sequence list.", nrow(b)), "message", 6)
    })
    # a GenBank file of the user (Design tab): annotated record, its genes / exons can be chosen as the target
    observeEvent(input$gb_file, {
      res <- tryCatch({
        txt <- readLines(input$gb_file$datapath[1], warn = FALSE)
        g <- pr_gb_parse(txt); g$raw <- paste(txt, collapse = "\n"); g
      }, error = function(e) { notify(paste("Cannot read the GenBank file:", conditionMessage(e))); NULL })
      req(!is.null(res)); rv$gb <- res
      r <- mk_recs(setNames(res$seq, paste(res$accession, res$definition)), "GenBank file")
      rv$pool <- if (!is.null(rv$pool)) rbind(rv$pool, r[, names(rv$pool)]) else r
      notify(sprintf("GenBank record %s loaded (%d feature(s)).", res$accession, nrow(res$features)), "message", 6)
    })
    output$pool_info <- renderUI({
      p <- rv$pool
      if (is.null(p)) return(helpText("Nothing yet. Search NCBI above (your own files are loaded in the Design tab)."))
      tags$p(sprintf("%d sequence(s); length %d-%d nt (median %d); %d species. %s", nrow(p), min(p$Length), max(p$Length), as.integer(stats::median(p$Length)),
                     length(unique(p$Species[nzchar(p$Species)])),
                     if (!is.null(rv$aln)) sprintf("Aligned: %d sequences x %d columns.", length(rv$aln), nchar(rv$aln[1])) else "Not aligned yet."))
    })
    output$pool_tbl <- DT::renderDT({
      p <- rv$pool; req(p)
      DT::datatable(p[, c("ID", "Species", "Marker", "BIN", "Length", "Source")], rownames = FALSE, selection = "multiple",
                    options = list(pageLength = 8, scrollX = TRUE, dom = "frtip"))
    })
    pool_idx <- reactive({ s <- input$pool_tbl_rows_selected; if (length(s)) s else seq_len(nrow(rv$pool)) })
    observeEvent(input$pool_rm, { s <- input$pool_tbl_rows_selected; req(length(s), rv$pool); rv$pool <- rv$pool[-s, , drop = FALSE]; rv$aln <- NULL })
    observeEvent(input$pool_clear, { rv$pool <- NULL; rv$aln <- NULL })
    do_align <- function() {
      req(rv$pool)
      p <- rv$pool[pool_idx(), , drop = FALSE]
      if (nrow(p) < 2L) { notify("At least 2 sequences are needed to align (step 1).", "warning", 8); return(invisible(NULL)) }
      if (nrow(p) > 400L) { p <- p[seq_len(400L), , drop = FALSE]; notify("Only the first 400 sequences are aligned here.", "warning", 8) }
      s <- setNames(p$Sequence, make.unique(trimws(paste(p$ID, gsub("[^A-Za-z0-9]+", "_", p$Species)))))
      a <- tryCatch(withProgress(message = "Aligning ...", value = 0.5, pr_align(s)), error = function(e) { notify(conditionMessage(e)); NULL })
      if (!is.null(a)) { rv$aln <- a; notify(sprintf("Aligned: %d sequences, %d columns.", length(a), nchar(a[1])), "message", 6) }
    }
    # choosing "Sequences of step 1" aligns them at once, so the alignment is visible below
    observeEvent(input$pool_tbl_rows_selected, { rv$aln <- NULL }, ignoreInit = TRUE, priority = 10)
    observeEvent(list(input$aln_src, input$mode, rv$pool, input$pool_tbl_rows_selected), {
      if (identical(input$mode, "degen") && identical(input$aln_src, "pool") && !is.null(rv$pool) && is.null(rv$aln)) do_align()
    }, ignoreInit = FALSE)
    output$pool_dl <- downloadHandler(filename = function() "sequences.fasta", content = function(file) {
      p <- rv$pool; req(p)
      writeLines(paste0(">", p$ID, ifelse(nzchar(p$Species), paste0(" ", p$Species), ""), "\n", p$Sequence), file)
    })

    lib <- reactive({
      L <- pr_primer_library()
      b <- rv$bpdb; if (is.null(b)) b <- pr_bold_builtin()
      if (!is.null(b)) {
        if (!"Marker" %in% names(b)) b$Marker <- b$Target
        if (!"Submitter" %in% names(b)) b$Submitter <- ""
        L <- rbind(data.frame(Name = b$Name, Sequence = b$Sequence, Direction = b$Direction, Marker = b$Marker, Taxa = "", Reference = ifelse(nzchar(b$Submitter), paste("BOLD,", b$Submitter), "BOLD primer database"), stringsAsFactors = FALSE), L)
      }
      bp <- pr_bold_primers(rv$bold)
      if (!is.null(bp)) {
        bp <- data.frame(Name = bp$Name, Sequence = bp$Sequence, Direction = "", Marker = paste(.pr_nz(input$bold_marker, ""), collapse = ","),
                         Taxa = .pr_nz(input$bold_taxon, ""), Reference = sprintf("named in %d BOLD record(s)", bp$Records), stringsAsFactors = FALSE)
        L <- rbind(bp, L)
      }
      p2 <- pr_pr2_builtin()
      if (!is.null(p2)) L <- rbind(data.frame(Name = p2$Name, Sequence = p2$Sequence, Direction = p2$Direction, Marker = p2$Marker, Taxa = p2$Taxa, Reference = p2$Reference, stringsAsFactors = FALSE), L)
      L$Source <- ifelse(grepl("^BOLD", L$Reference), "BOLD", ifelse(grepl("^PR2", L$Reference), "PR2", ifelse(grepl("named in", L$Reference), "BOLD records", "Published")))
      src <- input$lib_src; if (!is.null(src)) L <- L[L$Source %in% src | (L$Source == "BOLD records"), , drop = FALSE]
      q <- trimws(.pr_nz(input$lib_q, ""))
      if (nzchar(q)) L <- L[grepl(q, apply(L, 1, paste, collapse = " "), ignore.case = TRUE), , drop = FALSE]
      L
    })
    output$lib_info <- renderUI({
      bp <- pr_bold_primers(rv$bold)
      n <- nrow(lib())
      helpText(sprintf("%s primer(s) shown (BOLD primer database, PR2-primers and a reference list of widely used primers). Use 'Show ... entries' to list more per page. ", format(n, big.mark = ",")),
               "Check the original paper before ordering.")
    })
    output$lib_tbl <- DT::renderDT(DT::datatable(lib(), rownames = FALSE, selection = "multiple", options = list(pageLength = 10, lengthMenu = c(5, 10, 25, 50, 100, 500, 1000), scrollX = TRUE, dom = "lrtip")))
    observeEvent(input$lib_send, {
      s <- input$lib_tbl_rows_selected; req(length(s))
      L <- lib()[s, , drop = FALSE]
      add_primers(paste0(L$Name, "\t", L$Sequence))
      updateTabsetPanel(session, "pt", selected = "ana")
    })
    observeEvent(input$pdb_load, {
      res <- tryCatch({
        tabs <- list()
        if (!is.null(input$pdb_file) && nrow(input$pdb_file) > 0L) tabs[[length(tabs) + 1L]] <- pr_read_primer_table(input$pdb_file$datapath[1], input$pdb_file$name[1])
        txt <- input$pdb_txt
        if (!is.null(txt) && nzchar(trimws(txt))) {
          x <- pr_parse_primers(txt)
          if (!is.null(x) && nrow(x)) tabs[[length(tabs) + 1L]] <- data.frame(Name = x$Name, Sequence = x$Sequence, Direction = "", Target = "", Source = "pasted", stringsAsFactors = FALSE)
        }
        if (!length(tabs)) stop("Choose a table or paste primers (name  sequence, one per line).")
        d <- do.call(rbind, tabs); d$Name <- make.unique(d$Name); d
      }, error = function(e) { notify(conditionMessage(e), dur = 15); NULL })
      req(!is.null(res)); rv$pdb <- res
      notify(sprintf("%d primer(s) loaded.", nrow(res)), "message", 6)
    })
    output$pdb_tbl <- DT::renderDT({
      x <- rv$pdb
      if (is.null(x)) return(DT::datatable(data.frame(Note = "Nothing loaded yet."), rownames = FALSE, options = list(dom = "t")))
      DT::datatable(x, rownames = FALSE, selection = "multiple", filter = "top", options = list(pageLength = 8, scrollX = TRUE, dom = "frtip"))
    })
    observeEvent(input$pdb_send, {
      s <- input$pdb_tbl_rows_selected; req(length(s), rv$pdb)
      L <- rv$pdb[s, , drop = FALSE]
      add_primers(paste0(L$Name, "\t", L$Sequence))
      updateTabsetPanel(session, "pt", selected = "ana")
    })
    add_primers <- function(lines) {
      cur <- .pr_nz(input$prim_txt, "")
      updateTextAreaInput(session, "prim_txt", value = paste(c(if (nzchar(trimws(cur))) cur, lines), collapse = "\n"))
    }

    # ---------------------------------------------------------------------------------------------- GenBank annotation
    gb_idx <- reactive({
      g <- rv$gb; req(g); f <- g$features
      keep <- switch(.pr_nz(input$gb_type, "gene_cds"), gene_cds = f$Type %in% c("gene", "CDS", "rRNA", "tRNA", "ncRNA"),
                     mrna = f$Type %in% c("mRNA", "exon", "CDS", "intron"), all = f$Type != "source", f$Type != "source")
      which(keep)
    })
    output$gb_info <- renderUI({
      g <- rv$gb
      if (is.null(g) && !is.null(rv$gb_err)) return(tags$div(class = "alert alert-danger", strong("Annotation not loaded: "), rv$gb_err, " Press 'Get GenBank annotation' again."))
      if (is.null(g)) return(helpText("Step 1: click one sequence in the table above. Step 2: press 'Get GenBank annotation' (NCBI records) or load a GenBank file. Step 3: click a gene in the table that appears here, then 'Use selected gene in Design'."))
      tags$p(strong(g$accession), " - ", g$definition, sprintf(" (%s nt, %d features)", format(nchar(g$seq), big.mark = ","), nrow(g$features)))
    })
    output$gb_tbl <- DT::renderDT({
      g <- rv$gb; req(g); i <- gb_idx()
      DT::datatable(g$features[i, c("Type", "Name", "Product", "Start", "End", "Strand", "Length", "Parts"), drop = FALSE], rownames = FALSE, selection = "single",
                    options = list(pageLength = 8, scrollX = TRUE, dom = "frtip"))
    })
    observe({
      g <- rv$gb; f <- if (is.null(g)) NULL else g$features
      if (is.null(f) || !nrow(f)) { updateSelectizeInput(session, "gb_feat", choices = character(), selected = character()); return() }
      ok <- which(f$Type != "source")
      lab <- sprintf("%s [%s] %s-%s (%s)", ifelse(nzchar(f$Name), f$Name, f$Product), f$Type, f$Start, f$End, f$Strand)
      updateSelectizeInput(session, "gb_feat", choices = setNames(as.character(ok), lab[ok]), selected = character())
    })
    observeEvent(input$gb_get, {
      p <- rv$pool; i <- input$pool_tbl_rows_selected
      if (is.null(p) || !length(i)) { notify("Select one sequence in the table first."); return() }
      acc <- p$Accession[i[1]]; if (!nzchar(acc)) acc <- p$ID[i[1]]
      rv$gb_err <- NULL
      g <- tryCatch(withProgress(message = "Downloading the GenBank record ...", value = 0.4, pr_gb_fetch(acc, email())),
                    error = function(e) { rv$gb_err <- conditionMessage(e); notify(conditionMessage(e), dur = 20); NULL })
      req(!is.null(g)); rv$gb <- g
      notify(sprintf("%s: %d annotated features.", g$accession, nrow(g$features)), "message", 6)
    })
    observeEvent(input$gb_use, {
      g <- rv$gb; s <- input$gb_tbl_rows_selected
      if (is.null(g)) { notify("No annotation yet: click a sequence in the table above, press 'Get GenBank annotation', then choose the gene."); return() }
      if (!length(s)) { notify("Click a gene / feature row in the annotation table first."); return() }
      fi <- gb_idx()[s[1]]
      updateRadioButtons(session, "tpl_src", selected = "gb")
      updateSelectizeInput(session, "gb_feat", selected = as.character(fi))
      updateNumericInput(session, "tgt_from", value = NA); updateNumericInput(session, "tgt_to", value = NA)
      updateTabsetPanel(session, "pt", selected = "design")
    })
    observeEvent(input$gb_splice, {
      g <- rv$gb
      if (is.null(g)) { notify("No annotation yet: click a sequence in the table above and press 'Get GenBank annotation' first."); return() }
      fi <- tryCatch(gb_fi(), error = function(e) NA_integer_)
      if (is.na(fi)) { notify("Click a gene / mRNA / CDS row in the annotation table first."); return() }
      tx <- tryCatch(pr_gb_feature_tx(g, fi), error = function(e) NULL)
      tl <- pr_gb_transcripts(g)
      id <- if (is.null(tx)) "" else tx$id
      if (!nzchar(id) || !(id %in% tl$id)) {
        notify("This feature has a single exon, so it has no exon-exon junction. Choose a gene whose mRNA / CDS has 2 or more exons (column 'Parts'), or use the PCR mode.", "warning", 12); return()
      }
      rv$tx_sel <- id
      updateRadioButtons(session, "mode", selected = "junction"); updateRadioButtons(session, "tx_src", selected = "gb")
      updateTabsetPanel(session, "pt", selected = "design")
    })
    # the structure of the annotated record is drawn as soon as it is loaded (no position has to be typed)
    gb_fi <- reactive({
      g <- rv$gb; req(g); i <- gb_idx(); s <- input$gb_tbl_rows_selected
      if (length(s) && s[1] <= length(i)) return(i[s[1]])
      f <- g$features; m <- which(f$Type %in% c("mRNA", "CDS", "tRNA", "rRNA", "ncRNA") & f$Parts >= 2L)
      if (length(m)) m[1] else NA_integer_
    })
    output$gb_viz <- renderUI({
      req(rv$gb)
      tagList(plotOutput(ns("gb_map"), height = "210px"), uiOutput(ns("gb_gene_info")), plotOutput(ns("gb_gene"), height = "250px"))
    })
    output$gb_map <- renderPlot({
      g <- rv$gb; req(g); fi <- gb_fi()
      gg <- pr_plot_record(g, fi); req(!is.null(gg)); gg
    })
    gb_tx <- reactive({ g <- rv$gb; req(g); fi <- gb_fi(); req(!is.na(fi)); pr_gb_feature_tx(g, fi) })
    output$gb_gene <- renderPlot({ tx <- gb_tx(); pr_plot_tx(tx) })
    output$gb_gene_info <- renderUI({
      g <- rv$gb; req(g)
      tx <- tryCatch(gb_tx(), error = function(e) NULL)
      if (is.null(tx)) return(helpText("Click a gene / mRNA / CDS row in the table below to see its exons."))
      ne <- length(tx$exons)
      if (ne == 1L) return(tags$p(strong(tx$name), ": one exon (no intron), ", format(tx$exons, big.mark = ","), " nt. Exon-exon junction primers need a gene with 2 or more exons."))
      gs <- pmin(tx$gstart, tx$gend); ge <- pmax(tx$gstart, tx$gend)
      tags$p(strong(tx$name), sprintf(": %d exons. ", ne), paste(sprintf("exon %d = %s-%s (%d nt)", seq_len(ne), format(gs, big.mark = ","), format(ge, big.mark = ","), tx$exons), collapse = "; "),
             if (tx$strand == "-") " (gene on the - strand: exon 1 has the highest positions)")
    })
    output$gb_dl <- downloadHandler(filename = function() paste0(if (is.null(rv$gb)) "record" else rv$gb$accession, ".gb"),
                                    content = function(file) { g <- rv$gb; req(g); writeLines(if (is.null(g$raw)) "" else g$raw, file) })
    observeEvent(input$to_design, {
      p <- rv$pool; i <- input$pool_tbl_rows_selected
      if (is.null(p) || !length(i)) { notify("Select one sequence in the table first."); return() }
      ch <- setNames(seq_len(nrow(p)), paste0(p$ID, " ", p$Species, " (", p$Length, " nt)"))
      updateSelectInput(session, "tpl_pick", choices = ch, selected = i[1])
      g <- rv$gb
      updateRadioButtons(session, "tpl_src", selected = if (!is.null(g) && identical(g$accession, p$ID[i[1]])) "gb" else "pool")
      updateTabsetPanel(session, "pt", selected = "design")
    })

    # ---------------------------------------------------------------------------------------------- 2. design
    observe({
      p <- rv$pool
      ch <- if (is.null(p)) character() else setNames(seq_len(nrow(p)), paste0(p$ID, " ", p$Species, " (", p$Length, " nt)"))
      updateSelectInput(session, "tpl_pick", choices = ch)
    })
    aln_in <- reactive({
      s <- input$aln_src
      a <- switch(s, pool = rv$aln, trim = tryCatch(aln_trim(), error = function(e) NULL), all = tryCatch(aln_all(), error = function(e) NULL),
                  paste = { x <- pr_parse_fasta(.pr_read_txt(input$aln_txt, input$aln_file))
                    if (length(x) < 3L) NULL
                    else if (length(unique(nchar(x))) == 1L) x
                    else tryCatch(withProgress(message = "Aligning ...", value = 0.5, pr_align(setNames(gsub("[-. ]", "", x), make.unique(names(x))))), error = function(e) NULL) })
      .pr_aln_vec(a)
    })
    # ---- alignment view (degenerate mode): an interactive canvas (drag to move, wheel to zoom), drawn in the browser
    output$msa_ctl <- renderUI({
      req(identical(input$mode, "degen"))
      a <- aln_in()
      if (is.null(a)) return(helpText("No alignment yet. Choose the source above: sequences of step 1 are aligned at once, or take the alignment of the MSA tab, or paste / upload a FASTA."))
      helpText(strong(sprintf("%d sequences x %d columns.", length(a), nchar(a[1]))),
               " Drag to move, mouse wheel to zoom, Shift + wheel to scroll sideways, Alt + wheel up / down; buttons + / - / Fit all / Go to primers. Solid colour = differs from the most frequent base; bar = conservation.")
    })
    sent <- new.env()
    observe({
      req(identical(input$mode, "degen"))
      a <- aln_in(); req(a)
      d <- rv$design; mf <- .pr_nz(input$dg_minfreq, 0.1); sel <- input$design_tbl_rows_selected
      msg <- list(id = ns("msa_box"), minfreq = mf)
      if (!identical(sent$a, a)) { msg$names <- as.list(unname(names(a))); msg$seqs <- as.list(unname(toupper(a))); sent$a <- a; sent$key <- NULL }
      tr <- list()
      if (!is.null(d) && identical(d$mode, "degen") && !is.null(d$pairs) && nrow(d$pairs) > 0L) {
        p <- d$pairs[min(if (length(sel)) sel[1] else 1L, nrow(d$pairs)), ]
        tr[[1]] <- list(name = "Forward", label = "Forward primer ->", start = p$F_col, end = p$F_end_col, color = "#2c7be5")
        tr[[2]] <- list(name = "Reverse", label = "<- Reverse primer", start = p$R_col, end = p$R_end_col, color = "#e8503a")
        if ("Probe" %in% names(p) && !is.na(p$Probe)) tr[[3]] <- list(name = "Probe", label = "Probe", start = p$Probe_col, end = p$Probe_col + nchar(p$Probe) - 1L, color = "#22b07d")
      }
      msg$tracks <- tr
      key <- paste(unlist(lapply(tr, function(z) c(z$start, z$end))), collapse = ",")
      if (length(tr) && !identical(sent$key, key)) { msg$focus <- list(min(vapply(tr, function(z) z$start, 0)), max(vapply(tr, function(z) z$end, 0))); sent$key <- key }
      session$sendCustomMessage("prMsaSet", msg)
    })
    aln_any <- reactive({ a <- rv$aln; if (is.null(a)) a <- .pr_aln_vec(tryCatch(aln_trim(), error = function(e) NULL)); if (is.null(a)) a <- .pr_aln_vec(tryCatch(aln_all(), error = function(e) NULL)); a })
    output$aln_info <- renderUI({
      a <- aln_in()
      if (is.null(a)) return(helpText("No alignment: align the sequences in step 1, or use the MSA tab."))
      ok <- length(unique(nchar(a))) == 1L
      tags$p(class = if (ok) "text-success" else "text-danger", sprintf("%d sequences x %s columns%s", length(a), paste(unique(nchar(a)), collapse = "/"), if (ok) "" else " - NOT aligned"))
    })
    tpl <- reactive({
      src <- input$tpl_src; tf <- .pr_nz(input$tgt_from, NA); tt <- .pr_nz(input$tgt_to, NA)
      if (is.na(tf) != is.na(tt)) stop("Fill both 'from' and 'to' of the region, or leave both empty.")
      mk <- function(name, sq, exons = NULL, region = c(tf, tt)) {
        sq <- .pr_clean(sq)
        list(name = name, seq = sq, exons = exons, region = region, win = list(seq = sq, from = 1L, to = nchar(sq), strand = "+"), gb = FALSE)
      }
      if (src == "paste") {
        x <- pr_parse_fasta(.pr_read_txt(input$tpl_txt, input$tpl_file)); req(length(x) > 0L)
        x <- gsub("[-. ]", "", x)
        if (identical(input$jn_type, "fasta") && identical(input$mode, "junction") && length(x) >= 2L)
          return(mk("exons", paste(.pr_clean(x), collapse = ""), exons = nchar(.pr_clean(x)), region = c(NA, NA)))
        return(mk(names(x)[1], x[1]))
      }
      if (src == "pool") {
        p <- rv$pool; i <- as.integer(input$tpl_pick); req(!is.null(p), length(i) == 1L, !is.na(i))
        return(mk(p$ID[i], p$Sequence[i]))
      }
      if (src == "gb") {
        g <- rv$gb; req(!is.null(g)); n <- nchar(g$seq)
        fi <- suppressWarnings(as.integer(.pr_nz(input$gb_feat, NA)))
        fstr <- if (!is.na(fi)) g$features$Strand[fi] else "+"
        if (is.na(tf)) {
          if (is.na(fi)) stop("Choose the target gene / feature of the GenBank record (or type a region).")
          tf <- g$features$Start[fi]; tt <- g$features$End[fi]
        }
        tf <- max(1L, as.integer(tf)); tt <- min(n, as.integer(tt))
        strand <- if (isTRUE(input$gb_orient) && identical(fstr, "-")) "-" else "+"
        fl <- as.integer(.pr_nz(input$gb_flank, 300))
        w <- pr_gb_window(g, tf - fl, tt + fl, strand)
        r1 <- pr_win_pos(w, tf); r2 <- pr_win_pos(w, tt)
        return(list(name = paste0(g$accession, if (!is.na(fi)) paste0(" ", g$features$Name[fi]) else ""), seq = .pr_clean(w$seq), exons = NULL,
                    region = c(min(r1, r2), max(r1, r2)), win = w, gb = TRUE))
      }
      if (src == "spliced") { sp <- rv$spliced; req(!is.null(sp)); return(mk(sp$name, sp$seq, exons = sp$exons, region = c(NA, NA))) }
      if (src == "cons") { a <- aln_any(); req(!is.null(a)); return(mk("consensus of the alignment", .pr_consensus(a))) }
      cs <- tryCatch(consensus(), error = function(e) NULL); cs <- if (is.null(cs)) NULL else unlist(cs); req(length(cs) > 0L)
      mk("assembly consensus", as.character(cs[1]))
    })
    output$tpl_info <- renderUI({
      t <- tryCatch(tpl(), error = function(e) NULL)
      if (is.null(t)) return(NULL)
      extra <- if (isTRUE(t$gb)) sprintf("window %s-%s of the record%s", format(t$win$from, big.mark = ","), format(t$win$to, big.mark = ","),
                                         if (t$win$strand == "-") ", read as the reverse complement (gene direction)" else "") else ""
      tags$p(class = "text-muted", sprintf("%s: %s nt%s%s", t$name, format(nchar(t$seq), big.mark = ","),
                                           if (length(t$exons)) paste0("; exons: ", paste(t$exons, collapse = ", ")) else "", if (nzchar(extra)) paste0("; ", extra) else ""))
    })

    # ---- mRNA mode: the exons come from the annotation, nothing has to be typed
    output$tx_pick_ui <- renderUI({
      g <- rv$gb
      if (is.null(g)) return(helpText("No GenBank record yet. In step 1 find a sequence, select it and press 'Get GenBank annotation' (or load a GenBank file): ",
                                      "the exons of its mRNA / CDS are then read automatically and drawn here."))
      tl <- pr_gb_transcripts(g)
      if (!nrow(tl)) return(tags$div(class = "alert alert-warning", "This record has no annotated transcript with 2 or more exons (no spliced mRNA / CDS, no exon features), ",
                                     "so there is no exon-exon junction to design on. Use the PCR mode, or paste an mRNA."))
      sel <- rv$tx_sel; if (is.null(sel) || !(sel %in% tl$id)) sel <- tl$id[1]
      selectInput(ns("tx_pick"), "Transcript (spliced sequence)", choices = setNames(tl$id, tl$label), selected = sel, width = "100%")
    })
    tx_obj <- reactive({
      if (identical(input$tx_src, "paste")) {
        x <- pr_parse_fasta(.pr_read_txt(input$tx_txt, input$tx_file)); req(length(x) > 0L)
        nm <- .pr_nz(names(x)[1], "mRNA"); x <- .pr_clean(gsub("[-. ]", "", x))
        if (identical(input$jn_type, "fasta") && length(x) >= 2L) return(.pr_tx_from_exons(nchar(x), paste(x, collapse = ""), "exons"))
        sq <- x[1]
        v <- suppressWarnings(as.numeric(strsplit(gsub("[^0-9.,; ]", "", .pr_nz(input$jn_txt, "")), "[ ,;]+")[[1]])); v <- v[!is.na(v)]
        if (identical(input$jn_type, "len")) { if (length(v) < 2L) stop("Give at least two exon lengths."); v <- cumsum(v)[-length(v)] }
        if (!length(v)) stop("Give the exon boundaries (the last base of each exon, e.g. 120, 340) or choose 'Each FASTA record is one exon'.")
        return(.pr_tx_from_junc(v, nchar(sq), nm, sq))
      }
      g <- rv$gb; req(!is.null(g)); id <- .pr_nz(input$tx_pick, ""); req(nzchar(id))
      pr_gb_transcript(g, id)
    })
    output$tx_map <- renderPlot({
      tx <- tryCatch(tx_obj(), error = function(e) NULL); req(!is.null(tx))
      pr_plot_tx(tx)
    })
    output$tx_tbl <- DT::renderDT({
      tx <- tryCatch(tx_obj(), error = function(e) NULL); req(!is.null(tx))
      o <- data.frame(Exon = seq_along(tx$exons), `Transcript start` = tx$tstart, `Transcript end` = tx$tend, Length = tx$exons, check.names = FALSE)
      if (length(tx$gstart)) { o$`Record start` <- pmin(tx$gstart, tx$gend); o$`Record end` <- pmax(tx$gstart, tx$gend) }
      if (length(tx$exons) > 1L) o$`Junction with next exon (last base)` <- c(tx$tend[-length(tx$exons)], NA)
      DT::datatable(o, rownames = FALSE, selection = "none", options = list(dom = "t", pageLength = 100, scrollX = TRUE, scrollY = "170px", scrollCollapse = TRUE))
    })
    observeEvent(input$design_go, {
      md <- input$mode; p <- prm()
      res <- tryCatch(withProgress(message = "Designing primers ...", value = 0.3, {
        if (md == "degen") {
          a <- aln_in(); if (is.null(a)) stop("No alignment. Align the sequences in step 1 first, or choose the MSA tab alignment.")
          r <- pr_design_degen(a, p, max_deg = .pr_nz(input$dg_max, 64), minfreq = .pr_nz(input$dg_minfreq, 0.1), max_gap = .pr_nz(input$dg_gap, 0.2),
                               mm = as.integer(.pr_nz(input$dg_mm, 0)), min_cov = .pr_nz(input$dg_cov, 0.8), from = as.integer(.pr_nz(input$dg_from, 1)),
                               to = if (is.na(.pr_nz(input$dg_to, NA))) NA else as.integer(input$dg_to), probe = isTRUE(input$probe),
                               probe_delta = .pr_nz(input$probe_delta, 8))
          list(mode = "degen", pairs = r$pairs, ctx = r$ctx, mk = r$mk, nF = r$nF, nR = r$nR, n = r$ctx$w)
        } else if (md == "junction") {
          tx <- tx_obj()
          if (length(tx$exons) < 2L) stop("This transcript has a single exon, so it has no exon-exon junction. Choose a transcript with introns, or use the PCR mode.")
          j <- cumsum(tx$exons)[-length(tx$exons)]
          r <- pr_design_junction(tx$seq, j, p, span = input$jn_span, ov5 = as.integer(.pr_nz(input$jn_ov5, 7)), ov3 = as.integer(.pr_nz(input$jn_ov3, 4)),
                                  probe = isTRUE(input$probe), probe_delta = .pr_nz(input$probe_delta, 8))
          list(mode = "junction", pairs = r$pairs, n = r$n, junc = r$junc, tpl = tx$seq, win = list(seq = tx$seq, from = 1L, to = nchar(tx$seq), strand = "+"), feat = NULL, tx = tx)
        } else {
          t <- tpl(); tf <- t$region[1]; tt <- t$region[2]
          ft <- if (isTRUE(t$gb)) pr_gb_win_features(rv$gb, t$win) else NULL
          if (md == "pcr") {
            r <- pr_design_pcr(t$seq, tf, tt, p, probe = isTRUE(input$probe), probe_delta = .pr_nz(input$probe_delta, 8),
                               inside = identical(input$tgt_role, "inside"))
            list(mode = "pcr", pairs = r$pairs, nF = r$nF, nR = r$nR, n = r$n, region = c(tf, tt), tpl = t$seq, win = t$win, feat = ft)
          } else if (md == "seq") {
            r <- pr_design_seq(t$seq, tf, tt, p, offset = as.integer(.pr_nz(input$sq_offset, 50)), readlen = as.integer(.pr_nz(input$sq_read, 700)),
                               overlap = as.integer(.pr_nz(input$sq_overlap, 100)), both = isTRUE(input$sq_both))
            list(mode = "seq", prim = r, n = nchar(t$seq), region = c(tf, tt), tpl = t$seq, win = t$win, feat = ft)
          }
        }
      }), error = function(e) { notify(conditionMessage(e), dur = 15); NULL })
      req(!is.null(res)); rv$design <- res
    })

    rec_rng <- function(win, a, b) { x <- pr_rec_pos(win, a); y <- pr_rec_pos(win, b); cbind(pmin(x, y), pmax(x, y)) }
    win_of <- function(d) if (is.null(d$win)) list(from = 1L, to = d$n, strand = "+") else d$win
    design_disp <- reactive({
      d <- rv$design; if (is.null(d)) return(NULL)
      w <- win_of(d); fs <- if (w$strand == "-") "-" else "+"; rs <- if (w$strand == "-") "+" else "-"
      if (d$mode == "seq") {
        x <- d$prim; if (is.null(x)) return(NULL)
        r <- rec_rng(w, x$start, x$end); cv <- rec_rng(w, x$cover_from, x$cover_to)
        return(data.frame(Name = x$name, Strand = ifelse(x$strand == "F", fs, rs), `Sequence 5'->3'` = x$seq, Start = r[, 1], End = r[, 2], Length = x$len,
                          `Tm (C)` = x$tm, `GC %` = x$gc, `Read covers` = paste0(cv[, 1], "-", cv[, 2]), check.names = FALSE, stringsAsFactors = FALSE))
      }
      x <- d$pairs; if (is.null(x)) return(NULL)
      if (d$mode == "degen") {
        o <- data.frame(`#` = seq_len(nrow(x)), `Forward 5'->3'` = x$F_seq, `F deg` = x$F_deg, `F Tm` = round(x$F_tm, 1), `F Tm range` = x$F_tm_range, `F GC %` = round(x$F_gc),
                        `Reverse 5'->3'` = x$R_seq, `R deg` = x$R_deg, `R Tm` = round(x$R_tm, 1), `R Tm range` = x$R_tm_range, `R GC %` = round(x$R_gc),
                        `Product (median)` = round(x$product), `Product range` = x$prod_range, `Sequences covered` = paste0(round(100 * x$both_cov), "%"),
                        `dTm` = round(x$dTm, 1), `Hetero-dimer dG` = round(x$heterodimer, 1), `F col` = x$F_col, `R col (end)` = x$R_end_col,
                        check.names = FALSE, stringsAsFactors = FALSE)
        if ("Probe" %in% names(x)) {
          o$Probe <- ifelse(is.na(x$Probe), "none found", x$Probe); o$`Probe Tm` <- round(x$Probe_Tm, 1); o$`Probe col` <- x$Probe_col
          o$`Probe coverage` <- ifelse(is.na(x$Probe_cov), "", paste0(round(100 * x$Probe_cov), "%")); o$`Probe criteria` <- ifelse(is.na(x$Probe_note), "", x$Probe_note)
        }
        return(o)
      }
      fr <- rec_rng(w, x$F_start, x$F_end); rr <- rec_rng(w, x$R_start, x$R_end)
      o <- data.frame(`#` = seq_len(nrow(x)), `Forward 5'->3'` = x$F_seq, `F strand` = fs, `F start` = fr[, 1], `F end` = fr[, 2], `F Tm` = round(x$F_tm, 1), `F GC %` = round(x$F_gc),
                      `Reverse 5'->3'` = x$R_seq, `R strand` = rs, `R start` = rr[, 1], `R end` = rr[, 2], `R Tm` = round(x$R_tm, 1), `R GC %` = round(x$R_gc),
                      Product = x$product, dTm = round(x$dTm, 1), `Hetero-dimer dG` = round(x$heterodimer, 1), Penalty = round(x$penalty, 2),
                      check.names = FALSE, stringsAsFactors = FALSE)
      if (d$mode == "junction") {
        o$`Junction (exons)` <- x$junc_label; o$`F exon(s)` <- x$exons_F; o$`R exon(s)` <- x$exons_R; o$`Product exons` <- x$exons_prod; o$Type <- x$type
        first <- c("#", "Junction (exons)", "F exon(s)", "R exon(s)", "Forward 5'->3'", "Reverse 5'->3'")
        o <- o[, c(first, setdiff(names(o), first)), drop = FALSE]
      }
      if ("Probe" %in% names(x)) {
        pr <- rec_rng(w, x$Probe_start, x$Probe_start + nchar(ifelse(is.na(x$Probe), "", x$Probe)) - 1L)
        o$Probe <- x$Probe; o$`Probe Tm` <- round(x$Probe_Tm, 1); o$`Probe start` <- pr[, 1]; o$`Probe end` <- pr[, 2]
      }
      o
    })
    output$design_msg <- renderUI({
      d <- rv$design
      if (is.null(d)) return(helpText("Choose what to design, give the template, then press Design primers."))
      ok <- if (d$mode == "seq") !is.null(d$prim) else !is.null(d$pairs)
      if (ok && d$mode == "degen" && "Probe" %in% names(d$pairs) && all(is.na(d$pairs$Probe)))
        return(tags$div(class = "alert alert-warning", sprintf("%d result(s), but no probe was found between the primers, even with relaxed criteria (coverage down to 50 %%, 1 mismatch, up to 16 variants, Tm within 8 C of the primers + %s C): ",
                          nrow(design_disp()), .pr_nz(input$probe_delta, 8)), "there is no 18-30 nt window with 5' base not G and GC 30-80 % in the region between these primers. Try another primer pair (the probe is searched for the first 10 pairs), a wider product size, or a different 'Probe Tm above the primers'."))
      if (ok && d$mode == "degen" && "Probe_note" %in% names(d$pairs) && any(!is.na(d$pairs$Probe_note) & d$pairs$Probe_note != "strict"))
        return(tags$div(class = "alert alert-info", sprintf("%d result(s). For some pairs the probe needed relaxed criteria (see 'Probe criteria' and 'Probe coverage' in the table).", nrow(design_disp()))))
      if (ok) return(tags$p(class = "text-success", sprintf("%d result(s).", nrow(design_disp()))))
      tags$div(class = "alert alert-warning",
               if (d$mode == "seq") "No sequencing primer fits. Relax the Tm / GC limits or the read settings."
               else if (d$mode == "junction") "No primer pair was found around the exon-exon junctions with these settings. Widen the Tm / GC range or the product size, allow either primer to span the junction, or choose 'Primers in two different exons'."
               else sprintf("No primer pair found (%s forward and %s reverse candidate primers passed the filters). Widen the Tm / GC range, the product size or the allowed Tm difference%s.",
                            .pr_nz(d$nF, 0), .pr_nz(d$nR, 0),
                            if (d$mode == "degen") ", lower the coverage, or raise the degeneracy" else ""))
    })
    output$design_tbl <- DT::renderDT({
      x <- design_disp(); req(x)
      x <- as.data.frame(lapply(x, function(v) if (is.list(v)) vapply(v, function(z) paste(z, collapse = ", "), "") else if (is.factor(v)) as.character(v) else if (is.numeric(v)) { v[!is.finite(v)] <- NA; v } else v),
                         check.names = FALSE, stringsAsFactors = FALSE)
      names(x) <- make.unique(trimws(gsub("[^A-Za-z0-9 %()/'.+-]", "", names(x))))
      names(x)[names(x) == ""] <- "No"
      DT::datatable(x, rownames = FALSE, selection = "single", options = list(pageLength = 10, lengthMenu = c(5, 10, 25, 50), scrollX = TRUE, dom = "ltip"))
    }, server = FALSE)
    sel_row <- reactive({ i <- input$design_tbl_rows_selected; if (length(i)) i[1] else 1L })
    sel_primers <- function() {
      d <- rv$design; x <- design_disp(); req(d, x); i <- sel_row()
      if (d$mode == "seq") return(data.frame(Name = x$Name[i], Sequence = x$`Sequence 5'->3'`[i], stringsAsFactors = FALSE))
      p <- d$pairs[i, ]
      out <- data.frame(Name = c(sprintf("Pair%d_F", i), sprintf("Pair%d_R", i)), Sequence = c(p$F_seq, p$R_seq), stringsAsFactors = FALSE)
      if ("Probe" %in% names(p) && !is.na(p$Probe)) out <- rbind(out, data.frame(Name = sprintf("Pair%d_probe", i), Sequence = p$Probe, stringsAsFactors = FALSE))
      out
    }
    observeEvent(input$send_ana, {
      s <- tryCatch(sel_primers(), error = function(e) NULL); req(!is.null(s))
      add_primers(paste0(s$Name, "\t", s$Sequence))
      updateTabsetPanel(session, "pt", selected = "ana")
    })
    output$design_dl <- downloadHandler(filename = function() "primers_design.csv", content = function(file) utils::write.csv(design_disp(), file, row.names = FALSE))
    detail_rows <- reactive({
      d <- rv$design; x <- design_disp(); req(d, x); i <- sel_row()
      w <- win_of(d); fs <- if (w$strand == "-") "-" else "+"; rs <- if (w$strand == "-") "+" else "-"
      tmgc <- function(tm, gc) list(Tm = round(tm, 1), `GC %` = round(gc))
      if (d$mode == "seq") {
        q <- d$prim[i, ]; r <- rec_rng(w, q$start, q$end)
        return(list(rows = data.frame(Primer = q$name, `Sequence (5'->3')` = q$seq, `Template strand` = if (q$strand == "F") fs else rs, Start = r[1, 1], Stop = r[1, 2],
                                      Length = q$len, Tm = round(q$tm, 1), `GC %` = round(q$gc), check.names = FALSE, stringsAsFactors = FALSE), note = NULL))
      }
      p <- d$pairs[i, ]
      if (d$mode == "degen") {
        ctx <- d$ctx; refn <- if (length(ctx$names) && !is.na(ctx$names[1])) ctx$names[1] else "sequence 1"
        rp <- function(col) as.integer(ctx$cumng[1, col])                       # position in the first sequence of the alignment
        mkrow <- function(nm, sq, c1, c2, tm, gc, cov) data.frame(Primer = nm, `Sequence (5'->3')` = sq, `Alignment columns` = paste0(c1, "-", c2),
                          `Position in ref` = paste0(rp(c1), "-", rp(c2)), Length = nchar(sq), Tm = round(tm, 1), `GC %` = round(gc),
                          `Variants` = unname(vapply(sq, .pr_deg, 0)), `Sequences matched` = paste0(round(100 * cov), "%"), check.names = FALSE, stringsAsFactors = FALSE, row.names = NULL)
        rows <- rbind(mkrow("Forward", p$F_seq, p$F_col, p$F_end_col, p$F_tm, p$F_gc, p$F_cov), mkrow("Reverse", p$R_seq, p$R_col, p$R_end_col, p$R_tm, p$R_gc, p$R_cov))
        note <- sprintf("Position in ref = position in the first sequence of the alignment (%s). Product %s bp (median; range %s).", refn, round(p$product), p$prod_range)
        if ("Probe" %in% names(p)) {
          if (!is.na(p$Probe)) {
            pe <- p$Probe_col + nchar(p$Probe) - 1L
            rows <- rbind(rows, mkrow("Probe", p$Probe, p$Probe_col, pe, p$Probe_Tm, .pr_gc(p$Probe), p$Probe_cov))
            if (!identical(p$Probe_note, "strict")) note <- paste(note, sprintf("The probe needed relaxed criteria (%s): check its coverage.", p$Probe_note))
          } else note <- paste(note, "No probe was found for this pair (try another pair).")
        }
        return(list(rows = rows, note = note))
      }
      fr <- rec_rng(w, p$F_start, p$F_end); rr <- rec_rng(w, p$R_start, p$R_end)
      rows <- data.frame(Primer = c("Forward", "Reverse"), `Sequence (5'->3')` = c(p$F_seq, p$R_seq), `Template strand` = c(fs, rs), Start = c(fr[1, 1], rr[1, 1]),
                         Stop = c(fr[1, 2], rr[1, 2]), Length = nchar(c(p$F_seq, p$R_seq)), Tm = round(c(p$F_tm, p$R_tm), 1), `GC %` = round(c(p$F_gc, p$R_gc)),
                         check.names = FALSE, stringsAsFactors = FALSE)
      if ("Probe" %in% names(p) && !is.na(p$Probe)) {
        pb <- rec_rng(w, p$Probe_start, p$Probe_start + nchar(p$Probe) - 1L)
        rows <- rbind(rows, data.frame(Primer = "Probe", `Sequence (5'->3')` = p$Probe, `Template strand` = if (identical(p$Probe_strand, "+")) fs else rs, Start = pb[1, 1], Stop = pb[1, 2],
                                       Length = nchar(p$Probe), Tm = round(p$Probe_Tm, 1), `GC %` = round(.pr_gc(p$Probe)), check.names = FALSE, stringsAsFactors = FALSE))
      }
      if (d$mode == "junction" && !is.null(d$tx)) {
        ex <- c(p$exons_F, p$exons_R)
        if (nrow(rows) == 3L) ex <- c(ex, .pr_exon_span(d$tx, p$Probe_start, p$Probe_start + nchar(p$Probe) - 1L))
        rows$`Exon(s)` <- ex; rows$`Template strand` <- NULL
        return(list(rows = rows, note = sprintf("Positions are on the spliced transcript. Product %d bp over exons %s; %s (junction between exons %s).",
                                                p$product, p$exons_prod, p$type, p$junc_label)))
      }
      list(rows = rows, note = sprintf("Product length: %d bp, from %s to %s.", p$product, format(min(fr, rr), big.mark = ","), format(max(fr, rr), big.mark = ",")))
    })
    output$design_detail <- renderUI({
      dr <- tryCatch(detail_rows(), error = function(e) NULL); if (is.null(dr)) return(NULL)
      df <- dr$rows
      tagList(
        tags$table(class = "table table-sm",
          tags$thead(tags$tr(lapply(names(df), tags$th))),
          tags$tbody(lapply(seq_len(nrow(df)), function(k) tags$tr(lapply(seq_along(df), function(j) {
            v <- df[[j]][k]
            if (names(df)[j] == "Sequence (5'->3')") tags$td(class = "pr-seq", v) else tags$td(format(v, big.mark = ","))
          }))))),
        if (!is.null(dr$note)) tags$p(strong(dr$note)))
    })
    output$seq_view <- renderUI({
      d <- rv$design; dr <- tryCatch(detail_rows(), error = function(e) NULL)
      req(d, dr)
      if (identical(d$mode, "degen")) {
        ctx <- d$ctx; mk <- d$mk; req(!is.null(mk)); p <- d$pairs[sel_row(), ]
        cons <- paste(ifelse(mk > 0L, .pr_code_of_mask[pmax(mk, 1L)], "-"), collapse = "")
        marks <- data.frame(a = c(p$F_col, p$R_col), b = c(p$F_end_col, p$R_end_col), cls = c("pr-f", "pr-r"), stringsAsFactors = FALSE)
        hasp <- "Probe_col" %in% names(p) && !is.na(p$Probe_col)
        if (hasp) marks <- rbind(marks, data.frame(a = p$Probe_col, b = p$Probe_col + nchar(p$Probe) - 1L, cls = "pr-p", stringsAsFactors = FALSE))
        a <- max(1L, p$F_col - 30L); b <- min(ctx$w, p$R_end_col + 30L)
        if (b - a > 3000L) return(helpText("The product is longer than 3000 columns: the sequence view is skipped."))
        return(tagList(tags$div(class = "mb-2", tags$span(class = "pr-f px-2 me-2", "Forward primer"), tags$span(class = "pr-r px-2 me-2", "Reverse primer"),
                                if (hasp) tags$span(class = "pr-p px-2", "Probe")),
                       .pr_seq_html(cons, list(from = 1L, to = ctx$w, strand = "+"), a, b, marks),
                       helpText("IUPAC consensus of the alignment (every base present in at least 'Min base frequency' of the sequences; - = gap). Numbers = alignment columns. ",
                                "The forward primer reads as shown; the reverse primer is the reverse complement of the highlighted columns.")))
      }
      req(d$mode %in% c("pcr", "junction"), !is.null(d$tpl))
      p <- d$pairs[sel_row(), ]
      marks <- data.frame(a = c(p$F_start, p$R_start), b = c(p$F_end, p$R_end), cls = c("pr-f", "pr-r"), stringsAsFactors = FALSE)
      if ("Probe" %in% names(p) && !is.na(p$Probe))
        marks <- rbind(marks, data.frame(a = p$Probe_start, b = p$Probe_start + nchar(p$Probe) - 1L, cls = "pr-p", stringsAsFactors = FALSE))
      a <- max(1L, p$F_start - 30L); b <- min(nchar(d$tpl), p$R_end + 30L)
      if (b - a > 3000L) return(helpText("The product is longer than 3000 nt: the sequence view is skipped."))
      tagList(tags$div(class = "mb-2", tags$span(class = "pr-f px-2 me-2", "Left primer"), tags$span(class = "pr-r px-2 me-2", "Right primer"),
                       if (any(marks$cls == "pr-p")) tags$span(class = "pr-p px-2", "Probe")),
              .pr_seq_html(d$tpl, win_of(d), a, b, marks, junc = if (identical(d$mode, "junction")) d$junc else NULL),
              helpText("Lavender = forward (left) primer, yellow = reverse (right) primer (their footprint on the template), green = probe. Numbers = position of the first base of each line; 50 bases per line in blocks of 10.",
                       if (identical(d$mode, "junction")) " A purple bar after a base = the end of an exon (exon-exon junction); a primer crossing a bar spans two exons."))
    })
    output$map <- renderPlot({
      d <- rv$design; req(d, design_disp()); i <- sel_row(); w <- win_of(d)
      th <- ggplot2::theme_minimal(base_size = 13) + ggplot2::theme(axis.text.y = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(), legend.position = "top")
      cols <- c(Forward = "#2c7be5", Reverse = "#e8503a", Probe = "#22b07d", gene = "#9aa5b1", CDS = "#2e9d6a", exon = "#16a085", mRNA = "#7f8c8d", rRNA = "#b9770e", tRNA = "#8e44ad")
      xs <- ggplot2::scale_x_continuous(labels = function(b) format(round(pr_rec_pos(w, b)), big.mark = ",", trim = TRUE))
      ftr <- function() {
        if (is.null(d$feat)) return(list())
        f <- d$feat; f$ymin <- -1.5 - 0.4 * (f$lane - 1); f$ymax <- f$ymin + 0.28
        lab <- f[!duplicated(paste(f$label, f$lane)), ]
        list(ggplot2::geom_rect(data = f, ggplot2::aes(xmin = start, xmax = end, ymin = ymin, ymax = ymax, fill = type)),
             ggplot2::geom_text(data = lab, ggplot2::aes(x = (start + end) / 2, y = ymax + 0.1, label = label), size = 3, check_overlap = TRUE))
      }
      ylim <- if (is.null(d$feat)) c(-1.2, 1.1) else c(-2.4, 1.1)
      if (d$mode == "seq") {
        x <- d$prim; x$row <- seq_len(nrow(x)); x$y <- ifelse(x$strand == "F", 0.3, -0.3)
        x$what <- ifelse(x$strand == "F", "Forward", "Reverse"); x$sel <- x$row == i
        return(ggplot2::ggplot() +
          ggplot2::annotate("segment", x = 1, xend = d$n, y = 0, yend = 0, linewidth = 1.2, colour = "grey40") +
          ggplot2::annotate("rect", xmin = max(1, .pr_nz(d$region[1], 1)), xmax = .pr_nz(d$region[2], d$n), ymin = -1, ymax = 1, alpha = 0.08) +
          ggplot2::geom_segment(data = x, ggplot2::aes(x = cover_from, xend = cover_to, y = y * 2, yend = y * 2, colour = what), linewidth = ifelse(x$sel, 3, 1.2), alpha = 0.6,
                                arrow = grid::arrow(length = grid::unit(0.12, "inches"))) +
          ggplot2::geom_rect(data = x, ggplot2::aes(xmin = start, xmax = end, ymin = y - 0.12, ymax = y + 0.12, fill = what)) +
          ftr() + xs + ggplot2::scale_y_continuous(limits = c(min(ylim[1], -0.8), 0.8)) +
          ggplot2::scale_fill_manual(values = cols) + ggplot2::scale_colour_manual(values = cols) +
          ggplot2::labs(x = "position on the record (nt)", y = NULL, fill = NULL, colour = NULL, subtitle = "primer = box; arrow = part of the sequence read (after the lost bases)") + th)
      }
      if (d$mode == "junction" && !is.null(d$tx)) {
        p <- d$pairs[i, ]
        pb <- if ("Probe" %in% names(p) && !is.na(p$Probe)) c(p$Probe_start, p$Probe_start + nchar(p$Probe) - 1L) else NULL
        return(pr_plot_tx(d$tx, pairs = d$pairs, sel = i, zoom = isTRUE(input$map_zoom)))
      }
      if (d$mode == "degen") {
        ctx <- d$ctx; p <- d$pairs[i, ]
        a <- max(1L, p$F_col - 40L); b <- min(ctx$w, p$R_end_col + 40L)        # zoom on the amplicon (the alignment can have thousands of columns)
        df <- data.frame(col = a:b, ident = ctx$ident[a:b])
        rects <- data.frame(xmin = c(p$F_col, p$R_col), xmax = c(p$F_end_col, p$R_end_col), what = c("Forward", "Reverse"), ymin = 1.05, ymax = 1.25, stringsAsFactors = FALSE)
        hasp <- "Probe_col" %in% names(p) && !is.na(p$Probe_col)
        if (hasp) rects <- rbind(rects, data.frame(xmin = p$Probe_col, xmax = p$Probe_col + nchar(p$Probe) - 1, what = "Probe", ymin = 1.32, ymax = 1.52, stringsAsFactors = FALSE))
        rects$xmin <- rects$xmin - 0.5; rects$xmax <- rects$xmax + 0.5
        return(ggplot2::ggplot() +
          ggplot2::geom_col(data = df, ggplot2::aes(col, ident), fill = "grey75", width = 1) +
          ggplot2::geom_rect(data = rects, ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = what)) +
          ggplot2::geom_text(data = rects, ggplot2::aes(x = (xmin + xmax) / 2, y = ymax + 0.07, label = what, colour = what), size = 3.6, fontface = "bold", show.legend = FALSE) +
          ggplot2::annotate("segment", x = p$F_col, xend = p$R_end_col, y = -0.08, yend = -0.08, colour = "grey30") +
          ggplot2::annotate("text", x = (p$F_col + p$R_end_col) / 2, y = -0.2, size = 3.6, label = sprintf("product about %s bp (median over the sequences)", round(p$product))) +
          ggplot2::scale_fill_manual(values = cols) + ggplot2::scale_colour_manual(values = cols) +
          ggplot2::scale_y_continuous(limits = c(-0.3, if (hasp) 1.7 else 1.45), breaks = c(0, 0.5, 1)) +
          ggplot2::labs(x = "alignment column", y = "conservation", fill = NULL,
                        subtitle = sprintf("columns %s-%s of %s; grey bars = fraction of the sequences with the most common base", a, b, format(ctx$w, big.mark = ","))) +
          th + ggplot2::theme(axis.text.y = ggplot2::element_text()))
      }
      p <- d$pairs[i, ]
      rects <- data.frame(xmin = c(p$F_start, p$R_start), xmax = c(p$F_end, p$R_end), ymin = c(0.15, -0.55), ymax = c(0.55, -0.15), what = c("Forward", "Reverse"))
      if ("Probe" %in% names(p) && !is.na(p$Probe)) rects <- rbind(rects, data.frame(xmin = p$Probe_start, xmax = p$Probe_start + nchar(p$Probe) - 1, ymin = 0.65, ymax = 0.95, what = "Probe"))
      g <- ggplot2::ggplot() +
        ggplot2::annotate("segment", x = 1, xend = d$n, y = 0, yend = 0, linewidth = 1.2, colour = "grey40")
      if (d$mode == "pcr" && !any(is.na(d$region))) g <- g + ggplot2::annotate("rect", xmin = d$region[1], xmax = d$region[2], ymin = -1, ymax = 1, alpha = 0.1)
      if (d$mode == "junction") g <- g + ggplot2::geom_vline(xintercept = d$junc + 0.5, linetype = "dashed", colour = "#8e44ad")
      g + ggplot2::geom_rect(data = rects, ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = what)) +
        ggplot2::annotate("segment", x = p$F_start, xend = p$R_end, y = -0.85, yend = -0.85, colour = "grey30") +
        ggplot2::annotate("text", x = (p$F_start + p$R_end) / 2, y = -1.05, label = paste0(p$product, " bp"), size = 4) +
        ftr() + xs + ggplot2::scale_fill_manual(values = cols) + ggplot2::scale_y_continuous(limits = ylim) +
        ggplot2::labs(x = "position on the record (nt)", y = NULL, fill = NULL,
                      subtitle = if (d$mode == "junction") "dashed lines = exon-exon junctions" else if (!is.null(d$feat)) "gray / green boxes = genes and CDS of the GenBank record" else NULL) + th
    })

    # ---------------------------------------------------------------------------------------------- 3. quality
    observeEvent(input$ana_go, {
      prim <- pr_parse_primers(input$prim_txt)
      if (is.null(prim) || !nrow(prim)) { notify("Type or paste at least one primer (8 nt or longer)."); return() }
      p <- prm()
      res <- tryCatch(withProgress(message = "Analysing ...", value = 0.4, {
        list(ana = pr_analyse(prim, p), pairs = pr_pair_table(prim, p))
      }), error = function(e) { notify(conditionMessage(e)); NULL })
      req(!is.null(res))
      rv$prim <- prim; rv$ana <- res$ana; rv$pairs <- res$pairs; rv$bind <- NULL
      if (isTRUE(input$bind_on)) {
        seqs <- switch(input$bind_src,
                       pool = if (is.null(rv$pool)) NULL else setNames(rv$pool$Sequence, make.unique(rv$pool$ID)),
                       aln = { a <- aln_any(); if (is.null(a)) NULL else setNames(gsub("[-.]", "", a), names(a)) },
                       paste = { x <- pr_parse_fasta(input$bind_txt); if (length(x)) x else NULL })
        if (is.null(seqs)) { notify("No sequences to check the binding in (step 1, MSA tab or pasted FASTA).", "warning"); return() }
        mx <- as.integer(.pr_nz(input$bind_max, 300))
        b <- withProgress(message = "Scanning the sequences ...", value = 0, {
          lapply(seq_len(nrow(prim)), function(k) {
            incProgress(1 / nrow(prim), detail = prim$Name[k])
            pr_bind_scan(prim$Sequence[k], seqs, mx)
          })
        })
        names(b) <- prim$Name; rv$bind <- b
      }
    })
    observe({
      pr <- rv$prim; nm <- if (is.null(pr)) character() else pr$Name
      updateSelectizeInput(session, "bl_fpick", choices = c("", nm), selected = if (length(nm)) nm[1] else "")
      updateSelectizeInput(session, "bl_rpick", choices = c("", nm), selected = if (length(nm) > 1L) nm[2] else "")
    })
    output$ana_msg <- renderUI({
      if (is.null(rv$ana)) return(helpText("Paste primers on the left (or send them from Design / the primer list), then press Analyse."))
      v <- table(factor(rv$ana$Verdict, c("OK", "CHECK", "FAIL")))
      np <- sum(rv$ana$Role == "Probe")
      tags$p(sprintf("%d primer(s)%s: %d OK, %d to check, %d fail.", nrow(rv$ana) - np, if (np) sprintf(" and %d probe(s)", np) else "", v[["OK"]], v[["CHECK"]], v[["FAIL"]]))
    })
    output$ana_tbl <- DT::renderDT({
      x <- rv$ana; req(x)
      .pr_style_verdict(DT::datatable(x, rownames = FALSE, selection = "single", options = list(pageLength = 15, scrollX = TRUE, dom = "t")))
    })
    output$pair_tbl <- DT::renderDT({
      x <- rv$pairs; if (is.null(x)) return(DT::datatable(data.frame(Note = "Give at least two primers."), rownames = FALSE, options = list(dom = "t")))
      .pr_style_verdict(DT::datatable(x, rownames = FALSE, selection = list(mode = "single", selected = 1), options = list(pageLength = 15, scrollX = TRUE, dom = "t")))
    })
    output$qrep <- renderUI({
      pr <- rv$prim; req(pr, rv$ana)
      sel <- .pr_nz(input$rep_sel, ""); req(nzchar(sel))
      kind <- sub(":.*$", "", sel); i <- suppressWarnings(as.integer(sub("^.*:", "", sel))); req(!is.na(i))
      pm <- prm()
      blk <- function(k, label, cls) {
        st <- pr_struct_primer(pr$Sequence[k], pm)
        tagList(.pr_p3_head(paste(label, pr$Name[k]), cls, st$seq),
                .pr_p3_metrics(st, st$seq, NULL, nchar(st$seq), pr_tm(pr$Sequence[k], pm)[["mean"]]),
                tags$details(open = NA, tags$summary("Secondary structures"), .pr_struct_set(st)))
      }
      if (kind == "one") { req(i <= nrow(pr)); return(tagList(blk(i, "Primer", "p3-l"), helpText("Any / End / HP = Tm (C) of the most stable self-dimer, 3' self-dimer and hairpin (0 = none)."))) }
      x <- rv$pairs; req(x, i <= nrow(x))
      a <- match(x$`Primer 1`[i], pr$Name); b <- match(x$`Primer 2`[i], pr$Name); req(!is.na(a), !is.na(b))
      cr <- pr_struct_pair(pr$Sequence[a], pr$Sequence[b], pm)
      tagList(blk(a, "Left Primer", "p3-l"), blk(b, "Right Primer", "p3-r"),
              tags$div(class = "p3-head p3-pair", tags$b("Pair:"), .pr_metric("Delta Tm", sprintf("%.1f C", x$`Delta Tm`[i])),
                       .pr_metric("Any", sprintf("%.1f", .pr_tmv(cr$any))), .pr_metric("End", sprintf("%.1f", .pr_tmv(cr$end))),
                       .pr_metric("3' pairs", paste(x$`3' pairs (bp)`[i], "bp")), .pr_metric("Verdict", x$Verdict[i])),
              tags$details(open = NA, tags$summary("Cross dimer (pair) secondary structures"),
                           .pr_struct_block("Compl Any Secondary Structure", cr$any), .pr_struct_block("Compl End Secondary Structure", cr$end)),
              helpText("Any / End / HP = Tm (C) of the most stable self-dimer, 3' self-dimer and hairpin (0 = none); for the pair, of the cross dimer."))
    })
    observe({
      pr <- rv$prim; x <- rv$pairs
      if (is.null(pr) || !nrow(pr)) { updateSelectInput(session, "rep_sel", choices = character()); return() }
      ones <- setNames(paste0("one:", seq_len(nrow(pr))), pr$Name)
      prs <- if (is.null(x)) NULL else setNames(paste0("pair:", seq_len(nrow(x))), paste(x$`Primer 1`, "+", x$`Primer 2`))
      updateSelectInput(session, "rep_sel", choices = if (is.null(prs)) ones else list(`Pairs (with cross dimers)` = prs, `Single primers` = ones),
                        selected = if (is.null(prs)) "one:1" else "pair:1")
    })
    observeEvent(input$pair_tbl_rows_selected, { i <- input$pair_tbl_rows_selected; if (length(i)) updateSelectInput(session, "rep_sel", selected = paste0("pair:", i[1])) }, ignoreInit = TRUE)
    observeEvent(input$ana_tbl_rows_selected, { i <- input$ana_tbl_rows_selected; if (length(i)) updateSelectInput(session, "rep_sel", selected = paste0("one:", i[1])) }, ignoreInit = TRUE)
    bind_sum <- reactive({
      b <- rv$bind; req(b)
      s <- t(vapply(b, pr_bind_summary, numeric(6)))
      data.frame(Primer = rownames(s), s, check.names = FALSE, row.names = NULL, stringsAsFactors = FALSE)
    })
    output$bind_tbl <- DT::renderDT(DT::datatable(bind_sum(), rownames = FALSE, selection = list(mode = "single", selected = 1), options = list(dom = "t", scrollX = TRUE)))
    output$bind_detail <- DT::renderDT({
      b <- rv$bind; req(b); i <- input$bind_tbl_rows_selected; i <- if (length(i)) i[1] else 1L
      DT::datatable(b[[i]], rownames = FALSE, selection = "none", options = list(pageLength = 8, scrollX = TRUE, dom = "frtip"),
                    caption = paste("Best binding site per sequence:", names(b)[i]))
    })

    # ---------------------------------------------------------------------------------------------- 4. BLAST
    # the search itself runs on the NCBI website: the app only fills in the form with the primers, which the user can edit first
    pick_seq <- function(nm) { pr <- rv$prim; if (is.null(pr) || !nzchar(.pr_nz(nm, ""))) return(NULL); k <- match(nm, pr$Name); if (is.na(k)) NULL else .pr_clean(pr$Sequence[k]) }
    observeEvent(input$bl_fpick, { x <- pick_seq(input$bl_fpick); if (!is.null(x)) updateTextInput(session, "bl_fseq", value = x) }, ignoreInit = FALSE)
    observeEvent(input$bl_rpick, { x <- pick_seq(input$bl_rpick); if (!is.null(x)) updateTextInput(session, "bl_rseq", value = x) }, ignoreInit = FALSE)
    observeEvent(input$bl_fclear, { updateTextInput(session, "bl_fseq", value = ""); updateSelectizeInput(session, "bl_fpick", selected = "") })
    observeEvent(input$bl_rclear, { updateTextInput(session, "bl_rseq", value = ""); updateSelectizeInput(session, "bl_rpick", selected = "") })
    iu <- list(R = "AG", Y = "CT", S = "GC", W = "AT", K = "GT", M = "AC", B = "CGT", D = "AGT", H = "ACT", V = "ACG", N = "ACGT")
    clean_in <- function(x) toupper(gsub("[^A-Za-z]", "", .pr_nz(x, "")))
    variants_of <- function(sq, maxn = 64L) {
      out <- ""
      for (ch in strsplit(sq, "")[[1]]) { opt <- if (!is.null(iu[[ch]])) strsplit(iu[[ch]], "")[[1]] else ch; out <- as.vector(outer(out, opt, paste0)); if (length(out) > 4096L) break }
      list(v = utils::head(out, maxn), n = prod(vapply(strsplit(sq, "")[[1]], function(ch) if (!is.null(iu[[ch]])) nchar(iu[[ch]]) else 1L, 1L)))
    }
    amb_ui <- function(sq, id) {
      if (!nzchar(sq)) return(NULL)
      ch <- strsplit(sq, "")[[1]]; ab <- which(ch %in% names(iu)); bad <- which(!ch %in% c("A", "C", "G", "T", names(iu)))
      if (length(bad)) return(helpText(style = "color:#c0392b;", sprintf("Not a nucleotide letter: %s", paste(unique(ch[bad]), collapse = " "))))
      if (!length(ab)) return(helpText(sprintf("%d nt, no ambiguous base.", length(ch))))
      vr <- variants_of(sq)
      tagList(helpText(sprintf("%d nt; ambiguous: %s. It stands for %s sequence(s).", length(ch),
                               paste0(ch[ab], " (", vapply(ch[ab], function(z) paste(strsplit(iu[[z]], "")[[1]], collapse = "/"), ""), ") at ", ab, collapse = ", "), format(vr$n, big.mark = ","))),
              selectizeInput(ns(id), if (vr$n > 64) "Replace by one variant (first 64 shown; edit the letters above for others)" else "Replace by one variant",
                             choices = c("(keep as typed)" = "", vr$v), selected = ""))
    }
    output$bl_fvar <- renderUI(amb_ui(clean_in(input$bl_fseq), "bl_fv"))
    output$bl_rvar <- renderUI(amb_ui(clean_in(input$bl_rseq), "bl_rv"))
    observeEvent(input$bl_fv, { if (nzchar(.pr_nz(input$bl_fv, ""))) updateTextInput(session, "bl_fseq", value = input$bl_fv) }, ignoreInit = TRUE)
    observeEvent(input$bl_rv, { if (nzchar(.pr_nz(input$bl_rv, ""))) updateTextInput(session, "bl_rseq", value = input$bl_rv) }, ignoreInit = TRUE)
    output$bl_pb <- renderUI({
      fs <- clean_in(input$bl_fseq); rs <- clean_in(input$bl_rseq)
      if (!nzchar(fs) && !nzchar(rs)) return(helpText("Choose or type at least one primer."))
      enc <- function(x) utils::URLencode(as.character(x), reserved = TRUE)
      btn <- function(href, label) tags$a(href = href, target = "_blank", class = "btn btn-primary", icon("arrow-up-right-from-square"), paste0(" ", label))
      out <- list()
      if (nzchar(fs) && nzchar(rs)) {
        q <- c(PRIMER_LEFT_INPUT = fs, PRIMER_RIGHT_INPUT = rs)
        out[[1]] <- btn(paste0("https://www.ncbi.nlm.nih.gov/tools/primer-blast/index.cgi?", paste0(names(q), "=", vapply(q, enc, ""), collapse = "&")), "Open in NCBI Primer-BLAST (pair)")
      }
      nm <- c(if (nzchar(fs)) "forward", if (nzchar(rs)) "reverse"); sq <- c(fs, rs)[nzchar(c(fs, rs))]
      wq <- c(PAGE_TYPE = "BlastSearch", PROGRAM = "blastn", BLAST_PROGRAMS = "blastn", LINK_LOC = "blasthome", WORD_SIZE = 7,
              EXPECT = 1000, NUCL_REWARD = 1, NUCL_PENALTY = -3, QUERY = paste0(">", nm, "\n", sq, collapse = "\n"))
      out[[length(out) + 1L]] <- tags$a(href = paste0("https://blast.ncbi.nlm.nih.gov/Blast.cgi?", paste0(names(wq), "=", vapply(wq, enc, ""), collapse = "&")),
                                        target = "_blank", class = "btn btn-outline-primary", icon("arrow-up-right-from-square"), " Open in NCBI BLAST (each primer)")
      do.call(tagList, lapply(out, function(z) tagList(z, " ")))
    })
    output$bl_seq <- DT::renderDT({
      pr <- rv$prim
      if (is.null(pr) || !nrow(pr)) return(DT::datatable(data.frame(Note = "No primer yet: send primers from Design or the primer list to the Quality tab and press Analyse."), rownames = FALSE, options = list(dom = "t")))
      sq <- vapply(pr$Sequence, .pr_clean, "")
      ab <- vapply(strsplit(sq, ""), function(ch) paste(unique(ch[!ch %in% c("A", "C", "G", "T")]), collapse = " "), "")
      x <- data.frame(Name = pr$Name, `Sequence 5'->3'` = sq, `Reverse complement` = unname(.pr_rc(sq)), Length = nchar(sq),
                      `Ambiguous letters` = ab, `Degeneracy (sequences)` = vapply(sq, .pr_deg, 0), check.names = FALSE, stringsAsFactors = FALSE)
      DT::datatable(x, rownames = FALSE, selection = "none", options = list(pageLength = 10, scrollX = TRUE, dom = "tip"))
    })

    list(designed = reactive(rv$design), primers = reactive(rv$prim), pool = reactive(rv$pool))
  })
}
