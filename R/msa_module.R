# Sangerly - "MSA and trimming" tab
#   1. Multiple sequence alignment (DECIPHER, or the 'msa' package: ClustalOmega / ClustalW / Muscle)
#   2. Column trimming in the spirit of trimAl:
#        - built in (R):  gap threshold (-gt), conservation threshold (-st), both, ragged ends only,
#                         and a minimum percentage of columns to keep (-cons)
#        - the real trimAl program, if it is installed (automated1, gappyout, strict, strictplus, ...)
#
# Optional packages:   BiocManager::install(c("DECIPHER", "msa"))
# Optional program:    trimal (https://github.com/inab/trimal) on the PATH, or give its full path in the tab.

# ---------------------------------------------------------------------------------------------------
# Core functions (no Shiny)
# ---------------------------------------------------------------------------------------------------

# is a package installed? (does NOT load it: loading DECIPHER and msa takes several seconds at every start)
.pkg_installed <- function(p) nzchar(system.file(package = p))

.msa_methods <- function() {
  m <- character()
  if (.pkg_installed("DECIPHER")) m <- c(m, "DECIPHER (AlignSeqs)" = "decipher")
  if (.pkg_installed("msa"))
    m <- c(m, "Clustal Omega (msa)" = "ClustalOmega", "ClustalW (msa)" = "ClustalW", "MUSCLE (msa)" = "Muscle")
  m
}


# trimAl: a path typed by the user, else the copy shipped inside the app folder (trimAL/trimal.exe),
# else the PATH. The app is started from its own folder, so relative paths work.
.msa_find_trimal <- function(user = "") {
  user <- trimws(user)
  if (nzchar(user)) return(if (file.exists(user)) normalizePath(user, winslash = "/") else "")
  cand <- c(file.path("trimAL", "trimal.exe"), file.path("trimAL", "trimal"),
            file.path("trimal", "trimal.exe"), file.path("trimal", "trimal"),
            file.path("bin", "trimal.exe"), file.path("bin", "trimal"), "trimal.exe")
  cand <- cand[file.exists(cand)]
  if (length(cand)) return(normalizePath(cand[1], winslash = "/"))
  unname(Sys.which("trimal"))
}

# Numeric progress 0-100 (bar + "NN %"), closed with $close().
.msa_progress <- function(session, msg) {
  p <- shiny::Progress$new(session, min = 0, max = 100)
  p$set(0, message = msg, detail = "0 %")
  list(set = function(pct, detail = NULL)
         p$set(value = pct, message = msg,
               detail = paste0(round(pct), " %", if (!is.null(detail)) paste0(" - ", detail) else "")),
       close = function() p$close())
}


# Reverse-complement the sequences that share more 8-mers with the reverse complement of the
# reference (the longest sequence) than with the reference itself. Sequences longer than `maxlen` bp
# are left as they are (k-mer sets of whole genomes would be too heavy); that is only a convenience
# step, the alignment itself has no size limit here.
.msa_orient <- function(seqs, k = 8L, maxlen = 2e5, progress = NULL) {
  if (length(seqs) < 2) return(list(seqs = seqs, flipped = character()))
  rc <- function(s) as.character(Biostrings::reverseComplement(Biostrings::DNAString(s)))
  kmers <- function(s) {
    n <- nchar(s)
    if (n < k) return(character())
    unique(substring(s, 1:(n - k + 1L), k:n))
  }
  ok <- nchar(seqs) <= maxlen
  if (sum(ok) < 2) return(list(seqs = seqs, flipped = character()))
  ref <- which(ok)[which.max(nchar(seqs[ok]))]
  kr <- kmers(seqs[ref])
  flipped <- character()
  for (i in seq_along(seqs)) {
    if (!is.null(progress)) progress(i, length(seqs))
    if (i == ref || !ok[i]) next
    fw <- length(intersect(kmers(seqs[i]), kr))
    rv <- length(intersect(kmers(rc(seqs[i])), kr))
    if (rv > fw) { seqs[i] <- rc(seqs[i]); flipped <- c(flipped, names(seqs)[i]) }
  }
  list(seqs = seqs, flipped = flipped)
}

# seqs: named character vector (ACGT/IUPAC, no gaps). Returns a named character vector (aligned).
msa_align <- function(seqs, method = "decipher") {
  if (length(seqs) < 2) stop("At least 2 sequences are needed for an alignment.")
  dna <- Biostrings::DNAStringSet(seqs)
  names(dna) <- names(seqs)
  if (method == "decipher") {
    if (!requireNamespace("DECIPHER", quietly = TRUE))
      stop("The DECIPHER package is not installed: BiocManager::install(\"DECIPHER\")")
    res <- DECIPHER::AlignSeqs(dna, verbose = FALSE, processors = 1)
  } else {
    if (!requireNamespace("msa", quietly = TRUE))
      stop("The msa package is not installed: BiocManager::install(\"msa\")")
    res0 <- msa::msa(dna, method = method, type = "dna", verbose = FALSE)
    res <- Biostrings::unmasked(res0)
  }
  out <- as.character(res)
  names(out) <- names(seqs)
  out
}

# Character matrix [sequence x column] of an aligned named character vector.
msa_matrix <- function(aln) {
  m <- do.call(rbind, strsplit(toupper(unname(aln)), ""))
  rownames(m) <- names(aln)
  m
}

# Per-column statistics: gap fraction and conservation (frequency of the most common base among
# the sequences that have a base in that column).
msa_column_stats <- function(m) {
  isgap <- m == "-" | m == "." | m == "?"
  gapfrac <- colMeans(isgap)
  cons <- vapply(seq_len(ncol(m)), function(j) {
    b <- m[!isgap[, j], j]
    if (!length(b)) return(0)
    max(table(b)) / length(b)
  }, numeric(1))
  data.frame(Col = seq_len(ncol(m)), GapFrac = gapfrac, Cons = cons, stringsAsFactors = FALSE)
}

# Built-in trimming. mode: "none", "gt", "st", "both", "ends".
#   gt        keep columns where at least `gt` (0-1) of the sequences have a base       (trimAl -gt)
#   st        keep columns whose conservation is at least `st` (0-1)                   (like trimAl -st)
#   ends      remove the ragged ends: keep from the first to the last column where at least `ends_cov`
#             of the sequences have a base (inner gaps stay)
#   cons_min  minimum percentage (0-100) of the columns to keep; if the rule keeps fewer, the best
#             columns (most bases, then most conserved) are added back                  (trimAl -cons)
msa_trim_native <- function(m, mode = "gt", gt = 0.5, st = 0.5, ends_cov = 0.5, cons_min = 0) {
  s <- msa_column_stats(m)
  nc <- nrow(s)
  have <- 1 - s$GapFrac
  keep <- switch(mode,
    none = rep(TRUE, nc),
    gt   = have >= gt - 1e-12,
    st   = have > 0 & s$Cons >= st - 1e-12,
    both = have >= gt - 1e-12 & s$Cons >= st - 1e-12,
    ends = {
      w <- which(have >= ends_cov - 1e-12)
      if (length(w)) seq_len(nc) >= min(w) & seq_len(nc) <= max(w) else rep(FALSE, nc)
    },
    rep(TRUE, nc))
  keep <- keep & have > 0                     # a column with no base at all is always dropped
  if (!is.null(cons_min) && cons_min > 0) {
    need <- ceiling(nc * cons_min / 100)
    if (sum(keep) < need) {
      ord <- order(-have, -s$Cons, s$Col)
      keep[ord[seq_len(min(need, nc))]] <- TRUE
    }
  }
  list(kept = keep, stats = s)
}

# The real trimAl program. Returns list(kept = logical or NULL, trimmed = named character vector).
msa_trim_trimal <- function(aln, bin, flag = "-automated1", gt = NULL, st = NULL, cons = NULL) {
  fin <- tempfile(fileext = ".fasta"); fout <- tempfile(fileext = ".fasta")
  on.exit(unlink(c(fin, fout)), add = TRUE)
  writeLines(msa_fasta_text(aln), fin)
  args <- c("-in", shQuote(fin), "-out", shQuote(fout), "-fasta", "-colnumbering")
  if (identical(flag, "custom")) {
    if (!is.null(gt))   args <- c(args, "-gt", format(gt))
    if (!is.null(st))   args <- c(args, "-st", format(st))
    if (!is.null(cons) && cons > 0) args <- c(args, "-cons", format(cons))
  } else {
    args <- c(args, flag)
  }
  # system2() already quotes the command itself (do not shQuote it); the file arguments are quoted above
  out <- tryCatch(suppressWarnings(system2(bin, args, stdout = TRUE, stderr = TRUE)),
                  error = function(e) stop("trimAl could not be run: ", conditionMessage(e)))
  if (!file.exists(fout) || file.size(fout) == 0)
    stop("trimAl produced no output. Its message: ", paste(utils::head(out, 5), collapse = " "))
  recs <- blast_parse_records(paste(readLines(fout, warn = FALSE), collapse = "\n"))
  trimmed <- stats::setNames(recs$seq, recs$id)
  kept <- NULL
  i <- grep("ColumnsMap", out, fixed = TRUE)
  if (length(i)) {
    txt <- paste(out[i[1]:length(out)], collapse = " ")
    txt <- sub(".*ColumnsMap", "", txt)
    num <- regmatches(txt, gregexpr("[0-9]+", txt))[[1]]
    cols <- as.integer(num) + 1L                       # trimAl numbers the columns from 0
    nc <- nchar(aln[1])
    if (length(cols) && length(cols) == nchar(trimmed[1]) && all(cols >= 1L & cols <= nc))
      kept <- seq_len(nc) %in% cols
  }
  list(kept = kept, trimmed = trimmed)
}

msa_apply_mask <- function(aln, kept) {
  m <- msa_matrix(aln)
  out <- apply(m[, kept, drop = FALSE], 1, paste, collapse = "")
  stats::setNames(as.character(out), names(aln))
}

msa_fasta_text <- function(x, width = 60) {
  paste(vapply(seq_along(x), function(i) {
    s <- x[[i]]
    if (!nzchar(s)) return(paste0(">", names(x)[i]))
    ch <- substring(s, seq(1, nchar(s), by = width), pmin(seq(width, nchar(s) + width - 1, by = width), nchar(s)))
    paste(c(paste0(">", names(x)[i]), ch), collapse = "\n")
  }, ""), collapse = "\n")
}

.msa_style <- function() {
  shiny::tags$style(shiny::HTML("
    .msa-box{overflow-x:auto;background:#fff;border:1px solid #dfe3e8;border-radius:6px;padding:8px}
    .msa-box pre{margin:0 0 10px 0;font-family:Consolas,'Courier New',monospace;font-size:12px;line-height:1.35;background:transparent;border:0;padding:0}
    .msa-box .ruler{color:#7b8794}
    .mA{background:#c6f0d6}.mC{background:#cfe4fb}.mG{background:#e3e3e3}.mT{background:#ffd0d0}
    .mN{background:#f4f0c8}.mgap{color:#9aa5b1}.mrm{opacity:.28}
  "))
}

# Colour-coded alignment. kept = logical per column (columns removed by trimming are faded) or NULL.
msa_html <- function(m, kept = NULL, width = 80L, max_cols = 4000L) {
  nc <- ncol(m); shown <- min(nc, max_cols)
  nm <- rownames(m); nmw <- min(max(nchar(nm)), 28L)
  lab <- formatC(substr(nm, 1, nmw), width = -nmw)
  lab <- htmltools::htmlEscape(lab)
  cls <- function(b, k) {
    c1 <- ifelse(b %in% c("A", "C", "G", "T"), paste0("m", b), ifelse(b %in% c("-", ".", "?"), "mgap", "mN"))
    if (!is.null(k)) c1 <- ifelse(k, c1, paste(c1, "mrm"))
    c1
  }
  blocks <- character()
  for (s in seq(1L, shown, by = width)) {
    e <- min(s + width - 1L, shown)
    lines <- vapply(seq_len(nrow(m)), function(i) {
      b <- m[i, s:e]
      k <- if (is.null(kept)) NULL else kept[s:e]
      paste0(lab[i], " ", paste(sprintf("<span class='%s'>%s</span>", cls(b, k), b), collapse = ""))
    }, "")
    ruler <- sprintf("<span class='ruler'>%s %d</span>", strrep(" ", nmw), s)
    blocks <- c(blocks, paste0("<pre>", ruler, "\n", paste(lines, collapse = "\n"), "</pre>"))
  }
  note <- if (nc > shown) sprintf("<p class='text-muted'>Showing the first %d of %d columns (download the FASTA for the rest).</p>", shown, nc) else ""
  shiny::HTML(paste0("<div class='msa-box'>", paste(blocks, collapse = ""), note, "</div>"))
}

# ---------------------------------------------------------------------------------------------------
# Shiny module
# ---------------------------------------------------------------------------------------------------


# ---- sequence list: drag to reorder, right-click menu (rename, move, delete) -------------------------------
.msa_seqlist_head <- function() {
  shiny::tagList(
    shiny::tags$style(shiny::HTML("
      .sg-sl{max-height:560px;overflow:auto;border:1px solid #dfe3e8;border-radius:8px;background:#fff}
      .sg-sl-row{display:flex;gap:10px;align-items:center;padding:5px 10px;border-bottom:1px solid #eef1f5;
                 cursor:grab;-webkit-user-select:none;user-select:none;font-size:13px}
      .sg-sl-row:hover{background:#f4f8ff}
      .sg-sl-row.sel{background:#dcebff}
      .sg-sl-row.ov-t{box-shadow:inset 0 3px 0 #2c7be5}
      .sg-sl-row.ov-b{box-shadow:inset 0 -3px 0 #2c7be5}
      .sg-sl-n{width:36px;color:#7b8794;text-align:right}
      .sg-sl-name{flex:1;font-family:Consolas,'Courier New',monospace;overflow:hidden;text-overflow:ellipsis;white-space:nowrap}
      .sg-sl-bp{color:#52606d;width:90px;text-align:right}
      #sg-sl-menu{display:none;position:fixed;z-index:10000;background:#fff;border:1px solid #cfd6df;border-radius:8px;
                  box-shadow:0 6px 20px rgba(0,0,0,.2);padding:4px;min-width:180px}
      .sg-sl-mi{padding:6px 12px;cursor:pointer;border-radius:5px;font-size:13px}
      .sg-sl-mi:hover{background:#eaf2ff}.sg-sl-mi.danger{color:#c0392b}
    ")),
    shiny::tags$script(shiny::HTML("
(function(){
  if (window.__sgSeqList) return; window.__sgSeqList = true;
  var menu = null, drag = null, ctx = {box: null}, last = null;
  function rows(box){ return Array.prototype.slice.call(box.querySelectorAll('.sg-sl-row')); }
  function nm(list){ return list.map(function(r){ return r.getAttribute('data-name'); }); }
  function sel(box){ return rows(box).filter(function(r){ return r.classList.contains('sel'); }); }
  function up(el, cls){ return el && el.closest ? el.closest(cls) : null; }
  function send(id, obj){ obj.ts = Date.now(); Shiny.setInputValue(id, obj, {priority: 'event'}); }
  function hide(){ if (menu) menu.style.display = 'none'; }
  function build(){
    menu = document.createElement('div'); menu.id = 'sg-sl-menu';
    [['rename','Rename...'],['up','Move up'],['down','Move down'],['top','Move to top'],['bottom','Move to bottom'],['delete','Delete']]
      .forEach(function(d){
        var b = document.createElement('div'); b.className = 'sg-sl-mi' + (d[0] === 'delete' ? ' danger' : '');
        b.setAttribute('data-act', d[0]); b.textContent = d[1]; menu.appendChild(b); });
    document.body.appendChild(menu);
  }
  function renumber(box){ rows(box).forEach(function(r, i){ r.querySelector('.sg-sl-n').textContent = i + 1; }); }
  function clearMarks(){ document.querySelectorAll('.sg-sl-row.ov-t,.sg-sl-row.ov-b').forEach(function(r){ r.classList.remove('ov-t','ov-b'); }); }

  document.addEventListener('contextmenu', function(e){
    var row = up(e.target, '.sg-sl-row');
    if (!row) { hide(); return; }
    e.preventDefault();
    var box = up(row, '.sg-sl');
    if (!row.classList.contains('sel')) { rows(box).forEach(function(r){ r.classList.remove('sel'); }); row.classList.add('sel'); }
    ctx.box = box;
    if (!menu) build();
    menu.style.display = 'block';
    menu.style.left = Math.max(2, Math.min(e.clientX, window.innerWidth - menu.offsetWidth - 4)) + 'px';
    menu.style.top  = Math.max(2, Math.min(e.clientY, window.innerHeight - menu.offsetHeight - 4)) + 'px';
  });

  document.addEventListener('click', function(e){
    var mi = up(e.target, '.sg-sl-mi');
    if (mi && ctx.box) {
      var act = mi.getAttribute('data-act'), box = ctx.box, names = nm(sel(box));
      hide();
      if (!names.length) return;
      var id = box.getAttribute('data-act');
      if (act === 'rename') {
        var v = window.prompt('New name for \"' + names[0] + '\":', names[0]);
        if (v === null) return;
        send(id, {act: 'rename', names: [names[0]], value: v});
      } else if (act === 'delete') {
        if (!window.confirm('Delete ' + names.length + ' sequence(s)?')) return;
        send(id, {act: 'delete', names: names});
      } else send(id, {act: act, names: names});
      return;
    }
    hide();
    var row = up(e.target, '.sg-sl-row');
    if (!row) return;
    var box2 = up(row, '.sg-sl'), all = rows(box2);
    if (e.shiftKey && last && box2.contains(last)) {
      var a = all.indexOf(last), b = all.indexOf(row);
      all.forEach(function(r, i){ r.classList.toggle('sel', i >= Math.min(a, b) && i <= Math.max(a, b)); });
    } else if (e.ctrlKey || e.metaKey) { row.classList.toggle('sel'); last = row; }
    else { all.forEach(function(r){ r.classList.remove('sel'); }); row.classList.add('sel'); last = row; }
  });
  document.addEventListener('scroll', hide, true);

  document.addEventListener('dragstart', function(e){
    var row = up(e.target, '.sg-sl-row'); if (!row) return;
    drag = row;
    if (!row.classList.contains('sel')) { var b = up(row, '.sg-sl'); rows(b).forEach(function(r){ r.classList.remove('sel'); }); row.classList.add('sel'); }
    e.dataTransfer.effectAllowed = 'move';
    try { e.dataTransfer.setData('text/plain', row.getAttribute('data-name')); } catch (x) {}
  });
  document.addEventListener('dragover', function(e){
    if (!drag) return;
    var box = up(e.target, '.sg-sl'); if (!box) return;
    e.preventDefault(); clearMarks();
    var row = up(e.target, '.sg-sl-row');
    if (row) { var r = row.getBoundingClientRect(); row.classList.add(e.clientY > r.top + r.height / 2 ? 'ov-b' : 'ov-t'); }
  });
  document.addEventListener('drop', function(e){
    if (!drag) return;
    var box = up(e.target, '.sg-sl'); if (!box) return;
    e.preventDefault(); clearMarks();
    var row = up(e.target, '.sg-sl-row');
    var moving = sel(box); if (moving.indexOf(drag) < 0) moving = [drag];
    if (row && moving.indexOf(row) >= 0) { drag = null; return; }
    var ref = null;
    if (row) {
      var r = row.getBoundingClientRect();
      ref = e.clientY > r.top + r.height / 2 ? row.nextElementSibling : row;
      while (ref && moving.indexOf(ref) >= 0) ref = ref.nextElementSibling;
    }
    moving.forEach(function(m){ box.insertBefore(m, ref); });
    renumber(box);
    send(box.getAttribute('data-order'), {names: nm(rows(box))});
    drag = null;
  });
  document.addEventListener('dragend', function(){ clearMarks(); drag = null; });
})();
"))
  )
}

# ---- alignment editor (AliView-like): the JavaScript is read from www/sangerly_msa_editor.js and embedded inline ----
.msaed_js_tag <- function() {
  p <- c("www/sangerly_msa_editor.js", "sangerly_msa_editor.js", "R/sangerly_msa_editor.js", "../www/sangerly_msa_editor.js")
  p <- p[file.exists(p)]
  if (length(p)) {
    js <- paste(readLines(p[1], warn = FALSE, encoding = "UTF-8"), collapse = "\n")
    js <- gsub("</script", "<\\/script", js, fixed = TRUE)
    return(shiny::tags$script(shiny::HTML(js)))
  }
  message("Sangerly: sangerly_msa_editor.js not found (working directory: ", getwd(), "). Falling back to src=.")
  shiny::tags$script(src = "sangerly_msa_editor.js")
}

.msaed_style <- function() {
  shiny::tags$style(shiny::HTML("
    .msaed-root{background:#fff;border:1px solid #dfe3e8;border-radius:8px;margin-bottom:10px}
    .msaed-toolbar{display:flex;flex-wrap:wrap;gap:6px 10px;align-items:center;padding:6px 10px;border-bottom:1px solid #e3e8ee;font-size:12.5px;color:#1f3c5f}
    .msaed-toolbar label{display:flex;align-items:center;gap:6px;margin:0;font-weight:500}
    .msaed-toolbar input[type=range]{width:110px}
    .msaed-sep{width:1px;height:20px;background:#d5dbe3;margin:0 2px}
    .msaed-info{padding:4px 10px;font-size:12px;color:#52606d;border-bottom:1px solid #eef1f5;min-height:26px}
    .msaed-scroll{display:none;position:relative;height:480px;overflow:auto;outline:none}
    .msaed-scroll:focus{box-shadow:inset 0 0 0 2px rgba(44,123,229,.35)}
    .msaed-canvas{display:block;position:sticky;left:0;top:0}
    .msaed-spacer{position:absolute;left:0;top:0;height:1px;width:1px;pointer-events:none}
    .msaed-empty{padding:28px;text-align:center;color:#6b7785}
    .msaed-ovbox{display:none;border-top:1px solid #e3e8ee;background:#f5f7fa}
    .msaed-ov{display:block;cursor:pointer}
  "))
}

msaEditorUI <- function(id, readonly = FALSE, empty_text = "Align (or load an alignment) and it appears here.") {
  b <- function(act, label, title = NULL, cls = "btn-light")
    shiny::tags$button(type = "button", class = paste("btn btn-sm border", cls), `data-act` = act, title = title, label)
  sep <- function() shiny::span(class = "msaed-sep")
  bars <- shiny::tags$label("Bars",
    shiny::tags$select(class = "msaed-bars form-select form-select-sm", style = "width:auto",
      shiny::tags$option(value = "cons", "Conservation"), shiny::tags$option(value = "occ", "Bases (not gaps)"),
      shiny::tags$option(value = "none", "None")))
  top <- shiny::div(class = "msaed-toolbar",
    shiny::tags$label("Zoom", shiny::tags$input(type = "range", class = "msaed-zoom", min = 5, max = 24, step = 1, value = 12)),
    if (!readonly) shiny::tags$label(shiny::tags$input(type = "checkbox", class = "msaed-showrm", checked = NA), "Show removed columns (faded)"),
    shiny::tags$label(shiny::tags$input(type = "checkbox", class = "msaed-dots"), "Dots = same as reference (or consensus)"),
    bars,
    if (!readonly) shiny::tags$label("Delete key on bases",
      shiny::tags$select(class = "msaed-delmode form-select form-select-sm", style = "width:auto",
        shiny::tags$option(value = "col", "removes the column"), shiny::tags$option(value = "gap", "replaces by a gap"),
        shiny::tags$option(value = "shift", "deletes and shifts left"))))
  actions <- if (readonly)
    shiny::div(class = "msaed-toolbar", b("dl_kept", "FASTA (aligned)"), b("dl_ungapped", "FASTA (no gaps)"))
  else shiny::div(class = "msaed-toolbar",
      b("undo", "Undo", "Ctrl+Z"), b("redo", "Redo", "Ctrl+Y"), sep(),
      b("delcols", "Delete columns", "Leave the selected columns out of the result (they can be restored). Key: Delete on a column selection"),
      b("restore", "Restore columns"), b("restoreall", "Restore all"),
      b("allgap", "Remove all-gap columns"), sep(),
      b("gap", "Set to gap (-)", "Replace the selected bases by gaps. Key: Delete"),
      b("insgap", "Insert gap", "Insert gap(s) before the selection; the sequence shifts right. Key: Space"),
      b("shiftdel", "Delete + shift left", "Delete the selected bases and close the hole. Key: Backspace"),
      b("delrows", "Delete sequences"), sep(),
      b("setref", "Set reference", "Click a sequence name first, then press this: that sequence becomes the reference"),
      b("clearref", "Clear reference"),
      b("copyref", "Copy from reference (=)", "Copy the reference's bases onto the selected bases. Key: ="),
      b("nref", "N -> reference base", "Replace N by the reference base (A, C, G or T): in the selection, or in the whole alignment if nothing is selected"),
      sep(),
      b("commit", "Apply edits", "Send this alignment back to Sangerly (other tabs, downloads)", "btn-primary"),
      b("dl_all", "FASTA (all columns)"), b("dl_kept", "FASTA (kept columns)"), b("dl_ungapped", "FASTA (kept, no gaps)"))
  shiny::div(id = id, class = "msaed-root", `data-readonly` = if (readonly) "1" else "0",
    top, actions,
    shiny::div(class = "msaed-info"),
    shiny::div(class = "msaed-scroll", tabindex = 0, shiny::tags$canvas(class = "msaed-canvas"), shiny::div(class = "msaed-spacer")),
    shiny::div(class = "msaed-ovbox", title = "Overview: click or drag to move along the alignment",
               shiny::tags$canvas(class = "msaed-ov")),
    shiny::div(class = "msaed-empty", empty_text))
}

msaUI <- function(id) {
  ns <- shiny::NS(id)
  meth <- .msa_methods()
  shiny::div(class = "container-fluid",
    .msa_style(), .msaed_style(), .msaed_js_tag(), .msa_seqlist_head(),
    shiny::fluidRow(
      shiny::column(4,
        .card("1. Sequences",
          shiny::fileInput(ns("fa_files"), "Add FASTA files (select one or many files at once; repeat as often as you like)",
                           multiple = TRUE, buttonLabel = "Add files...", placeholder = "",
                           accept = c(".fa", ".fasta", ".fas", ".fna", ".ffn", ".txt")),
          if (.sv_local()) {
            if (.Platform$OS.type == "windows")
              shiny::div(style = "display:flex;gap:8px;align-items:center;flex-wrap:wrap;margin-bottom:6px",
                shiny::actionButton(ns("browse_files"), "Browse files...", icon = shiny::icon("file")),
                shiny::actionButton(ns("browse_dir"), "Browse folder...", icon = shiny::icon("folder-open")),
                shiny::checkboxInput(ns("subfolders"), "with subfolders", TRUE))
            else shiny::tagList(
              shiny::textInput(ns("folder"), "Or add a whole folder (path on this computer)", "",
                               placeholder = "/home/me/data"),
              shiny::div(style = "display:flex;gap:10px;align-items:center;flex-wrap:wrap",
                shiny::actionButton(ns("add_folder"), "Add folder", icon = shiny::icon("folder-open")),
                shiny::checkboxInput(ns("subfolders"), "with subfolders", TRUE)))
          },
          shiny::div(style = "margin-top:6px",
            shiny::actionButton(ns("use_asm"), "Add consensus from the Assembly", icon = shiny::icon("arrow-down"))),
          shiny::textAreaInput(ns("seqs"), "Or paste FASTA", rows = 3, width = "100%",
                               placeholder = ">sample1\nACGT...\n>sample2\nACGT..."),
          shiny::actionButton(ns("add_paste"), "Add pasted sequences", icon = shiny::icon("plus")),
          shiny::hr(),
          shiny::uiOutput(ns("pool_count")),
          DT::dataTableOutput(ns("pool_tbl")),
          shiny::div(style = "display:flex;gap:6px;margin-top:8px",
            shiny::actionButton(ns("pool_rm"), "Remove selected", icon = shiny::icon("minus")),
            shiny::actionButton(ns("pool_clear"), "Clear all")),
          shiny::checkboxInput(ns("blast_replace"), "'Send to MSA' from BLAST replaces the current list (untick to add to it)", TRUE),
          shiny::helpText("A file may hold one or many sequences. There is no size limit; every sequence in this list is aligned.")),
        .card("2. Multiple sequence alignment",
          if (length(meth)) shiny::selectInput(ns("method"), "Method", meth)
          else shiny::div(class = "alert alert-warning",
                          "No alignment package found. Install one: ",
                          shiny::tags$code("BiocManager::install(\"DECIPHER\")"), " (or ", shiny::tags$code("\"msa\""), ")."),
          shiny::checkboxInput(ns("orient"), "Reverse-complement sequences that fit better in the opposite orientation", TRUE),
          shiny::checkboxInput(ns("prealigned"), "My sequences are already aligned (skip the MSA, only trim)", FALSE),
          shiny::actionButton(ns("run"), "Align", icon = shiny::icon("play"), class = "btn-primary")),
        .card("3. Trimming",
          shiny::radioButtons(ns("engine"), "Trimming engine",
            c("Built in (R): trimAl-style rules" = "native", "trimAl program (if installed)" = "trimal"), "native"),
          shiny::conditionalPanel(sprintf("input['%s'] == 'native'", ns("engine")),
            shiny::selectInput(ns("mode"), "Rule",
              c("Gap threshold (-gt)" = "gt", "Conservation threshold (-st)" = "st",
                "Gap AND conservation" = "both", "Ragged ends only" = "ends", "No trimming" = "none"), "gt"),
            shiny::conditionalPanel(sprintf("input['%s'] == 'gt' || input['%s'] == 'both'", ns("mode"), ns("mode")),
              shiny::sliderInput(ns("gt"), "Minimum fraction of sequences with a base in the column", 0, 1, 0.5, 0.05)),
            shiny::conditionalPanel(sprintf("input['%s'] == 'st' || input['%s'] == 'both'", ns("mode"), ns("mode")),
              shiny::sliderInput(ns("st"), "Minimum conservation of the column (most common base)", 0, 1, 0.5, 0.05)),
            shiny::conditionalPanel(sprintf("input['%s'] == 'ends'", ns("mode")),
              shiny::sliderInput(ns("ends_cov"), "Ends: first/last column where this fraction of sequences has a base", 0, 1, 0.5, 0.05)),
            shiny::sliderInput(ns("cons"), "Always keep at least this % of the columns (-cons)", 0, 100, 0, 5)),
          shiny::conditionalPanel(sprintf("input['%s'] == 'trimal'", ns("engine")),
            shiny::selectInput(ns("tflag"), "trimAl mode",
              c("automated1 (recommended for ML trees)" = "-automated1", "gappyout" = "-gappyout", "strict" = "-strict",
                "strictplus" = "-strictplus", "nogaps (remove all columns with gaps)" = "-nogaps",
                "noallgaps (remove only all-gap columns)" = "-noallgaps", "custom (-gt / -st / -cons from the sliders)" = "custom")),
            shiny::conditionalPanel(sprintf("input['%s'] == 'custom'", ns("tflag")),
              shiny::sliderInput(ns("t_gt"), "-gt", 0, 1, 0.5, 0.05),
              shiny::sliderInput(ns("t_st"), "-st (0 = not used)", 0, 1, 0, 0.05),
              shiny::sliderInput(ns("t_cons"), "-cons (%)", 0, 100, 0, 5)),
            shiny::uiOutput(ns("tpath_ui")),     # only shown when trimAl is not found automatically
            shiny::uiOutput(ns("trimal_status"))),
          shiny::actionButton(ns("trim_go"), "Apply trimming", icon = shiny::icon("scissors")))),
      shiny::column(8,
        shiny::uiOutput(ns("msg")),
        shiny::uiOutput(ns("boxes")),
        shiny::tabsetPanel(id = ns("tabs"),
          shiny::tabPanel("Editor", shiny::br(),
            shiny::div(style = "display:flex;gap:10px;align-items:center;flex-wrap:wrap;margin-bottom:6px",
              shiny::actionButton(ns("ed_reset"), "Reload from alignment + trimming (discards manual edits)", icon = shiny::icon("rotate-left")),
              shiny::helpText(style = "margin:0",
                "Click a base, a column (ruler) or a name; type a letter to change bases; Delete / Backspace / Space as in the buttons. ",
                "Removed columns stay in the file, faded, and can be restored. Press 'Apply edits' to use the result in the other tabs.")),
            msaEditorUI(ns("editor"))),
          shiny::tabPanel("Sequences", shiny::br(),
            shiny::div(style = "display:flex;gap:8px;align-items:center;flex-wrap:wrap;margin-bottom:6px",
              shiny::actionButton(ns("sl_az"), "Sort A-Z", icon = shiny::icon("arrow-down-a-z")),
              shiny::actionButton(ns("sl_rev"), "Reverse order", icon = shiny::icon("arrows-up-down")),
              shiny::helpText(style = "margin:0",
                "Drag a row to move it (several rows: click, Ctrl+click or Shift+click first). Right-click for Rename, Move and Delete. ",
                "The alignment, the trimmed result and the Phylogeny tab use this order. Press 'Apply edits' in the Editor before changing the order, ",
                "otherwise unapplied edits are lost.")),
            shiny::uiOutput(ns("seqlist_ui"))),
          shiny::tabPanel("Trimmed result", shiny::br(),
            shiny::helpText("The alignment that remains after the trimming. You can edit it here too (change an N, delete a column ",
                            "or a sequence, ...): press 'Apply edits' to make it the result used by the Phylogeny tab and the downloads. ",
                            "Applying edits in the Editor tab, or trimming again, rebuilds this result and replaces edits made here."),
            msaEditorUI(ns("trimview"), empty_text = "No trimmed result yet: align, then trim.")),
          shiny::tabPanel("Site statistics", shiny::br(),
            shiny::fluidRow(
              shiny::column(5, shiny::radioButtons(ns("ss_src"), "Alignment", inline = TRUE,
                                                   c("Trimmed result" = "trim", "Whole alignment" = "all"))),
              shiny::column(7, shiny::checkboxInput(ns("ss_cd"),
                "Leave out the sites that have a gap or an ambiguous base (complete deletion, as in MEGA)", TRUE))),
            shiny::uiOutput(ns("ss_note")),
            DT::DTOutput(ns("ss_sum")),
            shiny::plotOutput(ns("ss_map"), height = "170px"),
            shiny::div(style = "display:flex;gap:6px;flex-wrap:wrap;margin-top:8px",
              shiny::downloadButton(ns("dl_ss_xlsx"), "Excel (summary + every site)"),
              shiny::downloadButton(ns("dl_ss_csv"), "Sites (CSV)")),
            shiny::helpText("Conserved: one base type in the site. Variable: two or more. Parsimony-informative: at least two base types, ",
                            "each present in at least two sequences. Singleton: variable but not parsimony-informative.")),
          shiny::tabPanel("Columns", shiny::br(), shiny::plotOutput(ns("col_plot"), height = "340px"),
                          shiny::helpText("Blue = kept, grey = removed. Conservation = frequency of the most common base ",
                                          "among the sequences that have a base in that column."))),
        shiny::br(),
        shiny::div(style = "display:flex;gap:6px;flex-wrap:wrap",
          shiny::downloadButton(ns("dl_aln"), "Alignment (FASTA)"),
          shiny::downloadButton(ns("dl_trim"), "Trimmed alignment (FASTA)"),
          shiny::downloadButton(ns("dl_ungap"), "Trimmed, gaps removed (FASTA)"),
          shiny::downloadButton(ns("dl_cols"), "Columns (CSV)")))))
}

msaServer <- function(id, consensus = shiny::reactive(NULL)) {
  shiny::moduleServer(id, function(input, output, session) {
    rv <- shiny::reactiveValues(aln = NULL, kept = NULL, trimmed = NULL, msg = NULL, info = NULL, ver = 0L,
                                pool = data.frame(id = character(), source = character(), bp = integer(),
                                                  seq = character(), stringsAsFactors = FALSE))

    set_msg <- function(txt, type = "info") rv$msg <- if (is.null(txt)) NULL else list(txt = txt, type = type)

    output$tpath_ui <- shiny::renderUI({
      if (nzchar(.msa_find_trimal(""))) NULL      # bundled or on the PATH: no need to ask for a path
      else shiny::textInput(session$ns("tpath"), "Path to trimal (not found automatically)", "", placeholder = "C:/tools/trimal.exe")
    })
    trimal_bin <- shiny::reactive(.msa_find_trimal(.nz(input$tpath, "")))
    output$trimal_status <- shiny::renderUI({
      if (nzchar(trimal_bin())) shiny::helpText(paste("trimAl found:", trimal_bin()))
      else shiny::div(class = "alert alert-warning",
                      "trimAl was not found. Put trimal.exe in the app folder (trimAL/trimal.exe), add it to the PATH, ",
                      "or type its full path above. Until then use the built-in rules.")
    })


    # ---- reset / reorder / rename / delete (used by the buttons and by the Sequences tab) -------------------
    reset_alignment <- function() {
      rv$aln <- NULL; rv$kept <- NULL; rv$trimmed <- NULL; rv$info <- NULL; rv$msg <- NULL
      session$sendCustomMessage("msaed_clear", list(id = session$ns("editor")))
      session$sendCustomMessage("msaed_clear", list(id = session$ns("trimview")))
      invisible()
    }
    bump <- function() if (!is.null(rv$aln)) rv$ver <- rv$ver + 1L     # reload the editors
    cur_names <- function() if (!is.null(rv$aln)) names(rv$aln) else rv$pool$id
    reorder_named <- function(x, ord) {
      if (is.null(x)) return(x)
      x[c(ord[ord %in% names(x)], setdiff(names(x), ord))]
    }
    apply_order <- function(ord) {
      ord <- as.character(ord)
      rv$aln <- reorder_named(rv$aln, ord); rv$trimmed <- reorder_named(rv$trimmed, ord)
      if (nrow(rv$pool)) {
        ix <- match(ord, rv$pool$id); ix <- ix[!is.na(ix)]
        rv$pool <- rv$pool[c(ix, setdiff(seq_len(nrow(rv$pool)), ix)), , drop = FALSE]
      }
      bump()
    }
    rename_one <- function(old, new) {
      new <- trimws(new)
      if (!nzchar(new) || identical(new, old)) return(invisible())
      if (new %in% setdiff(unique(c(names(rv$aln), rv$pool$id)), old)) {
        shiny::showNotification(sprintf("The name '%s' is already used by another sequence.", new), type = "warning"); return(invisible())
      }
      fix <- function(x) { if (!is.null(x)) names(x)[names(x) == old] <- new; x }
      rv$aln <- fix(rv$aln); rv$trimmed <- fix(rv$trimmed)
      p <- rv$pool; p$id[p$id == old] <- new; rv$pool <- p
      bump()
    }
    delete_names <- function(ids) {
      rv$pool <- rv$pool[!rv$pool$id %in% ids, , drop = FALSE]
      if (!is.null(rv$aln)) {
        k <- !names(rv$aln) %in% ids
        if (!any(k)) { reset_alignment(); return(invisible()) }
        if (!all(k)) {
          rv$aln <- rv$aln[k]
          if (!is.null(rv$trimmed)) { kt <- !names(rv$trimmed) %in% ids; rv$trimmed <- if (any(kt)) rv$trimmed[kt] else NULL }
          bump()
        }
      }
      invisible()
    }
    output$seqlist_ui <- shiny::renderUI({
      nm <- cur_names()
      if (!length(nm)) return(shiny::helpText("No sequences yet. Add sequences, then press Align."))
      bp <- if (!is.null(rv$aln)) nchar(gsub("-", "", rv$aln)) else rv$pool$bp
      shiny::tagList(
        shiny::helpText(if (!is.null(rv$aln)) sprintf("%d sequences of the current alignment.", length(nm))
                        else sprintf("%d sequences in the list (not aligned yet).", length(nm))),
        shiny::div(class = "sg-sl", `data-order` = session$ns("seq_order"), `data-act` = session$ns("seq_act"),
          lapply(seq_along(nm), function(i)
            shiny::div(class = "sg-sl-row", draggable = "true", `data-name` = nm[i],
              shiny::span(class = "sg-sl-n", i), shiny::span(class = "sg-sl-name", title = nm[i], nm[i]),
              shiny::span(class = "sg-sl-bp", paste(format(bp[i], big.mark = ","), "bp"))))))
    })
    shiny::observeEvent(input$seq_order, apply_order(unlist(input$seq_order$names)))
    shiny::observeEvent(input$sl_az, { cur <- cur_names(); if (length(cur)) apply_order(cur[order(tolower(cur))]) })
    shiny::observeEvent(input$sl_rev, { cur <- cur_names(); if (length(cur)) apply_order(rev(cur)) })
    shiny::observeEvent(input$seq_act, {
      a <- input$seq_act; cur <- cur_names()
      nm <- as.character(unlist(a$names)); nm <- nm[nm %in% cur]
      if (!length(nm)) return()
      switch(as.character(a$act),
        rename = rename_one(nm[1], as.character(a$value)),
        delete = delete_names(nm),
        top    = apply_order(c(nm, setdiff(cur, nm))),
        bottom = apply_order(c(setdiff(cur, nm), nm)),
        up     = { o <- cur; for (n in nm) { i <- match(n, o); if (i > 1 && !(o[i - 1] %in% nm)) o[c(i - 1, i)] <- o[c(i, i - 1)] }; apply_order(o) },
        down   = { o <- cur; for (n in rev(nm)) { i <- match(n, o); if (i < length(o) && !(o[i + 1] %in% nm)) o[c(i, i + 1)] <- o[c(i + 1, i)] }; apply_order(o) })
    })

    # ---- sequence list (accumulates; every row is aligned) ----------------------------------------
    add_records <- function(recs, source) {
      if (is.null(recs) || !nrow(recs)) return(0L)
      recs <- recs[nzchar(gsub("-", "", recs$seq)), , drop = FALSE]
      if (!nrow(recs)) return(0L)
      ids <- make.unique(c(rv$pool$id, recs$id), sep = "_")
      ids <- ids[length(rv$pool$id) + seq_len(nrow(recs))]    # (x[-integer(0)] would be empty)
      nd <- data.frame(id = ids, source = source, bp = nchar(gsub("-", "", recs$seq)),
                       seq = recs$seq, stringsAsFactors = FALSE)
      rv$pool <- rbind(rv$pool, nd)
      nrow(nd)
    }
    read_fasta_file <- function(path, label) {
      recs <- blast_parse_records(paste(readLines(path, warn = FALSE), collapse = "\n"))
      if (nrow(recs) == 1 && identical(recs$header, "query"))      # no ">" line: name it after the file
        recs$id <- gsub("[^A-Za-z0-9._-]", "_", tools::file_path_sans_ext(label))
      recs
    }
    add_files <- function(paths, labels) {
      prog <- .msa_progress(session, "Reading sequences")
      on.exit(prog$close(), add = TRUE)
      n <- 0L
      for (i in seq_along(paths)) {
        prog$set(100 * (i - 1) / length(paths), sprintf("%d/%d %s", i, length(paths), labels[i]))
        r <- tryCatch(read_fasta_file(paths[i], labels[i]), error = function(e) NULL)
        if (is.null(r)) { shiny::showNotification(paste("Could not read", labels[i]), type = "warning"); next }
        n <- n + add_records(r, labels[i])
      }
      prog$set(100)
      shiny::showNotification(sprintf("%d sequence(s) added from %d file(s). %d in the list.",
                                      n, length(paths), nrow(rv$pool)), type = "message")
    }
    shiny::observeEvent(input$fa_files, {
      f <- input$fa_files
      add_files(f$datapath, f$name)
    })
    add_dir <- function(d) {
      d <- trimws(gsub("^\"|\"$", "", d))
      if (!nzchar(d) || !dir.exists(d)) { shiny::showNotification("That folder does not exist.", type = "warning"); return() }
      fl <- list.files(d, pattern = "\\.(fa|fasta|fas|fna|ffn|txt)$", recursive = isTRUE(input$subfolders),
                       full.names = TRUE, ignore.case = TRUE)
      if (!length(fl)) { shiny::showNotification("No FASTA files (.fa .fasta .fas .fna .ffn .txt) in that folder.", type = "warning"); return() }
      add_files(fl, basename(fl))
    }
    shiny::observeEvent(input$add_folder, add_dir(.nz(input$folder, "")))
    shiny::observeEvent(input$browse_dir, {
      d <- tryCatch(utils::choose.dir(caption = "Choose the folder with the FASTA files"), error = function(e) NA_character_)
      if (is.na(d) || !nzchar(d)) return()
      add_dir(d)
    })
    shiny::observeEvent(input$browse_files, {
      p <- tryCatch(utils::choose.files(
        filters = matrix(c("FASTA files (*.fa, *.fasta, *.fas, *.fna, *.ffn, *.txt)", "*.fa;*.fasta;*.fas;*.fna;*.ffn;*.txt",
                           "All files (*.*)", "*.*"), 2, byrow = TRUE), multi = TRUE),
        error = function(e) character())
      if (!length(p)) return()
      p <- normalizePath(p, winslash = "/")
      add_files(p, basename(p))
    })
    shiny::observeEvent(input$add_paste, {
      txt <- .nz(input$seqs, "")
      if (!nzchar(trimws(txt))) { shiny::showNotification("Paste some FASTA first.", type = "warning"); return() }
      n <- add_records(blast_parse_records(txt), "pasted")
      shiny::updateTextAreaInput(session, "seqs", value = "")
      shiny::showNotification(sprintf("%d sequence(s) added. %d in the list.", n, nrow(rv$pool)), type = "message")
    })
    shiny::observeEvent(input$use_asm, {
      x <- consensus()
      if (is.null(x) || !length(x)) {
        shiny::showNotification("The Assembly tab has no consensus yet. Assemble first.", type = "warning"); return()
      }
      txt <- if (length(x) == 1 && grepl("^>", x[[1]])) x[[1]]
             else paste0(">", .nz(names(x), paste0("seq", seq_along(x))), "\n", gsub("-", "", x), collapse = "\n")
      n <- add_records(blast_parse_records(txt), "Assembly")
      shiny::showNotification(sprintf("%d consensus sequence(s) added. %d in the list.", n, nrow(rv$pool)), type = "message")
    })
    output$pool_count <- shiny::renderUI({
      n <- nrow(rv$pool)
      shiny::helpText(if (!n) "No sequences yet." else sprintf("%d sequence(s), %s bp in total.", n, format(sum(rv$pool$bp), big.mark = ",")))
    })
    output$pool_tbl <- DT::renderDataTable({
      d <- rv$pool[, c("id", "source", "bp"), drop = FALSE]
      DT::datatable(d, rownames = FALSE, selection = "multiple",
                    options = list(pageLength = 8, scrollX = TRUE, dom = "ftp"))
    })
    shiny::observeEvent(input$pool_rm, {
      s <- input$pool_tbl_rows_selected
      if (!length(s)) { shiny::showNotification("Select rows in the list first.", type = "message"); return() }
      delete_names(rv$pool$id[s])          # also leaves them out of the alignment, so nothing stays behind
    })
    shiny::observeEvent(input$pool_clear, {
      rv$pool <- rv$pool[0, , drop = FALSE]
      reset_alignment()                    # the alignment, the trimmed result and both editors are emptied too
      shiny::showNotification("List, alignment and trimmed result were cleared.", type = "message")
    })

    # ---- trimming ----------------------------------------------------------------------------------
    do_trim <- function(aln) {
      if (is.null(aln)) return(invisible())
      rv$ver <- rv$ver + 1L                       # the Editor tab reloads whenever the alignment / trimming changes
      if (identical(input$engine, "trimal")) {
        bin <- trimal_bin()
        if (!nzchar(bin)) { set_msg("trimAl was not found (see the status under the trimming settings).", "warning"); return(invisible()) }
        r <- tryCatch(msa_trim_trimal(aln, bin, .nz(input$tflag, "-automated1"),
                                      gt = input$t_gt, st = if (isTRUE(input$t_st > 0)) input$t_st else NULL,
                                      cons = input$t_cons),
                      error = function(e) { set_msg(conditionMessage(e), "danger"); NULL })
        if (is.null(r)) return(invisible())
        rv$trimmed <- r$trimmed
        rv$kept <- r$kept
        if (is.null(r$kept)) set_msg("trimAl ran, but its column list could not be read, so the removed columns are not highlighted.", "info")
        else set_msg(NULL)
      } else {
        m <- msa_matrix(aln)
        r <- msa_trim_native(m, .nz(input$mode, "gt"), .nz(input$gt, 0.5), .nz(input$st, 0.5),
                             .nz(input$ends_cov, 0.5), .nz(input$cons, 0))
        rv$kept <- r$kept
        rv$trimmed <- if (any(r$kept)) msa_apply_mask(aln, r$kept) else NULL
        if (!any(r$kept)) set_msg("These settings remove every column. Relax the thresholds.", "warning") else set_msg(NULL)
      }
      rv$trim_meta <- list(engine = .nz(input$engine, "native"), mode = .nz(input$mode, "gt"), gt = input$gt, st = input$st,
                           ends_cov = input$ends_cov, cons = input$cons, tflag = input$tflag, t_gt = input$t_gt,
                           t_st = input$t_st, t_cons = input$t_cons)
      invisible()
    }

    # ---- alignment (numeric progress 0-100) -------------------------------------------------------
    run_align <- function() {
      pool <- rv$pool
      if (nrow(pool) < 2) { set_msg("Add at least 2 sequences first.", "warning"); return(invisible()) }
      prog <- .msa_progress(session, "Aligning")
      on.exit(prog$close(), add = TRUE)
      seqs <- stats::setNames(pool$seq, pool$id)
      if (isTRUE(input$prealigned)) {
        prog$set(10, "checking the alignment")
        seqs <- gsub("[^A-Z.?-]", "", toupper(seqs)); seqs <- gsub("[.?]", "-", seqs)
        if (length(unique(nchar(seqs))) != 1) {
          set_msg("The sequences are not all the same length, so they are not an alignment. Untick the 'already aligned' box to align them.", "danger")
          return(invisible())
        }
        rv$meta <- list(prealigned = TRUE, n = length(seqs))
        rv$aln <- seqs; rv$info <- NULL; rv$kept <- NULL; rv$trimmed <- NULL
        prog$set(60, "trimming"); do_trim(rv$aln); prog$set(100, "done")
        return(invisible())
      }
      prog$set(3, "preparing")
      seqs <- gsub("U", "T", gsub("[^A-Z]", "", toupper(gsub("-", "", seqs))))
      seqs <- seqs[nzchar(seqs)]
      if (length(seqs) < 2) { set_msg("At least 2 non-empty sequences are needed.", "warning"); return(invisible()) }
      if (!length(.msa_methods())) { set_msg("No alignment package is installed (DECIPHER or msa).", "danger"); return(invisible()) }
      method <- .nz(input$method, .msa_methods()[[1]])
      flipped <- character()
      if (isTRUE(input$orient)) {
        o <- .msa_orient(seqs, progress = function(i, n) prog$set(5 + 25 * i / n, sprintf("orienting %d/%d", i, n)))
        seqs <- o$seqs; flipped <- o$flipped
      }
      mm <- .msa_methods()
      rv$meta <- list(prealigned = FALSE, method = method, label = names(mm)[match(method, mm)],
                      orient = isTRUE(input$orient), n = length(seqs), flipped = length(flipped))
      prog$set(30, sprintf("aligning %d sequences (this step reports no finer progress)", length(seqs)))
      res <- tryCatch(msa_align(seqs, method), error = function(e) { set_msg(conditionMessage(e), "danger"); NULL })
      if (is.null(res)) return(invisible())
      prog$set(90, "trimming")
      rv$aln <- res; rv$kept <- NULL; rv$trimmed <- NULL
      rv$info <- if (length(flipped)) sprintf("Reverse-complemented: %s.", paste(flipped, collapse = ", ")) else NULL
      set_msg(NULL)
      do_trim(rv$aln)
      prog$set(100, "done")
      invisible()
    }
    shiny::observeEvent(input$run, run_align())

    # ---- alignment editor (JavaScript, www/sangerly_msa_editor.js) --------------------------------
    push_editor <- function() {
      aln <- rv$aln
      if (is.null(aln)) return(invisible())
      kept <- rv$kept
      if (!is.null(kept) && length(kept) != nchar(aln[1])) kept <- NULL
      if (is.null(kept) && !is.null(rv$trimmed)) aln <- rv$trimmed     # trimAl gave no column list: edit its result
      session$sendCustomMessage("msaed_load", list(
        id = session$ns("editor"), names = unname(names(aln)), seqs = unname(as.character(aln)),
        keep = if (is.null(kept)) NULL else paste(as.integer(kept), collapse = "")))
      invisible()
    }
    push_trim <- function() {
      tr <- rv$trimmed
      if (is.null(tr)) session$sendCustomMessage("msaed_clear", list(id = session$ns("trimview")))
      else session$sendCustomMessage("msaed_load", list(id = session$ns("trimview"), names = unname(names(tr)),
                                                        seqs = unname(as.character(tr)), keep = NULL))
      invisible()
    }
    shiny::observeEvent(rv$ver, { push_editor(); push_trim() }, ignoreInit = TRUE)
    shiny::observeEvent(input$ed_reset, { push_editor(); push_trim() })
    shiny::observeEvent(input$trimview_state, {
      s <- input$trimview_state
      nm <- as.character(unlist(s$names)); sq <- as.character(unlist(s$seqs))
      if (!length(nm) || length(nm) != length(sq)) return()
      aln <- stats::setNames(sq, nm)
      kept <- strsplit(as.character(.nz(s$keep, "")), "")[[1]] == "1"
      if (length(kept) != nchar(aln[1])) kept <- rep(TRUE, nchar(aln[1]))
      if (!any(kept)) { shiny::showNotification("Every column was removed; nothing was applied.", type = "warning"); return() }
      rv$trimmed <- msa_apply_mask(aln, kept)
      shiny::showNotification(sprintf("Trimmed result updated: %d sequences, %d columns.", length(aln), sum(kept)), type = "message")
    })
    shiny::observeEvent(input$editor_state, {
      s <- input$editor_state
      nm <- as.character(unlist(s$names)); sq <- as.character(unlist(s$seqs))
      if (!length(nm) || length(nm) != length(sq)) return()
      aln <- stats::setNames(sq, nm)
      kept <- strsplit(as.character(.nz(s$keep, "")), "")[[1]] == "1"
      if (length(kept) != nchar(aln[1])) kept <- rep(TRUE, nchar(aln[1]))
      rv$aln <- aln; rv$kept <- kept
      rv$trimmed <- if (any(kept)) msa_apply_mask(aln, kept) else NULL
      set_msg(NULL)
      push_trim()
      shiny::showNotification(sprintf("Edited alignment applied: %d sequences, %d of %d columns kept.",
                                      length(aln), sum(kept), length(kept)), type = "message")
    })
    shiny::observeEvent(input$trim_go, {
      if (is.null(rv$aln)) { set_msg("Align (or load an alignment) first.", "warning"); return() }
      do_trim(rv$aln)
    })

    # ---- outputs -----------------------------------------------------------------------------------
    output$msg <- shiny::renderUI({
      parts <- list()
      if (!is.null(rv$msg)) parts <- c(parts, list(shiny::div(class = paste0("alert alert-", rv$msg$type), rv$msg$txt)))
      if (!is.null(rv$info)) parts <- c(parts, list(shiny::div(class = "alert alert-info", rv$info)))
      do.call(shiny::tagList, parts)
    })
    output$boxes <- shiny::renderUI({
      shiny::req(rv$aln)
      nc <- nchar(rv$aln[1]); nk <- if (is.null(rv$trimmed)) NA else nchar(rv$trimmed[1])
      shiny::fluidRow(
        shiny::column(3, .infobox("list", "Sequences", length(rv$aln), "#2c7be5")),
        shiny::column(3, .infobox("arrows-alt-h", "Alignment columns", nc, "#6f42c1")),
        shiny::column(3, .infobox("check", "Columns kept", if (is.na(nk)) "-" else sprintf("%d (%.0f%%)", nk, 100 * nk / nc), "#27ae60")),
        shiny::column(3, .infobox("scissors", "Columns removed", if (is.na(nk)) "-" else nc - nk, "#e67e22")))
    })
    col_df <- shiny::reactive({
      shiny::req(rv$aln)
      s <- msa_column_stats(msa_matrix(rv$aln))
      s$Kept <- if (is.null(rv$kept)) NA else rv$kept
      s
    })
    output$col_plot <- shiny::renderPlot({
      s <- col_df()
      d <- rbind(data.frame(Col = s$Col, Value = 1 - s$GapFrac, Kept = s$Kept, Metric = "Fraction of sequences with a base"),
                 data.frame(Col = s$Col, Value = s$Cons, Kept = s$Kept, Metric = "Conservation"))
      d$Fill <- ifelse(is.na(d$Kept), "n/a", ifelse(d$Kept, "kept", "removed"))
      ggplot2::ggplot(d, ggplot2::aes(Col, Value, fill = Fill)) +
        ggplot2::geom_col(width = 1) +
        ggplot2::facet_wrap(~Metric, ncol = 1) +
        ggplot2::scale_fill_manual(values = c(kept = "#2c7be5", removed = "#c9ced6", "n/a" = "#7b8794"), name = NULL) +
        ggplot2::labs(x = "Alignment column", y = NULL) +
        ggplot2::coord_cartesian(ylim = c(0, 1)) +
        ggplot2::theme_minimal(base_size = 12)
    })

    # ---- site statistics (conserved / variable / parsimony-informative / singleton) --------------------------------
    ss_aln <- shiny::reactive({
      a <- if (identical(input$ss_src, "all") || is.null(rv$trimmed)) rv$aln else rv$trimmed
      shiny::req(a, length(a) >= 2)
      a
    })
    ss_df <- shiny::reactive({
      a <- ss_aln()
      M <- do.call(rbind, strsplit(unname(a), ""))
      rownames(M) <- names(a)
      .site_stats(M, isTRUE(input$ss_cd))
    })
    output$ss_note <- shiny::renderUI({
      a <- ss_aln()
      shiny::helpText(sprintf("%d sequences, %d columns (%s).", length(a), nchar(a[[1]]),
        if (identical(input$ss_src, "all") || is.null(rv$trimmed)) "whole alignment" else "trimmed result"))
    })
    output$ss_sum <- DT::renderDT(
      DT::datatable(.site_summary(ss_df()), rownames = FALSE, options = list(dom = "t", paging = FALSE)))
    output$ss_map <- shiny::renderPlot({
      d <- ss_df()
      d$Class <- factor(d$Class, levels = names(.site_colours))
      ggplot2::ggplot(d, ggplot2::aes(x = Site, y = 1, fill = Class)) +
        ggplot2::geom_tile(height = 1) +
        ggplot2::scale_fill_manual(values = .site_colours, drop = TRUE, name = NULL) +
        ggplot2::scale_x_continuous(expand = c(0, 0)) +
        ggplot2::labs(x = "Alignment column", y = NULL) +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::theme(axis.text.y = ggplot2::element_blank(), panel.grid = ggplot2::element_blank(),
                       legend.position = "top")
    })
    output$dl_ss_xlsx <- shiny::downloadHandler(
      filename = function() "alignment_site_statistics.xlsx",
      content = function(file) {
        d <- ss_df()
        sheets <- list(Summary = .site_summary(d), Base_composition = .site_composition(d), Sites = d,
                       Variable_sites = d[grepl("^Variable", d$Class), , drop = FALSE])
        tryCatch(.sg_xlsx(sheets, file), error = function(e) {
          shiny::showNotification(paste("Excel file not written:", conditionMessage(e)), type = "error", duration = 15)
          writeLines(conditionMessage(e), file)
        })
      })
    output$dl_ss_csv <- shiny::downloadHandler(
      filename = function() "alignment_sites.csv",
      content = function(file) utils::write.csv(ss_df(), file, row.names = FALSE))

    output$dl_aln <- shiny::downloadHandler(
      filename = function() "alignment.fasta",
      content = function(file) { shiny::req(rv$aln); writeLines(msa_fasta_text(rv$aln), file) })
    output$dl_trim <- shiny::downloadHandler(
      filename = function() "alignment_trimmed.fasta",
      content = function(file) { shiny::req(rv$trimmed); writeLines(msa_fasta_text(rv$trimmed), file) })
    output$dl_ungap <- shiny::downloadHandler(
      filename = function() "trimmed_ungapped.fasta",
      content = function(file) { shiny::req(rv$trimmed); writeLines(msa_fasta_text(stats::setNames(gsub("-", "", rv$trimmed), names(rv$trimmed))), file) })
    output$dl_cols <- shiny::downloadHandler(
      filename = function() "alignment_columns.csv",
      content = function(file) { utils::write.csv(col_df(), file, row.names = FALSE) })

    # ---- project file: what is saved / restored (used by the Project tab) -----------------------------
    get_state <- function() shiny::isolate(list(pool = rv$pool, aln = rv$aln, kept = rv$kept, trimmed = rv$trimmed,
                                                info = rv$info, meta = rv$meta, trim_meta = rv$trim_meta))
    set_state <- function(s) {
      rv$pool <- if (is.data.frame(s$pool)) s$pool else rv$pool[0, , drop = FALSE]
      rv$aln <- s$aln; rv$kept <- s$kept; rv$trimmed <- s$trimmed; rv$info <- s$info
      rv$meta <- s$meta; rv$trim_meta <- s$trim_meta; rv$msg <- NULL
      if (is.null(s$aln)) reset_alignment() else rv$ver <- rv$ver + 1L     # the editors reload from rv
      invisible()
    }
    list(get_state = get_state, set_state = set_state,
         pool = shiny::reactive(rv$pool),
         meta = shiny::reactive(list(align = rv$meta, trim = rv$trim_meta)),
         clear = function() { rv$pool <- rv$pool[0, , drop = FALSE]; reset_alignment() },
         aln = shiny::reactive(rv$aln), trimmed = shiny::reactive(rv$trimmed),
         add = function(recs, source = "BLAST") {
           # a new 'Send to MSA' from BLAST starts a fresh list unless the box is unticked
           if (identical(source, "BLAST") && isTRUE(shiny::isolate(input$blast_replace))) {
             rv$pool <- rv$pool[0, , drop = FALSE]; reset_alignment()
           }
           add_records(recs, source)
         })
  })
}
