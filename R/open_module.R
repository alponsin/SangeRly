# ==============================================================================
# SangeRly - "Open" tab: open and look at files, no assembly needed
#   .ab1 chromatograms | FASTA | FASTQ | GenBank | Newick trees | tables (CSV / TSV) | text
# Needs (sourced first or at least before the app runs):
#   assembly_core.R   read_ab1, prep_read, default_params, revcomp_chr, .nz
#   extras_core.R     translate_six, six_frame_html, GENETIC_CODES
#   assembly_module.R .card, .infobox, plot_read_chromatogram, plot_phred, .sv_local
# Packages: shiny, ggplot2, DT.  Optional: ape (tree files), Biostrings (translation).
# ==============================================================================

suppressPackageStartupMessages({ library(shiny); library(ggplot2) })

# ------------------------------------------------------------------------------
# Core functions (no Shiny)
# ------------------------------------------------------------------------------

# extensions that are picked up when a whole folder is opened
.op_known_ext <- c("ab1", "abi", "fa", "fasta", "fas", "fna", "ffn", "faa", "fsa", "seq", "fastq", "fq",
                   "gb", "gbk", "gbff", "genbank", "nwk", "newick", "tre", "tree", "treefile", "contree",
                   "csv", "tsv", "tab", "txt")

.op_label <- function(type) {
  m <- c(ab1 = "AB1 chromatogram", fasta = "FASTA", plain = "Plain sequence", fastq = "FASTQ", genbank = "GenBank",
         newick = "Newick tree", table = "Table", text = "Text", scf = "SCF (not supported)",
         binary = "Binary (not supported)", gz = "Compressed (not supported)", empty = "Empty file")
  out <- unname(m[type]); out[is.na(out)] <- type[is.na(out)]; out
}

.op_kind <- function(type) {
  switch(type, ab1 = "ab1", fasta = , plain = "fasta", fastq = "fastq", genbank = "genbank",
         newick = "tree", table = "table", text = "text", "msg")
}

.op_size <- function(x)
  ifelse(x >= 1e6, sprintf("%.1f MB", x / 1e6), ifelse(x >= 1e3, sprintf("%.0f KB", x / 1e3), sprintf("%.0f B", x)))

.op_num <- function(x, default) if (is.null(x) || length(x) != 1L || is.na(x)) default else x

# Format of a file, from its first bytes and its extension
.op_detect <- function(path, name = basename(path)) {
  ext <- tolower(tools::file_ext(name))
  b <- readBin(path, "raw", n = 4096L)
  if (!length(b)) return("empty")
  if (length(b) >= 4 && identical(b[1:4], charToRaw("ABIF"))) return("ab1")
  if (length(b) >= 4 && identical(b[1:4], charToRaw(".scf"))) return("scf")
  if (length(b) >= 2 && identical(b[1:2], as.raw(c(0x1f, 0x8b)))) return("gz")
  if (any(b == as.raw(0))) return("binary")
  txt <- iconv(rawToChar(b), "UTF-8", "UTF-8", sub = "?")
  lines <- strsplit(txt, "\r?\n")[[1]]
  lines <- lines[nzchar(trimws(lines))]
  if (!length(lines)) return("empty")
  first <- sub("^\ufeff", "", trimws(lines[1]))
  if (grepl("^>", first)) return("fasta")
  if (grepl("^@", first) && length(lines) >= 3 && grepl("^\\+", lines[3])) return("fastq")
  if (grepl("^LOCUS", first)) return("genbank")
  if (grepl("^\\(", first) && (ext %in% c("nwk", "newick", "tre", "tree", "treefile", "contree") || grepl(";\\s*$", first)))
    return("newick")
  if (ext %in% c("csv", "tsv", "tab")) return("table")
  if (ext %in% c("seq", "txt", "fa", "fasta", "fas", "fna", "ffn", "fsa", "faa")) {
    body <- toupper(gsub("[^A-Za-z]", "", paste(lines, collapse = "")))
    if (nchar(body) >= 20) {
      nuc <- nchar(gsub("[^ACGTUNRYKMSWBDHV]", "", body))
      if (nuc / nchar(body) >= 0.9) return("plain")
    }
  }
  "text"
}

.op_wrap <- function(s, w = 60L) {
  if (!nzchar(s)) return("")
  st <- seq(1L, nchar(s), by = w)
  substring(s, st, pmin(st + w - 1L, nchar(s)))
}

.op_fasta_text <- function(ids, seqs, desc = NULL, w = 60L) {
  unlist(lapply(seq_along(seqs), function(i)
    c(paste0(">", ids[i], if (!is.null(desc) && !is.na(desc[i]) && nzchar(desc[i])) paste0(" ", desc[i]) else ""),
      .op_wrap(seqs[i], w))))
}

# ---- sequence files -----------------------------------------------------------
.op_recs <- function(id, desc, seq, qual = NA_character_)
  data.frame(id = id, desc = desc, seq = seq, qual = qual, stringsAsFactors = FALSE)

.op_read_fasta <- function(path, name = basename(path)) {
  lines <- sub("\r$", "", readLines(path, warn = FALSE))
  lines <- lines[nzchar(trimws(lines))]
  hdr <- grepl("^>", lines)
  if (!any(hdr)) {
    s <- gsub("[^A-Za-z*.-]", "", paste(lines, collapse = ""))
    if (!nzchar(s)) stop("No sequence found in this file.")
    return(.op_recs(gsub("[^A-Za-z0-9._-]", "_", tools::file_path_sans_ext(name)), "", s))
  }
  grp <- cumsum(hdr)
  h <- sub("^>\\s*", "", lines[hdr])
  id <- sub("\\s.*$", "", h)
  desc <- trimws(sub("^\\S+\\s*", "", h))
  id[!nzchar(id)] <- paste0("seq", seq_along(id))[!nzchar(id)]
  body <- split(lines[!hdr & grp > 0], factor(grp[!hdr & grp > 0], levels = seq_along(h)))
  s <- vapply(body, function(x) gsub("[^A-Za-z*.-]", "", paste(x, collapse = "")), "", USE.NAMES = FALSE)
  .op_recs(id, desc, s)
}

.op_read_fastq <- function(path) {
  lines <- sub("\r$", "", readLines(path, warn = FALSE))
  while (length(lines) && !nzchar(lines[length(lines)])) lines <- lines[-length(lines)]
  n <- length(lines) %/% 4L
  if (n < 1L) stop("This is not a valid FASTQ file (fewer than 4 lines).")
  i <- seq(1L, by = 4L, length.out = n)
  if (!all(grepl("^@", lines[i])) || !all(grepl("^\\+", lines[i + 2L])))
    stop("This is not a valid FASTQ file (records must be 4 lines: @name, sequence, +, quality).")
  h <- sub("^@", "", lines[i])
  .op_recs(sub("\\s.*$", "", h), trimws(sub("^\\S+\\s*", "", h)), toupper(lines[i + 1L]), lines[i + 3L])
}

# One GenBank record (the lines between LOCUS and //)
.op_parse_gb <- function(ch) {
  block <- function(key) {
    i <- grep(paste0("^", key, "(\\s|$)"), ch)[1]
    if (is.na(i)) return("")
    j <- i + 1L
    while (j <= length(ch) && grepl("^\\s", ch[j])) j <- j + 1L
    rows <- ch[i:(j - 1L)]
    rows[1] <- sub(paste0("^", key), "", rows[1])
    gsub("\\s+", " ", trimws(paste(rows, collapse = " ")))
  }
  loc <- block("LOCUS")
  org <- grep("^\\s+ORGANISM\\s+", ch, value = TRUE)
  org <- if (length(org)) trimws(sub("^\\s+ORGANISM\\s+", "", org[1])) else ""
  qget <- function(q, key) {
    m <- grep(paste0("^", key, "="), q, value = TRUE)
    if (!length(m)) "" else gsub("^\"|\"$", "", sub(paste0("^", key, "="), "", m[1]))
  }
  feat <- data.frame(Type = character(), Location = character(), Gene = character(), Product = character(),
                     Qualifiers = character(), stringsAsFactors = FALSE)
  fi <- grep("^FEATURES", ch)
  oi <- grep("^(ORIGIN|CONTIG|BASE COUNT)", ch)
  if (length(fi)) {
    endf <- if (any(oi > fi[1])) min(oi[oi > fi[1]]) - 1L else length(ch)
    fl <- if (endf > fi[1]) ch[(fi[1] + 1L):endf] else character()
    types <- character(); locs <- character(); quals <- list(); state <- "loc"
    for (ln in fl) {
      if (grepl("^ {5}\\S", ln)) {
        types <- c(types, trimws(substr(ln, 6L, 20L))); locs <- c(locs, trimws(substr(ln, 21L, nchar(ln))))
        quals[[length(types)]] <- character(); state <- "loc"
      } else if (grepl("^ {21}/", ln) && length(types)) {
        quals[[length(types)]] <- c(quals[[length(types)]], sub("^ {21}/", "", ln)); state <- "q"
      } else if (length(types)) {
        cont <- trimws(ln)
        if (state == "loc") locs[length(types)] <- paste0(locs[length(types)], cont)
        else {
          k <- length(quals[[length(types)]])
          if (k) quals[[length(types)]][k] <- paste(quals[[length(types)]][k], cont)
        }
      }
    }
    if (length(types))
      feat <- data.frame(Type = types, Location = locs,
                         Gene = vapply(quals, qget, "", key = "gene"),
                         Product = vapply(quals, qget, "", key = "product"),
                         Qualifiers = substr(vapply(quals, paste, "", collapse = "; "), 1L, 400L),
                         stringsAsFactors = FALSE)
  }
  o <- grep("^ORIGIN", ch)
  seq <- if (length(o) && o[1] < length(ch)) toupper(gsub("[^A-Za-z]", "", paste(ch[(o[1] + 1L):length(ch)], collapse = ""))) else ""
  list(id = strsplit(loc, " ")[[1]][1], locus = loc, definition = block("DEFINITION"), accession = block("ACCESSION"),
       version = block("VERSION"), organism = org, features = feat, seq = seq)
}

.op_read_genbank <- function(path) {
  lines <- sub("\r$", "", readLines(path, warn = FALSE))
  ends <- which(lines == "//")
  starts <- c(1L, ends + 1L)
  starts <- starts[starts <= length(lines)]
  stops <- c(ends, length(lines))[seq_along(starts)]
  gb <- list()
  for (k in seq_along(starts)) {
    ch <- lines[starts[k]:stops[k]]
    ch <- ch[ch != "//"]
    if (any(grepl("^LOCUS", ch))) gb[[length(gb) + 1L]] <- .op_parse_gb(ch)
  }
  if (!length(gb)) stop("No GenBank record (LOCUS line) found in this file.")
  gb
}

# ---- statistics ---------------------------------------------------------------
.op_seq_kind <- function(s) {                    # "nucleotide" or "protein" (vectorised)
  L <- toupper(gsub("[^A-Za-z]", "", s))
  n <- nchar(L)
  nuc <- nchar(gsub("[^ACGTUNRYKMSWBDHV]", "", L))
  ifelse(n == 0, "nucleotide", ifelse(nuc / pmax(n, 1L) >= 0.9, "nucleotide", "protein"))
}

.op_recs_stats <- function(recs) {
  L <- gsub("[-.*]", "", toupper(recs$seq))
  typ <- .op_seq_kind(L)
  cnt <- function(ch) nchar(gsub(paste0("[^", ch, "]"), "", L))
  gc <- cnt("GC"); at <- cnt("ATU"); nn <- cnt("N")
  nuc <- typ == "nucleotide"
  d <- data.frame(Name = recs$id, Length = nchar(L), Type = typ,
                  GC_pct = ifelse(nuc & gc + at > 0, round(100 * gc / pmax(gc + at, 1L), 1), NA_real_),
                  N_count = ifelse(nuc, nn, NA_integer_), stringsAsFactors = FALSE)
  if (any(!is.na(recs$qual)))
    d$Mean_Phred <- vapply(recs$qual, function(q) if (is.na(q) || !nzchar(q)) NA_real_ else round(mean(utf8ToInt(q) - 33L), 1), 0)
  d
}

# Read one file of a given type into a list(ok = TRUE, ...) or list(ok = FALSE, msg = ...)
.op_load <- function(path, name, type) {
  pack <- function(recs, gb = NULL) list(ok = TRUE, recs = recs, stats = .op_recs_stats(recs), gb = gb,
                                         aligned = nrow(recs) > 1L && length(unique(nchar(recs$seq))) == 1L)
  tryCatch(
    switch(type,
      ab1 = { raw <- read_ab1(path, name); list(ok = TRUE, raw = raw, meta = .op_abif_meta(path), orig = raw$seq, del = rep(FALSE, length(raw$qual))) },
      fasta = , plain = pack(.op_read_fasta(path, name)),
      fastq = pack(.op_read_fastq(path)),
      genbank = { gb <- .op_read_genbank(path)
                  pack(.op_recs(vapply(gb, function(g) g$id, ""), vapply(gb, function(g) g$definition, ""),
                                vapply(gb, function(g) g$seq, "")), gb) },
      newick = .op_read_tree(path),
      table = .op_read_table(path),
      text = list(ok = TRUE, lines = readLines(path, n = 3000L, warn = FALSE)),
      scf = list(ok = FALSE, msg = "SCF chromatograms cannot be opened yet. Export the sample as .ab1 or FASTA from the sequencer software."),
      gz = list(ok = FALSE, msg = "This file is compressed (gzip). Unzip it first, then open the unzipped file."),
      binary = list(ok = FALSE, msg = "This is a binary file in a format that SangeRly cannot show."),
      empty = list(ok = FALSE, msg = "This file is empty."),
      list(ok = FALSE, msg = paste("Unknown format:", type))),
    error = function(e) list(ok = FALSE, msg = conditionMessage(e)))
}

.op_read_tree <- function(path) {
  if (!requireNamespace("ape", quietly = TRUE))
    stop("The 'ape' package is needed to open tree files: install.packages(\"ape\")")
  tr <- tryCatch(suppressWarnings(ape::read.tree(path)), error = function(e) NULL)
  if (is.null(tr)) tr <- tryCatch(suppressWarnings(ape::read.nexus(path)), error = function(e) NULL)
  if (is.null(tr)) stop("This file is not a Newick (or Nexus) tree.")
  n <- 1L
  if (inherits(tr, "multiPhylo")) { n <- length(tr); tr <- tr[[1]] }
  list(ok = TRUE, tree = tr, ntrees = n)
}

.op_read_table <- function(path) {
  l1 <- readLines(path, n = 1L, warn = FALSE)
  seps <- c("\t", ",", ";")
  cnt <- vapply(seps, function(s) lengths(regmatches(l1, gregexpr(s, l1, fixed = TRUE))), 1L)
  sep <- if (max(cnt) == 0L) "\t" else seps[which.max(cnt)]
  d <- utils::read.table(path, sep = sep, header = TRUE, quote = "\"", comment.char = "", check.names = FALSE,
                         stringsAsFactors = FALSE, fill = TRUE, nrows = 100000L)
  list(ok = TRUE, df = d, sep = sep)
}

# ---- AB1 run information (strings, dates and times that read_ab1 does not need) --
.op_abif_meta <- function(path) {
  raw <- readBin(path, "raw", n = file.size(path))
  i32 <- function(o) readBin(raw[o + 0:3], "integer", n = 1, size = 4, endian = "big")
  u16 <- function(o) readBin(raw[o + 0:1], "integer", n = 1, size = 2, signed = FALSE, endian = "big")
  nent <- i32(19); doff <- i32(27)
  o <- doff + (seq_len(nent) - 1L) * 28L + 1L
  tag <- vapply(o, function(x) rawToChar(raw[x:(x + 3L)]), "")
  num <- vapply(o, function(x) i32(x + 4L), 1L)
  typ <- vapply(o, function(x) u16(x + 8L), 1L)
  cnt <- vapply(o, function(x) i32(x + 12L), 1L)
  size <- vapply(o, function(x) i32(x + 16L), 1L)
  clean <- function(x) iconv(x, "latin1", "UTF-8", sub = "?")
  val <- function(t, n = 1L) {
    k <- which(tag == t & num == n)[1]
    if (is.na(k) || size[k] <= 0) return(NA_character_)
    st <- if (size[k] <= 4) o[k] + 20L else i32(o[k] + 20L) + 1L
    b <- raw[st:(st + size[k] - 1L)]
    switch(as.character(typ[k]),
      "18" = { l <- as.integer(b[1]); b <- b[-1]; b <- b[seq_len(min(l, length(b)))]; clean(rawToChar(b[b != as.raw(0)])) },
      "19" = , "2" = clean(rawToChar(b[b != as.raw(0)])),
      "10" = sprintf("%04d-%02d-%02d", readBin(b[1:2], "integer", size = 2, endian = "big"), as.integer(b[3]), as.integer(b[4])),
      "11" = sprintf("%02d:%02d:%02d", as.integer(b[1]), as.integer(b[2]), as.integer(b[3])),
      "4" = as.character(readBin(b, "integer", n = cnt[k], size = 2, endian = "big")[1]),
      "5" = as.character(readBin(b, "integer", n = 1, size = 4, endian = "big")),
      NA_character_)
  }
  dt <- function(d, t) { a <- val(d[1], d[2]); b <- val(t[1], t[2]); if (is.na(a)) NA_character_ else paste(a, if (is.na(b)) "" else b) }
  sw <- stats::na.omit(c(val("SVER", 1L), val("SVER", 2L)))
  items <- list("Sample name" = val("SMPL"), "Instrument model" = val("MODL"), "Instrument name" = val("MCHN"),
                "Dye set" = val("DySN"), "Well" = val("TUBE"), "Lane / capillary" = val("LANE"),
                "Plate / container" = val("CTNM"), "Run start" = dt(c("RUND", 1L), c("RUNT", 1L)),
                "Run end" = dt(c("RUND", 2L), c("RUNT", 2L)), "Run module" = val("RMdN"),
                "Basecaller / mobility file" = val("PDMF", 2L), "Polymer" = trimws(val("GTyp")),
                "Software" = if (length(sw)) paste(sw, collapse = ", ") else NA_character_,
                "User" = val("User"), "Comment" = val("CMNT"))
  v <- unlist(items)
  v <- v[!is.na(v) & nzchar(trimws(v))]
  data.frame(Item = names(v), Value = unname(v), stringsAsFactors = FALSE)
}

# Write a copy of an .ab1 file with edited base calls (same number of bases; the traces are untouched).
# ABIF keeps the called bases in PBAS tag number 2 (number 1 is the instrument's original); like other
# editors, the edited calls go to number 2.
.op_ab1_patch <- function(path, seq, out, keep = NULL) {
  raw <- readBin(path, "raw", n = file.size(path))
  i32 <- function(o) readBin(raw[o + 0:3], "integer", n = 1, size = 4, endian = "big")
  put32 <- function(o, v) raw[o + 0:3] <<- writeBin(as.integer(v), raw(), size = 4, endian = "big")
  nent <- i32(19); doff <- i32(27)
  o <- doff + (seq_len(nent) - 1L) * 28L + 1L
  tag <- vapply(o, function(x) rawToChar(raw[x:(x + 3L)]), "")
  num <- vapply(o, function(x) i32(x + 4L), 1L)
  cnt <- vapply(o, function(x) i32(x + 12L), 1L)
  size <- vapply(o, function(x) i32(x + 16L), 1L)
  pick <- function(t) { k <- which(tag == t & num == 2L)[1]; if (is.na(k)) k <- which(tag == t & num == 1L)[1]; k }
  start <- function(k) if (size[k] <= 4) o[k] + 20L else i32(o[k] + 20L) + 1L
  kb <- pick("PBAS")
  if (is.na(kb)) stop("This file has no base calls to edit.")
  n <- cnt[kb]
  if (nchar(seq) != n) stop("The edited sequence and the file disagree on the number of bases.")
  if (is.null(keep)) keep <- rep(TRUE, n)
  if (length(keep) != n || !any(keep)) stop("Nothing would be left of this read.")
  chars <- charToRaw(seq)
  raw[start(kb) + seq_len(n) - 1L] <- chars          # edited base calls (same positions)
  if (!all(keep)) {                                   # deleted bases: shrink PBAS, PCON and PLOC together
    m <- sum(keep)
    for (t in c("PBAS", "PCON", "PLOC")) {
      k <- pick(t)
      if (is.na(k) || cnt[k] != n) stop("This file has no consistent ", t, " data to shorten.")
      w <- size[k] %/% n                              # bytes per base (1 for PBAS and PCON, 2 for PLOC)
      st <- start(k)
      old <- matrix(raw[st + seq_len(n * w) - 1L], nrow = w)
      new <- as.vector(old[, keep, drop = FALSE])
      raw[st + seq_len(m * w) - 1L] <- new
      put32(o[k] + 12L, m); put32(o[k] + 16L, m * w)
    }
  }
  writeBin(raw, out)
  invisible(out)
}

# A read as plot_read_chromatogram() / plot_phred() expect it, with trimming only for display
.op_ab1_rd <- function(raw, q_trim = 20, win = 10) {
  n <- length(raw$qual)
  P <- default_params(); P$q_trim <- q_trim; P$win <- win
  qc <- tryCatch(prep_read(raw, P)$qc, error = function(e) NULL)
  from <- if (is.null(qc)) 1L else as.integer(qc$from)
  to <- if (is.null(qc)) n else as.integer(qc$to)
  list(name = raw$name, seq = raw$seq, qual = raw$qual, ploc = raw$ploc, traces = raw$traces, rc = FALSE,
       raw = list(seq = strsplit(raw$seq, "")[[1]], qual = raw$qual, ploc = raw$ploc),
       qc = list(qual_full = raw$qual, from = from, to = to, n_full = n), trimmed = !is.null(qc))
}

# ---- HTML views -----------------------------------------------------------------
.op_col <- function(x) gsub("([ACGTUN])", "<i class='b\\1'>\\1</i>", x)

# One sequence in lines of `width`, in blocks of 10, with the position of the first base of each line
.op_seq_html <- function(seq, width = 60L, colour = TRUE, protein = FALSE, max_nt = 60000L) {
  seq <- toupper(seq); n <- nchar(seq)
  if (n == 0) return(HTML("<p><em>(empty sequence)</em></p>"))
  shown <- min(n, max_nt); s <- substr(seq, 1L, shown)
  width <- max(10L, as.integer(width))
  nd <- nchar(as.character(shown))
  starts <- seq(1L, shown, by = width)
  lines <- vapply(starts, function(a) {
    chunk <- substr(s, a, min(shown, a + width - 1L))
    st <- seq(1L, nchar(chunk), by = 10L)
    blocks <- substring(chunk, st, pmin(st + 9L, nchar(chunk)))
    if (colour && !protein) blocks <- .op_col(blocks)
    sprintf("%s  %s", formatC(a, width = nd), paste(blocks, collapse = " "))
  }, "")
  HTML(paste0("<pre class='op-seq'>", paste(lines, collapse = "\n"), "</pre>",
              if (n > shown) sprintf("<p><em>Showing the first %s of %s bases; download the FASTA file for all of it.</em></p>",
                                     format(shown, big.mark = ","), format(n, big.mark = ",")) else ""))
}

# Aligned sequences: names on the left, blocks of `width` columns, * under the identical columns
.op_aln_html <- function(seqs, width = 80L, dots = FALSE, max_seq = 200L, max_col = 8000L) {
  n_all <- length(seqs); k <- min(n_all, max_seq)
  nm <- names(seqs)[seq_len(k)]; s <- toupper(seqs[seq_len(k)])
  L_all <- max(nchar(s)); L <- min(L_all, max_col)
  s <- substr(s, 1L, L); s <- paste0(s, strrep("-", L - nchar(s)))
  M <- do.call(rbind, strsplit(s, ""))
  same <- if (k > 1) apply(M, 2, function(x) x[1] != "-" && all(x == x[1])) else rep(FALSE, L)
  if (dots && k > 1) {
    sub <- M[-1, , drop = FALSE]
    eq <- sub == matrix(M[1, ], k - 1L, L, byrow = TRUE) & sub != "-"
    sub[eq] <- "."; M[-1, ] <- sub
  }
  nw <- min(24L, max(nchar(nm)))
  lab <- formatC(substr(nm, 1L, nw), width = -nw)
  width <- max(20L, as.integer(width))
  rows <- character()
  for (a in seq(1L, L, by = width)) {
    b <- min(L, a + width - 1L)
    rows <- c(rows, paste0(strrep(" ", nw + 1L), a))
    for (r in seq_len(k))
      rows <- c(rows, paste0(htmltools::htmlEscape(lab[r]), " ", .op_col(paste(M[r, a:b], collapse = ""))))
    if (k > 1) rows <- c(rows, paste0(strrep(" ", nw + 1L), paste(ifelse(same[a:b], "*", " "), collapse = "")))
    rows <- c(rows, "")
  }
  note <- if (n_all > k || L_all > L)
    sprintf("<p><em>Showing the first %d of %d sequences and the first %d of %d columns.</em></p>", k, n_all, L, L_all) else ""
  HTML(paste0("<pre class='op-seq'>", paste(rows, collapse = "\n"), "</pre>", note))
}

# ------------------------------------------------------------------------------
# UI
# ------------------------------------------------------------------------------
# Pressing Save / Save as / Send to MSA first applies any pending edit of the AliView-like editor.
.op_shim <- tags$script(HTML("
(function(){
  if (window.__opShim) return; window.__opShim = true;
  document.addEventListener('click', function(e){
    var b = e.target.closest ? e.target.closest('button') : null;
    if (!b || !b.id || !/(save|save_as|to_msa)$/.test(b.id)) return;
    var rid = b.id.replace(/(save_as|save|to_msa)$/, 'ed');
    var reg = window.SangerlyMSA && window.SangerlyMSA._registry, ed = reg && reg[rid];
    if (ed && ed.model && ed.model.rev !== ed.model.savedRev && ed.root && ed.root.offsetParent) ed.commit();
  }, true);
})();"))

# Everything that can be shown is in the page from the start; the server only says which kind of file is open
# (output$kind_out) and the panels are shown or hidden. That keeps the viewer and the editor alive between files.
.op_main_ui <- function(ns) {
  cp <- function(kinds, ...)
    conditionalPanel(condition = sprintf("[%s].indexOf(output.kind_out) >= 0", paste0("'", kinds, "'", collapse = ",")), ns = ns, ...)
  formats_table <- tags$table(class = "table table-sm",
    tags$thead(tags$tr(tags$th("Format"), tags$th("Extensions"), tags$th("What you see"))),
    tags$tbody(
      tags$tr(tags$td("AB1 chromatogram"), tags$td(".ab1 .abi"), tags$td("The chromatogram viewer of the Assembly tab; delete or change bases, then save")),
      tags$tr(tags$td("FASTA / plain sequence"), tags$td(".fa .fasta .fas .fna .fsa .seq .txt"), tags$td("The alignment editor of the MSA tab: edit directly, rename, delete, then save")),
      tags$tr(tags$td("FASTQ"), tags$td(".fastq .fq"), tags$td("The same editor, plus Phred quality")),
      tags$tr(tags$td("GenBank"), tags$td(".gb .gbk .gbff"), tags$td("The same editor, plus record information and features")),
      tags$tr(tags$td("Newick tree"), tags$td(".nwk .tre .tree .treefile"), tags$td("The tree drawn, with support values")),
      tags$tr(tags$td("Table"), tags$td(".csv .tsv"), tags$td("The table, searchable")),
      tags$tr(tags$td("Other text"), tags$td("any"), tags$td("The first 3000 lines"))))
  bar <- div(class = "op-bar",
    actionButton(ns("save"), "Save", icon = icon("floppy-disk")),
    actionButton(ns("save_as"), "Save as...", icon = icon("file-export")),
    actionButton(ns("undo"), "Undo", icon = icon("rotate-left")),
    actionButton(ns("to_msa"), "Send to the MSA tab", icon = icon("arrow-right")),
    uiOutput(ns("dirty_note"), inline = TRUE))
  tagList(
    cp("none", .card("Open a file",
                     p("Use ", tags$b("Open files..."), " on the left (you can choose many files at once, from any folder). ",
                       "AB1 files open in the chromatogram viewer and FASTA / FASTQ / GenBank files in the alignment editor, ",
                       "both editable. Nothing is written to disk until you press Save or Save as."),
                     formats_table)),
    cp("msg", uiOutput(ns("msg"))),
    cp(c("ab1", "fasta", "fastq", "genbank", "tree", "table", "text"), bar),
    cp("ab1", tagList(
      tags$details(class = "op-trim", tags$summary("AB1 trimming (display only)"),
        div(class = "d-flex flex-wrap gap-3 align-items-end",
          numericInput(ns("q_trim"), "Window mean Phred >=", 20, 0, 60, width = "190px"),
          numericInput(ns("win"), "Window size", 10, 2, 40, width = "130px")),
        helpText("Only used to mark and export the trimmed part of an .ab1 read here. Nothing is changed in the Assembly tab.")),
      uiOutput(ns("ab1_boxes")),
      tabsetPanel(id = ns("ab1_tabs"),
        tabPanel("Chromatogram", br(), svViewerUI(ns("viewer")),
          helpText("Click the viewer first, then: click a base to select it (Shift+click or drag along the ruler for a range); ",
                   "Delete or Backspace removes the selected bases; type A C G T N (or an IUPAC letter) to change a base, 0 puts the original back; ",
                   "Ctrl+Z undoes; Ctrl+wheel zooms. Changed bases are green, deleted bases are struck out. ",
                   "Save writes the changes into the .ab1 file (the traces are kept), or Save as makes a copy.")),
        tabPanel("Quality", br(), plotOutput(ns("phred_plot"), height = "340px"),
                 helpText("Dotted lines are Q20 and Q30. Orange = the part kept by the trimming set on the left. Deleted bases are left out.")),
        tabPanel("Sequence", br(),
          fluidRow(column(4, radioButtons(ns("ab1_region"), NULL, c("Trimmed region" = "trim", "Whole read" = "all"))),
                   column(3, checkboxInput(ns("ab1_rc"), "Reverse complement", FALSE)),
                   column(3, checkboxInput(ns("ab1_col"), "Colour the bases", TRUE)),
                   column(2, numericInput(ns("ab1_w"), "Per line", 60, 10, 200, 10))),
          div(class = "op-bar", downloadButton(ns("dl_ab1_fa"), "FASTA"), downloadButton(ns("dl_ab1_fq"), "FASTQ"),
              downloadButton(ns("dl_ab1_csv"), "Bases + Phred (CSV)")),
          uiOutput(ns("ab1_seq_view"))),
        tabPanel("File info", br(), DT::DTOutput(ns("ab1_info")))))),
    cp(c("fasta", "fastq", "genbank"),
      tabsetPanel(id = ns("seq_tabs"),
        tabPanel("Editor", value = "ed", br(),
          sg_findbar(ns, paste("Finds the sequence and its reverse complement in every sequence of the file (strand + = as typed,",
                               "strand - = reverse complement). The sequences as last applied are searched (press Apply edits first to include your edits).",
                               "The hits are coloured in the editor; click a row of the result, or use the arrows, to go to a hit.",
                               "Enter = find (again = next hit), Shift+Enter = previous hit, Esc = close.")),
          msaEditorUI(ns("ed"), empty_text = "Open a FASTA, FASTQ or GenBank file and it appears here.",
                      extra = sg_hint(paste("Click a base or drag to select; type a letter to change it; Delete / Space / Backspace as described in the editing tools;",
                                            "double-click a name to rename it; right-click a name for more. Press Apply edits (Save and Save as do this for you),",
                                            "then Save or Save as to write the file. Ctrl+F opens the find bar.")))),
        tabPanel("Sequences", value = "recs", br(),
          div(class = "op-bar",
              actionButton(ns("del_recs"), "Delete selected", icon = icon("trash")),
              actionButton(ns("rc_recs"), "Reverse complement selected"),
              actionButton(ns("edit_seq"), "Edit sequence...", icon = icon("pen")),
              actionButton(ns("crop_seq"), "Crop..."),
              actionButton(ns("add_seq"), "Add sequences...", icon = icon("plus"))),
          DT::DTOutput(ns("recs_tbl")),
          helpText("Double-click a Name or Description to rename it. Click a row to select it (Ctrl/Shift-click for several).")),
        tabPanel("Record info", value = "gbinfo", br(), DT::DTOutput(ns("gb_info"))),
        tabPanel("Features", value = "gbfeat", br(), DT::DTOutput(ns("gb_feat"))),
        tabPanel("Quality", value = "fq", br(), plotOutput(ns("fq_plot"), height = "320px"),
                 helpText("Phred quality of every base of the selected read (first 5000 bases).")),
        tabPanel("Translation", value = "tr", br(),
          fluidRow(column(6, selectInput(ns("tr_code"), "Genetic code",
                                         if (exists("GENETIC_CODES")) GENETIC_CODES else c("Standard (1)" = "1"), "1")),
                   column(3, numericInput(ns("tr_min"), "Shortest ORF (aa)", 30, 10, 500, 10))),
          uiOutput(ns("tr_view")), DT::DTOutput(ns("tr_orfs"))))),
    cp("tree", tagList(
      fluidRow(column(3, selectInput(ns("tree_layout"), "Layout", c("Phylogram" = "phylogram", "Cladogram" = "cladogram",
                                                                   "Unrooted" = "unrooted", "Fan" = "fan", "Radial" = "radial"))),
               column(3, sliderInput(ns("tree_cex"), "Label size", 0.4, 2, 0.9, 0.1)),
               column(3, checkboxInput(ns("tree_sup"), "Show node labels (support values)", TRUE)),
               column(3, br(), downloadButton(ns("dl_tree"), "Newick"))),
      uiOutput(ns("tree_info")), plotOutput(ns("tree_plot"), height = "auto"))),
    cp("table", tagList(uiOutput(ns("table_info")), DT::DTOutput(ns("table_view")))),
    cp("text", tagList(uiOutput(ns("text_info")), verbatimTextOutput(ns("text_view")))))
}

openUI <- function(id) {
  ns <- NS(id)
  css <- tags$style(HTML("
    pre.op-seq{font-family:Consolas,Menlo,monospace;font-size:12.5px;line-height:1.5;background:#fff;border:1px solid #e3e8ee;
               border-radius:6px;padding:10px;overflow-x:auto;white-space:pre}
    pre.op-seq i{font-style:normal}
    .bA{background:#c6f0d6}.bC{background:#cfe4fb}.bG{background:#e3e3e3}.bT,.bU{background:#ffd0d0}.bN{background:#f4f0c8}
    .op-hdr{font-family:Consolas,Menlo,monospace;font-size:12.5px;color:#1f3c5f;margin:6px 0;word-break:break-all}
    .op-bar{display:flex;gap:8px;flex-wrap:wrap;margin-bottom:10px;align-items:center}
    .op-ed textarea{font-family:Consolas,Menlo,monospace;font-size:12.5px;white-space:pre}
  "))
  div(class = "op-root", sidebarLayout(
    sidebarPanel(
      width = 3,
      div(class = "sidebar-sec",
        tags$details(open = NA, tags$summary("Open files"),
          fileInput(ns("files"), "Choose one or more files (any format below)", multiple = TRUE,
                    buttonLabel = "Open files...", placeholder = "no file chosen"),
          if (.sv_local()) tagList(
            textInput(ns("path"), "Or type the path of a file or a folder", "", placeholder = "C:/data/read.ab1"),
            div(style = "display:flex;gap:8px;flex-wrap:wrap;margin-bottom:8px",
                actionButton(ns("path_add"), "Open path"),
                if (.Platform$OS.type == "windows") tagList(actionButton(ns("browse_files"), "Browse files..."),
                                                            actionButton(ns("browse_dir"), "Browse folder..."))),
            checkboxInput(ns("path_rec"), "Include subfolders when a folder is opened", TRUE)),
          div(style = "display:flex;gap:6px;flex-wrap:wrap",
              actionButton(ns("rm"), "Close this file", icon = icon("xmark")),
              actionButton(ns("clear"), "Close all"))))),
    mainPanel(width = 9, css, .op_shim, uiOutput(ns("files_card")), .op_main_ui(ns))))
}

# ------------------------------------------------------------------------------
# Server
# ------------------------------------------------------------------------------
# msa_add: function(recs, source) that adds sequences to the MSA tab (msa$add); go: function(tab) that opens a tab.
openServer <- function(id, msa_add = NULL, go = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    rv <- reactiveValues(files = data.frame(key = character(), name = character(), path = character(), type = character(),
                                            size = numeric(), where = character(), dirty = logical(), stringsAsFactors = FALSE),
                         sel = NA_integer_, ver = 0L, sver = 0L, rel = 0L)
    kind <- reactiveVal("none")
    cache <- new.env(parent = emptyenv())
    hist <- new.env(parent = emptyenv())     # undo history, one list per file
    up_dir <- file.path(tempdir(), paste0("sangerly_open_", paste(sample(c(letters, 0:9), 8, TRUE), collapse = "")))
    dir.create(up_dir, showWarnings = FALSE, recursive = TRUE)
    up_n <- 0L

    # ---- adding and removing files ---------------------------------------------------------------------
    add_items <- function(name, path, size, key, where) {
      keep <- !(key %in% rv$files$key) & !duplicated(key)
      if (!any(keep)) {
        showNotification("That file is already open.", type = "message", duration = 4)
        i <- match(key[1], rv$files$key); if (!is.na(i)) rv$sel <- i
        return(invisible())
      }
      type <- mapply(function(p, n) tryCatch(.op_detect(p, n), error = function(e) "text"),
                     path[keep], name[keep], USE.NAMES = FALSE)
      nd <- data.frame(key = key[keep], name = name[keep], path = path[keep], type = type, size = size[keep],
                       where = where[keep], dirty = FALSE, stringsAsFactors = FALSE)
      first <- nrow(rv$files) + 1L
      rv$files <- rbind(rv$files, nd)
      rv$sel <- first
      showNotification(sprintf("%d file(s) opened%s.", nrow(nd),
                               if (sum(!keep)) sprintf(", %d already open (skipped)", sum(!keep)) else ""),
                       type = "message", duration = 4)
    }
    observeEvent(input$files, {
      u <- input$files; req(nrow(u) > 0)
      dest <- character(nrow(u))
      for (i in seq_len(nrow(u))) {
        up_n <<- up_n + 1L
        dest[i] <- file.path(up_dir, sprintf("%05d_%s", up_n, u$name[i]))
        file.copy(u$datapath[i], dest[i], overwrite = TRUE)
      }
      add_items(u$name, dest, u$size, paste(u$name, u$size, sep = "|"), rep("(uploaded)", nrow(u)))
    })
    add_paths <- function(p) {
      p <- normalizePath(p, winslash = "/", mustWork = FALSE)
      add_items(basename(p), p, file.size(p), p, dirname(p))
    }
    add_dir <- function(d) {
      f <- list.files(d, recursive = isTRUE(input$path_rec), full.names = TRUE)
      f <- f[!dir.exists(f) & tolower(tools::file_ext(f)) %in% .op_known_ext]
      if (!length(f)) { showNotification("No file with a known extension in that folder.", type = "warning"); return() }
      if (length(f) > 500L) {
        showNotification(sprintf("%d files found; only the first 500 are opened.", length(f)), type = "warning")
        f <- f[1:500]
      }
      add_paths(f)
    }
    observeEvent(input$path_add, {
      p <- trimws(gsub("^\"|\"$", "", .nz(input$path, "")))
      if (!nzchar(p)) { showNotification("Type the path of a file or a folder first.", type = "warning"); return() }
      if (dir.exists(p)) add_dir(p)
      else if (file.exists(p)) add_paths(p)
      else showNotification("That path does not exist.", type = "error")
    })
    # Qt window, Navigation area: double-click on a file opens it here (sgQt.openPaths)
    observeEvent(input$qt_open, {
      p <- unlist(input$qt_open$paths); p <- p[file.exists(p) & !dir.exists(p)]
      if (!length(p)) { showNotification("That file was not found.", type = "warning"); return() }
      p <- normalizePath(p, winslash = "/", mustWork = FALSE)
      new <- p[!(p %in% rv$files$key)]
      if (length(new)) add_paths(new)                 # files that are open already are only selected, without a message
      i <- match(p[1], rv$files$key); if (!is.na(i)) rv$sel <- i
    })
    observeEvent(input$browse_files, {
      p <- tryCatch(utils::choose.files(multi = TRUE), error = function(e) character())
      if (length(p)) add_paths(p)
    })
    observeEvent(input$browse_dir, {
      d <- tryCatch(utils::choose.dir(caption = "Choose the folder with the files"), error = function(e) NA_character_)
      if (length(d) == 1L && !is.na(d) && nzchar(d)) add_dir(d)
    })
    drop_current <- function() {
      i <- rv$sel; if (is.na(i) || i > nrow(rv$files)) return()
      key <- rv$files$key[i]
      for (e in list(cache, hist)) if (exists(key, envir = e, inherits = FALSE)) rm(list = key, envir = e)
      rv$files <- rv$files[-i, , drop = FALSE]
      rv$sel <- if (nrow(rv$files)) min(i, nrow(rv$files)) else NA_integer_
    }
    drop_all <- function() {
      rv$files <- rv$files[0, , drop = FALSE]; rv$sel <- NA_integer_
      rm(list = ls(cache, all.names = TRUE), envir = cache); rm(list = ls(hist, all.names = TRUE), envir = hist)
    }
    confirm <- function(title, msg, id, label)
      showModal(modalDialog(title = title, msg, footer = tagList(modalButton("Cancel"), actionButton(ns(id), label)), easyClose = TRUE))
    observeEvent(input$clear, {
      if (any(rv$files$dirty)) confirm("Unsaved changes", "Some files have changes that are not saved. Close all anyway?", "clear_force", "Close without saving")
      else drop_all()
    })
    observeEvent(input$clear_force, { removeModal(); drop_all() })
    observeEvent(input$rm, {
      i <- rv$sel; if (is.na(i) || i > nrow(rv$files)) return()
      if (isTRUE(rv$files$dirty[i])) confirm("Unsaved changes", "This file has changes that are not saved. Close it anyway?", "rm_force", "Close without saving")
      else drop_current()
    })
    observeEvent(input$rm_force, { removeModal(); drop_current() })

    # ---- the list of open files ------------------------------------------------------------------------
    output$files_card <- renderUI({
      req(nrow(rv$files) > 0)
      .card("Opened files",
            DT::DTOutput(ns("files_tbl")),
            div(style = "display:flex;gap:6px;flex-wrap:wrap;margin-top:6px",
                actionButton(ns("prev"), "< Previous"), actionButton(ns("nxt"), "Next >"),
                actionButton(ns("rename"), "Rename file...", icon = icon("pen")),
                actionButton(ns("del_disk"), "Delete file from disk...", icon = icon("trash"))),
            helpText("Click a row to open that file. * = changed and not saved yet."))
    })
    output$files_tbl <- DT::renderDT({
      f <- rv$files; req(nrow(f) > 0)
      d <- data.frame(File = paste0(f$name, ifelse(f$dirty, "  *", "")), Format = .op_label(f$type), Size = .op_size(f$size), stringsAsFactors = FALSE)
      DT::datatable(d, rownames = FALSE, selection = list(mode = "single", selected = isolate(rv$sel), target = "row"),
                    options = list(pageLength = 6, scrollX = TRUE, dom = if (nrow(d) > 6) "ftp" else "t"))
    })
    observeEvent(input$files_tbl_rows_selected, {
      s <- input$files_tbl_rows_selected
      if (length(s) && s[1] <= nrow(rv$files)) rv$sel <- as.integer(s[1])
    })
    move <- function(d) {
      n <- nrow(rv$files); if (!n || is.na(rv$sel)) return()
      k <- min(n, max(1L, rv$sel + d)); rv$sel <- k
      p <- DT::dataTableProxy("files_tbl")
      DT::selectRows(p, k); DT::selectPage(p, (k - 1L) %/% 6L + 1L)
    }
    observeEvent(input$prev, move(-1L))
    observeEvent(input$nxt, move(1L))

    # ---- the current file ----------------------------------------------------------------------------------
    cur <- reactive({
      rv$ver
      i <- rv$sel; f <- rv$files
      req(!is.na(i), i <= nrow(f))
      row <- f[i, , drop = FALSE]
      if (!exists(row$key, envir = cache, inherits = FALSE))
        assign(row$key, .op_load(row$path, row$name, row$type), envir = cache)
      list(row = row, type = row$type, data = get(row$key, envir = cache, inherits = FALSE))
    })
    observe({
      k <- tryCatch({ x <- cur(); if (!isTRUE(x$data$ok)) "msg" else .op_kind(x$type) }, error = function(e) "none")
      if (!identical(k, isolate(kind()))) kind(k)
    })
    # ===== editing =========================================================================================
    # mutate(): change the data of the current file with fn(old) -> new; the old data goes on the undo stack.
    mutate <- function(fn, structural = TRUE) {
      i <- rv$sel
      if (is.na(i) || i > nrow(rv$files) || !exists(rv$files$key[i], envir = cache, inherits = FALSE)) return(invisible(FALSE))
      key <- rv$files$key[i]
      old <- get(key, envir = cache, inherits = FALSE)
      new <- tryCatch(fn(old), error = function(e) { showNotification(conditionMessage(e), type = "error", duration = 8); NULL })
      if (is.null(new)) return(invisible(FALSE))
      h <- if (exists(key, envir = hist, inherits = FALSE)) get(key, envir = hist) else list()
      assign(key, c(list(old), h)[seq_len(min(31L, length(h) + 1L))], envir = hist)
      assign(key, new, envir = cache)
      rv$files$dirty[i] <- TRUE
      rv$ver <- rv$ver + 1L
      if (structural) rv$sver <- rv$sver + 1L
      invisible(TRUE)
    }
    rebuild <- function(d, recs, gb = d$gb) {
      rownames(recs) <- NULL
      list(ok = TRUE, recs = recs, stats = .op_recs_stats(recs), gb = gb, aligned = d$aligned)
    }
    do_undo <- function() {
      i <- rv$sel; req(!is.na(i), i <= nrow(rv$files)); key <- rv$files$key[i]
      h <- if (exists(key, envir = hist, inherits = FALSE)) get(key, envir = hist) else list()
      if (!length(h)) { showNotification("Nothing to undo.", type = "message"); return() }
      assign(key, h[[1]], envir = cache); assign(key, h[-1], envir = hist)
      if (length(h) == 1L) rv$files$dirty[i] <- FALSE
      rv$ver <- rv$ver + 1L; rv$sver <- rv$sver + 1L
    }
    observeEvent(input$undo, do_undo())
    observeEvent(input$viewer_undo, do_undo())
    output$dirty_note <- renderUI({
      i <- rv$sel
      if (!is.na(i) && i <= nrow(rv$files) && isTRUE(rv$files$dirty[i]))
        span(style = "color:#c0392b;font-weight:600", "Unsaved changes")
    })
    sel_rows <- function() {
      s <- input$recs_tbl_rows_selected
      if (!length(s)) { showNotification("Select one or more rows in the table first.", type = "warning"); return(NULL) }
      as.integer(s)
    }

    # --- sequence files: rename, delete, reverse complement, edit, crop, add
    observeEvent(input$recs_tbl_cell_edit, {
      e <- input$recs_tbl_cell_edit
      field <- switch(as.character(as.integer(e$col) + 1L), "2" = "id", "3" = "desc", NULL)
      if (is.null(field)) return()
      v <- trimws(as.character(e$value))
      if (field == "id") v <- gsub("\\s+", "_", v)
      if (field == "id" && !nzchar(v)) { showNotification("A name cannot be empty.", type = "error"); rv$sver <- rv$sver + 1L; return() }
      mutate(function(d) {
        d$recs[e$row, field] <- v
        if (field == "id" && !is.null(d$gb)) d$gb[[e$row]]$id <- v
        rebuild(d, d$recs)
      })
    })
    observeEvent(input$del_recs, {
      sel <- sel_rows(); req(sel)
      mutate(function(d) {
        n <- nrow(d$recs)
        if (length(sel) >= n) stop("That would remove every sequence. Use \"Close this file\" instead.")
        keep <- setdiff(seq_len(n), sel)
        rebuild(d, d$recs[keep, , drop = FALSE], if (!is.null(d$gb)) d$gb[keep])
      })
    })
    observeEvent(input$rc_recs, {
      sel <- sel_rows(); req(sel)
      mutate(function(d) {
        r <- d$recs; gb <- d$gb
        for (j in sel) if (d$stats$Type[j] == "nucleotide") {
          r$seq[j] <- revcomp_chr(toupper(r$seq[j]))
          if (!is.na(r$qual[j])) r$qual[j] <- paste(rev(strsplit(r$qual[j], "")[[1]]), collapse = "")
          if (!is.null(gb)) gb[[j]]$seq <- r$seq[j]
        }
        rebuild(d, r, gb)
      })
    })
    open_edit_seq <- function() {
      d <- tryCatch(recs(), error = function(e) NULL); req(d); i <- rec_i()
      s <- d$recs$seq[i]
      if (nchar(s) > 200000L) { showNotification("This sequence is too long to edit here (over 200,000 bases).", type = "warning"); return() }
      showModal(modalDialog(title = paste("Edit sequence:", d$recs$id[i]), size = "l",
        div(class = "op-ed", textAreaInput(ns("ed_text"), NULL, value = paste(.op_wrap(s, 60L), collapse = "\n"), rows = 12, width = "100%")),
        helpText("Type or paste bases; spaces, line breaks and numbers are ignored. For a FASTQ read the number of bases must stay the same (use Crop to shorten)."),
        footer = tagList(modalButton("Cancel"), actionButton(ns("ed_ok"), "Apply")), easyClose = FALSE))
    }
    observeEvent(input$edit_seq, open_edit_seq())
    observeEvent(input$edit_seq2, open_edit_seq())
    observeEvent(input$ed_ok, {
      i <- rec_i(); s <- toupper(gsub("[^A-Za-z*.-]", "", .nz(input$ed_text, "")))
      if (!nzchar(s)) { showNotification("The sequence is empty.", type = "error"); return() }
      if (mutate(function(d) {
        q <- d$recs$qual[i]
        if (!is.na(q) && nchar(s) != nchar(q))
          stop("This is a FASTQ read: keep the same number of bases (you can replace bases, not add or remove them). Use Crop to shorten it.")
        d$recs$seq[i] <- s
        if (!is.null(d$gb)) d$gb[[i]]$seq <- s
        rebuild(d, d$recs)
      })) { removeModal(); if (!is.null(tryCatch(recs()$gb, error = function(e) NULL)))
        showNotification("GenBank feature positions are not updated by sequence edits.", type = "warning") }
    })
    observeEvent(input$crop_seq, {
      d <- tryCatch(recs(), error = function(e) NULL); req(d); i <- rec_i(); n <- nchar(d$recs$seq[i])
      showModal(modalDialog(title = paste("Crop:", d$recs$id[i]),
        fluidRow(column(6, numericInput(ns("crop_from"), "Keep from base", 1, 1, n, 1)), column(6, numericInput(ns("crop_to"), "to base", n, 1, n, 1))),
        helpText(sprintf("The sequence has %d bases. Only the chosen range is kept.", n)),
        footer = tagList(modalButton("Cancel"), actionButton(ns("crop_ok"), "Crop")), easyClose = TRUE))
    })
    observeEvent(input$crop_ok, {
      i <- rec_i(); a <- as.integer(.nz(input$crop_from, 1)); b <- as.integer(.nz(input$crop_to, 1))
      if (mutate(function(d) {
        n <- nchar(d$recs$seq[i])
        if (a < 1L || b > n || a > b) stop(sprintf("Choose a range inside 1-%d.", n))
        d$recs$seq[i] <- substr(d$recs$seq[i], a, b)
        if (!is.na(d$recs$qual[i])) d$recs$qual[i] <- substr(d$recs$qual[i], a, b)
        if (!is.null(d$gb)) d$gb[[i]]$seq <- d$recs$seq[i]
        rebuild(d, d$recs)
      })) removeModal()
    })
    observeEvent(input$add_seq, {
      showModal(modalDialog(title = "Add sequences", size = "l",
        div(class = "op-ed", textAreaInput(ns("add_text"), NULL, rows = 10, width = "100%", placeholder = ">name\nACGT...")),
        helpText("Paste FASTA (one or more sequences). Without a > line the text is added as one sequence named \"new_sequence\"."),
        footer = tagList(modalButton("Cancel"), actionButton(ns("add_ok"), "Add")), easyClose = FALSE))
    })
    observeEvent(input$add_ok, {
      txt <- .nz(input$add_text, ""); if (!nzchar(trimws(txt))) { showNotification("Nothing to add.", type = "warning"); return() }
      tf <- tempfile(fileext = ".fa"); on.exit(unlink(tf)); writeLines(txt, tf)
      nr <- tryCatch(.op_read_fasta(tf, "new_sequence.fa"), error = function(e) NULL)
      if (is.null(nr) || !nrow(nr)) { showNotification("No sequence found in the pasted text.", type = "error"); return() }
      if (add_recs(nr)) removeModal()
    })
    add_recs <- function(nr) {
      mutate(function(d) {
        if (any(!is.na(d$recs$qual))) stop("Sequences cannot be added to a FASTQ file.")
        if (!is.null(d$gb)) stop("Sequences cannot be added to a GenBank file. Save it as FASTA first (Save as), then open that file.")
        r <- rebuild(d, rbind(d$recs, nr)); r$aligned <- length(unique(nchar(r$recs$seq))) == 1L; r
      })
    }
    # Qt window, File > Import: sequences from files or from the clipboard are added to the alignment that is open
    observeEvent(input$qt_import, {
      x <- input$qt_import; nr <- list(); bad <- character()
      for (p in unlist(x$paths)) {
        if (!file.exists(p) || dir.exists(p)) next
        r <- tryCatch({ ty <- .op_detect(p)
          switch(ty, fasta = , plain = .op_read_fasta(p), fastq = .op_read_fastq(p),
                 genbank = { gb <- .op_read_genbank(p); .op_recs(vapply(gb, function(g) g$id, ""), vapply(gb, function(g) g$definition, ""), vapply(gb, function(g) g$seq, "")) },
                 NULL) }, error = function(e) NULL)
        if (is.null(r) || !nrow(r)) bad <- c(bad, basename(p)) else { r$qual <- NA_character_; nr[[length(nr) + 1L]] <- r }
      }
      if (nzchar(.nz(x$text, ""))) {
        tf <- tempfile(fileext = ".fa"); writeLines(x$text, tf)
        r <- tryCatch(.op_read_fasta(tf, "clipboard.fa"), error = function(e) NULL); unlink(tf)
        if (is.null(r) || !nrow(r)) bad <- c(bad, "clipboard") else { r$qual <- NA_character_; nr[[length(nr) + 1L]] <- r }
      }
      if (length(bad)) showNotification(paste("No sequences found in:", paste(bad, collapse = ", ")), type = "warning", duration = 7)
      if (!length(nr)) return()
      nr <- do.call(rbind, nr)
      k <- if (!is.na(rv$sel) && rv$sel <= nrow(rv$files)) rv$files$type[rv$sel] else ""
      if (k %in% c("fasta", "plain") && add_recs(nr))
        showNotification(sprintf("%d sequence(s) added to the open alignment.", nrow(nr)), type = "message", duration = 5)
      else if (!(k %in% c("fasta", "plain"))) {      # nothing suitable is open: the sequences become a new file
        tf <- file.path(up_dir, sprintf("imported_%d.fasta", as.integer(Sys.time()) %% 100000L))
        writeLines(.op_fasta_text(nr$id, nr$seq, nr$desc), tf); add_paths(tf)
      }
    })

    # ===== saving ==========================================================================================
    fmt_choices <- function(type) switch(type,
      ab1 = c("AB1 (original traces + your base calls)" = "ab1", "FASTA (the Sequence tab)" = "fasta", "FASTQ (the Sequence tab, with Phred)" = "fastq"),
      fasta = , plain = c("FASTA" = "fasta", "CSV" = "csv"),
      fastq = c("FASTQ" = "fastq", "FASTA" = "fasta", "CSV" = "csv"),
      genbank = c("FASTA" = "fasta", "CSV" = "csv"),
      c("Copy of the file" = "copy"))
    fmt_ext <- c(ab1 = ".ab1", fasta = ".fasta", fastq = ".fastq", csv = ".csv")
    write_out <- function(path, fmt) {
      x <- cur(); d <- x$data; id <- tools::file_path_sans_ext(x$row$name)
      if (fmt == "copy") { if (!file.copy(x$row$path, path, overwrite = TRUE)) stop("Could not copy the file."); return(invisible()) }
      if (x$type == "ab1") {
        if (fmt == "ab1") return(invisible(.op_ab1_patch(x$row$path, d$raw$seq, path, keep = !d$del)))
        sq <- ab1_seq()
        if (fmt == "fasta") writeLines(.op_fasta_text(id, sq$seq), path)
        else writeLines(c(paste0("@", id), sq$seq, "+", intToUtf8(pmin(pmax(sq$qual, 0L), 93L) + 33L)), path)
        return(invisible())
      }
      r <- d$recs
      switch(fmt,
        fasta = writeLines(.op_fasta_text(r$id, r$seq, r$desc), path),
        fastq = { if (anyNA(r$qual)) stop("These sequences have no quality values, so they cannot be saved as FASTQ.")
                  writeLines(as.vector(rbind(paste0("@", r$id, ifelse(nzchar(r$desc), " ", ""), r$desc), r$seq, "+", r$qual)), path) },
        csv = utils::write.csv(data.frame(Name = r$id, Description = r$desc, Length = d$stats$Length, Sequence = r$seq), path, row.names = FALSE),
        stop("Unknown format."))
      invisible()
    }
    is_local <- function(row) row$where != "(uploaded)"

    observeEvent(input$save, {
      x <- cur()
      if (!is_local(x$row)) { showNotification("This file came through the browser upload, so there is no file on your computer to overwrite. Use Save as, then Download.", type = "warning", duration = 10); return() }
      if (x$type == "genbank") { showNotification("A GenBank file cannot be overwritten here (the features would be lost). Use Save as to keep the sequences as FASTA.", type = "warning", duration = 10); return() }
      if (!isTRUE(x$row$dirty)) { showNotification("No changes to save.", type = "message"); return() }
      showModal(modalDialog(title = "Save", p("Overwrite ", tags$code(x$row$path), "?"),
        if (x$type == "plain") helpText("This file will be written as FASTA."),
        checkboxInput(ns("sv_bak"), "Keep a backup of the old file (name.bak)", TRUE),
        footer = tagList(modalButton("Cancel"), actionButton(ns("sv_ok"), "Save")), easyClose = TRUE))
    })
    observeEvent(input$sv_ok, {
      x <- cur(); i <- rv$sel; fmt <- switch(x$type, ab1 = "ab1", fastq = "fastq", "fasta")
      tmp <- tempfile(); on.exit(unlink(tmp))
      ok <- tryCatch({ write_out(tmp, fmt); TRUE }, error = function(e) { showNotification(conditionMessage(e), type = "error", duration = 10); FALSE })
      if (!ok) return()
      if (isTRUE(input$sv_bak)) file.copy(x$row$path, paste0(x$row$path, ".bak"), overwrite = TRUE)
      if (!file.copy(tmp, x$row$path, overwrite = TRUE)) { showNotification("Could not write the file (is it open in another program, or read-only?).", type = "error", duration = 10); return() }
      if (x$type == "ab1" && exists(x$row$key, envir = cache, inherits = FALSE)) rm(list = x$row$key, envir = cache)   # reload from the saved file
      if (exists(x$row$key, envir = hist, inherits = FALSE)) rm(list = x$row$key, envir = hist)
      rv$files$dirty[i] <- FALSE; rv$files$size[i] <- file.size(x$row$path); rv$ver <- rv$ver + 1L; rv$rel <- rv$rel + 1L
      removeModal(); showNotification("Saved.", type = "message")
    })

    observeEvent(input$save_as, {
      x <- cur(); ch <- fmt_choices(x$type); local <- .sv_local()
      def_dir <- if (is_local(x$row)) x$row$where else path.expand("~")
      showModal(modalDialog(title = "Save as",
        selectInput(ns("sa_fmt"), "Format", ch),
        textInput(ns("sa_name"), "File name", paste0(tools::file_path_sans_ext(x$row$name), if (x$type %in% c("ab1", seq_types)) fmt_ext[[ch[[1]]]] else paste0(".", tools::file_ext(x$row$name)))),
        if (local) div(style = "display:flex;gap:8px;align-items:flex-end",
                       div(style = "flex:1", textInput(ns("sa_dir"), "Folder on this computer", def_dir, width = "100%")),
                       if (.Platform$OS.type == "windows") div(style = "margin-bottom:15px", actionButton(ns("sa_browse"), "Browse...", icon = icon("folder-open")))),
        helpText(if (local) "\"Save to folder\" writes the file on this computer and opens the new file in the list. \"Download\" uses the browser's download." else "The file is downloaded by the browser."),
        footer = tagList(modalButton("Close"), downloadButton(ns("sa_dl"), "Download"),
                         if (local) actionButton(ns("sa_go"), "Save to folder", class = "btn-primary")), easyClose = TRUE))
    })
    observeEvent(input$sa_browse, {
      start <- trimws(.nz(input$sa_dir, ""))
      d <- tryCatch(utils::choose.dir(default = if (dir.exists(start)) gsub("/", "\\\\", start) else "", caption = "Choose the folder to save in"),
                    error = function(e) NA_character_)
      if (length(d) == 1L && !is.na(d) && nzchar(d)) updateTextInput(session, "sa_dir", value = gsub("\\\\", "/", d))
    })
    observeEvent(input$sa_fmt, {
      if (input$sa_fmt %in% names(fmt_ext) || input$sa_fmt %in% fmt_ext)
        updateTextInput(session, "sa_name", value = paste0(tools::file_path_sans_ext(.nz(input$sa_name, "file")), fmt_ext[[input$sa_fmt]]))
    }, ignoreInit = TRUE)
    output$sa_dl <- downloadHandler(function() .nz(input$sa_name, "file"), function(file) write_out(file, input$sa_fmt))
    sa_pending <- reactiveVal(NULL)
    do_save_as <- function(target) {
      tmp <- tempfile(); on.exit(unlink(tmp))
      ok <- tryCatch({ write_out(tmp, input$sa_fmt); TRUE }, error = function(e) { showNotification(conditionMessage(e), type = "error", duration = 10); FALSE })
      if (!ok) return()
      if (!file.copy(tmp, target, overwrite = TRUE)) { showNotification("Could not write the file there.", type = "error", duration = 10); return() }
      removeModal(); showNotification(paste("Saved:", target), type = "message", duration = 8)
      j <- match(target, rv$files$key)
      if (!is.na(j)) {                       # replaced a file that is open: reload it
        for (e in list(cache, hist)) if (exists(target, envir = e, inherits = FALSE)) rm(list = target, envir = e)
        rv$files$dirty[j] <- FALSE; rv$files$size[j] <- file.size(target); rv$ver <- rv$ver + 1L; rv$sver <- rv$sver + 1L; rv$rel <- rv$rel + 1L
      } else add_paths(target)
    }
    observeEvent(input$sa_go, {
      nm <- trimws(.nz(input$sa_name, "")); dir <- trimws(gsub("^\"|\"$", "", .nz(input$sa_dir, "")))
      if (!nzchar(nm) || grepl("[\\\\/:*?\"<>|]", nm)) { showNotification("Give a valid file name (no \\ / : * ? \" < > |).", type = "error"); return() }
      if (!nzchar(dir)) { showNotification("Give a folder.", type = "error"); return() }
      if (!dir.exists(dir) && !dir.create(dir, recursive = TRUE, showWarnings = FALSE)) { showNotification("Could not create that folder.", type = "error"); return() }
      target <- normalizePath(file.path(dir, nm), winslash = "/", mustWork = FALSE)
      if (file.exists(target)) {
        sa_pending(target)
        confirm("File exists", paste("Replace", target, "?"), "sa_force", "Replace")
      } else do_save_as(target)
    })
    observeEvent(input$sa_force, { t <- sa_pending(); req(t); do_save_as(t) })

    # ===== rename / delete a file ==========================================================================
    observeEvent(input$rename, {
      i <- rv$sel; req(!is.na(i), i <= nrow(rv$files))
      showModal(modalDialog(title = "Rename file", textInput(ns("rn_name"), "New file name", rv$files$name[i]),
        helpText(if (is_local(rv$files[i, ])) "The file is renamed on your computer." else "This file was uploaded through the browser: only its name in this list changes."),
        footer = tagList(modalButton("Cancel"), actionButton(ns("rn_ok"), "Rename")), easyClose = TRUE))
    })
    observeEvent(input$rn_ok, {
      i <- rv$sel; req(!is.na(i), i <= nrow(rv$files)); nm <- trimws(.nz(input$rn_name, ""))
      if (!nzchar(nm) || grepl("[\\\\/:*?\"<>|]", nm)) { showNotification("Give a valid file name (no \\ / : * ? \" < > |).", type = "error"); return() }
      f <- rv$files[i, ]; old <- f$key; new <- old
      if (is_local(f)) {
        np <- file.path(dirname(f$path), nm)
        if (file.exists(np)) { showNotification("A file with that name already exists.", type = "error"); return() }
        if (!file.rename(f$path, np)) { showNotification("Could not rename the file (is it open in another program?).", type = "error"); return() }
        new <- normalizePath(np, winslash = "/", mustWork = FALSE)
        rv$files$path[i] <- new; rv$files$key[i] <- new
        for (e in list(cache, hist)) if (exists(old, envir = e, inherits = FALSE)) { assign(new, get(old, envir = e), envir = e); rm(list = old, envir = e) }
      }
      rv$files$name[i] <- nm
      if (exists(new, envir = cache, inherits = FALSE)) {
        d <- get(new, envir = cache)
        if (!is.null(d$raw)) { d$raw$name <- nm; assign(new, d, envir = cache) }
      }
      rv$ver <- rv$ver + 1L; removeModal()
    })
    observeEvent(input$del_disk, {
      i <- rv$sel; req(!is.na(i), i <= nrow(rv$files)); f <- rv$files[i, ]
      if (!is_local(f)) { showNotification("This file was uploaded through the browser, so there is nothing on your computer to delete. Use \"Close this file\".", type = "warning", duration = 8); return() }
      confirm("Delete from disk", tagList(p("Delete this file from your computer?"), tags$code(f$path), p(tags$b("This cannot be undone."))),
              "del_ok", "Delete the file")
    })
    observeEvent(input$del_ok, {
      i <- rv$sel; req(!is.na(i), i <= nrow(rv$files)); f <- rv$files[i, ]
      removeModal()
      if (file.remove(f$path)) { drop_current(); showNotification(paste("Deleted:", f$name), type = "message") }
      else showNotification("Could not delete the file (is it open in another program?).", type = "error")
    })

    seq_types <- c("fasta", "plain", "fastq", "genbank")
    recs <- reactive({
      x <- cur(); req(x$type %in% seq_types, isTRUE(x$data$ok)); x$data
    })
    rec_i <- reactive({
      d <- recs(); i <- input$recs_tbl_rows_selected
      if (length(i) && i[1] <= nrow(d$recs)) as.integer(i[1]) else 1L
    })

    # ---- the main area depends on the kind of file -----------------------------------------------------------
    output$kind_out <- renderText(kind())
    outputOptions(output, "kind_out", suspendWhenHidden = FALSE)
    observe({                                  # GenBank / FASTQ only tabs
      k <- kind()
      for (t in c("gbinfo", "gbfeat")) (if (k == "genbank") showTab else hideTab)("seq_tabs", t, session = session)
      (if (k == "fastq") showTab else hideTab)("seq_tabs", "fq", session = session)
    })

    # ---- the AliView-like editor for sequence files -------------------------------------------------------
    cur_key <- reactive({ i <- rv$sel; f <- rv$files; if (is.na(i) || i > nrow(f)) NA_character_ else f$key[i] })
    skip_push <- FALSE
    push_ed <- function() {
      d <- tryCatch(cur()$data, error = function(e) NULL); req(d, isTRUE(d$ok), !is.null(d$recs))
      r <- d$recs; L <- nchar(r$seq); mx <- max(L)
      session$sendCustomMessage("msaed_load", list(id = ns("ed"), names = unname(r$id),
                                                   seqs = unname(paste0(r$seq, strrep("-", mx - L))), keep = NULL))
    }
    observe({
      k <- kind(); cur_key(); rv$sver; rv$rel
      if (k %in% c("fasta", "fastq", "genbank")) {
        if (skip_push) { skip_push <<- FALSE; return() }
        isolate(push_ed())
      }
    })
    observeEvent(input$ed_state, {           # "Apply edits" in the editor
      x <- tryCatch(cur(), error = function(e) NULL)
      if (is.null(x) || !(x$type %in% seq_types) || !isTRUE(x$data$ok)) return()
      s <- input$ed_state
      nm <- as.character(unlist(s$names)); sq <- as.character(unlist(s$seqs))
      if (!length(nm) || length(nm) != length(sq)) return()
      keep <- strsplit(as.character(.nz(s$keep, "")), "")[[1]] == "1"
      if (length(keep) == nchar(sq[1]) && !all(keep))
        sq <- vapply(strsplit(sq, ""), function(z) paste(z[keep], collapse = ""), "")
      if (!isTRUE(x$data$aligned)) sq <- sub("-+$", "", sq)         # padding added only to show unequal lengths
      if (any(!nzchar(sq))) { showNotification("A sequence became empty; delete that row in the editor instead.", type = "error"); return() }
      skip_push <<- TRUE
      ok <- mutate(function(d) {
        r <- d$recs; m <- match(nm, r$id)
        same <- length(nm) == nrow(r) && identical(nm, r$id)
        nr <- data.frame(id = nm, desc = ifelse(is.na(m), "", r$desc[m]), seq = sq, qual = NA_character_, stringsAsFactors = FALSE)
        if (any(!is.na(r$qual)) && same) {
          q <- r$qual; ok_q <- !is.na(q) & nchar(gsub("[-.]", "", sq)) == nchar(q)
          nr$qual[ok_q] <- q[ok_q]
        }
        gb <- d$gb
        if (!is.null(gb)) { if (same) { for (j in seq_along(gb)) gb[[j]]$seq <- gsub("[-.]", "", sq[j]) } else gb <- NULL }
        rebuild(d, nr, gb)
      })
      if (!isTRUE(ok)) skip_push <<- FALSE
      else {
        nd <- tryCatch(cur()$data, error = function(e) NULL)
        if (!is.null(nd) && any(!is.na(x$data$recs$qual)) && any(is.na(nd$recs$qual)))
          showNotification("Quality values were dropped for sequences whose bases were added or removed.", type = "warning")
        showNotification("Edits applied. Press Save or Save as to write the file.", type = "message", duration = 6)
      }
    })

    # ---- find a sequence in the open file, on both strands (the find bar opens with Ctrl+F; see R/ui_helpers.R) ------------
    find_hits <- sg_find_server(input, output, session, editor_id = ns("ed"), msg_none = "Open a sequence file first.",
      get_aln = function() {
        d <- tryCatch(cur()$data, error = function(e) NULL)
        if (is.null(d) || !isTRUE(d$ok) || is.null(d$recs) || !nrow(d$recs)) return(NULL)
        r <- d$recs; L <- nchar(r$seq); mx <- max(L)
        stats::setNames(paste0(r$seq, strrep("-", mx - L)), r$id)        # the same padded rows the editor shows
      })
    observeEvent(cur_key(), { find_hits(NULL) })

    # ---- the chromatogram viewer for .ab1 files (the same viewer as the Assembly tab) ---------------------------------
    vid <- ns("viewer")
    flagstr <- function(x) paste(as.integer(x), collapse = "")
    observe({
      k <- kind(); cur_key(); rv$rel
      if (identical(k, "ab1")) isolate({
        d <- tryCatch(cur()$data, error = function(e) NULL)
        if (!is.null(d) && isTRUE(d$ok)) {
          r <- d$raw; n <- length(r$qual)
          aln <- list(reads = list(list(traces = r$traces, ploc = r$ploc)), idx = matrix(seq_len(n), 1L))
          session$sendCustomMessage("sv_static", list(
            id = vid, nc = n, nr = 1L, has_ref = FALSE, ref = "", refpos = NULL,
            colors = list(base = as.list(BASE_COLS), trace = as.list(TRACE_COLS)),
            reads = list(list(label = tools::file_path_sans_ext(r$name), base = d$orig, q = I(as.integer(r$qual)), tr = sv_trace(aln, 1L)))))
        }
      })
    })
    observe({
      x <- tryCatch(cur(), error = function(e) NULL)
      req(x, x$type == "ab1", isTRUE(x$data$ok))
      d <- x$data; n <- length(d$raw$qual)
      man <- strsplit(d$raw$seq, "")[[1]] != strsplit(d$orig, "")[[1]]
      session$sendCustomMessage("sv_dyn", list(id = vid, final = d$raw$seq, review = strrep("0", n), manual = flagstr(man),
                                               delcol = flagstr(d$del), mask = I(flagstr(d$del))))
    })
    observeEvent(input$viewer_del, {
      x <- tryCatch(cur(), error = function(e) NULL); req(x, x$type == "ab1")
      a <- suppressWarnings(as.integer(input$viewer_del$from)); b <- suppressWarnings(as.integer(input$viewer_del$to))
      if (mutate(function(d) {
        n <- length(d$del); a <- max(1L, a); b <- min(n, b)
        if (is.na(a) || is.na(b) || a > b) stop("Nothing is selected.")
        ks <- (a:b)[!d$del[a:b]]
        if (!length(ks)) stop("Those bases are already deleted.")
        if (sum(!d$del) - length(ks) < 1L) stop("That would delete the whole read.")
        d$del[ks] <- TRUE; d
      }, structural = FALSE))
        showNotification("Bases deleted (shown struck out). Ctrl+Z undoes; Save writes the shorter read to the file.", type = "message", duration = 6)
    })
    observeEvent(input$viewer_set, {
      x <- tryCatch(cur(), error = function(e) NULL); req(x, x$type == "ab1")
      a <- suppressWarnings(as.integer(input$viewer_set$from)); b <- suppressWarnings(as.integer(input$viewer_set$to))
      tok <- toupper(as.character(input$viewer_set$base))
      mutate(function(d) {
        n <- length(d$del); a <- max(1L, a); b <- min(n, b)
        if (is.na(a) || is.na(b) || a > b) stop("Nothing is selected.")
        s <- strsplit(d$raw$seq, "")[[1]]
        if (identical(tok, "RULE")) { s[a:b] <- strsplit(d$orig, "")[[1]][a:b]; d$del[a:b] <- FALSE }
        else if (grepl("^[ACGTNRYKMSWBDHV]$", tok)) { ks <- (a:b)[!d$del[a:b]]; if (!length(ks)) stop("Those bases are deleted; press Ctrl+Z to restore them."); s[ks] <- tok }
        else stop("Type A C G T N or an IUPAC letter to change a base (0 = back to the original). Use Delete to remove a base.")
        d$raw$seq <- paste(s, collapse = ""); d
      }, structural = FALSE)
    })

    output$msg <- renderUI({
      x <- cur()
      div(class = "alert alert-danger", tags$b(x$row$name), tags$br(), x$data$msg)
    })

    # ===== AB1 =================================================================================================
    ab1 <- reactive({ x <- cur(); req(x$type == "ab1", isTRUE(x$data$ok)); c(x$data, list(row = x$row)) })
    ab1_rd <- reactive({                      # the read without its deleted bases
      d <- ab1(); k <- !d$del; req(any(k)); raw <- d$raw
      raw$seq <- paste(strsplit(raw$seq, "")[[1]][k], collapse = ""); raw$qual <- raw$qual[k]; raw$ploc <- raw$ploc[k]
      .op_ab1_rd(raw, .op_num(input$q_trim, 20), .op_num(input$win, 10))
    })
    ab1_seq <- reactive({
      rd <- ab1_rd(); n <- length(rd$qual)
      whole <- identical(input$ab1_region, "all")
      from <- if (whole) 1L else rd$qc$from; to <- if (whole) n else rd$qc$to
      s <- substr(rd$seq, from, to); q <- rd$qual[from:to]
      rc <- isTRUE(input$ab1_rc)
      if (rc) { s <- revcomp_chr(s); q <- rev(q) }
      list(seq = s, qual = q, from = from, to = to, rc = rc, whole = whole)
    })

    output$ab1_boxes <- renderUI({
      rd <- ab1_rd(); q <- rd$qual; f <- rd$qc$from; t <- rd$qc$to; kq <- q[f:t]
      fluidRow(
        column(3, .infobox("dna", "Bases called", format(length(q), big.mark = ","), "#2c7be5")),
        column(3, .infobox("chart-simple", "Mean Phred (whole read)", sprintf("%.1f", mean(q)), "#6f42c1")),
        column(3, .infobox("arrows-left-right", "Kept after trimming",
                           if (rd$trimmed) sprintf("%d-%d (%d bp)", f, t, t - f + 1L) else "no good region", "#e67e22")),
        column(3, .infobox("circle-check", "Q30 bases (kept part)", sprintf("%.0f%%", 100 * mean(kq >= 30)), "#27ae60")))
    })

    output$phred_plot <- renderPlot({ plot_phred(ab1_rd()) }, res = 96)

    output$ab1_seq_view <- renderUI({
      sq <- ab1_seq()
      tagList(div(class = "op-hdr", sprintf(">%s %s%s", tools::file_path_sans_ext(ab1()$row$name),
                                            if (sq$whole) "whole read" else sprintf("bases %d-%d", sq$from, sq$to),
                                            if (sq$rc) " (reverse complement)" else "")),
              .op_seq_html(sq$seq, .op_num(input$ab1_w, 60), !identical(input$ab1_col, FALSE)))
    })
    ab1_fname <- function(ext) {
      sq <- ab1_seq()
      paste0(tools::file_path_sans_ext(ab1()$row$name), if (!sq$whole) "_trimmed" else "", if (sq$rc) "_rc" else "", ext)
    }
    output$dl_ab1_fa <- downloadHandler(function() ab1_fname(".fasta"), function(file) {
      sq <- ab1_seq()
      writeLines(.op_fasta_text(tools::file_path_sans_ext(ab1()$row$name), sq$seq,
                                if (sq$whole) "whole read" else sprintf("bases %d-%d%s", sq$from, sq$to, if (sq$rc) " reverse complement" else "")), file)
    })
    output$dl_ab1_fq <- downloadHandler(function() ab1_fname(".fastq"), function(file) {
      sq <- ab1_seq()
      writeLines(c(paste0("@", tools::file_path_sans_ext(ab1()$row$name)), sq$seq, "+",
                   intToUtf8(pmin(pmax(sq$qual, 0L), 93L) + 33L)), file)
    })
    output$dl_ab1_csv <- downloadHandler(function() ab1_fname("_bases.csv"), function(file) {
      sq <- ab1_seq(); rd <- ab1_rd()
      idx <- if (sq$rc) rev(sq$from:sq$to) else sq$from:sq$to
      utils::write.csv(data.frame(Base_no = idx, Base = strsplit(sq$seq, "")[[1]], Phred = sq$qual,
                                  Peak_position = rd$ploc[idx]), file, row.names = FALSE)
    })
    output$ab1_info <- DT::renderDT({
      x <- ab1(); rd <- ab1_rd()
      base <- data.frame(
        Item = c("File", "Size", "Bases called (kept)", "Bases deleted (not saved yet)", "Bases changed from the original", "Trace points per channel", "Channel order (dyes)", "Median peak spacing (points)"),
        Value = c(x$row$name, .op_size(x$row$size), as.character(length(rd$qual)), as.character(sum(x$del)),
                  as.character(sum(strsplit(x$raw$seq, "")[[1]] != strsplit(x$orig, "")[[1]])), as.character(nrow(rd$traces)),
                  paste(colnames(rd$traces), collapse = " "), sprintf("%.1f", stats::median(diff(rd$ploc)))),
        stringsAsFactors = FALSE)
      DT::datatable(rbind(base, x$meta), rownames = FALSE, selection = "none",
                    options = list(dom = "t", paging = FALSE, ordering = FALSE))
    })

    # ===== FASTA / FASTQ / GenBank =============================================================================
    output$recs_tbl <- DT::renderDT({
      rv$sel; rv$sver                      # redraw when another file is chosen or the list itself changes (not on a rename)
      d <- isolate(recs())
      tb <- cbind(No = seq_len(nrow(d$stats)), Name = d$recs$id, Description = d$recs$desc,
                  d$stats[, setdiff(names(d$stats), "Name"), drop = FALSE])
      DT::datatable(tb, rownames = FALSE,
                    editable = list(target = "cell", disable = list(columns = setdiff(seq_len(ncol(tb)) - 1L, c(1L, 2L)))),
                    selection = list(mode = "multiple", selected = 1L, target = "row"),
                    options = list(pageLength = 10, scrollX = TRUE, dom = if (nrow(tb) > 10) "ftip" else "t"))
    })
    output$seq_boxes <- renderUI({
      d <- recs(); st <- d$stats[rec_i(), ]
      fluidRow(
        column(3, .infobox("ruler", "Length", format(st$Length, big.mark = ","), "#2c7be5")),
        column(3, .infobox("dna", "GC content", if (is.na(st$GC_pct)) "-" else paste0(st$GC_pct, "%"), "#27ae60")),
        column(3, if (!is.null(st$Mean_Phred)) .infobox("chart-simple", "Mean Phred", sprintf("%.1f", st$Mean_Phred), "#6f42c1")
                  else .infobox("circle-question", "N (unknown) bases", if (is.na(st$N_count)) "-" else st$N_count, "#6f42c1")),
        column(3, .infobox("layer-group", "Type", st$Type, "#e67e22")))
    })
    seq_now <- reactive({           # sequence of the selected record, as the controls ask
      d <- recs(); i <- rec_i(); r <- d$recs[i, ]
      s <- r$seq
      if (!identical(input$seq_nogap, FALSE)) s <- gsub("[-.]", "", s)
      prot <- d$stats$Type[i] == "protein"
      rc <- isTRUE(input$seq_rc) && !prot
      if (rc && nchar(s) <= 300000L) s <- revcomp_chr(toupper(s)) else rc <- FALSE
      list(rec = r, seq = s, prot = prot, rc = rc)
    })
    output$seq_view <- renderUI({
      x <- seq_now(); r <- x$rec
      tagList(div(class = "op-hdr", paste0(">", r$id, if (nzchar(r$desc)) paste0(" ", r$desc) else "",
                                           if (x$rc) " (reverse complement)" else "")),
              .op_seq_html(x$seq, .op_num(input$seq_w, 60), !identical(input$seq_col, FALSE), x$prot))
    })
    output$dl_fa_sel <- downloadHandler(
      function() paste0(gsub("[^A-Za-z0-9._-]", "_", seq_now()$rec$id), ".fasta"),
      function(file) { x <- seq_now(); writeLines(.op_fasta_text(x$rec$id, x$seq, x$rec$desc), file) })
    output$dl_fa_all <- downloadHandler(
      function() paste0(tools::file_path_sans_ext(cur()$row$name), "_all.fasta"),
      function(file) { r <- recs()$recs; writeLines(.op_fasta_text(r$id, r$seq, r$desc), file) })

    output$fq_plot <- renderPlot({
      d <- recs(); q <- d$recs$qual[rec_i()]
      validate(need(!is.na(q) && nzchar(q), "This file has no quality values."))
      v <- utils::head(utf8ToInt(q) - 33L, 5000L)
      df <- data.frame(i = seq_along(v), q = v)
      ggplot(df, aes(i, q)) + geom_col(width = 1, fill = "#2a9d8f") +
        geom_hline(yintercept = c(20, 30), linetype = 3, colour = "grey40") +
        labs(x = "base", y = "Phred", subtitle = sprintf("%s: dotted lines = Q20 and Q30", d$recs$id[rec_i()])) +
        theme_minimal(base_size = 11)
    }, res = 96)

    output$aln_view <- renderUI({
      r <- recs()$recs
      validate(need(nrow(r) >= 2, "Open a file with at least two sequences to see an alignment."))
      validate(need(length(unique(nchar(r$seq))) == 1,
                    "The sequences have different lengths, so this is not an alignment. Send them to the MSA tab to align them."))
      .op_aln_html(stats::setNames(r$seq, r$id), .op_num(input$aln_w, 80), isTRUE(input$aln_dots))
    })

    six <- reactive({
      d <- recs(); i <- rec_i()
      validate(need(d$stats$Type[i] == "nucleotide", "Translation is for nucleotide sequences."))
      validate(need(requireNamespace("Biostrings", quietly = TRUE), "The Biostrings package is needed for translation."))
      nt <- gsub("[-.*]", "", d$recs$seq[i])
      validate(need(nchar(nt) >= 3, "The sequence is too short to translate."))
      translate_six(substr(nt, 1L, 30000L), .nz(input$tr_code, "1"), as.integer(.op_num(input$tr_min, 30)))
    })
    output$tr_view <- renderUI({
      s <- six()
      tagList(helpText("ORFs of at least the chosen length are shaded; green = start (Met), blue = stop (*)."),
              six_frame_html(s, width = 60L))
    })
    output$tr_orfs <- DT::renderDT({
      o <- six()$orfs; req(nrow(o) > 0)
      o$Protein <- ifelse(nchar(o$Protein) > 60, paste0(substr(o$Protein, 1, 57), "..."), o$Protein)
      DT::datatable(o, rownames = FALSE, selection = "none", options = list(dom = "tp", pageLength = 8, scrollX = TRUE))
    })

    gb_rec <- reactive({ d <- recs(); req(!is.null(d$gb)); d$gb[[rec_i()]] })
    output$gb_info <- DT::renderDT({
      g <- gb_rec()
      DT::datatable(data.frame(Item = c("LOCUS", "Definition", "Accession", "Version", "Organism", "Length (bp)", "Features"),
                               Value = c(g$locus, g$definition, g$accession, g$version, g$organism, nchar(g$seq), nrow(g$features)),
                               stringsAsFactors = FALSE),
                    rownames = FALSE, selection = "none", options = list(dom = "t", paging = FALSE, ordering = FALSE))
    })
    output$gb_feat <- DT::renderDT({
      f <- gb_rec()$features
      validate(need(nrow(f) > 0, "This record has no features."))
      DT::datatable(f, rownames = FALSE, selection = "none", options = list(pageLength = 15, scrollX = TRUE, dom = "ftp"))
    })

    # send what is open to the MSA tab
    observeEvent(input$to_msa, {
      if (!is.function(msa_add)) { showNotification("The MSA tab is not connected.", type = "warning"); return() }
      x <- tryCatch(cur(), error = function(e) NULL); req(x)
      r <- if (x$type == "ab1") {
        data.frame(id = tools::file_path_sans_ext(x$row$name), seq = ab1_seq()$seq, stringsAsFactors = FALSE)
      } else if (x$type %in% seq_types) {
        k <- x$data$stats$Type == "nucleotide"
        data.frame(id = x$data$recs$id[k], seq = gsub("[-.]", "", toupper(x$data$recs$seq[k])), stringsAsFactors = FALSE)
      } else NULL
      if (is.null(r) || !nrow(r)) { showNotification("There is no nucleotide sequence to send.", type = "warning"); return() }
      if (nrow(r) > 1000L) { r <- r[1:1000, ]; showNotification("Only the first 1000 sequences were sent.", type = "warning") }
      msa_add(r, "Open"); n <- nrow(r)
      showNotification(sprintf("%d sequence(s) sent to the MSA tab. Press Align there.", n), type = "message", duration = 8)
      if (is.function(go)) go("msa")
    })

    # ===== Newick tree =========================================================================================
    tree_obj <- reactive({ x <- cur(); req(x$type == "newick", isTRUE(x$data$ok)); x$data })
    output$tree_info <- renderUI({
      d <- tree_obj(); t <- d$tree
      helpText(sprintf("%d tips, %d internal nodes; %s; %s%s.", ape::Ntip(t), ape::Nnode(t),
                       if (is.null(t$edge.length)) "no branch lengths" else "with branch lengths",
                       if (ape::is.rooted(t)) "rooted" else "unrooted",
                       if (d$ntrees > 1) sprintf(" (the file has %d trees; the first is shown)", d$ntrees) else ""))
    })
    tree_h <- function() {
      t <- tryCatch(tree_obj()$tree, error = function(e) NULL)
      if (is.null(t)) return(60L)
      if (.nz(input$tree_layout, "phylogram") %in% c("phylogram", "cladogram"))
        as.integer(max(380, 22 * ape::Ntip(t) * max(0.6, .op_num(input$tree_cex, 0.9)) + 80)) else 640L
    }
    output$tree_plot <- renderPlot({
      t <- tree_obj()$tree; cex <- .op_num(input$tree_cex, 0.9); type <- .nz(input$tree_layout, "phylogram")
      op <- graphics::par(mar = c(2, 1, 1, 1)); on.exit(graphics::par(op))
      ape::plot.phylo(t, type = type, use.edge.length = !identical(type, "cladogram"), cex = cex, edge.width = 1.6,
                      label.offset = if (type %in% c("phylogram", "cladogram")) 0.002 else 0)
      if (isTRUE(input$tree_sup) && !is.null(t$node.label)) {
        lab <- as.character(t$node.label); lab[is.na(lab)] <- ""
        if (any(nzchar(lab))) ape::nodelabels(lab, frame = "none", adj = c(1.15, -0.35), cex = cex * 0.75, col = "#c0392b")
      }
      if (identical(type, "phylogram") && !is.null(t$edge.length)) ape::add.scale.bar(cex = cex * 0.8)
    }, height = function() tree_h(), res = 96)
    output$dl_tree <- downloadHandler(function() paste0(tools::file_path_sans_ext(cur()$row$name), "_copy.nwk"),
                                      function(file) ape::write.tree(tree_obj()$tree, file))

    # ===== table and text ===================================================================================
    tbl <- reactive({ x <- cur(); req(x$type == "table", isTRUE(x$data$ok)); x$data })
    output$table_info <- renderUI({
      d <- tbl(); helpText(sprintf("%s rows, %d columns; separator: %s.", format(nrow(d$df), big.mark = ","), ncol(d$df),
                                   switch(d$sep, "\t" = "tab", "," = "comma", ";" = "semicolon")))
    })
    output$table_view <- DT::renderDT({
      DT::datatable(tbl()$df, rownames = FALSE, filter = "top", options = list(pageLength = 15, scrollX = TRUE))
    })
    txt <- reactive({ x <- cur(); req(x$type == "text", isTRUE(x$data$ok)); x$data$lines })
    output$text_info <- renderUI({
      helpText(sprintf("%s lines shown%s.", format(length(txt()), big.mark = ","),
                       if (length(txt()) >= 3000L) " (the first 3000 only)" else ""))
    })
    output$text_view <- renderText(paste(txt(), collapse = "\n"))

    # For other tabs (Barcode gap): every open sequence file whose records are aligned (3 or more, all of the same length),
    # as list(file name = named character vector). Files are read once and kept in the cache.
    alns <- reactive({
      rv$ver; rv$sver
      f <- rv$files; out <- list()
      for (i in seq_len(nrow(f))) {
        if (!(f$type[i] %in% seq_types)) next
        key <- f$key[i]
        d <- tryCatch({
          if (!exists(key, envir = cache, inherits = FALSE)) assign(key, .op_load(f$path[i], f$name[i], f$type[i]), envir = cache)
          get(key, envir = cache, inherits = FALSE)
        }, error = function(e) NULL)
        if (is.null(d) || !isTRUE(d$ok) || is.null(d$recs) || nrow(d$recs) < 3L) next
        s <- gsub("U", "T", gsub("[.?]", "-", toupper(as.character(d$recs$seq))))
        if (length(unique(nchar(s))) != 1L) next
        out[[f$name[i]]] <- stats::setNames(s, make.unique(as.character(d$recs$id)))
      }
      out
    })
    list(aln = alns)
  })
}
