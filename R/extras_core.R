# ==============================================================================
# Sangerly - sangeranalyseR-style extras (no Shiny)
#   * genetic distance between reads + dendrogram
#   * translation of the consensus and stop-codon detection
#   * stand-alone HTML report
# Packages: Biostrings (getGeneticCode only), htmltools (installed with shiny)
# ==============================================================================

GENETIC_CODES <- c("Standard (1)" = "1",
                   "Vertebrate mitochondrial (2)" = "2",
                   "Yeast mitochondrial (3)" = "3",
                   "Mold, protozoan, coelenterate mitochondrial; mycoplasma/spiroplasma (4)" = "4",
                   "Invertebrate mitochondrial (5)" = "5",
                   "Ciliate, dasycladacean, hexamita nuclear (6)" = "6",
                   "Echinoderm and flatworm mitochondrial (9)" = "9",
                   "Euplotid nuclear (10)" = "10",
                   "Bacterial, archaeal, plant plastid (11)" = "11",
                   "Alternative yeast nuclear (12)" = "12",
                   "Ascidian mitochondrial (13)" = "13",
                   "Alternative flatworm mitochondrial (14)" = "14",
                   "Blepharisma nuclear (15)" = "15",
                   "Chlorophycean mitochondrial (16)" = "16",
                   "Trematode mitochondrial (21)" = "21",
                   "Scenedesmus obliquus mitochondrial (22)" = "22",
                   "Pterobranchia mitochondrial (24)" = "24")

# ------------------------------------------------------------------------------
# Genetic distance (p-distance) between reads over commonly covered columns
# ------------------------------------------------------------------------------
aln_dist <- function(aln, mask = NULL, min_overlap = 30) {
  B <- aln_eff_base(aln, mask); nr <- nrow(B)
  D <- matrix(NA_real_, nr, nr, dimnames = list(aln$labels, aln$labels))
  ok <- B %in% c("A", "C", "G", "T")
  dim(ok) <- dim(B)
  for (i in seq_len(nr)) for (j in i:nr) {
    if (i == j) { D[i, j] <- 0; next }
    both <- ok[i, ] & ok[j, ]
    if (sum(both) >= min_overlap) D[i, j] <- D[j, i] <- mean(B[i, both] != B[j, both])
  }
  D
}

dist_to_hclust <- function(D) {
  if (nrow(D) < 3) return(NULL)
  Dm <- D
  mx <- suppressWarnings(max(Dm, na.rm = TRUE)); if (!is.finite(mx) || mx <= 0) mx <- 0.1
  Dm[is.na(Dm)] <- mx * 1.5
  stats::hclust(stats::as.dist(Dm), method = "average")
}

# ------------------------------------------------------------------------------
# Translation and stop codons
# ------------------------------------------------------------------------------
translate_cf <- function(cf, code_id = "1", frame = 1L) {
  code <- Biostrings::getGeneticCode(as.character(code_id))
  cols <- which(cf$Final != "-")
  s <- paste(cf$Final[cols], collapse = "")
  n <- nchar(s)
  empty_stops <- data.frame(Aa_pos = integer(), Nt_start = integer(), Col = integer(),
                            RefPos = integer(), Codon = character(), stringsAsFactors = FALSE)
  if (n - frame + 1L < 3L) return(list(aa = "", stops = empty_stops, n_x = 0L, n_codons = 0L))
  starts <- seq(frame, n - 2L, by = 3L)
  cod <- substring(s, starts, starts + 2L)
  aa <- unname(code[cod]); aa[is.na(aa)] <- "X"
  st <- which(aa == "*")
  stops <- if (length(st)) {
    data.frame(Aa_pos = st, Nt_start = starts[st], Col = cols[starts[st]],
               RefPos = cf$RefPos[cols[starts[st]]], Codon = cod[st], stringsAsFactors = FALSE)
  } else empty_stops
  list(aa = paste(aa, collapse = ""), stops = stops, n_x = sum(aa == "X"), n_codons = length(aa))
}

best_frame <- function(cf, code_id = "1") {
  ns <- vapply(1:3, function(f) nrow(translate_cf(cf, code_id, f)$stops), 1L)
  which.min(ns)
}

# ------------------------------------------------------------------------------
# Stand-alone HTML report (counterpart of generateReport in sangeranalyseR)
# ------------------------------------------------------------------------------
.html_table <- function(df, max_rows = 500) {
  if (is.null(df) || !nrow(df)) return("<p><em>(none)</em></p>")
  extra <- if (nrow(df) > max_rows) sprintf("<p><em>Showing %d of %d rows.</em></p>", max_rows, nrow(df)) else ""
  df <- utils::head(df, max_rows)
  df[] <- lapply(df, function(x) ifelse(is.na(x), "", as.character(x)))
  hdr <- paste0("<tr>", paste0("<th>", htmltools::htmlEscape(names(df)), "</th>", collapse = ""), "</tr>")
  body <- vapply(seq_len(nrow(df)), function(i)
    paste0("<tr>", paste0("<td>", htmltools::htmlEscape(unlist(df[i, ], use.names = FALSE)), "</td>", collapse = ""), "</tr>"), "")
  paste0("<table>", hdr, paste(body, collapse = ""), "</table>", extra)
}

build_report_html <- function(aln, cf, contigs, qc, snps, indels, diffs, stops = NULL,
                              n_deleted_bases = 0L, n_deleted_cols = 0L, excluded = character(),
                              title = "Sangerly report") {
  esc <- htmltools::htmlEscape
  P <- aln$P
  pm <- data.frame(
    Parameter = c("Trimming method", "Trimming threshold", "Good Phred", "Poor Phred",
                  "Signal ratio cutoff", "Type 4 handling", "Deleted bases", "Deleted consensus positions", "Excluded reads"),
    Value = c(if (identical(P$trim_method, "M1")) "M1 modified Mott" else "M2 sliding window",
              if (identical(P$trim_method, "M1")) sprintf("cutoff %g", P$mott_cutoff) else sprintf("Q>=%g, window %g", P$q_trim, P$win),
              P$q_good, P$q_poor, P$ratio, P$type4, n_deleted_bases, n_deleted_cols,
              if (length(excluded)) paste(excluded, collapse = ", ") else "-"),
    stringsAsFactors = FALSE)
  types <- as.data.frame(table(Type = cf$Type), stringsAsFactors = FALSE); names(types)[2] <- "Columns"
  dec <- cf[cf$Review | cf$Manual | cf$Type == "ColDel",
            c("Col", "RefPos", "Cov", "Type", "Auto", "Final", "Source", "Recommendation", "Reason")]
  fasta <- paste0(">", names(contigs), "\n", unlist(contigs), collapse = "\n")
  css <- paste("body{font-family:system-ui,Segoe UI,Arial,sans-serif;max-width:1100px;margin:24px auto;padding:0 16px;color:#222}",
               "h1{color:#1f4e79}h2{border-bottom:2px solid #d6e4f0;padding-bottom:4px;margin-top:32px;color:#2c5f8a}",
               "table{border-collapse:collapse;font-size:12.5px;margin:8px 0}th,td{border:1px solid #ccd;padding:3px 7px;text-align:left}",
               "th{background:#eef3f8}pre{background:#f6f8fa;padding:10px;overflow:auto;white-space:pre-wrap;word-break:break-all}",
               ".note{background:#fff6dd;padding:8px 12px;border-radius:4px}", sep = "\n")
  parts <- c(
    sprintf("<!doctype html><html><head><meta charset='utf-8'><title>%s</title><style>%s</style></head><body>", esc(title), css),
    sprintf("<h1>%s</h1><p>Generated %s by Sangerly. Mode: %s. %d reads, %d alignment columns.</p>",
            esc(title), format(Sys.time(), "%Y-%m-%d %H:%M"),
            if (identical(aln$mode, "ref")) paste0("mapped to reference ", esc(substr(aln$ref$name, 1, 80))) else "Forward-Reverse overlap, no reference",
            nrow(aln$base), ncol(aln$base)),
    if (length(aln$notes)) paste0("<div class='note'><strong>Notes:</strong><br>", paste(esc(aln$notes), collapse = "<br>"), "</div>"),
    "<h2>Consensus</h2>", sprintf("<pre>%s</pre>", esc(fasta)),
    "<h2>Parameters</h2>", .html_table(pm),
    "<h2>Per-read QC</h2>", .html_table(qc),
    "<h2>Decision types</h2>", .html_table(types),
    "<h2>Reviewed, edited and deleted positions</h2>",
    "<p><code>Source</code> is rule, manual or deleted. Type ColDel = position deleted from the consensus.</p>",
    .html_table(dec),
    "<h2>SNPs relative to the reference</h2>", .html_table(snps),
    "<h2>Indels</h2>", .html_table(indels),
    "<h2>Read differences from the consensus</h2>", .html_table(diffs),
    if (!is.null(stops)) c("<h2>Translation</h2>", sprintf("<p>%s</p>", esc(stops$desc)), sprintf("<pre>%s</pre>", esc(stops$aa)),
                           "<h3>Open reading frames (six-frame translation)</h3>", .html_table(stops$table)),
    "</body></html>")
  paste(parts, collapse = "\n")
}

# ==============================================================================
# Six-frame translation (like ExPASy Translate): 5'3' frames 1-3 and 3'5' frames 1-3
# ==============================================================================
translate_six <- function(nt, code_id = "1", min_aa = 30L) {
  nt <- gsub("U", "T", toupper(gsub("[^A-Za-z]", "", nt)))
  n <- nchar(nt)
  code <- tryCatch(Biostrings::getGeneticCode(as.character(code_id)),
                   error = function(e) stop("Genetic code ", code_id,
                     " is not available in the installed Biostrings; update Bioconductor."))
  rc <- if (n > 0) revcomp_chr(nt) else ""
  frames <- list(); orfs <- list()
  for (strand in c("+", "-")) for (f in 1:3) {
    s <- if (strand == "+") nt else rc
    nm <- sprintf("%s Frame %d", if (strand == "+") "5'3'" else "3'5'", f)
    fr <- list(name = nm, aa = character(), left = integer(), mid = integer(), right = integer(),
               cls = integer(), inorf = logical())
    if (n - f + 1L >= 3L) {
      st <- seq(f, n - 2L, by = 3L)
      aa <- unname(code[substring(s, st, st + 2L)]); aa[is.na(aa)] <- "X"
      # nucleotide span of every codon in FORWARD coordinates of the input sequence
      if (strand == "+") { left <- st; right <- st + 2L } else { left <- n - st - 1L; right <- n - st + 1L }
      cls <- ifelse(aa == "*", 3L, 0L)          # 0 plain, 1 inside ORF, 2 start Met, 3 stop codon
      inorf <- rep(FALSE, length(aa))
      na <- length(aa); i <- 1L
      while (i <= na) {
        if (aa[i] != "M") { i <- i + 1L; next }
        j <- i
        while (j <= na && aa[j] != "*") j <- j + 1L
        stop_found <- j <= na
        e <- if (stop_found) j else na
        len <- e - i + 1L - as.integer(stop_found)
        if (len >= min_aa) {
          inorf[i:e] <- TRUE
          idx <- i:e; idx <- idx[cls[idx] != 3L]; cls[idx] <- 1L; cls[i] <- 2L
          a <- if (strand == "+") left[i] else left[e]
          b <- if (strand == "+") right[e] else right[i]
          orfs[[length(orfs) + 1L]] <- data.frame(
            Frame = nm, Start_nt = a, End_nt = b, Length_aa = len,
            Stop = if (stop_found) "yes" else "no (runs off the end)",
            Protein = paste(aa[i:(e - as.integer(stop_found))], collapse = ""), stringsAsFactors = FALSE)
        }
        i <- e + 1L
      }
      fr <- list(name = nm, aa = aa, left = left, mid = left + 1L, right = right, cls = cls, inorf = inorf)
    }
    frames[[length(frames) + 1L]] <- fr
  }
  cols <- c("ORF", "Frame", "Start_nt", "End_nt", "Length_aa", "Stop", "Protein")
  o <- if (length(orfs)) do.call(rbind, orfs) else
    data.frame(Frame = character(), Start_nt = integer(), End_nt = integer(), Length_aa = integer(),
               Stop = character(), Protein = character(), stringsAsFactors = FALSE)
  o <- o[order(-o$Length_aa), , drop = FALSE]
  o$ORF <- if (nrow(o)) paste0("ORF", seq_len(nrow(o))) else character()
  list(n = n, nt = nt, frames = frames, orfs = o[, cols, drop = FALSE], code_id = code_id, min_aa = min_aa)
}

# Lines of the six-frame layout: three forward frames, the sequence, three reverse frames
.six_lines <- function(six) {
  n <- six$n; nt <- strsplit(six$nt, "")[[1]]
  L <- lapply(six$frames, function(fr) {
    ch <- rep(" ", n); cl <- rep(0L, n)
    if (length(fr$aa)) {
      w <- which(fr$inorf)
      if (length(w)) cl[unlist(lapply(w, function(i) fr$left[i]:fr$right[i]))] <- 1L
      ch[fr$mid] <- fr$aa
      cl[fr$mid] <- fr$cls
    }
    list(label = fr$name, ch = ch, cl = cl)
  })
  c(L[1:3], list(list(label = "Sequence", ch = nt, cl = rep(0L, n))), L[4:6])
}

six_frame_text <- function(six, width = 60L) {
  if (six$n == 0) return("(no sequence)")
  Ls <- .six_lines(six); n <- six$n; out <- character()
  for (a in seq(1L, n, by = width)) {
    b <- min(n, a + width - 1L)
    out <- c(out, sprintf("%s %d", formatC("", width = 13), a))
    for (L in Ls) {
      line <- sprintf("%s %s", formatC(L$label, width = -13), paste(L$ch[a:b], collapse = ""))
      if (identical(L$label, "Sequence")) line <- sprintf("%s %d", line, b)
      out <- c(out, line)
    }
    out <- c(out, "")
  }
  c(sprintf("Six-frame translation, genetic code %s; ORFs >= %d aa start with Met and end at a stop (*).",
            six$code_id, six$min_aa), "", out)
}

# ExPASy-like HTML: ORFs shaded, start Met green, stop codons (*) blue
six_frame_html <- function(six, width = 60L, max_nt = 6000L) {
  if (six$n == 0) return(htmltools::HTML("<p><em>(no sequence)</em></p>"))
  Ls <- .six_lines(six); n <- min(six$n, max_nt)
  seg <- function(ch, cl) {
    r <- rle(cl); ends <- cumsum(r$lengths); starts <- ends - r$lengths + 1L
    paste(vapply(seq_along(r$values), function(k) {
      txt <- paste(ch[starts[k]:ends[k]], collapse = "")
      if (r$values[k] == 0L) txt else sprintf("<span class='sx%d'>%s</span>", r$values[k], txt)
    }, ""), collapse = "")
  }
  out <- character()
  for (a in seq(1L, n, by = width)) {
    b <- min(n, a + width - 1L)
    out <- c(out, sprintf("%s %d", formatC("", width = 13), a))
    for (L in Ls) {
      line <- sprintf("%s %s", htmltools::htmlEscape(formatC(L$label, width = -13)), seg(L$ch[a:b], L$cl[a:b]))
      if (identical(L$label, "Sequence")) line <- sprintf("%s %d", line, b)
      out <- c(out, line)
    }
    out <- c(out, "")
  }
  htmltools::HTML(paste0("<pre class='sx'>", paste(out, collapse = "\n"), "</pre>",
                         if (six$n > max_nt) sprintf("<p><em>Showing the first %d nt only.</em></p>", max_nt) else ""))
}

# ==============================================================================
# Cropping the consensus by primers
# ==============================================================================
.mask4 <- function(x) {
  m <- c(A = 1L, C = 2L, G = 4L, T = 8L, R = 5L, Y = 10L, S = 6L, W = 9L, K = 12L, M = 3L,
         B = 14L, D = 13L, H = 11L, V = 7L, N = 15L)
  r <- unname(m[x]); r[is.na(r)] <- 15L; r
}
.clean_primer <- function(x) {
  if (is.null(x) || is.na(x)) return("")
  gsub("U", "T", toupper(gsub("[^A-Za-z]", "", x)))
}

# Best occurrence of a primer (IUPAC aware) allowing up to max_mm mismatches; NULL if none
find_primer <- function(nt, primer, max_mm = 2L) {
  s <- .mask4(strsplit(nt, "")[[1]]); p <- .mask4(strsplit(primer, "")[[1]])
  n <- length(s); L <- length(p)
  if (L == 0L || n < L) return(NULL)
  mm <- vapply(seq_len(n - L + 1L), function(i) sum(bitwAnd(s[i:(i + L - 1L)], p) == 0L), 1L)
  best <- which(mm == min(mm))[1]
  if (mm[best] > max_mm) return(NULL)
  list(start = best, end = best + L - 1L, mismatches = mm[best])
}

# Crop to the region defined by the forward primer and the reverse primer (given 5'->3').
# If the primers are found on the opposite strand the sequence is reverse-complemented.
crop_by_primers <- function(nt, fwd = "", rev = "", remove = TRUE, max_mm = 2L) {
  fwd <- .clean_primer(fwd); rev <- .clean_primer(rev)
  out <- list(seq = nt, notes = character(), changed = FALSE)
  if (!nzchar(fwd) && !nzchar(rev)) return(out)
  rc_rev <- if (nzchar(rev)) revcomp_chr(rev) else ""
  look <- function(s) list(f = if (nzchar(fwd)) find_primer(s, fwd, max_mm),
                           r = if (nzchar(rc_rev)) find_primer(s, rc_rev, max_mm))
  nf <- function(h) as.integer(!is.null(h$f)) + as.integer(!is.null(h$r))
  want <- as.integer(nzchar(fwd)) + as.integer(nzchar(rev))
  hit <- look(nt); notes <- character()
  if (nf(hit) < want) {
    rcs <- revcomp_chr(nt); h2 <- look(rcs)
    if (nf(h2) > nf(hit)) {
      nt <- rcs; hit <- h2
      notes <- c(notes, "Primers matched the opposite strand: the consensus was reverse-complemented.")
    }
  }
  from <- 1L; to <- nchar(nt)
  if (!is.null(hit$f)) {
    from <- if (remove) hit$f$end + 1L else hit$f$start
    notes <- c(notes, sprintf("Forward primer found at %d-%d (%d mismatch%s).", hit$f$start, hit$f$end,
                              hit$f$mismatches, if (hit$f$mismatches == 1) "" else "es"))
  } else if (nzchar(fwd)) notes <- c(notes, sprintf("Forward primer not found (max %d mismatches): start not cropped.", max_mm))
  if (!is.null(hit$r)) {
    to <- if (remove) hit$r$start - 1L else hit$r$end
    notes <- c(notes, sprintf("Reverse primer (as reverse complement) found at %d-%d (%d mismatch%s).",
                              hit$r$start, hit$r$end, hit$r$mismatches, if (hit$r$mismatches == 1) "" else "es"))
  } else if (nzchar(rev)) notes <- c(notes, sprintf("Reverse primer not found (max %d mismatches): end not cropped.", max_mm))
  if (from > to) {
    out$notes <- c(notes, "The primers overlap or are in the wrong order: no cropping applied.")
    return(out)
  }
  out$seq <- substr(nt, from, to); out$notes <- notes
  out$changed <- from > 1L || to < nchar(nt)
  out
}
