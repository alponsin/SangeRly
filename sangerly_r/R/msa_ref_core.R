# SangeRly - reference-guided multiple alignment (like ViralMSA) and the list of mutations (no Shiny).
#
# msa_refalign(): works like ViralMSA. With minimap2 installed it follows ViralMSA.py step by step (minimap2 + CIGAR, primary
#                 alignments only, unmapped sequences left out, optional 2-pass insertion strategy). Without minimap2 a built-in
#                 aligner does the same job. Every sequence is aligned to ONE reference (both strands are tried, the better one is kept) and the
#                 result is an alignment in the numbering of the reference. Insertions relative to the reference are
#                 left out (like ViralMSA) or kept as extra columns. Up to about 7 kb the alignment is exact; longer sequences
#                 (viral genomes, 30 kb or more) are placed with shared k-mers and only the stretches between them are aligned.
# .aln_variants(): reads ANY alignment against one of its sequences and lists substitutions, insertions, deletions
#                 and ambiguous bases. The ends of a sequence that does not reach the ends of the reference are
#                 not counted as deletions.

# ---- anchored alignment for long sequences (viral genomes of 30 kb or more) ----------------------------------------------------
# A full pairwise table (length x length) is too big for a genome, so a long sequence is placed on the reference with exact
# k-mers that occur ONCE in both sequences ("anchors"), the anchors are chained in order, and only the short stretches between
# the anchors (and the two ends) are aligned with the exact aligner.

# k-mers (ACGT only) that occur once in s, and their start positions
.kmer_index <- function(s, k) {
  n <- nchar(s)
  if (n < k) return(list(km = character(), pos = integer()))
  st <- seq_len(n - k + 1L)
  km <- substring(s, st, st + k - 1L)
  keep <- grepl("^[ACGT]+$", km) & !(duplicated(km) | duplicated(km, fromLast = TRUE))
  list(km = km[keep], pos = which(keep))
}

# indices of a longest strictly increasing subsequence of x
.lis_idx <- function(x) {
  n <- length(x)
  if (!n) return(integer())
  tv <- numeric(n); ti <- integer(n); prev <- integer(n); len <- 0L
  for (i in seq_len(n)) {
    lo <- 1L; hi <- len + 1L
    while (lo < hi) { mid <- (lo + hi) %/% 2L; if (tv[mid] < x[i]) lo <- mid + 1L else hi <- mid }
    tv[lo] <- x[i]; ti[lo] <- i; prev[i] <- if (lo > 1L) ti[lo - 1L] else 0L
    if (lo > len) len <- lo
  }
  out <- integer(len); j <- ti[len]
  for (p in len:1) { out[p] <- j; j <- prev[j] }
  out
}

# chained anchors, merged into blocks that do not overlap: list(r = start on the reference, q = start on the query, len)
.anchor_blocks <- function(ri, qi, k) {
  common <- intersect(ri$km, qi$km)
  if (!length(common)) return(NULL)
  rp <- ri$pos[match(common, ri$km)]; qp <- qi$pos[match(common, qi$km)]
  o <- order(rp); rp <- rp[o]; qp <- qp[o]
  keep <- .lis_idx(qp); rp <- rp[keep]; qp <- qp[keep]
  R0 <- integer(length(rp)); Q0 <- R0; Ln <- R0; nb <- 0L
  for (i in seq_along(rp)) {
    r <- rp[i]; q <- qp[i]
    if (nb > 0L) {
      pr <- R0[nb]; pq <- Q0[nb]; pl <- Ln[nb]
      if (r - q == pr - pq && r <= pr + pl) { Ln[nb] <- max(pl, r + k - pr); next }   # same diagonal: extend the block
      if (r < pr + pl || q < pq + pl) next                                          # overlaps the block on another diagonal
    }
    nb <- nb + 1L; R0[nb] <- r; Q0[nb] <- q; Ln[nb] <- k
  }
  list(r = R0[seq_len(nb)], q = Q0[seq_len(nb)], len = Ln[seq_len(nb)])
}

# copy an alignment (pa = query row, sa = reference row) into base (one base per reference position) and ins (insertions
# found just before reference position i). refpos = the reference position before the first reference base of the alignment.
.fill_aln <- function(pa, sa, refpos, base, ins) {
  pending <- ""
  for (k in seq_along(sa)) {
    if (sa[k] == "-") pending <- paste0(pending, pa[k])
    else {
      refpos <- refpos + 1L
      if (nzchar(pending)) { ins[refpos] <- paste0(ins[refpos], pending); pending <- "" }
      base[refpos] <- pa[k]
    }
  }
  if (nzchar(pending) && refpos + 1L <= length(ins)) ins[refpos + 1L] <- paste0(ins[refpos + 1L], pending)
  list(base = base, ins = ins)
}

# ---- ViralMSA's own way: minimap2 + the CIGAR of the SAM file ---------------------------------------------------------------------
# Same steps as ViralMSA.py (v1.2.1): minimap2 -a --score-N=0 --secondary=no --sam-hit-only; only primary alignments (flag 0 / 16);
# the row of a sequence is written from the CIGAR (M = bases, D = gaps, I and S = left out); sequences that minimap2 cannot map
# are left out. With keep_ins the 2-pass "insertion strategy" is used: the insertions found in pass 1 are put into an augmented
# reference (leaf words of a trie per position, as in ViralMSA) and everything is mapped again to the augmented reference.

# minimap2 has no Windows build. Where it is found:
#   - a path typed by the user (type  wsl  to force WSL), else a copy inside the app folder, else the PATH (Linux, macOS);
#   - on Windows: minimap2 installed inside WSL (Windows Subsystem for Linux). It is then called as  wsl -e minimap2 ...  and the
#     result is "WSL:minimap2".
.mm2_state <- new.env(parent = emptyenv())

.mm2_wsl_path <- function() {      # path of minimap2 inside WSL, or "" (re-checked at most every 30 s: starting WSL takes a moment)
  if (!is.null(.mm2_state$t) && as.numeric(difftime(Sys.time(), .mm2_state$t, units = "secs")) < 30) return(.mm2_state$path)
  p <- ""
  if (.Platform$OS.type == "windows" && nzchar(Sys.which("wsl"))) {
    out <- tryCatch(suppressWarnings(system2("wsl", c("-e", "sh", "-c", shQuote("command -v minimap2")), stdout = TRUE, stderr = FALSE)),
                    error = function(e) character())
    st <- attr(out, "status")
    if ((is.null(st) || identical(as.integer(st), 0L)) && length(out) && startsWith(out[1], "/")) p <- trimws(out[1])
  }
  .mm2_state$path <- p; .mm2_state$t <- Sys.time()
  p
}

.msa_find_minimap2 <- function(user = "") {
  user <- trimws(user)
  if (tolower(user) == "wsl") return(if (nzchar(.mm2_wsl_path())) "WSL:minimap2" else "")
  if (nzchar(user)) return(if (file.exists(user)) normalizePath(user, winslash = "/") else "")
  cand <- c(file.path("minimap2", "minimap2.exe"), file.path("minimap2", "minimap2"),
            file.path("bin", "minimap2.exe"), file.path("bin", "minimap2"), "minimap2.exe")
  cand <- cand[file.exists(cand)]
  if (length(cand)) return(normalizePath(cand[1], winslash = "/"))
  p <- unname(Sys.which("minimap2"))
  if (nzchar(p)) return(p)
  if (nzchar(.mm2_wsl_path())) return("WSL:minimap2")
  ""
}

# C:/Users/x/file  ->  /mnt/c/Users/x/file  (the path of a Windows file as WSL sees it)
.wsl_path <- function(p) {
  p <- gsub("\\\\", "/", normalizePath(p, winslash = "/", mustWork = FALSE))
  m <- regmatches(p, regexec("^([A-Za-z]):(.*)$", p))[[1]]
  if (length(m) == 3L) paste0("/mnt/", tolower(m[2]), m[3]) else p
}

# runs minimap2 on a reference and named query sequences; returns the lines of the SAM output
.mm2_sam <- function(bin, ref_seq, qseqs) {
  d <- tempfile("mm2_"); dir.create(d); on.exit(unlink(d, recursive = TRUE), add = TRUE)
  rf <- file.path(d, "ref.fa"); qf <- file.path(d, "q.fa"); of <- file.path(d, "out.sam"); lf <- file.path(d, "log.txt")
  writeLines(c(">ref", ref_seq), rf)
  writeLines(as.vector(rbind(paste0(">", names(qseqs)), unname(qseqs))), qf)
  thr <- tryCatch(max(1L, min(8L, parallel::detectCores())), error = function(e) 1L)
  opts <- c("-t", as.character(thr), "--score-N=0", "--secondary=no", "--sam-hit-only", "-a")
  if (startsWith(bin, "WSL:")) {
    cmd <- "wsl"; args <- c("-e", .mm2_wsl_path(), opts, shQuote(.wsl_path(rf)), shQuote(.wsl_path(qf)))
  } else {
    cmd <- bin; args <- c(opts, shQuote(rf), shQuote(qf))
  }
  st <- tryCatch(suppressWarnings(system2(cmd, args, stdout = of, stderr = lf)),
                 error = function(e) stop("minimap2 could not be run: ", conditionMessage(e)))
  if (!identical(as.integer(st), 0L) || !file.exists(of))
    stop("minimap2 failed. Its message: ", paste(utils::tail(if (file.exists(lf)) readLines(lf, warn = FALSE) else "", 5), collapse = " "))
  readLines(of, warn = FALSE)
}

# the primary alignments of a SAM file: id, flag, pos (1-based), cigar, seq
.sam_records <- function(sam) {
  sam <- sam[nzchar(sam) & !startsWith(sam, "@")]
  empty <- data.frame(id = character(), flag = integer(), pos = integer(), cigar = character(), seq = character(), stringsAsFactors = FALSE)
  if (!length(sam)) return(empty)
  f <- strsplit(sam, "\t", fixed = TRUE)
  f <- f[vapply(f, length, 1L) >= 11L]
  if (!length(f)) return(empty)
  g <- function(i) vapply(f, `[`, character(1), i)
  d <- data.frame(id = g(1L), flag = suppressWarnings(as.integer(g(2L))), pos = suppressWarnings(as.integer(g(4L))),
                  cigar = g(6L), seq = toupper(g(10L)), stringsAsFactors = FALSE)
  d[d$flag %in% c(0L, 16L) & d$cigar != "*" & d$seq != "*" & !is.na(d$pos) & d$pos >= 1L, , drop = FALSE]
}

.cigar_ops <- function(cigar) {
  ops <- regmatches(cigar, gregexpr("[0-9]+[MIDNSHP=X]", cigar))[[1]]
  list(n = as.integer(substr(ops, 1L, nchar(ops) - 1L)), op = substr(ops, nchar(ops), nchar(ops)))
}

# one row of length L per primary alignment (ViralMSA: gaps before, M/=/X bases, D gaps, I and S left out, gaps after)
.sam_rows <- function(sam, L) {
  d <- .sam_records(sam)
  rows <- character(nrow(d))
  for (i in seq_len(nrow(d))) {
    co <- .cigar_ops(d$cigar[i]); parts <- character(length(co$op)); qi <- 1L
    for (j in seq_along(co$op)) {
      o <- co$op[j]; n <- co$n[j]
      if (o %in% c("M", "=", "X")) { parts[j] <- substr(d$seq[i], qi, qi + n - 1L); qi <- qi + n }
      else if (o %in% c("D", "N")) parts[j] <- strrep("-", n)
      else if (o %in% c("I", "S")) qi <- qi + n
    }
    row <- paste0(strrep("-", d$pos[i] - 1L), paste(parts, collapse = ""))
    if (nchar(row) < L) row <- paste0(row, strrep("-", L - nchar(row)))
    rows[i] <- row
  }
  list(rows = stats::setNames(rows, d$id), rev = d$id[d$flag == 16L])
}

# insertions[[k]] = the insertion strings found just before reference base k (k = L + 1: after the last base)
.sam_insertions <- function(sam, L) {
  d <- .sam_records(sam); ins <- vector("list", L + 1L)
  for (i in seq_len(nrow(d))) {
    co <- .cigar_ops(d$cigar[i]); rp <- d$pos[i] - 1L; qi <- 1L
    for (j in seq_along(co$op)) {
      o <- co$op[j]; n <- co$n[j]
      if (o %in% c("M", "=", "X")) { rp <- rp + n; qi <- qi + n }
      else if (o == "I") {
        k <- rp + 1L
        if (k <= L + 1L) ins[[k]] <- c(ins[[k]], substr(d$seq[i], qi, qi + n - 1L))
        qi <- qi + n
      }
      else if (o %in% c("D", "N")) rp <- rp + n
      else if (o == "S") qi <- qi + n
    }
  }
  ins
}

# the augmented reference: before every base, the insertion words that are not the beginning of another word, in sorted order
.vmsa_augment <- function(ref, ins) {
  L <- nchar(ref); ref_ch <- strsplit(ref, "")[[1]]; pre <- character(L + 1L)
  for (k in which(!vapply(ins, is.null, logical(1)))) {
    w <- sort(unique(ins[[k]]), method = "radix")
    leaf <- vapply(w, function(x) !any(startsWith(w, x) & w != x), logical(1))
    pre[k] <- paste(w[leaf], collapse = "")
  }
  paste0(paste0(pre[seq_len(L)], ref_ch, collapse = ""), pre[L + 1L])
}

.vmsa_minimap2 <- function(q, ref_id, ref, keep_ins, bin, progress = NULL) {
  tick <- function(i) if (!is.null(progress)) progress(i, 4L)
  ids <- paste0("q", seq_along(q))
  qs <- stats::setNames(vapply(q, function(s) toupper(gsub("[^A-Z]", "", s)), ""), ids)
  qs <- qs[nzchar(qs)]
  ref_use <- ref; ref_row <- NULL
  tick(1L)
  if (keep_ins) {                                         # pass 1: find the insertions
    ins <- .sam_insertions(.mm2_sam(bin, ref, qs), nchar(ref))
    ref_use <- .vmsa_augment(ref, ins)
    qs <- c(stats::setNames(ref, "q0"), qs)               # the original reference is mapped to the augmented one too
    tick(2L)
  }
  sam <- .mm2_sam(bin, ref_use, qs)                       # (pass 2 when keep_ins)
  tick(3L)
  r <- .sam_rows(sam, nchar(ref_use))
  got <- r$rows
  if (!length(got) || !any(names(got) %in% ids)) stop("minimap2 could not map any sequence to the reference.")
  ok <- ids[ids %in% names(got)]
  rows <- c(if (keep_ins) stats::setNames(ref_use, paste0("AUG_REF_", ref_id)) else stats::setNames(ref, ref_id),
            if (keep_ins && "q0" %in% names(got)) stats::setNames(got[["q0"]], paste0("ORIG_REF_", ref_id)),
            stats::setNames(unname(got[ok]), names(q)[match(ok, ids)]))
  tick(4L)
  list(aln = rows,
       flipped = names(q)[match(intersect(ids, r$rev), ids)],
       reference = if (keep_ins && "q0" %in% names(got)) paste0("ORIG_REF_", ref_id) else if (keep_ins) paste0("AUG_REF_", ref_id) else ref_id,
       skipped = names(q)[!(ids %in% names(got))],
       engine = if (keep_ins) "minimap2, 2-pass insertion strategy (as ViralMSA)" else "minimap2 (as ViralMSA)")
}

msa_refalign <- function(seqs, ref_id, keep_ins = FALSE, progress = NULL, minimap2 = "") {
  if (!(ref_id %in% names(seqs))) stop("The reference is not among the sequences.")
  ref <- toupper(gsub("[^A-Z]", "", seqs[[ref_id]])); L <- nchar(ref)
  q <- seqs[names(seqs) != ref_id]
  if (!length(q)) stop("Add at least one sequence besides the reference.")
  if (L < 20) stop("The reference is too short.")
  if (nzchar(minimap2)) return(.vmsa_minimap2(q, ref_id, ref, keep_ins, minimap2, progress))   # ViralMSA's own way
  mat <- .pw("nucleotideSubstitutionMatrix")(match = 2, mismatch = -3, baseOnly = FALSE)
  pa_fn <- .pw("pairwiseAlignment"); sc_fn <- .pw("score")
  ap_fn <- .pw("alignedPattern"); as_fn <- .pw("alignedSubject"); su_fn <- .pw("subject")

  max_direct <- 6e7      # sequence length x reference length: up to this the exact aligner is used for the whole sequence
  W <- 3000L             # anchored method: window aligned exactly at each end of the sequence
  ref_ch <- strsplit(ref, "")[[1]]
  ri_cache <- list()
  ref_idx <- function(k) {
    key <- as.character(k)
    if (is.null(ri_cache[[key]])) ri_cache[[key]] <<- .kmer_index(ref, k)
    ri_cache[[key]]
  }
  aln_chars <- function(a) list(pa = strsplit(as.character(ap_fn(a)), "")[[1]], sa = strsplit(as.character(as_fn(a)), "")[[1]])
  pair_aln <- function(Qs, Rs, type) pa_fn(pattern = Qs, subject = Rs, type = type, substitutionMatrix = mat, gapOpening = 8, gapExtension = 4)

  # exact: the whole sequence against the whole reference (both strands tried)
  one_direct <- function(s) {
    best <- NULL
    for (rc in c(FALSE, TRUE)) {
      qs <- if (rc) revcomp_chr(s) else s
      a <- pair_aln(qs, ref, "overlap")
      sc <- as.numeric(sc_fn(a))
      if (is.null(best) || sc > best$score) best <- list(score = sc, a = a, rc = rc)
    }
    ac <- aln_chars(best$a)
    s0 <- as.integer(BiocGenerics::start(su_fn(best$a)))
    f <- .fill_aln(ac$pa, ac$sa, s0 - 1L, rep("-", L), character(L + 1L))
    list(base = f$base, ins = f$ins, rc = best$rc)
  }

  # anchored: for sequences too long for the exact method (viral genomes)
  one_anchored <- function(s, nm) {
    best <- NULL
    for (k in c(15L, 12L)) {                       # 12-mers only if 15-mers find too little (divergent sequences)
      ri <- ref_idx(k)
      for (rc in c(FALSE, TRUE)) {
        qs <- if (rc) revcomp_chr(s) else s
        bl <- .anchor_blocks(ri, .kmer_index(qs, k), k)
        sc <- if (is.null(bl)) 0 else sum(bl$len)
        if (is.null(best) || sc > best$score) best <- list(score = sc, b = bl, rc = rc, qs = qs)
      }
      if (!is.null(best$b) && length(best$b$len) >= 3L) break
    }
    if (is.null(best$b) || length(best$b$len) < 3L)
      stop(sprintf("No similar region was found between the reference and %s (too different, or too short).",
                   if (nzchar(nm)) paste0("'", nm, "'") else "a sequence"))
    b <- best$b; nbk <- length(b$len)
    qc <- strsplit(best$qs, "")[[1]]
    base <- rep("-", L); ins <- character(L + 1L)
    for (i in seq_len(nbk)) {
      r <- b$r[i]; q <- b$q[i]; len <- b$len[i]
      base[r:(r + len - 1L)] <- qc[q:(q + len - 1L)]
      if (i > 1L) {                                # the stretch between block i-1 and block i
        gr0 <- b$r[i - 1L] + b$len[i - 1L]; gq0 <- b$q[i - 1L] + b$len[i - 1L]
        dr <- r - gr0; dq <- q - gq0
        if (dr > 0L && dq > 0L) {
          Rg <- ref_ch[gr0:(r - 1L)]; Qg <- qc[gq0:(q - 1L)]
          if (dr == dq && sum(Rg != Qg) <= 7L) base[gr0:(r - 1L)] <- Qg      # substitutions only (cheaper than two gaps)
          else if (as.numeric(dr) * as.numeric(dq) <= 4e7) {
            a <- pair_aln(paste(Qg, collapse = ""), paste(Rg, collapse = ""), "global")
            ac <- aln_chars(a)
            f <- .fill_aln(ac$pa, ac$sa, gr0 - 1L, base, ins); base <- f$base; ins <- f$ins
          } else ins[r] <- paste0(ins[r], paste(Qg, collapse = ""))          # too big to align: kept as an insertion
        } else if (dq > 0L) ins[r] <- paste0(ins[r], paste(qc[gq0:(q - 1L)], collapse = ""))   # pure insertion
        # dr > 0 and dq == 0: a deletion; the reference positions stay "-"
      }
    }
    # the two ends: exact alignment (free end gaps) of a window next to the first / last block
    end_aln <- function(Rs, Qs, ref_off) {
      a <- tryCatch(pair_aln(Qs, Rs, "overlap"), error = function(e) NULL)
      if (is.null(a)) return(NULL)
      ac <- aln_chars(a)
      if (!length(ac$sa)) return(NULL)
      f <- .fill_aln(ac$pa, ac$sa, ref_off + as.integer(BiocGenerics::start(su_fn(a))) - 1L, base, ins)
      f
    }
    r1 <- b$r[1]; q1 <- b$q[1]
    if (r1 > 1L && q1 > 1L) {
      ra <- max(1L, r1 - W); qa <- max(1L, q1 - W)
      f <- end_aln(paste(ref_ch[ra:(r1 - 1L)], collapse = ""), paste(qc[qa:(q1 - 1L)], collapse = ""), ra - 1L)
      if (!is.null(f)) { base <- f$base; ins <- f$ins }
    }
    r2 <- b$r[nbk] + b$len[nbk]; q2 <- b$q[nbk] + b$len[nbk]
    if (r2 <= L && q2 <= length(qc)) {
      rb <- min(L, r2 + W - 1L); qb <- min(length(qc), q2 + W - 1L)
      f <- end_aln(paste(ref_ch[r2:rb], collapse = ""), paste(qc[q2:qb], collapse = ""), r2 - 1L)
      if (!is.null(f)) { base <- f$base; ins <- f$ins }
    }
    list(base = base, ins = ins, rc = best$rc)
  }

  one <- function(s, nm = "") {
    s <- toupper(gsub("[^A-Z]", "", s))
    if (as.numeric(nchar(s)) * as.numeric(L) <= max_direct) one_direct(s) else one_anchored(s, nm)
  }

  res <- vector("list", length(q)); err <- character(length(q))
  for (i in seq_along(q)) {
    if (!is.null(progress)) progress(i, length(q))
    res[[i]] <- tryCatch(one(q[[i]], names(q)[i]), error = function(e) { err[i] <<- conditionMessage(e); NULL })
  }
  bad <- vapply(res, is.null, logical(1))           # like ViralMSA: a sequence that cannot be placed is left out (and reported)
  if (all(bad)) stop("No sequence could be placed on the reference. ", err[1])
  skipped <- names(q)[bad]; q <- q[!bad]; res <- res[!bad]
  flipped <- names(q)[vapply(res, function(r) isTRUE(r$rc), logical(1))]
  # insertion columns (only when asked): the longest insertion at every reference position
  mx <- if (keep_ins) Reduce(pmax, lapply(res, function(r) nchar(r$ins)), accumulate = FALSE) else rep(0L, L + 1L)
  make_row <- function(base, ins) {
    if (!keep_ins || !any(mx[seq_len(L)] > 0)) return(paste(base, collapse = ""))
    ex <- character(L); ix <- which(mx[seq_len(L)] > 0L)
    ex[ix] <- paste0(ins[ix], strrep("-", mx[ix] - nchar(ins[ix])))
    paste0(ex, base, collapse = "")
  }
  ref_ch <- strsplit(ref, "")[[1]]
  rows <- c(make_row(ref_ch, character(L + 1L)),
            vapply(res, function(r) make_row(r$base, r$ins), ""))
  names(rows) <- c(ref_id, names(q))
  list(aln = rows, flipped = flipped, reference = ref_id, skipped = skipped, engine = "built-in R aligner")
}

# ---- mutations -------------------------------------------------------------------------------------------------------
.aln_variants <- function(aln, ref_name, refpos = NULL) {
  empty <- data.frame(Sequence = character(), Type = character(), Position = integer(), End = integer(),
                      Length = integer(), Ref = character(), Alt = character(), Notation = character(),
                      Column = integer(), stringsAsFactors = FALSE)
  if (is.null(aln) || !(ref_name %in% names(aln))) return(empty)
  M <- do.call(rbind, strsplit(toupper(unname(aln)), "")); rownames(M) <- names(aln)
  r <- M[ref_name, ]; isb <- r != "-"
  rpos <- if (is.null(refpos)) cumsum(isb) else refpos     # refpos: reference position of every column (trimmed alignments keep the original numbering)
  acgt <- c("A", "C", "G", "T")
  runs <- function(ix) {                                    # consecutive column numbers -> first and last column of every run
    if (!length(ix)) return(list(a = integer(), b = integer()))
    brk <- c(TRUE, diff(ix) != 1L)
    list(a = ix[brk], b = ix[c(brk[-1L], TRUE)])
  }
  others <- setdiff(rownames(M), ref_name)
  out <- vector("list", length(others)); n <- 0L; col <- seq_len(ncol(M))
  for (nm in others) {
    s <- M[nm, ]; ng <- which(s != "-"); if (!length(ng)) next
    inr <- col >= min(ng) & col <= max(ng)                  # only between the first and the last base of the sequence
    sub <- which(inr & isb & s != "-" & s != r & s != "N" & r != "N")
    del <- runs(which(inr & isb & s == "-"))
    ins <- runs(which(inr & !isb & s != "-"))
    parts <- list()
    if (length(sub)) parts[[length(parts) + 1L]] <- data.frame(
      Sequence = nm, Type = ifelse(s[sub] %in% acgt, "Substitution", "Ambiguous"), Position = rpos[sub], End = rpos[sub], Length = 1L,
      Ref = r[sub], Alt = s[sub], Notation = sprintf("%d%s>%s", rpos[sub], r[sub], s[sub]), Column = sub, stringsAsFactors = FALSE)
    if (length(del$a)) {
      len <- del$b - del$a + 1L
      parts[[length(parts) + 1L]] <- data.frame(
        Sequence = nm, Type = "Deletion", Position = rpos[del$a], End = rpos[del$b], Length = len,
        Ref = vapply(seq_along(len), function(i) paste(r[del$a[i]:del$b[i]], collapse = ""), ""), Alt = "-",
        Notation = ifelse(len == 1L, sprintf("%ddel", rpos[del$a]), sprintf("%d_%ddel", rpos[del$a], rpos[del$b])),
        Column = del$a, stringsAsFactors = FALSE)
    }
    if (length(ins$a)) {
      len <- ins$b - ins$a + 1L
      alt <- vapply(seq_along(len), function(i) paste(s[ins$a[i]:ins$b[i]], collapse = ""), "")
      parts[[length(parts) + 1L]] <- data.frame(
        Sequence = nm, Type = "Insertion", Position = rpos[ins$a], End = rpos[ins$a] + 1L, Length = len, Ref = "-", Alt = alt,
        Notation = sprintf("%d_%dins%s", rpos[ins$a], rpos[ins$a] + 1L, alt), Column = ins$a, stringsAsFactors = FALSE)
    }
    if (length(parts)) { n <- n + 1L; out[[n]] <- do.call(rbind, parts) }
  }
  if (!n) return(empty)
  d <- do.call(rbind, out[seq_len(n)])
  d[order(d$Position, d$Sequence), , drop = FALSE]
}

# which columns of the whole alignment are in a trimmed copy (same columns in the same order, some removed)? NULL when that cannot be told
.match_columns <- function(aln, tr) {
  nm <- intersect(names(tr), names(aln))
  if (length(nm) < 2L) return(NULL)
  A <- do.call(rbind, strsplit(toupper(unname(aln[nm])), "")); B <- do.call(rbind, strsplit(toupper(unname(tr[nm])), ""))
  ca <- apply(A, 2L, paste, collapse = ""); cb <- apply(B, 2L, paste, collapse = "")
  keep <- logical(length(ca)); j <- 1L; nb <- length(cb)
  for (i in seq_along(ca)) if (j <= nb && ca[i] == cb[j]) { keep[i] <- TRUE; j <- j + 1L }
  if (j <= nb) NULL else keep
}

.aln_var_summary <- function(aln, ref_name, v, refpos = NULL) {
  M <- do.call(rbind, strsplit(toupper(unname(aln)), "")); rownames(M) <- names(aln)
  r <- M[ref_name, ]; isb <- r != "-"; rpos <- if (is.null(refpos)) cumsum(isb) else refpos
  rows <- lapply(setdiff(rownames(M), ref_name), function(nm) {
    s <- M[nm, ]; ng <- which(s != "-")
    if (!length(ng)) return(data.frame(Sequence = nm, Covers_from = NA, Covers_to = NA, Covered_positions = 0L, Identity_pct = NA_real_,
                                        Substitutions = 0L, Ambiguous = 0L, Insertions = 0L, Deletions = 0L, Inserted_bases = 0L,
                                        Deleted_bases = 0L, N_bases = 0L, stringsAsFactors = FALSE))
    a <- min(ng); b <- max(ng); cols <- a:b; rc <- cols[isb[cols]]
    cov <- length(rc); same <- sum(s[rc] == r[rc] & s[rc] != "-" & r[rc] != "N")
    x <- v[v$Sequence == nm, , drop = FALSE]
    data.frame(Sequence = nm, Covers_from = rpos[a] + as.integer(!isb[a]), Covers_to = rpos[b], Covered_positions = cov,
               Identity_pct = if (cov) round(100 * same / cov, 2) else NA_real_,
               Substitutions = sum(x$Type == "Substitution"), Ambiguous = sum(x$Type == "Ambiguous"),
               Insertions = sum(x$Type == "Insertion"), Deletions = sum(x$Type == "Deletion"),
               Inserted_bases = sum(x$Length[x$Type == "Insertion"]), Deleted_bases = sum(x$Length[x$Type == "Deletion"]),
               N_bases = sum(s[cols] == "N"), stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

# which sequences carry which mutation (1 / 0); the same mutation in several sequences is one row
.aln_var_matrix <- function(v, seq_names) {
  if (!nrow(v)) return(data.frame(Mutation = character(), Type = character(), Position = integer(), Count = integer(), stringsAsFactors = FALSE))
  key <- paste(v$Type, v$Notation)
  u <- v[!duplicated(key), c("Type", "Notation", "Position"), drop = FALSE]
  m <- vapply(seq_along(seq_names), function(j) as.integer(paste(u$Type, u$Notation) %in% key[v$Sequence == seq_names[j]]), integer(nrow(u)))
  if (is.null(dim(m))) m <- matrix(m, nrow = nrow(u))
  colnames(m) <- seq_names
  d <- data.frame(Mutation = u$Notation, Type = u$Type, Position = u$Position, Count = rowSums(m), m, check.names = FALSE, stringsAsFactors = FALSE)
  d[order(d$Position, d$Type), , drop = FALSE]
}

.var_colours <- c(Substitution = "#d62728", Ambiguous = "#ff7f0e", Insertion = "#2ca02c", Deletion = "#222222")

.aln_var_svg <- function(v, seq_names, L, width = 1000L, max_rows = 400L) {
  esc <- htmltools::htmlEscape
  rows <- utils::head(seq_names, max_rows); h <- 16L; top <- 24L; lab <- 190L; W <- width + lab + 20L
  H <- top + h * length(rows) + 22L; sc <- width / max(L, 1)
  step <- 10^floor(log10(max(L, 10) / 5)); step <- step * c(1, 2, 5, 10)[which(L / (step * c(1, 2, 5, 10)) <= 12)[1]]
  ticks <- seq(0, L, by = step)
  parts <- c(sprintf('<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d" font-family="Segoe UI,Arial,sans-serif" font-size="11">', W, H),
             sprintf('<text x="%d" y="14" font-size="11" fill="#555">position on the reference</text>', lab),
             sprintf('<line x1="%d" y1="%d" x2="%d" y2="%d" stroke="#999"/>', lab, top - 4L, lab + width, top - 4L),
             sprintf('<line x1="%.1f" y1="%d" x2="%.1f" y2="%d" stroke="#999"/><text x="%.1f" y="%d" fill="#555" text-anchor="middle">%d</text>',
                     lab + ticks * sc, top - 8L, lab + ticks * sc, top - 4L, lab + ticks * sc, top - 10L, as.integer(ticks)))
  if (length(rows)) {
    y <- top + (seq_along(rows) - 1L) * h
    short <- .short_labels(rows); short <- ifelse(nchar(short) > 30L, paste0(substr(short, 1L, 28L), ".."), short)
    parts <- c(parts, sprintf('<text x="%d" y="%d" text-anchor="end" fill="#222"><title>%s</title>%s</text>', lab - 6L, y + 11L, esc(rows), esc(short)),
               sprintf('<rect x="%d" y="%d" width="%d" height="%d" fill="#eef1f5"/>', lab, y + 2L, width, h - 4L))
    ri <- match(v$Sequence, rows); ok <- !is.na(ri); x <- v[ok, , drop = FALSE]; ri <- ri[ok]
    if (nrow(x)) {
      key <- paste(ri, x$Type, floor((x$Position - 1) * sc / 1.5)); first <- !duplicated(key)
      cnt <- tabulate(match(key, key[first])); x <- x[first, , drop = FALSE]; ri <- ri[first]
      tip <- paste0(x$Type, " ", x$Notation, ifelse(cnt > 1L, sprintf(" (+%d nearby)", cnt - 1L), ""))
      parts <- c(parts, sprintf('<rect x="%.2f" y="%d" width="%.2f" height="%d" fill="%s"><title>%s</title></rect>',
                                lab + (x$Position - 1) * sc, top + (ri - 1L) * h + 1L, pmax(1.6, (x$End - x$Position + 1) * sc), h - 2L,
                                unname(.var_colours[x$Type]), esc(tip)))
    }
  }
  paste(c(parts, "</svg>"), collapse = "\n")
}

# ---- zoomable map of the mutations (used in the Mutations tab) -----------------------------------------------------------------------
# Names like "MW308135.1_Yellow_fever_virus_isolate_X__complete_genome": the words that most names share are left out of the
# label (the full name stays in the tooltip), so the part that tells the sequences apart is visible.
.short_labels <- function(nm) {
  if (length(nm) < 3L) return(nm)
  tk <- strsplit(nm, "_", fixed = TRUE)
  freq <- table(unlist(lapply(tk, unique)))
  common <- names(freq)[freq >= 0.7 * length(nm)]
  vapply(tk, function(x) { y <- x[!(x %in% common)]; if (length(y)) paste(y, collapse = "_") else paste(x, collapse = "_") }, "")
}


.aln_var_html <- function(v, summ, ref_name, ref_len, seq_names, title = "Mutations relative to the reference") {
  esc <- htmltools::htmlEscape
  tbl <- function(d) {
    if (!nrow(d)) return("<p><i>none</i></p>")
    th <- paste0("<tr>", paste0("<th>", esc(names(d)), "</th>", collapse = ""), "</tr>")
    td <- apply(d, 1, function(r) paste0("<tr>", paste0("<td>", esc(ifelse(is.na(r), "", as.character(r))), "</td>", collapse = ""), "</tr>"))
    paste0("<table>", th, paste(td, collapse = "\n"), "</table>")
  }
  vc <- v; if (nrow(vc)) vc$Type <- sprintf('<span class="t t%s">%s</span>', substr(vc$Type, 1, 3), vc$Type)
  cap <- 20000L
  vtbl <- if (!nrow(v)) "<p><i>No mutation found.</i></p>" else {
    th <- paste0("<tr>", paste0("<th>", esc(c("Sequence", "Type", "Notation", "Position", "End", "Length", "Ref", "Alt")), "</th>", collapse = ""), "</tr>")
    vv <- utils::head(v, cap)
    td <- sprintf("<tr><td>%s</td><td><span class='t t%s'>%s</span></td><td>%s</td><td>%d</td><td>%d</td><td>%d</td><td>%s</td><td>%s</td></tr>",
                  esc(vv$Sequence), substr(vv$Type, 1, 3), esc(vv$Type), esc(vv$Notation), as.integer(vv$Position), as.integer(vv$End),
                  as.integer(vv$Length), esc(vv$Ref), esc(vv$Alt))
    paste0(if (nrow(v) > cap) sprintf("<p>Only the first %d of %d mutations are listed here; the complete list is in the CSV file.</p>", cap, nrow(v)) else "",
           "<table>", th, paste(td, collapse = "\n"), "</table>")
  }
  tot <- table(factor(v$Type, levels = names(.var_colours)))
  paste0('<!DOCTYPE html><html><head><meta charset="utf-8"><title>', esc(title), '</title><style>',
    'body{font-family:"Segoe UI",Arial,sans-serif;margin:24px;color:#1b2a3a}h1{font-size:22px}h2{font-size:17px;margin-top:28px;color:#1c4c8c}',
    'table{border-collapse:collapse;font-size:13px;margin:8px 0}th,td{border:1px solid #d5dbe3;padding:4px 8px;text-align:left}',
    'th{background:#eef3fa}tr:nth-child(even) td{background:#fafbfd}.t{padding:1px 6px;border-radius:8px;color:#fff;font-size:12px}',
    '.tSub{background:#d62728}.tAmb{background:#ff7f0e}.tIns{background:#2ca02c}.tDel{background:#222}.box{background:#f4f7fb;border:1px solid #d5dbe3;border-radius:8px;padding:10px 14px;display:inline-block;margin-right:10px}',
    '.leg span{margin-right:14px}</style></head><body>',
    '<h1>', esc(title), '</h1>',
    '<p>Reference: <b>', esc(ref_name), '</b> (', ref_len, ' positions). Sequences compared: ', length(seq_names) - 1L, '. Created by SangeRly, ',
    format(Sys.time(), "%Y-%m-%d %H:%M"), '.</p>',
    '<div><span class="box"><b>', tot[["Substitution"]], '</b> substitutions</span><span class="box"><b>', tot[["Insertion"]],
    '</b> insertions</span><span class="box"><b>', tot[["Deletion"]], '</b> deletions</span><span class="box"><b>', tot[["Ambiguous"]], '</b> ambiguous bases</span></div>',
    '<h2>Map of the mutations</h2><p class="leg"><span style="color:#d62728">&#9632; substitution</span><span style="color:#ff7f0e">&#9632; ambiguous (IUPAC)</span>',
    '<span style="color:#2ca02c">&#9632; insertion</span><span style="color:#222">&#9632; deletion</span></p>',
    .aln_var_svg(v, setdiff(seq_names, ref_name), ref_len),
    '<h2>Per sequence</h2>', tbl(summ),
    '<h2>All mutations</h2>', vtbl,
    '<p style="color:#6c757d;font-size:12px">Notation: 123A&gt;G substitution at reference position 123; 150_152del deletion; 200_201insACT insertion between positions 200 and 201. ',
    'The ends of a sequence that does not reach the ends of the reference are not counted as deletions.</p>',
    '</body></html>')
}
