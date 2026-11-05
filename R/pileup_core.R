# ==============================================================================
# Sangerly - Multi-read pileup (v4)
# Mode 1: Forward + Reverse pair, no reference.
# Mode 2: many reads / tiled regions mapped to a reference.
#
# aln object:
#   ok, msg, mode ("denovo"/"ref"), notes, P
#   reads   : list of oriented reads (seq, qual, ploc, idx, traces, rc, raw, qc, pk)
#   labels, strand ("Forward"/"Reverse")
#   ref     : NULL or list(name, seq); refbase; refpos (reference position per column)
#   base    : character matrix [read x column]; NA = not covered, "-" = deletion
#   idx     : integer matrix, base position within the oriented read
#   q, ratio, top, sec : Phred, secondary/primary peak ratio, primary/secondary peak letter
#   cols    : automatic per-column decisions (see aln_decide)
#   insertions, map_info
#
# Deletion model:
#   mask   : logical matrix [read x column]; TRUE = that read's base is deleted.
#   delcol : logical vector [column]; TRUE = the whole position is deleted, i.e. the
#            consensus base is removed too (not replaced by N) and no read contributes.
# Deleted bases are ignored by decisions, consensus, SNPs and genetic distances.
#
# Needs assembly_core.R (read_ab1, prep_read, rc_read, pair_alignment, .iupac, .pw,
# final_cols, .nz) to be sourced first.
# ==============================================================================

# Primary/secondary peak for EVERY base of one read (computed once) ----------------
.peak_table <- function(rd) {
  n <- length(rd$ploc)
  top <- character(n); sec <- character(n); ratio <- numeric(n)
  tr <- rd$traces; nr <- nrow(tr)
  for (i in seq_len(n)) {
    sc <- min(max(rd$ploc[i], 1L), nr)
    rng <- max(1, sc - 1):min(nr, sc + 1)
    v <- apply(tr[rng, , drop = FALSE], 2, max)
    o <- order(v, decreasing = TRUE)
    top[i] <- names(v)[o[1]]; sec[i] <- names(v)[o[2]]
    ratio[i] <- if (v[o[1]] > 0) unname(v[o[2]] / v[o[1]]) else NA_real_
  }
  data.frame(top = top, second = sec, ratio = ratio, stringsAsFactors = FALSE)
}

# Supporting matrices (Phred, ratio, peaks) from the base and idx matrices ----------
aln_matrices <- function(reads, base, idx) {
  nr <- nrow(base); nc <- ncol(base)
  q <- ratio <- matrix(NA_real_, nr, nc)
  top <- sec <- matrix(NA_character_, nr, nc)
  for (i in seq_len(nr)) {
    reads[[i]]$pk <- .peak_table(reads[[i]])
    ok <- which(!is.na(idx[i, ]))
    ii <- idx[i, ok]
    q[i, ok] <- reads[[i]]$qual[ii]
    ratio[i, ok] <- reads[[i]]$pk$ratio[ii]
    top[i, ok] <- reads[[i]]$pk$top[ii]
    sec[i, ok] <- reads[[i]]$pk$second[ii]
    # Phred of a deletion cell = minimum of the neighbouring bases of the same read
    dl <- which(!is.na(base[i, ]) & base[i, ] == "-")
    for (k in dl) {
      cand <- c(if (any(ok < k)) q[i, max(ok[ok < k])],
                if (any(ok > k)) q[i, min(ok[ok > k])])
      q[i, k] <- if (length(cand)) min(cand) else NA_real_
    }
  }
  list(reads = reads, base = base, idx = idx, q = q, ratio = ratio, top = top, sec = sec)
}

# ------------------------------------------------------------------------------
# Decision for one column covered by >= 1 read
# ------------------------------------------------------------------------------
decide_multi <- function(nm, b, q, rt, sec, top, P) {
  mk <- function(type, base, rec, reason, review)
    list(type = type, base = base, rec = rec, reason = reason, review = review)
  q[is.na(q)] <- 0
  q[b == "N"] <- 0
  ev <- paste(sprintf("%s=%s(Q%d)", nm, b, as.integer(round(q))), collapse = "; ")
  n <- length(b)

  if (n == 1) {
    if (b == "N" || q < P$q_poor)
      return(mk("S-low", "N", "Mark as N (single read, low quality)", paste("Single-read coverage;", ev), TRUE))
    if (b == "-")
      return(mk("Del", if (q >= P$q_good) "-" else "N",
                "Deletion supported by one read only: verify", paste("Single-read coverage;", ev), TRUE))
    return(mk("S", b, "Only one read covers this position", paste("Single-read coverage;", ev), FALSE))
  }

  good <- q >= P$q_poor & b != "N"
  if (!any(good))
    return(mk("5", "N", "Mark as N, delete the position, or resequence", paste("All reads are low quality;", ev), TRUE))

  # gap in some reads vs base in others (possible indel / homopolymer)
  if (any(good & b == "-") && any(good & b != "-")) {
    w <- tapply(q[good], b[good], sum)
    best <- names(w)[which.max(w)]
    return(mk("Indel", best, "Gap in some reads vs base in others: check for an indel / homopolymer",
              paste("Phred support:", paste(sprintf("%s=%d", names(w), as.integer(w)), collapse = " "), "|", ev), TRUE))
  }

  # consistent double peak in >= 2 good reads
  hi <- good & !is.na(rt) & rt >= P$ratio & b != "-"
  if (sum(hi) >= 2) {
    s <- unique(c(top[hi], sec[hi]))
    base <- if (length(s) == 2) .iupac(s) else "N"
    return(mk("6", base, "Consistent double peak across reads: consider heterozygosity / mixed template",
              paste0(ev, "; highest secondary ratio ", sprintf("%.2f", max(rt[hi], na.rm = TRUE))), TRUE))
  }

  w <- tapply(q[good], b[good], sum)
  if (length(w) == 1) {
    a <- names(w)
    if (all(b == a))
      return(mk(if (a == "-") "Del" else "1", a,
                if (a == "-") "Deletion supported by all reads: verify" else "Accept", ev,
                a == "-"))
    return(mk("2", a, "Accept the base from the good-quality read(s) (poor read differs)", ev,
              max(q[good]) < P$q_good))
  }

  ord <- sort(w, decreasing = TRUE)
  best <- names(ord)[1]
  if (ord[1] - ord[2] >= P$dq_pref)
    return(mk("2", best, "Base with the largest Phred support; verify visually",
              paste("Phred support:", paste(sprintf("%s=%d", names(ord), as.integer(ord)), collapse = " "), "|", ev),
              TRUE))
  base <- switch(P$type4, iupac = .iupac(names(ord)[1:2]), F = best, "N")
  mk("4", base, "Good reads disagree: verify manually on the chromatogram",
     paste("Balanced Phred support:", paste(sprintf("%s=%d", names(ord), as.integer(ord)), collapse = " "), "|", ev),
     TRUE)
}

# Mask including whole deleted positions
aln_mask_all <- function(mask, delcol = NULL) {
  if (!is.null(delcol) && any(delcol)) mask[, delcol] <- TRUE
  mask
}

# Effective bases after deleted cells are removed
aln_eff_base <- function(aln, mask = NULL) {
  B <- aln$base
  if (!is.null(mask) && any(mask)) B[mask] <- NA_character_
  B
}

# Automatic decisions for all columns; recomputed whenever deletions change
aln_decide <- function(aln, mask = NULL, delcol = NULL, P = aln$P) {
  B <- aln$base; Q <- aln$q; R <- aln$ratio; S <- aln$sec; T <- aln$top
  if (!is.null(mask) && any(mask)) {
    B[mask] <- NA_character_; Q[mask] <- NA_real_; R[mask] <- NA_real_
    S[mask] <- NA_character_; T[mask] <- NA_character_
  }
  nc <- ncol(B)
  if (is.null(delcol)) delcol <- rep(FALSE, nc)
  dec <- vector("list", nc)
  for (k in seq_len(nc)) {
    cv <- which(!is.na(B[, k]))
    dec[[k]] <- if (delcol[k]) {
      list(type = "ColDel", base = "-", rec = "Position deleted by user (removed from the consensus)",
           reason = "Deleted", review = FALSE)
    } else if (!length(cv)) {
      list(type = "NoCov", base = "N",
           rec = "No read left at this position (coverage gap or all bases deleted): fill with N, delete the position, or resequence",
           reason = "No coverage", review = TRUE)
    } else {
      decide_multi(aln$labels[cv], B[cv, k], Q[cv, k], R[cv, k], S[cv, k], T[cv, k], P)
    }
  }
  g <- function(f, type) vapply(dec, function(x) x[[f]], type)
  data.frame(Col = seq_len(nc), RefPos = aln$refpos,
             RefBase = if (is.null(aln$refbase)) rep(NA_character_, nc) else aln$refbase,
             Cov = as.integer(colSums(!is.na(B))), Type = g("type", ""), Auto = g("base", ""),
             Recommendation = g("rec", ""), Reason = g("reason", ""), Review = g("review", TRUE),
             stringsAsFactors = FALSE)
}

# Delete low-Phred bases.
#   scope "reads"   : delete the bases only; consensus is recomputed from the remaining reads
#   scope "orphans" : also delete the consensus position where no usable base remains
#   scope "columns" : also delete every consensus position that contained a low base
# Returns the updated mask/delcol and counts.
aln_delete_low <- function(aln, mask, delcol, thr, scope = "orphans") {
  nr <- nrow(aln$base); nc <- ncol(aln$base)
  colmat <- matrix(delcol, nr, nc, byrow = TRUE)
  covered <- !is.na(aln$base)
  low <- covered & !is.na(aln$q) & aln$q < thr & !mask & !colmat
  m2 <- mask | low
  newcols <- rep(FALSE, nc)
  if (identical(scope, "orphans")) {
    left <- colSums(covered & !m2 & !colmat)
    newcols <- !delcol & colSums(covered) > 0 & left == 0 & colSums(low) > 0
  } else if (identical(scope, "columns")) {
    newcols <- !delcol & colSums(low) > 0
  }
  list(mask = m2, delcol = delcol | newcols, n_bases = sum(low), n_cols = sum(newcols))
}

# ------------------------------------------------------------------------------
# Mapping reads to a reference
# ------------------------------------------------------------------------------
walk_alignment <- function(a1, a2, q_start, r_start) {
  qp <- q_start - 1L; rp <- r_start - 1L
  refcol <- integer(0); base <- character(0); idx <- integer(0)
  ins <- list(); cur <- NULL
  for (j in seq_along(a1)) {
    x <- a1[j]; y <- a2[j]
    if (x != "-") qp <- qp + 1L
    if (y == "-") {                      # insertion relative to the reference
      if (is.null(cur)) cur <- list(after = rp, seq = x) else cur$seq <- paste0(cur$seq, x)
      next
    }
    if (!is.null(cur)) { ins[[length(ins) + 1]] <- cur; cur <- NULL }
    rp <- rp + 1L
    refcol <- c(refcol, rp); base <- c(base, x)
    idx <- c(idx, if (x == "-") NA_integer_ else qp)
  }
  if (!is.null(cur)) ins[[length(ins) + 1]] <- cur
  list(refcol = refcol, base = base, idx = idx, ins = ins)
}

.map_run <- function(s, ref_seq, mat)
  .pw("pairwiseAlignment")(pattern = s, subject = ref_seq, type = "local",
                           substitutionMatrix = mat, gapOpening = 8, gapExtension = 4)

.empty_ins <- function()
  data.frame(Read = character(), After_ref_pos = integer(), Inserted = character(), stringsAsFactors = FALSE)

build_aln_ref <- function(raws, ref, P, labels = NULL) {
  if (nchar(ref$seq) > 30000)
    return(list(ok = FALSE, msg = paste0("The reference is too long (", nchar(ref$seq),
      " bp). Use the target region sequence (< 30 kb), not a whole genome.")))
  if (is.null(labels)) labels <- vapply(raws, function(r) r$name, "")
  mat <- .pw("nucleotideSubstitutionMatrix")(match = 2, mismatch = -3, baseOnly = FALSE)
  maps <- list(); notes <- character()

  for (j in seq_along(raws)) {
    lab <- labels[j]
    rd0 <- tryCatch(prep_read(raws[[j]], P), error = function(e) {
      notes <<- c(notes, paste0(lab, ": ", conditionMessage(e), " (skipped)")); NULL })
    if (is.null(rd0)) next
    rd1 <- rc_read(rd0)
    a0 <- .map_run(rd0$seq, ref$seq, mat); a1 <- .map_run(rd1$seq, ref$seq, mat)
    use1 <- as.numeric(.pw("score")(a1)) > as.numeric(.pw("score")(a0))
    rd <- if (use1) rd1 else rd0
    aln <- if (use1) a1 else a0
    s1 <- strsplit(as.character(.pw("alignedPattern")(aln)), "")[[1]]
    s2 <- strsplit(as.character(.pw("alignedSubject")(aln)), "")[[1]]
    qs <- BiocGenerics::start(.pw("pattern")(aln)); qe <- BiocGenerics::end(.pw("pattern")(aln))
    rs <- BiocGenerics::start(.pw("subject")(aln))
    both <- s1 != "-" & s2 != "-"
    ident <- if (any(both)) mean(s1[both] == s2[both]) else 0
    if (sum(both) < 50 || ident < 0.70) {
      notes <- c(notes, sprintf("%s: does not map to the reference (%d matching bp, identity %.0f%%); skipped.",
                                lab, sum(both), 100 * ident))
      next
    }
    w <- walk_alignment(s1, s2, qs, rs)
    unm5 <- qs - 1L; unm3 <- nchar(rd$seq) - qe
    if (unm5 > 20 || unm3 > 20)
      notes <- c(notes, sprintf("%s: %d bp (5' end) and %d bp (3' end) do not map to the reference.",
                                lab, unm5, unm3))
    maps[[length(maps) + 1]] <- c(list(read = rd, label = lab, rc = use1, ident = ident), w)
  }

  if (!length(maps))
    return(list(ok = FALSE, msg = paste("No read maps to the reference.", paste(notes, collapse = " "))))

  starts <- vapply(maps, function(m) min(m$refcol), 1L)
  maps <- maps[order(starts)]
  lo <- min(vapply(maps, function(m) min(m$refcol), 1L))
  hi <- max(vapply(maps, function(m) max(m$refcol), 1L))
  nc <- hi - lo + 1L; nr <- length(maps)

  base <- matrix(NA_character_, nr, nc)
  idx <- matrix(NA_integer_, nr, nc)
  for (i in seq_len(nr)) {
    kk <- maps[[i]]$refcol - lo + 1L
    base[i, kk] <- maps[[i]]$base
    idx[i, kk] <- maps[[i]]$idx
  }
  lbl <- make.unique(vapply(maps, function(m) m$label, ""))
  reads <- lapply(maps, function(m) m$read)
  M <- aln_matrices(reads, base, idx)

  ins <- do.call(rbind, lapply(seq_along(maps), function(i) {
    if (!length(maps[[i]]$ins)) return(NULL)
    data.frame(Read = lbl[i],
               After_ref_pos = vapply(maps[[i]]$ins, function(z) z$after, 1L),
               Inserted = vapply(maps[[i]]$ins, function(z) z$seq, ""), stringsAsFactors = FALSE)
  }))
  if (is.null(ins)) ins <- .empty_ins()
  if (nrow(ins))
    notes <- c(notes, sprintf("%d insertion(s) relative to the reference detected; they are NOT included in the consensus (see the SNP & variants tab).", nrow(ins)))

  aln <- list(ok = TRUE, mode = "ref", notes = notes, P = P, ref = ref,
              refbase = strsplit(ref$seq, "")[[1]][lo:hi], refpos = lo:hi,
              reads = M$reads, labels = lbl,
              strand = vapply(maps, function(m) if (m$rc) "Reverse" else "Forward", ""),
              base = M$base, idx = M$idx, q = M$q, ratio = M$ratio, top = M$top, sec = M$sec,
              insertions = ins,
              map_info = data.frame(Read = lbl,
                                    Strand = vapply(maps, function(m) if (m$rc) "reverse-complemented" else "as is", ""),
                                    Identity_to_ref = round(100 * vapply(maps, function(m) m$ident, 1), 1),
                                    stringsAsFactors = FALSE))
  aln$cols <- aln_decide(aln, NULL, NULL)
  aln
}

# ------------------------------------------------------------------------------
# Mode without a reference: Forward + Reverse pair
# ------------------------------------------------------------------------------
build_aln_denovo <- function(f_raw, r_raw, P, labels = c("Forward", "Reverse")) {
  res <- pair_alignment(f_raw, r_raw, P)
  if (!isTRUE(res$ok)) return(list(ok = FALSE, msg = res$msg))
  cols0 <- res$cols; n <- nrow(cols0)
  base <- rbind(cols0$Fwd, cols0$Rev)
  idx <- rbind(as.integer(cols0$Fwd_idx), as.integer(cols0$Rev_idx))
  for (i in 1:2) {
    nz <- which(!is.na(idx[i, ]))
    out <- seq_len(n) < min(nz) | seq_len(n) > max(nz)
    base[i, out] <- NA_character_
  }
  M <- aln_matrices(list(res$f, res$r), base, idx)
  notes <- sprintf("Overlap %d bp (identity %.1f%%); R %s.", res$overlap_len, 100 * res$overlap_identity,
                   if (res$r_was_rc) "was reverse-complemented" else "was used as is")
  aln <- list(ok = TRUE, mode = "denovo", notes = notes, P = P, ref = NULL, refbase = NULL,
              refpos = rep(NA_integer_, n), reads = M$reads, labels = labels,
              strand = c("Forward", "Reverse"),
              base = M$base, idx = M$idx, q = M$q, ratio = M$ratio, top = M$top, sec = M$sec,
              insertions = .empty_ins(), map_info = NULL)
  aln$cols <- aln_decide(aln, NULL, NULL)
  aln
}

# ------------------------------------------------------------------------------
# Final result (automatic + manual decisions) and derived tables
# ------------------------------------------------------------------------------
aln_final <- function(aln, cols, ov = character(), mask = NULL) {
  cf <- final_cols(cols, ov)
  cf$Source <- ifelse(cf$Type == "ColDel", "deleted", ifelse(cf$Manual, "manual", "rule"))
  B <- aln_eff_base(aln, mask)
  agree <- colSums(sweep(B, 2, cf$Final, "==") & !is.na(B))
  cf$Supported <- sprintf("%d/%d reads", as.integer(agree), as.integer(cf$Cov))
  cf
}

aln_variants <- function(aln, cf) {
  empty <- data.frame(RefPos = integer(), Reference = character(), Consensus = character(),
                      Kind = character(), Type = character(), Supported = character(),
                      Needs_review = logical(), stringsAsFactors = FALSE)
  if (is.null(aln$ref)) return(empty)
  d <- cf[!cf$Type %in% c("NoCov", "ColDel") & cf$Final != cf$RefBase, , drop = FALSE]
  if (!nrow(d)) return(empty)
  data.frame(RefPos = d$RefPos, Reference = d$RefBase, Consensus = d$Final,
             Kind = ifelse(d$Final == "-", "Deletion",
                           ifelse(d$Final %in% c("A", "C", "G", "T"), "Substitution", "Ambiguous (N/IUPAC)")),
             Type = d$Type, Supported = d$Supported, Needs_review = d$Review & !d$Manual,
             stringsAsFactors = FALSE)
}

# SNPs relative to the reference (substitutions + heterozygous/ambiguous candidates)
aln_snps <- function(aln, cf, mask = NULL) {
  empty <- data.frame(RefPos = integer(), Col = integer(), Ref = character(), Alt = character(),
                      Class = character(), Cov = integer(), Reads_alt = integer(), Reads_ref = integer(),
                      Phred_alt_min = numeric(), Phred_alt_mean = numeric(), Type = character(),
                      Source = character(), Needs_review = logical(), stringsAsFactors = FALSE)
  if (is.null(aln$ref)) return(empty)
  d <- which(!cf$Type %in% c("NoCov", "ColDel") & cf$Final != "-" & cf$Final != cf$RefBase)
  if (!length(d)) return(empty)
  B <- aln_eff_base(aln, mask)
  rows <- lapply(d, function(k) {
    alt <- cf$Final[k]; rf <- cf$RefBase[k]
    isb <- alt %in% c("A", "C", "G", "T")
    a_rd <- which(!is.na(B[, k]) & B[, k] == alt)
    r_rd <- which(!is.na(B[, k]) & B[, k] == rf)
    qa <- aln$q[a_rd, k]
    n_cov <- cf$Cov[k]
    cls <- if (isb) {
      if (length(a_rd) == n_cov) "SNP (all reads agree)"
      else if (length(a_rd) >= 2 || length(a_rd) / n_cov >= 0.5) "SNP (majority of reads)"
      else "SNP (single read; verify)"
    } else if (alt == "N") "Ambiguous (N)" else "Heterozygous / mixed-template candidate (IUPAC)"
    data.frame(RefPos = cf$RefPos[k], Col = k, Ref = rf, Alt = alt, Class = cls, Cov = n_cov,
               Reads_alt = length(a_rd), Reads_ref = length(r_rd),
               Phred_alt_min = if (length(qa)) round(min(qa, na.rm = TRUE)) else NA_real_,
               Phred_alt_mean = if (length(qa)) round(mean(qa, na.rm = TRUE), 1) else NA_real_,
               Type = cf$Type[k], Source = cf$Source[k], Needs_review = cf$Review[k] && !cf$Manual[k],
               stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

# Indels: deletions in the consensus, Indel-type columns, insertions relative to the reference
aln_indels <- function(aln, cf) {
  out <- list()
  runs <- function(ix) if (length(ix)) split(ix, cumsum(c(1, diff(ix) != 1))) else list()
  ref_txt <- function(g) if (all(is.na(cf$RefBase[g]))) "-" else paste(cf$RefBase[g], collapse = "")
  for (g in runs(which(cf$Final == "-" & !cf$Type %in% c("NoCov", "ColDel"))))
    out[[length(out) + 1]] <- data.frame(Kind = "Deletion in consensus", Col_start = min(g), Col_end = max(g),
      Ref_start = cf$RefPos[min(g)], Ref_end = cf$RefPos[max(g)], Length = length(g),
      Detail = paste("Reference bases:", ref_txt(g)), stringsAsFactors = FALSE)
  for (g in runs(which(cf$Type == "Indel")))
    out[[length(out) + 1]] <- data.frame(Kind = "Gap in some reads (possible indel)", Col_start = min(g),
      Col_end = max(g), Ref_start = cf$RefPos[min(g)], Ref_end = cf$RefPos[max(g)], Length = length(g),
      Detail = "Check the chromatogram and homopolymers", stringsAsFactors = FALSE)
  ins <- aln$insertions
  if (!is.null(ins) && nrow(ins))
    for (i in seq_len(nrow(ins)))
      out[[length(out) + 1]] <- data.frame(Kind = "Insertion relative to reference", Col_start = NA_integer_,
        Col_end = NA_integer_, Ref_start = ins$After_ref_pos[i], Ref_end = ins$After_ref_pos[i],
        Length = nchar(ins$Inserted[i]), Detail = paste0(ins$Read[i], ": ", ins$Inserted[i]), stringsAsFactors = FALSE)
  if (!length(out))
    return(data.frame(Kind = character(), Col_start = integer(), Col_end = integer(), Ref_start = integer(),
                      Ref_end = integer(), Length = integer(), Detail = character(), stringsAsFactors = FALSE))
  do.call(rbind, out)
}

# Differences of each read from the consensus ("difference data frame")
aln_read_diffs <- function(aln, cf, mask = NULL, max_rows = 5000) {
  B <- aln_eff_base(aln, mask)
  rows <- list()
  for (i in seq_len(nrow(B))) {
    d <- which(!is.na(B[i, ]) & B[i, ] != cf$Final & !cf$Type %in% c("NoCov", "ColDel"))
    if (length(d))
      rows[[length(rows) + 1]] <- data.frame(Read = aln$labels[i], Col = d, RefPos = cf$RefPos[d],
        Consensus = cf$Final[d], Read_base = B[i, d], Phred = round(aln$q[i, d]), stringsAsFactors = FALSE)
  }
  if (!length(rows))
    return(data.frame(Read = character(), Col = integer(), RefPos = integer(), Consensus = character(),
                      Read_base = character(), Phred = numeric(), stringsAsFactors = FALSE))
  utils::head(do.call(rbind, rows), max_rows)
}

aln_evidence_table <- function(aln, k, mask = NULL) {
  rows <- lapply(seq_len(nrow(aln$base)), function(i) {
    if (is.na(aln$base[i, k])) return(NULL)
    deleted <- !is.null(mask) && isTRUE(mask[i, k])
    data.frame(Read = aln$labels[i], Base = aln$base[i, k], Phred = round(aln$q[i, k]),
               Top_peak = aln$top[i, k], Second_peak = aln$sec[i, k],
               Ratio = round(aln$ratio[i, k], 2),
               Pos_in_read = if (is.na(aln$idx[i, k])) NA_integer_ else as.integer(aln$idx[i, k]),
               Status = if (deleted) "DELETED from alignment" else "used",
               stringsAsFactors = FALSE)
  })
  rows <- do.call(rbind, rows)
  if (is.null(rows)) data.frame(Info = "No read covers this position") else rows
}

aln_alignment_fasta <- function(aln, cf, mask = NULL) {
  B <- aln_eff_base(aln, mask)
  rowstr <- function(m) { m[is.na(m)] <- "-"; paste(m, collapse = "") }
  nm <- c(if (!is.null(aln$ref)) "Reference", aln$labels, "Consensus")
  sq <- c(if (!is.null(aln$ref)) paste(aln$refbase, collapse = ""),
          vapply(seq_len(nrow(B)), function(i) rowstr(B[i, ]), ""),
          paste(cf$Final, collapse = ""))
  paste0(">", nm, "\n", sq, collapse = "\n")
}
