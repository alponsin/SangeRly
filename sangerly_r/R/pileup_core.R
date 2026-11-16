# ==============================================================================
# SangeRly - Multi-read pileup (v4)
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
#   ins_cols: TRUE when the alignment has extra columns for insertions / read ends that run past the reference.
#             Those columns have refpos = NA and refbase = "-" (like the mapped-reads view of UGENE), so the
#             numbering of the reference positions never changes.
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
  isIns <- if (is.null(aln$refbase)) rep(FALSE, nc) else aln$refbase == "-"      # column the reference does not have
  dec <- vector("list", nc)
  for (k in seq_len(nc)) {
    cv <- which(!is.na(B[, k]))
    dec[[k]] <- if (delcol[k]) {
      list(type = "ColDel", base = "-", rec = "Position deleted by user (removed from the consensus)",
           reason = "Deleted", review = FALSE)
    } else if (isIns[k] && length(cv)) {
      b <- B[cv, k]; ng <- b != "-"
      ins_rec <- "Bases the reference does not have (insertion, or a read end past the reference): kept so that you decide - keep this position or delete it"
      if (!any(ng)) {
        list(type = "Ins", base = "-", rec = "No read has a base here any more", reason = "No inserted base", review = FALSE)
      } else if (sum(ng) * 2 < length(cv)) {
        list(type = "Ins", base = "-", rec = paste("Only a minority of the reads have a base here: left out of the consensus.", ins_rec),
             reason = paste(sprintf("%s=%s", aln$labels[cv], b), collapse = "; "), review = TRUE)
      } else {
        d <- decide_multi(aln$labels[cv][ng], b[ng], Q[cv, k][ng], R[cv, k][ng], S[cv, k][ng], T[cv, k][ng], P)
        list(type = "Ins", base = d$base, rec = ins_rec, reason = d$reason, review = TRUE)
      }
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
      if (is.null(cur)) cur <- list(after = rp, seq = x, idx = qp)
      else { cur$seq <- paste0(cur$seq, x); cur$idx <- c(cur$idx, qp) }
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

.map_run <- function(s, ref_seq, mat, type = "local")
  .pw("pairwiseAlignment")(pattern = s, subject = ref_seq, type = type,
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
  # "whole": every read is placed end to end on the reference ("overlap" alignment: only the ends that run past
  #          the reference are free). Poor or messy stretches stay in the alignment as mismatches / gaps / N,
  #          so the user decides what to delete or change (like the mapped-reads view of UGENE).
  # "local": only the best matching part of each read is kept (the older behaviour).
  whole <- !identical(P$map_mode, "local")
  atype <- if (whole) "overlap" else "local"
  # limits for accepting a read: set by the user (Reference-guided panel); without a setting 55% in "whole" mode, 70% in "local" mode
  min_id <- if (!is.null(P$map_min_id) && is.finite(P$map_min_id)) P$map_min_id / 100 else if (whole) 0.55 else 0.70
  min_len <- if (!is.null(P$map_min_len) && is.finite(P$map_min_len)) P$map_min_len else 50

  for (j in seq_along(raws)) {
    lab <- labels[j]
    rd0 <- tryCatch(prep_read(raws[[j]], P), error = function(e) {
      notes <<- c(notes, paste0(lab, ": ", conditionMessage(e), " (skipped)")); NULL })
    if (is.null(rd0)) next
    rd1 <- rc_read(rd0)
    a0 <- .map_run(rd0$seq, ref$seq, mat, atype); a1 <- .map_run(rd1$seq, ref$seq, mat, atype)
    use1 <- as.numeric(.pw("score")(a1)) > as.numeric(.pw("score")(a0))
    rd <- if (use1) rd1 else rd0
    aln <- if (use1) a1 else a0
    s1 <- strsplit(as.character(.pw("alignedPattern")(aln)), "")[[1]]
    s2 <- strsplit(as.character(.pw("alignedSubject")(aln)), "")[[1]]
    qs <- BiocGenerics::start(.pw("pattern")(aln)); qe <- BiocGenerics::end(.pw("pattern")(aln))
    rs <- BiocGenerics::start(.pw("subject")(aln))
    both <- s1 != "-" & s2 != "-"
    ident <- if (any(both)) mean(s1[both] == s2[both]) else 0
    if (sum(both) < min_len || ident < min_id) {
      notes <- c(notes, sprintf("%s: does not map to the reference (%d aligned bp, identity %.0f%%; the limits are %.0f%% and %d bp, lower them in the Reference-guided panel); skipped.",
                                lab, sum(both), 100 * ident, 100 * min_id, as.integer(min_len)))
      next
    }
    w <- walk_alignment(s1, s2, qs, rs)
    unm5 <- qs - 1L; unm3 <- nchar(rd$seq) - qe
    if (unm5 > 20 || unm3 > 20)
      notes <- c(notes, sprintf("%s: %d bp (5' end) and %d bp (3' end) do not map to the reference.",
                                lab, unm5, unm3))
    if (whole && isTRUE(P$keep_ins) && length(w$refcol)) {
      # the ends of the read that run past the reference are not in the alignment: keep them as insertions
      sq <- strsplit(toupper(rd$seq), "")[[1]]
      if (unm5 > 0)
        w$ins <- c(list(list(after = rs - 1L, seq = paste(sq[seq_len(unm5)], collapse = ""), idx = seq_len(unm5))), w$ins)
      if (unm3 > 0)
        w$ins[[length(w$ins) + 1L]] <- list(after = max(w$refcol),
                                            seq = paste(sq[(qe + 1L):length(sq)], collapse = ""),
                                            idx = (qe + 1L):length(sq))
    }
    maps[[length(maps) + 1]] <- c(list(read = rd, label = lab, rc = use1, ident = ident), w)
  }

  if (!length(maps))
    return(list(ok = FALSE, msg = paste("No read maps to the reference.", paste(notes, collapse = " "))))

  starts <- vapply(maps, function(m) min(m$refcol), 1L)
  maps <- maps[order(starts)]
  lo <- min(vapply(maps, function(m) min(m$refcol), 1L))
  hi <- max(vapply(maps, function(m) max(m$refcol), 1L))
  if (isTRUE(P$full_ref)) { lo <- 1L; hi <- nchar(ref$seq) }      # whole reference; positions without a read become N
  nr <- length(maps)

  # Extra columns for insertions and for read ends past the reference (the reference shows "-" there).
  # They hang after reference position `site`; the numbering of the reference positions does not change.
  keep_ins <- isTRUE(P$keep_ins)
  ins_all <- lapply(maps, function(m) if (keep_ins) Filter(function(z) z$after >= lo - 1L && z$after <= hi, m$ins) else list())
  after_v <- unlist(lapply(ins_all, function(l) vapply(l, function(z) as.integer(z$after), 1L)))
  len_v <- unlist(lapply(ins_all, function(l) vapply(l, function(z) nchar(z$seq), 1L)))
  if (length(after_v)) {
    Kt <- tapply(len_v, after_v, max)
    sites <- as.integer(names(Kt)); K <- as.integer(Kt)
    o <- order(sites); sites <- sites[o]; K <- K[o]
  } else { sites <- integer(0); K <- integer(0) }
  cumK <- cumsum(K)
  cumKprev <- if (length(K)) c(0L, cumK)[seq_along(K)] else integer(0)
  shift <- function(p) {
    ci <- findInterval(p - 1L, sites)
    if (!length(cumK)) return(rep(0L, length(p)))
    ifelse(ci == 0L, 0L, cumK[pmax(ci, 1L)])
  }
  cidx <- function(p) (p - lo + 1L) + shift(p)               # alignment column of reference position p
  site_start <- (sites - lo + 1L) + cumKprev + 1L            # first extra column of each site
  nc <- (hi - lo + 1L) + sum(K)

  base <- matrix(NA_character_, nr, nc)
  idx <- matrix(NA_integer_, nr, nc)
  for (i in seq_len(nr)) {
    kk <- cidx(maps[[i]]$refcol)
    base[i, kk] <- maps[[i]]$base
    idx[i, kk] <- maps[[i]]$idx
    for (z in ins_all[[i]]) {                                 # bases inserted by this read
      j <- match(as.integer(z$after), sites)
      L <- nchar(z$seq)
      cc <- if (z$after == lo - 1L) site_start[j] + (K[j] - L) + seq_len(L) - 1L   # left end: against the reference start
            else site_start[j] + seq_len(L) - 1L
      base[i, cc] <- strsplit(z$seq, "")[[1]]
      idx[i, cc] <- z$idx
    }
    if (length(sites)) {                                      # a read that spans an insertion site but has no base there: "-"
      rmin <- min(maps[[i]]$refcol); rmax <- max(maps[[i]]$refcol)
      for (j in which(sites >= rmin & sites < rmax)) {
        blk <- site_start[j] + seq_len(K[j]) - 1L
        cur <- base[i, blk]; cur[is.na(cur)] <- "-"; base[i, blk] <- cur
      }
    }
  }
  refpos_v <- rep(NA_integer_, nc); refpos_v[cidx(lo:hi)] <- lo:hi
  refbase_v <- rep("-", nc); refbase_v[cidx(lo:hi)] <- strsplit(ref$seq, "")[[1]][lo:hi]
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
    notes <- c(notes, if (length(sites))
      sprintf("%d insertion(s) / read end(s) past the reference are kept as extra columns (the reference shows '-'; the reference numbering does not change). They are flagged for review: keep them or delete those positions.", nrow(ins))
    else sprintf("%d insertion(s) relative to the reference detected; they are NOT included in the consensus (see the SNP & variants tab).", nrow(ins)))

  notes <- c(notes, sprintf("Alignment columns: reference positions %d to %d of %d (%s); reads placed %s.",
                            lo, hi, nchar(ref$seq),
                            if (isTRUE(P$full_ref)) "whole reference region" else "only the part covered by reads",
                            if (whole) "end to end" else "locally"))
  aln <- list(ok = TRUE, mode = "ref", notes = notes, P = P, ref = ref,
              refbase = refbase_v, refpos = refpos_v, ins_cols = length(sites) > 0,
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
             Kind = ifelse(d$RefBase == "-", "Insertion / read end past the reference",
                    ifelse(d$Final == "-", "Deletion",
                           ifelse(d$Final %in% c("A", "C", "G", "T"), "Substitution", "Ambiguous (N/IUPAC)"))),
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
  d <- which(!cf$Type %in% c("NoCov", "ColDel") & cf$Final != "-" & cf$Final != cf$RefBase & cf$RefBase != "-")
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
  not_ins <- is.na(cf$RefBase) | cf$RefBase != "-"                 # columns that exist in the reference
  for (g in runs(which(cf$Final == "-" & !cf$Type %in% c("NoCov", "ColDel") & not_ins)))
    out[[length(out) + 1]] <- data.frame(Kind = "Deletion in consensus", Col_start = min(g), Col_end = max(g),
      Ref_start = cf$RefPos[min(g)], Ref_end = cf$RefPos[max(g)], Length = length(g),
      Detail = paste("Reference bases:", ref_txt(g)), stringsAsFactors = FALSE)
  for (g in runs(which(cf$Type == "Indel")))
    out[[length(out) + 1]] <- data.frame(Kind = "Gap in some reads (possible indel)", Col_start = min(g),
      Col_end = max(g), Ref_start = cf$RefPos[min(g)], Ref_end = cf$RefPos[max(g)], Length = length(g),
      Detail = "Check the chromatogram and homopolymers", stringsAsFactors = FALSE)
  if (!is.null(aln$refbase)) {                                       # insertions kept as extra columns
    prev_ref <- function(k) { v <- cf$RefPos[seq_len(k - 1L)]; v <- v[!is.na(v)]; if (length(v)) max(v) else 0L }
    for (g in runs(which(!not_ins & cf$Final != "-" & !cf$Type %in% c("NoCov", "ColDel"))))
      out[[length(out) + 1]] <- data.frame(Kind = "Insertion / read end past the reference (kept in the consensus)",
        Col_start = min(g), Col_end = max(g), Ref_start = prev_ref(min(g)), Ref_end = prev_ref(min(g)), Length = length(g),
        Detail = paste0("after reference position ", prev_ref(min(g)), ": ", paste(cf$Final[g], collapse = "")), stringsAsFactors = FALSE)
  }
  ins <- if (isTRUE(aln$ins_cols)) NULL else aln$insertions
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

# Alignment (reference, reads, consensus) as FASTA. Positions deleted by the user (Type "ColDel") are removed from
# EVERY row (reference, reads and consensus), so all rows stay the same length and shift together.
aln_alignment_fasta <- function(aln, cf, mask = NULL, drop_deleted = TRUE) {
  B <- aln_eff_base(aln, mask)
  keep <- if (drop_deleted && !is.null(cf$Type) && length(cf$Type) == ncol(B)) cf$Type != "ColDel" else rep(TRUE, ncol(B))
  rowstr <- function(m) { m[is.na(m)] <- "-"; paste(m[keep], collapse = "") }
  nm <- c(if (!is.null(aln$ref)) "Reference", aln$labels, "Consensus")
  sq <- c(if (!is.null(aln$ref)) paste(aln$refbase[keep], collapse = ""),
          vapply(seq_len(nrow(B)), function(i) rowstr(B[i, ]), ""),
          paste(cf$Final[keep], collapse = ""))
  paste0(">", nm, "\n", sq, collapse = "\n")
}
