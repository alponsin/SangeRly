# primer_core.R - functions (no Shiny) for the "Primer mining and design" tab of SangeRly
#
#   * sequences:  BOLD (v5 portal API), NCBI nucleotide (Entrez), FASTA, MSA tab, assembly consensus
#   * design:     conventional PCR pair, Sanger sequencing primers (with primer walking), mRNA exon-exon junction,
#                 degenerate primers from an alignment; optional probe (hydrolysis probe, qPCR)
#   * quality:    Tm (SantaLucia 1998 nearest-neighbour), GC, runs, 3' stability, hairpin, self- / hetero-dimer (estimates),
#                 binding in a set of sequences (mismatches, 3' end)
#   * BLAST:      NCBI web BLAST with short-query settings, plus predicted amplicons
#
# Packages: httr, jsonlite (network); Biostrings + DECIPHER or msa (only to align mined sequences here).

# ---- IUPAC helpers --------------------------------------------------------------------------------------
.pr_iupac <- c(A = 1L, C = 2L, G = 4L, T = 8L, M = 3L, R = 5L, S = 6L, V = 7L, W = 9L, Y = 10L,
               H = 11L, K = 12L, D = 13L, B = 14L, N = 15L)                      # bit mask: A=1 C=2 G=4 T=8
.pr_code_of_mask <- c("A", "C", "M", "G", "R", "S", "V", "T", "W", "Y", "H", "K", "D", "B", "N")  # index = mask 1..15
.pr_nbits <- c(1L, 1L, 2L, 1L, 2L, 2L, 3L, 1L, 2L, 2L, 3L, 2L, 3L, 3L, 4L)      # number of bases of each mask
.pr_comp <- c(A = "T", C = "G", G = "C", T = "A", M = "K", R = "Y", S = "S", V = "B", W = "W", Y = "R",
              H = "D", K = "M", D = "H", B = "V", N = "N")
.pr_pc <- function(k) c(0L, .pr_nbits)[k + 1L]                                   # number of bases of a mask (0 = none)

.pr_clean <- function(x) {                                                       # upper case, U -> T, only IUPAC letters
  x <- gsub("U", "T", toupper(as.character(x)), fixed = TRUE)
  gsub("[^ACGTMRSVWYHKDBN]", "", x)
}
.pr_rc <- function(x) vapply(as.character(x), function(s) {
  if (is.na(s) || !nzchar(s)) return("")
  paste(rev(.pr_comp[strsplit(s, "")[[1]]]), collapse = "")
}, "", USE.NAMES = FALSE)
.pr_masks <- function(s) unname(.pr_iupac[strsplit(s, "")[[1]]])
.pr_deg <- function(s) prod(.pr_nbits[.pr_masks(s)])
.pr_gc <- function(s) { m <- .pr_masks(s); 100 * mean(.pr_pc(bitwAnd(m, 6L)) / .pr_pc(m)) }
.pr_ngc <- function(s) sum(.pr_pc(bitwAnd(.pr_masks(s), 6L)) / .pr_pc(.pr_masks(s)))
.pr_maxrun <- function(s) max(rle(strsplit(s, "")[[1]])$lengths)

# up to max_out spread-out variants of a (degenerate) primer
.pr_variants <- function(s, max_out = 16L) {
  m <- .pr_masks(s)
  opts <- lapply(m, function(k) c("A", "C", "G", "T")[bitwAnd(k, c(1L, 2L, 4L, 8L)) > 0L])
  ln <- lengths(opts); tot <- prod(ln)
  idx <- if (tot <= max_out) seq_len(tot) - 1 else unique(floor(seq(0, tot - 1, length.out = max_out)))
  vapply(idx, function(ix) {
    out <- character(length(opts)); r <- ix
    for (p in rev(seq_along(opts))) { out[p] <- opts[[p]][(r %% ln[p]) + 1]; r <- r %/% ln[p] }
    paste(out, collapse = "")
  }, "")
}

# ---- thermodynamics (SantaLucia 1998, unified parameters) -------------------------------------------------------
# dH in kcal/mol, dS in cal/(K mol), for the 5'->3' dinucleotide of the top strand
.pr_nn <- list(
  AA = c(-7.9, -22.2), TT = c(-7.9, -22.2), AT = c(-7.2, -20.4), TA = c(-7.2, -21.3),
  CA = c(-8.5, -22.7), TG = c(-8.5, -22.7), GT = c(-8.4, -22.4), AC = c(-8.4, -22.4),
  CT = c(-7.8, -21.0), AG = c(-7.8, -21.0), GA = c(-8.2, -22.2), TC = c(-8.2, -22.2),
  CG = c(-10.6, -27.2), GC = c(-9.8, -24.4), GG = c(-8.0, -19.9), CC = c(-8.0, -19.9))

pr_params <- function(...) {
  p <- list(len_min = 18L, len_opt = 20L, len_max = 25L, tm_min = 57, tm_opt = 60, tm_max = 63,
            gc_min = 40, gc_opt = 50, gc_max = 60, max_run = 4L, max_gc3 = 3L, stab3_max = 9,
            hp_min = -3.5, sd_min = -7, sd3_min = -4, hd_min = -7, hd3_min = -4, max_dtm = 3,
            prod_min = 100L, prod_max = 700L, na = 50, mg = 1.5, dntp = 0.2, conc = 200,
            n3_max = 2L, run_flag = 4L, struct_tm_max = 47, n_out = 5L)   # quality checks: 3' pairing (bp), run length, Tm of a structure
  a <- list(...)
  for (k in names(a)) p[[k]] <- a[[k]]
  p
}

# Tm of an unambiguous A/C/G/T oligo (degC): nearest neighbour, SantaLucia 1998 salt correction
# (monovalent + 120*sqrt(Mg - dNTP) as sodium equivalent), strand concentration conc (nM), C/4 for non-self-complementary
.pr_tm1 <- function(s, prm) {
  n <- nchar(s); if (n < 2L) return(NA_real_)
  b <- strsplit(s, "")[[1]]
  dh <- 0; ds <- 0
  for (i in seq_len(n - 1L)) { v <- .pr_nn[[paste0(b[i], b[i + 1L])]]; dh <- dh + v[1]; ds <- ds + v[2] }
  for (e in c(b[1], b[n])) {
    if (e %in% c("G", "C")) { dh <- dh + 0.1; ds <- ds - 2.8 } else { dh <- dh + 2.3; ds <- ds + 4.1 }
  }
  selfc <- identical(s, .pr_rc(s))
  if (selfc) ds <- ds - 1.4
  na_eq <- (prm$na + 120 * sqrt(max(prm$mg - prm$dntp, 0))) / 1000
  ds <- ds + 0.368 * (n - 1L) * log(na_eq)
  x <- if (selfc) 1 else 4
  dh * 1000 / (ds + 1.987 * log(prm$conc * 1e-9 / x)) - 273.15
}
# mean / min / max Tm over the variants of a degenerate primer
pr_tm <- function(s, prm = pr_params()) {
  v <- .pr_variants(s, 16L)
  t <- vapply(v, function(z) .pr_tm1(z, prm), 0)
  c(mean = mean(t), min = min(t), max = max(t))
}

# free energy (kcal/mol, 37 degC) of the nearest-neighbour stacks of a perfect duplex
.pr_dg_nn <- function(s) {
  b <- strsplit(s, "")[[1]]; n <- length(b); g <- 0
  if (n < 2L) return(0)
  for (i in seq_len(n - 1L)) { v <- .pr_nn[[paste0(b[i], b[i + 1L])]]; g <- g + v[1] - 310.15 * v[2] / 1000 }
  g
}
.pr_dg_seq <- function(s) {                      # stacks + initiation of both ends
  b <- strsplit(s, "")[[1]]; n <- length(b)
  .pr_dg_nn(s) + sum(vapply(c(b[1], b[n]), function(e) if (e %in% c("G", "C")) 0.98 else 1.03, 0))
}
.pr_stab3 <- function(s) { n <- nchar(s); -.pr_dg_nn(substr(s, max(1L, n - 4L), n)) }   # 3' end stability (positive; primer3 style)

.pr_loop_pen <- function(n) {
  tab <- c(3.5, 3.5, 3.3, 4.0, 4.2, 4.3, 4.5)    # loops of 3..9 nt
  if (n <= 9L) tab[n - 2L] else 4.5 + 1.08 * log(n / 9)
}
# most stable hairpin (stem >= 3 bp, loop >= 3 nt); any = whole primer, end = stem contains the 3' end base
.pr_hairpin <- function(s) {
  b <- strsplit(s, "")[[1]]; n <- length(b)
  best <- 0; best3 <- 0
  if (n < 9L) return(list(any = 0, end = 0))
  for (i in seq_len(n)) {
    if (i + 8L > n) break
    for (j in (i + 8L):n) {
      k <- 0L
      while (j - i - 2L * k >= 4L && .pr_comp[[b[i + k]]] == b[j - k]) k <- k + 1L
      if (k >= 3L) {
        for (kk in 3L:k) {
          g <- .pr_dg_nn(substr(s, i, i + kk - 1L)) + .pr_loop_pen(j - i - 2L * kk + 1L)
          if (g < best) best <- g
          if (j == n && g < best3) best3 <- g
        }
      }
    }
  }
  list(any = best, end = best3)
}
# most stable ungapped duplex between a and b (a with b); any = best, end = best one that contains a 3' end
.pr_dimer <- function(a, b = a) {
  A <- strsplit(a, "")[[1]]; B <- strsplit(.pr_rc(b), "")[[1]]
  na <- length(A); nb <- length(B)
  best <- 0; best3 <- 0; n3 <- 0L; nany <- 0L          # n3 = longest perfect run of base pairs that includes a 3' end base
  for (d in seq(-(nb - 1L), na - 1L)) {
    lo <- max(1L, 1L + d); hi <- min(na, nb + d)
    if (hi - lo + 1L < 3L) next
    i <- lo:hi; j <- i - d
    r <- rle(A[i] == B[j])
    ends <- cumsum(r$lengths); starts <- ends - r$lengths + 1L
    for (k in which(r$values & r$lengths >= 3L)) {
      s <- i[starts[k]]; e <- i[ends[k]]
      g <- .pr_dg_seq(paste(A[s:e], collapse = ""))
      if (g < best) best <- g
      if (r$lengths[k] > nany) nany <- r$lengths[k]
      if (e == na || j[starts[k]] == 1L) { if (g < best3) best3 <- g; if (r$lengths[k] > n3) n3 <- r$lengths[k] }
    }
  }
  list(any = best, end = best3, n3 = n3, nany = nany)
}

# ---- quality of one primer ------------------------------------------------------------------------------------
pr_eval <- function(s, prm = pr_params()) {
  s <- .pr_clean(s); n <- nchar(s)
  vars <- .pr_variants(s, 8L)
  tm <- pr_tm(s, prm)
  hp <- c(0, 0); sd <- c(0, 0)
  for (v in vars) {
    h <- .pr_hairpin(v); d <- .pr_dimer(v, v)
    hp <- pmin(hp, c(h$any, h$end)); sd <- pmin(sd, c(d$any, d$end))
  }
  tail5 <- substr(s, max(1L, n - 4L), n)
  m <- .pr_masks(s)
  list(seq = s, len = n, gc = .pr_gc(s), tm = tm[["mean"]], tm_min = tm[["min"]], tm_max = tm[["max"]],
       deg = .pr_deg(s), run = .pr_maxrun(s), gc3 = round(.pr_ngc(tail5), 1),
       stab3 = .pr_stab3(vars[1]), hp = hp[1], hp3 = hp[2], sd = sd[1], sd3 = sd[2],
       dinuc = grepl("([ACGT]{2})\\1{4,}", vars[1]), deg3 = .pr_pc(m[n]) > 1L)
}
# status OK / CHECK / FAIL with the reasons
pr_verdict <- function(e, prm = pr_params(), probe = FALSE, ref_tm = NA_real_) {
  notes <- character(); sev <- 0L
  add <- function(txt, s) { notes <<- c(notes, txt); sev <<- max(sev, s) }
  if (probe) {                       # hydrolysis (TaqMan-type) probe: other rules than a primer (no GC clamp, 5' G quenches the dye, Tm above the primers)
    if (e$len < 18L) add("probe short (<18 nt)", 1L)
    if (e$len > 30L) add("probe long (>30 nt: quenching is poor, consider an internal quencher)", 1L)
    if (substr(e$seq, 1L, 1L) == "G") add("5' base is G (quenches the reporter dye)", 1L)
    if (e$gc < 30 || e$gc > 80) add(sprintf("GC %.0f%% outside 30-80", e$gc), 2L)
    else if (e$gc < 35 || e$gc > 65) add(sprintf("GC %.0f%% (aim 35-65)", e$gc), 1L)
    if (!is.na(ref_tm)) {
      d <- e$tm - ref_tm
      if (d < 3) add(sprintf("Tm only %.1f C above the primers (aim +6 to +10 so that it binds first)", d), 2L)
      else if (d < 5) add(sprintf("Tm %.1f C above the primers (aim +6 to +10)", d), 1L)
      else if (d > 15) add(sprintf("Tm %.1f C above the primers (very high)", d), 1L)
    }
    if (e$run >= 4L && grepl("G{4,}", e$seq)) add("run of 4+ G (poor synthesis / G-quadruplex)", 1L)
    else if (e$run >= 6L) add(sprintf("run of %d identical bases", e$run), 2L)
    else if (e$run >= 5L) add(sprintf("run of %d identical bases", e$run), 1L)
    if (e$hp < -4) add(sprintf("stable hairpin (%.1f kcal/mol)", e$hp), 2L)
    else if (e$hp < -2) add(sprintf("hairpin %.1f kcal/mol", e$hp), 1L)
    if (e$sd < -9) add(sprintf("strong self-dimer (%.1f)", e$sd), 2L)
    else if (e$sd < -6) add(sprintf("self-dimer %.1f", e$sd), 1L)
    if (e$deg > 16) add(sprintf("degeneracy %g (a probe should have few variants)", e$deg), 1L)
    notes <- c(notes, "probe rules apply: the primer rules on the GC clamp and the 3' end are not used")
    return(list(status = c("OK", "CHECK", "FAIL")[sev + 1L], notes = paste(notes, collapse = "; ")))
  }
  if (e$len < 15L) add("short (<15 nt)", 1L)
  if (e$len > 35L) add("long (>35 nt)", 1L)
  if (e$tm < prm$tm_min || e$tm > prm$tm_max) add(sprintf("Tm %.1f outside %g-%g", e$tm, prm$tm_min, prm$tm_max), 1L)
  if (e$gc < 30 || e$gc > 70) add(sprintf("GC %.0f%% extreme", e$gc), 2L)
  else if (e$gc < prm$gc_min || e$gc > prm$gc_max) add(sprintf("GC %.0f%% outside %g-%g", e$gc, prm$gc_min, prm$gc_max), 1L)
  if (e$gc3 < 1) add("no G/C in the last 5 nt (weak 3' end)", 1L)
  if (e$gc3 > prm$max_gc3) add("GC-rich 3' end (mispriming risk)", 1L)
  if (e$stab3 > prm$stab3_max) add(sprintf("3' end too stable (%.1f)", e$stab3), 1L)
  if (isTRUE(e$gc2 == 0)) add("no G/C in the last 2 nt (no GC clamp)", 1L)
  if (e$run >= 6L) add(sprintf("run of %d identical bases", e$run), 2L)
  else if (e$run >= prm$run_flag) add(sprintf("run of %d identical bases", e$run), 1L)
  if (isTRUE(e$dinuc)) add("dinucleotide repeat", 1L)
  if (e$hp < -4) add(sprintf("stable hairpin (%.1f kcal/mol)", e$hp), 2L)
  else if (e$hp < -2) add(sprintf("hairpin %.1f kcal/mol", e$hp), 1L)
  if (e$hp3 < -2) add(sprintf("hairpin at the 3' end (%.1f)", e$hp3), 1L)
  if (e$sd < -9) add(sprintf("strong self-dimer (%.1f)", e$sd), 2L)
  else if (e$sd < -6) add(sprintf("self-dimer %.1f", e$sd), 1L)
  if (e$sd3 < -5) add(sprintf("3' self-dimer (%.1f)", e$sd3), 2L)
  else if (e$sd3 < -3) add(sprintf("3' self-dimer %.1f", e$sd3), 1L)
  if (isTRUE(e$n3 > prm$n3_max)) add(sprintf("3' end pairs with the primer itself over %d bases (self-dimer / primer-dimer risk)", e$n3), 1L)
  if (isTRUE(e$sd_tm > prm$struct_tm_max)) add(sprintf("self-dimer Tm %.0f C (> %g)", e$sd_tm, prm$struct_tm_max), 1L)
  if (e$deg > 256) add(sprintf("degeneracy %g (very high)", e$deg), 1L)
  if (e$deg > 1 && (e$tm_max - e$tm_min) > 8) add(sprintf("Tm spread %.1f degC between variants", e$tm_max - e$tm_min), 1L)
  if (isTRUE(e$deg3)) add("degenerate base at the 3' end", 1L)
  list(status = c("OK", "CHECK", "FAIL")[sev + 1L], notes = paste(notes, collapse = "; "))
}

# ---- primer candidates on one template ---------------------------------------------------------------------------
# footprint [a, b] on the template; strand "F": primer = footprint; "R": primer = reverse complement of the footprint
.pr_make_cand <- function(s, a, b, strand, prm) {
  gc <- .pr_gc(s); if (gc < prm$gc_min || gc > prm$gc_max) return(NULL)
  if (.pr_maxrun(s) > prm$max_run) return(NULL)
  tm <- .pr_tm1(s, prm); if (is.na(tm) || tm < prm$tm_min || tm > prm$tm_max) return(NULL)
  n <- nchar(s)
  gc3 <- .pr_ngc(substr(s, n - 4L, n)); if (gc3 > prm$max_gc3) return(NULL)
  st3 <- .pr_stab3(s); if (st3 > prm$stab3_max) return(NULL)
  if (grepl("([ACGT]{2})\\1{4,}", s)) return(NULL)
  hp <- .pr_hairpin(s); if (hp$any < prm$hp_min || hp$end < prm$hp_min + 1) return(NULL)
  sd <- .pr_dimer(s, s); if (sd$any < prm$sd_min || sd$end < prm$sd3_min) return(NULL)
  pen <- abs(tm - prm$tm_opt) + 0.2 * abs(n - prm$len_opt) + 0.05 * abs(gc - prm$gc_opt) +
         0.3 * (gc3 < 1) + 0.1 * max(0, -hp$any) + 0.05 * max(0, -sd$any)
  data.frame(strand = strand, start = a, end = b, seq = s, len = n, tm = tm, gc = gc, stab3 = st3,
             hairpin = hp$any, selfdimer = sd$any, pen = pen, stringsAsFactors = FALSE)
}
# all acceptable primers whose footprint lies inside [from, to]; filt(a, b) may reject footprints
pr_primer_cands <- function(tpl, prm, from = 1L, to = nchar(tpl), strand = "F", filt = NULL) {
  n <- nchar(tpl); from <- max(1L, as.integer(from)); to <- min(n, as.integer(to))
  out <- list(); k <- 0L
  for (L in prm$len_min:prm$len_max) {
    if (to - from + 1L < L) next
    for (a in from:(to - L + 1L)) {
      b <- a + L - 1L
      if (!is.null(filt) && !isTRUE(filt(a, b))) next
      fp <- substr(tpl, a, b)
      if (grepl("[^ACGT]", fp)) next
      s <- if (strand == "F") fp else .pr_rc(fp)
      z <- .pr_make_cand(s, a, b, strand, prm)
      if (!is.null(z)) { k <- k + 1L; out[[k]] <- z }
    }
  }
  if (!k) return(NULL)
  do.call(rbind, out)
}

# F x R pairs; product = from the start of F to the end of R
pr_pair_cands <- function(Fc, Rc, prm, top = 60L, keep = 30L) {
  if (is.null(Fc) || is.null(Rc) || !nrow(Fc) || !nrow(Rc)) return(NULL)
  Fc <- Fc[order(Fc$pen), , drop = FALSE]; Fc <- Fc[seq_len(min(nrow(Fc), top)), , drop = FALSE]
  Rc <- Rc[order(Rc$pen), , drop = FALSE]; Rc <- Rc[seq_len(min(nrow(Rc), top)), , drop = FALSE]
  popt <- (prm$prod_min + prm$prod_max) / 2
  rows <- list(); k <- 0L
  for (i in seq_len(nrow(Fc))) for (j in seq_len(nrow(Rc))) {
    prod <- Rc$end[j] - Fc$start[i] + 1L
    if (prod < prm$prod_min || prod > prm$prod_max) next
    if (Rc$start[j] <= Fc$end[i]) next
    dtm <- abs(Fc$tm[i] - Rc$tm[j]); if (dtm > prm$max_dtm) next
    hd <- .pr_dimer(Fc$seq[i], Rc$seq[j]); if (hd$any < prm$hd_min || hd$end < prm$hd3_min) next
    k <- k + 1L
    rows[[k]] <- data.frame(F_seq = Fc$seq[i], F_start = Fc$start[i], F_end = Fc$end[i], F_tm = Fc$tm[i], F_gc = Fc$gc[i],
                            R_seq = Rc$seq[j], R_start = Rc$start[j], R_end = Rc$end[j], R_tm = Rc$tm[j], R_gc = Rc$gc[j],
                            product = prod, dTm = dtm, heterodimer = hd$any,
                            penalty = Fc$pen[i] + Rc$pen[j] + 0.5 * dtm + 0.002 * abs(prod - popt) + 0.05 * max(0, -hd$any),
                            stringsAsFactors = FALSE)
  }
  if (!k) return(NULL)
  d <- do.call(rbind, rows)
  d <- d[order(d$penalty), , drop = FALSE]
  d <- d[seq_len(min(nrow(d), keep)), , drop = FALSE]
  rownames(d) <- NULL; d
}

# hydrolysis probe inside [from, to] (Tm near tm_target, 5' base not G, GC 30-80 %, no GGGG; strand with more C than G preferred)
pr_probe <- function(tpl, from, to, prm, tm_target, top = 3L) {
  from <- max(1L, as.integer(from)); to <- min(nchar(tpl), as.integer(to))
  rows <- list(); k <- 0L
  for (L in 18L:30L) {
    if (to - from + 1L < L) next
    for (a in from:(to - L + 1L)) {
      fp <- substr(tpl, a, a + L - 1L)
      if (grepl("[^ACGT]", fp)) next
      for (st in c("+", "-")) {
        s <- if (st == "+") fp else .pr_rc(fp)
        if (substr(s, 1L, 1L) == "G") next
        gc <- .pr_gc(s); if (gc < 30 || gc > 80) next
        if (grepl("GGGG", s)) next
        tm <- .pr_tm1(s, prm); if (is.na(tm) || tm < tm_target - 3 || tm > tm_target + 3) next
        b <- strsplit(s, "")[[1]]
        pen <- abs(tm - tm_target) + 0.05 * abs(gc - 50) + 0.5 * (sum(b == "C") < sum(b == "G"))
        k <- k + 1L
        rows[[k]] <- data.frame(start = a, end = a + L - 1L, strand = st, seq = s, tm = tm, gc = gc, pen = pen, stringsAsFactors = FALSE)
      }
    }
  }
  if (!k) return(NULL)
  d <- do.call(rbind, rows); d <- d[order(d$pen), , drop = FALSE]
  d[seq_len(min(top, nrow(d))), , drop = FALSE]
}
.pr_add_probe <- function(pairs, tpl, prm, delta, n_pairs = 10L) {
  pairs$Probe <- NA_character_; pairs$Probe_Tm <- NA_real_; pairs$Probe_start <- NA_integer_; pairs$Probe_strand <- NA_character_
  for (i in seq_len(min(n_pairs, nrow(pairs)))) {
    tmt <- mean(c(pairs$F_tm[i], pairs$R_tm[i])) + delta
    pb <- pr_probe(tpl, pairs$F_end[i] + 2L, pairs$R_start[i] - 2L, prm, tmt, top = 1L)
    if (!is.null(pb) && nrow(pb)) {
      pairs$Probe[i] <- pb$seq[1]; pairs$Probe_Tm[i] <- pb$tm[1]; pairs$Probe_start[i] <- pb$start[1]; pairs$Probe_strand[i] <- pb$strand[1]
    }
  }
  pairs
}

# ---- conventional PCR -------------------------------------------------------------------------------------------
# tf, tt: region that must be inside the product (NA = no constraint)
pr_design_pcr <- function(tpl, tf = NA, tt = NA, prm = pr_params(), probe = FALSE, probe_delta = 8, inside = FALSE) {
  tpl <- .pr_clean(tpl); n <- nchar(tpl)
  if (n < prm$prod_min) stop("The sequence (", n, " nt) is shorter than the minimum product size.")
  if ((is.na(tf) || is.na(tt)) && n > 20000L)
    stop("The template is long (", n, " nt): give the region (From / To) or choose a gene of an annotated GenBank record.")
  if (!is.na(tf) && !is.na(tt) && inside) {            # primers INSIDE the region: amplify part of it
    tf <- max(1L, as.integer(tf)); tt <- min(n, as.integer(tt))
    if (tt - tf + 1L < prm$prod_min) stop("The region (", tt - tf + 1L, " nt) is shorter than the minimum product size.")
    Ff <- tf; Ft <- tt - prm$len_min; Rf <- tf + prm$len_min; Rt <- tt
  } else if (is.na(tf) || is.na(tt)) {
    Ff <- 1L; Ft <- n - prm$len_min; Rf <- prm$len_min; Rt <- n
  } else {
    tf <- max(1L, as.integer(tf)); tt <- min(n, as.integer(tt))
    if (tt < tf) stop("'To' is smaller than 'From' in the target region.")
    if (tt - tf + 1L > prm$prod_max)
      stop("The region (", tt - tf + 1L, " nt) is longer than the maximum product size (", prm$prod_max,
           "): raise the maximum product size, or choose 'Primers inside the region' to amplify part of it.")
    Ff <- max(1L, tt - prm$prod_max + 1L); Ft <- tf - 1L
    Rf <- tt + 1L; Rt <- min(n, tf + prm$prod_max - 1L)
    if (Ft < prm$len_min) stop("The target starts too close to the 5' end: no room for a forward primer.")
    if (Rf > n - prm$len_min + 1L) stop("The target ends too close to the 3' end: no room for a reverse primer.")
  }
  Fc <- pr_primer_cands(tpl, prm, Ff, Ft, "F"); Rc <- pr_primer_cands(tpl, prm, Rf, Rt, "R")
  nF <- if (is.null(Fc)) 0L else nrow(Fc); nR <- if (is.null(Rc)) 0L else nrow(Rc)
  pairs <- pr_pair_cands(Fc, Rc, prm, keep = prm$n_out)
  if (!is.null(pairs) && probe) pairs <- .pr_add_probe(pairs, tpl, prm, probe_delta)
  list(pairs = pairs, nF = nF, nR = nR, n = n)
}

# ---- Sanger sequencing primers (with primer walking) ------------------------------------------------------------
# first primer 'offset' nt before the region, next ones every (readlen - overlap) nt; returns footprints on tpl
.pr_walk <- function(tpl, s, e, prm, offset, readlen, overlap, strand) {
  n <- nchar(tpl)
  if (strand == "R") { tpl <- .pr_rc(tpl); ss <- n - e + 1L; e <- n - s + 1L; s <- ss }
  out <- list(); k <- 0L
  target <- s - offset
  repeat {
    hi <- min(n, max(target, prm$len_min))
    got <- NULL
    for (w in c(80L, 200L, 400L)) {
      lo <- max(1L, hi - w)
      c1 <- pr_primer_cands(tpl, prm, lo - prm$len_max + 1L, hi, "F")
      if (!is.null(c1)) c1 <- c1[c1$end >= lo & c1$end <= hi, , drop = FALSE]
      if (!is.null(c1) && nrow(c1)) { got <- c1[which.min(c1$pen), , drop = FALSE]; break }
    }
    if (is.null(got)) break
    k <- k + 1L
    cf <- got$end + offset; ct <- min(n, got$end + readlen)
    out[[k]] <- data.frame(got, cover_from = cf, cover_to = ct, stringsAsFactors = FALSE)
    if (ct >= e || k >= 30L) break
    target <- ct - overlap - offset
    if (target <= got$end) break
  }
  if (!k) return(NULL)
  d <- do.call(rbind, out)
  if (strand == "R") {
    a <- n - d$end + 1L; b <- n - d$start + 1L; cf <- n - d$cover_to + 1L; ct <- n - d$cover_from + 1L
    d$start <- a; d$end <- b; d$cover_from <- cf; d$cover_to <- ct
  }
  d$strand <- strand
  d
}
pr_design_seq <- function(tpl, tf = NA, tt = NA, prm = pr_params(), offset = 50L, readlen = 700L, overlap = 100L, both = TRUE) {
  tpl <- .pr_clean(tpl); n <- nchar(tpl)
  if ((is.na(tf) || is.na(tt)) && n > 20000L) stop("The template is long (", n, " nt): give the region (From / To).")
  s <- if (is.na(tf)) 1L else max(1L, as.integer(tf)); e <- if (is.na(tt)) n else min(n, as.integer(tt))
  if (e < s) stop("'To' is smaller than 'From'.")
  if (readlen <= offset + overlap + 50L) stop("Read length must be larger than offset + overlap + 50.")
  parts <- list(.pr_walk(tpl, s, e, prm, offset, readlen, overlap, "F"))
  if (both) parts[[2]] <- .pr_walk(tpl, s, e, prm, offset, readlen, overlap, "R")
  parts <- Filter(Negate(is.null), parts)
  if (!length(parts)) return(NULL)
  d <- do.call(rbind, parts)
  d$name <- ave(seq_len(nrow(d)), d$strand, FUN = seq_along)
  d$name <- paste0("Seq", d$strand, d$name)
  d <- d[order(d$strand, d$start), , drop = FALSE]
  d$gc <- round(d$gc, 1); d$tm <- round(d$tm, 1)
  rownames(d) <- NULL
  d[, c("name", "strand", "start", "end", "seq", "len", "tm", "gc", "cover_from", "cover_to")]
}

# ---- mRNA / cDNA: primers across an exon-exon junction --------------------------------------------------------------
# junc: positions of the last base of an exon (the junction lies between j and j + 1)
# span: "forward", "reverse", "either" = one primer spans a junction; "flank" = primers in two different exons
pr_design_junction <- function(tpl, junc, prm = pr_params(), span = "either", ov5 = 7L, ov3 = 4L, probe = FALSE, probe_delta = 8) {
  tpl <- .pr_clean(tpl); n <- nchar(tpl)
  junc <- sort(unique(as.integer(junc))); junc <- junc[!is.na(junc) & junc >= 1L & junc < n]
  if (!length(junc)) stop("No valid junction position (they must be between 1 and ", n - 1L, ").")
  inside <- function(a, b) !any(junc >= a & junc < b)          # footprint lies in ONE exon
  rows <- list(); k <- 0L
  add <- function(p, type, j) if (!is.null(p) && nrow(p)) { p$junction <- j; p$type <- type; k <<- k + 1L; rows[[k]] <<- p }
  for (j in junc) {
    if (span %in% c("forward", "either")) {
      Fc <- pr_primer_cands(tpl, prm, j - prm$len_max + 1L, j + prm$len_max, "F",
                            filt = function(a, b) a <= j - ov5 + 1L && b >= j + ov3)
      Rc <- pr_primer_cands(tpl, prm, j + ov3 + 1L, min(n, j + prm$prod_max), "R", filt = inside)
      add(pr_pair_cands(Fc, Rc, prm), "forward primer spans the junction", j)
    }
    if (span %in% c("reverse", "either")) {
      Rc <- pr_primer_cands(tpl, prm, j - prm$len_max + 1L, j + prm$len_max, "R",
                            filt = function(a, b) a <= j - ov3 + 1L && b >= j + ov5)
      Fc <- pr_primer_cands(tpl, prm, max(1L, j - prm$prod_max), j, "F", filt = inside)
      add(pr_pair_cands(Fc, Rc, prm), "reverse primer spans the junction", j)
    }
    if (identical(span, "flank")) {
      Fc <- pr_primer_cands(tpl, prm, max(1L, j - prm$prod_max + 1L), j, "F", filt = inside)
      Rc <- pr_primer_cands(tpl, prm, j + 1L, min(n, j + prm$prod_max - 1L), "R", filt = inside)
      add(pr_pair_cands(Fc, Rc, prm), "primers in two different exons", j)
    }
  }
  if (!k) return(list(pairs = NULL, n = n, junc = junc))
  d <- do.call(rbind, rows)
  d <- d[!duplicated(paste(d$F_seq, d$R_seq)), , drop = FALSE]
  d <- d[order(d$penalty), , drop = FALSE]
  d <- d[seq_len(min(prm$n_out, nrow(d))), , drop = FALSE]
  d$exon_F <- vapply(d$F_end, function(x) 1L + sum(junc < x), 1L)
  d$exon_R <- vapply(d$R_start, function(x) 1L + sum(junc < x), 1L)
  tx <- .pr_tx_from_junc(junc, n)
  d$exons_F <- vapply(seq_len(nrow(d)), function(i) .pr_exon_span(tx, d$F_start[i], d$F_end[i]), "")
  d$exons_R <- vapply(seq_len(nrow(d)), function(i) .pr_exon_span(tx, d$R_start[i], d$R_end[i]), "")
  d$exons_prod <- vapply(seq_len(nrow(d)), function(i) .pr_exon_span(tx, d$F_start[i], d$R_end[i], "-"), "")
  d$junc_label <- vapply(d$junction, function(j) { k <- match(j, junc); sprintf("%d | %d", k, k + 1L) }, "")
  rownames(d) <- NULL
  if (probe) d <- .pr_add_probe(d, tpl, prm, probe_delta)
  list(pairs = d, n = n, junc = junc)
}

# ---- degenerate primers from an alignment -------------------------------------------------------------------------------
# aln: named character vector of aligned sequences (same width, '-' = gap)
pr_aln_ctx <- function(aln) {
  aln <- toupper(gsub("U", "T", as.character(aln), fixed = TRUE)); nm <- names(aln)
  w <- unique(nchar(aln)); if (length(w) != 1L) stop("The sequences are not aligned (different lengths).")
  M <- do.call(rbind, strsplit(aln, ""))
  n <- nrow(M)
  bases <- c("A", "C", "G", "T")
  cnt <- vapply(bases, function(b) colSums(M == b), numeric(w))               # w x 4
  if (w == 1L) cnt <- matrix(cnt, nrow = 1L)
  tot <- rowSums(cnt)
  SM <- matrix(0L, n, w)                                                       # IUPAC mask of every base (gap = 0)
  for (ch in names(.pr_iupac)) SM[M == ch] <- .pr_iupac[[ch]]
  ng <- SM > 0L
  list(M = M, SM = SM, cnt = cnt, tot = tot, n = n, w = w, names = nm, gapfrac = 1 - tot / n,
       cumng = t(apply(ng, 1, cumsum)), ident = ifelse(tot > 0, apply(cnt, 1, max) / pmax(tot, 1), 0),
       major = bases[apply(cnt, 1, which.max)])
}
# mask of each column: the bases with frequency >= minfreq
.pr_col_masks <- function(ctx, minfreq) {
  w <- ctx$w; mk <- integer(w)
  for (c in seq_len(w)) {
    if (ctx$tot[c] < 1) next
    f <- ctx$cnt[c, ] / ctx$tot[c]
    keep <- f >= minfreq & ctx$cnt[c, ] > 0
    if (!any(keep)) keep[which.max(ctx$cnt[c, ])] <- TRUE
    mk[c] <- sum(c(1L, 2L, 4L, 8L)[keep])
  }
  mk
}
# scan windows; strand "F" (3' end on the right) or "R" (3' end on the left); returns candidate windows with coverage
.pr_win_scan <- function(ctx, mk, Ls, strand, prm, max_deg, max_per_col, max_gap, mm, three = 3L, from = 1L, to = ctx$w,
                         fix3 = TRUE, max_deg3 = 1L, keep_top = 150L) {
  nb <- .pr_pc(mk); w <- ctx$w
  okcol <- ctx$gapfrac <= max_gap & nb >= 1L & nb <= max_per_col
  badc <- c(0L, cumsum(!okcol)); lg <- c(0, cumsum(log2(pmax(nb, 1L))))
  out <- list(); k <- 0L
  for (L in Ls) {
    if (min(w, to) - L + 1L < max(1L, from)) next
    for (a0 in seq.int(max(1L, from), min(w, to) - L + 1L)) {
      e <- a0 + L - 1L
      if (badc[e + 1L] - badc[a0] > 0L) next
      deg <- 2^(lg[e + 1L] - lg[a0]); if (deg > max_deg) next
      cols3 <- if (strand == "F") (e - 4L):e else a0:(a0 + 4L)
      last <- if (strand == "F") e else a0
      if (fix3 && nb[last] > 1L) next
      if (sum(nb[cols3] > 1L) > max_deg3) next
      X <- ctx$SM[, a0:e, drop = FALSE]
      P <- matrix(mk[a0:e], ctx$n, L, byrow = TRUE)
      bad <- matrix(bitwAnd(X, P) == 0L, nrow = ctx$n)          # bitwAnd drops the dimensions
      present <- rowSums(X == 0L) <= 0.1 * L
      if (sum(present) < 2L) next
      i3 <- if (strand == "F") (L - three + 1L):L else 1L:three
      okseq <- present & rowSums(bad) <= mm & rowSums(bad[, i3, drop = FALSE]) == 0L
      cov <- sum(okseq) / sum(present)
      k <- k + 1L
      out[[k]] <- list(start = a0, end = e, len = L, deg = deg, cov = cov, npres = sum(present), ok = okseq, present = present)
    }
  }
  if (!k) return(NULL)
  sc <- vapply(out, function(z) z$cov - 0.02 * log2(z$deg), 0)
  out <- out[order(-sc)]
  out[seq_len(min(keep_top, length(out)))]
}
.pr_win_seq <- function(mk, st, en, strand) {
  s <- paste(.pr_code_of_mask[mk[st:en]], collapse = "")
  if (strand == "R") .pr_rc(s) else s
}
# degenerate pairs (and optional probe) from an alignment
pr_design_degen <- function(aln, prm = pr_params(), max_deg = 64, minfreq = 0.1, max_gap = 0.2, mm = 0L, min_cov = 0.8,
                            from = 1L, to = NA, probe = FALSE, probe_delta = 8, ref = 1L, max_per_col = 3L, top_win = 40L) {
  ctx <- pr_aln_ctx(aln)
  if (ctx$n < 3L) stop("Give an alignment of at least 3 sequences.")
  if (is.na(to)) to <- ctx$w
  mk <- .pr_col_masks(ctx, minfreq)
  Ls <- prm$len_min:prm$len_max
  evalwin <- function(win, strand) {
    rows <- list(); k <- 0L
    for (z in win) {
      if (z$cov < min_cov) next
      s <- .pr_win_seq(mk, z$start, z$end, strand)
      e <- tryCatch(pr_eval(s, prm), error = function(err) NULL)
      if (is.null(e)) next
      if (e$tm < prm$tm_min || e$tm > prm$tm_max) next
      if (e$gc < prm$gc_min - 5 || e$gc > prm$gc_max + 5) next
      if (e$run > prm$max_run || e$hp < -4 || e$sd < -9 || e$sd3 < -5) next
      if (e$deg > 1 && (e$tm_max - e$tm_min) > 8) next
      k <- k + 1L
      rows[[k]] <- list(z = z, e = e, seq = s, strand = strand,
                        score = z$cov - 0.02 * log2(z$deg) - 0.05 * abs(e$tm - prm$tm_opt))
    }
    rows <- rows[order(-vapply(rows, function(r) r$score, 0))]
    rows[seq_len(min(top_win, length(rows)))]
  }
  wF <- .pr_win_scan(ctx, mk, Ls, "F", prm, max_deg, max_per_col, max_gap, mm, from = from, to = to)
  wR <- .pr_win_scan(ctx, mk, Ls, "R", prm, max_deg, max_per_col, max_gap, mm, from = from, to = to)
  Fw <- evalwin(wF, "F"); Rw <- evalwin(wR, "R")
  if (!length(Fw) || !length(Rw)) {
    return(list(pairs = NULL, ctx = ctx, mk = mk, nF = length(Fw), nR = length(Rw)))
  }
  refpos <- function(col) as.integer(ctx$cumng[ref, col])
  rows <- list(); k <- 0L
  for (f in Fw) for (r in Rw) {
    if (r$z$start <= f$z$end) next
    pl <- ctx$cumng[, r$z$end] - (if (f$z$start > 1L) ctx$cumng[, f$z$start - 1L] else 0L)
    okp <- ctx$SM[, f$z$start] > 0L & ctx$SM[, r$z$end] > 0L
    if (!any(okp)) next
    prod <- stats::median(pl[okp])
    if (prod < prm$prod_min || prod > prm$prod_max) next
    dtm <- abs(f$e$tm - r$e$tm); if (dtm > prm$max_dtm + 1) next
    hd <- .pr_dimer(.pr_variants(f$seq, 1L), .pr_variants(r$seq, 1L))
    if (hd$any < prm$hd_min - 1 || hd$end < prm$hd3_min - 1) next
    both <- f$z$ok & r$z$ok
    bc <- mean(both[f$z$present | r$z$present])
    k <- k + 1L
    rows[[k]] <- data.frame(F_seq = f$seq, F_col = f$z$start, F_end_col = f$z$end, F_deg = f$z$deg, F_tm = f$e$tm,
                            F_tm_range = sprintf("%.1f-%.1f", f$e$tm_min, f$e$tm_max), F_gc = f$e$gc, F_cov = f$z$cov,
                            R_seq = r$seq, R_col = r$z$start, R_end_col = r$z$end, R_deg = r$z$deg, R_tm = r$e$tm,
                            R_tm_range = sprintf("%.1f-%.1f", r$e$tm_min, r$e$tm_max), R_gc = r$e$gc, R_cov = r$z$cov,
                            product = prod, prod_range = sprintf("%d-%d", min(pl[okp]), max(pl[okp])),
                            both_cov = bc, dTm = dtm, heterodimer = hd$any, ref_F = refpos(f$z$start), ref_R = refpos(r$z$end),
                            score = f$score + r$score - 0.3 * dtm, stringsAsFactors = FALSE)
  }
  if (!k) return(list(pairs = NULL, ctx = ctx, mk = mk, nF = length(Fw), nR = length(Rw)))
  d <- do.call(rbind, rows)
  d <- d[order(-d$score), , drop = FALSE]
  d <- d[seq_len(min(prm$n_out, nrow(d))), , drop = FALSE]; rownames(d) <- NULL
  if (probe) {
    d$Probe <- NA_character_; d$Probe_Tm <- NA_real_; d$Probe_col <- NA_integer_; d$Probe_cov <- NA_real_; d$Probe_note <- NA_character_
    np <- min(10L, nrow(d))
    lo <- min(d$F_end_col[seq_len(np)]) + 2L; hi <- max(d$R_col[seq_len(np)]) - 2L
    # three rounds, from strict to relaxed; every window between the primers is scanned (all are kept), then the best one for each pair is chosen.
    # filters of a probe: 5' base not G (it quenches the dye), no GGGG, GC 30-80 %, Tm near the primers + probe_delta
    tiers <- list(list(note = "strict", cov = min_cov, mm = mm, deg = 4, per = 2L, d3 = 4L, tol = 3),
                  list(note = "relaxed", cov = max(0.5, min_cov - 0.2), mm = max(1L, mm), deg = 8, per = 3L, d3 = 5L, tol = 5),
                  list(note = "relaxed (low coverage)", cov = 0.5, mm = max(1L, mm), deg = 16, per = 3L, d3 = 5L, tol = 8))
    for (tr in tiers) {
      todo <- which(is.na(d$Probe[seq_len(np)])); if (!length(todo)) break
      cand <- list()
      if (hi - lo >= 17L) for (st in c("F", "R")) {
        wp <- .pr_win_scan(ctx, mk, 18L:30L, st, prm, tr$deg, tr$per, max_gap, tr$mm, from = lo, to = hi, fix3 = FALSE, max_deg3 = tr$d3, keep_top = 100000L)
        for (z in wp) {
          if (z$cov < tr$cov) next
          sq <- .pr_win_seq(mk, z$start, z$end, st)
          if (substr(sq, 1L, 1L) == "G" || grepl("GGGG", sq, fixed = TRUE)) next
          gc <- .pr_gc(sq); if (gc < 30 || gc > 80) next
          cand[[length(cand) + 1L]] <- list(z = z, seq = sq, tm = pr_tm(sq, prm)[["mean"]])
        }
      }
      if (!length(cand)) next
      for (i in todo) {
        tmt <- mean(c(d$F_tm[i], d$R_tm[i])) + probe_delta
        ok <- Filter(function(cd) cd$z$start >= d$F_end_col[i] + 2L && cd$z$end <= d$R_col[i] - 2L && abs(cd$tm - tmt) <= tr$tol, cand)
        if (!length(ok)) next
        sc <- vapply(ok, function(cd) cd$z$cov - 0.02 * log2(cd$z$deg) - 0.02 * abs(cd$tm - tmt), 0)
        bst <- ok[[which.max(sc)]]
        d$Probe[i] <- bst$seq; d$Probe_Tm[i] <- bst$tm; d$Probe_col[i] <- bst$z$start; d$Probe_cov[i] <- bst$z$cov; d$Probe_note[i] <- tr$note
      }
    }
  }
  list(pairs = d, ctx = ctx, mk = mk, nF = length(Fw), nR = length(Rw))
}

# align (unaligned) sequences: DECIPHER, otherwise msa
pr_align <- function(seqs) {
  seqs <- seqs[nchar(seqs) > 0]
  if (length(seqs) < 2L) stop("At least 2 sequences are needed to align.")
  if (!requireNamespace("Biostrings", quietly = TRUE)) stop("The 'Biostrings' package is required (BiocManager::install(\"Biostrings\")).")
  ss <- Biostrings::DNAStringSet(toupper(gsub("[^A-Za-z]", "", seqs)))
  names(ss) <- make.unique(names(seqs))
  if (requireNamespace("DECIPHER", quietly = TRUE)) {
    a <- DECIPHER::AlignSeqs(ss, verbose = FALSE)
  } else if (requireNamespace("msa", quietly = TRUE)) {
    a <- as(msa::msa(ss, method = "ClustalW", type = "dna"), "DNAStringSet")
  } else stop("Install DECIPHER or msa (BiocManager::install(\"DECIPHER\")) to align the sequences.")
  out <- as.character(a); names(out) <- names(ss); out
}

# ---- sequences: FASTA, BOLD, NCBI ---------------------------------------------------------------------------------------
pr_parse_fasta <- function(txt) {
  ln <- unlist(strsplit(paste(txt, collapse = "\n"), "\r?\n")); ln <- trimws(ln); ln <- ln[nzchar(ln)]
  if (!length(ln)) return(character())
  if (!any(startsWith(ln, ">"))) return(c(sequence = toupper(gsub("[^A-Za-z.-]", "", paste(ln, collapse = "")))))
  hdr <- which(startsWith(ln, ">")); out <- character(length(hdr))
  for (k in seq_along(hdr)) {
    a <- hdr[k] + 1L; b <- if (k < length(hdr)) hdr[k + 1L] - 1L else length(ln)
    out[k] <- if (a <= b) toupper(gsub("[^A-Za-z.-]", "", paste(ln[a:b], collapse = ""))) else ""
  }
  names(out) <- sub("^>\\s*", "", ln[hdr])
  out[nzchar(out)]
}
# every NCBI request goes through .pr_net: HTTP/1.1 (avoids many "Connection was reset" errors) + automatic retries with back-off
.pr_log <- function(...) { m <- paste0(format(Sys.time(), "%H:%M:%S"), "  ", paste0(..., collapse = "")); .pr_state$log <- c(utils::tail(.pr_state$log, 40L), m); invisible(m) }
.pr_state <- new.env()                                          # last transport errors, for the error messages
.pr_cfg <- function() httr::config(http_version = 2L)          # 2L = CURL_HTTP_VERSION_1_1
.pr_net <- function(verb, url, ..., tries = 4L) {
  r <- try(httr::RETRY(verb, url, ..., httr::user_agent("SangeRly (primer tab)"), config = .pr_cfg(),
                       times = tries, pause_base = 3, pause_cap = 20, terminate_on = c(400L, 401L, 403L, 404L), quiet = TRUE), silent = TRUE)
  if (inherits(r, "try-error")) {
    msg <- trimws(gsub("\\s+", " ", sub("^Error[^:]*:\\s*", "", as.character(r))))
    stop(sprintf("Cannot reach NCBI (%s). Tried %d times. Check the internet connection / VPN / proxy / firewall and try again in a minute.", msg, tries), call. = FALSE)
  }
  r
}
# fallback when httr / libcurl cannot connect (e.g. Windows + proxy / antivirus resetting the connection): the system download method
# (wininet = the Windows / Internet Explorer connection settings), which often works where libcurl does not
.pr_fallback_get <- function(url, query = list(), timeout = 60) {
  qs <- if (length(query)) paste0("?", paste0(names(query), "=", vapply(query, function(v) utils::URLencode(as.character(v)[1], reserved = TRUE), ""), collapse = "&")) else ""
  tf <- tempfile(); on.exit(unlink(tf), add = TRUE)
  old <- options(timeout = max(timeout, 60)); on.exit(options(old), add = TRUE)
  errs <- character()
  for (m in c(if (.Platform$OS.type == "windows") "wininet", "libcurl")) {
    ok <- tryCatch({ suppressWarnings(utils::download.file(paste0(url, qs), tf, method = m, quiet = TRUE, mode = "wb")); file.exists(tf) && file.size(tf) > 0 },
                   error = function(e) { errs <<- c(errs, sprintf("%s: %s", m, trimws(conditionMessage(e)))); FALSE })
    if (isTRUE(ok)) return(paste(readLines(tf, warn = FALSE, encoding = "UTF-8"), collapse = "\n"))
  }
  .pr_state$fallback_err <- paste(errs, collapse = " | ")
  NULL
}
.pr_http <- function(url, query = list(), timeout = 60) {
  r <- try(.pr_net("GET", url, query = query, httr::timeout(timeout)), silent = TRUE)
  if (inherits(r, "try-error")) {
    x <- .pr_fallback_get(url, query, timeout)
    if (is.null(x)) stop(conditionMessage(attr(r, "condition")), call. = FALSE)
    return(structure(list(text = x), class = "pr_resp"))
  }
  if (httr::status_code(r) >= 400L) stop(sprintf("HTTP %d from %s", httr::status_code(r), sub("\\?.*$", "", url)), call. = FALSE)
  r
}
.pr_text <- function(r) if (inherits(r, "pr_resp")) r$text else httr::content(r, as = "text", encoding = "UTF-8")
.pr_collapse <- function(x) {            # a column that may hold lists / vectors -> one string per row
  if (is.null(x)) return(character())
  vapply(x, function(v) if (length(v) == 0L || all(is.na(v))) "" else paste(unlist(v), collapse = "; "), "", USE.NAMES = FALSE)
}

# BOLD v5 portal API: tax:<name> -> query_id -> documents (JSON, 1000 per page). Marker filter is applied here.
pr_bold_fetch <- function(taxon, markers = character(), max_n = 300L, per_species = 0L, progress = function(m) NULL) {
  if (!requireNamespace("httr", quietly = TRUE) || !requireNamespace("jsonlite", quietly = TRUE)) stop("Install httr and jsonlite.")
  base <- "https://portal.boldsystems.org/api"
  taxon <- trimws(taxon); if (!nzchar(taxon)) stop("Type a taxon name (for example Danaus plexippus, Coleoptera, Aves).")
  qid_of <- function(q) {
    r <- try(.pr_http(paste0(base, "/query"), list(query = q, extent = "full"), 60), silent = TRUE)
    if (inherits(r, "try-error")) return(NULL)
    txt <- .pr_text(r)
    x <- try(jsonlite::fromJSON(txt), silent = TRUE)
    id <- if (!inherits(x, "try-error") && is.list(x) && !is.null(x$query_id)) x$query_id
          else if (!inherits(x, "try-error") && is.character(x) && length(x) == 1L) x else NULL
    if (is.null(id) || !nzchar(id)) NULL else id
  }
  progress("BOLD: preparing the query ...")
  qid <- qid_of(paste0("tax:", taxon))
  tc <- paste0(toupper(substr(taxon, 1L, 1L)), tolower(substring(taxon, 2L)))      # BOLD names are case-sensitive: "zingiberaceae" -> "Zingiberaceae"
  if (is.null(qid) && tc != taxon) { taxon <- tc; qid <- qid_of(paste0("tax:", taxon)) }
  if (is.null(qid)) {                    # let BOLD resolve the name to a rank (preprocessor), then retry
    r <- try(.pr_http(paste0(base, "/query/preprocessor"), list(query = paste0("tax:", taxon)), 60), silent = TRUE)
    if (!inherits(r, "try-error")) {
      m <- regmatches(.pr_text(r), regexpr("tax:[a-z]+:[A-Za-z0-9 ._-]+", .pr_text(r)))
      if (length(m)) qid <- qid_of(m[1])
    }
  }
  if (is.null(qid)) stop("BOLD did not accept '", taxon, "'. Check the spelling (scientific name, e.g. Danaus plexippus) or try a higher rank.")
  recs <- list(); got <- 0L; start <- 0L; scanned <- 0L; last_first <- ""
  cap <- 20000L
  repeat {
    progress(sprintf("BOLD: reading records %d-%d (kept %d) ...", start + 1L, start + 1000L, got))
    r <- try(.pr_http(sprintf("%s/documents/%s", base, qid), list(start = start, length = 1000L), 120), silent = TRUE)
    if (inherits(r, "try-error")) { if (start == 0L) stop("BOLD could not deliver the records: ", conditionMessage(attr(r, "condition"))) else break }
    x <- try(jsonlite::fromJSON(.pr_text(r), simplifyVector = TRUE), silent = TRUE)
    d <- if (inherits(x, "try-error")) NULL else if (is.data.frame(x)) x else if (is.list(x) && is.data.frame(x$data)) x$data else NULL
    if (is.null(d) || !NROW(d)) break
    first <- as.character(d[[1]][1]); if (identical(first, last_first)) break
    last_first <- first
    scanned <- scanned + NROW(d)
    nuc <- .pr_collapse(d$nuc)
    mk <- if (!is.null(d$marker_code)) .pr_collapse(d$marker_code) else rep("", NROW(d))
    sp <- if (!is.null(d$species)) .pr_collapse(d$species) else rep("", NROW(d))
    sel <- nzchar(nuc) & (length(markers) == 0L | mk %in% markers)
    if (any(sel)) {
      pid <- if (!is.null(d$processid)) .pr_collapse(d$processid) else as.character(seq_len(NROW(d)) + start)
      bin <- if (!is.null(d$bin_uri)) .pr_collapse(d$bin_uri) else rep("", NROW(d))
      acc <- if (!is.null(d$insdc_acs)) .pr_collapse(d$insdc_acs) else rep("", NROW(d))
      prim <- rep("", NROW(d))
      for (cn in grep("primer", names(d), ignore.case = TRUE, value = TRUE)) prim <- paste(prim, .pr_collapse(d[[cn]]))
      z <- data.frame(ID = pid, Species = sp, Marker = mk, BIN = bin, Accession = acc, Primers = trimws(prim),
                      Sequence = gsub("[^A-Za-z]", "", toupper(nuc)), stringsAsFactors = FALSE)[sel, , drop = FALSE]
      recs[[length(recs) + 1L]] <- z; got <- got + nrow(z)
    }
    if (got >= max_n * 3L || scanned >= cap || NROW(d) < 1000L) break
    start <- start + 1000L
  }
  if (!length(recs)) stop("BOLD returned no sequence", if (length(markers)) paste0(" for marker(s) ", paste(markers, collapse = ", ")), " for '", taxon, "'.")
  z <- do.call(rbind, recs)
  z$Sequence[nchar(z$Sequence) < 50L] <- ""; z <- z[nzchar(z$Sequence), , drop = FALSE]
  if (per_species > 0L) {                                  # at most N sequences per species (the longest ones)
    z <- z[order(-nchar(z$Sequence)), , drop = FALSE]
    z <- z[ave(seq_len(nrow(z)), ifelse(nzchar(z$Species), z$Species, z$ID), FUN = seq_along) <= per_species, , drop = FALSE]
  }
  z <- z[seq_len(min(nrow(z), max_n)), , drop = FALSE]
  z$Length <- nchar(z$Sequence); z$Source <- "BOLD"
  rownames(z) <- NULL
  z
}
# Entrez nucleotide: term in NCBI syntax, e.g. Danaus plexippus[Organism] AND COI[Gene]
pr_ncbi_fetch <- function(term, max_n = 200L, email = "", progress = function(m) NULL) {
  if (!requireNamespace("httr", quietly = TRUE) || !requireNamespace("jsonlite", quietly = TRUE)) stop("Install httr and jsonlite.")
  term <- trimws(term); if (!nzchar(term)) stop("Type an NCBI search, e.g. Danaus plexippus[Organism] AND COI[Gene].")
  base <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
  com <- list(tool = "SangeRly"); if (nzchar(trimws(email))) com$email <- trimws(email)
  progress("NCBI: searching ...")
  r <- .pr_http(paste0(base, "esearch.fcgi"), c(list(db = "nuccore", term = term, retmax = max_n, retmode = "json"), com))
  ids <- jsonlite::fromJSON(.pr_text(r))$esearchresult$idlist
  if (!length(ids)) stop("No NCBI nucleotide record matches this search.")
  recs <- character()
  for (chunk in split(ids, ceiling(seq_along(ids) / 100L))) {
    progress(sprintf("NCBI: downloading %d of %d sequences ...", length(recs) + length(chunk), length(ids)))
    r <- .pr_http(paste0(base, "efetch.fcgi"), c(list(db = "nuccore", id = paste(chunk, collapse = ","), rettype = "fasta", retmode = "text"), com), 120)
    recs <- c(recs, pr_parse_fasta(.pr_text(r)))
    Sys.sleep(0.4)
  }
  nm <- unname(names(recs)); recs <- unname(recs)
  z <- data.frame(ID = sub("\\s.*$", "", nm), Species = sub("^\\S+\\s+", "", nm), Marker = "", BIN = "", Accession = sub("\\s.*$", "", nm),
                  Primers = "", Sequence = gsub("[^A-Za-z]", "", recs), stringsAsFactors = FALSE)
  z$Length <- nchar(z$Sequence); z$Source <- "NCBI"
  rownames(z) <- NULL; z
}
# primer names + sequences written inside BOLD records, e.g. "LepF1 (ATTCAACCAATCATAAAGATATTGG)"
pr_bold_primers <- function(recs) {
  if (is.null(recs) || !"Primers" %in% names(recs)) return(NULL)
  tx <- recs$Primers[nzchar(recs$Primers)]; if (!length(tx)) return(NULL)
  hits <- unlist(regmatches(tx, gregexpr("[A-Za-z0-9_.-]+\\s*\\(([ACGTRYSWKMBDHVNI]{12,})\\)", tx)))
  if (!length(hits)) return(NULL)
  nm <- trimws(sub("\\s*\\(.*$", "", hits)); sq <- sub("^.*\\(([A-Z]+)\\).*$", "\\1", hits)
  key <- paste(nm, sq); tb <- table(key)
  u <- !duplicated(key)
  data.frame(Name = nm[u], Sequence = sq[u], Records = as.integer(tb[key[u]]), Source = "BOLD records", stringsAsFactors = FALSE)
}
# a small reference list of widely used primers (check the original publication before ordering)
pr_primer_library <- function() {
  L <- function(name, seq, dir, marker, taxa, ref) data.frame(Name = name, Sequence = seq, Direction = dir, Marker = marker, Taxa = taxa, Reference = ref, stringsAsFactors = FALSE)
  rbind(
    L("LCO1490", "GGTCAACAAATCATAAAGATATTGG", "F", "COI", "Metazoa (invertebrates)", "Folmer et al. 1994"),
    L("HCO2198", "TAAACTTCAGGGTGACCAAAAAATCA", "R", "COI", "Metazoa (invertebrates)", "Folmer et al. 1994"),
    L("dgLCO1490", "GGTCAACAAATCATAAAGAYATYGG", "F", "COI", "Metazoa (invertebrates)", "Meyer 2003"),
    L("dgHCO2198", "TAAACTTCAGGGTGACCAAARAAYCA", "R", "COI", "Metazoa (invertebrates)", "Meyer 2003"),
    L("LepF1", "ATTCAACCAATCATAAAGATATTGG", "F", "COI", "Lepidoptera", "Hebert et al. 2004"),
    L("LepR1", "TAAACTTCTGGATGTCCAAAAAATCA", "R", "COI", "Lepidoptera", "Hebert et al. 2004"),
    L("mlCOIintF", "GGWACWGGWTGAACWGTWTAYCCYCC", "F", "COI", "Metazoa (metabarcoding)", "Leray et al. 2013"),
    L("BF3", "CCHGAYATRGCHTTYCCHCG", "F", "COI", "Arthropoda / freshwater invertebrates", "Elbrecht & Leese 2017"),
    L("BR2", "TCDGGRTGNCCRAARAAYCA", "R", "COI", "Arthropoda / freshwater invertebrates", "Elbrecht & Leese 2017"),
    L("FishF1", "TCAACCAACCACAAAGACATTGGCAC", "F", "COI", "Fish", "Ward et al. 2005"),
    L("FishR1", "TAGACTTCTGGGTGGCCAAAGAATCA", "R", "COI", "Fish", "Ward et al. 2005"),
    L("FishF2", "TCGACTAATCATAAAGATATCGGCAC", "F", "COI", "Fish", "Ward et al. 2005"),
    L("FishR2", "ACTTCAGGGTGACCGAAGAATCAGAA", "R", "COI", "Fish", "Ward et al. 2005"),
    L("MiFish-U-F", "GTCGGTAAAACTCGTGCCAGC", "F", "12S", "Fish (eDNA)", "Miya et al. 2015"),
    L("MiFish-U-R", "CATAGTGGGGTATCTAATCCCAGTTTG", "R", "12S", "Fish (eDNA)", "Miya et al. 2015"),
    L("16Sar-L", "CGCCTGTTTATCAAAAACAT", "F", "16S", "Metazoa", "Palumbi et al. 1991"),
    L("16Sbr-H", "CCGGTCTGAACTCAGATCACGT", "R", "16S", "Metazoa", "Palumbi et al. 1991"),
    L("L14841", "AAAAAGCTTCCATCCAACATCTCAGCATGATGAAA", "F", "Cytb", "Vertebrates", "Kocher et al. 1989"),
    L("H15149", "AAACTGCAGCCCCTCAGAATGATATTTGTCCTCA", "R", "Cytb", "Vertebrates", "Kocher et al. 1989"),
    L("rbcLa-F", "ATGTCACCACAAACAGAGACTAAAGC", "F", "rbcL", "Plants", "CBOL Plant Working Group 2009"),
    L("rbcLa-R", "GTAAAATCAAGTCCACCRCG", "R", "rbcL", "Plants", "CBOL Plant Working Group 2009"),
    L("MatK_390f", "CGATCTATTCATTCAATATTTC", "F", "matK", "Plants", "Cuenoud et al. 2002"),
    L("MatK_1326r", "TCTAGCACACGAAAGTCGAAGT", "R", "matK", "Plants", "Cuenoud et al. 2002"),
    L("psbA3_f", "GTTATGCATGAACGTAATGCTC", "F", "trnH-psbA", "Plants", "Sang et al. 1997"),
    L("trnHf_05", "CGCGCATGGTGGATTCACAATCC", "R", "trnH-psbA", "Plants", "Tate & Simpson 2003"),
    L("ITS1", "TCCGTAGGTGAACCTGCGG", "F", "ITS", "Fungi / plants", "White et al. 1990"),
    L("ITS4", "TCCTCCGCTTATTGATATGC", "R", "ITS", "Fungi / plants", "White et al. 1990"),
    L("ITS5", "GGAAGTAAAAGTCGTAACAAGG", "F", "ITS", "Fungi / plants", "White et al. 1990"),
    L("ITS2", "GCTGCGTTCTTCATCGATGC", "R", "ITS", "Fungi / plants", "White et al. 1990"),
    L("ITS1F", "CTTGGTCATTTAGAGGAAGTAA", "F", "ITS", "Fungi", "Gardes & Bruns 1993"),
    L("27F", "AGAGTTTGATCMTGGCTCAG", "F", "16S rRNA", "Bacteria", "Lane 1991"),
    L("1492R", "GGTTACCTTGTTACGACTT", "R", "16S rRNA", "Bacteria", "Lane 1991"),
    L("515F", "GTGYCAGCMGCCGCGGTAA", "F", "16S rRNA (V4)", "Bacteria / Archaea", "Parada et al. 2016"),
    L("806R", "GGACTACNVGGGTWTCTAAT", "R", "16S rRNA (V4)", "Bacteria / Archaea", "Apprill et al. 2015"),
    L("341F", "CCTACGGGNGGCWGCAG", "F", "16S rRNA (V3-V4)", "Bacteria", "Klindworth et al. 2013"),
    L("805R", "GACTACHVGGGTATCTAATCC", "R", "16S rRNA (V3-V4)", "Bacteria", "Klindworth et al. 2013"),
    L("M13F(-20)", "GTAAAACGACGGCCAGT", "F", "sequencing (vector)", "universal", "Messing 1983"),
    L("M13R", "CAGGAAACAGCTATGAC", "R", "sequencing (vector)", "universal", "Messing 1983"),
    L("T7", "TAATACGACTCACTATAGGG", "F", "sequencing (vector)", "universal", "Studier & Moffatt 1986"),
    L("T3", "ATTAACCCTCACTAAAGGGA", "F", "sequencing (vector)", "universal", "-"),
    L("SP6", "ATTTAGGTGACACTATAG", "F", "sequencing (vector)", "universal", "-"))
}

# ---- primers typed by the user -----------------------------------------------------------------------------------------------
pr_parse_primers <- function(txt) {
  ln <- trimws(unlist(strsplit(paste(txt, collapse = "\n"), "\r?\n"))); ln <- ln[nzchar(ln)]
  if (!length(ln)) return(NULL)
  if (any(startsWith(ln, ">"))) {
    f <- pr_parse_fasta(ln); sq <- .pr_clean(f)
    return(data.frame(Name = names(f), Sequence = sq, stringsAsFactors = FALSE)[nzchar(sq), , drop = FALSE])
  }
  nm <- character(length(ln)); sq <- character(length(ln))
  for (i in seq_along(ln)) {
    tok <- strsplit(ln[i], "[ \t,;]+")[[1]]
    sq[i] <- .pr_clean(tail(tok, 1L))
    nm[i] <- if (length(tok) > 1L) paste(head(tok, -1L), collapse = "_") else sprintf("P%d", i)
  }
  d <- data.frame(Name = make.unique(nm), Sequence = sq, stringsAsFactors = FALSE)
  d[nchar(d$Sequence) >= 8L, , drop = FALSE]
}
pr_analyse <- function(prim, prm = pr_params()) {
  g0 <- function(x, f) if (is.null(x)) 0 else x[[f]]
  isp <- grepl("probe|taqman|(^|[_-])P$", prim$Name, ignore.case = TRUE)
  ref <- if (any(!isp)) stats::median(vapply(prim$Sequence[!isp], function(z) pr_tm(.pr_clean(z), prm)[["mean"]], 0)) else NA_real_
  rows <- lapply(seq_len(nrow(prim)), function(i) {
    e <- pr_eval(prim$Sequence[i], prm)
    st <- pr_struct_primer(e$seq, prm)                       # structures with mismatches, salt correction and drawings
    e$hp <- g0(st$hp, "dg"); e$hp3 <- g0(st$hp3, "dg"); e$sd <- g0(st$sd, "dg"); e$sd3 <- g0(st$sd3, "dg")
    e$hp_tm <- max(0, g0(st$hp, "tm")); e$sd_tm <- max(0, g0(st$sd, "tm"))
    e$n3 <- .pr_dimer(st$seq, st$seq)$n3
    ch <- strsplit(e$seq, "")[[1]]
    e$gc2 <- sum(ch[max(1L, e$len - 1L):e$len] %in% c("G", "C", "S"))
    v <- pr_verdict(e, prm, probe = isp[i], ref_tm = ref)
    data.frame(Name = prim$Name[i], Role = if (isp[i]) "Probe" else "Primer", Sequence = e$seq, Length = e$len, `GC %` = round(e$gc, 1), `Tm (C)` = round(e$tm, 1),
               `Tm range` = if (e$deg > 1) sprintf("%.1f-%.1f", e$tm_min, e$tm_max) else "",
               `Ambiguous bases` = sum(!ch %in% c("A", "C", "G", "T")), `Degeneracy (sequences)` = e$deg,
               `GC clamp (last 2)` = e$gc2, `GC in last 5` = e$gc3, `Longest run` = e$run, `3' stability` = round(e$stab3, 1),
               `Hairpin dG` = round(e$hp, 1), `Hairpin Tm` = round(e$hp_tm, 1), `Self-dimer dG` = round(e$sd, 1), `Self-dimer Tm` = round(e$sd_tm, 1),
               `3' self-dimer dG` = round(e$sd3, 1), `3' pairs (bp)` = e$n3,
               Verdict = v$status, Notes = v$notes, check.names = FALSE, stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}
pr_pair_table <- function(prim, prm = pr_params()) {
  n <- nrow(prim); if (n < 2L) return(NULL)
  rows <- list(); k <- 0L
  full <- n * (n - 1L) / 2L <= 36L                            # many pairs: fast estimate only
  for (i in 1:(n - 1L)) for (j in (i + 1L):n) {
    a <- .pr_variants(prim$Sequence[i], 1L); b <- .pr_variants(prim$Sequence[j], 1L)
    d <- .pr_dimer(a, b); gany <- d$any; gend <- d$end; tmx <- 0
    if (full) {
      sa <- .pr_dup_struct(a, b, prm, FALSE, FALSE); se <- .pr_dup_struct(a, b, prm, TRUE, FALSE)
      gany <- if (is.null(sa)) 0 else sa$dg; gend <- if (is.null(se)) 0 else se$dg
      tmx <- max(0, if (is.null(sa)) 0 else sa$tm, if (is.null(se)) 0 else se$tm)
    }
    k <- k + 1L
    rows[[k]] <- data.frame(`Primer 1` = prim$Name[i], `Primer 2` = prim$Name[j],
                            `Delta Tm` = round(abs(pr_tm(prim$Sequence[i], prm)[["mean"]] - pr_tm(prim$Sequence[j], prm)[["mean"]]), 1),
                            `Hetero-dimer dG` = round(gany, 1), `3' hetero-dimer dG` = round(gend, 1), `3' pairs (bp)` = d$n3,
                            Verdict = if (gany < -9 || gend < -5) "FAIL" else if (gany < -6 || gend < -3 || d$n3 > prm$n3_max || tmx > prm$struct_tm_max) "CHECK" else "OK",
                            check.names = FALSE, stringsAsFactors = FALSE)
  }
  do.call(rbind, rows)
}

# ---- binding of a primer in a set of sequences (both strands, IUPAC aware) ----------------------------------------------------
pr_bind_scan <- function(primer, seqs, max_n = 300L) {
  primer <- .pr_clean(primer); L <- nchar(primer)
  seqs <- seqs[seq_len(min(length(seqs), max_n))]
  pF <- .pr_masks(primer); pR <- .pr_masks(.pr_rc(primer))
  rows <- vector("list", length(seqs))
  for (k in seq_along(seqs)) {
    s <- .pr_clean(gsub("[-. ]", "", seqs[[k]])); n <- nchar(s)
    nm <- if (!is.null(names(seqs))) names(seqs)[k] else as.character(k)
    if (n < L) { rows[[k]] <- data.frame(Sequence = nm, Strand = NA_character_, Pos = NA_integer_, Mismatches = NA_integer_, End3 = NA_integer_, stringsAsFactors = FALSE); next }
    sm <- .pr_masks(s); best <- NULL
    for (st in c("+", "-")) {
      pm <- if (st == "+") pF else pR
      mmv <- vapply(seq_len(n - L + 1L), function(i) sum(bitwAnd(sm[i:(i + L - 1L)], pm) == 0L), 0L)
      i <- which.min(mmv)
      i3 <- if (st == "+") (L - 2L):L else 1L:3L
      e3 <- sum(bitwAnd(sm[(i - 1L) + i3], pm[i3]) == 0L)
      if (is.null(best) || mmv[i] < best$Mismatches) best <- data.frame(Sequence = nm, Strand = st, Pos = i, Mismatches = mmv[i], End3 = e3, stringsAsFactors = FALSE)
    }
    rows[[k]] <- best
  }
  do.call(rbind, rows)
}
pr_bind_summary <- function(sc) {
  ok <- sc[!is.na(sc$Mismatches), , drop = FALSE]; n <- nrow(ok)
  if (!n) return(c(`Sequences` = 0, `0 mismatches` = 0, `<= 1` = 0, `<= 2` = 0, `<= 3` = 0, `3' end intact (<= 2 mm)` = 0))
  c(`Sequences` = n, `0 mismatches` = sum(ok$Mismatches == 0L), `<= 1` = sum(ok$Mismatches <= 1L), `<= 2` = sum(ok$Mismatches <= 2L),
    `<= 3` = sum(ok$Mismatches <= 3L), `3' end intact (<= 2 mm)` = sum(ok$Mismatches <= 2L & ok$End3 == 0L))
}

# ---- NCBI web BLAST for primers (short-query settings) + predicted amplicons ---------------------------------------------------
.pr_blast_url <- "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
.pr_blast_get <- function(query, timeout = 60, tries = 2L) {            # text answer of Blast.cgi (GET), with the fallback transport
  r <- try(.pr_net("GET", .pr_blast_url, query = query, httr::timeout(timeout), tries = tries), silent = TRUE)
  if (!inherits(r, "try-error")) {
    if (httr::status_code(r) >= 400L) stop(sprintf("HTTP %d from NCBI BLAST", httr::status_code(r)), call. = FALSE)
    return(.pr_text(r))
  }
  x <- .pr_fallback_get(.pr_blast_url, query, timeout)
  if (is.null(x)) stop(paste0(conditionMessage(attr(r, "condition")), if (nzchar(.pr_nz0(.pr_state$fallback_err))) paste0(" (system download: ", .pr_state$fallback_err, ")") else ""), call. = FALSE)
  x
}
.pr_nz0 <- function(x) if (is.null(x) || !length(x) || is.na(x[1])) "" else as.character(x[1])
# 1. send the primers: returns the job (RID ...) at once, the search runs at NCBI
pr_blast_submit <- function(queries, db = "core_nt", organism = "", hits = 100L, expect = 1000, word = 7L, email = "") {
  if (!requireNamespace("httr", quietly = TRUE)) stop("Install httr.")
  nm <- gsub("[^A-Za-z0-9_.-]", "_", names(queries))
  fa <- paste0(">", nm, "\n", unname(queries), collapse = "\n")
  # BLAST_PROGRAMS = "blastn" is essential: without it NCBI runs megablast (word size >= 16), which can never hit a 20-mer primer
  body <- list(CMD = "Put", PROGRAM = "blastn", BLAST_PROGRAMS = "blastn", DATABASE = db, QUERY = fa, EXPECT = expect, WORD_SIZE = word,
               NUCL_REWARD = 1, NUCL_PENALTY = -3, GAPCOSTS = "5 2", FILTER = "F", HITLIST_SIZE = hits, TOOL = "SangeRly")
  if (nzchar(trimws(email))) body$EMAIL <- trimws(email)
  org <- trimws(organism)
  if (nzchar(org)) body$ENTREZ_QUERY <- if (grepl("\\[", org)) org else paste0(org, "[Organism]")
  .pr_state$log <- NULL; .pr_log("submit: ", length(queries), " primer(s), db=", db, ", word=", word)
  r <- try(.pr_net("POST", .pr_blast_url, body = body, encode = "form", httr::timeout(25), tries = 2L), silent = TRUE)
  if (inherits(r, "try-error")) .pr_log("POST failed: ", conditionMessage(attr(r, "condition"))) else .pr_log("POST HTTP ", httr::status_code(r))
  txt <- if (!inherits(r, "try-error") && httr::status_code(r) < 400L) .pr_text(r) else {
    .pr_log("retrying the same request as GET"); tg <- try(.pr_blast_get(body, 30, 1L), silent = TRUE)
    if (inherits(tg, "try-error")) { .pr_log("GET failed: ", conditionMessage(attr(tg, "condition"))); stop(conditionMessage(attr(tg, "condition")), call. = FALSE) }
    .pr_log("GET ok"); tg }
  rid <- sub(".*RID = ([A-Z0-9-]+).*", "\\1", regmatches(txt, regexpr("RID = [A-Z0-9-]+", txt)))
  if (!length(rid) || !nzchar(rid)) {
    why <- trimws(gsub("\\s+", " ", gsub("<[^>]*>", " ", substr(txt, 1L, 1500L))))
    stop("NCBI BLAST did not accept the request (no RID returned). ", if (nzchar(why)) paste0("NCBI answered: ", substr(why, 1L, 220L), " ") else "", "Try again in a few minutes.", call. = FALSE)
  }
  .pr_log("RID ", if (length(rid)) rid else "(none)")
  rtoe <- suppressWarnings(as.integer(sub(".*RTOE = ([0-9]+).*", "\\1", regmatches(txt, regexpr("RTOE = [0-9]+", txt)))))
  list(rid = rid, rtoe = if (length(rtoe) && !is.na(rtoe)) rtoe else 10L, nm = nm, hits = hits, t0 = Sys.time())
}
# 2. ask NCBI whether it is finished
pr_blast_status <- function(rid) {
  s <- .pr_blast_get(list(CMD = "Get", FORMAT_OBJECT = "SearchInfo", RID = rid), 30, 2L)
  st <- regmatches(s, regexpr("Status=[A-Z]+", s)); th <- regmatches(s, regexpr("ThereAreHits=[a-z]+", s))
  list(status = if (length(st)) sub("Status=", "", st) else "WAITING", there = if (length(th)) grepl("=yes", th) else NA)
}
# 3. download and read the hit table
pr_blast_fetch <- function(rid, nm, hits = 100L) {
  tb <- .pr_blast_get(list(CMD = "Get", FORMAT_TYPE = "Tabular", ALIGNMENT_VIEW = "Tabular", RID = rid, ALIGNMENTS = hits, DESCRIPTIONS = hits), 120, 3L)
  ln <- strsplit(tb, "\r?\n")[[1]]
  keep <- character(); lab <- character(); ql <- ""
  for (l in ln) {                                   # data lines + the "# Query: ..." comment that precedes each block
    if (startsWith(l, "# Query:")) { ql <- trimws(sub("^# Query:\\s*", "", l)); next }
    if (nzchar(l) && !startsWith(l, "#") && grepl("\t", l, fixed = TRUE)) { keep <- c(keep, l); lab <- c(lab, ql) }
  }
  if (!length(keep)) return(list(hits = NULL, rid = rid))
  f <- strsplit(keep, "\t", fixed = TRUE)
  ok <- vapply(f, length, 1L) >= 12L; f <- f[ok]; lab <- lab[ok]
  if (!length(f)) return(list(hits = NULL, rid = rid))
  idc <- function(x) vapply(strsplit(x, "\\|"), function(p) { p <- p[nzchar(p)]; if (length(p)) tail(p, 1L) else "" }, "")
  q1 <- idc(vapply(f, `[`, "", 1L))
  if (!all(q1 %in% nm)) {                           # NCBI may rename queries to Query_<n>: recover the primer name
    byl <- vapply(lab, function(z) { tk <- strsplit(z, "[ |]+")[[1]]; k <- which(nm %in% tk); if (length(k)) nm[k[1]] else NA_character_ }, "", USE.NAMES = FALSE)
    if (!anyNA(byl)) q1 <- byl else {
      ids <- unique(q1); num <- suppressWarnings(as.numeric(sub("^\\D*", "", ids)))
      if (length(ids) == length(nm) && !anyNA(num)) q1 <- nm[match(q1, ids[order(num)])]
    }
  }
  d <- data.frame(Query = q1, Subject = idc(vapply(f, `[`, "", 2L)),
                  Identity = as.numeric(vapply(f, `[`, "", 3L)), Length = as.integer(vapply(f, `[`, "", 4L)),
                  Mismatches = as.integer(vapply(f, `[`, "", 5L)), Gaps = as.integer(vapply(f, `[`, "", 6L)),
                  Qstart = as.integer(vapply(f, `[`, "", 7L)), Qend = as.integer(vapply(f, `[`, "", 8L)),
                  Sstart = as.integer(vapply(f, `[`, "", 9L)), Send = as.integer(vapply(f, `[`, "", 10L)),
                  Evalue = as.numeric(vapply(f, `[`, "", 11L)), Bits = as.numeric(vapply(f, `[`, "", 12L)), stringsAsFactors = FALSE)
  list(hits = d, rid = rid)
}
# all in one (blocking) - the app itself uses the three steps above so that it stays responsive and BLAST can be cancelled
pr_blast_run <- function(queries, db = "core_nt", organism = "", hits = 100L, expect = 1000, word = 7L, email = "",
                         timeout = 600, progress = function(m) NULL) {
  j <- pr_blast_submit(queries, db, organism, hits, expect, word, email)
  Sys.sleep(min(j$rtoe, 20L))
  repeat {
    el <- as.numeric(difftime(Sys.time(), j$t0, units = "secs"))
    if (el > timeout) stop("BLAST took longer than ", timeout, " s. RID ", j$rid, " can still be fetched at blast.ncbi.nlm.nih.gov.")
    progress(sprintf("BLAST: waiting for NCBI (%d s, RID %s) ...", as.integer(el), j$rid))
    st <- try(pr_blast_status(j$rid), silent = TRUE)
    if (!inherits(st, "try-error")) {
      if (st$status == "READY") break
      if (st$status %in% c("FAILED", "UNKNOWN")) stop("NCBI BLAST reported status ", st$status, " for RID ", j$rid, ".")
    }
    Sys.sleep(15)
  }
  if (isFALSE(st$there)) return(list(hits = NULL, rid = j$rid))
  pr_blast_fetch(j$rid, j$nm, hits)
}
# titles of subject accessions (best effort)
pr_ncbi_titles <- function(acc, email = "") {
  acc <- unique(acc[nzchar(acc)]); if (!length(acc)) return(character())
  acc <- acc[seq_len(min(length(acc), 150L))]
  com <- list(tool = "SangeRly"); if (nzchar(trimws(email))) com$email <- trimws(email)
  r <- try(.pr_http("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi",
                    c(list(db = "nuccore", id = paste(acc, collapse = ","), retmode = "json"), com), 60), silent = TRUE)
  if (inherits(r, "try-error")) return(character())
  x <- try(jsonlite::fromJSON(.pr_text(r), simplifyVector = FALSE)$result, silent = TRUE)
  if (inherits(x, "try-error") || is.null(x)) return(character())
  out <- character()
  for (u in x$uids) { z <- x[[u]]; if (!is.null(z$accessionversion)) out[z$accessionversion] <- as.character(z$title) }
  out
}
# a primer pair on the same subject, opposite strands, product <= max_prod; each primer: <= max_mm differences and the 3' end aligned
pr_amplicons <- function(hits, f_name, r_name, f_len, r_len, max_mm = 3L, max_prod = 3000L) {
  if (is.null(hits) || !nrow(hits)) return(NULL)
  prep <- function(nm, len) {
    h <- hits[hits$Query == gsub("[^A-Za-z0-9_.-]", "_", nm), , drop = FALSE]
    h <- h[h$Qend == len & (h$Mismatches + h$Gaps + (len - h$Length)) <= max_mm, , drop = FALSE]
    if (!nrow(h)) return(h)
    h$plus <- h$Send > h$Sstart
    h$five <- ifelse(h$plus, h$Sstart - (h$Qstart - 1L), h$Sstart + (h$Qstart - 1L))
    h$mm <- h$Mismatches + h$Gaps + (len - h$Length)
    h
  }
  Fh <- prep(f_name, f_len); Rh <- prep(r_name, r_len)
  if (!nrow(Fh) || !nrow(Rh)) return(NULL)
  rows <- list(); k <- 0L
  for (s in intersect(Fh$Subject, Rh$Subject)) {
    a <- Fh[Fh$Subject == s, , drop = FALSE]; b <- Rh[Rh$Subject == s, , drop = FALSE]
    for (i in seq_len(nrow(a))) for (j in seq_len(nrow(b))) {
      if (a$plus[i] == b$plus[j]) next
      p <- if (a$plus[i]) a[i, ] else b[j, ]; m <- if (a$plus[i]) b[j, ] else a[i, ]
      size <- m$five - p$five + 1L
      if (size < 1L || size > max_prod) next
      k <- k + 1L
      rows[[k]] <- data.frame(Subject = s, Product = size, `Forward mismatches` = a$mm[i], `Reverse mismatches` = b$mm[j],
                              Orientation = if (a$plus[i]) "F+ / R-" else "F- / R+", check.names = FALSE, stringsAsFactors = FALSE)
    }
  }
  if (!k) return(NULL)
  d <- do.call(rbind, rows); d <- d[order(d$`Forward mismatches` + d$`Reverse mismatches`, d$Product), , drop = FALSE]
  rownames(d) <- NULL; d
}

# ---- GenBank flat file (annotated record): features, windows, spliced sequences ---------------------------------------------
# location -> parts (start, end) in transcript order + strand; join(), complement(), <, > are understood
.pr_loc_parts <- function(loc) {
  loc <- gsub("\\s+", "", loc)
  minus <- grepl("^complement\\(", loc)
  m <- regmatches(loc, gregexpr("[0-9]+(\\.\\.[0-9]+)?", gsub("[<>]", "", loc)))[[1]]
  if (!length(m)) return(NULL)
  st <- as.integer(sub("\\.\\..*$", "", m))
  en <- as.integer(ifelse(grepl("\\.\\.", m), sub("^.*\\.\\.", "", m), m))
  p <- data.frame(start = pmin(st, en), end = pmax(st, en))
  if (minus) p <- p[rev(seq_len(nrow(p))), , drop = FALSE]
  list(parts = p, strand = if (minus) "-" else "+")
}
# first record of a GenBank flat file -> list(accession, definition, organism, seq, features, parts)
pr_gb_parse <- function(txt) {
  ln <- unlist(strsplit(paste(txt, collapse = "\n"), "\r?\n")); ln <- ln[!is.na(ln)]
  e <- which(ln == "//"); if (length(e)) ln <- ln[seq_len(e[1])]
  fi <- grep("^FEATURES", ln); oi <- grep("^ORIGIN", ln)
  if (!length(oi) || oi[1] >= length(ln)) stop("This is not a GenBank record (no ORIGIN / sequence section).")
  hdr <- ln[seq_len((if (length(fi)) fi[1] else oi[1]) - 1L)]
  v <- hdr[grep("^VERSION", hdr)]
  acc <- if (length(v)) sub("^VERSION\\s+(\\S+).*", "\\1", v[1]) else sub("^LOCUS\\s+(\\S+).*", "\\1", hdr[1])
  di <- grep("^DEFINITION", hdr); def <- ""
  if (length(di)) {
    k <- di[1]; def <- sub("^DEFINITION\\s+", "", hdr[k]); k <- k + 1L
    while (k <= length(hdr) && grepl("^ {12}", hdr[k])) { def <- paste(def, trimws(hdr[k])); k <- k + 1L }
  }
  og <- hdr[grep("^\\s*ORGANISM", hdr)]; org <- if (length(og)) sub("^\\s*ORGANISM\\s+", "", og[1]) else ""
  sq <- toupper(gsub("[^A-Za-z]", "", paste(ln[(oi[1] + 1L):length(ln)], collapse = "")))
  fl <- if (length(fi)) ln[(fi[1] + 1L):(oi[1] - 1L)] else character()
  feat <- list(); cur <- NULL
  for (l in fl) {
    if (grepl("^ {5}[^ ]", l)) {
      if (!is.null(cur)) feat[[length(feat) + 1L]] <- cur
      cur <- list(type = trimws(substr(l, 6L, 20L)), loc = trimws(substr(l, 22L, nchar(l))), q = list(), last = NULL)
    } else if (!is.null(cur)) {
      t <- trimws(substr(l, 22L, nchar(l)))
      if (startsWith(t, "/")) {
        nm <- sub("^/([^=]+).*$", "\\1", t)
        val <- if (grepl("=", t, fixed = TRUE)) sub("^/[^=]+=", "", t) else "TRUE"
        if (is.null(cur$q[[nm]])) cur$q[[nm]] <- val
        cur$last <- nm
      } else if (is.null(cur$last)) {
        cur$loc <- paste0(cur$loc, t)
      } else {
        cur$q[[cur$last]] <- paste(cur$q[[cur$last]], t)
      }
    }
  }
  if (!is.null(cur)) feat[[length(feat) + 1L]] <- cur
  rows <- list(); parts <- list()
  for (f in feat) {
    if (f$type == "source") next
    lp <- .pr_loc_parts(f$loc); if (is.null(lp)) next
    q <- lapply(f$q, function(z) gsub("^\"|\"$", "", trimws(z)))
    nm <- q[["gene"]]; if (is.null(nm)) nm <- q[["locus_tag"]]; if (is.null(nm)) nm <- q[["product"]]; if (is.null(nm)) nm <- ""
    pr <- q[["product"]]; if (is.null(pr)) pr <- ""
    rows[[length(rows) + 1L]] <- data.frame(Type = f$type, Name = nm, Product = pr, Start = min(lp$parts$start), End = max(lp$parts$end),
                                           Strand = lp$strand, Length = sum(lp$parts$end - lp$parts$start + 1L), Parts = nrow(lp$parts),
                                           stringsAsFactors = FALSE)
    parts[[length(parts) + 1L]] <- lp$parts
  }
  feats <- if (length(rows)) do.call(rbind, rows) else data.frame(Type = character(), Name = character(), Product = character(), Start = integer(),
                                                                  End = integer(), Strand = character(), Length = integer(), Parts = integer(), stringsAsFactors = FALSE)
  list(accession = acc, definition = def, organism = org, seq = sq, features = feats, parts = parts)
}
pr_gb_fetch <- function(acc, email = "") {
  if (!requireNamespace("httr", quietly = TRUE)) stop("Install httr.")
  acc <- trimws(acc); if (!nzchar(acc)) stop("This record has no accession to download.")
  com <- list(tool = "SangeRly"); if (nzchar(trimws(email))) com$email <- trimws(email)
  txt <- ""; err <- NULL
  for (rt in c("gbwithparts", "gb")) {                      # a big genome can fail with one rettype: try the other
    r <- tryCatch(.pr_http("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
                           c(list(db = "nuccore", id = acc, rettype = rt, retmode = "text"), com), 240),
                  error = function(e) { err <<- conditionMessage(e); NULL })
    if (is.null(r)) next
    txt <- .pr_text(r)
    if (grepl("(^|\\n)ORIGIN", txt) && grepl("(^|\\n)FEATURES", txt)) { err <- NULL; break }
    err <- sprintf("NCBI returned no annotated GenBank record for %s.", acc)
  }
  if (!is.null(err) || !nzchar(txt)) stop(if (is.null(err)) "Empty answer from NCBI." else err, call. = FALSE)
  g <- pr_gb_parse(txt); g$raw <- txt
  g
}
# window [from, to] of the record, optionally read as the reverse complement; win$strand tells which
pr_gb_window <- function(g, from, to, strand = "+") {
  n <- nchar(g$seq); from <- max(1L, as.integer(from)); to <- min(n, as.integer(to))
  s <- substr(g$seq, from, to); if (strand == "-") s <- .pr_rc(s)
  list(seq = s, from = from, to = to, strand = strand)
}
pr_win_pos <- function(win, rec) if (win$strand == "-") win$to - rec + 1L else rec - win$from + 1L    # record -> window position
pr_rec_pos <- function(win, pos) if (win$strand == "-") win$to - pos + 1L else win$from + pos - 1L    # window -> record position
# spliced (cDNA-like) sequence of feature i, read in the direction of the gene, and the exon lengths
pr_gb_spliced <- function(g, i) {
  p <- g$parts[[i]]; f <- g$features[i, ]
  ex <- vapply(seq_len(nrow(p)), function(k) { s <- substr(g$seq, p$start[k], p$end[k]); if (f$Strand == "-") .pr_rc(s) else s }, "")
  list(name = paste0(g$accession, " ", f$Type, " ", if (nzchar(f$Name)) f$Name else f$Product), seq = paste(ex, collapse = ""), exons = nchar(ex))
}
# features overlapping the window, as drawn parts in window coordinates (lane 1 = gene / RNA, lane 2 = CDS / exon / mRNA)
pr_gb_win_features <- function(g, win, types = c("gene", "CDS", "exon", "mRNA", "rRNA", "tRNA")) {
  f <- g$features; if (!nrow(f)) return(NULL)
  idx <- which(f$Type %in% types & f$End >= win$from & f$Start <= win$to)
  rows <- list()
  for (i in idx) {
    p <- g$parts[[i]]; a <- pmax(p$start, win$from); b <- pmin(p$end, win$to); ok <- a <= b
    if (!any(ok)) next
    x <- pr_win_pos(win, a[ok]); y <- pr_win_pos(win, b[ok])
    rows[[length(rows) + 1L]] <- data.frame(start = pmin(x, y), end = pmax(x, y),
                                           label = if (nzchar(f$Name[i])) f$Name[i] else f$Product[i], type = f$Type[i],
                                           lane = if (f$Type[i] %in% c("gene", "rRNA", "tRNA")) 1L else 2L, stringsAsFactors = FALSE)
  }
  if (!length(rows)) NULL else do.call(rbind, rows)
}


# ---- secondary structures with drawings (Quality tab and the Primer3Plus-style report) ---------------------------------------
# Dimers: ungapped alignment of a with b (antiparallel), internal mismatches allowed (<= max_mm, each costs .pr_mm_dg); nearest-neighbour
# stacks (SantaLucia 1998), terminal initiation, salt correction of the reaction settings. All dG at 37 degC (kcal/mol); Tm at prm$conc.
.pr_mm_dg <- 1.7
.pr_end_par <- function(b) if (b %in% c("G", "C")) c(0.1, -2.8) else c(2.3, 4.1)       # dH, dS of a terminal base pair
.pr_naeq <- function(prm) (prm$na + 120 * sqrt(max(prm$mg - prm$dntp, 0))) / 1000
.pr_struct_tm <- function(dh, ds, prm, x = 4) dh * 1000 / (ds + 1.987 * log(prm$conc * 1e-9 / x)) - 273.15

.pr_art_dup <- function(a, b, st) {
  Bn <- paste(rev(strsplit(b, "")[[1]]), collapse = "")
  left <- min(1L, 1L + st$d); pt <- 1L - left; pb <- (1L + st$d) - left
  m <- rep(" ", nchar(a)); m[st$ii[st$mt]] <- "|"
  paste(c(paste0("5' ", strrep(" ", pt), a, " 3'"),
          paste0("   ", strrep(" ", pt), paste(m, collapse = "")),
          paste0("3' ", strrep(" ", pb), Bn, " 5'")), collapse = "\n")
}
# need3 = TRUE: only structures that contain a 3' terminal base of a or b
.pr_dup_struct <- function(a, b, prm = pr_params(), need3 = FALSE, self = FALSE, max_mm = 2L) {
  A <- strsplit(a, "")[[1]]; Bn <- rev(strsplit(b, "")[[1]])
  na <- length(A); nb <- length(Bn)
  if (na < 3L || nb < 3L) return(NULL)
  cmp <- unname(.pr_comp[A])
  mm_ds <- -.pr_mm_dg * 1000 / 310.15
  lnna <- log(.pr_naeq(prm))
  best <- NULL
  for (d in seq(-(nb - 1L), na - 1L)) {
    lo <- max(1L, 1L + d); hi <- min(na, nb + d)
    if (hi - lo + 1L < 3L) next
    ii <- lo:hi; jj <- ii - d
    mt <- cmp[ii] == Bn[jj]
    L <- length(ii)
    for (s in which(mt)) {
      if (s >= L) next
      sdh <- 0; sds <- 0; mm <- 0L
      for (e in (s + 1L):L) {
        v <- .pr_nn[[paste0(A[ii[e - 1L]], A[ii[e]])]]
        sdh <- sdh + v[1]; sds <- sds + v[2]
        if (!mt[e]) { mm <- mm + 1L; if (mm > max_mm) break; sds <- sds + mm_ds; next }
        if ((e - s + 1L) - mm < 3L) next
        if (need3 && !(ii[e] == na || jj[s] == 1L)) next
        t1 <- .pr_end_par(A[ii[s]]); t2 <- .pr_end_par(A[ii[e]])
        dh <- sdh + t1[1] + t2[1]
        ds <- sds + t1[2] + t2[2] + 0.368 * (e - s) * lnna
        dg <- dh - 310.15 * ds / 1000
        if (dg < 0 && (is.null(best) || dg < best$dg)) best <- list(dg = dg, dh = dh, ds = ds, d = d, ii = ii[s:e], mt = mt[s:e])
      }
    }
  }
  if (is.null(best)) return(NULL)
  best$tm <- .pr_struct_tm(best$dh, best$ds, prm, if (self) 1 else 4)
  best$art <- .pr_art_dup(a, b, best)
  best
}
.pr_art_hp <- function(h) {
  b <- strsplit(h$s, "")[[1]]; n <- h$n; i <- h$i; j <- h$j; k <- h$k
  up <- paste(b[1:(i + k - 1L)], collapse = "")
  low <- paste(b[n:(j - k + 1L)], collapse = "")                # 3' tail, read towards the stem
  c0 <- i - (n - j); pad <- max(0L, 1L - c0)
  paste(c(paste0("5' ", strrep(" ", pad), up, "\u2510"),
          paste0("   ", strrep(" ", pad + i - 1L), strrep("|", k)),
          paste0("3' ", strrep(" ", pad + c0 - 1L), low, "\u2518")), collapse = "\n")
}
# hairpins: perfect stem >= 3 bp, loop >= 3 nt; any = best, end = best one whose stem contains the 3' terminal base
.pr_hairpin_struct <- function(s, prm = pr_params()) {
  b <- strsplit(s, "")[[1]]; n <- length(b)
  if (n < 9L) return(list(any = NULL, end = NULL))
  cp <- unname(.pr_comp[b])
  best <- NULL; best3 <- NULL
  for (i in seq_len(n - 8L)) for (j in (i + 8L):n) {
    k <- 0L
    while (j - i - 2L * k >= 4L && cp[i + k] == b[j - k]) k <- k + 1L
    if (k < 3L) next
    for (kk in 3L:k) {
      stem <- b[i:(i + kk - 1L)]
      dh <- 0; ds <- 0
      for (u in seq_len(kk - 1L)) { v <- .pr_nn[[paste0(stem[u], stem[u + 1L])]]; dh <- dh + v[1]; ds <- ds + v[2] }
      loop <- j - i - 2L * kk + 1L                              # simplified: stem stacks + loop penalty (no end initiation)
      ds <- ds - .pr_loop_pen(loop) * 1000 / 310.15
      dg <- dh - 310.15 * ds / 1000
      if (dg >= 0) next
      cand <- list(dg = dg, dh = dh, ds = ds, tm = dh * 1000 / ds - 273.15, i = i, j = j, k = kk, n = n, s = s)
      if (is.null(best) || dg < best$dg) best <- cand
      if (j == n && (is.null(best3) || dg < best3$dg)) best3 <- cand
    }
  }
  if (!is.null(best)) best$art <- .pr_art_hp(best)
  if (!is.null(best3)) best3$art <- .pr_art_hp(best3)
  list(any = best, end = best3)
}
# all structures of one primer: hairpin (hp, hp3 = at the 3' end), self-dimer (sd, sd3 = at the 3' end); degenerate primers: first variant
pr_struct_primer <- function(s, prm = pr_params()) {
  s <- .pr_variants(.pr_clean(s), 1L)
  h <- .pr_hairpin_struct(s, prm)
  list(seq = s, hp = h$any, hp3 = h$end, sd = .pr_dup_struct(s, s, prm, FALSE, TRUE), sd3 = .pr_dup_struct(s, s, prm, TRUE, TRUE))
}
# cross structures of a primer pair
pr_struct_pair <- function(f, r, prm = pr_params()) {
  f <- .pr_variants(.pr_clean(f), 1L); r <- .pr_variants(.pr_clean(r), 1L)
  list(any = .pr_dup_struct(f, r, prm, FALSE, FALSE), end = .pr_dup_struct(f, r, prm, TRUE, FALSE))
}
# Tm of a PCR product from its sequence (GC / length / salt formula, as primer3 does)
pr_product_tm <- function(x, prm = pr_params()) {
  x <- .pr_clean(x); n <- nchar(x); if (n < 20L) return(NA_real_)
  81.5 + 16.6 * log10(.pr_naeq(prm)) + 41 * (.pr_ngc(x) / n) - 675 / n
}


# quick connection test (BLAST and Entrez servers, with both ways of connecting) for the BLAST tab
pr_net_test <- function() {
  urls <- c(`BLAST server (blast.ncbi.nlm.nih.gov)` = "https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get&FORMAT_OBJECT=SearchInfo&RID=TEST",
            `Entrez server (eutils.ncbi.nlm.nih.gov)` = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/einfo.fcgi?db=nuccore")
  do.call(rbind, lapply(names(urls), function(n) {
    u <- urls[[n]]
    a <- tryCatch({ r <- httr::GET(u, httr::timeout(20), httr::user_agent("SangeRly (test)"), .pr_cfg()); paste("OK (HTTP", httr::status_code(r), ")") },
                  error = function(e) paste("FAILED:", trimws(conditionMessage(e))))
    b <- tryCatch({ x <- .pr_fallback_get(u, list(), 30); if (is.null(x)) paste("FAILED:", .pr_nz0(.pr_state$fallback_err)) else "OK" },
                  error = function(e) paste("FAILED:", trimws(conditionMessage(e))))
    data.frame(Server = n, `httr / libcurl` = a, `System download (Windows connection)` = b, check.names = FALSE, stringsAsFactors = FALSE)
  }))
}


# ---- primer tables from other databases (PR2-primers, PrimerBank, qPrimerDB ... : saved / copied tables) --------------------------
# the table may have any column names: the columns that hold primer sequences are recognised by their content
.pr_primers_from_df <- function(df) {
  df[] <- lapply(df, function(x) trimws(as.character(x)))
  looks <- vapply(df, function(x) {
    x <- x[!is.na(x) & nzchar(x)]; if (!length(x)) return(FALSE)
    mean(grepl("^[ACGTURYSWKMBDHVN]{8,80}$", gsub("[ -]", "", toupper(x)))) >= 0.8
  }, TRUE)
  sc <- names(df)[looks]
  if (!length(sc)) stop("No column with primer sequences (A, C, G, T / IUPAC letters, 8-80 nt) was found in the table.")
  rest <- setdiff(names(df), sc)
  find <- function(pat) { k <- grep(pat, rest, ignore.case = TRUE, value = TRUE); if (length(k)) k[1] else NA_character_ }
  ncl <- find("name|label|primer|set|id"); dcl <- find("dir|orient|strand|sense"); tcl <- find("gene|target|marker|region|locus|taxon|species|organism")
  get <- function(cl) if (is.na(cl)) rep("", nrow(df)) else df[[cl]]
  out <- lapply(sc, function(cn) {
    nm <- get(ncl); nm <- ifelse(is.na(nm) | !nzchar(nm), paste0("primer", seq_len(nrow(df))), nm)
    if (length(sc) > 1L) nm <- paste0(nm, "_", cn)
    data.frame(Name = gsub("\\s+", "_", nm), Sequence = .pr_clean(gsub("[ -]", "", df[[cn]])),
               Direction = if (length(sc) > 1L) cn else get(dcl), Target = get(tcl), stringsAsFactors = FALSE)
  })
  d <- do.call(rbind, out); d <- d[nchar(d$Sequence) >= 8L, , drop = FALSE]
  if (!nrow(d)) stop("The table has no primer sequence.")
  d$Name <- make.unique(d$Name); rownames(d) <- NULL; d
}
.pr_primers_from_fasta <- function(path) {
  ln <- readLines(path, warn = FALSE); ln <- ln[nzchar(trimws(ln))]
  h <- which(startsWith(ln, ">")); if (!length(h)) stop("No FASTA record found in the file.")
  e <- c(h[-1] - 1L, length(ln))
  sq <- vapply(seq_along(h), function(i) if (e[i] > h[i]) paste(ln[(h[i] + 1L):e[i]], collapse = "") else "", "")
  hd <- sub("^>\\s*", "", ln[h]); f <- strsplit(hd, "\\s*\\|\\s*|\\s+")
  nm <- vapply(f, `[`, "", 1L)
  dr <- vapply(f, function(z) { k <- z[toupper(z) %in% c("F", "R", "FORWARD", "REVERSE")]; if (length(k)) toupper(substr(k[1], 1L, 1L)) else "" }, "")
  tg <- vapply(f, function(z) paste(z[-1L][!toupper(z[-1L]) %in% c("F", "R", "FORWARD", "REVERSE")], collapse = " "), "")
  df <- data.frame(Name = nm, Sequence = sq, Direction = dr, Target = tg, stringsAsFactors = FALSE)
  df
}
# BOLD primer database (portal.boldsystems.org/primer): the portal has no documented primer API, so a few likely addresses are tried;
# when they do not answer, the FASTA / table exported from the primer page can be loaded instead
pr_bold_primer_db <- function(query = "") {
  if (!requireNamespace("jsonlite", quietly = TRUE)) stop("Install jsonlite.")
  d <- NULL
  for (u in c("https://portal.boldsystems.org/api/primers", "https://portal.boldsystems.org/api/primer", "https://portal.boldsystems.org/api/primers/all")) {
    r <- try(.pr_http(u, list(), 40), silent = TRUE); if (inherits(r, "try-error")) next
    x <- try(jsonlite::fromJSON(.pr_text(r), simplifyVector = TRUE), silent = TRUE); if (inherits(x, "try-error")) next
    x <- if (is.data.frame(x)) x else if (is.list(x) && is.data.frame(x$data)) x$data else NULL
    if (!is.null(x) && nrow(x)) { d <- tryCatch(.pr_primers_from_df(as.data.frame(lapply(x, function(v) if (is.list(v)) vapply(v, function(z) paste(z, collapse = " "), "") else v), stringsAsFactors = FALSE)), error = function(e) NULL); if (!is.null(d)) break }
  }
  if (is.null(d)) stop("The BOLD portal did not answer a primer request from the app (it has no documented primer interface). Open portal.boldsystems.org/primer, search there, press the FASTA button and load the file in 'Load primers from another database' below.")
  q <- trimws(query)
  if (nzchar(q)) d <- d[grepl(q, apply(d, 1L, paste, collapse = " "), ignore.case = TRUE), , drop = FALSE]
  if (!nrow(d)) stop("No BOLD primer matches '", q, "'.")
  d$Source <- "BOLD primer database"; d
}
# the BOLD primer database (portal.boldsystems.org/primer), built into the app: 3,788 primers (Name, Sequence, Marker, Direction, Submitter) from the
# BOLD export; entries with the same name AND sequence were merged, the same name with another sequence was numbered _2, _3 ...
.pr_bold_z <- c(
  "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"
)
pr_bold_builtin <- function() {
  if (!is.null(.pr_state$boldb)) return(.pr_state$boldb)
  d <- tryCatch({
    r <- jsonlite::base64_dec(paste(.pr_bold_z, collapse = ""))
    utils::read.csv(text = rawToChar(memDecompress(r, type = "gzip")), stringsAsFactors = FALSE, na.strings = character())
  }, error = function(e) NULL)
  if (is.null(d) || !nrow(d)) return(NULL)
  .pr_state$boldb <- d; d
}
# PR2-primers (app.pr2-primers.org): the table exported on 2026-11-15 (primers, probes and blocking primers of the 18S / 28S / ITS / 16S / 23S rRNA genes), built into the app
.pr_pr2_z <- c(
  "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"
)
pr_pr2_builtin <- function() {
  if (!is.null(.pr_state$pr2b)) return(.pr_state$pr2b)
  d <- tryCatch({
    r <- jsonlite::base64_dec(paste(.pr_pr2_z, collapse = ""))
    utils::read.csv(text = rawToChar(memDecompress(r, type = "gzip")), stringsAsFactors = FALSE, na.strings = character())
  }, error = function(e) NULL)
  if (is.null(d) || !nrow(d)) return(NULL)
  .pr_state$pr2b <- d; d
}
pr_read_primer_table <- function(path, fname = path) {
  ext <- tolower(tools::file_ext(fname))
  if (!ext %in% c("xlsx", "xls") && startsWith(trimws(readLines(path, n = 1L, warn = FALSE)[1]), ">")) {
    d <- .pr_primers_from_fasta(path); d <- .pr_primers_from_df(d); d$Source <- basename(fname); return(d)
  }
  df <- if (ext %in% c("xlsx", "xls")) {
    if (!requireNamespace("readxl", quietly = TRUE)) stop("Reading an Excel file needs the 'readxl' package; install it or save the table as CSV.")
    as.data.frame(readxl::read_excel(path, sheet = 1), stringsAsFactors = FALSE)
  } else {
    first <- readLines(path, n = 1L, warn = FALSE)
    sep <- if (grepl("\t", first)) "\t" else if (grepl(";", first) && !grepl(",", first)) ";" else ","
    utils::read.table(path, header = TRUE, sep = sep, quote = "\"", comment.char = "", fill = TRUE, stringsAsFactors = FALSE,
                      check.names = FALSE, strip.white = TRUE)
  }
  if (!nrow(df)) stop("The table is empty.")
  d <- .pr_primers_from_df(df); d$Source <- basename(fname); d
}


# ---- transcripts: exon structure read from a GenBank record ---------------------------------------------------------------
# A transcript = list(id, name, gene, strand, seq (spliced, 5'->3'), exons (lengths), tstart / tend (exon limits on the transcript),
#                     gstart / gend (exon limits on the record, in transcript order), cds = c(first, last base on the transcript) or NULL)
.pr_g2t <- function(tx, pos) {                                   # record position -> transcript position (NA when not in an exon)
  for (k in seq_along(tx$exons)) if (pos >= tx$gstart[k] && pos <= tx$gend[k])
    return(if (tx$strand == "-") tx$tstart[k] + (tx$gend[k] - pos) else tx$tstart[k] + (pos - tx$gstart[k]))
  NA_integer_
}
.pr_tx_from_exons <- function(len, seq = NULL, name = "transcript") {
  len <- as.integer(len); ce <- cumsum(len)
  list(id = "x", name = name, gene = "", strand = "+", seq = seq, exons = len, tstart = ce - len + 1L, tend = ce, gstart = NULL, gend = NULL, cds = NULL)
}
.pr_tx_from_junc <- function(junc, n, name = "transcript", seq = NULL) {   # junc = last base of every exon but the last
  junc <- sort(unique(as.integer(junc))); junc <- junc[!is.na(junc) & junc >= 1L & junc < n]
  .pr_tx_from_exons(diff(c(0L, junc, as.integer(n))), seq, name)
}
# exons touched by the footprint [a, b]: "3" or "2/3" (sep = "-" for a range: "2-4")
.pr_exon_span <- function(tx, a, b, sep = "/") {
  k <- which(tx$tend >= a & tx$tstart <= b)
  if (!length(k)) "" else if (length(k) == 1L) as.character(k) else if (identical(sep, "-")) paste0(min(k), "-", max(k)) else paste(k, collapse = sep)
}
.pr_tx_build <- function(g, parts, strand, name, id, gene = "") {
  ex <- vapply(seq_len(nrow(parts)), function(k) { s <- substr(g$seq, parts$start[k], parts$end[k]); if (strand == "-") .pr_rc(s) else s }, "")
  len <- nchar(ex); ce <- cumsum(len)
  tx <- list(id = id, name = name, gene = gene, strand = strand, seq = paste(ex, collapse = ""), exons = len,
             tstart = ce - len + 1L, tend = ce, gstart = parts$start, gend = parts$end, cds = NULL)
  f <- g$features
  for (ci in which(f$Type == "CDS" & f$Strand == strand & f$Start >= min(parts$start) & f$End <= max(parts$end))) {   # coding region on the transcript
    cp <- g$parts[[ci]]
    a <- .pr_g2t(tx, if (strand == "-") cp$end[1] else cp$start[1]); b <- .pr_g2t(tx, if (strand == "-") cp$start[nrow(cp)] else cp$end[nrow(cp)])
    if (!is.na(a) && !is.na(b) && b > a) { tx$cds <- c(a, b); break }
  }
  tx
}
# transcripts that can be chosen: mRNA / CDS / RNA features with 2 or more exons, else the exon features (e.g. RefSeq mRNA records), grouped by gene
pr_gb_transcripts <- function(g) {
  f <- g$features; ids <- character(); lab <- character()
  if (nrow(f)) {
    ty <- c("mRNA", "ncRNA", "misc_RNA", "transcript", "CDS", "tRNA", "rRNA")
    i <- which(f$Type %in% ty & f$Parts >= 2L); i <- i[order(match(f$Type[i], ty), f$Start[i])]
    for (k in i) {
      ids <- c(ids, paste0("f:", k))
      lab <- c(lab, sprintf("%s [%s]: %d exons, %s nt", if (nzchar(f$Name[k])) f$Name[k] else f$Product[k], f$Type[k], f$Parts[k], format(f$Length[k], big.mark = ",")))
    }
    ex <- which(f$Type == "exon")
    if (length(ex) >= 2L) {
      nm <- ifelse(nzchar(f$Name[ex]), f$Name[ex], "(unnamed)")
      for (gn in unique(nm)) {
        e <- ex[nm == gn]; if (length(e) < 2L || gn %in% f$Name[i]) next
        ids <- c(ids, paste0("e:", gn)); lab <- c(lab, sprintf("%s [exon features]: %d exons", gn, length(e)))
      }
    }
  }
  data.frame(id = ids, label = lab, stringsAsFactors = FALSE)
}
pr_gb_transcript <- function(g, id) {
  f <- g$features; kind <- sub(":.*$", "", id); key <- sub("^[a-z]:", "", id)
  if (kind == "f") {
    i <- as.integer(key); if (is.na(i) || i < 1L || i > nrow(f)) stop("Unknown feature.")
    nm <- if (nzchar(f$Name[i])) f$Name[i] else f$Product[i]
    return(.pr_tx_build(g, g$parts[[i]], f$Strand[i], paste(g$accession, nm, f$Type[i]), id, f$Name[i]))
  }
  ex <- which(f$Type == "exon" & (f$Name == key | (key == "(unnamed)" & !nzchar(f$Name))))
  if (!length(ex)) stop("No exon feature for this gene.")
  strand <- f$Strand[ex[1]]
  parts <- do.call(rbind, lapply(ex, function(k) g$parts[[k]]))
  parts <- parts[!duplicated(paste(parts$start, parts$end)), , drop = FALSE]
  parts <- parts[order(parts$start), , drop = FALSE]
  if (strand == "-") parts <- parts[rev(seq_len(nrow(parts))), , drop = FALSE]
  .pr_tx_build(g, parts, strand, paste(g$accession, key, "exons"), id, key)
}
# the transcript that belongs to feature i: itself (mRNA / CDS / RNA), or for a gene its mRNA (else CDS) with the same strand inside it
pr_gb_feature_tx <- function(g, i) {
  f <- g$features
  if (f$Type[i] == "gene") {
    ch <- which(f$Type %in% c("mRNA", "CDS") & f$Strand == f$Strand[i] & f$Start >= f$Start[i] & f$End <= f$End[i])
    if (length(ch)) { ch <- ch[order(match(f$Type[ch], c("mRNA", "CDS")))]; i <- ch[1] }
  }
  pr_gb_transcript(g, paste0("f:", i))
}

# exon structure of a transcript, optionally with the primers (F, R, P = c(first, last) base on the transcript)
# Transcript map in the style of the NCBI Primer-BLAST graphic: exon track (numbered) on top, one row per designed pair
# (forward = blue arrow ->, reverse = red arrow <-, thin line = amplicon). pairs = data.frame with F_start/F_end/R_start/R_end (+ Probe, Probe_start);
# sel = highlighted row; zoom = TRUE shows only the selected amplicon (+ margin).
pr_plot_tx <- function(tx, F = NULL, R = NULL, P = NULL, pairs = NULL, sel = 1L, zoom = FALSE) {
  L <- sum(tx$exons); ne <- length(tx$exons)
  if (is.null(pairs) && !is.null(F) && !is.null(R)) {
    pairs <- data.frame(F_start = F[1], F_end = F[2], R_start = R[1], R_end = R[2], stringsAsFactors = FALSE)
    if (!is.null(P)) { pairs$Probe_start <- P[1]; pairs$Probe <- strrep("N", P[2] - P[1] + 1L) }
  }
  np <- if (is.null(pairs)) 0L else nrow(pairs)
  sel <- if (np) max(1L, min(as.integer(sel), np)) else 1L
  lo <- 0.5; hi <- L + 0.5
  if (np && zoom) {
    a <- min(pairs$F_start[sel], pairs$R_start[sel]); b <- max(pairs$F_end[sel], pairs$R_end[sel]); m <- max(40, 0.12 * (b - a))
    lo <- max(0.5, a - m); hi <- min(L + 0.5, b + m)
  }
  ex <- data.frame(k = seq_len(ne), xmin = tx$tstart - 0.5, xmax = tx$tend + 0.5, len = tx$exons, stringsAsFactors = FALSE)
  ex <- ex[ex$xmax >= lo & ex$xmin <= hi, , drop = FALSE]
  ex$fill <- ifelse(ex$k %% 2L == 1L, "Exon (odd)", "Exon (even)")
  ex$mid <- (pmax(ex$xmin, lo) + pmin(ex$xmax, hi)) / 2
  ex$lab <- sprintf("exon %d\n%d nt", ex$k, ex$len)
  cols <- c(`Exon (odd)` = "#5b6b7c", `Exon (even)` = "#2b3947", Forward = "#1f6fe0", Reverse = "#e8503a", Probe = "#22b07d", CDS = "#2e9d6a")
  span <- hi - lo
  g <- ggplot2::ggplot() +
    ggplot2::annotate("segment", x = lo, xend = hi, y = 0, yend = 0, colour = "grey60") +
    ggplot2::geom_rect(data = ex, ggplot2::aes(xmin = xmin, xmax = xmax, ymin = -0.2, ymax = 0.2, fill = fill), colour = "white", linewidth = 0.5) +
    ggplot2::geom_text(data = ex, ggplot2::aes(x = mid, y = 0, label = k), size = if (nrow(ex) > 15L) 3 else 4, colour = "white", fontface = "bold", check_overlap = TRUE)
  if (nrow(ex) <= 14L) g <- g + ggplot2::geom_text(data = ex, ggplot2::aes(x = mid, y = 0.5, label = lab), size = 3, colour = "grey25", lineheight = 0.9, check_overlap = TRUE)
  if (!is.null(tx$cds)) g <- g + ggplot2::annotate("rect", xmin = tx$cds[1] - 0.5, xmax = tx$cds[2] + 0.5, ymin = -0.32, ymax = -0.26, fill = cols[["CDS"]]) +
    ggplot2::annotate("text", x = max(lo, tx$cds[1]), y = -0.42, label = "CDS", size = 2.8, colour = cols[["CDS"]], hjust = 0)
  lab_of <- function(ab) .pr_exon_span(tx, ab[1], ab[2])
  ymin_all <- -0.5
  if (np) {
    pairs$row <- seq_len(np); pairs$y <- -0.85 - 0.55 * (pairs$row - 1L)
    pairs$sel <- pairs$row == sel
    pairs$prod <- pmax(pairs$F_end, pairs$R_end) - pmin(pairs$F_start, pairs$R_start) + 1L
    pairs <- pairs[pmax(pairs$F_end, pairs$R_end) >= lo & pmin(pairs$F_start, pairs$R_start) <= hi, , drop = FALSE]
    ymin_all <- min(pairs$y) - 0.45
    pairs$eF <- vapply(seq_len(nrow(pairs)), function(i) lab_of(c(pairs$F_start[i], pairs$F_end[i])), "")
    pairs$eR <- vapply(seq_len(nrow(pairs)), function(i) lab_of(c(pairs$R_start[i], pairs$R_end[i])), "")
    pairs$label <- sprintf("Pair %d  |  F: exon %s,  R: exon %s,  %d bp", pairs$row, pairs$eF, pairs$eR, pairs$prod)
    pairs$lx <- ifelse(pmin(pairs$F_start, pairs$R_start) - lo > 0.45 * span, pmin(pairs$F_start, pairs$R_start) - 0.01 * span, pmax(pairs$F_end, pairs$R_end) + 0.01 * span)
    pairs$hj <- ifelse(pmin(pairs$F_start, pairs$R_start) - lo > 0.45 * span, 1, 0)
    # a primer is drawn at least ~1% of the view wide so that it stays visible on long transcripts
    mw <- 0.01 * span
    ar <- grid::arrow(length = grid::unit(0.11, "inches"), type = "closed")
    pairs$Fx2 <- pmax(pairs$F_end, pairs$F_start + mw); pairs$Rx1 <- pmin(pairs$R_start, pairs$R_end - mw)
    g <- g +
      ggplot2::geom_segment(data = pairs, ggplot2::aes(x = F_end, xend = R_start, y = y, yend = y, linewidth = sel), colour = "grey55", show.legend = FALSE) +
      ggplot2::geom_segment(data = pairs, ggplot2::aes(x = F_start, xend = Fx2, y = y, yend = y, colour = "Forward", alpha = sel), linewidth = 2.6, arrow = ar, show.legend = FALSE) +
      ggplot2::geom_segment(data = pairs, ggplot2::aes(x = R_end, xend = Rx1, y = y, yend = y, colour = "Reverse", alpha = sel), linewidth = 2.6, arrow = ar, show.legend = FALSE) +
      ggplot2::geom_text(data = pairs, ggplot2::aes(x = lx, y = y + 0.2, label = label, hjust = hj, fontface = ifelse(sel, "bold", "plain")), size = 3.1, colour = "grey20") +
      ggplot2::scale_linewidth_manual(values = c(`TRUE` = 0.9, `FALSE` = 0.4), guide = "none") +
      ggplot2::scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.55), guide = "none") +
      ggplot2::scale_colour_manual(values = cols, guide = "none")
    ps <- pairs[pairs$sel & "Probe" %in% names(pairs) & !is.na(pairs$Probe), , drop = FALSE]
    if (nrow(ps)) g <- g + ggplot2::geom_rect(data = ps, ggplot2::aes(xmin = Probe_start - 0.5, xmax = Probe_start + nchar(Probe) - 0.5, ymin = y - 0.13, ymax = y + 0.13, fill = "Probe"), colour = "white")
  }
  g + ggplot2::scale_fill_manual(values = cols, breaks = c("Exon (odd)", "Probe")[c(TRUE, np > 0L && "Probe" %in% names(pairs) && any(!is.na(pairs$Probe)))], name = NULL) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.02), labels = function(b) format(b, big.mark = ",", trim = TRUE)) +
    ggplot2::scale_y_continuous(limits = c(ymin_all, 0.95)) +
    ggplot2::coord_cartesian(xlim = c(lo, hi)) +
    ggplot2::labs(x = "position on the transcript (nt)", y = NULL,
                  subtitle = sprintf("%s: %d exon%s, %s nt%s%s", tx$name, ne, if (ne == 1L) "" else "s", format(L, big.mark = ","),
                                     if (np) "; blue arrow = forward, red arrow = reverse, thin line = amplicon" else "",
                                     if (np && zoom) " (zoomed to the selected pair)" else "")) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(), panel.grid.major.y = ggplot2::element_blank(),
                   legend.position = "top")
}
# overview of the whole record (genes, CDS, mRNA, RNAs; + strand above, - strand below); sel = feature to highlight
pr_plot_record <- function(g, sel = NULL) {
  f <- g$features; n <- nchar(g$seq)
  tys <- c("gene", "mRNA", "CDS", "rRNA", "tRNA", "ncRNA"); if (!any(f$Type %in% c("mRNA", "CDS"))) tys <- c(tys, "exon")
  keep <- which(f$Type %in% tys); if (!length(keep)) return(NULL)
  rows <- lapply(keep, function(i) {
    p <- g$parts[[i]]
    data.frame(xmin = p$start - 0.5, xmax = p$end + 0.5, type = f$Type[i], y = if (f$Strand[i] == "-") -0.7 else 0.7, idx = i, stringsAsFactors = FALSE)
  })
  d <- do.call(rbind, rows)
  d$h <- c(gene = 0.07, mRNA = 0.16, exon = 0.2, CDS = 0.24, rRNA = 0.24, tRNA = 0.24, ncRNA = 0.24)[d$type]
  d <- d[order(match(d$type, tys)), , drop = FALSE]
  cols <- c(gene = "#9aa5b1", mRNA = "#7f8c8d", exon = "#16a085", CDS = "#2e9d6a", rRNA = "#b9770e", tRNA = "#8e44ad", ncRNA = "#8e44ad")
  gg <- ggplot2::ggplot() + ggplot2::annotate("segment", x = 0, xend = n, y = 0, yend = 0, colour = "grey70") +
    ggplot2::geom_rect(data = d, ggplot2::aes(xmin = xmin, xmax = xmax, ymin = y - h, ymax = y + h, fill = type))
  gn <- f[f$Type == "gene", , drop = FALSE]
  if (nrow(gn) && nrow(gn) <= 60L) {
    gn$y <- ifelse(gn$Strand == "-", -1.05, 1.05); gn$mid <- (gn$Start + gn$End) / 2
    gg <- gg + ggplot2::geom_text(data = gn, ggplot2::aes(x = mid, y = y, label = Name), size = 2.8, check_overlap = TRUE)
  }
  if (!is.null(sel) && !is.na(sel)) {
    p <- g$parts[[sel]]; y <- if (f$Strand[sel] == "-") -0.7 else 0.7
    gg <- gg + ggplot2::annotate("rect", xmin = p$start - 0.5, xmax = p$end + 0.5, ymin = y - 0.3, ymax = y + 0.3, fill = "#e8503a") +
      ggplot2::annotate("text", x = mean(c(min(p$start), max(p$end))), y = if (y > 0) 1.45 else -1.45, size = 3.6, colour = "#e8503a", fontface = "bold",
                        label = paste0(if (nzchar(f$Name[sel])) f$Name[sel] else f$Product[sel], " (", f$Type[sel], ")"))
  }
  gg + ggplot2::scale_fill_manual(values = cols, name = NULL) +
    ggplot2::scale_x_continuous(limits = c(0, n), labels = function(b) format(b, big.mark = ",", trim = TRUE)) +
    ggplot2::scale_y_continuous(limits = c(-1.7, 1.7)) +
    ggplot2::labs(x = "position on the record (nt)", y = NULL, subtitle = "+ strand above the line, - strand below; red = selected feature") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(), panel.grid.major.y = ggplot2::element_blank(), legend.position = "top")
}

# ---- alignment view (Jalview-like): sequences, IUPAC consensus, conservation bar, and the primers / probe on top ------------------------------
# aln: aligned sequences (named); from/to: columns shown; F, R, P: c(first, last) alignment column of the forward primer, reverse primer, probe.
# Bases equal to the most frequent base of the column are pale, different bases (the variation) are solid.
pr_plot_msa <- function(aln, from = 1L, to = NULL, F = NULL, R = NULL, P = NULL, maxrows = 25L, minfreq = 0.1) {
  aln <- toupper(gsub("U", "T", as.character(aln))); nm <- names(aln); if (is.null(nm)) nm <- paste0("seq", seq_along(aln))
  n <- length(aln); L <- nchar(aln[1])
  from <- max(1L, as.integer(from)); to <- min(L, if (is.null(to) || is.na(to)) from + 99L else as.integer(to)); if (to < from) to <- min(L, from + 99L)
  cols <- from:to; nc <- length(cols)
  M <- do.call(rbind, strsplit(substring(aln, from, to), ""))                 # n x nc
  bases <- c("A", "C", "G", "T")
  iu <- c(A = "A", C = "C", G = "G", T = "T", AG = "R", CT = "Y", CG = "S", AT = "W", GT = "K", AC = "M", CGT = "B", AGT = "D", ACT = "H", ACG = "V", ACGT = "N")
  st <- lapply(seq_len(nc), function(j) {
    x <- M[, j]; cnt <- vapply(bases, function(b) sum(x == b), 0); gap <- sum(x %in% c("-", ".", " "))
    nb <- sum(cnt); maj <- if (nb) bases[which.max(cnt)] else "-"
    set <- bases[nb > 0 & cnt / max(nb, 1) >= minfreq]
    cons <- if (gap >= n / 2 || !nb) "-" else if (!length(set)) maj else unname(iu[paste(set, collapse = "")])
    list(maj = maj, cons = cons, cons_freq = max(cnt) / n)
  })
  maj <- vapply(st, `[[`, "", "maj"); cons <- vapply(st, `[[`, "", "cons"); cf <- vapply(st, `[[`, 0, "cons_freq")
  ns <- min(n, maxrows); top <- 1.9 + ns - 1
  d <- data.frame(x = rep(cols, each = ns), i = rep(seq_len(ns), times = nc), stringsAsFactors = FALSE)
  d$b <- M[cbind(d$i, rep(seq_len(nc), each = ns))]; d$y <- 1.9 + (ns - d$i)
  d$m <- maj[rep(seq_len(nc), each = ns)]
  d$cls <- ifelse(d$b %in% c("-", ".", " "), "gap", ifelse(d$b == d$m, paste0("same", d$b), paste0("diff", ifelse(d$b %in% bases, d$b, "N"))))
  fills <- c(sameA = "#d9efd4", sameC = "#d3e3f6", sameG = "#fbe6cb", sameT = "#f7d6d6", diffA = "#2f9e57", diffC = "#2f6fc4", diffG = "#e58a2c", diffT = "#d64545", diffN = "#8e8e8e", gap = "#ececec")
  d$tc <- ifelse(startsWith(d$cls, "diff"), "white", "grey25")
  cdf <- data.frame(x = cols, cons = cons, cf = cf, stringsAsFactors = FALSE)
  sz <- if (nc <= 60) 3.6 else if (nc <= 120) 2.9 else if (nc <= 200) 2.2 else 0
  g <- ggplot2::ggplot() +
    ggplot2::geom_tile(data = d, ggplot2::aes(x = x, y = y, fill = cls), width = 0.96, height = 0.9) +
    ggplot2::geom_tile(data = cdf, ggplot2::aes(x = x, y = 0.5), width = 0.96, height = 0.9, fill = "#fff6c8") +
    ggplot2::geom_rect(data = cdf, ggplot2::aes(xmin = x - 0.48, xmax = x + 0.48, ymin = -1.3, ymax = -1.3 + cf), fill = "#e07b73", colour = NA)
  if (sz > 0) g <- g +
    ggplot2::geom_text(data = d, ggplot2::aes(x = x, y = y, label = ifelse(cls == "gap", "-", b), colour = tc), size = sz, show.legend = FALSE) +
    ggplot2::geom_text(data = cdf, ggplot2::aes(x = x, y = 0.5, label = cons), size = sz, fontface = "bold", colour = "grey15")
  g <- g + ggplot2::scale_colour_identity()
  track <- list(list(F, "Forward primer", "#2c7be5", 0), list(R, "Reverse primer", "#e8503a", 0), list(P, "Probe", "#22b07d", 1))
  yv <- top + 1.7
  for (t in track) {
    r <- t[[1]]; if (is.null(r) || anyNA(r)) next
    a <- max(r[1], from); b <- min(r[2], to); if (b < a) next
    yy <- yv + t[[4]] * 0.95
    g <- g +
      ggplot2::annotate("rect", xmin = a - 0.5, xmax = b + 0.5, ymin = yy - 0.38, ymax = yy + 0.38, fill = t[[3]]) +
      ggplot2::annotate("rect", xmin = a - 0.5, xmax = b + 0.5, ymin = -0.5, ymax = top + 0.5, fill = t[[3]], alpha = 0.10) +
      ggplot2::annotate("text", x = (a + b) / 2, y = yy, label = if (t[[2]] == "Reverse primer") "<- Reverse primer" else if (t[[2]] == "Forward primer") "Forward primer ->" else "Probe", colour = "white", size = 3.2, fontface = "bold")
  }
  ymax <- if (!is.null(P)) yv + 1.5 else yv + 0.7
  lab <- substr(nm[seq_len(ns)], 1L, 24L)
  g + ggplot2::scale_fill_manual(values = fills, guide = "none") +
    ggplot2::scale_x_continuous(breaks = seq(ceiling(from / 10) * 10, to, by = 10), expand = ggplot2::expansion(add = 0.6)) +
    ggplot2::scale_y_continuous(breaks = c(-0.8, 0.5, 1.9 + (ns - seq_len(ns))), labels = c("Conservation", "Consensus", lab), limits = c(-1.5, ymax), expand = c(0, 0)) +
    ggplot2::labs(x = "alignment column", y = NULL,
                  subtitle = sprintf("Columns %d-%d of %d. %s Solid colour = differs from the most frequent base; consensus letters are IUPAC codes (bases at >= %g%% of the sequences); the bar is the share of sequences with the most frequent base.",
                                     from, to, L, if (n > ns) sprintf("First %d of %d sequences shown (consensus and conservation use all).", ns, n) else sprintf("%d sequences.", n), 100 * minfreq)) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(), axis.text.y = ggplot2::element_text(size = 8, hjust = 1), plot.subtitle = ggplot2::element_text(size = 8.5, colour = "grey30"))
}
