# SangeRly - population genetics, core functions (no Shiny)
#   * clean alignment (complete deletion of sites with gaps / N / IUPAC, like the DnaSP default)
#   * haplotypes, diversity (S, eta, Nh, Hd, pi, k, thetaW), Tajima's D, Fu's Fs (log scale, no overflow)
#   * coalescent simulation for the p-value of Fu's Fs
#   * mismatch distribution + raggedness
#   * AMOVA (2 or 3 levels) with permutation tests, pairwise PhiST
#   * export text for DnaSP (FASTA / NEXUS / PHYLIP / MEGA), Arlequin (.arp) and PopART (NEXUS + TRAITS)
#
# Packages: ape (distances, DNAbin).  The haplotype networks use pegas (see popgen_module.R).
# The formulas are the published ones; compare a first dataset with DnaSP / Arlequin before publishing.

.pg_bases <- c("A", "C", "G", "T")

# ---------------------------------------------------------------------------------------------------
# data
# ---------------------------------------------------------------------------------------------------

# aln: named character vector (aligned).
#   gap = "excluded" (default, complete deletion): every column with a gap, N or IUPAC code in any sequence is dropped.
#   gap = "fifth": columns with N / IUPAC codes are dropped, gaps are kept and counted as a fifth state.
#   gap = "level": a column is kept when at most the fraction `miss` of its sequences has a gap / N / IUPAC code (Arlequin's AllowedLevelOfMissingData);
#     the distance functions must then skip the missing entries inside each pairwise comparison.
#   gap = "pairwise": every column is kept; sites that are not A/C/G/T are skipped inside each pairwise comparison.
# Returns list(M = matrix n x L, total = columns before, kept = L, cols = original column numbers of the kept columns).
pg_prepare <- function(aln, gap = "excluded", miss = 0) {
  x <- toupper(gsub("U", "T", unname(aln)))
  M <- do.call(rbind, strsplit(x, ""))
  rownames(M) <- names(aln)
  acgt <- M == "A" | M == "C" | M == "G" | M == "T"
  ok <- switch(gap,
    pairwise = rep(TRUE, ncol(M)),
    level = colMeans(!acgt) <= miss + 1e-12,
    fifth = colSums(!(acgt | M == "-")) == 0,
    colSums(!acgt) == 0)
  list(M = M[, ok, drop = FALSE], total = ncol(M), kept = sum(ok), cols = which(ok))
}

.pg_dnabin <- function(M) {
  x <- M; x[] <- tolower(M)
  ape::as.DNAbin(x)
}

# Population labels from the sequence names.
#   prefix: the first `nfield` fields separated by `sep`;  chars: the first `nchars` characters;
#   regex: the first () group of the expression;  all: one population.
pg_pop_from_names <- function(nm, method = "prefix", sep = "_", nfield = 1L, nchars = 3L, regex = "") {
  if (!nzchar(sep)) sep <- "_"
  out <- switch(method,
    prefix = vapply(strsplit(nm, sep, fixed = TRUE),
                    function(p) paste(p[seq_len(max(1L, min(nfield, length(p))))], collapse = sep), ""),
    chars = substr(nm, 1, nchars),
    regex = {
      r <- tryCatch(regmatches(nm, regexec(regex, nm)), error = function(e) NULL)
      if (is.null(r)) rep(NA_character_, length(nm))
      else vapply(r, function(m) if (length(m) >= 2) m[2] else NA_character_, "")
    },
    rep("all", length(nm)))
  out[is.na(out) | !nzchar(out)] <- "unassigned"
  out
}

# Haplotypes of the rows of M: H1, H2... ordered by frequency (ties: order of first appearance).
pg_haplotypes <- function(M) {
  seqs <- apply(M, 1, paste, collapse = "")
  fr <- table(seqs)
  ord <- order(-as.integer(fr), match(names(fr), seqs))
  uniq <- names(fr)[ord]
  id <- paste0("Hap_", seq_along(uniq))
  list(id = id, seq = uniq, n = as.integer(fr)[ord], hap = factor(id[match(seqs, uniq)], levels = id))
}

# ---------------------------------------------------------------------------------------------------
# neutrality tests
# ---------------------------------------------------------------------------------------------------

# Tajima (1989). The p-value follows the beta distribution of the paper (as DnaSP); p_normal is the normal approximation.
pg_tajima <- function(n, S, k) {
  na <- c(D = NA_real_, p_beta = NA_real_, p_normal = NA_real_)
  if (n < 4 || S < 1 || is.na(k)) return(na)
  i <- seq_len(n - 1)
  a1 <- sum(1 / i); a2 <- sum(1 / i^2)
  b1 <- (n + 1) / (3 * (n - 1)); b2 <- 2 * (n^2 + n + 3) / (9 * n * (n - 1))
  c1 <- b1 - 1 / a1; c2 <- b2 - (n + 2) / (a1 * n) + a2 / a1^2
  e1 <- c1 / a1; e2 <- c2 / (a1^2 + a2)
  D <- (k - S / a1) / sqrt(e1 * S + e2 * S * (S - 1))
  Dmin <- (2 / n - 1 / a1) / sqrt(e2); Dmax <- ((n + 1) / (2 * n) - 1 / a1) / sqrt(e2)
  tmp1 <- 1 + Dmin * Dmax; tmp2 <- Dmax - Dmin
  a <- -tmp1 * Dmax / tmp2; b <- tmp1 * Dmin / tmp2
  x <- (D - Dmin) / tmp2
  pb <- if (is.finite(x) && x >= 0 && x <= 1) { q <- stats::pbeta(x, b, a); 2 * min(q, 1 - q) } else NA_real_
  c(D = D, p_beta = pb, p_normal = 2 * stats::pnorm(-abs(D)))
}

# log of the unsigned Stirling numbers of the first kind |s(n, j)|, j = 1..n (recursion in log space)
pg_fs_prep <- function(n) {
  row <- 0                                       # n = 1: |s(1, 1)| = 1
  if (n >= 2) for (m in 1:(n - 1)) {
    a <- c(log(m) + row, -Inf); b <- c(-Inf, row)
    row <- pmax(a, b) + log1p(exp(-abs(a - b)))
  }
  row
}

# Fu (1997): Fs = ln(S' / (1 - S')), S' = P(K >= k | theta_pi) = sum_{j=k}^{n} |s(n,j)| theta^j / (theta (theta+1) ... (theta+n-1)).
# n = number of sequences, k = number of haplotypes, theta = mean pairwise differences (theta_pi).
pg_fs <- function(n, k, theta, logc = NULL) {
  if (is.na(theta) || theta <= 0 || k < 1 || n < 2) return(NA_real_)
  if (is.null(logc)) logc <- pg_fs_prep(n)
  j <- k:n
  lt <- logc[j] + j * log(theta)
  mx <- max(lt)
  lS <- mx + log(sum(exp(lt - mx))) - sum(log(theta + 0:(n - 1)))
  if (lS >= 0) return(Inf)
  lS - log1p(-exp(lS))
}

# One neutral coalescent sample (constant size, infinite sites, no recombination). Time is in units of 2N generations;
# theta is the per-locus value (4 N mu); mutations ~ Poisson(theta * total branch length / 2).
.pg_sim_one <- function(n, theta) {
  lineages <- as.list(seq_len(n)); birth <- rep(0, n)
  len <- numeric(0); members <- list(); tnow <- 0; m <- n
  while (m > 1) {
    tnow <- tnow + stats::rexp(1, m * (m - 1) / 2)
    pick <- sample.int(m, 2)
    for (q in pick) { len <- c(len, tnow - birth[q]); members[[length(members) + 1L]] <- lineages[[q]] }
    newl <- c(lineages[[pick[1]]], lineages[[pick[2]]])
    keep <- setdiff(seq_len(m), pick)
    lineages <- c(lineages[keep], list(newl)); birth <- c(birth[keep], tnow)
    m <- m - 1L
  }
  tot <- sum(len)
  S <- stats::rpois(1, theta * tot / 2)
  if (S == 0) return(c(S = 0, k = 1, pi = 0))
  cnt <- as.vector(stats::rmultinom(1, S, len / tot))
  H <- matrix(0L, n, S); col <- 0L
  for (b in which(cnt > 0)) {
    cc <- col + seq_len(cnt[b]); H[members[[b]], cc] <- 1L; col <- col + cnt[b]
  }
  cs <- colSums(H)
  c(S = S, k = length(unique(apply(H, 1, paste, collapse = ""))), pi = sum(2 * cs * (n - cs)) / (n * (n - 1)))
}

# p-value of Fu's Fs by coalescent simulation given theta (= theta_pi): P(Fs_sim <= Fs_obs).
pg_fs_pvalue <- function(n, theta, fs_obs, nrep = 1000, progress = NULL) {
  if (is.na(fs_obs) || is.na(theta) || theta <= 0 || n < 3) return(NA_real_)
  logc <- pg_fs_prep(n)
  sim <- rep(NA_real_, nrep)
  for (r in seq_len(nrep)) {
    s <- .pg_sim_one(n, theta)
    sim[r] <- pg_fs(n, s[["k"]], s[["pi"]], logc)
    if (!is.null(progress) && r %% 50 == 0) progress(r / nrep)
  }
  ok <- !is.na(sim)
  if (!any(ok)) return(NA_real_)
  mean(sim[ok] <= fs_obs + 1e-9)
}

# ---------------------------------------------------------------------------------------------------
# diversity statistics of the rows of M
# ---------------------------------------------------------------------------------------------------
# M may contain gaps / N / IUPAC codes (then k is computed with pairwise deletion: proportion of differences among the comparable sites of each pair x L).
# Lden = number of sites used to turn k, S into per-site values (default: ncol(M)).
pg_stats <- function(M, Lden = NULL) {
  n <- nrow(M); L <- ncol(M); if (is.null(Lden)) Lden <- L
  cnt <- rbind(colSums(M == "A"), colSums(M == "C"), colSums(M == "G"), colSums(M == "T"))
  nal <- colSums(cnt > 0)
  S <- sum(nal > 1); eta <- sum(nal - 1)
  k <- if (n > 1) sum(n^2 - colSums(cnt^2)) / (n * (n - 1)) else NA_real_
  if (n > 1 && L > 0 && any(!(M %in% c("A", "C", "G", "T")))) {
    dd <- suppressWarnings(ape::dist.dna(.pg_dnabin(M), model = "raw", pairwise.deletion = TRUE, as.matrix = FALSE))
    k <- mean(dd, na.rm = TRUE) * L; if (!is.finite(k)) k <- 0
  }
  hp <- pg_haplotypes(M); p <- hp$n / n; nh <- length(hp$n)
  p2 <- sum(p^2); p3 <- sum(p^3)
  Hd <- if (n > 1) n / (n - 1) * (1 - p2) else NA_real_
  vHd <- if (n > 2) 2 / (n * (n - 1)) * (2 * (n - 2) * (p3 - p2^2) + p2 - p2^2) else NA_real_
  a1 <- if (n > 1) sum(1 / seq_len(n - 1)) else NA_real_
  taj <- pg_tajima(n, S, k)
  list(n = n, L = L, S = S, eta = eta, Nh = nh, Hd = Hd, SD_Hd = sqrt(max(vHd, 0)),
       k = k, pi = if (Lden > 0) k / Lden else NA_real_, thetaW = if (Lden > 0 && !is.na(a1)) S / a1 / Lden else NA_real_,
       D = taj[["D"]], p_D = taj[["p_beta"]], Fs = if (n >= 2) pg_fs(n, nh, k) else NA_real_)
}

# ---------------------------------------------------------------------------------------------------
# mismatch distribution
# ---------------------------------------------------------------------------------------------------
pg_mismatch <- function(M) {
  n <- nrow(M)
  if (n < 3) return(NULL)
  d <- as.vector(stats::as.dist(ape::dist.dna(.pg_dnabin(M), model = "N", as.matrix = TRUE)))
  tb <- table(factor(d, levels = 0:max(d)))
  x <- as.numeric(tb) / sum(tb)
  list(diff = as.integer(names(tb)), freq = x, count = as.integer(tb), mean = mean(d),
       raggedness = sum(diff(c(x, 0))^2))                   # Harpending (1994)
}

# ---------------------------------------------------------------------------------------------------
# AMOVA. D = matrix of squared distances (number of differences, or a distance used as squared distance as Arlequin does)
# ---------------------------------------------------------------------------------------------------
.pg_within <- function(D, lab) {                                # sum over groups of (1 / 2 n_g) * sum of D inside the group
  idx <- split(seq_len(nrow(D)), lab)
  sum(vapply(idx, function(i) sum(D[i, i, drop = FALSE]) / (2 * length(i)), 0))
}

# two levels: individuals within populations
.pg_amova2 <- function(D, pop) {
  N <- nrow(D); pop <- factor(pop); P <- nlevels(pop); np <- as.vector(table(pop))
  ssd_t <- sum(D) / (2 * N); ssd_w <- .pg_within(D, pop); ssd_a <- ssd_t - ssd_w
  df_a <- P - 1; df_w <- N - P
  msd_a <- ssd_a / df_a; msd_w <- ssd_w / df_w
  n0 <- (N - sum(np^2) / N) / (P - 1)
  s_a <- (msd_a - msd_w) / n0; s_w <- msd_w
  list(ssd = c(a = ssd_a, w = ssd_w, t = ssd_t), df = c(a = df_a, w = df_w, t = N - 1),
       sigma = c(a = s_a, w = s_w), phiST = if ((s_a + s_w) > 0) s_a / (s_a + s_w) else NA_real_)
}

# permutation p-value: (number of permuted values >= observed + 1) / (permutations + 1); NA when the statistic is undefined
.pg_p <- function(perm, obs) {
  if (is.na(obs)) return(NA_real_)
  (sum(perm >= obs - 1e-12, na.rm = TRUE) + 1) / (sum(!is.na(perm)) + 1)
}

# three levels: individuals within populations within groups
.pg_amova3 <- function(D, pop, grp) {
  N <- nrow(D); pop <- factor(pop); grp <- factor(grp)
  P <- nlevels(pop); G <- nlevels(grp)
  np <- table(pop)
  gof <- tapply(as.character(grp), pop, function(x) x[1])           # group of each population
  Ng <- tapply(as.vector(np), gof[names(np)], sum)                  # sequences per group
  ssd_t <- sum(D) / (2 * N)
  ssd_wp <- .pg_within(D, pop); ssd_wg <- .pg_within(D, grp)
  ssd_ag <- ssd_t - ssd_wg; ssd_ap <- ssd_wg - ssd_wp
  df_ag <- G - 1; df_ap <- P - G; df_wp <- N - P
  msd_ag <- ssd_ag / df_ag; msd_ap <- ssd_ap / df_ap; msd_wp <- ssd_wp / df_wp
  sq_by_g <- tapply(as.vector(np)^2, gof[names(np)], sum)           # sum of n_p^2 inside each group
  n1 <- (N - sum(sq_by_g / Ng[names(sq_by_g)])) / (P - G)
  n2 <- (sum(sq_by_g / Ng[names(sq_by_g)]) - sum(as.vector(np)^2) / N) / (G - 1)
  n3 <- (N - sum(Ng^2) / N) / (G - 1)
  s_c <- msd_wp
  s_b <- (msd_ap - msd_wp) / n1
  s_a <- (msd_ag - msd_wp - n2 * s_b) / n3
  tot <- s_a + s_b + s_c
  list(ssd = c(ag = ssd_ag, ap = ssd_ap, wp = ssd_wp, t = ssd_t), df = c(ag = df_ag, ap = df_ap, wp = df_wp, t = N - 1),
       sigma = c(a = s_a, b = s_b, c = s_c), n = c(n1 = n1, n2 = n2, n3 = n3),
       phiST = (s_a + s_b) / tot, phiSC = s_b / (s_b + s_c), phiCT = s_a / tot)
}

# Full analysis with permutation tests. grp = NULL (or a single group) gives the two-level AMOVA.
pg_amova <- function(D, pop, grp = NULL, nperm = 1000, progress = NULL) {
  D <- as.matrix(D); pop <- factor(pop); N <- nrow(D)
  three <- !is.null(grp) && nlevels(factor(grp)) > 1 && nlevels(factor(grp)) < nlevels(pop)
  if (nlevels(pop) < 2) stop("At least 2 populations are needed for an AMOVA.")
  if (N <= nlevels(pop)) stop("The AMOVA needs more sequences than populations.")
  tick <- function(i) if (!is.null(progress) && i %% 25 == 0) progress(i / nperm)
  if (!three) {
    obs <- .pg_amova2(D, pop)
    perm <- vapply(seq_len(nperm), function(i) { tick(i); .pg_amova2(D, sample(pop))$phiST }, 0)
    return(list(levels = 2, obs = obs, nperm = nperm, p = c(phiST = .pg_p(perm, obs$phiST))))
  }
  grp <- factor(grp)
  gof <- tapply(as.character(grp), pop, function(x) unique(x))
  if (any(vapply(gof, length, 1L) > 1)) stop("A population is in more than one group. Give every population one group.")
  obs <- .pg_amova3(D, pop, grp)
  popv <- as.character(pop); grpv <- as.character(grp); pops <- levels(pop)
  gpop <- vapply(pops, function(p) as.character(gof[[p]]), "")       # group of each population
  pm <- matrix(NA_real_, nperm, 3, dimnames = list(NULL, c("phiST", "phiSC", "phiCT")))
  for (i in seq_len(nperm)) {
    tick(i)
    # PhiST: sequences permuted among all populations (the group follows the new population)
    p1 <- sample(popv); a <- .pg_amova3(D, p1, gpop[p1])
    # PhiSC: sequences permuted among populations inside each group
    p2 <- popv; for (g in unique(grpv)) { ii <- which(grpv == g); p2[ii] <- sample(popv[ii]) }
    b <- .pg_amova3(D, p2, grpv)
    # PhiCT: whole populations permuted among groups
    gnew <- sample(gpop); names(gnew) <- pops
    cc <- .pg_amova3(D, popv, gnew[popv])
    pm[i, ] <- c(a$phiST, b$phiSC, cc$phiCT)
  }
  p <- c(phiST = .pg_p(pm[, "phiST"], obs$phiST), phiSC = .pg_p(pm[, "phiSC"], obs$phiSC),
         phiCT = .pg_p(pm[, "phiCT"], obs$phiCT))
  list(levels = 3, obs = obs, nperm = nperm, p = p)
}

# Pairwise PhiST (two-level AMOVA of every pair of populations) and its p-values.
pg_pairwise_phist <- function(D, pop, nperm = 200, progress = NULL) {
  D <- as.matrix(D); pop <- factor(pop); pops <- levels(pop); P <- length(pops)
  phi <- pv <- matrix(NA_real_, P, P, dimnames = list(pops, pops))
  tot <- P * (P - 1) / 2; done <- 0
  for (i in seq_len(P - 1)) for (j in (i + 1):P) {
    ii <- which(pop %in% pops[c(i, j)])
    sub <- D[ii, ii, drop = FALSE]; lab <- droplevels(pop[ii])
    if (length(ii) > 2 && all(table(lab) >= 1)) {
      o <- .pg_amova2(sub, lab)$phiST
      pr <- vapply(seq_len(nperm), function(r) .pg_amova2(sub, sample(lab))$phiST, 0)
      phi[i, j] <- phi[j, i] <- o
      pv[i, j] <- pv[j, i] <- .pg_p(pr, o)
    }
    done <- done + 1; if (!is.null(progress)) progress(done / tot)
  }
  list(phi = phi, p = pv)
}

# ---------------------------------------------------------------------------------------------------
# export text
# ---------------------------------------------------------------------------------------------------
# Arlequin "Compute pairwise differences (pi)": mean distance between populations (PiXY), within (PiX) and the corrected value PiXY - (PiX + PiY) / 2
pg_pop_pi <- function(D, pop) {
  D <- as.matrix(D); pop <- factor(pop); pops <- levels(pop); P <- length(pops)
  xy <- matrix(NA_real_, P, P, dimnames = list(pops, pops)); w <- stats::setNames(rep(NA_real_, P), pops)
  for (i in seq_len(P)) {
    ii <- which(pop == pops[i])
    if (length(ii) > 1) w[i] <- sum(D[ii, ii][upper.tri(D[ii, ii])]) / choose(length(ii), 2)
    for (j in seq_len(P)) { jj <- which(pop == pops[j]); if (i != j) xy[i, j] <- mean(D[ii, jj, drop = FALSE]) }
  }
  net <- xy - (outer(w, rep(1, P)) + outer(rep(1, P), w)) / 2; diag(net) <- 0
  xy[row(xy) == col(xy)] <- w
  list(xy = xy, within = w, net = net)
}
pg_reynolds <- function(phi) { r <- -log(1 - phi); r[!is.finite(r)] <- NA_real_; r }

pg_safe_names <- function(nm) make.unique(gsub("[^A-Za-z0-9_]", "_", nm), sep = "_")

pg_fasta <- function(seqs, nms, width = 60) {
  out <- character()
  for (i in seq_along(seqs)) {
    s <- seqs[i]
    out <- c(out, paste0(">", nms[i]),
             if (nchar(s)) substring(s, seq(1, nchar(s), by = width), pmin(seq(width, nchar(s) + width - 1, by = width), nchar(s))) else "")
  }
  paste(out, collapse = "\n")
}

# NEXUS for DnaSP, MEGA, PAUP...; with a SETS block (one TAXSET per population) when pops is given.
pg_nexus <- function(seqs, nms, pops = NULL) {
  w <- max(nchar(nms))
  out <- c("#NEXUS", "", "BEGIN DATA;", sprintf("  DIMENSIONS NTAX=%d NCHAR=%d;", length(seqs), nchar(seqs[1])),
           "  FORMAT DATATYPE=DNA MISSING=? GAP=-;", "  MATRIX",
           sprintf("  %s %s", formatC(nms, width = -w), seqs), "  ;", "END;")
  if (!is.null(pops)) {
    out <- c(out, "", "BEGIN SETS;")
    for (p in unique(pops)) out <- c(out, sprintf("  TAXSET %s = %s;", pg_safe_names(p)[1], paste(nms[pops == p], collapse = " ")))
    out <- c(out, "END;")
  }
  paste(out, collapse = "\n")
}

# PopART: haplotypes + TRAITS block (counts of each haplotype in each population).
pg_nexus_popart <- function(hap_seq, hap_id, counts) {
  w <- max(nchar(hap_id)); pops <- pg_safe_names(colnames(counts))
  out <- c("#NEXUS", "", "BEGIN TAXA;", sprintf("  DIMENSIONS NTAX=%d;", length(hap_id)),
           paste0("  TAXLABELS ", paste(hap_id, collapse = " "), ";"), "END;", "",
           "BEGIN CHARACTERS;", sprintf("  DIMENSIONS NCHAR=%d;", nchar(hap_seq[1])),
           "  FORMAT DATATYPE=DNA MISSING=? GAP=-;", "  MATRIX",
           sprintf("  %s %s", formatC(hap_id, width = -w), hap_seq), "  ;", "END;", "",
           "BEGIN TRAITS;", sprintf("  DIMENSIONS NTRAITS=%d;", ncol(counts)),
           "  FORMAT LABELS=YES MISSING=? SEPARATOR=COMMA;",
           paste0("  TRAITLABELS ", paste(pops, collapse = " "), ";"), "  MATRIX",
           sprintf("  %s %s", formatC(hap_id, width = -w), apply(counts, 1, paste, collapse = ",")), "  ;", "END;")
  paste(out, collapse = "\n")
}

# PHYLIP (names of at most 10 characters, as the strict format wants). Returns list(text, map).
pg_phylip <- function(seqs, nms) {
  short <- substr(gsub("[^A-Za-z0-9_]", "_", nms), 1, 10)
  short <- make.unique(short, sep = "")
  short <- substr(short, 1, 10)
  if (anyDuplicated(short)) short <- sprintf("S%09d", seq_along(nms))
  text <- c(sprintf(" %d %d", length(seqs), nchar(seqs[1])), paste0(formatC(short, width = -10), seqs))
  list(text = paste(text, collapse = "\n"), map = data.frame(Phylip_name = short, Sequence = nms, stringsAsFactors = FALSE))
}

pg_mega <- function(seqs, nms, title = "SangeRly alignment") {
  out <- c("#mega", paste0("!Title ", title, ";"), "!Format DataType=DNA indel=-;", "")
  for (i in seq_along(seqs)) out <- c(out, paste0("#", nms[i]), seqs[i])
  paste(out, collapse = "\n")
}

# Arlequin project file in the haplotype-frequency format that Arlequin AND DnaSP read:
# one sample per population, each line of SampleData = haplotype id, frequency, sequence.
# groups = NULL or a named vector (names = populations, values = group) for the structure (AMOVA).
pg_arlequin <- function(hap_seq, hap_id, counts, groups = NULL, title = "SangeRly") {
  pops <- colnames(counts)
  out <- c("[Profile]", sprintf("\tTitle=\"%s\"", title), sprintf("\tNbSamples=%d", length(pops)),
           "\tGenotypicData=0", "\tGameticPhase=0", "\tDataType=DNA", "\tLocusSeparator=NONE", "\tMissingData=\"?\"", "",
           "[Data]", "", "[[Samples]]")
  for (p in pops) {
    ii <- which(counts[, p] > 0)
    out <- c(out, "", sprintf("\tSampleName=\"%s\"", p), sprintf("\tSampleSize=%d", sum(counts[, p])), "\tSampleData={",
             sprintf("\t\t%s\t%d\t%s", hap_id[ii], counts[ii, p], hap_seq[ii]), "\t}")
  }
  if (!is.null(groups) && length(unique(groups)) >= 1) {
    ug <- unique(groups)
    out <- c(out, "", "[[Structure]]", "", "\tStructureName=\"SangeRly groups\"", sprintf("\tNbGroups=%d", length(ug)))
    for (g in ug) out <- c(out, sprintf("\t#%s", g), "\tGroup={", sprintf("\t\t\"%s\"", names(groups)[groups == g]), "\t}")
  }
  paste(out, collapse = "\n")
}

# ---------------------------------------------------------------------------------------------------
# matching the sequence names with the keys of a CSV (accession numbers, with or without the version)
# ---------------------------------------------------------------------------------------------------
.pg_acc <- function(x) {
  x <- toupper(x)
  lapply(regmatches(x, gregexpr("[A-Z]{1,2}_?[0-9]{5,9}(\\.[0-9]+)?", x)), function(v) unique(sub("\\.[0-9]+$", "", v)))
}

# Returns list(index = position of the CSV key for every sequence name (NA = no match), how = "exact" / "prefix" / "accession").
#   exact: same text;  prefix: the name starts with the key ("MT298407.1" in "MT298407.1_Orthetrum_sabina_Italy");
#   accession: the accession number of the key is found anywhere in the name (version ignored).
pg_match_keys <- function(seq_names, keys) {
  norm <- function(x) gsub("[^A-Z0-9]+", "_", toupper(trimws(as.character(x))))
  sn <- norm(seq_names); kn <- norm(keys)
  idx <- rep(NA_integer_, length(sn)); how <- rep("", length(sn))
  m <- match(sn, kn); ok <- !is.na(m); idx[ok] <- m[ok]; how[ok] <- "exact"
  for (i in which(is.na(idx))) {
    hit <- which(vapply(kn, function(k) nzchar(k) && startsWith(sn[i], paste0(k, "_")), TRUE))
    if (length(hit)) { idx[i] <- hit[1]; how[i] <- "prefix" }
  }
  acc_n <- .pg_acc(seq_names)
  acc_k <- vapply(.pg_acc(keys), function(v) if (length(v)) v[1] else NA_character_, "")
  for (i in which(is.na(idx))) {
    hit <- which(acc_k %in% acc_n[[i]])
    if (length(hit)) { idx[i] <- hit[1]; how[i] <- "accession" }
  }
  list(index = idx, how = how)
}



# Example data: ITS alignment of Durio zibethinus (60 sequences, 654 sites) and the table of locations (columns ID, Species, Location, as in an Excel sheet).
.pg_durio_hap <- c(
  "CGCGGACGCGTTAACGAAACACCGAGGGGGCGGGGGTGCCCCCGCGCCCCGAGCCCCCTCCGGCCCAGGGAGGGCGCGACGGTCGCTTCGTCCCCCCTACGGGGCGGGAGGAACGATCGGCGTCGTCGCCCTCCGCCGAACAACGAACCCCCGGCGCGAATCGCGCCAAGGAACCGAAATGAAAACGGAGCGCGCGCGCCGTCGCCATCCCGTCCGCGGTGCGTGCGCGGCGGCGCCGCGACTCCTTTTGTCGCGAAAAAACGAAACGACTCTCGGCAACGGATATCTTGGCTCTCGCATCGATGAAGAACGTAGCGAAATGCGATACTTGGTGTGAATTGCAGAATCCCGTGAACCATCGAGTCTTTGAACGCAAGTTGCGCCCCAAGCCATTAGGCCAGGGCACGCCTGCCTGGGTGTCATGCACCGTCGCCCCCTCCGACCCCCCGCCCTCCGGGCTCCGGGGGAAGGGGGCGGAAGATGGCCTCCCGCGATCCCGTGGGTCGCGGTTGGCCCAAAAGCGAGTCCCCGGCAACGACAGCGCCGTGACCGTCGGTGGTGGAGCTCCCGGCGTGCCTCGTCCCCGGTGACGCGCGCTCGTTGTCCATCCGGACCCCTTCGGACCCCCCGCATCGCTCGATCGATGCTCGCGTC",
  "CGCGGACGCGTTAACGAAACACCGAGGGGGCGAGGGTGCCCCCGCGCCCCGAGCCCCCTCCGGCCCCGGGAGGGCGCGACGGTCGCTTCGTCCCCCCTACGGGGCGGGCGAAACGATCGGCGTCGTCGCCCTCCGCCGAACAACGAACCCCCGGCGCGAATCGCGCCAAGGAACCGAAACGAAAACGGAGCGCGCGCGCCGTCGCCATCCCGTCCGCGGTGCGTGCGCGGCGGCGCCGCGACTCCTTTTGTCGCGAAAAAACGAAACGACTCTCGGCAACGGATATCTCGGCTCTCGCATCGATGAAGAAAGTAGCGGAATGCGAGACTTGGTGTGAATTGCAGAATCCCGTGAACCATCAAGTCTTTGAACGCAAGTTGCACCCCCAGCCATTCGGCCAGGGCACGCCTGCCTGGGTGTCATGCACCGTCGCCCCCTCCGACCCCCCGCCCTCCGGGCTCCGGGGGAAGGGGGCGGAAGATGGCCTCCCGCGATCCCGTGGGTCGCGGTTGGCCCAAAAGCGAGTCCCCGGCAACGACAGCGCCGTGACCGTCGGTGGTGGAGCTCCCGACGTGCCTCGTCCCCGGTGACGCGCGCTCGTTGTCCATCCGGACCCCTTCGGACCCCCCGCATCGCTCGATCGATGCTCGCGTC",
  "CGCGGACGCGTTAACGAAACACCGAGGGGGCGAGGGTGCCCCCGCGCCCCGAGCCCCCTCCGGCCCCGGGAGGGCGCGACGGTCGCTTCGTCCCCCCTACGGGGAGGGCGGAACGATCGGCGTCGTCGCCCTCCGCCGAACAACGAACCCCCGGCGCGAATCGCGCCAAGGAACCGAAACGAAAACGGAGCGCGCGCGCCGTCGCCATCCCGTCCGCGGTGCGTGCGCGGCGGCGCCGTGACTCCTTTTGTCGCGAAAAAACGAAACGACTCTCGGCAACGGATATCTTGGCTCTCGCATCGATGAAGAACGTAGCGAAATGCGATACTTGGTGTAAATTGCAGAATCCCATGAACCGTCGAGTTTTTGAATTCAAGTTACGCCCCAAGCCGTTAGGTTAGGGAACTCTTGCCTGGGTGTCATGCACCGTCGCCCCTTCCGACCCCCCGCCCTTCGGGCTCCGGGGGAAGGGGGCGGAAGATGGCCTCCCGCGATCCCGTGGGTCACGGTTGGCCCAAAAGTGAGTCCTTGGCAATGACAGCGCCGTGACCGTCGGTGGTAGAGCTCCCGACGTGCCTCGTCCTCGGTGACGCGCGCTCGTTGTCCATCCGGACCCCTTCGGACCCCCCGCATCGCTCGATCGATGCTCGCGTC",
  "CGAGGGTGTGTTAAAAACACACCGAGGGGGCGGGGGTGGCCCCGCGCCCCGAGCCCCCTCCGGCCCCGGGAGGGCGCGACGGTCGCTTTGTCCCCCCTACGGGGCGGGCGGAACGATTGGGGTTGTCGCCCTCCGCCGAACAACGAACCCCCGGCGCGAATCGCGCCAAGGAACCGAAACGAAAACGGAGCGCGCGCGCCGTTGCCATCCCGTTCGCGGTGTGTGCGCGGCGGCGCCGCGACTCCTTTTGTCGCGAAAAAACGAAACGACTCTCGGCAACGGATATCTCGGCTCTCGCATCGATGAAGAACGTAGCGAAATGCGATACTTGGTGTGAATTGCAGAATCCCGTGAACCATCGAGTCTTTGAACGCAAGTTGCGCCCCAAGCCATTAGGCCAGGGCACGCCTGCCTGGGTGTCATGCACCGTCGCCCCCTCCGACCCCCCGCCCTCCGGGCTCCGGGGGAAGGGGGCGGAAGATGGCCTCCCGCGATCCCGTGGGTCGCGGTTGGCCCAAAAGCGAGTCCCCGGCAACGACAGCGCCGTGACCGTCGGTGGTGGAGCTCCCGGCGTGCCTCGTCCCCGGTGACGCGCGCTCGTTGTCCATCCGGACCCCTTCGGACCCCCCGCATCGCTCGATCGATGCTCGCGTC",
  "CGCGGACGCGTTAAAGAAACACCGAGGGGGCGGGGGTGCCCCCGCGCCCCGAGCCCCCTCCGGCCCAGGGAGGGCGCGACGGTCGCTTCGTCCCCCCTACGGGGCGGGCGGAACGATCGGCGTCGTCGCCCTCCGCCGAACAACGAACCCCCGGCGCGAATCGCGCCAAGGAACCGAAACGAAAACGGAGCGCGCGCGCCGTCGCCATCCCGTCCGCGGTGCGTGCGCGGCGGCGCCGCGACTCCTTTTGTCGCGAAAAAACGAAACGACTCTCGGCAACGGATATCTCGGCTCTCGCATCGATGAAGAACGTAGCGAAATGCGATACTTGGTGTGAATTGCAGAATCCCGTGAACCATCGAGTCTTTGAACGCAAGTTGCGCCCCAAGCCATTAGGCCAGGGCACGCCTGCCTGGGTGTCATGCACCGTCGCCCCCTCCGACCCCCCGCCCTCCGGGCTCCGGGGGAAGGGGGCGGAAGATGGCCCCCCGCGATCCCGGGGGTCGCGGTGGGCCCAAAAGCGAGTCCCCGGCAACGACAGCGCCGTGACCGTCGGGGGTGGAGCCCCCGACGTGCCTCGCCCCCGGTGACGCGCGCTCGTTGTCCATCCGGACCCCTTCGGACCCCCCGCATCGCTCGATCGATGCTCGCGTC",
  "CGCGGACGCGTTAACGAAACACCGAGGGGGCGGGGGTGCCCCCGCGCCCCGAGCCCCCTCCGGCCCAGGGAGGGCGCGACGGTCGCTTCGTCCCCCCTACGGGGCGGGCGGAACGATCGGCGTCTTCGCCCTCCGCCGAACAACGAACCCCCGGCGCGAATCGCGCCAAGGAACCGAAACAAAAACGGAGCGCGCGCGCCGTCGCCATCCCGTCCGCGGTGCGTGCGCGGCGGCGCCGCGACTCCTTTTGTCGCGAAAAAACGAAACGACTCTCGGCAACGGATATCTCGGCTCTCGCATCGATGAAAAACGTAGCGAAATGCGATACTTGGTGTGAATTGCAGAATCCCGTGAACCATCGAGTCTTTGAACGCAAGTTGCGCCCCAAGCCATTAGGCCAGGGCACGCCTGCCTGGGTGTCATGCACCGTCGCCCCCTCCGACCCCCCGCCCTCCGGGCTCCGGGGGAAGGGGGCGGAAGATGGCCTCCCGCGATCCCGTGGGTCGCGGTTGGCCCAAAAGCGAGTCCCCGGCAACGACAGCGCCGTGACCGTCGGTGGTGGAGCTCCCGGCGTGCCTCGTCCCCGGTGACGCGCGCTCGTTGTCCATCCGGACCCCTTCGGACCCCCCGCATCGCTCGATCGATGCTCGCGTC",
  "CGCGGACGCGTTAACGAAACACCGAGGGGGCGGGGGTGCCCCCGCGCCCCGAGCCCCTTCCGGCCCAGGGAGGGCGCGACGGTCGCTTTGTCCCCCCTACGGGGCGGGCGGAACGATCGGCGTCGTCGCCCTCCGCCGAACAACGAACCCCCGGCGCGAATCGCGCCAAGGAACCGAAACGAAAACGGAGCGCGCGCGCCGTCGCCATCCCGTCCGCGGTGCGTGGGCGGCGGCGCCGCGACTCCTTTTGTCGCGAAAAAACGAAACGACTCTCGGCAACGGATATCTCGGCTCTCGCATCGATGAAGAACGTAGCGAAATGCGATACTTGGTGTGAATTGCAGAATCCCGTGAACCATCGAGTCTTTGAACGCAAGTTGCGCCCCAAGCCATTAGGCCAGGGCACGCCTGCCTGGGTGTCATGCACCGTCGCCCCCTCCGACCCCCCGCCCTCCGGGCTCCGGGGGAAGGGGGCGGAAGATGGCCTCCCGCGATCCCGTGGGTCGCGGTTGGCCCAAAAGCGAGTCCCCGGCAACGACAGCGCCGTGACCGTCGGTGGTGGAGCTCCCGGCGTGCCTCGTCCCCGGTGACGCGCGCTCGTTGTCCATCCGGACCCCTTCGGACCCCCCGCATCGCTCGATCGATGCTCGCGTC",
  "CGCGGACGCGTTAACGAAACACCGAGGGGGCGGGGGTGCCCCCGCGCCCCGAGCCCCCTCCGGCCCAGGGAGGGCGCGACGGTCGCTTCGTCCCCCCTACGGGGCGGGCGGAACGATCGGCGTCGTCGCCCTCCGCCGAACAACGAACCCCCGGCGCGAATCGCGCCAAGGAACCGAAACGAAAACGGAGCGCGCGCGCCGTCGCCATCCCGTCCGCGGTGCGTGCGCGGCGGCGCCGCGACTCCTTTTGTCGCGAAAAAACGAAACGACTCTCGGCAACGGATATCTCGGCTCTCGCATCGATGAAGAACGTAGCGAAATGCGATACTTGGTGTGAATTGCAGAATCCCGTGAACCATCGAGTCTTTGAACGCAAGTTGCGCCCCAAGCCATTAGGCCAGGGCACGCCTGCCTGGGTGTCATGCACCGTCGCCCCCTCCGACCCCCCGCCCTCCGGGCTCCGGGGGAAGGGGGGGGAAGATGGCCTCCCGCGATCCCGTGGGTCGCGGTTGGCCCAAAAGCGAGTCCCCGGCAACGACAGCGCCGTGACCGTCGGTGGGGGAGCTCCCGACGTGCCTCGTCCCCGGTGACGCGCGCTCGTTGTCCATCCGGACCCCTTCGGACCCCCCGCATCGCTCGATCGATGCTCGCGTC",
  "CGCGGACGCGTTAACGAAAGACCGAGGGGGCGGGGGTGCCCCCGCGCCCCGAGCCCCCTCCGGCCCAGGGAGGGCGCGACGGTCGCTTCGTCCCCCCTACGGGGCGGGCGGAACGATCGGCGTCGTCGCCCTCCGCCGAACAACGAACCCCCGGCGCGAATCGCGCCAAGGAACCGAAACGAAAACGGAGCGCGCGCGCCGTCGCCATCCCGTCCGCGGTGCGTGCGCGGCGGCGCCGCGACTCCTTTTGTCGCGAAAAAACGAAACGACTCTCGGCAACGGATATCTCGGCTCTCGCATCGATGAAGAACGTAGCGAAATGCGATACTTGGTGTGAATTGCAGAATCCCGTGAACCATCGAGTCTTTGAACGCAAGTTGCGCCCCAAGCCATTAGGCCAGGGCACGCCTGCCTGGGTGTCATGCACCGTCGCCCCCTCCGACCCCCCGCCCTCCGGGCTCCGGGGGAAGGGGGCGGAAGATGGCCTCCCGCGATCCCGTGGGTCGCGGTTGGCCCAAAAGCGAGTCCCCGGCAACGACAGCGCCGTGACCGTCGGTGGTGGAGCTCCCGGCGTGCCTCGTCCCCGGTGACGCGCGCTCGTTGTCCATCCGGACCCCTTCGGACCCCCCGCATCGCTCGATCGATGCTCGCGTC",
  "CGCGGACGCGTTAACGAAACACCGAGGGGGCGGGGGTGCCCCCGCGCCCCGAGCCCCCTCCGGCCCAGGGAGGGCGCGACGGTCGCTTCGTCCCCCCTACGGGGCGGGCGGAACGATCGGCGTCGTCGCCCTCCGCCGAACAACGAACCCCCGGCGCGAATCGCGCCAAGGAACCGAAACGAAAACGGAGCGCGCGCGCCGTCGCCATCCCGTCCGCGGTGCGTGCGCGGCGGCGCCGCGACTCCTTTTGTCGCGAAAAAACGAAACGACTCTCGGCAACGGATATCTCGGCTCTCGCATCGATGAAGAACGTAGCGAAATGCGATACTTGGTGTGAATTGCAGAATCCCGTGAACCATCGAGTCTTTGAACGCAAGTTGCGCCCCAAGCCATTAGGCCAGGGCACGCCTGCCTGGGTGTCATGCACCGTCGCCCCCTCCGACCCCCCGCCCTCCGGGCTCCGGGGGAAGGGGGCGGAAGATGGCCTCCCGCGATCCCGTGGGTCGCGGTTGGCCCAAAAGCGAGTCCCCGGCAACGACAGCGCCGTGACCGTCGGTGGTGGAGCTCCCGGCGTGCCTCGTCCCCGGTGACCCGCGCTCGTTGTCCATCCGGACCCCTTCGGACCCCCCGCATCGCTCGATCGATGCTCGCGTC",
  "CGCGGACGCGTTAACGAAACACCGAGGGGGCGGGGGTGCCCCCGCGCCCCGAGCCCCCTCCGGCCCAGGGAGGGCGCGACGGTCGCTTTGTCCCCCCTACGGGGCGGGCGGAACGATCGGCGTCGTCGCCCTCCGCCGAACAACGAACCCCCGGCGCGAATCGCGCCAAGGAACCGAAACGAAAACGGAGCGCGCGCGCCGTCGCCATCCCGTCCGCGGTGCGTGCGCGGCGGCGCCGCGACTCCTTTTGTCGCGAAAAAACGAAACGACTCTCGGCAACGGATATCTCGGCTCTCGCATCGATGAAGAACGTAGCGAAATGCGATACTTGGTGTGAATTGCAGAATCCCGTGAACCATCGAGTCTTTGAACGCAAGTTGCGCCCCAAGCCATTAGGCCAGGGCACGCCTGCCTGGGTGTCATGCACCGTCGCCCCCTCCGACCCCCCGCCCTCCGGGCTCCGGGGGAAGGGGGCGGAAGATGGCCTCCCGCGATCCCGTGGGTCGCGGTTGGCCCAAAAGCGAGTCCCCGGCAACGACAGCGCCGTGACCGTCGGTGGTGGAGCTCCCGGCGTGCCTCGTCCCCGGTGACGCGCGCTCGTTGTCCATCCGGACCCCTTCGGACCCCCCGCATCGCTCGATCGATGCTCGCGTC",
  "CGCGGACGCGTTAACGAAACACCGAGGGGGCGGGGGTGCCCCCGCGCCCCGAGCCCCCTCCGGCCCAGGGAGGGCGCGACGGTCGCTTCGTCCCCCCTACGGGGCGGGCGGAACGATCGGCGTCGTCGCCCTCCGCCGAACAACGAACCCCCGGCGCGAATCGCGCCAAGGAACCGAAACGAAAACGGAGCGCGCGCGCCGTCGCCATCCCGTCCGCGGTGCGTGCGCGGCGGCGCCGCGACTCCTTTTGTCGCGAAAAAACGAAACGACTCTCGGCAACGGATATCTCGGCTCTCGCATCGATGAAGAACGTAGCGAAATGCGATACTTGGTGTGAATTGCAGAATCCCGTGAACCATCGAGTCTTTGAACGCAAGTTGCGCCCCAAGCCATTAGGCCAGGGCACGCCTGCCTGGGTGTCATGCACCGTCGCCCCCTCCGACCCCCCGCCCTCCGGGCTCCGGGGGAAGGGGGCGGAAGATGGCCTCCCGCGATCCCGTGGGTCGCGGTTGGCCCAAAAGCGAGTCCCCGGCAACGACAGCGCCGTGACCGTCGGTGGTGGAGCTCCCGGCGTGCCTCGTCCCCGGTGACGCGCGCTCGTTGTCCATCCGGACCCCTTCGGACCCCCCGCATCGCTCGATCGATGCTCGCGTC",
  "CGCGGACGCGTTAACGAAACACCGAGGGGGCGGGGGTGCCCCCGCGCCCCGAGCCCCCTCCGGCCCAGGGAGGGCGCGACGGTCGCTTCGTCCCCCCTACGGGGCGGGCGGAACGATCGGCGTCGTCGCCCTCCGCCGAACAACGAACCCCCGGCGCGAATCGCGCCAAGGAACCGAAACGAAAACGGAGCGCGCGCGCCGTCGCCATCCCGTCCGCGGTGCGTGCGCGGCGGCGCCGCGACTCCTTTTGTCGCGAAAAAACGAAACGACTTTCGGCAACGGATATCTCGGCTCTCGCATCGATGAAGAACGTAGCGAAATGCGATACTTGGTGTGAATTGCAGAATCCCGTGAACCATCGAGTCTTTGAACGCAAGTTGCGCCCCAAGCCATTAGGCCAGGGCACGCCTGCCTGGGTGTCATGCACCGTCGCCCCCTCCGACCCCCCGCCCTCCGGGCTCCGGGGGAAGGGGGCGGAAGATGGCCTCCCGCGATCCCGTGGGTCGCGGTTGGCCCAAAAGCGAGTCCCCGGCAACGACAGCGCCGTGACCGTCGGTGGTGGAGCTCCCGACGTGCCTCGTCCCCGGTGACGCGCGCTCGTTGTCCATCCGGACCCCTTCGGACCCCCCGCATCGCTCGATCGATGCTCGCGTC",
  "CGCGGACGCGTTAACGAAACACCGAGGGGGCGGGGATGCCCCCGCGCCCCGAGCCCCCTCCGGCCCAGGGAGGGCGCGACGGTCGCTTCGTCCCCCCTACGGGGCGGGCGGAACGATCGGCGTCGTCGCCCTCCGCCGAACAACGAACCCCCGGCGCGAATCGCGCCAAGGAACCGAAACGAAAACGGAGCGCGCGCGCCGTCGCCATCCCGTCCGCGGTGCGTGCGCGGCGGCGCCGCGACTCCTTTTGTCGCGAAAAAACGAAACGACTCTCGGCAACGGATATCTCGGCTCTCGCATCGATGAAGAACGTAGCGAAATGCGATACTTGGTGTGAATTGCAGAATCCCGTGAACCATCGAGTCTTTGAACGCAAGTTGCGCCCCAAGCCATTAGGCCAGGGCACGCCTGCCTGGGTGTCATGCACCGTCGCCCCCTCCGACCCCCCGCCCTCCGGGCTCCGGGGGAAGGGGGCGGAAGATGGCCTCCCGCGATCCCGTGGGTCGCGGTTGGCCCAAAAGCGAGTCCCCGGCAACGACAGCGCCGTGACCGTCGGTGGTGGAGCTCCCGACGTGCCTCGTCCCCGGTGACGCGCGCTCGTTGTCCATCCGGACCCCTTCGGACCCCCCGCATCGCTCGATCGATGCTCGCGTC",
  "CGCGGACGCGTTAACGAAACACCGAGGGGGCGGGGGTGCCCCCGCGCCCCGAGCCCCCTCCGGCCCAGGGAGGGCGCGACGGTCGCTTCGTCCCCCCTACGGGGCGGGCGGAACGATCGGCGTCGTCGCCCTCCGCCGAACAACGAACCCCCGGCGCGAATCGCGCCAAGGAACCGAAACGAAAACGGAGCGCGCGCGCCGTCGCCATCCCGTCCGCGGTGCGTGCGCGGCGGCGCCGCGACTCCTTTTGTCGCGAAAAAACGAAACGACTCTCGGCAACGGATATCTCGGCTCTCGCATCGATGAAGAACGTAGCGAAATGCGATACTTGGTGTGAATTGCAGAATCCCGTGAACCATCGAGTCTTTGAACGCAAGTTGCGCCCCAAGCCATTAGGCCAGGGCACGCCTGCCTGGGTGTCATGCACCGTCGCCCCCTCCGACCCCCCGCCCTCCGGGCTCCGGGGGAAGGGGGCGGAAGATGGCCTCCCGCGATCCCGTGGGTCGCGGTTGGCCCAAAAGCGAGTCCCCGGCAACGACAGCGCCGTGACCGTCGGTGGTGGAGCTCCCGACGTGCCTCGTCCCCGGTGACGCGCGCTCGTTGTCCATCCGGACCCCTTCGGACCCCCCGCATCGCTCGATCGATGCTCGCGTC"
)
.pg_durio_map <- c(
  "Durio_zibethinus_cultivar_Kinoso1" = 1,
  "Durio_zibethinus_cultivar_Kinoso3" = 1,
  "Durio_zibethinus_cultivar_Durian_Kunyit1" = 2,
  "Durio_zibethinus_cultivar_Durian_Kunyit7" = 2,
  "Durio_zibethinus_cultivar_Duriat1" = 3,
  "Durio_zibethinus_cultivar_Duriat10" = 3,
  "ON815633.1_Durio_zibethinus_cultivar_Duriat" = 3,
  "MF629759.1_Durio_zibethinus_isolate_J" = 4,
  "ON817266.1_Durio_zibethinus_cultivar_Durian_Kunyit" = 2,
  "MT230652.1_Durio_zibethinus_voucher_DUR17" = 5,
  "ON817274.1_Durio_zibethinus_cultivar_Kinoso" = 1,
  "MF629757.1_Durio_zibethinus_isolate_H" = 6,
  "MF629769.1_Durio_zibethinus_isolate_T" = 6,
  "MW316655.1_Durio_zibethinus_isolate_Chuong_Bo_TG" = 7,
  "MT230651.1_Durio_zibethinus_voucher_DUR12" = 8,
  "OK052902.1_Durio_zibethinus_voucher_A4716_5.8S" = 9,
  "MF629772.1_Durio_zibethinus_isolate_W" = 10,
  "MW186242.1_Durio_zibethinus_isolate_Ri-6" = 11,
  "OR775700.1_Durio_zibethinus_voucher_ITS-Ri6-BT" = 12,
  "OR775701.1_Durio_zibethinus_voucher_ITS-MT-BT" = 12,
  "MW318976.1_Durio_zibethinus_isolate_Chin_Hoa_BT" = 12,
  "MW318979.1_Durio_zibethinus_isolate_Bi_BT" = 12,
  "MW318978.1_Durio_zibethinus_isolate_Kho_Qua_Xanh_VL" = 12,
  "MW318974.1_Durio_zibethinus_isolate_Monthong_CT" = 12,
  "MW318975.1_Durio_zibethinus_isolate_Monthong_VL" = 12,
  "MW318973.1_Durio_zibethinus_isolate_Chin_Hoa_CT" = 12,
  "AF233304.1_Durio_zibethinus" = 12,
  "AF287710.1_Durio_zibethinus" = 12,
  "MT230656.1_Durio_zibethinus_voucher_DUR33" = 12,
  "MT230655.1_Durio_zibethinus_voucher_DUR24" = 12,
  "MF629753.1_Durio_zibethinus_isolate_D" = 12,
  "MF629754.1_Durio_zibethinus_isolate_E" = 12,
  "MF629755.1_Durio_zibethinus_isolate_F" = 12,
  "MF629756.1_Durio_zibethinus_isolate_G" = 12,
  "MF629761.1_Durio_zibethinus_isolate_L" = 12,
  "MF629763.1_Durio_zibethinus_isolate_N" = 12,
  "MF629764.1_Durio_zibethinus_isolate_O" = 12,
  "MF629766.1_Durio_zibethinus_isolate_Q" = 12,
  "MF629767.1_Durio_zibethinus_isolate_R" = 12,
  "MF629771.1_Durio_zibethinus_isolate_V" = 12,
  "MF629774.1_Durio_zibethinus_isolate_Y" = 12,
  "MF629775.1_Durio_zibethinus_isolate_Z" = 12,
  "MF629777.1_Durio_zibethinus_isolate_BB" = 12,
  "MF629778.1_Durio_zibethinus_isolate_CC" = 12,
  "MF629779.1_Durio_zibethinus_isolate_DD" = 12,
  "MW316477.1_Durio_zibethinus_isolate_Ri-6_VCAQMN" = 12,
  "MW316599.1_Durio_zibethinus_isolate_Monthong_TG" = 12,
  "MW316621.1_Durio_zibethinus_isolate_Ri-6_TG" = 12,
  "MW316622.1_Durio_zibethinus_isolate_Ri-6_CT" = 12,
  "MW316637.1_Durio_zibethinus_isolate_Sua_Hat_Lep_CT" = 12,
  "MW316650.1_Durio_zibethinus_isolate_Sau_Huu" = 12,
  "AF233305.1_Durio_zibethinus" = 13,
  "MF629762.1_Durio_zibethinus_isolate_M" = 14,
  "MT230663.1_Durio_zibethinus_voucher_DUR37" = 15,
  "MF629752.1_Durio_zibethinus_isolate_C" = 15,
  "MT230662.1_Durio_zibethinus_voucher_DUR36" = 15,
  "MT230654.1_Durio_zibethinus_voucher_DUR23" = 15,
  "MT230664.1_Durio_zibethinus_voucher_DUR39" = 15,
  "MT230650.1_Durio_zibethinus_voucher_DUR10" = 15,
  "MF629776.1_Durio_zibethinus_isolate_AA" = 15
)
pg_example <- function() stats::setNames(.pg_durio_hap[.pg_durio_map], names(.pg_durio_map))
pg_example_sets <- function() data.frame(ID = c("AF287710.1_Durio_zibethinus", "MT230650.1_Durio_zibethinus_voucher_DUR10", "MT230651.1_Durio_zibethinus_voucher_DUR12", "MT230652.1_Durio_zibethinus_voucher_DUR17", "MT230654.1_Durio_zibethinus_voucher_DUR23", "MT230655.1_Durio_zibethinus_voucher_DUR24", "MT230656.1_Durio_zibethinus_voucher_DUR33", "MT230662.1_Durio_zibethinus_voucher_DUR36", "MT230663.1_Durio_zibethinus_voucher_DUR37", "MT230664.1_Durio_zibethinus_voucher_DUR39", "Durio_zibethinus_cultivar_Duriat1", "Durio_zibethinus_cultivar_Duriat10", "Durio_zibethinus_cultivar_Kinoso1", "Durio_zibethinus_cultivar_Kinoso3", "ON815633.1_Durio_zibethinus_cultivar_Duriat", "ON817274.1_Durio_zibethinus_cultivar_Kinoso", "MF629752.1_Durio_zibethinus_isolate_C", "MF629753.1_Durio_zibethinus_isolate_D", "MF629754.1_Durio_zibethinus_isolate_E", "MF629755.1_Durio_zibethinus_isolate_F", "MF629756.1_Durio_zibethinus_isolate_G", "MF629757.1_Durio_zibethinus_isolate_H", "MF629759.1_Durio_zibethinus_isolate_J", "MF629761.1_Durio_zibethinus_isolate_L", "MF629762.1_Durio_zibethinus_isolate_M", "MF629763.1_Durio_zibethinus_isolate_N", "MF629764.1_Durio_zibethinus_isolate_O", "MF629766.1_Durio_zibethinus_isolate_Q", "MF629767.1_Durio_zibethinus_isolate_R", "MF629769.1_Durio_zibethinus_isolate_T", "MF629771.1_Durio_zibethinus_isolate_V", "MF629772.1_Durio_zibethinus_isolate_W", "MF629774.1_Durio_zibethinus_isolate_Y", "MF629775.1_Durio_zibethinus_isolate_Z", "MF629776.1_Durio_zibethinus_isolate_AA", "MF629777.1_Durio_zibethinus_isolate_BB", "MF629778.1_Durio_zibethinus_isolate_CC", "MF629779.1_Durio_zibethinus_isolate_DD", "OK052902.1_Durio_zibethinus_voucher_A4716_5.8S", "AF233305.1_Durio_zibethinus", "AF233304.1_Durio_zibethinus", "Durio_zibethinus_cultivar_Durian_Kunyit1", "Durio_zibethinus_cultivar_Durian_Kunyit7", "ON817266.1_Durio_zibethinus_cultivar_Durian_Kunyit", "MW186242.1_Durio_zibethinus_isolate_Ri-6", "MW316477.1_Durio_zibethinus_isolate_Ri-6_VCAQMN", "MW316599.1_Durio_zibethinus_isolate_Monthong_TG", "MW316621.1_Durio_zibethinus_isolate_Ri-6_TG", "MW316622.1_Durio_zibethinus_isolate_Ri-6_CT", "MW316637.1_Durio_zibethinus_isolate_Sua_Hat_Lep_CT", "MW316650.1_Durio_zibethinus_isolate_Sau_Huu", "MW316655.1_Durio_zibethinus_isolate_Chuong_Bo_TG", "MW318973.1_Durio_zibethinus_isolate_Chin_Hoa_CT", "MW318974.1_Durio_zibethinus_isolate_Monthong_CT", "MW318975.1_Durio_zibethinus_isolate_Monthong_VL", "MW318976.1_Durio_zibethinus_isolate_Chin_Hoa_BT", "MW318978.1_Durio_zibethinus_isolate_Kho_Qua_Xanh_VL", "MW318979.1_Durio_zibethinus_isolate_Bi_BT", "OR775700.1_Durio_zibethinus_voucher_ITS-Ri6-BT", "OR775701.1_Durio_zibethinus_voucher_ITS-MT-BT"),
  Species = c("Durio_zibethinus", "Durio_zibethinus", "Durio_zibethinus", "Durio_zibethinus", "Durio_zibethinus", "Durio_zibethinus", "Durio_zibethinus", "Durio_zibethinus", "Durio_zibethinus", "Durio_zibethinus", "Durio_zibethinus", "Durio_zibethinus", "Durio_zibethinus", "Durio_zibethinus", "Durio_zibethinus", "Durio_zibethinus", "Durio_zibethinus", "Durio_zibethinus", "Durio_zibethinus", "Durio_zibethinus", "Durio_zibethinus", "Durio_zibethinus", "Durio_zibethinus", "Durio_zibethinus", "Durio_zibethinus", "Durio_zibethinus", "Durio_zibethinus", "Durio_zibethinus", "Durio_zibethinus", "Durio_zibethinus", "Durio_zibethinus", "Durio_zibethinus", "Durio_zibethinus", "Durio_zibethinus", "Durio_zibethinus", "Durio_zibethinus", "Durio_zibethinus", "Durio_zibethinus", "Durio_zibethinus", "Durio_zibethinus", "Durio_zibethinus", "Durio_zibethinus", "Durio_zibethinus", "Durio_zibethinus", "Durio_zibethinus", "Durio_zibethinus", "Durio_zibethinus", "Durio_zibethinus", "Durio_zibethinus", "Durio_zibethinus", "Durio_zibethinus", "Durio_zibethinus", "Durio_zibethinus", "Durio_zibethinus", "Durio_zibethinus", "Durio_zibethinus", "Durio_zibethinus", "Durio_zibethinus", "Durio_zibethinus", "Durio_zibethinus"),
  Location = c("Borneo", "Borneo", "Borneo", "Borneo", "Borneo", "Borneo", "Borneo", "Borneo", "Borneo", "Borneo", "Siberut", "Siberut", "Siberut", "Siberut", "Siberut", "Siberut", "North Maluku", "North Maluku", "North Maluku", "North Maluku", "North Maluku", "North Maluku", "North Maluku", "North Maluku", "North Maluku", "North Maluku", "North Maluku", "North Maluku", "North Maluku", "North Maluku", "North Maluku", "North Maluku", "North Maluku", "North Maluku", "North Maluku", "North Maluku", "North Maluku", "North Maluku", "Peninsula Malaysia", "Peninsula Malaysia", "Sulawesi", "Sumatra", "Sumatra", "Sumatra", "Viet Nam", "Viet Nam", "Viet Nam", "Viet Nam", "Viet Nam", "Viet Nam", "Viet Nam", "Viet Nam", "Viet Nam", "Viet Nam", "Viet Nam", "Viet Nam", "Viet Nam", "Viet Nam", "Viet Nam", "Viet Nam"), stringsAsFactors = FALSE)

# which columns of a table of sequence sets are the sequence name, the population and the group (by the header; the first match wins)
pg_pick_set_cols <- function(nm) {
  nm <- as.character(nm); lc <- gsub("[^a-z0-9]+", "", tolower(trimws(nm)))
  find <- function(keys, skip = character()) { for (k in keys) { i <- which(lc == k & !nm %in% skip); if (length(i)) return(nm[i[1]]) }; NA_character_ }
  key <- find(c("sequence", "sequences", "sequencename", "seq", "seqid", "id", "sampleid", "sample", "accession", "accessionnumber", "name", "taxon", "label", "isolate")); if (is.na(key)) key <- nm[1]
  pop <- find(c("population", "pop", "location", "locality", "site", "samplingsite", "sampling", "area", "origin", "deme", "province", "place", "region", "country", "island"), key)
  if (is.na(pop)) { rest <- nm[!nm %in% key & !lc %in% c("species", "taxon", "latitude", "longitude", "lat", "lon", "long", "x", "y")]; pop <- if (length(rest)) rest[1] else setdiff(nm, key)[1] }
  grp <- find(c("group", "region", "island", "country", "cluster", "clade", "lineage", "area"), c(key, pop))
  list(key = key, pop = pop, grp = if (is.na(grp)) NA_character_ else grp)
}

# ---------------------------------------------------------------------------------------------------
# reading the sequence-set table: csv / tsv / txt / xlsx (no extra package needed for xlsx; readxl is used when it is installed)
# ---------------------------------------------------------------------------------------------------
.pg_xml_dec <- function(x) { x <- gsub("&lt;", "<", x, fixed = TRUE); x <- gsub("&gt;", ">", x, fixed = TRUE); x <- gsub("&quot;", "\"", x, fixed = TRUE)
  x <- gsub("&apos;", "'", x, fixed = TRUE); x <- gsub("&#10;", "\n", x, fixed = TRUE); gsub("&amp;", "&", x, fixed = TRUE) }
.pg_xml_text <- function(chunk) {                       # all <t>..</t> texts of one string item
  m <- regmatches(chunk, gregexpr("<t(?:\\s[^>]*)?>[^<]*</t>|<t\\s*/>", chunk, perl = TRUE))[[1]]
  .pg_xml_dec(paste(sub("</t>$", "", sub("^<t(?:\\s[^>]*)?>", "", m, perl = TRUE)), collapse = ""))
}
.pg_col_idx <- function(ref) { l <- strsplit(sub("[0-9]+$", "", toupper(ref)), "")[[1]]; sum((utf8ToInt(paste(l, collapse = "")) - 64L) * 26L^(rev(seq_along(l)) - 1L)) }
# Returns a list of data frames (one per sheet, first row = header, all columns as text)
pg_read_xlsx <- function(path) {
  td <- tempfile("xlsx"); dir.create(td); on.exit(unlink(td, recursive = TRUE), add = TRUE)
  utils::unzip(path, exdir = td)
  rd <- function(f) { p <- file.path(td, f); if (!file.exists(p)) return(""); x <- readBin(p, "raw", file.info(p)$size); iconv(rawToChar(x), "UTF-8", "UTF-8", sub = "") }
  wb <- rd("xl/workbook.xml"); rels <- rd("xl/_rels/workbook.xml.rels")
  sh <- regmatches(wb, gregexpr("<sheet\\b[^>]*>", wb))[[1]]
  attr <- function(tag, a) { m <- regmatches(tag, regexpr(sprintf("\\b%s=\"[^\"]*\"", a), tag)); if (length(m)) .pg_xml_dec(sub(sprintf("^%s=\"", a), "", sub("\"$", "", m))) else NA_character_ }
  rl <- regmatches(rels, gregexpr("<Relationship\\b[^>]*>", rels))[[1]]
  rid_target <- stats::setNames(vapply(rl, attr, "", "Target", USE.NAMES = FALSE), vapply(rl, attr, "", "Id", USE.NAMES = FALSE))
  ss <- rd("xl/sharedStrings.xml")
  sst <- if (nzchar(ss)) vapply(regmatches(ss, gregexpr("<si\\b.*?</si>|<si\\s*/>", ss, perl = TRUE))[[1]], .pg_xml_text, "", USE.NAMES = FALSE) else character()
  out <- list()
  for (k in seq_along(sh)) {
    nm <- attr(sh[k], "name"); rid <- attr(sh[k], "r:id"); tg <- rid_target[[rid]]
    f <- if (is.null(tg) || is.na(tg)) sprintf("xl/worksheets/sheet%d.xml", k) else sub("^/", "", if (startsWith(tg, "/")) tg else paste0("xl/", tg))
    x <- rd(f); if (!nzchar(x)) next
    cells <- regmatches(x, gregexpr("<c\\b[^>]*?/>|<c\\b[^>]*?>.*?</c>", x, perl = TRUE))[[1]]
    if (!length(cells)) { out[[nm]] <- data.frame(); next }
    ref <- vapply(cells, attr, "", "r", USE.NAMES = FALSE); typ <- vapply(cells, attr, "", "t", USE.NAMES = FALSE)
    val <- vapply(seq_along(cells), function(i) {
      ce <- cells[i]; if (grepl("/>$", ce) && !grepl("</c>", ce, fixed = TRUE)) return(NA_character_)
      if (identical(typ[i], "inlineStr")) return(.pg_xml_text(ce))
      v <- regmatches(ce, regexpr("<v>[^<]*</v>", ce)); if (!length(v)) return(NA_character_)
      v <- .pg_xml_dec(sub("</v>$", "", sub("^<v>", "", v)))
      if (identical(typ[i], "s")) sst[as.integer(v) + 1L] else if (identical(typ[i], "b")) c("FALSE", "TRUE")[as.integer(v) + 1L] else v
    }, "")
    ci <- vapply(ref, .pg_col_idx, 1); ri <- as.integer(sub("^[A-Za-z]+", "", ref))
    m <- matrix(NA_character_, max(ri), max(ci)); m[cbind(ri, ci)] <- val
    m <- m[rowSums(!is.na(m) & nzchar(trimws(m))) > 0, , drop = FALSE]
    if (nrow(m) < 1) { out[[nm]] <- data.frame(); next }
    hd <- trimws(ifelse(is.na(m[1, ]), "", m[1, ])); keep <- nzchar(hd)
    out[[nm]] <- as.data.frame(m[-1, keep, drop = FALSE], stringsAsFactors = FALSE, optional = TRUE) |> stats::setNames(hd[keep])
  }
  out
}
# Reads the table of sequence sets from csv / tsv / txt / xlsx / xls; returns a data frame or a character message (the reason it failed)
pg_read_sets <- function(path, name = path) {
  ext <- tolower(tools::file_ext(name))
  pick <- function(L) {                                   # the first sheet that has a Sequence or Population column, otherwise the first one with two columns
    L <- L[vapply(L, function(d) ncol(d) >= 2 && nrow(d) >= 1, TRUE)]; if (!length(L)) return(NULL)
    has <- vapply(L, function(d) any(tolower(trimws(names(d))) %in% c("sequence", "population")), TRUE)
    if (any(has)) L[[which(has)[1]]] else L[[1]]
  }
  if (ext %in% c("xlsx", "xlsm", "xls")) {
    if (requireNamespace("readxl", quietly = TRUE)) {
      L <- tryCatch({ sn <- readxl::excel_sheets(path); stats::setNames(lapply(sn, function(s) as.data.frame(readxl::read_excel(path, sheet = s, col_types = "text", .name_repair = "minimal"), stringsAsFactors = FALSE)), sn) }, error = function(e) NULL)
    } else L <- NULL
    if (is.null(L) && ext != "xls") L <- tryCatch(pg_read_xlsx(path), error = function(e) NULL)
    if (is.null(L)) return(if (ext == "xls") "An old .xls file cannot be read here. In Excel use File > Save As > Excel Workbook (.xlsx) or CSV and load that file." else "This Excel file could not be read. Save it again as .xlsx or CSV.")
    d <- pick(L); if (is.null(d)) return("The Excel file has no sheet with at least two columns (Sequence, Population).")
    return(d)
  }
  rd <- list(function() utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE, fileEncoding = "UTF-8-BOM", colClasses = "character"),
             function() utils::read.csv2(path, stringsAsFactors = FALSE, check.names = FALSE, fileEncoding = "UTF-8-BOM", colClasses = "character"),
             function() utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE, fileEncoding = "UTF-8-BOM", colClasses = "character"))
  for (f in rd) { d <- tryCatch(f(), error = function(e) NULL); if (!is.null(d) && ncol(d) >= 2) return(d) }
  "The file needs at least two columns: Sequence and Population (and optionally Group). Separator: comma, semicolon or tab."
}

# ---------------------------------------------------------------------------------------------------
# pairwise differences with the gap options of DnaSP
# ---------------------------------------------------------------------------------------------------
# gap = "excluded" / "pairwise" / "fifth" as in pg_prepare. Returns d (differences), nv (sites compared), p (differences per site).
pg_pairdiff <- function(M, gap = "excluded") {
  states <- if (identical(gap, "fifth")) c("A", "C", "G", "T", "-") else c("A", "C", "G", "T")
  valid <- matrix(M %in% states, nrow(M))
  nv <- tcrossprod(valid * 1)
  same <- 0
  for (s in states) same <- same + tcrossprod(((M == s) & valid) * 1)
  d <- nv - same
  p <- ifelse(nv > 0, d / nv, NA_real_)
  dimnames(d) <- dimnames(p) <- dimnames(nv) <- list(rownames(M), rownames(M))
  diag(d) <- 0; diag(p) <- 0
  list(d = d, p = p, nv = nv)
}

.pg_jc <- function(p) -0.75 * log(1 - 4 / 3 * pmin(p, 0.7499))
.pg_hd <- function(h) { n <- length(h); if (n < 2) return(NA_real_); q <- as.integer(table(h)) / n; n / (n - 1) * (1 - sum(q^2)) }
.pg_pair_mean <- function(Dm, ia, ib = NULL) {
  if (is.null(ib)) { if (length(ia) < 2) return(NA_real_); x <- Dm[ia, ia, drop = FALSE]; mean(x[upper.tri(x)]) }
  else mean(Dm[ia, ib, drop = FALSE])
}
# Nm from an Fst-like statistic; genome: haploid (2), diploid autosome (4), X-linked (3), Y-linked (1)
pg_nm <- function(f, genome = "diploid") {
  mult <- switch(genome, haploid = 2, diploid = 4, X = 3, Y = 1, 4)
  ifelse(is.na(f), NA_real_, ifelse(f > 0, (1 - f) / (mult * f), Inf))
}

# ---------------------------------------------------------------------------------------------------
# "Gene flow and genetic differentiation" in the style of DnaSP
#   Hs / Hst  Hudson, Boos & Kaplan 1992 (weights n - 2 when every sample has more than 2 sequences)
#   Ks / Kst  weights n;  Fst  Hudson, Slatkin & Maddison 1992 = 1 - mean(K within) / mean(K between)
#   Nst       the same with Jukes-Cantor distances (Lynch & Crease 1990);  Snn  Hudson 2000
#   Chi-square of the haplotype x population table
# ---------------------------------------------------------------------------------------------------
pg_geneflow <- function(M, pop, gap = "excluded", genome = "diploid", nperm = 0, seed = NULL, progress = NULL) {
  pop <- as.character(pop); pops <- unique(pop); P <- length(pops); n <- nrow(M)
  if (P < 2) stop("At least 2 populations are needed.")
  if (!is.null(seed) && !is.na(seed)) set.seed(seed)
  pdm <- pg_pairdiff(M, gap); D <- pdm$d; Pm <- pdm$p; J <- .pg_jc(Pm); diag(J) <- 0
  acgt <- M == "A" | M == "C" | M == "G" | M == "T"
  hc <- colSums(!acgt) == 0
  if (!any(hc)) stop("No site is A/C/G/T in every sequence, so haplotypes cannot be defined.")
  hap <- pg_haplotypes(M[, hc, drop = FALSE])$hap
  # nearest neighbours of every sequence (fixed by the data; only the labels are permuted)
  nnl <- lapply(seq_len(n), function(i) { d <- Pm[i, ]; d[i] <- Inf; which(d == min(d, na.rm = TRUE)) })
  snn <- function(lab) mean(vapply(seq_len(n), function(i) mean(lab[nnl[[i]]] == lab[i]), 0))
  core <- function(lab) {
    ids <- split(seq_len(n), factor(lab, levels = unique(lab))); ni <- vapply(ids, length, 1L)
    Hd <- vapply(ids, function(i) .pg_hd(hap[i]), 0); K <- vapply(ids, function(i) .pg_pair_mean(D, i), 0)
    okH <- !is.na(Hd); okK <- !is.na(K)
    wH <- if (all(ni > 2)) ni - 2 else ni
    Hs <- sum((wH * Hd)[okH]) / sum(wH[okH]); Ht <- .pg_hd(hap)
    Ks <- sum((ni * K)[okK]) / sum(ni[okK]); Kt <- .pg_pair_mean(D, seq_len(n))
    tb <- table(hap, factor(lab, levels = names(ids))); E <- outer(rowSums(tb), colSums(tb)) / n
    c(chi2 = sum((tb - E)^2 / E), Hs = Hs, Ht = Ht, Hst = 1 - Hs / Ht, Ks = Ks, Kt = Kt, Kst = 1 - Ks / Kt, Snn = snn(lab))
  }
  obs <- core(pop)
  idx <- lapply(pops, function(p) which(pop == p)); names(idx) <- pops; ni <- vapply(idx, length, 1L)
  nh <- nlevels(droplevels(hap)); df <- (nh - 1) * (P - 1)
  perm <- NULL
  if (nperm > 0) {
    perm <- matrix(NA_real_, nperm, 4, dimnames = list(NULL, c("chi2", "Hst", "Kst", "Snn")))
    for (r in seq_len(nperm)) {
      perm[r, ] <- core(sample(pop))[c("chi2", "Hst", "Kst", "Snn")]
      if (!is.null(progress) && r %% 25 == 0) progress(r / nperm)
    }
  }
  pp <- function(nm) if (is.null(perm) || is.na(obs[[nm]])) NA_real_ else (sum(perm[, nm] >= obs[[nm]] - 1e-12, na.rm = TRUE) + 1) / (sum(!is.na(perm[, nm])) + 1)
  # F-statistics from pairwise means
  Kw <- vapply(idx, function(i) .pg_pair_mean(D, i), 0); Jw <- vapply(idx, function(i) .pg_pair_mean(J, i), 0)
  Hdw <- vapply(idx, function(i) .pg_hd(hap[i]), 0)
  pr <- t(utils::combn(P, 2))
  Kb <- apply(pr, 1, function(q) .pg_pair_mean(D, idx[[q[1]]], idx[[q[2]]]))
  Jb <- apply(pr, 1, function(q) .pg_pair_mean(J, idx[[q[1]]], idx[[q[2]]]))
  Hb <- apply(pr, 1, function(q) { a <- table(factor(hap[idx[[q[1]]]], levels = levels(hap))); b <- table(factor(hap[idx[[q[2]]]], levels = levels(hap)))
                                  1 - sum(a * b) / (length(idx[[q[1]]]) * length(idx[[q[2]]])) })
  Fst <- 1 - mean(Kw) / mean(Kb); Nst <- 1 - mean(Jw) / mean(Jb); Gst <- 1 - mean(Hdw) / mean(Hb)
  # per population and total
  states <- if (identical(gap, "fifth")) c("A", "C", "G", "T", "-") else c("A", "C", "G", "T")
  Sof <- function(i) { x <- M[i, , drop = FALSE]; sum(vapply(seq_len(ncol(x)), function(j) { v <- x[, j]; length(unique(v[v %in% states])) > 1 }, TRUE)) }
  one <- function(i, label) {
    data.frame(Population = label, n = length(i), S = Sof(i), h = nlevels(droplevels(hap[i])), Hd = .pg_hd(hap[i]),
               K = .pg_pair_mean(D, i), Pi = .pg_pair_mean(Pm, i), PiJC = .pg_pair_mean(J, i), stringsAsFactors = FALSE)
  }
  popt <- do.call(rbind, c(lapply(pops, function(p) one(idx[[p]], p)), list(one(seq_len(n), "Total"))))
  # pairwise table
  prt <- do.call(rbind, lapply(seq_len(nrow(pr)), function(r) {
    a <- idx[[pr[r, 1]]]; b <- idx[[pr[r, 2]]]; ab <- c(a, b); na <- length(a); nb <- length(b)
    w <- if (na > 2 && nb > 2) c(na - 2, nb - 2) else c(na, nb)
    hs <- sum(w * c(.pg_hd(hap[a]), .pg_hd(hap[b]))) / sum(w); ht <- .pg_hd(hap[ab])
    ka <- .pg_pair_mean(D, a); kb <- .pg_pair_mean(D, b); ks <- (na * ka + nb * kb) / (na + nb); kt <- .pg_pair_mean(D, ab)
    kxy <- Kb[r]; fst <- 1 - mean(c(ka, kb)) / kxy; nst <- 1 - mean(c(.pg_pair_mean(J, a), .pg_pair_mean(J, b))) / Jb[r]
    gst <- 1 - mean(c(.pg_hd(hap[a]), .pg_hd(hap[b]))) / Hb[r]
    data.frame(Population1 = pops[pr[r, 1]], Population2 = pops[pr[r, 2]], Hs = hs, Hst = 1 - hs / ht, Ks = ks, Kst = 1 - ks / kt, Kxy = kxy,
               Gst = gst, Fst = fst, Nm_Fst = pg_nm(fst, genome), Nst = nst, stringsAsFactors = FALSE)
  }))
  tests <- data.frame(
    Statistic = c("Chi-square (haplotype table)", "Hs (Hudson et al. 1992)", "Hst", "Ks", "Kst", "Snn (Hudson 2000)", "Gst (Hudson et al. 1992 form)",
                  "Fst (Hudson, Slatkin & Maddison 1992)", "Nst (Lynch & Crease 1990, JC)"),
    Value = c(obs[["chi2"]], obs[["Hs"]], obs[["Hst"]], obs[["Ks"]], obs[["Kst"]], obs[["Snn"]], Gst, Fst, Nst),
    P_chi2_distribution = c(stats::pchisq(obs[["chi2"]], df, lower.tail = FALSE), rep(NA_real_, 8)),
    P_permutation = c(pp("chi2"), NA, pp("Hst"), NA, pp("Kst"), pp("Snn"), NA, NA, NA),
    Nm = c(NA, NA, NA, NA, NA, NA, pg_nm(Gst, genome), pg_nm(Fst, genome), pg_nm(Nst, genome)), stringsAsFactors = FALSE)
  list(pops = popt, tests = tests, pairwise = prt, df = df, Ht = obs[["Ht"]], Kt = obs[["Kt"]], haplotypes = nh, nperm = nperm,
       sites = ncol(M), gap = gap, genome = genome)
}

# Arlequin-like estimates of theta from one sample: Hom (homozygosity), S, k (number of haplotypes) and Pi.
pg_thetas <- function(n, S, hap_n, K) {
  a1 <- sum(1 / seq_len(n - 1)); p <- hap_n / n; k <- length(hap_n)
  Fh <- sum(p^2)
  # Theta(Hom) as Arlequin (Zouros 1979; Chakraborty & Weiss 1991): H = (n sum p^2 - 1) / (n - 1) is the unbiased homozygosity and theta solves
  # (1 - H) / H = theta * (1 + 2 (1 + theta) / ((2 + theta) (3 + theta))). Not defined when all the genes are different (H <= 0).
  Hu <- if (n > 1) (n * Fh - 1) / (n - 1) else NA_real_
  thH <- if (!is.finite(Hu) || Hu <= 1e-12) NA_real_ else if (Hu >= 1) 0 else {
    tg <- (1 - Hu) / Hu
    tryCatch(stats::uniroot(function(th) th * (1 + 2 * (1 + th) / ((2 + th) * (3 + th))) - tg, c(0, 1e6), tol = 1e-12)$root, error = function(e) NA_real_)
  }
  c(theta_Hom = thH, theta_S = S / a1,
    theta_k = if (k > 1 && k < n) tryCatch(stats::uniroot(function(th) sum(th / (th + 0:(n - 1))) - k, c(1e-8, 1e5))$root, error = function(e) NA_real_) else NA_real_,
    theta_Pi = K)
}

# ---------------------------------------------------------------------------------------------------
# haplotype distribution in the style of DnaSP ("Hap_8: 2  [20 24]")
# ---------------------------------------------------------------------------------------------------
pg_hap_text <- function(hap, nms) {
  rng <- function(v) {
    v <- sort(v); g <- split(v, cumsum(c(1, diff(v) != 1)))
    paste(vapply(g, function(x) if (length(x) > 1) paste0(min(x), "-", max(x)) else as.character(x), ""), collapse = " ")
  }
  idx <- split(seq_along(hap), hap)
  idx <- idx[vapply(idx, length, 1L) > 0]
  tab <- as.integer(table(hap)); n <- length(hap); p <- tab[tab > 0] / n
  Hd <- n / (n - 1) * (1 - sum(p^2))
  p2 <- sum(p^2); p3 <- sum(p^3)
  vH <- 2 / (n * (n - 1)) * (2 * (n - 2) * (p3 - p2^2) + p2 - p2^2)
  c(sprintf("Number of sequences: %d", n), sprintf("Number of haplotypes, h: %d", length(idx)),
    sprintf("Haplotype (gene) diversity, Hd: %.5f", Hd), sprintf("Variance of Hd: %.5f", vH),
    sprintf("Standard deviation of Hd: %.5f", sqrt(max(vH, 0))), "",
    "Haplotype data (sequence numbers)", vapply(names(idx), function(h) sprintf("%s: %d  [%s]", h, length(idx[[h]]), rng(idx[[h]])), ""),
    "", "Haplotype data (sequence names)",
    vapply(names(idx), function(h) sprintf("%s: %d  [%s]", h, length(idx[[h]]), paste(nms[idx[[h]]], collapse = " ")), ""))
}

# ---------------------------------------------------------------------------------------------------
# minimum spanning network among haplotypes (all edges that belong to some minimum spanning tree: Excoffier & Smouse 1994)
# ---------------------------------------------------------------------------------------------------
pg_hamming <- function(S) {
  k <- nrow(S); D <- matrix(0L, k, k)
  for (i in seq_len(k)) D[i, ] <- as.integer(rowSums(sweep(S, 2, S[i, ], "!=")))
  D
}
pg_net_msn <- function(hap_seq) {
  k <- length(hap_seq)
  if (k < 2) return(data.frame(a = integer(0), b = integer(0), steps = integer(0)))
  S <- do.call(rbind, strsplit(hap_seq, "")); D <- pg_hamming(S)
  inT <- c(TRUE, rep(FALSE, k - 1)); best <- D[1, ]; parent <- rep(1L, k); tree <- matrix(0L, 0, 2)
  for (s in seq_len(k - 1)) {
    cand <- which(!inT); j <- cand[which.min(best[cand])]
    tree <- rbind(tree, c(parent[j], j)); inT[j] <- TRUE
    upd <- which(!inT & D[j, ] < best); best[upd] <- D[j, upd]; parent[upd] <- j
  }
  adj <- vector("list", k)
  for (e in seq_len(nrow(tree))) { a <- tree[e, 1]; b <- tree[e, 2]; adj[[a]] <- c(adj[[a]], b); adj[[b]] <- c(adj[[b]], a) }
  B <- matrix(0L, k, k)                                   # largest edge on the tree path between two haplotypes
  for (r in seq_len(k)) {
    seen <- rep(FALSE, k); seen[r] <- TRUE; stack <- r
    while (length(stack)) {
      u <- stack[length(stack)]; stack <- stack[-length(stack)]
      for (v in adj[[u]]) if (!seen[v]) { seen[v] <- TRUE; B[r, v] <- max(B[r, u], D[u, v]); stack <- c(stack, v) }
    }
  }
  ii <- which(upper.tri(D) & D <= B, arr.ind = TRUE)
  data.frame(a = ii[, 1], b = ii[, 2], steps = D[ii], stringsAsFactors = FALSE)
}

# ---------------------------------------------------------------------------------------------------
# reports: self-contained HTML (Arlequin-like) and plain text (DnaSP-like)
# ---------------------------------------------------------------------------------------------------
.pg_esc <- function(x) { x <- gsub("&", "&amp;", x, fixed = TRUE); x <- gsub("<", "&lt;", x, fixed = TRUE); gsub(">", "&gt;", x, fixed = TRUE) }
.pg_html_table <- function(d, digits = 5) {
  if (is.null(d) || !nrow(d)) return("<p><i>(no data)</i></p>")
  fm <- function(x) if (is.numeric(x)) ifelse(is.na(x), "", ifelse(is.finite(x), format(round(x, digits), trim = TRUE, scientific = FALSE), as.character(x))) else .pg_esc(as.character(x))
  hd <- paste0("<th>", .pg_esc(names(d)), "</th>", collapse = "")
  rows <- vapply(seq_len(nrow(d)), function(i) paste0("<tr>", paste0("<td>", vapply(d, function(col) fm(col)[i], ""), "</td>", collapse = ""), "</tr>"), "")
  paste0("<table><thead><tr>", hd, "</tr></thead><tbody>", paste(rows, collapse = ""), "</tbody></table>")
}
# sections: list of list(h = heading, tbl = data.frame, pre = character lines, p = paragraph)
pg_html <- function(title, sections, subtitle = "") {
  css <- paste("body{font-family:Segoe UI,Arial,sans-serif;margin:24px auto;max-width:1100px;color:#1f2d3d;background:#fff}",
               "h1{background:#8a8f95;color:#fff;padding:10px 16px;margin:0 0 6px 0;font-size:22px}",
               "h2{border-bottom:2px solid #8a8f95;margin-top:28px;padding-bottom:4px;font-size:17px;color:#33414f}",
               "table{border-collapse:collapse;margin:8px 0;font-size:13px}th,td{border:1px solid #c9ced6;padding:3px 9px;text-align:right}",
               "th{background:#eef1f5;text-align:center}td:first-child,th:first-child{text-align:left}",
               "pre{background:#f6f7f9;border:1px solid #dfe3e8;padding:10px;font-size:12px;overflow:auto}",
               "p.note{color:#52606d;font-size:12.5px}", sep = "\n")
  body <- vapply(sections, function(s) paste0(
    if (!is.null(s$h)) paste0("<h2>", .pg_esc(s$h), "</h2>") else "",
    if (!is.null(s[["p"]])) paste0("<p class='note'>", .pg_esc(s[["p"]]), "</p>") else "",       # [[ ]]: $p would partly match $pre
    if (!is.null(s$tbl)) .pg_html_table(s$tbl, if (is.null(s$digits)) 5 else s$digits) else "",
    if (!is.null(s$pre)) paste0("<pre>", .pg_esc(paste(s$pre, collapse = "\n")), "</pre>") else ""), "")
  paste0("<!DOCTYPE html><html><head><meta charset='utf-8'><title>", .pg_esc(title), "</title><style>", css, "</style></head><body><h1>",
         .pg_esc(title), "</h1><p class='note'>", .pg_esc(subtitle), "</p>", paste(body, collapse = "\n"), "</body></html>")
}

# ---------------------------------------------------------------------------------------------------
# Ewens-Watterson neutrality test (infinite alleles model; haplotypes are the alleles)
#   F = homozygosity = sum of the squared haplotype frequencies. The null distribution is sampled from the Ewens sampling
#   formula conditional on the observed number of haplotypes k (theta is chosen so that E[K] = k; samples with K != k are discarded).
#   p_lower = P(F_sim <= F_obs) (too even: balancing selection), p_upper = P(F_sim >= F_obs) (too uneven: directional selection / expansion).
# ---------------------------------------------------------------------------------------------------
pg_ew_test <- function(counts, nsim = 1000, max_try = 300) {
  counts <- counts[counts > 0]; n <- sum(counts); k <- length(counts)
  F_obs <- sum((counts / n)^2)
  out <- list(n = n, k = k, F = F_obs, mean = NA_real_, p_lower = NA_real_, p_upper = NA_real_, nsim = 0L)
  if (nsim <= 0 || k < 2 || k >= n) return(out)
  th <- tryCatch(stats::uniroot(function(th) sum(th / (th + 0:(n - 1))) - k, c(1e-8, 1e5))$root, error = function(e) NA_real_)
  if (is.na(th)) return(out)
  sims <- numeric(nsim); got <- 0L; tries <- 0
  while (got < nsim && tries < nsim * max_try) {
    tries <- tries + 1
    cnt <- integer(0); bad <- FALSE
    for (i in 0:(n - 1)) {
      if (stats::runif(1) < th / (th + i)) {
        cnt <- c(cnt, 1L)
        if (length(cnt) > k) { bad <- TRUE; break }
      } else { j <- sample.int(length(cnt), 1, prob = cnt); cnt[j] <- cnt[j] + 1L }
    }
    if (!bad && length(cnt) == k) { got <- got + 1L; sims[got] <- sum((cnt / n)^2) }
  }
  if (got == 0L) return(out)
  sims <- sims[seq_len(got)]
  out$mean <- mean(sims); out$p_lower <- mean(sims <= F_obs + 1e-12); out$p_upper <- mean(sims >= F_obs - 1e-12); out$nsim <- got
  out
}

# ---------------------------------------------------------------------------------------------------
# Exact test of population differentiation (Raymond & Rousset 1995): Markov chain on the haplotype x population table
# with fixed margins. p = proportion of visited tables whose probability is <= that of the observed table.
# ---------------------------------------------------------------------------------------------------
pg_exact_test <- function(tab, steps = 100000, demem = 10000, seed = NULL) {
  tab <- as.matrix(tab); tab <- tab[rowSums(tab) > 0, colSums(tab) > 0, drop = FALSE]
  nr <- nrow(tab); nc <- ncol(tab)
  if (nr < 2 || nc < 2 || steps <= 0) return(list(p = NA_real_, se = NA_real_, steps = 0))
  if (!is.null(seed) && !is.na(seed)) set.seed(seed)
  cur <- tab; lp_obs <- -sum(lgamma(tab + 1)); lp <- lp_obs
  le <- integer(steps); last <- 0L
  for (s in seq_len(demem + steps)) {
    r <- sample.int(nr, 2); cc <- sample.int(nc, 2); a <- r[1]; b <- r[2]; j <- cc[1]; k <- cc[2]
    if (stats::runif(1) < 0.5) {                      # add one to (a,j) and (b,k), remove one from (a,k) and (b,j)
      if (cur[a, k] >= 1 && cur[b, j] >= 1) {
        ratio <- cur[a, k] * cur[b, j] / ((cur[a, j] + 1) * (cur[b, k] + 1))
        if (stats::runif(1) < ratio) { cur[a, j] <- cur[a, j] + 1; cur[b, k] <- cur[b, k] + 1; cur[a, k] <- cur[a, k] - 1; cur[b, j] <- cur[b, j] - 1; lp <- lp + log(ratio) }
      }
    } else {                                           # the opposite move
      if (cur[a, j] >= 1 && cur[b, k] >= 1) {
        ratio <- cur[a, j] * cur[b, k] / ((cur[a, k] + 1) * (cur[b, j] + 1))
        if (stats::runif(1) < ratio) { cur[a, k] <- cur[a, k] + 1; cur[b, j] <- cur[b, j] + 1; cur[a, j] <- cur[a, j] - 1; cur[b, k] <- cur[b, k] - 1; lp <- lp + log(ratio) }
      }
    }
    if (s > demem) le[s - demem] <- as.integer(lp <= lp_obs + 1e-8)
  }
  bm <- vapply(split(le, cut(seq_along(le), 10, labels = FALSE)), mean, 0)       # 10 batches for the standard error
  list(p = mean(le), se = stats::sd(bm) / sqrt(length(bm)), steps = steps)
}

# p-values of the exact test for every pair of populations (columns of the table)
pg_exact_pairwise <- function(tab, steps = 20000, demem = 2000, seed = NULL, progress = NULL) {
  P <- ncol(tab); pv <- matrix(NA_real_, P, P, dimnames = list(colnames(tab), colnames(tab)))
  tot <- P * (P - 1) / 2; done <- 0
  for (i in seq_len(P - 1)) for (j in (i + 1):P) {
    pv[i, j] <- pv[j, i] <- pg_exact_test(tab[, c(i, j), drop = FALSE], steps, demem, seed)$p
    done <- done + 1; if (!is.null(progress)) progress(done / tot)
  }
  pv
}

# ---------------------------------------------------------------------------------------------------
# Structure: model-based clustering of the sequences into K clusters (STRUCTURE / ADMIXTURE style bar plots)
#   Every variable site is coded as an allele (A/C/G/T) and the sites are treated as independent, as the "no linkage" model
#   of STRUCTURE does. Two models, both fitted by EM with several random starts:
#     admixture : every sequence has its own mixture of the K clusters, Q[i, k]  (bars with several colours)
#     mixture   : every sequence belongs to ONE cluster (Q[i, k] = probability of the cluster; usually almost 0 or 1);
#                 this suits one non-recombining locus (mtDNA, one gene), where the sites are not independent.
#   With linked sites the likelihood is over-confident: read the result as a description of the groups, not as a test.
# ---------------------------------------------------------------------------------------------------

# Allele indicator matrix of the variable sites of M (n x sites, letters A/C/G/T).
# Returns list(G = n x T matrix of 0/1, sid = site number of every column, S = sites used, cols = column numbers in M).
pg_str_data <- function(M) {
  var <- which(apply(M, 2, function(x) length(unique(x)) > 1))
  if (length(var) < 2) stop("Fewer than 2 variable sites: there is nothing to cluster. Use more sequences or a longer alignment.")
  blocks <- lapply(var, function(j) { a <- sort(unique(M[, j])); outer(M[, j], a, "==") * 1 })
  nb <- vapply(blocks, ncol, 1L); off <- c(0L, cumsum(nb))[seq_along(nb)]
  C <- vapply(seq_along(blocks), function(s) as.integer(max.col(blocks[[s]], ties.method = "first")) + off[s], integer(nrow(M)))   # n x S: column of G of the allele seen
  list(G = do.call(cbind, blocks), C = matrix(C, nrow(M)), sid = rep(seq_along(var), nb), S = length(var), cols = var)
}

.pg_str_normP <- function(P, sid) {                  # every site of every cluster sums to 1 (and never reaches 0)
  P <- pmax(P, 1e-8)
  s <- rowsum(t(P), sid, reorder = FALSE)            # sites x K
  P / t(s[sid, , drop = FALSE])
}

# One EM run. Returns list(Q, P, pi, lnL, iter).
.pg_str_em <- function(G, sid, K, model = "admix", maxit = 300, tol = 1e-7) {
  n <- nrow(G); Tn <- ncol(G)
  P <- .pg_str_normP(matrix(stats::runif(K * Tn), K, Tn), sid)
  Q <- matrix(stats::rgamma(n * K, 1), n, K); Q <- Q / rowSums(Q)
  pi_k <- rep(1 / K, K); ll_old <- -Inf; it <- 0L
  for (it in seq_len(maxit)) {
    if (model == "admix") {
      mix <- pmax(Q %*% P, 1e-300)
      ll <- sum(G * log(mix))
      W <- G / mix
      Qn <- Q * (W %*% t(P)); Qn <- Qn / rowSums(Qn)
      P <- .pg_str_normP(P * (t(Q) %*% W), sid); Q <- Qn
    } else {
      lp <- G %*% t(log(P)) + matrix(log(pmax(pi_k, 1e-300)), n, K, byrow = TRUE)
      m <- apply(lp, 1, max); R <- exp(lp - m); sR <- rowSums(R)
      ll <- sum(m + log(sR)); R <- R / sR
      P <- .pg_str_normP(t(R) %*% G + 1e-3, sid); pi_k <- colMeans(R); Q <- R
    }
    if (is.finite(ll_old) && abs(ll - ll_old) < tol * abs(ll)) break
    ll_old <- ll
  }
  ll <- if (model == "admix") sum(G * log(pmax(Q %*% P, 1e-300))) else {
    lp <- G %*% t(log(P)) + matrix(log(pmax(pi_k, 1e-300)), n, K, byrow = TRUE); m <- apply(lp, 1, max); sum(m + log(rowSums(exp(lp - m))))
  }
  list(Q = Q, P = P, pi = pi_k, lnL = ll, iter = it)
}


# One run of the Gibbs sampler of Pritchard, Stephens & Donnelly (2000), the algorithm of the program STRUCTURE, for haploid sequences:
# every variable site is a locus, each sequence has one allele per locus.
#   admixture model: z(i,l) ~ q(i); q(i) ~ Dirichlet(alpha, ..., alpha); p(k,l) ~ Dirichlet(1, ..., 1); alpha is fixed or updated by Metropolis-Hastings
#   no-admixture model: z(i) uniform over 1..K; p(k,l) ~ Dirichlet(1, ..., 1)
# Returns the posterior mean of Q (average of the sweeps after the burn-in), the mean ln likelihood of those sweeps and the mean alpha.
.pg_str_gibbs <- function(C, Tn, sid, K, model = "admix", burnin = 2000, nrep = 3000, alpha = 1, infer_alpha = TRUE, tick = NULL) {
  n <- nrow(C); S <- ncol(C); No <- n * S; cv <- as.vector(C)
  P <- .pg_str_normP(matrix(stats::rgamma(K * Tn, 1), K, Tn), sid)
  Q <- matrix(1 / K, n, K); zi <- sample.int(K, n, TRUE)
  Qsum <- matrix(0, n, K); ns <- 0L; llv <- numeric(nrep); alv <- numeric(nrep)
  rdir <- function(shape) { g <- matrix(pmax(stats::rgamma(length(shape), shape = as.vector(shape)), 1e-300), nrow(shape)); g / rowSums(g) }
  for (sw in seq_len(burnin + nrep)) {
    if (K == 1L) { Q[] <- 1; z <- rep(1L, No); ll <- sum(log(P[1, cv])) }
    if (model == "admix" && K > 1L) {
      cum <- matrix(0, No, K); run <- 0
      for (k in seq_len(K)) { run <- run + rep(Q[, k], S) * P[k, cv]; cum[, k] <- run }
      ll <- sum(log(cum[, K])); u <- stats::runif(No) * cum[, K]
      z <- rep(1L, No); for (k in seq_len(K - 1L)) z <- z + (cum[, k] < u)
      Zm <- matrix(z, n, S)
      Cnt <- vapply(seq_len(K), function(k) rowSums(Zm == k), numeric(n))
      Q <- rdir(alpha + Cnt)
      if (infer_alpha) {
        lq <- sum(log(pmax(Q, 1e-300))); f <- function(a) n * (lgamma(K * a) - K * lgamma(a)) + (a - 1) * lq
        a2 <- alpha + stats::rnorm(1, 0, 0.05)
        if (a2 > 0.001 && a2 < 10 && log(stats::runif(1)) < f(a2) - f(alpha)) alpha <- a2
      }
    } else if (model != "admix" && K > 1L) {
      lP <- log(P); LP <- vapply(seq_len(K), function(k) rowSums(matrix(lP[k, cv], n, S)), numeric(n))
      m <- apply(LP, 1, max); R <- exp(LP - m); cs <- rowSums(R); ll <- sum(m + log(cs)) - n * log(K)
      cum <- R; for (k in seq_len(K)[-1]) cum[, k] <- cum[, k - 1L] + R[, k]
      zi <- 1L + rowSums(cum[, -K, drop = FALSE] < stats::runif(n) * cs)
      Q <- matrix(0, n, K); Q[cbind(seq_len(n), zi)] <- 1; z <- rep(zi, S)
    }
    Nk <- vapply(seq_len(K), function(k) tabulate(cv[z == k], Tn), numeric(Tn))            # Tn x K
    P <- .pg_str_normP(matrix(stats::rgamma(K * Tn, shape = 1 + as.vector(t(Nk))), K, Tn), sid)
    if (sw > burnin) { j <- sw - burnin; Qsum <- Qsum + Q; llv[j] <- ll; alv[j] <- alpha }
    if (!is.null(tick) && sw %% 50L == 0L) tick(sw / (burnin + nrep))
  }
  list(Q = Qsum / nrep, lnL = mean(llv), lnL_sd = stats::sd(llv), alpha = mean(alv), iter = burnin + nrep)
}

# Several K, several runs each. progress(fraction, text) is optional. Returns list(Q = list by K (best run, clusters ordered by size),
# ll = table by K, runs = every run, S, Tn, model, n).
pg_structure <- function(M, Ks, reps = 5, model = "admix", maxit = 300, seed = NULL, progress = NULL, cv_frac = 0.05,
                         method = "em", burnin = 2000, nrep = 3000, alpha = 1, infer_alpha = TRUE) {
  if (!is.null(seed) && !is.na(seed)) set.seed(seed)
  d <- pg_str_data(M); n <- nrow(M)
  Ks <- sort(unique(as.integer(Ks)))
  if (any(Ks < 1) || any(Ks >= n)) stop("K must be at least 1 and smaller than the number of sequences.")
  # cross-entropy (as LEA / ADMIXTURE): a fraction of the genotypes (sequence x site) is hidden, the model is fitted on the rest,
  # and the hidden alleles are predicted; a smaller cross-entropy = a better K.
  mc <- identical(method, "mcmc"); cv <- cv_frac > 0 && !mc
  if (cv) { msk <- matrix(stats::runif(n * d$S) < cv_frac, n, d$S); Gm <- d$G * (!msk)[, d$sid, drop = FALSE]; Hm <- d$G * msk[, d$sid, drop = FALSE] }
  runs <- list(); best <- list(); tot <- length(Ks) * reps; done <- 0L
  for (K in Ks) {
    bk <- NULL; rl <- numeric(reps); ri <- integer(reps); rce <- rep(NA_real_, reps); ral <- rep(NA_real_, reps)
    for (r in seq_len(reps)) {
      if (mc) {
        fit <- .pg_str_gibbs(d$C, ncol(d$G), d$sid, K, model, burnin, nrep, alpha, infer_alpha,
                             tick = function(f) if (!is.null(progress)) progress((done + f) / tot, sprintf("K = %d, run %d of %d (MCMC %d%%)", K, r, reps, round(100 * f))))
        ral[r] <- fit$alpha
      } else fit <- .pg_str_em(d$G, d$sid, K, model, maxit)
      rl[r] <- fit$lnL; ri[r] <- fit$iter
      if (cv && sum(Hm) > 0) { fm <- .pg_str_em(Gm, d$sid, K, model, maxit); rce[r] <- -sum(Hm * log(pmax(fm$Q %*% fm$P, 1e-300))) / sum(Hm) }
      if (is.null(bk) || fit$lnL > bk$lnL) bk <- fit
      done <- done + 1L; if (!is.null(progress)) progress(done / tot, sprintf("K = %d, run %d of %d", K, r, reps))
    }
    o <- order(colSums(bk$Q), decreasing = TRUE)          # biggest cluster first, so the colours stay in a stable order
    Qb <- bk$Q[, o, drop = FALSE]; dimnames(Qb) <- list(rownames(M), paste0("Q", seq_len(K)))
    best[[as.character(K)]] <- Qb
    runs[[as.character(K)]] <- data.frame(K = K, run = seq_len(reps), lnL = rl, iterations = ri, cross_entropy = rce, alpha = ral)
  }
  runs_df <- do.call(rbind, runs)
  ll <- do.call(rbind, lapply(Ks, function(K) {
    rr <- runs[[as.character(K)]]; x <- rr$lnL; ce <- rr$cross_entropy; p <- if (model == "admix") K * (ncol(d$G) - d$S) + n * (K - 1) else K * (ncol(d$G) - d$S) + K - 1
    data.frame(K = K, runs = length(x), mean_lnL = mean(x), sd_lnL = if (length(x) > 1) stats::sd(x) else NA_real_, best_lnL = max(x),
               runs_at_best = sum(x >= max(x) - 1), BIC = -2 * max(x) + p * log(n * d$S),
               cross_entropy = if (all(is.na(ce))) NA_real_ else mean(ce, na.rm = TRUE),
               sd_cross_entropy = if (sum(!is.na(ce)) > 1) stats::sd(ce, na.rm = TRUE) else NA_real_, stringsAsFactors = FALSE)
  }))
  ll$deltaK <- NA_real_
  if (length(Ks) >= 3 && reps >= 2 && all(diff(Ks) == 1)) {
    L <- vapply(Ks, function(K) runs[[as.character(K)]]$lnL, numeric(reps))      # reps x K
    L <- matrix(L, nrow = reps)
    Lp <- L[, -1, drop = FALSE] - L[, -ncol(L), drop = FALSE]                    # L'(K) for K = Ks[2..]
    Lpp <- abs(Lp[, -1, drop = FALSE] - Lp[, -ncol(Lp), drop = FALSE])            # |L''(K)| for K = Ks[2..(n-1)]
    for (j in seq_len(ncol(Lpp))) { s <- stats::sd(L[, j + 1]); if (is.finite(s) && s > 0) ll$deltaK[j + 1] <- mean(Lpp[, j]) / s }
  }
  list(Q = best, ll = ll, runs = runs_df, S = d$S, Tn = ncol(d$G), model = model, n = n, reps = reps, method = if (mc) "mcmc" else "em", burnin = if (mc) burnin else NA, nrep = if (mc) nrep else NA)
}

# mean ancestry of every population (rows = populations in the order of `pops`)
pg_str_bypop <- function(Q, pop, pops = unique(pop)) {
  out <- t(vapply(pops, function(p) colMeans(Q[pop == p, , drop = FALSE]), numeric(ncol(Q))))
  dimnames(out) <- list(pops, colnames(Q)); out
}

pg_str_palette <- function(K, palette = "figure", custom = "") {
  cu <- trimws(strsplit(custom, "[,; ]+")[[1]]); cu <- cu[grepl("^#?[0-9A-Fa-f]{6}$", cu)]
  cu <- ifelse(substr(cu, 1, 1) == "#", cu, paste0("#", cu))
  base <- switch(palette,
    figure = c("#EE1D52", "#19F0F0", "#FF9EE0", "#2E9132", "#B04E0F", "#0A0AF0", "#FFE29A", "#7F7F7F", "#8E44AD", "#F39C12"),
    bright = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFD700", "#A65628", "#F781BF", "#999999", "#00CED1"),
    soft   = c("#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD"),
    grDevices::hcl.colors(10, "Dark 3"))
  col <- c(cu, base)
  if (K > length(col)) col <- c(col, grDevices::hcl.colors(K - length(col), "Set 3"))
  col[seq_len(K)]
}

# Bar plot as in population-structure papers: K = by population (bar width = number of sequences) and/or K = by individual.
#   mode: "both", "pop" or "ind". pop_order: populations from left to right (others follow). Draws on the current device.
pg_structure_plot <- function(Q, pop, colors, mode = "both", pop_order = NULL, title = "", cex = 0.85, sort_ind = FALSE, ind_names = FALSE) {
  K <- ncol(Q); n <- nrow(Q)
  pops <- unique(pop)
  if (!is.null(pop_order) && length(pop_order)) pops <- c(pop_order[pop_order %in% pops], setdiff(pops, pop_order))
  ord <- unlist(lapply(pops, function(p) {
    i <- which(pop == p)
    if (sort_ind && length(i) > 1) i <- i[order(max.col(Q[i, , drop = FALSE], ties.method = "first"), -apply(Q[i, , drop = FALSE], 1, max))]
    i
  }))
  Qo <- Q[ord, , drop = FALSE]; nm <- rownames(Q)[ord]
  cnt <- as.integer(table(factor(pop, levels = pops))); ends <- cumsum(cnt); starts <- ends - cnt
  PQ <- pg_str_bypop(Q, pop, pops)
  bars <- function(Qm, x0, x1) {
    cum <- t(apply(Qm, 1, cumsum)); if (K == 1) cum <- t(cum); lo <- cbind(0, cum[, -K, drop = FALSE])
    for (k in seq_len(K)) graphics::rect(x0, lo[, k], x1, cum[, k], col = colors[k], border = NA)
  }
  lab_in <- max(graphics::strwidth(pops, units = "inches", cex = cex)) * 0.72 + 0.15
  if (isTRUE(ind_names)) lab_in <- max(lab_in, max(graphics::strwidth(nm, units = "inches", cex = cex * 0.6)) * 0.72 + 0.1)
  both <- identical(mode, "both")
  if (both) graphics::layout(matrix(1:2, 2), heights = c(1, 1.15 + lab_in / 3)) else graphics::layout(1)
  op <- graphics::par(no.readonly = TRUE); on.exit({ graphics::par(op[intersect(names(op), c("mai", "xpd", "xaxs", "yaxs"))]); graphics::layout(1) }, add = TRUE)
  panel <- function(kind, bottom) {
    graphics::par(mai = c(if (bottom) lab_in else 0.05, 1.15, if (kind == "pop" || !both) 0.3 else 0.12, 0.15), xaxs = "i", yaxs = "i", xpd = FALSE)
    graphics::plot(NA, xlim = c(0, n), ylim = c(0, 1), axes = FALSE, xlab = "", ylab = "")
    if (kind == "pop") { bars(PQ, starts, ends) } else { bars(Qo, seq_len(n) - 1, seq_len(n)) }
    graphics::rect(starts, 0, ends, 1, border = "black", col = NA, lwd = 1)
    graphics::mtext(sprintf("K=%d (by %s)", K, if (kind == "pop") "pop" else "ind"), side = 2, line = 0.5, las = 1, font = 2, cex = cex)
    if (bottom) {
      graphics::par(xpd = NA)
      graphics::text((starts + ends) / 2, -0.03, pops, srt = 45, adj = c(1, 1), cex = cex)
      if (isTRUE(ind_names) && kind == "ind") graphics::text(seq_len(n) - 0.5, -0.03, nm, srt = 60, adj = c(1, 1), cex = cex * 0.6)
    }
  }
  if (nzchar(title)) { first <- TRUE } else first <- FALSE
  if (both) { panel("pop", FALSE); if (first) graphics::mtext(title, side = 3, adj = 0, line = 0.4, font = 2, cex = cex * 1.15); panel("ind", TRUE) }
  else { panel(mode, TRUE); if (first) graphics::mtext(title, side = 3, adj = 0, line = 0.4, font = 2, cex = cex * 1.15) }
  invisible(list(pops = pops, order = ord))
}


# Pie charts of the mean ancestry of every population on a plain longitude / latitude plot (add coastlines when the 'maps' package is installed).
#   PQ: populations x K; coords: data.frame(Population, lon, lat) ; size: radius as % of the width of the map; by_n: radius grows with sqrt(n)
pg_structure_map <- function(PQ, n, coords, colors, size = 3, by_n = FALSE, labels = TRUE, cex = 0.8, title = "") {
  m <- match(tolower(rownames(PQ)), tolower(coords$Population)); ok <- !is.na(m)
  if (!any(ok)) stop("None of the populations has coordinates.")
  PQ <- PQ[ok, , drop = FALSE]; n <- n[ok]; x <- coords$lon[m[ok]]; y <- coords$lat[m[ok]]
  px <- diff(range(x)); py <- diff(range(y)); pad <- max(px, py, 1) * 0.12
  xl <- range(x) + c(-pad, pad); yl <- range(y) + c(-pad, pad)
  if (requireNamespace("maps", quietly = TRUE)) {
    ok_map <- tryCatch({ maps::map("world", xlim = xl, ylim = yl, fill = TRUE, col = "grey93", border = "grey70", mar = c(4, 4, 2, 1)); graphics::axis(1); graphics::axis(2); graphics::box(); TRUE }, error = function(e) FALSE)
  } else ok_map <- FALSE
  if (!ok_map) { graphics::par(mar = c(4, 4, 2, 1)); graphics::plot(NA, xlim = xl, ylim = yl, asp = 1, xlab = "", ylab = ""); graphics::grid(col = "grey90", lty = 1) }
  graphics::title(main = title, xlab = "Longitude", ylab = "Latitude", adj = 0, font.main = 2)
  r0 <- size / 100 * diff(xl); K <- ncol(PQ)
  for (i in seq_len(nrow(PQ))) {
    r <- if (by_n) r0 * sqrt(n[i] / max(n)) * 1.4 else r0
    a <- -pi / 2; q <- PQ[i, ]
    for (k in seq_len(K)) if (q[k] > 0.001) {
      a1 <- a + 2 * pi * q[k]; th <- seq(a, a1, length.out = max(3, ceiling(60 * q[k])))
      graphics::polygon(x[i] + c(0, r * cos(th)), y[i] + c(0, r * sin(th)), col = colors[k], border = if (q[k] > 0.999) NA else "grey20", lwd = 0.5); a <- a1
    }
    graphics::symbols(x[i], y[i], circles = r, inches = FALSE, add = TRUE, fg = "grey20", lwd = 0.8)
    if (labels) graphics::text(x[i], y[i] + r * 1.15, rownames(PQ)[i], pos = 3, cex = cex, offset = 0.2, xpd = NA)
  }
  invisible(NULL)
}


# ---- Arlequin settings file (.ars) --------------------------------------------------------------------------------------
# Reads the settings that Arlequin saves (Settings > "Save settings") and turns the ones that this program has into the values
# of its own form; the others are listed, so nothing is silently dropped. pg_ars_write() does the opposite: it puts the values of
# the form back into the .ars file that was opened (every other line is kept as it was), so the same settings can be run in Arlequin.
pg_ars_parse <- function(txt) {
  ln <- trimws(unlist(strsplit(txt, "\r?\n"))); ln <- ln[!grepl("^#", ln)]
  m <- regmatches(ln, regexec("^([A-Za-z][A-Za-z0-9_]*)\\s*=\\s*(.*)$", ln)); m <- m[lengths(m) == 3]
  stats::setNames(lapply(m, function(x) trimws(x[3])), vapply(m, function(x) x[2], ""))
}

.pg_ars_dist <- c(`0` = "N", `1` = "JC69", `2` = "K80", `3` = "N", `6` = "TN93")      # Arlequin DistanceMethod code -> distance of this program

pg_ars_to_settings <- function(kv) {
  g <- function(k) { v <- kv[[k]]; if (is.null(v)) NA_real_ else suppressWarnings(as.numeric(v)) }
  task <- g("TaskNumber"); has_task <- is.finite(task)
  bit <- function(b) !has_task || bitwAnd(as.integer(task), as.integer(b)) > 0L
  vals <- list(); app <- list(); notes <- character()
  put <- function(key, id, val, label) { vals[[id]] <<- val; app[[length(app) + 1L]] <<- data.frame(Setting = label, Arlequin = key, Value = if (is.logical(val)) (if (val) "yes" else "no") else as.character(val), stringsAsFactors = FALSE) }
  clamp <- function(x, lo, hi, what) { y <- min(max(x, lo), hi); if (y != x) notes <<- c(notes, sprintf("%s = %s is outside the range of this program; %s was used.", what, format(x), format(y))); y }

  if (is.finite(g("NumPermutationsAMOVA")))     put("NumPermutationsAMOVA", "nperm", as.integer(clamp(g("NumPermutationsAMOVA"), 100, 1e7, "NumPermutationsAMOVA")), "AMOVA: permutations")
  if (is.finite(g("NumPermutationsFST")))       put("NumPermutationsFST", "nperm_pair", as.integer(clamp(g("NumPermutationsFST"), 50, 1e7, "NumPermutationsFST")), "Pairwise FST: permutations")
  if (is.finite(g("FSTSignifLevel")))           put("FSTSignifLevel", "alpha", clamp(g("FSTSignifLevel"), 0.001, 0.5, "FSTSignifLevel"), "Significance level")
  dm <- g("DistanceMethodAMOVA")
  if (is.finite(dm)) {
    k <- as.character(as.integer(dm))
    if (k %in% names(.pg_ars_dist)) {
      put("DistanceMethodAMOVA", "dmodel", unname(.pg_ars_dist[k]), "Distance method (AMOVA, pairwise FST)")
      if (k == "3") notes <- c(notes, "DistanceMethodAMOVA = 3 (p-distance): the number of pairwise differences was used (the same FST, only the scale of the variance components differs).")
    } else notes <- c(notes, sprintf("DistanceMethodAMOVA = %s (Tajima-Nei or Tamura) is not available here; the distance method of the form was kept.", k))
  }
  if (is.finite(g("GammaAValueAMOVA")))         put("GammaAValueAMOVA", "gamma", max(0, g("GammaAValueAMOVA")), "Gamma a value")
  # neutrality tests and the exact test: Arlequin runs them only when their bit of TaskNumber is on
  if (is.finite(g("NumSimulFuTest")) || has_task)
    put("NumSimulFuTest", "nsim", if (bit(16384L)) as.integer(nz_(g("NumSimulFuTest"), 1000)) else 0L, "Fu's Fs: simulated samples")
  ew_on <- bit(128L) && !identical(g("EwensWattersonHomozygosityTest"), 0)
  if (is.finite(g("NumIterationsNeutralityTests")) || has_task)
    put("NumIterationsNeutralityTests", "ew_nsim", if (ew_on) as.integer(nz_(g("NumIterationsNeutralityTests"), 1000)) else 0L, "Ewens-Watterson: simulated samples")
  if (is.finite(g("NumPermutPopDiff")) || has_task)
    put("NumPermutPopDiff", "ex_steps", if (bit(8192L)) as.integer(nz_(g("NumPermutPopDiff"), 100000)) else 0L, "Exact test of differentiation: Markov chain steps")
  if (is.finite(g("NumDememoPopDiff")))         put("NumDememoPopDiff", "ex_demem", as.integer(g("NumDememoPopDiff")), "Exact test: dememorization steps")
  if (is.finite(g("PrintSlatkinsDist")))        put("PrintSlatkinsDist", "slatkin", g("PrintSlatkinsDist") == 1, "Show Slatkin's linearized FST")
  if (is.finite(g("PrintConcestryCoeff")))      put("PrintConcestryCoeff", "reyn", g("PrintConcestryCoeff") == 1, "Show Reynolds' distance")
  if (is.finite(g("PrintMissIntermatchs")))     put("PrintMissIntermatchs", "pipop", g("PrintMissIntermatchs") == 1, "Show pairwise differences (pi) between populations")
  if (is.finite(g("AllowedLevelOfMissingData"))) put("AllowedLevelOfMissingData", "miss", clamp(g("AllowedLevelOfMissingData"), 0, 0.5, "AllowedLevelOfMissingData"), "Allowed level of missing data per site (AMOVA, pairwise FST, pi)")

  if (has_task) {
    off <- c(`512` = "AMOVA", `1024` = "Pairwise FST")
    for (b in names(off)) if (!bit(as.integer(b))) notes <- c(notes, sprintf("%s is switched off in this .ars file; here it is run with the button 'AMOVA and pairwise FST'.", off[[b]]))
    no <- c(`4` = "Mismatch distribution", `8` = "Haplotype frequency estimation (EM)", `16` = "Linkage disequilibrium", `32` = "Hardy-Weinberg test",
            `256` = "Chakraborty's test", `2048` = "Pairwise genetic distances (as a separate table)", `32768` = "Genotype assignment test", `65536` = "Mantel test")
    for (b in names(no)) if (bit(as.integer(b))) notes <- c(notes, sprintf("%s is switched on in this .ars file, but it is not available in this program (run it in Arlequin).", no[[b]]))
  }
  w <- c(DeletionWeight = g("DeletionWeight"), TransitionWeight = g("TransitionWeight"), TranversionWeight = g("TranversionWeight"))
  if (any(is.finite(w) & w != 1)) notes <- c(notes, "Weights of deletions, transitions or transversions different from 1 are not supported: all were taken as 1.")
  if (identical(g("ComputeConventionalFST"), 1)) notes <- c(notes, "Conventional FST (haplotype frequencies only) is not supported: FST uses the distances between haplotypes.")
  if (is.finite(g("AllowedLevelOfMissingData")) && g("AllowedLevelOfMissingData") > 0)
    notes <- c(notes, "AllowedLevelOfMissingData is used for AMOVA, pairwise FST and pairwise pi (sites above the level are dropped; the missing entries in the other sites are skipped pair by pair). The diversity indices, neutrality tests and the exact test still use only sites that are complete in every sequence.")
  # settings of Arlequin that this program does not have: every one that is switched on (or not at its neutral value) is listed
  chk <- list(
    list("DetectSelectedLoci", 0, "Detection of loci under selection (FDIST2, Beaumont & Nichols) is switched on: not available here (it needs many independent loci; with one linked sequence it is not meaningful)."),
    list("ComputeSumStatsWithinGroups", 0, "ComputeSumStatsWithinGroups = 1: summary statistics are computed per population here, not for pooled groups."),
    list("InferHaplotypeDefinitionFromDistanceMatrix", 0, "'Infer haplotype definition from distance matrix' is on: here haplotypes are the distinct sequences."),
    list("EliminateRedondHaplodefs", 0, "EliminateRedondHaplodefs is on: not available (identical haplotypes are always merged here)."),
    list("LocByLocAMOVA", 0, "Locus-by-locus AMOVA is on: not available (a sequence is analysed as one haplotypic locus)."),
    list("IncludeIndividualLevel", 0, "Individual level (FIS) in AMOVA: not available."),
    list("PopSpecificFST", 0, "Population-specific FST: not available."),
    list("PopSpecificFIS", 0, "Population-specific FIS: not available."),
    list("UnequalPopSizeDiv", 0, "Distances for unequal population sizes: not available."),
    list("PrintMinSpannNetworkPop", 0, "Minimum spanning tree/network per population: not available (use the Network tab)."),
    list("PrintMinSpannNetworkGlob", 0, "Minimum spanning network among haplotypes: not available (use the Network tab)."),
    list("PrintHaplDistMat", 0, "Print the distance matrix between haplotypes: not available."),
    list("PrintDistanceMatrix", 0, "Print the distance matrix between haplotypes (AMOVA): not available."),
    list("ComputeAllAlleleFreqs", 0, "Output of sample allele frequencies at all loci: not available."),
    list("ComputeMissDataFreq", 0, "Frequencies of missing data: not available."),
    list("ComputePopSFS", 0, "Site frequency spectrum: not available."), list("ComputeJointSFS", 0, "Joint site frequency spectrum: not available."),
    list("ComputeMultiSFS", 0, "Multidimensional site frequency spectrum: not available."),
    list("KeepNullDistrib", 0, "KeepNullDistrib: the null distributions of AMOVA are not saved to files."),
    list("MakeHWExactTest", 0, "Hardy-Weinberg exact test: not available."), list("MakeExactTestLD", 0, "Linkage disequilibrium test: not available."))
  for (cc in chk) { v <- g(cc[[1]]); if (is.finite(v) && v != cc[[2]]) notes <- c(notes, cc[[3]]) }
  if (is.finite(g("DistanceMethod")) && g("DistanceMethod") != 0) notes <- c(notes, "DistanceMethod (molecular diversity) is not 'pairwise difference': the mean number of pairwise differences is always used for the diversity indices.")
  if (is.finite(g("GammaAValue")) && g("GammaAValue") > 0) notes <- c(notes, "GammaAValue (molecular diversity) > 0: gamma-distributed rates are not used for the diversity indices (only for the AMOVA distance).")
  if (is.finite(g("TestSignificanceAMOVA")) && g("TestSignificanceAMOVA") == 0) notes <- c(notes, "Significance testing of AMOVA is off in this file; here the permutations are always run (set a small number to make it quick).")
  if (is.finite(g("TestSignificancePairewiseFST")) && g("TestSignificancePairewiseFST") == 0) notes <- c(notes, "Significance testing of pairwise FST is off in this file; here the permutations are always run.")
  if (is.finite(g("SignLevelPopDiff")) && is.finite(g("FSTSignifLevel")) && abs(g("SignLevelPopDiff") - g("FSTSignifLevel")) > 1e-9)
    notes <- c(notes, sprintf("The significance level of the exact test (%s) differs from that of pairwise FST (%s); one level (%s) is used for both.", format(g("SignLevelPopDiff")), format(g("FSTSignifLevel")), format(g("FSTSignifLevel"))))
  th <- g("ComputeTheta"); if (is.finite(th) && th != 15) notes <- c(notes, sprintf("ComputeTheta = %d: all four theta estimates are always printed.", as.integer(th)))
  if (has_task && !bit(64L)) notes <- c(notes, "Tajima's D is switched off in this .ars file; it is always computed here.")
  if (has_task && !bit(2L)) notes <- c(notes, "Molecular diversity is switched off in this .ars file; it is always computed here.")
  if (has_task && is.finite(g("NumSimulFuTest")) && is.finite(g("NumIterationsNeutralityTests")))
    notes <- c(notes, "The simulated Tajima's D p-value uses the same number of samples as Fu's Fs (NumSimulFuTest), as in Arlequin.")
  list(values = vals, applied = if (length(app)) do.call(rbind, app) else data.frame(Setting = character(), Arlequin = character(), Value = character()), notes = notes)
}
nz_ <- function(x, d) if (is.finite(x)) x else d

# values: list(nperm, nperm_pair, alpha, dmodel, gamma, nsim, ew_nsim, ex_steps, ex_demem, slatkin) as in the form
pg_ars_write <- function(txt, v) {
  setk <- function(t, key, val) { pat <- sprintf("(?m)^(%s[ \t]*=[ \t]*)[^\r\n]*", key); if (!grepl(pat, t, perl = TRUE)) t else sub(pat, paste0("\\1", val), t, perl = TRUE) }
  inv <- c(N = 0L, JC69 = 1L, K80 = 2L, TN93 = 6L)
  kv <- pg_ars_parse(txt); task <- suppressWarnings(as.integer(kv[["TaskNumber"]]))
  if (!is.na(task)) {
    sw <- function(t, b, on) if (on) bitwOr(t, b) else bitwAnd(t, bitwNot(b))
    task <- sw(task, 16384L, isTRUE(v$nsim > 0)); task <- sw(task, 128L, isTRUE(v$ew_nsim > 0)); task <- sw(task, 8192L, isTRUE(v$ex_steps > 0))
    txt <- setk(txt, "TaskNumber", task)
  }
  txt <- setk(txt, "NumPermutationsAMOVA", as.integer(v$nperm)); txt <- setk(txt, "NumPermutationsFST", as.integer(v$nperm_pair))
  txt <- setk(txt, "FSTSignifLevel", sprintf("%.6f", v$alpha)); txt <- setk(txt, "SignLevelPopDiff", sprintf("%.6f", v$alpha))
  if (!is.na(inv[v$dmodel])) txt <- setk(txt, "DistanceMethodAMOVA", inv[[v$dmodel]])
  txt <- setk(txt, "GammaAValueAMOVA", sprintf("%.6f", v$gamma))
  if (isTRUE(v$nsim > 0)) txt <- setk(txt, "NumSimulFuTest", as.integer(v$nsim))
  if (isTRUE(v$ew_nsim > 0)) { txt <- setk(txt, "NumIterationsNeutralityTests", as.integer(v$ew_nsim)); txt <- setk(txt, "EwensWattersonHomozygosityTest", 1L) }
  else txt <- setk(txt, "EwensWattersonHomozygosityTest", 0L)
  if (isTRUE(v$ex_steps > 0)) txt <- setk(txt, "NumPermutPopDiff", as.integer(v$ex_steps))
  txt <- setk(txt, "NumDememoPopDiff", as.integer(v$ex_demem)); txt <- setk(txt, "PrintSlatkinsDist", as.integer(isTRUE(v$slatkin)))
  txt <- setk(txt, "PrintConcestryCoeff", as.integer(isTRUE(v$reyn))); txt <- setk(txt, "PrintMissIntermatchs", as.integer(isTRUE(v$pipop)))
  if (is.finite(v$miss)) txt <- setk(txt, "AllowedLevelOfMissingData", sprintf("%.6f", v$miss))
  txt
}


# ---- Arlequin-style details of the diversity analysis --------------------------------------------------------------------
# Everything below is what the Arlequin result file prints in addition to the DnaSP-style table: the standard deviations of the mean
# number of pairwise differences (Tajima 1993) and of theta_Pi (Tajima 1983), transitions / transversions, private polymorphic
# sites, nucleotide composition, expected heterozygosity of every polymorphic site, and the simulated p of Tajima's D.

pg_arl_sd <- function(n, k) {
  if (n < 2 || !is.finite(k)) return(c(SD_k = NA_real_, SD_theta_Pi = NA_real_))
  v93 <- (3 * n * (n + 1) * k + 2 * (n^2 + n + 3) * k^2) / (11 * n^2 - 7 * n + 6)
  v83 <- (n + 1) * k / (3 * (n - 1)) + 2 * (n^2 + n + 3) * k^2 / (9 * n * (n - 1))
  c(SD_k = sqrt(max(v93, 0)), SD_theta_Pi = sqrt(max(v83, 0)))
}

# Tajima's D p-value as Arlequin: coalescent samples of n sequences with theta = theta_S of the data (S / a1); the D of every
# simulated sample is compared with the observed one: p = P(D_sim < D_obs). D is not defined for n < 4 or S = 0 (those samples are skipped).
pg_tajima_sim <- function(n, S, k, nsim = 1000, progress = NULL) {
  nm <- c("p", "theta_S", "mean_theta_S", "sd_theta_S", "mean_D", "sd_D", "nsim")
  out <- stats::setNames(rep(NA_real_, 7), nm)
  D <- pg_tajima(n, S, k)[["D"]]
  if (nsim <= 0 || !is.finite(D)) return(out)
  a1 <- sum(1 / seq_len(n - 1)); th <- S / a1
  sm <- matrix(NA_real_, nsim, 2)
  for (r in seq_len(nsim)) {
    s <- .pg_sim_one(n, th)
    sm[r, ] <- c(s[["S"]] / a1, pg_tajima(n, s[["S"]], s[["pi"]])[["D"]])
    if (!is.null(progress) && r %% 100 == 0) progress(r / nsim)
  }
  ok <- is.finite(sm[, 2])
  if (!any(ok)) return(out)
  out[] <- c(mean(sm[ok, 2] < D - 1e-9), th, mean(sm[, 1]), stats::sd(sm[, 1]), mean(sm[ok, 2]), stats::sd(sm[ok, 2]), nsim)
  out
}

.pg_ts <- function(a, b) paste(sort(c(a, b)), collapse = "") %in% c("AG", "CT")          # transition (purine-purine or pyrimidine-pyrimidine)

# details of every unit: "All" (all sequences) and each population. M = sequences x sites (A, C, G, T, - allowed), pop = population of each row.
pg_arl_details <- function(M, pop) {
  units <- c(list(All = seq_len(nrow(M))), split(seq_len(nrow(M)), pop))
  isvar <- function(X) if (ncol(X) == 0) integer(0) else which(apply(X, 2, function(x) length(unique(x)) > 1))
  poly <- lapply(units, function(ii) isvar(M[ii, , drop = FALSE]))
  out <- lapply(names(units), function(u) {
    ii <- units[[u]]; X <- M[ii, , drop = FALSE]; n <- nrow(X); pl <- poly[[u]]
    ts <- tv <- ts_s <- tv_s <- indel <- 0L
    for (j in pl) {
      x <- X[, j]; if (any(x == "-")) indel <- indel + 1L
      tb <- table(x[x %in% c("A", "C", "G", "T")]); if (length(tb) < 2) next
      ref <- names(tb)[which.max(tb)]; oth <- setdiff(names(tb), ref)
      t1 <- sum(vapply(oth, function(o) .pg_ts(ref, o), TRUE)); ts <- ts + t1; tv <- tv + (length(oth) - t1)
      ts_s <- ts_s + (t1 > 0); tv_s <- tv_s + (t1 < length(oth))
    }
    others <- if (u == "All") integer(0) else unique(unlist(poly[setdiff(names(units), c("All", u))]))
    comp <- table(factor(as.vector(X), levels = c("A", "C", "G", "T"))); comp <- 100 * as.numeric(comp) / max(1, sum(comp))
    het <- function(j) { cn <- table(X[, j]); if (n < 2) 0 else n / (n - 1) * (1 - sum((cn / n)^2)) }
    list(unit = u, n = n, L = ncol(X), poly = pl, ts = ts, tv = tv, subst = ts + tv, indel = indel, ts_sites = ts_s, tv_sites = tv_s,
         subst_sites = length(pl), private = if (u == "All") NA_integer_ else length(setdiff(pl, others)),
         comp = stats::setNames(comp, c("A", "C", "G", "T")),
         site_het = vapply(pl, het, 0), site_all = vapply(pl, function(j) length(unique(X[, j])), 1L))
  })
  stats::setNames(out, names(units))
}

.pg_f <- function(x, d = 5) ifelse(is.finite(x), formatC(x, digits = d, format = "f"), "N.A.")

# The block that Arlequin prints for one population sample. st = row of the diversity table (data.frame, one row), ew = row of the
# Ewens-Watterson table, dt = pg_arl_details()[[pop]], tj = pg_tajima_sim(), sites = original alignment column of every site of M.
pg_arl_pop_text <- function(st, ew, dt, tj, nsim_fs, sites = NULL) {
  f <- .pg_f; nm <- st$Population; n <- st$n
  ln <- function(lab, val) sprintf("%-45s: %s", lab, val)
  out <- c(strrep("=", 70), paste("== Sample :", nm), strrep("=", 70), "",
           sprintf("== Standard diversity indices : (%s)", nm), "",
           ln("No. of gene copies", n), ln("No. of haplotypes (sequences)", st$h), ln("No. of usable sites", st$Sites), ln("No. of polymorphic sites", dt$subst_sites), "")
  if (dt$subst_sites > 0) {
    s <- if (is.null(sites)) dt$poly else sites[dt$poly]
    out <- c(out, "Results are only shown for polymorphic sites", strrep("-", 49), "           Num.", "           gene     Num.        Exp.", "Site#     copies   alleles       Het", strrep("-", 49),
             sprintf("%6s %9d %9d %10s", s, n, dt$site_all, f(dt$site_het)), strrep("-", 49),
             sprintf("  Mean %9.3f %9.3f %10s", n, mean(dt$site_all), f(mean(dt$site_het))), sprintf("  s.d. %9.3f %9.3f %10s", 0, stats::sd(dt$site_all), f(stats::sd(dt$site_het))),
             strrep("-", 49), "")
  }
  out <- c(out, "Haplotype-level computations", "----------------------------",
           ln("Gene diversity", sprintf("%s +/- %s", f(st$Hd, 4), f(st$SD_Hd, 4))), "",
           sprintf("== Molecular diversity indices : (%s)", nm), "",
           ln("Sample size", f(n, 4)), ln("No. of haplotypes", st$h), "",
           ln("Number of observed transitions", dt$ts), ln("Number of observed transversions", dt$tv), ln("Number of substitutions", dt$subst),
           ln("Number of observed indels", dt$indel), ln("Number of polymorphic sites", dt$subst_sites),
           ln("Number of observed sites with transitions", dt$ts_sites), ln("Number of observed sites with transversions", dt$tv_sites),
           ln("Number of observed sites with substitutions", dt$subst_sites), ln("Number of private polymorphic sites", if (is.na(dt$private)) "N.A." else dt$private),
           ln("Number of usable nucleotide sites", dt$L), "",
           "Nucleotide composition (relative values)", sprintf("   %s : %5.2f%%", names(dt$comp), dt$comp), "",
           ln("Mean number of pairwise differences", sprintf("%s +/- %s", f(st$k, 6), f(st$SD_k, 6))),
           ln("Nucleotide diversity (average over sites)", sprintf("%s +/- %s", f(st$Pi, 6), f(st$SD_Pi, 6))),
           "(Standard deviations: Tajima 1993 for the mean number of pairwise differences, Tajima 1983 for theta_Pi)", "",
           ln("Theta(Hom)", f(st$theta_Hom, 6)), ln("Theta(k)", f(st$theta_k, 6)), ln("Theta(S)", f(st$theta_S, 6)),
           ln("Theta(Pi)", f(st$theta_Pi, 6)), ln("S.D. Theta(Pi)", f(st$SD_theta_Pi, 6)), "",
           sprintf("== Tajima's test of selective neutrality : (%s)", nm), "",
           ln("Sample size", n), ln("No. of sites with substitutions (S)", st$S), ln("Mean No. of pairwise differences (Pi)", f(st$k, 5)), "",
           if (is.finite(st$Tajima_D)) c(ln("Tajima's D", f(st$Tajima_D, 5)), ln("P-value, beta distribution (two-sided; Tajima 1989)", f(st$p_D, 5))) else "Tajima's D is not defined (it needs at least 4 sequences and one polymorphic site)",
           if (is.finite(tj[["p"]])) c(ln("No. of simulations", tj[["nsim"]]), ln("Obs. Theta(S)", f(tj[["theta_S"]], 5)), ln("Mean Theta(S)", f(tj[["mean_theta_S"]], 5)),
                                       ln("S.D. Theta(S)", f(tj[["sd_theta_S"]], 5)), ln("Mean D", f(tj[["mean_D"]], 5)), ln("S.D. D", f(tj[["sd_D"]], 5)),
                                       ln("P(D simul < D obs) (as Arlequin)", f(tj[["p"]], 5))), "",
           sprintf("== Ewens-Watterson tests of selective neutrality : (%s)", nm), "")
  if (is.finite(ew$F_mean_simulated)) out <- c(out, ln("No. of genes in sample", ew$n), ln("No. of haplotypes in sample", ew$h), ln("Observed F value", f(ew$F_obs)),
                                               ln("Expected F value", f(ew$F_mean_simulated)), ln("No. of simulated samples", ew$simulated),
                                               ln("Pr(rand F <= obs F)", f(ew$p_F_sim_le_F_obs)), ln("Pr(rand F >= obs F)", f(ew$p_F_sim_ge_F_obs)))
  else out <- c(out, ln("No. of genes in sample", ew$n), ln("No. of alleles in sample", ew$h),
                if (isTRUE(ew$h == ew$n)) "The test is impossible because all gene copies are different" else if (isTRUE(ew$h < 2)) "The test is impossible with a single haplotype" else "Not computed (the number of simulated samples is 0)")
  out <- c(out, "", sprintf("== Fu's Fs test of selective neutrality : (%s)", nm), "")
  if (is.finite(st$Fs)) out <- c(out, ln("Original No. of alleles (k)", st$h), ln("Theta(Pi)", f(st$k, 5)), ln("Fs", f(st$Fs, 5)),
                                 if (is.finite(st$p_Fs)) c(ln("No. of simulations", nsim_fs), ln("Prob(sim_Fs <= obs_Fs)", f(st$p_Fs, 5))))
  else out <- c(out, "Fu's Fs cannot be computed (one haplotype, or no difference between the sequences)")
  c(out, "")
}

# The tables of "Summary of computations done within populations". units = names of the populations (without "All").
pg_arl_summary <- function(stats, ew, det, units, sites = NULL) {
  f <- .pg_f; stu <- stats[match(units, stats$Population), , drop = FALSE]; ewu <- ew[match(units, ew$Population), , drop = FALSE]
  ms <- function(M) { M <- as.matrix(M); cbind(M, Mean = apply(M, 1, function(x) mean(as.numeric(x), na.rm = TRUE)), s.d. = apply(M, 1, function(x) stats::sd(as.numeric(x), na.rm = TRUE))) }
  dd <- det[units]
  basic <- cbind(Statistics = c("No. of gene copies", "No. of usable sites", "No. of polymorphic sites"),
                 as.data.frame(round(ms(rbind(vapply(dd, function(x) x$n, 0), vapply(dd, function(x) x$L, 0), vapply(dd, function(x) x$subst_sites, 0))), 3), row.names = NULL))
  names(basic)[2:(1 + length(units))] <- units
  allsites <- sort(unique(unlist(lapply(dd, function(x) x$poly))))
  mat <- function(g) { M <- matrix(0, length(allsites), length(units)); for (u in seq_along(units)) { x <- dd[[u]]; M[match(x$poly, allsites), u] <- g(x) }; M }
  lab <- if (is.null(sites)) allsites else sites[allsites]
  tabs <- function(M, d) { colnames(M) <- units; Mm <- ms(M); data.frame(Site = c(lab, "Mean", "s.d."), round(rbind(Mm[seq_along(allsites), , drop = FALSE], colMeans(Mm[seq_along(allsites), , drop = FALSE]),
                                                       apply(Mm[seq_along(allsites), , drop = FALSE], 2, stats::sd)), d), check.names = FALSE, row.names = NULL) }
  het <- if (length(allsites)) tabs(mat(function(x) x$site_het), 5) else NULL
  nal <- if (length(allsites)) tabs({ M <- mat(function(x) x$site_all); M[M == 0] <- 1; M }, 3) else NULL
  mol <- data.frame(Statistics = c("No. of transitions", "No. of transversions", "No. of substitutions", "No. of indels", "No. of private polymorphic sites",
                                   "Mean number of pairwise differences (k)", "S.D. of k", "Nucleotide diversity per site (Pi)", "Theta_Hom", "Theta_k", "Theta_S", "Theta_Pi", "S.D. Theta_Pi"),
                    stats::setNames(lapply(seq_along(units), function(i) { x <- dd[[i]]; s <- stu[i, ]
                      c(x$ts, x$tv, x$subst, x$indel, x$private, round(s$k, 5), round(s$SD_k, 5), round(s$Pi, 6), round(s$theta_Hom, 5), round(s$theta_k, 5), round(s$theta_S, 5), round(s$theta_Pi, 5), round(s$SD_theta_Pi, 5)) }), units),
                    check.names = FALSE)
  neu <- data.frame(Statistics = c("Ewens-Watterson: sample size", "Ewens-Watterson: no. of alleles", "Observed F", "Expected F", "Pr(rand F <= obs F)",
                                   "Tajima's D", "Tajima's D p (beta)", "Tajima's D p (simulated, as Arlequin)", "Fu's Fs", "Fu's Fs p"),
                    stats::setNames(lapply(seq_along(units), function(i) { e <- ewu[i, ]; s <- stu[i, ]
                      c(e$n, e$h, round(e$F_obs, 5), round(e$F_mean_simulated, 5), round(e$p_F_sim_le_F_obs, 5), round(s$Tajima_D, 5), round(s$p_D, 5), round(s$p_D_sim, 5), round(s$Fs, 5), round(s$p_Fs, 5)) }), units),
                    check.names = FALSE)
  list(basic = basic, het = het, nal = nal, mol = mol, neu = neu)
}

# AMOVA written as Arlequin prints it. am = pg_amova() result; the names of the variance components follow Arlequin (Va, Vb, Vc).
pg_arl_amova_text <- function(am, n_perm = am$nperm) {
  o <- am$obs; f <- .pg_f
  row <- function(lab, df, ssd, v, pct, tag) sprintf(" %-30s %5d %12s %14s %-3s %8s", lab, as.integer(df), f(ssd, 3), f(v, 5), tag, if (is.na(pct)) "" else f(pct, 2))
  hdr <- c(strrep("-", 80), sprintf(" %-30s %5s %12s %14s %12s", "Source of variation", "d.f.", "Sum of squares", "Variance comp.", "% variation"), strrep("-", 80))
  if (am$levels == 2) {
    sg <- unname(o$sigma[c("a", "w")]); tot <- sum(sg)
    body <- c(row("Among populations", o$df[["a"]], o$ssd[["a"]], sg[1], 100 * sg[1] / tot, "Va"), row("Within populations", o$df[["w"]], o$ssd[["w"]], sg[2], 100 * sg[2] / tot, "Vb"),
              row("Total", o$df[["t"]], o$ssd[["t"]], tot, NA, ""), strrep("-", 80), sprintf(" Fixation index   FST : %s", f(o$phiST, 5)), strrep("-", 80), "",
              sprintf("Significance tests (%d permutations)", n_perm), "Va and FST : P(rand. value >= obs. value) = %s", sprintf("             P-value = %s", f(am$p[["phiST"]], 5)))
    body[length(body) - 1] <- sprintf("Va and FST : P(rand. value >= obs. value) = %s", f(am$p[["phiST"]], 5)); body <- body[-length(body)]
  } else {
    sg <- unname(o$sigma[c("a", "b", "c")]); tot <- sum(sg)
    body <- c(row("Among groups", o$df[["ag"]], o$ssd[["ag"]], sg[1], 100 * sg[1] / tot, "Va"), row("Among populations within groups", o$df[["ap"]], o$ssd[["ap"]], sg[2], 100 * sg[2] / tot, "Vb"),
              row("Within populations", o$df[["wp"]], o$ssd[["wp"]], sg[3], 100 * sg[3] / tot, "Vc"), row("Total", o$df[["t"]], o$ssd[["t"]], tot, NA, ""), strrep("-", 80),
              "Fixation indices", sprintf("   FSC : %s", f(o$phiSC, 5)), sprintf("   FST : %s", f(o$phiST, 5)), sprintf("   FCT : %s", f(o$phiCT, 5)), strrep("-", 80), "",
              sprintf("Significance tests (%d permutations)", n_perm), sprintf("Vc and FST  : P(rand. value >= obs. value) = %s", f(am$p[["phiST"]], 5)),
              sprintf("Vb and FSC  : P(rand. value >= obs. value) = %s", f(am$p[["phiSC"]], 5)), sprintf("Va and FCT  : P(rand. value >= obs. value) = %s", f(am$p[["phiCT"]], 5)))
  }
  c(hdr, body)
}


# ---- the result file in the layout of the Arlequin result browser ------------------------------------------------------------
# pg_arl_report() returns the lines of the result (already escaped for HTML, with <A NAME> anchors) and the entries of the tree;
# pg_arl_files() writes them as the three files that Arlequin writes (NAME_main.htm = frames, NAME_tree.htm = menu, NAME.htm = results);
# pg_arl_single() puts the menu and the results in ONE html file (nothing else is needed to open it).
.pg_tj_na <- c(p = NA_real_, theta_S = NA_real_, mean_theta_S = NA_real_, sd_theta_S = NA_real_, mean_D = NA_real_, sd_D = NA_real_, nsim = NA_real_)

.pg_txt_table <- function(d, first_w = 40, w = 10) {          # fixed-width text of a data.frame (first column = row names)
  d <- as.data.frame(d, check.names = FALSE, stringsAsFactors = FALSE); if (!ncol(d)) return(character())
  cell <- function(x) { x <- as.character(x); x[is.na(x) | x == "NA" | x == "NaN"] <- "N.A."; x }
  hw <- pmax(w, nchar(names(d)[-1]) + 1L)
  hd <- paste0(formatC(names(d)[1], width = -first_w), paste0(vapply(seq_along(hw), function(j) formatC(names(d)[j + 1L], width = hw[j]), ""), collapse = " "))
  rows <- vapply(seq_len(nrow(d)), function(i) paste0(formatC(cell(d[[1]])[i], width = -first_w), paste0(vapply(seq_along(hw), function(j) formatC(cell(d[[j + 1]])[i], width = hw[j]), ""), collapse = " ")), "")
  ln <- strrep("-", nchar(hd)); c(ln, hd, ln, rows, ln)
}
.pg_lower_txt <- function(m, f = .pg_f, diag = "0.00000", star = FALSE) {      # lower triangle, columns numbered as in Arlequin
  P <- nrow(m); out <- paste0(strrep(" ", 12), paste(formatC(seq_len(P), width = 9), collapse = " "))
  for (i in seq_len(P)) out <- c(out, paste0(formatC(i, width = 12), "  ", paste(vapply(seq_len(i), function(j) formatC(if (i == j) (if (star) "*" else diag) else f(m[i, j]), width = 9), ""), collapse = " ")))
  out
}
.pg_pi_txt <- function(pp) {            # Arlequin layout: PiX on the diagonal, PiXY above, corrected PiXY - (PiX + PiY) / 2 below
  P <- nrow(pp$xy); m <- matrix("", P, P)
  for (i in seq_len(P)) for (j in seq_len(P)) m[i, j] <- .pg_f(if (i == j) pp$xy[i, i] else if (i < j) pp$xy[i, j] else pp$net[i, j], 5)
  out <- paste0(strrep(" ", 12), paste(formatC(seq_len(P), width = 11), collapse = " "))
  for (i in seq_len(P)) out <- c(out, paste0(formatC(i, width = 12), paste(formatC(m[i, ], width = 11), collapse = " ")))
  out
}
.pg_sig_txt <- function(m, alpha) {
  P <- nrow(m); out <- paste0(strrep(" ", 12), paste(formatC(seq_len(P), width = 9), collapse = " "))
  for (i in seq_len(P)) out <- c(out, paste0(formatC(i, width = 12), "  ", paste(vapply(seq_len(P), function(j) formatC(if (i == j) "" else if (is.finite(m[i, j]) && m[i, j] < alpha) "+" else "-", width = 9), ""), collapse = " ")))
  out
}

pg_arl_report <- function(d, a = NULL, info = list()) {
  esc <- .pg_esc; f <- .pg_f
  lines <- character(); tree <- data.frame(id = character(), label = character(), level = integer(), stringsAsFactors = FALSE)
  anchor <- function(id, label, level = 0L) { lines <<- c(lines, sprintf("<A NAME=\"%s\"></A>", id)); tree[nrow(tree) + 1L, ] <<- list(id, label, level) }
  add <- function(x) lines <<- c(lines, esc(x))
  title <- if (is.null(info$title)) "SangeRly" else info$title
  st <- d$stats; for (cc in c("SD_k", "SD_Pi", "SD_theta_Pi", "p_D_sim")) if (is.null(st[[cc]])) st[[cc]] <- NA_real_
  det <- if (!is.null(d$det)) d$det else pg_arl_details(d$M, d$pop)
  sites <- d$cols; units <- setdiff(st$Population, "All")
  anchor("run_information", "Run information", 0L)
  add(c(strrep("/", 66), sprintf("%s (%s)", title, format(Sys.time(), "%d/%m/%y at %H:%M:%S")), strrep("/", 66), "",
        "SangeRly results in the layout of the Arlequin result file. These are SangeRly's own computations, not Arlequin output.", "",
        "Project information:", "--------------------", sprintf("\tNbSamples      =  %d", length(units)), "\tDataType       =  DNA", "\tGenotypicData  =  0",
        sprintf("\tNo. of sequences = %d; usable sites (no gap / N / IUPAC code) = %d of %d", nrow(d$M), ncol(d$M), d$total), "",
        "Notes for comparing with Arlequin:",
        " - Nucleotide diversity is per usable site: if Arlequin was given only the variable sites its value is larger (same mean number of pairwise differences).",
        " - Theta(Hom) here is (1 - F) / F; Arlequin applies the iterative correction of Chakraborty & Weiss (1991) and gives a larger value.",
        " - Probabilities from simulations (Tajima's D, Fu's Fs, Ewens-Watterson) and the Markov chain of the exact test change a little from run to run; use more samples / steps for a final value.", ""))
  if (length(info$settings)) {
    anchor("settings", "Settings used for calculations", 1L)
    add(c(strrep("=", 30), "Settings used for Calculations", strrep("=", 30), "", info$settings, ""))
  }
  anchor("intra", "Analyses at the intra-population level", 0L)
  add(c(strrep("=", 78), "== ANALYSES AT THE INTRA-POPULATION LEVEL", strrep("=", 78), ""))
  for (u in units) {
    i <- match(u, st$Population); j <- match(u, d$ew$Population)
    anchor(paste0("pop_", make.names(u)), paste("Sample:", u), 1L)
    tj <- if (!is.null(d$tj[[u]])) d$tj[[u]] else .pg_tj_na
    lines <- c(lines, esc(pg_arl_pop_text(st[i, , drop = FALSE], d$ew[j, , drop = FALSE], det[[u]], tj, info$nsim_fs, sites)))
  }
  if (length(units) > 1) {
    sm <- pg_arl_summary(st, d$ew, det, units, sites)
    anchor("summary", "Summary of computations done within populations", 0L)
    add(c(strrep("=", 80), "== Summary of computations done within populations", strrep("=", 80), ""))
    for (k in c("basic", "het", "nal", "mol", "neu")) if (!is.null(sm[[k]])) {
      lab <- c(basic = "Basic properties", het = "Expected heterozygosity of every polymorphic site", nal = "Number of alleles of every polymorphic site", mol = "Molecular diversity indices", neu = "Neutrality tests")[[k]]
      anchor(paste0("sum_", k), lab, 1L); add(c(lab, .pg_txt_table(sm[[k]], if (k %in% c("het", "nal")) 12 else 42, 12), ""))
    }
  }
  if (!is.null(a)) {
    anchor("structure", "Genetic structure analysis", 0L)
    add(c(strrep("=", 80), "== GENETIC STRUCTURE ANALYSIS", strrep("=", 80), "", sprintf("Number of usable sites for distance computation : %d", a$sites), ""))
    anchor("amova", "AMOVA", 1L)
    grp <- if (is.null(info$groups)) NULL else info$groups
    add(c(strrep("=", 65), "  AMOVA ANALYSIS", strrep("=", 65), "", sprintf("Distance method: %s%s", c(N = "Pairwise difference", K80 = "Kimura 2-parameter", TN93 = "Tamura-Nei", JC69 = "Jukes-Cantor")[[a$model]],
                                                                         if (a$gamma > 0 && a$model != "N") sprintf(" (gamma a = %g)", a$gamma) else ""), "",
          if (length(grp)) c("Genetic structure to test:", unlist(lapply(names(grp), function(g) c(sprintf("  Group %s: %s", g, paste(grp[[g]], collapse = ", ")))))) , "",
          "AMOVA design and results:", "  Weir, B.S. and Cockerham, C.C. 1984; Excoffier, L., Smouse, P., and Quattro, J. 1992.", pg_arl_amova_text(a$am), ""))
    pw <- a$pw; ph <- pw$phi; P <- nrow(ph)
    anchor("pairwise", "Comparisons of pairs of population samples", 1L)
    add(c(strrep("=", 80), "== Comparisons of pairs of population samples", strrep("=", 80), "", "List of labels for population samples used below:", strrep("-", 49), "Label  \tPopulation name", "-----  \t---------------",
          sprintf("%3d:\t%s", seq_len(P), rownames(ph)), "", strrep("-", 24), "Population pairwise FSTs", strrep("-", 24), ""))
    add(c(sprintf("Distance method: %s", c(N = "Pairwise difference", K80 = "Kimura 2-parameter", TN93 = "Tamura-Nei", JC69 = "Jukes-Cantor")[[a$model]]), .pg_lower_txt(ph), ""))
    if (isTRUE(a$slatkin)) add(c("Slatkin's linearized FST = FST / (1 - FST)", .pg_lower_txt(ph / (1 - ph)), ""))
    if (isTRUE(a$reyn)) add(c(strrep("-", 20), "Reynolds' distances", strrep("-", 20), "", "D = -ln(1 - FST)", .pg_lower_txt(pg_reynolds(ph)), ""))
    if (!is.null(a$pipop)) {
      pp <- a$pipop; pn <- rownames(pp$xy)
      add(c(strrep("-", 60), "Pairwise differences between populations", strrep("-", 60), "",
            "Distance method: same as above. Diagonal: mean number of pairwise differences within populations (PiX);",
            "above the diagonal: average pairwise differences between populations (PiXY); below the diagonal: corrected PiXY - (PiX + PiY) / 2", "",
            .pg_pi_txt(pp), ""))
    }
    add(c(strrep("-", 12), "FST P values", strrep("-", 12), "", sprintf("Number of permutations : %d", info$nperm_pair), "", .pg_lower_txt(pw$p, function(x) formatC(x, digits = 5, format = "f"), star = TRUE), "",
          strrep("-", 12), "Matrix of significant Fst P values", sprintf("Significance Level=%.4f", a$alpha), strrep("-", 12), "", .pg_sig_txt(pw$p, a$alpha), ""))
  }
  ex <- d$exact
  if (!is.null(ex) && !is.na(ex$p)) {
    pops <- colnames(d$counts); P <- length(pops)
    anchor("exact", "Exact test of sample differentiation", 0L)
    add(c(strrep("=", 80), "== Exact Test of Sample Differentiation Based on Haplotype Frequencies", strrep("=", 80), "  Raymond M. and F. Rousset. 1995; Goudet, J., M. Raymond, T. de Meeus and F. Rousset, 1996.", "",
          "List of labels for population samples used below:", strrep("-", 49), "Label  \tPopulation name", "-----  \t---------------", sprintf("%3d\t%s", seq_len(P), pops), "",
          strrep("-", 74), "Global test of differentiation among samples:", strrep("-", 44), "",
          sprintf(" Non-differentiation: Exact P value = %s +- %s  (%d Markov steps done)", f(ex$p), f(ex$se), ex$steps), ""))
    if (!is.null(ex$pairwise)) {
      pv <- ex$pairwise
      add(c(strrep("-", 74), "Differentiation test between all pairs of samples:", strrep("-", 50), "", "Non-differentiation exact P values :", "", .pg_lower_txt(pv, star = TRUE), "",
            strrep("-", 32), sprintf("Table of significant differences (significance level=%.4f):", if (is.null(info$alpha)) 0.05 else info$alpha), strrep("-", 32), "", .pg_sig_txt(pv, if (is.null(info$alpha)) 0.05 else info$alpha), ""))
      sg <- pv < (if (is.null(info$alpha)) 0.05 else info$alpha); diag(sg) <- FALSE; sg[is.na(sg)] <- FALSE
      add(c(strrep("-", 68), "|" %+% paste(formatC(seq_len(P), width = 10), collapse = "|") %+% "|", strrep("-", 68), "|" %+% paste(formatC(colSums(sg), width = 10), collapse = "|") %+% "|", strrep("-", 68), ""))
    }
  }
  if (length(info$dnasp)) { anchor("geneflow", "Genetic differentiation and gene flow (DnaSP layout)", 0L); add(c(strrep("=", 80), "== GENETIC DIFFERENTIATION AND GENE FLOW", strrep("=", 80), "", info$dnasp, "")) }
  add(c(strrep("/", 66), sprintf("END OF RUN (%s)", format(Sys.time(), "%d/%m/%y at %H:%M:%S")), strrep("/", 66)))
  list(lines = lines, tree = tree)
}
`%+%` <- function(a, b) paste0(a, b)

.pg_arl_page_css <- "BODY {background-color: white; font-size: 8pt; font-family: Courier New}"
pg_arl_files <- function(name, rep) {
  main <- sprintf("<html>\n<!--SangeRly results frame file (layout of the Arlequin result browser)-->\n<head>\n<title> SangeRly Result Browser - %s\n</title>\n</head>\n<FRAMESET cols=\"230,*\">\n  <FRAME src=\"%s_tree.htm\" name=\"treeframe\">\n  <FRAME SRC=\"%s.htm\" name=\"basefrm\">\n</FRAMESET>\n</html>\n", .pg_esc(name), name, name)
  tr <- rep$tree
  tree <- paste0("<html><head><meta charset=\"utf-8\"><style>BODY{background:white;font-family:Arial,sans-serif;font-size:9pt}a{text-decoration:none;color:#0b3d91}a:hover{text-decoration:underline}</style></head><body>\n<b>", .pg_esc(name), "</b><br>\n",
                 paste0("<div style=\"margin-left:", 12 * tr$level, "px\"><a href=\"", name, ".htm#", tr$id, "\" target=\"basefrm\">", .pg_esc(tr$label), "</a></div>", collapse = "\n"), "\n</body></html>\n")
  res <- paste0("<html>\n<!--SangeRly results file-->\n<head>\n<meta charset=\"utf-8\">\n<style>\n", .pg_arl_page_css, "\n</style>\n</head>\n<body bgcolor=white>\n<PRE>\n", paste(rep$lines, collapse = "\n"), "\n</PRE>\n</body>\n</html>\n")
  stats::setNames(list(main, tree, res), paste0(name, c("_main.htm", "_tree.htm", ".htm")))
}
pg_arl_single <- function(name, rep) {
  tr <- rep$tree
  paste0("<!DOCTYPE html><html><head><meta charset=\"utf-8\"><title>", .pg_esc(name), "</title><style>", .pg_arl_page_css,
         "\n#menu{position:fixed;left:0;top:0;bottom:0;width:230px;overflow:auto;background:#f3f5f8;border-right:1px solid #cfd6df;padding:8px;font-family:Arial,sans-serif;font-size:9pt}",
         "\n#menu a{text-decoration:none;color:#0b3d91;display:block}#menu a:hover{text-decoration:underline}#res{margin-left:244px}</style></head><body>\n<div id=\"menu\"><b>", .pg_esc(name), "</b><br>\n",
         paste0("<a style=\"margin-left:", 12 * tr$level, "px\" href=\"#", tr$id, "\">", .pg_esc(tr$label), "</a>", collapse = "\n"), "\n</div>\n<div id=\"res\"><PRE>\n", paste(rep$lines, collapse = "\n"), "\n</PRE></div></body></html>\n")
}
# writes <dir>/<name>.res/ with the three files; returns the path of the _main.htm file
pg_arl_save <- function(dir, name, rep) {
  name <- gsub("[^A-Za-z0-9_.-]", "_", name); if (!nzchar(name)) name <- "sangerly"
  out <- file.path(dir, paste0(name, ".res")); dir.create(out, showWarnings = FALSE, recursive = TRUE)
  if (!dir.exists(out)) stop("The folder could not be created: ", out)
  fl <- pg_arl_files(name, rep)
  for (n in names(fl)) writeLines(fl[[n]], file.path(out, n), useBytes = TRUE)
  normalizePath(file.path(out, paste0(name, "_main.htm")), mustWork = FALSE)
}

.pg_default_dir <- function() {
  h <- Sys.getenv(if (.Platform$OS.type == "windows") "USERPROFILE" else "HOME"); d <- file.path(h, "Documents")
  gsub("\\\\", "/", if (nzchar(h) && dir.exists(d)) d else if (nzchar(h)) h else getwd())
}

# ---------------------------------------------------------------------------------------------------
# Fu & Li's D* and F* (1993, without an outgroup), Ramos-Onsins & Rozas' R2 (2002)
# Only sites with two alleles are used (as DnaSP 6). Significance: neutral coalescent samples with the same n and the same
# number of segregating sites S (the mutations are placed on the branches in proportion to their length).
# ---------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------
# Achaz (2009) general neutrality tests on the site frequency spectrum: Tajima's D, Fu and Li's D* and F*, Achaz's Y*
# Statistic = (Theta1 - Theta2) / sqrt(alpha * theta + beta * theta^2); Theta = sum_i w_i xi_i with xi_i = number of sites with i copies of the
# derived (here: any) variant; variances and covariances of xi_i from Fu (1995). theta = S / a_n, theta^2 = S (S - 1) / (a_n^2 + b_n).
# ---------------------------------------------------------------------------------------------------
.pg_fu95 <- function(n) {                      # matrix sigma2 (theta^2 part of Cov(xi_i, xi_j)); the theta part is diag(1 / i)
  i <- seq_len(n - 1); ai <- c(0, cumsum(1 / seq_len(n - 1)))[i]            # a_i = sum_{k < i} 1 / k
  an <- sum(1 / seq_len(n - 1)); an1 <- an + 1 / n
  bet <- function(k) 2 * n / ((n - k + 1) * (n - k)) * (an1 - c(0, cumsum(1 / seq_len(n)))[k]) - 2 / (n - k)     # beta_n(k), k = 1..n
  S2 <- matrix(0, n - 1, n - 1)
  for (a in i) for (b in i) {
    if (a == b) S2[a, b] <- if (2 * a < n) bet(a + 1) else if (2 * a == n) 2 * (an - ai[a]) / (n - a) - 1 / a^2 else bet(a) - 1 / a^2
    else { u <- max(a, b); v <- min(a, b)
      S2[a, b] <- if (u + v < n) (bet(u + 1) - bet(u)) / 2
                  else if (u + v == n) (an - ai[u]) / (n - u) + (an - ai[v]) / (n - v) - (bet(u) + bet(v + 1)) / 2 - 1 / (u * v)
                  else (bet(v) - bet(v + 1)) / 2 - 1 / (u * v) }
  }
  S2
}
.pg_w <- function(n, what) {                    # weights of the estimators of theta (symmetric in i and n - i, unbiased under the neutral model)
  i <- seq_len(n - 1); an <- sum(1 / i)
  switch(what,
    pi = i * (n - i) / choose(n, 2),
    W = rep(1 / an, n - 1),
    s = { w <- numeric(n - 1); w[c(1, n - 1)] <- (n - 1) / n; w },                              # singletons (folded: both ends)
    Wns = { w <- numeric(n - 1); k <- i >= 2 & i <= n - 2; w[k] <- 1 / sum(1 / i[k]); w },     # without singletons
    pins = { w <- numeric(n - 1); k <- i >= 2 & i <= n - 2; w[k] <- i[k] * (n - i[k]) / sum((n - i[k])); w })
}
.pg_achaz_const <- function(n) {
  an <- sum(1 / seq_len(n - 1)); bn <- sum(1 / seq_len(n - 1)^2); S2 <- .pg_fu95(n); i <- seq_len(n - 1)
  mk <- function(a, b) { d <- .pg_w(n, a) - .pg_w(n, b); c(alpha = sum(d^2 / i), beta = as.numeric(t(d) %*% S2 %*% d), w1 = list(.pg_w(n, a)), w2 = list(.pg_w(n, b))) }
  tests <- list(D = c("pi", "W"), Dstar = c("W", "s"), Fstar = c("pi", "s"), Ystar = c("Wns", "pins"))
  ct <- lapply(tests, function(p) { d <- .pg_w(n, p[1]) - .pg_w(n, p[2]); list(w1 = .pg_w(n, p[1]), w2 = .pg_w(n, p[2]), alpha = sum(d^2 / i), beta = as.numeric(t(d) %*% S2 %*% d)) })
  list(n = n, an = an, bn = bn, ct = ct, S2 = S2)
}
# xi: vector (length n - 1) with the number of sites having i copies of one of the two variants (any variant: weights are symmetric)
.pg_achaz_xi <- function(xi, cst) {
  S <- sum(xi); th <- S / cst$an; th2 <- S * (S - 1) / (cst$an^2 + cst$bn)
  vapply(cst$ct, function(z) { v <- z$alpha * th + z$beta * th2; if (!is.finite(v) || v <= 0) NA_real_ else (sum(z$w1 * xi) - sum(z$w2 * xi)) / sqrt(v) }, numeric(1))
}
.pg_xi <- function(H) { n <- nrow(H); cs <- colSums(H); tabulate(cs, n - 1) }

.pg_fuli_const <- function(n) {
  an <- sum(1 / seq_len(n - 1)); bn <- sum(1 / seq_len(n - 1)^2); an1 <- an + 1 / n
  cn <- if (n == 2) 1 else 2 * (n * an - 2 * (n - 1)) / ((n - 1) * (n - 2))
  dn <- cn + (n - 2) / (n - 1)^2 + 2 / (n - 1) * (1.5 - (2 * an1 - 3) / (n - 2) - 1 / n)
  vD <- ((n / (n - 1))^2 * bn + an^2 * dn - 2 * n * an * (dn + 1) / ((n - 1) * (n - 2))) / (an^2 + bn)
  uD <- n / (n - 1) * (an - n / (n - 1)) - vD
  vF <- ((2 * n^3 + 110 * n^2 - 255 * n + 153) / (9 * n^2 * (n - 1)) + 2 * (n - 1) * an / n^2 - 8 * bn / n) / (an^2 + bn)
  uF <- ((4 * n^2 + 19 * n + 3 - 12 * (n + 1) * an1) / (3 * n * (n - 1))) / an - vF
  list(an = an, vD = vD, uD = uD, vF = vF, uF = uF, ach = .pg_achaz_const(n))
}

# H: n x S matrix of 0/1 (one column per biallelic site). Returns D*, F*, R2 (and the numerators used for the checks).
.pg_fuli_H <- function(H, cst = NULL) {
  n <- nrow(H); S <- ncol(H)
  if (n < 4 || S < 1) return(c(Dstar = NA_real_, Fstar = NA_real_, Ystar = NA_real_, R2 = NA_real_, eta_s = NA_real_, k = NA_real_, Dstar_achaz = NA_real_, Fstar_classic = NA_real_))
  if (is.null(cst)) cst <- .pg_fuli_const(n)
  cs <- colSums(H); one <- cs == 1 | cs == n - 1
  eta <- sum(one); k <- sum(2 * cs * (n - cs)) / (n * (n - 1))
  Dn <- n / (n - 1) * S - cst$an * eta; Dd <- sqrt(cst$uD * S + cst$vD * S^2)
  Fn <- k - (n - 1) / n * eta;          Fd <- sqrt(cst$uF * S + cst$vF * S^2)
  U <- numeric(n)
  if (eta > 0) for (j in which(one)) { r <- if (cs[j] == 1) which(H[, j] == 1) else which(H[, j] == 0); U[r] <- U[r] + 1 }
  ac <- .pg_achaz_xi(tabulate(cs, n - 1), cst$ach)           # DnaSP 6: D* as Fu and Li (1993), F* with the general variances of Achaz (2009), Y* of Achaz (2008)
  c(Dstar = if (Dd > 0) Dn / Dd else NA_real_, Fstar = ac[["Fstar"]], Ystar = ac[["Ystar"]],
    R2 = sqrt(sum((U - k / 2)^2) / n) / S, eta_s = eta, k = k,
    Dstar_achaz = ac[["Dstar"]], Fstar_classic = if (Fd > 0) Fn / Fd else NA_real_)
}

# coalescent sample of n sequences with exactly S (infinite-sites) mutations; returns the n x S matrix of 0/1
.pg_sim_H <- function(n, S, theta = NULL) {
  lineages <- as.list(seq_len(n)); birth <- rep(0, n)
  len <- numeric(0); members <- list(); tnow <- 0; m <- n
  while (m > 1) {
    tnow <- tnow + stats::rexp(1, m * (m - 1) / 2)
    pick <- sample.int(m, 2)
    for (q in pick) { len <- c(len, tnow - birth[q]); members[[length(members) + 1L]] <- lineages[[q]] }
    newl <- c(lineages[[pick[1]]], lineages[[pick[2]]])
    keep <- setdiff(seq_len(m), pick)
    lineages <- c(lineages[keep], list(newl)); birth <- c(birth[keep], tnow)
    m <- m - 1L
  }
  if (!is.null(theta)) S <- stats::rpois(1, theta * sum(len) / 2)      # theta given: S is random (used to check the formulas)
  cnt <- as.vector(stats::rmultinom(1, S, len / sum(len)))
  H <- matrix(0L, n, S); col <- 0L
  for (b in which(cnt > 0)) { cc <- col + seq_len(cnt[b]); H[members[[b]], cc] <- 1L; col <- col + cnt[b] }
  H
}

.pg_fuli_na <- c(n = NA, S_biallelic = NA, singletons = NA, D_star = NA, p_D_star_le = NA, p_D_star_ge = NA,
                 F_star = NA, p_F_star_le = NA, p_F_star_ge = NA, Y_star = NA, p_Y_star_le = NA, p_Y_star_ge = NA, R2 = NA, p_R2_le = NA, p_R2_ge = NA,
                 D_star_Achaz_variance = NA, F_star_classic_variance = NA, nsim = NA)

# M: sequences x sites (A, C, G, T). nsim = number of simulated samples (0 = statistics only).
pg_fuli <- function(M, nsim = 1000, progress = NULL) {
  n <- nrow(M); out <- .pg_fuli_na; out[["n"]] <- n
  if (ncol(M) > 0) M <- M[, colSums(!(M == "A" | M == "C" | M == "G" | M == "T")) == 0, drop = FALSE]   # only sites complete in this set
  if (n < 4 || ncol(M) == 0) return(out)
  nal <- apply(M, 2, function(x) length(unique(x)))
  bi <- which(nal == 2); S <- length(bi)
  out[["S_biallelic"]] <- S
  if (S == 0) return(out)
  H <- vapply(bi, function(j) as.integer(M[, j] != M[1, j]), integer(n))
  if (S == 1) H <- matrix(H, n, 1)
  ob <- .pg_fuli_H(H)
  out[c("singletons", "D_star", "F_star", "Y_star", "R2", "D_star_Achaz_variance", "F_star_classic_variance")] <- ob[c("eta_s", "Dstar", "Fstar", "Ystar", "R2", "Dstar_achaz", "Fstar_classic")]
  if (nsim > 0) {
    cst <- .pg_fuli_const(n); sm <- matrix(NA_real_, nsim, 4)
    for (r in seq_len(nsim)) {
      sm[r, ] <- .pg_fuli_H(.pg_sim_H(n, S), cst)[c("Dstar", "Fstar", "R2", "Ystar")]
      if (!is.null(progress) && r %% 100 == 0) progress(r / nsim)
    }
    tail <- function(x, o, lo) { ok <- is.finite(x); if (!is.finite(o) || !any(ok)) NA_real_ else if (lo) mean(x[ok] <= o + 1e-9) else mean(x[ok] >= o - 1e-9) }
    out[c("p_D_star_le", "p_D_star_ge")] <- c(tail(sm[, 1], ob[["Dstar"]], TRUE), tail(sm[, 1], ob[["Dstar"]], FALSE))
    out[c("p_F_star_le", "p_F_star_ge")] <- c(tail(sm[, 2], ob[["Fstar"]], TRUE), tail(sm[, 2], ob[["Fstar"]], FALSE))
    out[c("p_Y_star_le", "p_Y_star_ge")] <- c(tail(sm[, 4], ob[["Ystar"]], TRUE), tail(sm[, 4], ob[["Ystar"]], FALSE))
    out[c("p_R2_le", "p_R2_ge")] <- c(tail(sm[, 3], ob[["R2"]], TRUE), tail(sm[, 3], ob[["R2"]], FALSE))
    out[["nsim"]] <- nsim
  }
  out
}

# ---------------------------------------------------------------------------------------------------
# sliding window. M = sites kept (complete deletion), cols = their column numbers in the alignment, total = alignment length,
# pop = population of every row; units = "All" and every population. Windows are laid on the ALIGNMENT columns; the per-site values
# use the sites kept inside the window.
# ---------------------------------------------------------------------------------------------------
pg_sliding <- function(M, cols, total, pop, window = 100, step = 25, min_sites = 10, by_pop = TRUE) {
  window <- max(2L, as.integer(window)); step <- max(1L, as.integer(step))
  if (total < window) window <- total
  starts <- seq(1L, max(1L, total - window + 1L), by = step)
  units <- c(list(All = seq_len(nrow(M))), if (by_pop) split(seq_len(nrow(M)), pop))
  res <- list()
  for (u in names(units)) {
    X <- M[units[[u]], , drop = FALSE]; n <- nrow(X)
    if (n < 2) next
    cnt <- rbind(colSums(X == "A"), colSums(X == "C"), colSums(X == "G"), colSums(X == "T"))
    dif <- n^2 - colSums(cnt^2); var <- colSums(cnt > 0) > 1
    a1 <- sum(1 / seq_len(n - 1))
    for (s in starts) {
      e <- s + window - 1L; ii <- which(cols >= s & cols <= e); L <- length(ii)
      if (L < min_sites) { res[[length(res) + 1L]] <- data.frame(Population = u, Start = s, End = e, Mid = (s + e) / 2, Sites = L, S = NA_real_, Pi = NA_real_, ThetaW = NA_real_, TajimaD = NA_real_, h = NA_real_, Hd = NA_real_); next }
      S <- sum(var[ii]); k <- sum(dif[ii]) / (n * (n - 1))
      hp <- pg_haplotypes(X[, ii, drop = FALSE]); p <- hp$n / n
      res[[length(res) + 1L]] <- data.frame(Population = u, Start = s, End = e, Mid = (s + e) / 2, Sites = L, S = S, Pi = k / L, ThetaW = S / a1 / L,
        TajimaD = if (n >= 4) pg_tajima(n, S, k)[["D"]] else NA_real_, h = length(hp$n), Hd = n / (n - 1) * (1 - sum(p^2)), stringsAsFactors = FALSE)
    }
  }
  if (!length(res)) return(NULL)
  do.call(rbind, res)
}

# ---------------------------------------------------------------------------------------------------
# median-joining network (Bandelt, Forster & Roehl 1999), as the program NETWORK (Fluxus) and DnaSP > Network:
#   1. minimum spanning network (union of all minimum spanning trees; links that exceed the tree by at most `eps` are kept too);
#   2. for every vertex and every two of its neighbours, the median vector (majority state of every site) is added when it lowers the connection
#      cost (distance to the three vertices < the cost of joining the three without it; <= with epsilon > 0). A site with three different states
#      has no median, the triplet is skipped;
#   3. steps 1-2 are repeated with the new vertices until no median vector is added;
#   4. median vectors that are left with two links or fewer are removed (the superfluous ones), then the network is drawn again from the remaining vertices.
# Checked against the .out file of NETWORK 10 (Durio zibethinus, epsilon 0): the same 15 links with the same lengths.
# hap_seq: character vector of haplotype sequences (equal length). Returns list(seq = sequences of the haplotypes then the median vectors,
# med = logical, edges = data.frame(a, b) indices into seq).
# ---------------------------------------------------------------------------------------------------
.pg_msn <- function(D, eps = 0) {
  k <- nrow(D); if (k < 2) return(matrix(0L, 0, 2))
  inT <- c(TRUE, rep(FALSE, k - 1)); best <- D[1, ]; parent <- rep(1L, k); tree <- matrix(0L, 0, 2)
  for (s in seq_len(k - 1)) { cand <- which(!inT); j <- cand[which.min(best[cand])]; tree <- rbind(tree, c(parent[j], j)); inT[j] <- TRUE
    upd <- which(!inT & D[j, ] < best); best[upd] <- D[j, upd]; parent[upd] <- j }
  adj <- vector("list", k); for (e in seq_len(nrow(tree))) { a <- tree[e, 1]; b <- tree[e, 2]; adj[[a]] <- c(adj[[a]], b); adj[[b]] <- c(adj[[b]], a) }
  B <- matrix(0L, k, k)                                   # largest edge on the tree path between two vertices
  for (r in seq_len(k)) { seen <- rep(FALSE, k); seen[r] <- TRUE; st <- r
    while (length(st)) { u <- st[length(st)]; st <- st[-length(st)]; for (v in adj[[u]]) if (!seen[v]) { seen[v] <- TRUE; B[r, v] <- max(B[r, u], D[u, v]); st <- c(st, v) } } }
  ii <- which(upper.tri(D) & D <= B + eps, arr.ind = TRUE)
  cbind(ii[, 1], ii[, 2])
}
pg_net_mj <- function(hap_seq, eps = 0, max_vertices = 600L) {
  k0 <- length(hap_seq); S0 <- do.call(rbind, strsplit(hap_seq, ""))
  keep <- which(apply(S0, 2, function(x) length(unique(x)) > 1))
  if (k0 < 2 || !length(keep)) return(list(seq = hap_seq, med = rep(FALSE, k0), edges = data.frame(a = integer(0), b = integer(0))))
  S <- S0[, keep, drop = FALSE]; isMed <- rep(FALSE, k0)
  key <- function(S) apply(S, 1, paste, collapse = "")
  for (iter in 1:60) {
    D <- pg_hamming(S); n <- nrow(S); E <- .pg_msn(D, eps)
    adjm <- matrix(FALSE, n, n); if (nrow(E)) { adjm[E] <- TRUE; adjm[E[, 2:1, drop = FALSE]] <- TRUE }
    ks <- key(S); newS <- list()
    for (u in seq_len(n)) {
      nb <- which(adjm[u, ]); if (length(nb) < 2) next
      pr <- utils::combn(nb, 2)
      for (r in seq_len(ncol(pr))) {
        b <- pr[1, r]; cc <- pr[2, r]; x <- S[u, ]; y <- S[b, ]; z <- S[cc, ]
        m <- ifelse(x == y | x == z, x, ifelse(y == z, y, NA_character_)); if (anyNA(m)) next
        kk <- paste(m, collapse = ""); if (kk %in% ks || !is.null(newS[[kk]])) next
        cost <- sum(m != x) + sum(m != y) + sum(m != z); mst <- D[u, b] + D[u, cc] + D[b, cc] - max(D[u, b], D[u, cc], D[b, cc])
        if (cost < mst || (eps > 0 && cost <= mst + eps)) newS[[kk]] <- m
      }
    }
    if (!length(newS) || n + length(newS) > max_vertices) break
    S <- rbind(S, do.call(rbind, unname(newS))); isMed <- c(isMed, rep(TRUE, length(newS)))
  }
  repeat {                                                # superfluous median vectors
    D <- pg_hamming(S); E <- .pg_msn(D, eps); deg <- tabulate(c(E), nrow(S))
    bad <- which(isMed & deg <= 2); if (!length(bad)) break
    S <- S[-bad[1], , drop = FALSE]; isMed <- isMed[-bad[1]]
  }
  E <- .pg_msn(pg_hamming(S), eps)
  full <- matrix(S0[1, ], nrow(S), ncol(S0), byrow = TRUE); full[, keep] <- S
  list(seq = apply(full, 1, paste, collapse = ""), med = isMed, edges = data.frame(a = E[, 1], b = E[, 2]))
}

# ---------------------------------------------------------------------------------------------------
# map of the haplotypes (as PopART): one pie per population at its coordinates, slices = haplotypes, area = number of sequences
#   counts: haplotypes x populations (integers); coords: data.frame(Population, lon, lat); colors: one per haplotype
# ---------------------------------------------------------------------------------------------------
.pg_hapcols <- function(n) {
  b <- c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffd92f", "#a65628", "#f781bf", "#00ced1", "#999999", "#1b9e77", "#7570b3", "#66a61e", "#e6ab02", "#a6761d", "#666666")
  if (n <= length(b)) b[seq_len(n)] else grDevices::hcl(h = seq(15, 375, length.out = n + 1)[-(n + 1)], c = 80, l = 60)
}
pg_hap_map <- function(counts, coords, colors, size = 4, labels = TRUE, cex = 0.8, title = "", land = "#e8e4d8", sea = "#dbe9f4",
                       hap_legend = TRUE, size_legend = TRUE, pad = 0.15, hires = TRUE) {
  pops <- colnames(counts); m <- match(tolower(pops), tolower(coords$Population)); ok <- !is.na(m)
  if (!any(ok)) stop("None of the populations has coordinates (names must be the same as in section 2).")
  counts <- counts[, ok, drop = FALSE]; pops <- pops[ok]; x <- coords$lon[m[ok]]; y <- coords$lat[m[ok]]
  keep <- rowSums(counts) > 0; counts <- counts[keep, , drop = FALSE]; colors <- colors[keep]
  n <- colSums(counts); nmax <- max(n)
  px <- diff(range(x)); py <- diff(range(y)); span <- max(px, py, 2); padd <- span * pad
  xl <- range(x) + c(-padd, padd) * 1.3; yl <- range(y) + c(-padd, padd)
  # the plot is kept in proportion (1 degree of longitude = 1 degree of latitude)
  if (diff(xl) < diff(yl) * 1.2) { c0 <- mean(xl); xl <- c0 + c(-1, 1) * diff(yl) * 0.6 }
  r0 <- size / 100 * diff(xl); rad <- function(k) r0 * pmax(0.28, sqrt(k / nmax))
  op <- graphics::par(mar = c(4, 4, 2.2, if (hap_legend || size_legend) 9 else 1), xpd = FALSE); on.exit(graphics::par(op), add = TRUE)
  graphics::plot(NA, xlim = xl, ylim = yl, asp = 1, xlab = "Longitude", ylab = "Latitude"); graphics::title(main = title, adj = 0, font.main = 2)
  graphics::rect(graphics::par("usr")[1], graphics::par("usr")[3], graphics::par("usr")[2], graphics::par("usr")[4], col = sea, border = NA)
  if (requireNamespace("maps", quietly = TRUE)) {
    db <- if (hires && diff(xl) < 60 && requireNamespace("mapdata", quietly = TRUE)) "worldHires" else "world"
    tryCatch(maps::map(db, xlim = xl, ylim = yl, fill = TRUE, col = land, border = "grey55", add = TRUE, lwd = 0.5, resolution = 0), error = function(e) tryCatch(maps::map("world", xlim = xl, ylim = yl, fill = TRUE, col = land, border = "grey55", add = TRUE), error = function(e2) NULL))
  }
  graphics::box(); graphics::grid(col = "#ffffff55", lty = 1)
  for (i in order(-n)) {
    r <- rad(n[i]); q <- counts[, i]; q <- q / sum(q); a <- pi / 2
    for (k in which(q > 0)) {
      a1 <- a - 2 * pi * q[k]; th <- seq(a, a1, length.out = max(3, ceiling(70 * q[k])))
      graphics::polygon(x[i] + c(0, r * cos(th)), y[i] + c(0, r * sin(th)), col = colors[k], border = if (q[k] > 0.999) NA else "grey15", lwd = 0.6); a <- a1
    }
    graphics::symbols(x[i], y[i], circles = r, inches = FALSE, add = TRUE, fg = "grey15", lwd = 0.9)
  }
  if (labels) for (i in seq_along(pops)) graphics::text(x[i], y[i] + rad(n[i]), pops[i], pos = 3, cex = cex, offset = 0.25, font = 2, col = "#8b1a1a", xpd = NA)
  usr <- graphics::par("usr"); lx <- usr[2] + diff(usr[1:2]) * 0.02; yy <- usr[4]
  if (size_legend && nmax > 1) {
    nb <- if (nmax >= 100) 100 else if (nmax >= 10) 10 else nmax; rb <- rad(nb); rs <- rad(1)
    cx <- lx + rb; cy <- yy - rb
    graphics::symbols(cx, cy, circles = rb, inches = FALSE, add = TRUE, fg = "grey15", bg = "white", xpd = NA)
    graphics::symbols(cx, cy - rb + rs, circles = rs, inches = FALSE, add = TRUE, fg = "grey15", bg = "white", xpd = NA)
    graphics::segments(cx, cy + rb, cx + rb * 1.3, cy + rb, xpd = NA); graphics::segments(cx, cy - rb + 2 * rs, cx + rb * 1.3, cy - rb + 2 * rs, xpd = NA)
    graphics::text(cx + rb * 1.4, cy + rb, sprintf("%d samples", nb), adj = c(0, 0.5), cex = 0.75, xpd = NA)
    graphics::text(cx + rb * 1.4, cy - rb + 2 * rs, "1 sample", adj = c(0, 0.5), cex = 0.75, xpd = NA)
    yy <- cy - rb - diff(usr[3:4]) * 0.04
  }
  if (hap_legend) graphics::legend(lx, yy, rownames(counts), pch = 21, pt.bg = colors, pt.cex = 1.3, bty = "n", cex = 0.78, xpd = NA, y.intersp = 1.05, ncol = if (nrow(counts) > 22) 2 else 1)
  invisible(NULL)
}

# ---------------------------------------------------------------------------------------------------
# Ewens (1972) sampling formula: Strobeck's S (DnaSP) and Chakraborty's test of population amalgamation (Arlequin)
# P(K = k | n, theta) = |s(n, k)| theta^k / (theta (theta + 1) ... (theta + n - 1))
# ---------------------------------------------------------------------------------------------------
pg_ewens_K <- function(n, theta) {
  if (!is.finite(theta) || theta <= 0 || n < 2) return(rep(NA_real_, max(n, 1)))
  lp <- pg_fs_prep(n) + (1:n) * log(theta) - sum(log(theta + 0:(n - 1))); exp(lp)
}
# n genes, k haplotypes; theta_pi = mean pairwise differences, theta_hom = Theta(Hom) (see pg_thetas)
pg_ewens_tests <- function(n, k, theta_pi, theta_hom) {
  out <- c(n = n, k = k, theta_Pi = theta_pi, S_P_K_le_k = NA, S_P_K_eq_k = NA, theta_Hom = theta_hom, Chakraborty_P_K_ge_k = NA)
  if (n >= 2 && is.finite(theta_pi) && theta_pi > 0 && k >= 1) { p <- pg_ewens_K(n, theta_pi); out[["S_P_K_le_k"]] <- sum(p[seq_len(k)]); out[["S_P_K_eq_k"]] <- p[k] }
  if (n >= 2 && is.finite(theta_hom) && theta_hom > 0 && k >= 1) { p <- pg_ewens_K(n, theta_hom); out[["Chakraborty_P_K_ge_k"]] <- sum(p[k:n]) }
  out
}

# ---------------------------------------------------------------------------------------------------
# folded site frequency spectrum of the biallelic sites (sites without gap / N / IUPAC), with the neutral expectation (constant size)
# ---------------------------------------------------------------------------------------------------
pg_sfs <- function(M) {
  n <- nrow(M); if (ncol(M) > 0) M <- M[, colSums(!(M == "A" | M == "C" | M == "G" | M == "T")) == 0, drop = FALSE]
  if (n < 4 || ncol(M) == 0) return(NULL)
  nal <- apply(M, 2, function(x) length(unique(x))); bi <- which(nal == 2)
  cs <- vapply(bi, function(j) sum(M[, j] != M[1, j]), 1L); cs <- pmin(cs, n - cs)
  cl <- seq_len(floor(n / 2)); S <- length(bi); an <- sum(1 / seq_len(n - 1))
  ex <- (S / an) * (1 / cl + 1 / (n - cl)) / ifelse(2 * cl == n, 2, 1)
  data.frame(Copies_of_minor_variant = cl, Observed_sites = tabulate(cs, length(cl)), Expected_neutral = ex,
             Observed_proportion = tabulate(cs, length(cl)) / max(S, 1), Expected_proportion = ex / max(sum(ex), 1e-12))
}

# ---------------------------------------------------------------------------------------------------
# linkage disequilibrium between pairs of biallelic sites: D, D', r^2, Fisher's exact test; ZnS (Kelly 1997), Za, ZZ (Rozas et al. 2001),
# Rm = minimum number of recombination events (Hudson and Kaplan 1985, four-gamete test)
# ---------------------------------------------------------------------------------------------------
pg_ld <- function(M, cols = NULL, maxsites = 400, fisher = TRUE) {
  n <- nrow(M); ok <- if (ncol(M)) colSums(!(M == "A" | M == "C" | M == "G" | M == "T")) == 0 else logical(0)
  if (is.null(cols)) cols <- seq_len(ncol(M))
  M <- M[, ok, drop = FALSE]; cols <- cols[ok]
  if (n < 4 || ncol(M) < 2) return(NULL)
  nal <- apply(M, 2, function(x) length(unique(x))); bi <- which(nal == 2); S <- length(bi)
  if (S < 2) return(NULL)
  trunc <- FALSE; if (S > maxsites) { bi <- bi[round(seq(1, S, length.out = maxsites))]; trunc <- TRUE }
  H <- vapply(bi, function(j) as.integer(M[, j] != M[1, j]), integer(n)); pos <- cols[bi]; s <- length(bi)
  rows <- vector("list", s * (s - 1) / 2); q <- 0L; inc <- matrix(FALSE, s, s)
  for (a in 1:(s - 1)) for (b in (a + 1):s) {
    n11 <- sum(H[, a] == 1 & H[, b] == 1); n10 <- sum(H[, a] == 1 & H[, b] == 0); n01 <- sum(H[, a] == 0 & H[, b] == 1); n00 <- n - n11 - n10 - n01
    pa <- (n11 + n10) / n; pb <- (n11 + n01) / n; Dv <- n11 / n - pa * pb
    Dmax <- if (Dv >= 0) min(pa * (1 - pb), (1 - pa) * pb) else min(pa * pb, (1 - pa) * (1 - pb))
    r2 <- Dv^2 / (pa * (1 - pa) * pb * (1 - pb))
    inc[a, b] <- n11 > 0 && n10 > 0 && n01 > 0 && n00 > 0
    pf <- if (fisher) stats::fisher.test(matrix(c(n11, n10, n01, n00), 2))$p.value else NA_real_
    q <- q + 1L; rows[[q]] <- c(a, b, pos[a], pos[b], pos[b] - pos[a], Dv, if (Dmax > 0) abs(Dv) / Dmax else NA_real_, r2, pf, inc[a, b])
  }
  tb <- as.data.frame(do.call(rbind, rows)); names(tb) <- c("ia", "ib", "Column_1", "Column_2", "Distance", "D", "Dprime", "r2", "P_Fisher", "Four_gametes")
  tb$Four_gametes <- tb$Four_gametes == 1; tb$P_Bonferroni <- pmin(1, tb$P_Fisher * nrow(tb))
  # Rm: the largest set of disjoint intervals (i, j) between incompatible sites (greedy on the right end)
  iv <- which(inc, arr.ind = TRUE); Rm <- 0L
  if (nrow(iv)) { iv <- iv[order(iv[, 2], iv[, 1]), , drop = FALSE]; last <- -Inf; for (r in seq_len(nrow(iv))) if (iv[r, 1] >= last) { Rm <- Rm + 1L; last <- iv[r, 2] } }
  adj <- tb$ib == tb$ia + 1L; ZnS <- mean(tb$r2, na.rm = TRUE); Za <- mean(tb$r2[adj], na.rm = TRUE)
  sm <- data.frame(Sites = s, Pairs = nrow(tb), ZnS = ZnS, Za = Za, ZZ = Za - ZnS, Rm = Rm, Pairs_with_four_gametes = sum(tb$Four_gametes),
                   Pairs_P_lt_0.05 = sum(tb$P_Fisher < 0.05, na.rm = TRUE), Pairs_P_Bonferroni_lt_0.05 = sum(tb$P_Bonferroni < 0.05, na.rm = TRUE),
                   Sites_in_alignment = S, Subsampled = trunc, check.names = FALSE)
  list(pairs = tb[, setdiff(names(tb), c("ia", "ib"))], summary = sm)
}

# ---------------------------------------------------------------------------------------------------
# Mantel test: correlation between two distance matrices (lower triangles), permuting the rows and columns of the first one
# ---------------------------------------------------------------------------------------------------
pg_mantel <- function(A, B, nperm = 9999, seed = NULL) {
  A <- as.matrix(A); B <- as.matrix(B); m <- nrow(A); if (m < 4) return(NULL)
  lt <- lower.tri(A); ok <- is.finite(A[lt]) & is.finite(B[lt]); if (sum(ok) < 3) return(NULL)
  if (!is.null(seed)) set.seed(seed)
  r0 <- stats::cor(A[lt][ok], B[lt][ok]); rp <- numeric(nperm)
  for (i in seq_len(nperm)) { p <- sample.int(m); Ap <- A[p, p]; x <- is.finite(Ap[lt]) & is.finite(B[lt]); rp[i] <- stats::cor(Ap[lt][x], B[lt][x]) }
  list(r = r0, p_right = (sum(rp >= r0 - 1e-12) + 1) / (nperm + 1), p_left = (sum(rp <= r0 + 1e-12) + 1) / (nperm + 1),
       p_two = (sum(abs(rp) >= abs(r0) - 1e-12) + 1) / (nperm + 1), nperm = nperm, npop = m, npairs = sum(ok))
}
pg_haversine <- function(lon, lat) {
  R <- 6371.0088; r <- pi / 180; m <- length(lon); D <- matrix(0, m, m)
  for (i in seq_len(m)) { dl <- (lon - lon[i]) * r; a <- sin((lat - lat[i]) * r / 2)^2 + cos(lat[i] * r) * cos(lat * r) * sin(dl / 2)^2; D[i, ] <- 2 * R * asin(pmin(1, sqrt(a))) }
  D
}

# ---------------------------------------------------------------------------------------------------
# Hardy-Weinberg test, site by site, for diploid sequences with IUPAC codes at the heterozygous sites (exact test of Wigginton, Cutler and Abecasis 2005)
# ---------------------------------------------------------------------------------------------------
.pg_hwe_exact <- function(nAB, nAA, nBB) {
  n <- nAB + nAA + nBB; nA <- 2 * nAA + nAB; nB <- 2 * nBB + nAB
  rare <- min(nA, nB); mid <- as.integer(rare * (2L * n - rare) / (2L * n)); if (mid %% 2 != rare %% 2) mid <- mid + 1L
  pr <- numeric(rare + 1); pr[mid + 1] <- 1
  h <- mid; hr <- (rare - mid) / 2; hc <- n - mid - hr
  while (h >= 2) { pr[h - 1] <- pr[h + 1] * h * (h - 1) / (4 * (hr + 1) * (hc + 1)); hr <- hr + 1; hc <- hc + 1; h <- h - 2 }
  h <- mid; hr <- (rare - mid) / 2; hc <- n - mid - hr
  while (h <= rare - 2) { pr[h + 3] <- pr[h + 1] * 4 * hr * hc / ((h + 2) * (h + 1)); hr <- hr - 1; hc <- hc - 1; h <- h + 2 }
  pr <- pr / sum(pr); min(1, sum(pr[pr <= pr[nAB + 1] + 1e-12]))
}
.pg_iupac2 <- c(R = "AG", Y = "CT", S = "CG", W = "AT", K = "GT", M = "AC")
# M: sequences x alignment columns (upper case, IUPAC kept). pop: optional vector (one test per population and for all)
pg_hw_sites <- function(M, pop = NULL, cols = NULL, minn = 4) {
  if (is.null(cols)) cols <- seq_len(ncol(M))
  units <- c(list(All = seq_len(nrow(M))), if (!is.null(pop)) split(seq_len(nrow(M)), pop))
  out <- list()
  for (u in names(units)) {
    X <- M[units[[u]], , drop = FALSE]
    for (j in seq_len(ncol(X))) {
      x <- X[, j]; het <- x %in% names(.pg_iupac2); hom <- x %in% c("A", "C", "G", "T")
      if (!any(het)) next
      hc <- unique(x[het]); if (length(hc) != 1) next                       # all heterozygotes must be of the same kind (the two bases of the code)
      b <- strsplit(.pg_iupac2[[hc]], "")[[1]]; if (any(!(x[hom] %in% b))) next
      nAB <- sum(het); nAA <- sum(x == b[1]); nBB <- sum(x == b[2]); n <- nAB + nAA + nBB
      if (n < minn || nAA + nBB == 0 && FALSE) next
      pA <- (2 * nAA + nAB) / (2 * n); He <- 2 * pA * (1 - pA) * 2 * n / (2 * n - 1); Ho <- nAB / n
      out[[length(out) + 1]] <- data.frame(Population = u, Column = cols[j], Alleles = paste(b, collapse = "/"), AA = nAA, AB = nAB, BB = nBB, n = n,
                                           Freq_first = pA, Ho = Ho, He = He, Fis = if (He > 0) 1 - Ho / He else NA_real_, P_exact = .pg_hwe_exact(nAB, nAA, nBB), stringsAsFactors = FALSE)
    }
  }
  if (!length(out)) return(NULL)
  d <- do.call(rbind, out); d$P_Bonferroni <- pmin(1, d$P_exact * ave(d$P_exact, d$Population, FUN = length)); d
}

# ---------------------------------------------------------------------------------------------------
# haplotype phasing by EM (Excoffier and Slatkin 1995) of diploid sequences whose heterozygous sites are IUPAC codes (R, Y, S, W, K, M)
# ---------------------------------------------------------------------------------------------------
pg_phase_em <- function(M, ids = NULL, maxhet = 12, maxit = 1000, tol = 1e-9) {
  n <- nrow(M); if (is.null(ids)) ids <- paste0("Ind", seq_len(n))
  good <- apply(M, 2, function(x) all(x %in% c("A", "C", "G", "T", names(.pg_iupac2))))
  X <- M[, good, drop = FALSE]; kept <- which(good)
  if (ncol(X) == 0) return(NULL)
  nhet <- rowSums(matrix(X %in% names(.pg_iupac2), nrow(X)))
  use <- which(nhet <= maxhet); skipped <- setdiff(seq_len(n), use)
  if (!length(use)) return(NULL)
  cand <- vector("list", length(use))                                         # per individual: matrix of candidate haplotype strings (2 columns)
  for (u in seq_along(use)) {
    x <- X[use[u], ]; hj <- which(x %in% names(.pg_iupac2)); a1 <- a2 <- x
    if (length(hj) == 0) { cand[[u]] <- cbind(paste(x, collapse = ""), paste(x, collapse = "")); next }
    bs <- lapply(hj, function(j) strsplit(.pg_iupac2[[x[j]]], "")[[1]])
    comb <- as.matrix(expand.grid(rep(list(1:2), length(hj)))); comb <- comb[comb[, 1] == 1, , drop = FALSE]   # unordered pairs: first het site fixed
    cand[[u]] <- t(vapply(seq_len(nrow(comb)), function(r) {
      h1 <- x; h2 <- x; for (q in seq_along(hj)) { h1[hj[q]] <- bs[[q]][comb[r, q]]; h2[hj[q]] <- bs[[q]][3 - comb[r, q]] }
      c(paste(h1, collapse = ""), paste(h2, collapse = "")) }, character(2)))
    if (length(hj) == 1) cand[[u]] <- matrix(cand[[u]], 1)
  }
  hap <- sort(unique(unlist(cand))); hi <- lapply(cand, function(m) cbind(match(m[, 1], hap), match(m[, 2], hap)))
  p <- rep(1 / length(hap), length(hap)); N <- length(use); it <- 0; ll <- -Inf
  repeat {
    it <- it + 1; cnt <- numeric(length(hap)); ll_new <- 0
    for (u in seq_len(N)) { m <- hi[[u]]; w <- p[m[, 1]] * p[m[, 2]] * ifelse(m[, 1] == m[, 2], 1, 2); sw <- sum(w); ll_new <- ll_new + log(sw); w <- w / sw
      for (cc in 1:2) { sw2 <- rowsum(w, m[, cc]); ix <- as.integer(rownames(sw2)); cnt[ix] <- cnt[ix] + sw2[, 1] } }
    pn <- cnt / (2 * N); d <- max(abs(pn - p)); p <- pn; ll <- ll_new
    if (d < tol || it >= maxit) break
  }
  best <- t(vapply(seq_len(N), function(u) { m <- hi[[u]]; w <- p[m[, 1]] * p[m[, 2]] * ifelse(m[, 1] == m[, 2], 1, 2); w <- w / sum(w); k <- which.max(w); c(m[k, 1], m[k, 2], w[k], nrow(m)) }, numeric(4)))
  names(hap) <- NULL
  freq <- data.frame(Haplotype = paste0("PH_", seq_along(hap)), Frequency = p, Copies = p * 2 * N, Sequence = hap, stringsAsFactors = FALSE)
  freq <- freq[order(-freq$Frequency), ]; freq <- freq[freq$Frequency > 1e-6 | seq_len(nrow(freq)) %in% match(unlist(lapply(seq_len(N), function(u) hap[best[u, 1:2]])), freq$Sequence), ]; freq$Haplotype <- paste0("PH_", seq_len(nrow(freq))); rk <- match(hap, freq$Sequence)
  ind <- data.frame(Individual = ids[use], Heterozygous_sites = nhet[use], Candidate_pairs = best[, 4], Haplotype_1 = freq$Haplotype[rk[best[, 1]]], Haplotype_2 = freq$Haplotype[rk[best[, 2]]],
                    Posterior = best[, 3], stringsAsFactors = FALSE)
  full <- function(seqs) { out <- matrix(rep(M[1, ], length(seqs)), length(seqs), ncol(M), byrow = TRUE); out[, kept] <- do.call(rbind, strsplit(seqs, "")); out }
  ph <- rbind(full(freq$Sequence[rk[best[, 1]]]), full(freq$Sequence[rk[best[, 2]]]))
  nm <- c(paste0(ids[use], "_a"), paste0(ids[use], "_b")); o <- order(rep(seq_len(N), 2), rep(1:2, each = N)); ph <- ph[o, , drop = FALSE]; rownames(ph) <- nm[o]
  list(freq = freq, ind = ind, phased = ph[, kept, drop = FALSE], kept_columns = kept, skipped = ids[skipped], iterations = it, loglik = ll,
       dropped_columns = ncol(M) - length(kept), n_ind = N)
}

# ---------------------------------------------------------------------------------------------------
# detection of loci under selection by Fst-He outliers (Beaumont and Nichols 1996, FDIST idea) for biallelic sites (SNPs).
# Fst = ratio of averages (Hudson et al. 1992; Bhatia et al. 2013) of haploid samples; the neutral expectation comes from coalescent simulations of an
# island model (one mutation per locus on a random branch, conditional on the observed sample sizes) with the migration rate tuned to the observed mean Fst.
# ---------------------------------------------------------------------------------------------------
.pg_fst_locus <- function(cnt, nn) {                    # cnt: allele counts (of variant 1) per population, nn: sample sizes
  r <- length(nn); p <- cnt / nn; hw <- 2 * p * (1 - p) * nn / (nn - 1)
  num <- 0; den <- 0
  for (i in 1:(r - 1)) for (j in (i + 1):r) { hb <- p[i] * (1 - p[j]) + p[j] * (1 - p[i]); num <- num + (hb - (hw[i] + hw[j]) / 2); den <- den + hb }
  c(num, den)
}
.pg_island_sim <- function(nn, Mig, nloc) {             # returns num, den, He per simulated locus
  D <- length(nn); res <- matrix(NA_real_, nloc, 3); ntot <- sum(nn); N <- ntot
  for (l in seq_len(nloc)) {
    deme <- rep(seq_len(D), nn); cn <- diag(D)[deme, , drop = FALSE]; k <- ntot                 # lineage descendant counts per deme
    bl <- numeric(0); bc <- list(); t0 <- rep(0, k); tt <- 0; alive <- rep(TRUE, k); ds <- deme
    cnts <- lapply(seq_len(k), function(i) cn[i, ]); tbirth <- rep(0, k); nlin <- k
    lin_d <- deme; lin_c <- cnts; lin_t <- tbirth
    while (length(lin_d) > 1) {
      kd <- tabulate(lin_d, D); rc <- kd * (kd - 1) / 2; rm <- length(lin_d) * Mig / 2; tot <- sum(rc) + rm
      tt <- tt + stats::rexp(1, tot)
      if (stats::runif(1) * tot < rm) {                                                          # migration of a random lineage
        i <- sample.int(length(lin_d), 1); lin_d[i] <- sample(setdiff(seq_len(D), lin_d[i]), 1)
      } else {                                                                                   # coalescence within a deme
        d <- if (sum(rc) > 0) sample.int(D, 1, prob = rc) else 1; ii <- which(lin_d == d); ab <- sample(ii, 2)
        for (q in ab) { bl <- c(bl, tt - lin_t[q]); bc[[length(bc) + 1]] <- lin_c[[q]] }
        lin_c[[ab[1]]] <- lin_c[[ab[1]]] + lin_c[[ab[2]]]; lin_t[ab[1]] <- tt
        lin_d <- lin_d[-ab[2]]; lin_c <- lin_c[-ab[2]]; lin_t <- lin_t[-ab[2]]
      }
    }
    j <- sample.int(length(bl), 1, prob = bl); cnt <- bc[[j]]
    fl <- .pg_fst_locus(cnt, nn); pt <- sum(cnt) / ntot
    res[l, ] <- c(fl, 2 * pt * (1 - pt) * ntot / (ntot - 1))
  }
  res
}
pg_fst_outlier <- function(M, pop, cols = NULL, nsim = 10000, mincount = 1, progress = NULL, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (is.null(cols)) cols <- seq_len(ncol(M))
  ok <- colSums(!(M == "A" | M == "C" | M == "G" | M == "T")) == 0; M <- M[, ok, drop = FALSE]; cols <- cols[ok]
  pop <- as.character(pop); tp <- table(pop); keep <- names(tp)[tp >= 2]; if (length(keep) < 2) return(NULL)
  sel <- pop %in% keep; M <- M[sel, , drop = FALSE]; pop <- pop[sel]; nn <- as.integer(table(factor(pop, keep)))
  nal <- apply(M, 2, function(x) length(unique(x))); bi <- which(nal == 2); if (length(bi) < 1) return(NULL)
  cn <- t(vapply(bi, function(j) { x <- M[, j]; v <- x != x[1]; vapply(keep, function(p) sum(v[pop == p]), 1) }, numeric(length(keep))))
  tot <- rowSums(cn); mac <- pmin(tot, sum(nn) - tot); sn <- which(mac >= mincount); if (!length(sn)) return(NULL)
  cn <- cn[sn, , drop = FALSE]; bi <- bi[sn]
  ob <- t(apply(cn, 1, function(c) { f <- .pg_fst_locus(c, nn); pt <- sum(c) / sum(nn); c(f, 2 * pt * (1 - pt) * sum(nn) / (sum(nn) - 1)) }))
  Fobs <- ob[, 1] / ob[, 2]; Fmean <- sum(ob[, 1]) / sum(ob[, 2]); He <- ob[, 3]
  # tune the migration rate: mean simulated Fst = observed mean Fst
  mf <- function(lm) { s <- .pg_island_sim(nn, exp(lm), 300); sum(s[, 1]) / sum(s[, 2]) }
  lo <- log(0.01); hi <- log(1000); target <- max(Fmean, 0.005)
  if (!is.null(progress)) progress(0.02, "calibrating the migration rate")
  for (i in 1:9) { mid <- (lo + hi) / 2; if (mf(mid) > target) lo <- mid else hi <- mid; if (!is.null(progress)) progress(0.02 + 0.1 * i / 9, "calibrating the migration rate") }
  Mig <- exp((lo + hi) / 2); sim <- NULL; chunk <- 500
  while (is.null(sim) || nrow(sim) < nsim) {
    sim <- rbind(sim, .pg_island_sim(nn, Mig, min(chunk, nsim - if (is.null(sim)) 0 else nrow(sim))))
    if (!is.null(progress)) progress(0.12 + 0.86 * nrow(sim) / nsim, "simulating the neutral island model")
  }
  Fs <- sim[, 1] / sim[, 2]; Fs[!is.finite(Fs)] <- NA; Hs <- sim[, 3]
  pv <- vapply(seq_along(Fobs), function(i) { w <- order(abs(Hs - He[i]))[seq_len(max(50, ceiling(0.05 * nsim)))]; (sum(Fs[w] >= Fobs[i], na.rm = TRUE) + 1) / (sum(!is.na(Fs[w])) + 1) }, 0)
  pl <- vapply(seq_along(Fobs), function(i) { w <- order(abs(Hs - He[i]))[seq_len(max(50, ceiling(0.05 * nsim)))]; (sum(Fs[w] <= Fobs[i], na.rm = TRUE) + 1) / (sum(!is.na(Fs[w])) + 1) }, 0)
  bins <- seq(0, 0.5, length.out = 11); mid <- (bins[-1] + bins[-11]) / 2
  env <- t(vapply(1:10, function(b) { w <- Hs >= bins[b] & Hs < bins[b + 1] + (b == 10) * 1e-9; if (sum(w) < 20) c(NA, NA, NA) else stats::quantile(Fs[w], c(0.5, 0.95, 0.99), na.rm = TRUE) }, numeric(3)))
  list(table = data.frame(Column = cols[bi], Fst = Fobs, He = He, P_high = pv, P_low = pl, Outlier = ifelse(pv < 0.05, "positive (divergent selection?)", ifelse(pl < 0.05, "negative (balancing selection?)", "")), stringsAsFactors = FALSE),
       envelope = data.frame(He = mid, median = env[, 1], q95 = env[, 2], q99 = env[, 3]), sim = data.frame(He = Hs, Fst = Fs), mean_fst = Fmean, migration = Mig, npop = length(keep), pops = keep, n = nn, nsim = nsim)
}
