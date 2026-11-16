# ==============================================================================
# SangeRly - "Barcode gap" tab
#
# Groups: every sequence has a SPECIES, a GENUS and a ROLE (Ingroup / Outgroup). They are taken from the
# sequence names (Genus_species after the accession) and can be edited in the table.
#
# Pair categories (every pairwise distance falls in exactly one):
#   Intraspecific : same species
#   Congeneric    : same genus, different species
#   Intergeneric  : different genus (an outgroup of another genus therefore counts here)
# The pairs of an outgroup sequence are used only when "Include the outgroup" is ticked.
#
# Barcode gap: smallest distance to another species > largest distance inside the species.
# By default the outgroup is left out of the gap (switch in the options).
#
# Statistics: Kruskal-Wallis + pairwise Wilcoxon (Holm) on the pair distances, Mann-Whitney intra vs inter
# (Cliff's delta, AUC), PERMANOVA and ANOSIM on the distance matrix (sequence level, base R, no extra package).
# Figures: histogram / density, box or violin plot, PCoA, distance heat map, species mean matrix,
# threshold curves, nearest-neighbour tree. Also: identification success (nearest neighbour and best close match).
# Needs: ape (already used by the app), ggplot2 >= 3.3, DT. App helpers used: .card(), .infobox(), .sg_xlsx().
# ==============================================================================

.bc_cat_cols <- c(Intraspecific = "#2c7be5", Congeneric = "#f5a623", Intergeneric = "#e5484d")

.bc_nz <- function(x, d) {
  if (is.null(x) || length(x) == 0L || (length(x) == 1L && (is.na(x) || identical(x, "")))) d else x
}

# ---- groups ------------------------------------------------------------------------------------------------
# Name of the group of each sequence. Leading accession numbers (MG272989_, NC_012920.1_, ...) are removed first.
.bc_group <- function(nm, rule = "species") {
  x <- sub("^[A-Z]{1,6}_?[0-9]{3,}(\\.[0-9]+)?[_ ]+", "", nm)
  tok <- strsplit(x, "[_ ]+")
  out <- vapply(tok, function(t) {
    t <- t[nzchar(t)]
    if (!length(t) || !grepl("^[A-Za-z]", t[1])) return("")   # no taxon name (e.g. "1392-2958"): the full sequence name is used
    switch(rule, genus = t[1], first = t[1], paste(utils::head(t, 2), collapse = "_"))
  }, "")
  out[!nzchar(out)] <- nm[!nzchar(out)]
  out
}

.bc_role <- function(x) {
  x <- as.character(x)
  ifelse(!is.na(x) & grepl("^\\s*o", x, ignore.case = TRUE), "Outgroup", "Ingroup")
}

# g: data.frame(Sequence, Species, Genus, Role). txt: character vector of choices. A choice equal to a sequence name marks
# exactly that sequence; any other text is a species, a genus or part of a sequence name (commas separate several).
.bc_set_outgroup <- function(g, txt) {
  tok <- trimws(as.character(txt)); tok <- tok[!is.na(tok) & nzchar(tok)]
  if (!length(tok)) return(g)
  hit <- g$Sequence %in% tok
  free <- trimws(unlist(strsplit(setdiff(tok, g$Sequence), "[,;]")))
  free <- tolower(gsub("[ ]+", "_", free[nzchar(free)]))
  if (length(free)) {
    sp <- tolower(gsub("[ ]+", "_", g$Species)); ge <- tolower(gsub("[ ]+", "_", g$Genus)); nm <- tolower(gsub("[ ]+", "_", g$Sequence))
    hit <- hit | Reduce(`|`, lapply(free, function(t) sp == t | ge == t | grepl(t, nm, fixed = TRUE)))
  }
  g$Role[hit] <- "Outgroup"
  g
}

# Palettes the user can choose for the groups (PCoA, tree)
.bc_pal_choices <- c("Dark 3", "Set 2", "Paired", "Pastel 1", "Dynamic", "Viridis", "Zissou 1", "Okabe-Ito", "Tableau 10")
.bc_cols <- function(lv, pal = "Dark 3") {
  lv <- unique(as.character(lv)); n <- max(length(lv), 1L)
  v <- if (pal %in% c("Okabe-Ito", "Tableau 10")) rep_len(unname(grDevices::palette.colors(palette = pal)), n)
       else grDevices::hcl.colors(n, pal)
  stats::setNames(v[seq_along(lv)], lv)
}

# Colours the user can pick for every figure: key, label on the colour box, default.
.bc_cc <- function(key, label, default) list(key = key, label = label, default = default)
.bc_cat_cc <- lapply(names(.bc_cat_cols), function(k) .bc_cc(k, k, unname(.bc_cat_cols[k])))
.bc_figcols <- list(
  dens = .bc_cat_cc, box = .bc_cat_cc,
  scat = list(.bc_cc("gap", "Barcode gap", "#27ae60"), .bc_cc("overlap", "Overlap", "#c0392b")),
  thr = list(.bc_cc("split", "False split", "#2c7be5"), .bc_cc("lump", "False lump", "#e5484d")),
  mat = list(.bc_cc("low", "Low", "#FDE725"), .bc_cc("high", "High", "#440154")))

# A colour box (browser colour picker): the input stays NULL until the user picks a colour, then the server uses it.
.bc_colorbox <- function(id, label, value) {
  shiny::div(style = "display:flex;flex-direction:column;align-items:center;font-size:12px;max-width:92px;text-align:center",
    shiny::tags$input(type = "color", id = id, value = value, title = label,
                      style = "width:46px;height:30px;padding:1px;border:1px solid #b8c2cc;border-radius:6px;cursor:pointer;background:#fff",
                      onchange = sprintf("Shiny.setInputValue('%s', this.value)", id)),
    shiny::span(style = "margin-top:2px;color:#52606d", label))
}

# The bar above each figure: its colours, palette or scale, font size and a Reset button.
.bc_figbar <- function(ns, fig) {
  boxes <- lapply(.bc_figcols[[fig]], function(cc) .bc_colorbox(ns(paste0("c_", fig, "_", cc$key)), cc$label, cc$default))
  extra <- switch(fig,
    pcoa = shiny::selectInput(ns("pal_pcoa"), "Palette", .bc_pal_choices, "Dark 3", width = "130px"),
    tree = shiny::tagList(
      shiny::selectInput(ns("pal_tree"), "Palette", .bc_pal_choices, "Dark 3", width = "130px"),
      shiny::div(style = "width:150px", shiny::sliderInput(ns("tree_cex"), "Label size", 0.4, 1.6, 0.9, step = 0.1, ticks = FALSE))),
    mat = shiny::selectInput(ns("mat_scale"), "Colour scale", c("Viridis" = "viridis", "Magma" = "magma", "Plasma" = "plasma", "Two colours (pick)" = "two"),
                             "viridis", width = "150px"),
    NULL)
  shiny::div(style = "display:flex;gap:14px;flex-wrap:wrap;align-items:flex-end;padding:8px 10px 2px;margin-bottom:8px;background:#f6f8fa;border-radius:8px",
    boxes, extra,
    if (!identical(fig, "tree"))
      shiny::div(style = "width:150px", shiny::sliderInput(ns(paste0("fs_", fig)), "Font size", 8, 28, if (identical(fig, "mat")) 12 else 13, step = 1, ticks = FALSE)),
    shiny::div(style = "margin-bottom:15px", shiny::actionButton(ns(paste0("reset_", fig)), "Reset", icon = shiny::icon("rotate-left"), class = "btn-sm")))
}

# ---- distances ---------------------------------------------------------------------------------------------
# Distance matrix of an aligned named character vector. model: "K80" (K2P), "raw" (p-distance), "JC69", "TN93".
.bc_dist <- function(aln, model = "K80") {
  if (!requireNamespace("ape", quietly = TRUE)) stop("The 'ape' package is required: install.packages(\"ape\")")
  m <- do.call(rbind, strsplit(tolower(unname(toupper(aln))), ""))
  rownames(m) <- names(aln)
  D <- ape::dist.dna(ape::as.DNAbin(m), model = model, pairwise.deletion = TRUE, as.matrix = TRUE)
  D[!is.finite(D)] <- NA_real_
  D[!is.na(D) & D <= 0] <- 0        # K2P can give -0 / tiny negative values for identical sequences
  D
}

# Rows to keep so that the matrix has no NA (PCoA, PERMANOVA and the tree need a complete matrix):
# the sequence with most NA is dropped, repeatedly.
.bc_complete <- function(D) {
  keep <- seq_len(nrow(D))
  repeat {
    if (length(keep) < 3L) break
    sub <- D[keep, keep, drop = FALSE]
    bad <- rowSums(is.na(sub))
    if (!any(bad > 0)) break
    keep <- keep[-which.max(bad)]
  }
  keep
}

# ---- pairs, per-species gap, summary -----------------------------------------------------------------------
# D: distance matrix; sp, ge, ro: species, genus, role of each row of D.
.bc_analyse <- function(D, sp, ge, ro, model = "K80", gap_out = FALSE) {
  nm <- rownames(D)
  ix <- which(lower.tri(D), arr.ind = TRUE)
  i <- ix[, 2]; j <- ix[, 1]
  pairs <- data.frame(Sequence_1 = nm[i], Sequence_2 = nm[j], Species_1 = sp[i], Species_2 = sp[j],
                      Genus_1 = ge[i], Genus_2 = ge[j], Role_1 = ro[i], Role_2 = ro[j],
                      Distance = D[ix], stringsAsFactors = FALSE)
  pairs$Category <- ifelse(pairs$Species_1 == pairs$Species_2, "Intraspecific",
                      ifelse(pairs$Genus_1 == pairs$Genus_2, "Congeneric", "Intergeneric"))
  use <- if (gap_out) rep(TRUE, length(sp)) else ro != "Outgroup"
  pairs$Use <- use[i] & use[j]

  pu <- pairs[pairs$Use & !is.na(pairs$Distance), , drop = FALSE]
  sl <- unique(sp[use])
  rows <- lapply(sl, function(g) {
    n <- sum(sp == g & use)
    gen <- ge[which(sp == g & use)[1]]
    intra <- pu$Distance[pu$Species_1 == g & pu$Species_2 == g]
    inter <- pu[xor(pu$Species_1 == g, pu$Species_2 == g), , drop = FALSE]
    other <- ifelse(inter$Species_1 == g, inter$Species_2, inter$Species_1)
    k <- if (nrow(inter)) which.min(inter$Distance) else NA_integer_
    cong <- inter[inter$Genus_1 == inter$Genus_2, , drop = FALSE]
    mx <- if (length(intra)) max(intra) else NA_real_
    mn <- if (!is.na(k)) inter$Distance[k] else NA_real_
    data.frame(Species = g, Genus = gen, N = n,
               Max_intra = mx,
               Mean_intra = if (length(intra)) mean(intra) else NA_real_,
               Nearest_species = if (!is.na(k)) other[k] else NA_character_,
               Min_inter = mn,
               Min_congeneric = if (nrow(cong)) min(cong$Distance) else NA_real_,
               Gap = if (is.na(mx) || is.na(mn)) NA_real_ else mn - mx,
               Verdict = if (is.na(mn)) "no other species"
                         else if (is.na(mx)) "single sequence (no intraspecific distance)"
                         else if (mn > mx) "barcode gap" else "overlap (no gap)",
               stringsAsFactors = FALSE)
  })
  species <- do.call(rbind, rows)

  same <- pairs$Species_1 == pairs$Species_2
  ok <- pairs$Use & !is.na(pairs$Distance)
  intra_all <- pairs$Distance[ok & same]; inter_all <- pairs$Distance[ok & !same]
  rng <- function(f, x) if (length(x)) f(x) else NA_real_
  cat_mean <- vapply(names(.bc_cat_cols), function(k) {
    x <- pairs$Distance[pairs$Use & pairs$Category == k & !is.na(pairs$Distance)]
    if (length(x)) mean(x) else NA_real_
  }, 0)
  summ <- list(model = model, n_seq = length(nm), n_species = length(sl),
               n_genera = length(unique(ge[use])), n_out = sum(ro == "Outgroup"),
               n_intra = length(intra_all), n_inter = length(inter_all),
               mean_intra = rng(mean, intra_all), max_intra = rng(max, intra_all),
               mean_inter = rng(mean, inter_all), min_inter = rng(min, inter_all),
               cat_mean = cat_mean, gap_out = gap_out)
  summ$gap <- if (is.na(summ$max_intra) || is.na(summ$min_inter)) NA_real_ else summ$min_inter - summ$max_intra
  summ$n_with_gap <- sum(species$Verdict == "barcode gap")
  summ$n_testable <- sum(species$Verdict %in% c("barcode gap", "overlap (no gap)"))
  list(pairs = pairs, species = species, summary = summ)
}

# ---- statistics --------------------------------------------------------------------------------------------
# Cliff's delta of a against b (-1 .. 1; positive: a tends to be larger), from ranks.
.bc_cliff <- function(a, b) {
  n1 <- length(a); n2 <- length(b)
  r <- rank(c(a, b))
  U <- sum(r[seq_len(n1)]) - n1 * (n1 + 1) / 2
  2 * U / (n1 * n2) - 1
}

.bc_catstats <- function(pairs) {
  rows <- lapply(names(.bc_cat_cols), function(k) {
    x <- pairs$Distance[pairs$Use & pairs$Category == k & !is.na(pairs$Distance)]
    if (!length(x)) return(NULL)
    data.frame(Category = k, N = length(x), Mean = mean(x), SD = if (length(x) > 1) stats::sd(x) else NA_real_,
               Median = stats::median(x), Min = min(x), Max = max(x), stringsAsFactors = FALSE)
  })
  rows <- rows[!vapply(rows, is.null, TRUE)]
  if (length(rows)) do.call(rbind, rows) else NULL
}

# Tests on the pair distances. NOTE: pairs share sequences, so they are not independent; the P values are
# optimistic. The sequence-level tests (PERMANOVA, ANOSIM) below are the safer ones.
.bc_tests <- function(pairs) {
  out <- list(kw = NULL, pw = NULL, ivi = NULL)
  d <- pairs[pairs$Use & !is.na(pairs$Distance), , drop = FALSE]
  d$Category <- factor(d$Category, levels = names(.bc_cat_cols))
  tab <- table(d$Category)
  d <- d[d$Category %in% names(tab)[tab >= 2], , drop = FALSE]
  d$Category <- droplevels(d$Category)
  if (nlevels(d$Category) >= 2) {
    k <- tryCatch(stats::kruskal.test(Distance ~ Category, data = d), error = function(e) NULL)
    if (!is.null(k))
      out$kw <- data.frame(Test = "Kruskal-Wallis (all categories)", Statistic = unname(k$statistic),
                           df = unname(k$parameter), P = k$p.value, stringsAsFactors = FALSE)
    pw <- tryCatch(suppressWarnings(stats::pairwise.wilcox.test(d$Distance, d$Category, p.adjust.method = "holm", exact = FALSE)),
                   error = function(e) NULL)
    if (!is.null(pw)) {
      pm <- pw$p.value; rows <- list()
      for (a in seq_len(nrow(pm))) for (b in seq_len(ncol(pm))) if (!is.na(pm[a, b])) {
        x <- d$Distance[d$Category == rownames(pm)[a]]; y <- d$Distance[d$Category == colnames(pm)[b]]
        rows[[length(rows) + 1L]] <- data.frame(
          Comparison = paste(rownames(pm)[a], "vs", colnames(pm)[b]), n1 = length(x), n2 = length(y),
          Median_1 = stats::median(x), Median_2 = stats::median(y), Cliffs_delta = .bc_cliff(x, y),
          P_holm = pm[a, b], stringsAsFactors = FALSE)
      }
      if (length(rows)) out$pw <- do.call(rbind, rows)
    }
  }
  ok <- pairs$Use & !is.na(pairs$Distance)
  same <- pairs$Species_1 == pairs$Species_2
  x <- pairs$Distance[ok & same]; y <- pairs$Distance[ok & !same]
  if (length(x) >= 2 && length(y) >= 2) {
    w <- suppressWarnings(stats::wilcox.test(x, y, exact = FALSE))
    n1 <- length(x); n2 <- length(y)
    out$ivi <- data.frame(Test = "Mann-Whitney: intraspecific vs interspecific", n_intra = n1, n_inter = n2,
                          Median_intra = stats::median(x), Median_inter = stats::median(y),
                          W = unname(w$statistic), P = w$p.value, Cliffs_delta = .bc_cliff(x, y),
                          AUC = 1 - unname(w$statistic) / (n1 * n2), stringsAsFactors = FALSE)
  }
  out
}

# PERMANOVA (Anderson 2001) on a complete distance matrix; grp: group of every row.
.bc_permanova <- function(D, grp, nperm = 999, seed = 1) {
  grp <- factor(grp); N <- nrow(D); a <- nlevels(grp)
  if (a < 2L || N <= a) return(NULL)
  A <- D^2
  sst <- sum(A[lower.tri(A)]) / N
  fstat <- function(g) {
    ssw <- 0
    for (l in levels(g)) {
      k <- which(g == l)
      if (length(k) > 1L) { s <- A[k, k, drop = FALSE]; ssw <- ssw + sum(s[lower.tri(s)]) / length(k) }
    }
    ssa <- sst - ssw
    c(F = (ssa / (a - 1)) / (ssw / (N - a)), R2 = ssa / sst)
  }
  obs <- fstat(grp)
  set.seed(seed)
  perm <- vapply(seq_len(nperm), function(i) fstat(sample(grp))[["F"]], 0)
  list(stat = obs[["F"]], R2 = obs[["R2"]], P = (sum(perm >= obs[["F"]]) + 1) / (nperm + 1))
}

# ANOSIM (Clarke 1993): R close to 1 = groups separated, close to 0 = no difference.
.bc_anosim <- function(D, grp, nperm = 999, seed = 1) {
  gi <- as.integer(factor(grp)); N <- nrow(D)
  if (max(gi) < 2L || N <= max(gi)) return(NULL)
  lt <- lower.tri(D)
  r <- rank(D[lt]); M <- length(r)
  ii <- row(D)[lt]; jj <- col(D)[lt]
  rstat <- function(g) {
    w <- g[ii] == g[jj]
    if (!any(w) || all(w)) return(NA_real_)
    (mean(r[!w]) - mean(r[w])) / (M / 2)
  }
  obs <- rstat(gi)
  if (is.na(obs)) return(NULL)
  set.seed(seed)
  perm <- vapply(seq_len(nperm), function(i) rstat(sample(gi)), 0)
  list(stat = obs, P = (sum(perm >= obs, na.rm = TRUE) + 1) / (nperm + 1))
}

.bc_permtests <- function(Dc, sp, ge, ro, nperm = 999, seed = 1) {
  rows <- list()
  add <- function(label, grp) {
    grp <- factor(grp)
    if (nlevels(grp) < 2L) return(invisible())
    pm <- .bc_permanova(Dc, grp, nperm, seed); an <- .bc_anosim(Dc, grp, nperm, seed)
    if (!is.null(pm))
      rows[[length(rows) + 1L]] <<- data.frame(Test = "PERMANOVA", Grouping = label, N = nrow(Dc), Groups = nlevels(grp),
        Statistic = pm$stat, R2 = pm$R2, P = pm$P, Permutations = nperm, stringsAsFactors = FALSE)
    if (!is.null(an))
      rows[[length(rows) + 1L]] <<- data.frame(Test = "ANOSIM", Grouping = label, N = nrow(Dc), Groups = nlevels(grp),
        Statistic = an$stat, R2 = NA_real_, P = an$P, Permutations = nperm, stringsAsFactors = FALSE)
  }
  add("Species", sp); add("Genus", ge); add("Ingroup vs outgroup", ro)
  if (length(rows)) do.call(rbind, rows) else NULL
}

# ---- PCoA, threshold, identification, tree ------------------------------------------------------------------
.bc_pcoa <- function(Dc) {
  n <- nrow(Dc); k <- min(3L, n - 1L)
  if (k < 2L) return(NULL)
  cm <- tryCatch(suppressWarnings(stats::cmdscale(stats::as.dist(Dc), k = k, eig = TRUE)), error = function(e) NULL)
  if (is.null(cm)) return(NULL)
  pos <- sum(cm$eig[cm$eig > 0]); neg <- sum(abs(cm$eig[cm$eig < 0]))
  co <- cm$points; colnames(co) <- paste0("PCoA", seq_len(ncol(co)))
  list(coords = co, var = cm$eig[seq_len(ncol(co))] / pos * 100, neg_pct = neg / (pos + neg) * 100)
}

# Distance threshold that best separates conspecific from heterospecific pairs.
# False split = conspecific pairs farther apart than the threshold; false lump = heterospecific pairs at or below it.
.bc_threshold <- function(pairs) {
  ok <- pairs$Use & !is.na(pairs$Distance)
  same <- pairs$Species_1 == pairs$Species_2
  intra <- pairs$Distance[ok & same]; inter <- pairs$Distance[ok & !same]
  if (length(intra) < 1L || length(inter) < 1L) return(NULL)
  allv <- sort(unique(c(0, intra, inter)))
  if (length(allv) > 400L) allv <- unique(stats::quantile(allv, seq(0, 1, length.out = 400), names = FALSE))
  fn <- vapply(allv, function(t) mean(intra > t), 0)
  fp <- vapply(allv, function(t) mean(inter <= t), 0)
  tot <- fn + fp
  b <- which.min(tot)
  list(curve = data.frame(Threshold = allv, False_split = fn, False_lump = fp, Total = tot),
       best = allv[b], fn = fn[b], fp = fp[b])
}

# Nearest neighbour (NN) and best close match (BCM, with a threshold) identification, ingroup only.
.bc_ident <- function(D, sp, ro, thr = NA_real_) {
  idx <- which(ro != "Outgroup")
  if (length(idx) < 3L) return(NULL)
  D2 <- D[idx, idx, drop = FALSE]; s2 <- sp[idx]; nm <- rownames(D2)
  rows <- lapply(seq_along(idx), function(q) {
    d <- D2[q, ]; d[q] <- NA_real_
    if (all(is.na(d)))
      return(data.frame(Sequence = nm[q], Species = s2[q], Nearest_distance = NA_real_, Nearest_species = "",
                        NN_result = "No comparison", BCM_result = "No comparison", stringsAsFactors = FALSE))
    mn <- min(d, na.rm = TRUE)
    nn <- which(!is.na(d) & d <= mn + 1e-10)
    sps <- unique(s2[nn])
    cls <- if (sum(s2 == s2[q]) == 1L) "Singleton (not testable)"
           else if (!(s2[q] %in% sps)) "Incorrect"
           else if (length(sps) > 1L) "Ambiguous" else "Correct"
    bcm <- if (cls == "Singleton (not testable)") cls else if (!is.na(thr) && mn > thr) "No match" else cls
    data.frame(Sequence = nm[q], Species = s2[q], Nearest_distance = mn, Nearest_species = paste(sps, collapse = "; "),
               NN_result = cls, BCM_result = bcm, stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

.bc_ident_summary <- function(id) {
  lv <- c("Correct", "Ambiguous", "Incorrect", "No match", "Singleton (not testable)")
  one <- function(col, label) {
    t <- table(factor(id[[col]], levels = lv))
    testable <- sum(t[c("Correct", "Ambiguous", "Incorrect", "No match")])
    data.frame(Method = label, Correct = unname(t["Correct"]), Ambiguous = unname(t["Ambiguous"]),
               Incorrect = unname(t["Incorrect"]), `No match` = unname(t["No match"]),
               `Singleton (not testable)` = unname(t["Singleton (not testable)"]),
               Success_pct = if (testable > 0) round(100 * unname(t["Correct"]) / testable, 1) else NA_real_,
               check.names = FALSE, stringsAsFactors = FALSE)
  }
  rbind(one("NN_result", "Nearest neighbour"), one("BCM_result", "Best close match (with threshold)"))
}

# Neighbour-joining tree of the complete distance matrix (as in the Phylogeny tab: negative branches set to 0).
# root: "outgroup" (when the matrix has outgroup sequences), "midpoint" or "none".
.bc_tree <- function(res, root = "outgroup") {
  Dc <- res$Dc
  if (nrow(Dc) < 3L) return(NULL)
  tr <- ape::nj(stats::as.dist(Dc))
  tr$edge.length[tr$edge.length < 0] <- 0
  out <- rownames(Dc)[res$role[res$keep] == "Outgroup"]
  if (identical(root, "outgroup") && length(out) && length(out) < nrow(Dc)) {
    tr <- tryCatch(ape::root(tr, outgroup = out, resolve.root = TRUE), error = function(e) tr)
  } else if (identical(root, "midpoint") && requireNamespace("phangorn", quietly = TRUE)) {
    tr <- tryCatch(phangorn::midpoint(tr), error = function(e) tr)
  }
  tr
}

# Height (pixels) of the tree plot: same rule as the Phylogeny tab (22 px per tip).
.bc_tree_h <- function(tr, cex = 0.9) if (is.null(tr)) 120 else max(380, 22 * ape::Ntip(tr) * max(0.6, cex) + 80)

# Is every species (or genus) monophyletic in the tree? lab: named vector (tip label -> group).
.bc_mono <- function(tr, lab) {
  rows <- lapply(unique(lab), function(g) {
    tips <- names(lab)[lab == g]
    if (length(tips) < 2L) return(NULL)
    mono <- tryCatch(ape::is.monophyletic(tr, tips), error = function(e) NA)
    data.frame(Group = g, N = length(tips), Monophyletic = if (is.na(mono)) "?" else if (mono) "yes" else "no",
               stringsAsFactors = FALSE)
  })
  rows <- rows[!vapply(rows, is.null, TRUE)]
  if (length(rows)) do.call(rbind, rows) else NULL
}

# ---- everything at once ------------------------------------------------------------------------------------
.bc_run <- function(D, sp, ge, ro, model = "K80", gap_out = FALSE, nperm = 999, seed = 1) {
  an <- .bc_analyse(D, sp, ge, ro, model, gap_out)
  kp <- .bc_complete(D)
  Dc <- D[kp, kp, drop = FALSE]
  c(an, list(D = D, Dc = Dc, sp = sp, ge = ge, role = ro, keep = kp,
             dropped = rownames(D)[setdiff(seq_len(nrow(D)), kp)],
             catstats = .bc_catstats(an$pairs),
             tests = .bc_tests(an$pairs),
             perm = if (length(kp) >= 4L) .bc_permtests(Dc, sp[kp], ge[kp], ro[kp], nperm, seed) else NULL,
             pcoa = if (length(kp) >= 3L) .bc_pcoa(Dc) else NULL,
             thr = .bc_threshold(an$pairs),
             nperm = nperm, model = model))
}

# ---- figures -----------------------------------------------------------------------------------------------
.bc_theme <- function(fs = 13) list(ggplot2::theme_minimal(base_size = fs), ggplot2::theme(legend.position = "top"))

.bc_plot_dens <- function(pairs, mode = "dens", adjust = 1, model = "K80", thr = NULL, cols = .bc_cat_cols, fs = 13) {
  d <- pairs[pairs$Use & !is.na(pairs$Distance), , drop = FALSE]
  tab <- table(d$Category)
  d$Category <- factor(d$Category, levels = names(.bc_cat_cols))
  if (mode != "hist") d <- d[d$Category %in% names(tab)[tab >= 2], , drop = FALSE]
  p <- ggplot2::ggplot(d, ggplot2::aes(Distance, fill = Category))
  if (identical(mode, "hist")) {
    p <- p + ggplot2::geom_histogram(bins = 40, alpha = 0.65, position = "identity", colour = "white", linewidth = 0.2) +
      ggplot2::labs(y = "Number of pairs")
  } else if (identical(mode, "dens")) {
    p <- p + ggplot2::geom_density(ggplot2::aes(colour = Category), alpha = 0.35, adjust = adjust, linewidth = 0.7) +
      ggplot2::labs(y = "Density")
  } else {
    p <- p + ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)), bins = 40, alpha = 0.3,
                                     position = "identity", colour = NA) +
      ggplot2::geom_density(ggplot2::aes(colour = Category), fill = NA, adjust = adjust, linewidth = 0.8) +
      ggplot2::labs(y = "Density")
  }
  p <- p + ggplot2::scale_fill_manual(values = cols, name = NULL, drop = TRUE) +
    ggplot2::scale_colour_manual(values = cols, name = NULL, drop = TRUE) +
    ggplot2::labs(x = paste0("Distance (", model, ")")) + .bc_theme(fs)
  if (!is.null(thr) && is.finite(thr))
    p <- p + ggplot2::geom_vline(xintercept = thr, linetype = 2, colour = "grey25") +
      ggplot2::annotate("text", x = thr, y = Inf, label = sprintf(" threshold %.4f", thr), hjust = 0, vjust = 1.6, size = fs * 0.8 / ggplot2::.pt)
  p
}

.bc_plot_box <- function(pairs, kw_p = NULL, violin = FALSE, pts = TRUE, model = "K80", cols = .bc_cat_cols, fs = 13) {
  d <- pairs[pairs$Use & !is.na(pairs$Distance), , drop = FALSE]
  lev <- names(.bc_cat_cols)[names(.bc_cat_cols) %in% d$Category]
  d$Category <- factor(d$Category, levels = lev)
  n <- table(d$Category)
  p <- ggplot2::ggplot(d, ggplot2::aes(Category, Distance, fill = Category))
  if (violin) {
    big <- names(n)[n >= 2]
    p <- p + ggplot2::geom_violin(data = d[d$Category %in% big, , drop = FALSE], alpha = 0.6, scale = "width", trim = TRUE) +
      ggplot2::geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", alpha = 0.9)
  } else {
    p <- p + ggplot2::geom_boxplot(alpha = 0.7, outlier.shape = NA)
  }
  if (pts) p <- p + ggplot2::geom_jitter(width = 0.18, size = 0.8, alpha = 0.35)
  p + ggplot2::scale_fill_manual(values = cols[lev], guide = "none") +
    ggplot2::scale_x_discrete(labels = paste0(lev, "\n(n = ", as.integer(n[lev]), ")")) +
    ggplot2::labs(x = NULL, y = paste0("Distance (", model, ")"),
                  subtitle = if (!is.null(kw_p) && is.finite(kw_p)) paste("Kruskal-Wallis P =", format.pval(kw_p, digits = 3)) else NULL) +
    ggplot2::theme_minimal(base_size = fs)
}

.bc_plot_pcoa <- function(pc, sp, ge, ro, colour_by = "species", axes = c(1, 2), ellipse = TRUE, labels = FALSE, pal = "Dark 3", fs = 13) {
  co <- pc$coords
  axes <- as.integer(axes)
  if (max(axes) > ncol(co)) axes <- c(1L, 2L)
  df <- data.frame(Sequence = rownames(co), X = co[, axes[1]], Y = co[, axes[2]], Role = ro, stringsAsFactors = FALSE)
  df$Group <- switch(colour_by, genus = ge, role = ro, sp)
  cols <- .bc_cols(df$Group, pal)
  p <- ggplot2::ggplot(df, ggplot2::aes(X, Y, colour = Group, shape = Role)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey85") + ggplot2::geom_vline(xintercept = 0, colour = "grey85") +
    ggplot2::geom_point(size = 3, alpha = 0.85) +
    ggplot2::scale_colour_manual(values = cols, name = NULL) +
    ggplot2::scale_shape_manual(values = c(Ingroup = 16, Outgroup = 17), name = NULL)
  if (length(unique(df$Role)) < 2L) p <- p + ggplot2::guides(shape = "none")   # only one role: the shape legend says nothing
  if (isTRUE(ellipse)) {
    tab <- table(df$Group); big <- names(tab)[tab >= 4]
    if (length(big))
      p <- p + ggplot2::stat_ellipse(data = df[df$Group %in% big, , drop = FALSE],
                                     ggplot2::aes(X, Y, colour = Group, group = Group), inherit.aes = FALSE,
                                     level = 0.95, type = "norm", linetype = 2, linewidth = 0.5, show.legend = FALSE)
  }
  if (isTRUE(labels))
    p <- p + ggplot2::geom_text(ggplot2::aes(label = Sequence), size = fs * 0.62 / ggplot2::.pt, vjust = -0.9, check_overlap = TRUE, show.legend = FALSE)
  p + ggplot2::labs(x = sprintf("PCoA %d (%.1f%%)", axes[1], pc$var[axes[1]]), y = sprintf("PCoA %d (%.1f%%)", axes[2], pc$var[axes[2]])) +
    ggplot2::coord_equal() + ggplot2::theme_minimal(base_size = fs) + ggplot2::theme(legend.position = "right")
}

# fill scale of the distance matrices: a viridis option or two colours picked by the user
.bc_fill <- function(ms, na, name) {
  if (identical(ms$scale, "two")) ggplot2::scale_fill_gradient(low = ms$low, high = ms$high, na.value = na, name = name)
  else ggplot2::scale_fill_viridis_c(option = .bc_nz(ms$scale, "viridis"), direction = -1, na.value = na, name = name)
}

.bc_plot_heat <- function(D, sp, ro, labels = TRUE, order = "group", model = "K80", ms = list(scale = "viridis"), fs = 12) {
  n <- nrow(D); nm <- rownames(D)
  ord <- if (identical(order, "cluster") && n >= 3L) stats::hclust(stats::as.dist(D), "average")$order
         else order(ro == "Outgroup", sp, seq_len(n))
  M <- D[ord, ord, drop = FALSE]
  df <- data.frame(x = rep(seq_len(n), each = n), y = n + 1 - rep(seq_len(n), times = n), Distance = as.vector(M))
  show <- isTRUE(labels) && n <= 60L
  p <- ggplot2::ggplot(df, ggplot2::aes(x, y, fill = Distance)) + ggplot2::geom_raster() +
    .bc_fill(ms, "grey85", model) +
    ggplot2::scale_x_continuous(breaks = if (show) seq_len(n) else NULL, labels = if (show) nm[ord] else NULL, expand = c(0, 0)) +
    ggplot2::scale_y_continuous(breaks = if (show) n + 1 - seq_len(n) else NULL, labels = if (show) nm[ord] else NULL, expand = c(0, 0)) +
    ggplot2::coord_equal() + ggplot2::labs(x = NULL, y = NULL) + ggplot2::theme_minimal(base_size = fs) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, size = max(3, fs * 0.5)),
                   axis.text.y = ggplot2::element_text(size = max(3, fs * 0.5)))
  if (identical(order, "group")) {
    spo <- sp[ord]; b <- which(spo[-1] != spo[-n]) + 0.5
    if (length(b)) p <- p + ggplot2::geom_vline(xintercept = b, colour = "white", linewidth = 0.3) +
      ggplot2::geom_hline(yintercept = n + 1 - b, colour = "white", linewidth = 0.3)
  }
  p
}

# Mean distance between (and within, on the diagonal) species.
.bc_spmatrix <- function(D, sp, ro, include_out = TRUE) {
  keep <- if (include_out) seq_along(sp) else which(ro != "Outgroup")
  D <- D[keep, keep, drop = FALSE]; sp <- sp[keep]
  lv <- unique(sp[order(ro[keep] == "Outgroup", sp)])
  M <- matrix(NA_real_, length(lv), length(lv), dimnames = list(lv, lv))
  for (a in seq_along(lv)) for (b in a:length(lv)) {
    ia <- which(sp == lv[a]); ib <- which(sp == lv[b])
    v <- if (a == b) { s <- D[ia, ia, drop = FALSE]; s[lower.tri(s)] } else as.vector(D[ia, ib, drop = FALSE])
    v <- v[!is.na(v)]
    if (length(v)) M[a, b] <- M[b, a] <- mean(v)
  }
  M
}

.bc_plot_spheat <- function(M, model = "K80", ms = list(scale = "viridis"), fs = 12) {
  n <- nrow(M); lv <- rownames(M)
  df <- expand.grid(i = seq_len(n), j = seq_len(n))
  df$Distance <- M[cbind(df$i, df$j)]
  df$Row <- factor(lv[df$i], levels = rev(lv)); df$Col <- factor(lv[df$j], levels = lv)
  p <- ggplot2::ggplot(df, ggplot2::aes(Col, Row, fill = Distance)) + ggplot2::geom_tile(colour = "white") +
    .bc_fill(ms, "grey90", model) +
    ggplot2::labs(x = NULL, y = NULL, subtitle = "Mean distance between species (diagonal: mean within the species)") +
    ggplot2::theme_minimal(base_size = fs) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(), axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  if (n <= 20L) {
    mid <- mean(range(df$Distance, na.rm = TRUE))
    df$lab <- ifelse(is.na(df$Distance), "", sprintf("%.3f", df$Distance))
    df$tc <- ifelse(!is.na(df$Distance) & df$Distance > mid, "white", "black")
    p <- p + ggplot2::geom_text(data = df, ggplot2::aes(label = lab, colour = tc), size = fs * 0.65 / ggplot2::.pt, show.legend = FALSE) +
      ggplot2::scale_colour_identity()
  }
  p
}

.bc_plot_scatter <- function(sp_tab, cols = c("barcode gap" = "#27ae60", "overlap" = "#c0392b"), fs = 13) {
  d <- sp_tab[!is.na(sp_tab$Max_intra) & !is.na(sp_tab$Min_inter), , drop = FALSE]
  d$Result <- ifelse(d$Gap > 0, "barcode gap", "overlap")
  lim <- c(0, max(d$Max_intra, d$Min_inter) * 1.08)
  ggplot2::ggplot(d, ggplot2::aes(Max_intra, Min_inter, colour = Result, label = Species)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey40") +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_text(vjust = -0.8, size = fs * 0.75 / ggplot2::.pt, check_overlap = TRUE, show.legend = FALSE) +
    ggplot2::scale_colour_manual(values = cols, name = NULL) +
    ggplot2::coord_cartesian(xlim = lim, ylim = lim) +
    ggplot2::labs(x = "Largest intraspecific distance", y = "Distance to the nearest other species") + .bc_theme(fs)
}

.bc_plot_thr <- function(th, user = NA_real_, model = "K80", cols = c("#2c7be5", "#e5484d"), fs = 13) {
  cv <- th$curve
  long <- rbind(data.frame(Threshold = cv$Threshold, Rate = cv$False_split, Error = "False split (same species farther than threshold)"),
                data.frame(Threshold = cv$Threshold, Rate = cv$False_lump, Error = "False lump (different species at or below threshold)"))
  p <- ggplot2::ggplot(long, ggplot2::aes(Threshold, Rate, colour = Error)) + ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::geom_vline(xintercept = th$best, linetype = 2, colour = "grey25") +
    ggplot2::annotate("text", x = th$best, y = Inf, label = sprintf(" suggested %.4f", th$best), hjust = 0, vjust = 1.6, size = fs * 0.8 / ggplot2::.pt) +
    ggplot2::scale_colour_manual(values = cols, name = NULL) +
    ggplot2::labs(x = paste0("Distance threshold (", model, ")"), y = "Proportion of pairs") + .bc_theme(fs) +
    ggplot2::guides(colour = ggplot2::guide_legend(nrow = 2))
  if (!is.na(user)) p <- p + ggplot2::geom_vline(xintercept = user, linetype = 3, colour = "#27ae60")
  p
}

# Same look as the trees of the Phylogeny tab (ape::plot.phylo, thick branches, scale bar); tips coloured by species or genus,
# the legend sits in its own panel on the right so that it never covers the tree.
.bc_draw_tree <- function(res, tr, colour_by = "species", type = "phylogram", cex = 0.9, tip_label = "sequence", title = NULL, pal = "Dark 3") {
  if (is.null(tr)) { graphics::plot.new(); graphics::text(0.5, 0.5, "At least 3 sequences with complete distances are needed."); return(invisible()) }
  Dc <- res$Dc
  grp <- if (identical(colour_by, "genus")) res$ge[res$keep] else res$sp[res$keep]
  names(grp) <- rownames(Dc)
  cols <- .bc_cols(grp, pal)
  tipcol <- unname(cols[grp[tr$tip.label]])
  if (identical(tip_label, "species")) tr$tip.label <- unname(res$sp[res$keep][match(tr$tip.label, rownames(Dc))])
  graphics::layout(matrix(1:2, 1), widths = c(4, 1)); on.exit(graphics::layout(1), add = TRUE)
  op <- graphics::par(mar = c(2, 1, 2.2, 0.5)); on.exit(graphics::par(op), add = TRUE)
  ape::plot.phylo(tr, type = type, use.edge.length = !identical(type, "cladogram"), cex = cex, edge.width = 1.6,
                  tip.color = tipcol, label.offset = if (type %in% c("phylogram", "cladogram")) 0.002 else 0, main = title)
  if (identical(type, "phylogram")) ape::add.scale.bar(cex = cex * 0.8)
  graphics::par(mar = c(2, 0, 2.2, 0))
  graphics::plot.new()
  graphics::legend("topleft", legend = names(cols), text.col = cols, pch = 15, col = cols, bty = "n", cex = max(0.6, cex * 0.85),
                   title = if (identical(colour_by, "genus")) "Genus" else "Species", title.adj = 0)
  invisible()
}

# ---- input files and group sources --------------------------------------------------------------------------
# Aligned or plain FASTA -> named character vector (U -> T, "." and "?" -> "-"). Names = the whole header line.
.bc_read_fasta <- function(path) {
  x <- readLines(path, warn = FALSE)
  x <- x[nzchar(trimws(x))]
  h <- grep("^>", x)
  if (!length(h)) return(NULL)
  id <- trimws(sub("^>", "", x[h]))
  end <- c(h[-1] - 1L, length(x))
  sq <- vapply(seq_along(h), function(i)
    if (end[i] > h[i]) paste(gsub("\\s", "", x[(h[i] + 1L):end[i]]), collapse = "") else "", "")
  stats::setNames(gsub("[.?]", "-", gsub("U", "T", toupper(sq))), make.unique(id))
}

# CSV / TSV / semicolon table with a header row.
.bc_read_table <- function(path) {
  l <- readLines(path, n = 5, warn = FALSE); l <- l[nzchar(l)]
  if (!length(l)) return(NULL)
  cnt <- function(ch) sum(gregexpr(ch, l[1], fixed = TRUE)[[1]] > 0)
  sep <- if (any(grepl("\t", l))) "\t" else if (cnt(";") > cnt(",")) ";" else ","
  tryCatch(utils::read.table(path, header = TRUE, sep = sep, quote = "\"", comment.char = "", stringsAsFactors = FALSE,
                             check.names = FALSE, fill = TRUE, strip.white = TRUE), error = function(e) NULL)
}

# Names are compared after replacing every symbol by "_" (the tree and population tabs clean the names this way).
.bc_key <- function(x) gsub("[^A-Za-z0-9]+", "_", trimws(as.character(x)))

# A data frame (sequence column + group column, optionally Genus and Role) or a named vector (names = sequences)
# -> data.frame(Sequence, Species[, Genus][, Role]).
.bc_as_groups <- function(x) {
  if (is.null(x)) return(NULL)
  if (is.atomic(x) && !is.null(names(x)))
    return(data.frame(Sequence = names(x), Species = as.character(x), stringsAsFactors = FALSE))
  if (!is.data.frame(x) || ncol(x) < 2L) return(NULL)
  nmx <- tolower(trimws(names(x)))
  pick <- function(cands, excl = integer()) which(nmx %in% cands & !(seq_along(nmx) %in% excl))[1]
  i <- pick(c("sequence", "id", "name", "label", "sample", "tip", "seq", "taxon", "accession", "individual"))
  if (is.na(i)) i <- 1L
  j <- pick(c("species", "group", "population", "pop", "cluster", "clade", "deme", "lineage", "haplogroup"), i)
  if (is.na(j)) j <- setdiff(seq_len(ncol(x)), i)[1]
  out <- data.frame(Sequence = as.character(x[[i]]), Species = as.character(x[[j]]), stringsAsFactors = FALSE)
  g <- pick("genus", c(i, j)); if (!is.na(g)) out$Genus <- as.character(x[[g]])
  r <- pick(c("role", "outgroup"), c(i, j)); if (!is.na(r)) out$Role <- .bc_role(x[[r]])
  out
}

# Clades of a tree: tips linked by branches with a patristic distance <= thr (single linkage) form one group.
.bc_tree_groups <- function(tr, thr) {
  if (is.null(tr) || ape::Ntip(tr) < 2L) return(NULL)
  cl <- stats::cutree(stats::hclust(stats::as.dist(ape::cophenetic.phylo(tr)), "single"), h = thr)
  data.frame(Sequence = names(cl), Species = sprintf("Clade_%02d", as.integer(cl)), stringsAsFactors = FALSE)
}

# data.frame display helper: P columns with 3 significant digits, other numbers rounded
.bc_fmt <- function(d, digits = 4) {
  if (is.null(d)) return(d)
  for (k in names(d)) if (is.numeric(d[[k]])) d[[k]] <- if (grepl("^P", k)) signif(d[[k]], 3) else round(d[[k]], digits)
  d
}

# ==============================================================================
# UI
# ==============================================================================
barcodeUI <- function(id) {
  ns <- shiny::NS(id)
  shiny::div(class = "container-fluid",
    # lets the Reset buttons put a colour back into a colour box
    shiny::tags$script(shiny::HTML("(function w(){ if (!(window.Shiny && Shiny.addCustomMessageHandler)) return setTimeout(w, 30);
      Shiny.addCustomMessageHandler('bc_setcolor', function(m){ var e = document.getElementById(m.id);
        if (e) { e.value = m.value; Shiny.setInputValue(m.id, m.value); } }); })();")),
    shiny::fluidRow(
      shiny::column(4,
        .card("1. Alignment and groups",
          shiny::selectInput(ns("src"), "Alignment from", c("An aligned FASTA file (open it here)" = "file")),
          shiny::conditionalPanel(sprintf("input['%s'] == 'file'", ns("src")),
            shiny::fileInput(ns("fa"), NULL, buttonLabel = "Open aligned FASTA...", placeholder = "",
                             accept = c(".fa", ".fasta", ".fas", ".fna", ".aln", ".txt"))),
          shiny::uiOutput(ns("src_info")),
          shiny::selectInput(ns("gsrc"), "Groups (species) from", c("The sequence names" = "names", "A group table (CSV / TSV)" = "file")),
          shiny::conditionalPanel(sprintf("input['%s'] == 'names'", ns("gsrc")),
            shiny::selectInput(ns("rule"), "Species = ",
              c("Genus_species (the name after the accession)" = "species", "Genus only" = "genus", "First word of the name" = "first"))),
          shiny::conditionalPanel(sprintf("input['%s'] == 'file'", ns("gsrc")),
            shiny::fileInput(ns("gfile"), NULL, buttonLabel = "Open group table...", placeholder = "",
                             accept = c(".csv", ".tsv", ".txt")),
            shiny::helpText("A header row, a column with the sequence names and a column with the group (Species, Group, Population, Cluster ...). ",
                            "Optional columns: Genus, Role. Names are matched ignoring spaces and symbols.")),
          shiny::conditionalPanel(sprintf("input['%s'].indexOf('t:') === 0", ns("gsrc")),
            shiny::numericInput(ns("gthr"), "Same clade when the distance along the tree is at most", 0.02, min = 0, step = 0.005)),
          shiny::actionButton(ns("detect"), "Apply groups", icon = shiny::icon("wand-magic-sparkles")),
          shiny::hr(),
          shiny::selectizeInput(ns("outgroup"), "Outgroup: choose sequences, species or genera", choices = NULL, multiple = TRUE,
            options = list(placeholder = "Choose from the list (or type)...", create = TRUE, plugins = list("remove_button"))),
          shiny::div(style = "display:flex;gap:6px;flex-wrap:wrap;margin-bottom:6px",
            shiny::actionButton(ns("apply_out"), "Set as outgroup", icon = shiny::icon("arrow-right-from-bracket")),
            shiny::actionButton(ns("clear_out"), "Clear outgroup")),
          shiny::helpText("Pick as many sequences as you like (a species or a genus marks all its sequences). ",
                          "Edit Species, Genus or Role (Ingroup / Outgroup) by double-clicking a cell. ",
                          "Sequences with the same Species are one species; the same Genus makes them congeneric."),
          DT::DTOutput(ns("tbl"))),
        .card("2. Distance and tests",
          shiny::selectInput(ns("model"), "Model",
            c("Kimura 2-parameter (K2P, standard in barcoding)" = "K80", "p-distance (raw)" = "raw", "Jukes-Cantor (JC69)" = "JC69", "Tamura-Nei (TN93)" = "TN93")),
          shiny::helpText("Pairwise deletion of gaps and missing data."),
          shiny::checkboxInput(ns("gap_out"), "Include the outgroup in the barcode gap and species statistics", FALSE),
          shiny::fluidRow(
            shiny::column(6, shiny::numericInput(ns("nperm"), "Permutations", 999, min = 99, max = 9999, step = 100)),
            shiny::column(6, shiny::numericInput(ns("seed"), "Random seed", 1, min = 1, step = 1))),
          shiny::actionButton(ns("run"), "Calculate", icon = shiny::icon("play"), class = "btn-primary")),
        .card("3. Display",
          shiny::radioButtons(ns("dens_mode"), "Distribution plot", c("Density curves" = "dens", "Histogram" = "hist", "Both" = "both"), inline = TRUE),
          shiny::sliderInput(ns("dens_adj"), "Density smoothing", 0.3, 3, 1, step = 0.1),
          shiny::checkboxInput(ns("dens_thr"), "Mark the suggested threshold", TRUE),
          shiny::radioButtons(ns("box_style"), "Box plot", c("Box" = "box", "Violin + box" = "violin"), inline = TRUE),
          shiny::checkboxInput(ns("box_pts"), "Show the pairs as points", TRUE),
          shiny::selectInput(ns("pc_col"), "PCoA colour", c("Species" = "species", "Genus" = "genus", "Ingroup / outgroup" = "role")),
          shiny::selectInput(ns("pc_axes"), "PCoA axes", c("1 and 2" = "1,2", "1 and 3" = "1,3", "2 and 3" = "2,3")),
          shiny::checkboxInput(ns("pc_ell"), "95% ellipses (groups with 4+ sequences)", TRUE),
          shiny::checkboxInput(ns("pc_lab"), "Sequence labels", FALSE),
          shiny::radioButtons(ns("heat_level"), "Distance matrix", c("Sequences" = "seq", "Species means" = "sp"), inline = TRUE),
          shiny::radioButtons(ns("heat_order"), "Order of the sequences", c("By group" = "group", "By clustering" = "cluster"), inline = TRUE),
          shiny::selectInput(ns("tree_col"), "Tree colour", c("Species" = "species", "Genus" = "genus")),
          shiny::selectInput(ns("tree_layout"), "Tree layout", c("Phylogram" = "phylogram", "Cladogram" = "cladogram", "Unrooted" = "unrooted", "Fan" = "fan", "Radial" = "radial")),
          shiny::selectInput(ns("tree_root"), "Tree root", c("Outgroup (if any)" = "outgroup", "Midpoint" = "midpoint", "None" = "none")),
          shiny::selectInput(ns("tree_lab"), "Tip label", c("Sequence name" = "sequence", "Species" = "species")),
          shiny::numericInput(ns("thr_user"), "Threshold for best close match (empty = suggested)", NA, min = 0, step = 0.005))),
      shiny::column(8,
        shiny::uiOutput(ns("msg")),
        shiny::uiOutput(ns("boxes")),
        shiny::tabsetPanel(id = ns("tabs"),
          shiny::tabPanel("Distribution", value = "dens", shiny::br(), .bc_figbar(ns, "dens"), shiny::plotOutput(ns("p_dens"), height = "400px"),
            shiny::helpText("Distances between sequences of the same species (intraspecific), of different species of the same genus (congeneric), ",
                            "and of different genera (intergeneric; an outgroup of another genus is counted here when it is included). ",
                            "A barcode gap exists when the intraspecific curve ends before the other curves begin.")),
          shiny::tabPanel("Box plot", value = "box", shiny::br(), .bc_figbar(ns, "box"), shiny::plotOutput(ns("p_box"), height = "430px")),
          shiny::tabPanel("Statistics", value = "tests", shiny::br(),
            shiny::h5("Distances per category"), DT::DTOutput(ns("tbl_cat")),
            shiny::h5("Kruskal-Wallis and pairwise Wilcoxon (Holm-adjusted P)"), DT::DTOutput(ns("tbl_kw")), DT::DTOutput(ns("tbl_pw")),
            shiny::h5("Intraspecific vs interspecific"), DT::DTOutput(ns("tbl_ivi")),
            shiny::helpText("Cliff's delta: -1 to 1 (positive: the first group is larger). AUC: probability that a random intraspecific distance is smaller than a random interspecific one (1 = perfect separation)."),
            shiny::h5("PERMANOVA and ANOSIM (sequence level)"), DT::DTOutput(ns("tbl_perm")),
            shiny::helpText("The pair distances share sequences, so they are not independent and the P values of the first tests are optimistic: use them to describe, ",
                            "and rely on PERMANOVA / ANOSIM (permutation of whole sequences) for inference. R2 = share of the variation explained by the grouping; ",
                            "ANOSIM R near 1 means the groups are well separated. Needs a complete matrix: sequences with missing distances are left out of these two.")),
          shiny::tabPanel("PCoA", value = "pcoa", shiny::br(), .bc_figbar(ns, "pcoa"), shiny::plotOutput(ns("p_pcoa"), height = "520px"), shiny::uiOutput(ns("pcoa_note"))),
          shiny::tabPanel("Distance matrix", value = "mat", shiny::br(), .bc_figbar(ns, "mat"), shiny::plotOutput(ns("p_mat"), height = "560px"), shiny::br(), DT::DTOutput(ns("tbl_mat"))),
          shiny::tabPanel("Per species", value = "sp", shiny::br(), DT::DTOutput(ns("tbl_species")),
            shiny::helpText("Max_intra = largest distance inside the species. Min_inter = distance to the nearest sequence of another species. ",
                            "Min_congeneric = the same, but only within the genus. Gap = Min_inter - Max_intra (positive: a barcode gap).")),
          shiny::tabPanel("Max intra vs nearest inter", value = "scat", shiny::br(), .bc_figbar(ns, "scat"), shiny::plotOutput(ns("p_scat"), height = "460px"),
            shiny::helpText("Each point is a species. Points above the diagonal have a barcode gap; points on or below it overlap with another species. ",
                            "Species with one sequence have no intraspecific distance and are not drawn.")),
          shiny::tabPanel("Threshold", value = "thr", shiny::br(), .bc_figbar(ns, "thr"), shiny::plotOutput(ns("p_thr"), height = "400px"), shiny::uiOutput(ns("thr_note"))),
          shiny::tabPanel("Identification", value = "ident", shiny::br(), DT::DTOutput(ns("tbl_ident_sum")), shiny::br(), DT::DTOutput(ns("tbl_ident")),
            shiny::helpText("Each ingroup sequence is identified with the others as the reference. Nearest neighbour: the closest sequence has the same species. ",
                            "Best close match: the same, but 'No match' when the closest sequence is farther than the threshold. ",
                            "'Ambiguous' = the closest sequences belong to several species. Species with one sequence cannot be tested.")),
          shiny::tabPanel("Tree", value = "tree", shiny::br(), .bc_figbar(ns, "tree"), shiny::plotOutput(ns("p_tree"), height = "auto"), shiny::br(), DT::DTOutput(ns("tbl_mono")),
            shiny::helpText("Neighbour-joining tree of the distance matrix, rooted on the outgroup when it is monophyletic. Species that are not monophyletic cannot be told apart by this marker.")),
          shiny::tabPanel("Pairs", value = "pairs", shiny::br(), DT::DTOutput(ns("tbl_pairs")))),
        shiny::br(),
        shiny::div(style = "display:flex;gap:6px;flex-wrap:wrap;align-items:flex-end",
          shiny::downloadButton(ns("dl_xlsx"), "Excel (all tables)"),
          shiny::downloadButton(ns("dl_species"), "Species (CSV)"),
          shiny::downloadButton(ns("dl_pairs"), "Pairs (CSV)"),
          shiny::div(style = "width:90px", shiny::selectInput(ns("fig_fmt"), "Figure", c("png", "pdf", "svg", "jpg"), width = "100%")),
          shiny::div(style = "width:70px", shiny::numericInput(ns("fig_w"), "Width", 8, min = 2, max = 30, width = "100%")),
          shiny::div(style = "width:70px", shiny::numericInput(ns("fig_h"), "Height", 6, min = 2, max = 30, width = "100%")),
          shiny::downloadButton(ns("dl_fig"), "Download the figure of the open tab")))))
}

# ==============================================================================
# Server
# ==============================================================================
# aln_trim / aln_all : alignments of the MSA tab (named character vectors).
# aln_extra          : optional reactive, named list of more alignments (e.g. files of "Open and view files"): list("name" = alignment).
# group_sources      : optional reactive, named list of group tables, e.g. list("Population structure" = data.frame(Sequence, Group)).
# tree_sources       : optional reactive, named list of ape trees, e.g. list("Phylogeny: NJ" = nj_tree); clades give the groups.
barcodeServer <- function(id, aln_trim = shiny::reactive(NULL), aln_all = shiny::reactive(NULL),
                          aln_extra = shiny::reactive(NULL), group_sources = shiny::reactive(NULL),
                          tree_sources = shiny::reactive(NULL)) {
  shiny::moduleServer(id, function(input, output, session) {
    rv <- shiny::reactiveValues(groups = NULL, gver = 0L, res = NULL, msg = NULL, file_aln = NULL)
    set_msg <- function(txt, type = "info") rv$msg <- if (is.null(txt)) NULL else list(txt = txt, type = type)
    dtab <- function(d, dom = "ftp", n = 15) DT::datatable(.bc_fmt(d), rownames = FALSE, options = list(pageLength = n, scrollX = TRUE, dom = dom))

    # ---- where the alignment comes from: a drop-down that lists whatever is available --------------------
    ok_aln <- function(a) !is.null(a) && length(a) >= 3
    shiny::observe({
      ch <- character()
      if (ok_aln(aln_trim())) ch <- c(ch, "MSA tab: kept columns (trimmed)" = "trim")
      if (ok_aln(aln_all())) ch <- c(ch, "MSA tab: whole alignment" = "all")
      ex <- aln_extra()
      for (n in names(ex)) if (ok_aln(ex[[n]])) ch <- c(ch, stats::setNames(paste0("x:", n), n))
      ch <- c(ch, "An aligned FASTA file (open it here)" = "file")
      cur <- shiny::isolate(input$src)
      sel <- if (!is.null(cur) && cur %in% ch && !(identical(cur, "file") && is.null(shiny::isolate(rv$file_aln)) && length(ch) > 1L)) cur else ch[[1]]
      shiny::updateSelectInput(session, "src", choices = ch, selected = sel)
    })
    shiny::observeEvent(input$fa, {
      a <- tryCatch(.bc_read_fasta(input$fa$datapath), error = function(e) NULL)
      if (!ok_aln(a)) { set_msg("The file needs at least 3 sequences in FASTA format.", "warning"); rv$file_aln <- NULL; return() }
      if (length(unique(nchar(a))) != 1L) {
        set_msg("The sequences in this file do not all have the same length, so it is not an alignment. Align it first in the MSA tab.", "danger")
        rv$file_aln <- NULL; return()
      }
      rv$file_aln <- a; set_msg(NULL)
    })
    src <- shiny::reactive({
      s <- .bc_nz(input$src, "file")
      a <- if (identical(s, "trim")) aln_trim()
           else if (identical(s, "all")) aln_all()
           else if (identical(s, "file")) rv$file_aln
           else if (startsWith(s, "x:")) aln_extra()[[sub("^x:", "", s)]]
           else NULL
      if (ok_aln(a)) a else NULL
    })
    output$src_info <- shiny::renderUI({
      a <- src()
      if (is.null(a)) shiny::helpText("No alignment yet (at least 3 sequences). Align in the MSA tab, or open an aligned FASTA file.")
      else shiny::helpText(sprintf("%d sequences, %d columns.", length(a), nchar(a[[1]])))
    })

    # ---- groups -------------------------------------------------------------------------------------------
    detect <- function(keep = FALSE) {
      a <- src(); if (is.null(a)) return(invisible())
      nm <- names(a)
      sp <- .bc_group(nm, .bc_nz(input$rule, "species"))
      ge <- .bc_group(nm, "genus")
      ro <- rep("Ingroup", length(nm))
      old <- rv$groups
      if (keep && !is.null(old)) {                       # keep what the user already edited
        m <- match(nm, old$Sequence); h <- !is.na(m)
        sp[h] <- old$Species[m[h]]; ge[h] <- old$Genus[m[h]]; ro[h] <- old$Role[m[h]]
      }
      rv$groups <- data.frame(Sequence = nm, Species = sp, Genus = ge, Role = ro, stringsAsFactors = FALSE)
      rv$gver <- rv$gver + 1L
    }
    # drop-down of group sources: names, a table, and whatever the app passes (population structure, phylogeny)
    shiny::observe({
      ch <- c("The sequence names" = "names", "A group table (CSV / TSV)" = "file")
      gs <- group_sources()
      for (n in names(gs)) if (!is.null(gs[[n]])) ch <- c(ch, stats::setNames(paste0("g:", n), n))
      ts <- tree_sources()
      for (n in names(ts)) if (!is.null(ts[[n]])) ch <- c(ch, stats::setNames(paste0("t:", n), paste0(n, " (clades)")))
      cur <- shiny::isolate(input$gsrc)
      shiny::updateSelectInput(session, "gsrc", choices = ch, selected = if (!is.null(cur) && cur %in% ch) cur else "names")
    })
    apply_groups <- function() {
      if (is.null(src())) { set_msg("Choose or open an alignment first.", "warning"); return(invisible()) }
      mode <- .bc_nz(input$gsrc, "names")
      if (identical(mode, "names")) { detect(FALSE); set_msg(NULL); return(invisible()) }
      gt <- NULL; label <- ""
      if (identical(mode, "file")) {
        f <- input$gfile
        if (is.null(f)) { set_msg("Open a group table first.", "warning"); return(invisible()) }
        gt <- .bc_as_groups(.bc_read_table(f$datapath)); label <- "the group table"
      } else if (startsWith(mode, "g:")) {
        label <- sub("^g:", "", mode); gt <- .bc_as_groups(group_sources()[[label]])
      } else if (startsWith(mode, "t:")) {
        nm <- sub("^t:", "", mode); thr <- .bc_nz(input$gthr, 0.02)
        gt <- .bc_tree_groups(tree_sources()[[nm]], thr); label <- sprintf("%s, clades at distance <= %s", nm, thr)
      }
      if (is.null(gt) || !nrow(gt)) {
        set_msg("No usable groups there: a table needs a column with the sequence names and a column with the group.", "warning")
        return(invisible())
      }
      detect(FALSE)                                      # start from the names, then overwrite what matches
      g <- rv$groups
      m <- match(.bc_key(g$Sequence), .bc_key(gt$Sequence))
      idx <- which(!is.na(m))
      val <- trimws(gt$Species[m[idx]]); good <- !is.na(val) & nzchar(val)
      idx <- idx[good]; val <- val[good]
      if (!length(idx)) { set_msg(sprintf("None of the sequence names of the alignment was found in %s.", label), "warning"); return(invisible()) }
      if (!any(grepl("[A-Za-z]", val))) val <- paste0("Group_", val)       # groups numbered 1, 2, 3 ...
      g$Species[idx] <- val
      if ("Genus" %in% names(gt)) { x <- gt$Genus[m[idx]]; g$Genus[idx] <- ifelse(is.na(x) | !nzchar(x), g$Genus[idx], x) }
      if ("Role" %in% names(gt)) g$Role[idx] <- gt$Role[m[idx]]
      rv$groups <- g; rv$gver <- rv$gver + 1L
      set_msg(sprintf("Groups taken from %s: %d of %d sequences found, %d groups.%s", label, length(idx), nrow(g), length(unique(val)),
                      if (length(idx) < nrow(g)) " The other sequences keep the group read from their names (edit them in the table)." else ""),
              if (length(idx) < nrow(g)) "warning" else "info")
    }
    shiny::observeEvent(input$detect, apply_groups())
    shiny::observeEvent(src(), {
      if (is.null(src())) { rv$groups <- NULL; rv$res <- NULL; rv$gver <- rv$gver + 1L }
      else detect(TRUE)
    }, ignoreNULL = FALSE)
    # the outgroup list: every sequence, then the species and genera of the group table
    shiny::observe({
      a <- src(); rv$gver; g <- rv$groups
      ch <- if (is.null(a)) character() else {
        l <- list(Sequences = unname(names(a)))
        if (!is.null(g)) { l$Species <- sort(unique(g$Species)); l$Genera <- sort(unique(g$Genus)) }
        l
      }
      shiny::updateSelectizeInput(session, "outgroup", choices = ch, selected = shiny::isolate(input$outgroup), server = FALSE)
    })
    shiny::observeEvent(input$apply_out, {
      g <- rv$groups; shiny::req(g)
      g <- .bc_set_outgroup(g, input$outgroup)
      rv$groups <- g; rv$gver <- rv$gver + 1L
      set_msg(sprintf("%d sequence(s) are now the outgroup.", sum(g$Role == "Outgroup")), "info")
    })
    shiny::observeEvent(input$clear_out, {
      g <- rv$groups; shiny::req(g)
      g$Role <- "Ingroup"; rv$groups <- g; rv$gver <- rv$gver + 1L
      shiny::updateSelectizeInput(session, "outgroup", selected = character(0))
      set_msg(NULL)
    })

    output$tbl <- DT::renderDT({
      rv$gver
      g <- shiny::isolate(rv$groups); shiny::req(g)
      DT::datatable(g, rownames = FALSE, selection = "none",
                    editable = list(target = "cell", disable = list(columns = 0)),
                    options = list(pageLength = 10, scrollX = TRUE, dom = "ftp"))
    })
    shiny::observeEvent(input$tbl_cell_edit, {
      info <- input$tbl_cell_edit
      g <- DT::editData(rv$groups, info, rownames = FALSE)
      if (identical(as.integer(info$col), 3L)) { g$Role <- .bc_role(g$Role); rv$gver <- rv$gver + 1L }
      rv$groups <- g
    })

    # ---- calculation --------------------------------------------------------------------------------------
    shiny::observeEvent(input$run, {
      a <- src(); g <- rv$groups
      if (is.null(a) || is.null(g)) { set_msg("Provide an alignment with at least 3 sequences first.", "warning"); return() }
      nm0 <- names(a); m <- match(nm0, g$Sequence)
      sp <- trimws(g$Species[m]); ge <- trimws(g$Genus[m]); ro <- .bc_role(g$Role[m])
      sp[is.na(sp) | !nzchar(sp)] <- "(unassigned)"; ge[is.na(ge) | !nzchar(ge)] <- "(unassigned)"
      gap_out <- isTRUE(input$gap_out)
      use <- if (gap_out) rep(TRUE, length(sp)) else ro != "Outgroup"
      if (length(unique(sp[use])) < 2) {
        set_msg("At least two species are needed in the ingroup (or tick 'Include the outgroup').", "warning"); return()
      }
      a2 <- a; names(a2) <- make.unique(nm0)
      model <- .bc_nz(input$model, "K80")
      r <- shiny::withProgress(message = "Calculating distances and tests ...", value = 0.5, {
        tryCatch({
          D <- .bc_dist(a2, model)
          .bc_run(D, sp, ge, ro, model, gap_out, as.integer(.bc_nz(input$nperm, 999)), as.integer(.bc_nz(input$seed, 1)))
        }, error = function(e) { set_msg(paste("Analysis failed:", conditionMessage(e)), "danger"); NULL })
      })
      if (is.null(r)) return()
      rv$res <- r
      notes <- character()
      if (length(r$dropped))
        notes <- c(notes, sprintf("%d sequence(s) with missing distances were left out of the PCoA, PERMANOVA/ANOSIM and the tree: %s.",
                                  length(r$dropped), paste(utils::head(r$dropped, 5), collapse = ", ")))
      bad <- unique(r$sp[!grepl("^[A-Za-z][A-Za-z.-]*_[A-Za-z]", r$sp)])
      if (length(bad))
        notes <- c(notes, sprintf("These species names do not look like Genus_species: %s. Edit them in the table (a sequence without a species name counts as its own species).",
                                  paste(utils::head(bad, 4), collapse = ", ")))
      z <- r$pairs[r$pairs$Use & !is.na(r$pairs$Distance) & r$pairs$Distance == 0 & r$pairs$Species_1 != r$pairs$Species_2, , drop = FALSE]
      if (nrow(z))
        notes <- c(notes, sprintf("%d pair(s) of different species have distance 0 (identical sequences), which removes the barcode gap, for example %s vs %s (%s / %s). Check the labels.",
                                  nrow(z), z$Sequence_1[1], z$Sequence_2[1], z$Species_1[1], z$Species_2[1]))
      set_msg(if (length(notes)) paste(notes, collapse = " ") else NULL, "warning")
    })

    output$msg <- shiny::renderUI({
      if (is.null(rv$msg)) return(NULL)
      shiny::div(class = paste0("alert alert-", rv$msg$type), rv$msg$txt)
    })
    res <- shiny::reactive({ r <- rv$res; shiny::req(r); r })
    f3 <- function(x) if (is.null(x) || is.na(x)) "-" else formatC(round(x, 4) + 0, format = "f", digits = 4)

    output$boxes <- shiny::renderUI({
      r <- res(); s <- r$summary
      verdict <- if (is.na(s$gap)) "-" else if (s$gap > 0) "Yes" else "No"
      cm <- s$cat_mean
      shiny::tagList(
        shiny::fluidRow(
          shiny::column(3, .infobox("layer-group", "Species | genera | outgroup", sprintf("%d | %d | %d seq.", s$n_species, s$n_genera, s$n_out), "#2c7be5")),
          shiny::column(3, .infobox("arrows-left-right", "Mean intra | congen. | intergen.",
                                    paste(f3(cm["Intraspecific"]), "|", f3(cm["Congeneric"]), "|", f3(cm["Intergeneric"])), "#6f42c1")),
          shiny::column(3, .infobox("ruler", "Max intra | Min inter", paste(f3(s$max_intra), "|", f3(s$min_inter)), "#e67e22")),
          shiny::column(3, .infobox("circle-check", "Barcode gap overall",
                                    if (identical(verdict, "-")) "-" else sprintf("%s (%d/%d spp)", verdict, s$n_with_gap, s$n_testable), if (identical(verdict, "Yes")) "#27ae60" else "#c0392b"))),
        shiny::helpText(sprintf("Model %s. Species with a gap: %d of %d that can be tested.%s", s$model, s$n_with_gap, s$n_testable,
                                if (s$n_out == 0) "" else if (isTRUE(s$gap_out)) " The outgroup is included." else " The outgroup is not part of the gap.")))
    })

    # ---- figures ------------------------------------------------------------------------------------------
    thr_eff <- function(r) {
      u <- input$thr_user
      if (!is.null(u) && length(u) == 1L && !is.na(u)) u else if (!is.null(r$thr)) r$thr$best else NA_real_
    }
    # colours and font size chosen above each figure (the defaults apply until the user touches a control)
    fig_cols <- function(fig) {
      d <- .bc_figcols[[fig]]
      stats::setNames(vapply(d, function(cc) .bc_nz(input[[paste0("c_", fig, "_", cc$key)]], cc$default), ""), vapply(d, function(cc) cc$key, ""))
    }
    fig_fs <- function(fig, d = 13) { v <- input[[paste0("fs_", fig)]]; if (is.null(v) || length(v) != 1L || is.na(v)) d else v }
    mat_ms <- function() { cc <- fig_cols("mat"); list(scale = .bc_nz(input$mat_scale, "viridis"), low = cc[["low"]], high = cc[["high"]]) }
    plot_dens <- shiny::reactive({
      r <- res()
      mode <- .bc_nz(input$dens_mode, "dens")
      tab <- table(r$pairs$Category[r$pairs$Use & !is.na(r$pairs$Distance)])
      shiny::validate(shiny::need(sum(tab >= (if (mode == "hist") 1 else 2)) > 0, "Not enough distances to draw this plot."))
      .bc_plot_dens(r$pairs, mode, .bc_nz(input$dens_adj, 1), r$model,
                    if (isTRUE(input$dens_thr) && !is.null(r$thr)) r$thr$best else NULL, fig_cols("dens"), fig_fs("dens"))
    })
    plot_box <- shiny::reactive({
      r <- res()
      .bc_plot_box(r$pairs, if (!is.null(r$tests$kw)) r$tests$kw$P else NULL, identical(input$box_style, "violin"),
                   isTRUE(input$box_pts), r$model, fig_cols("box"), fig_fs("box"))
    })
    plot_pcoa <- shiny::reactive({
      r <- res()
      shiny::validate(shiny::need(!is.null(r$pcoa), "The PCoA needs at least 3 sequences with a complete distance matrix."))
      ax <- as.integer(strsplit(.bc_nz(input$pc_axes, "1,2"), ",")[[1]])
      k <- r$keep
      .bc_plot_pcoa(r$pcoa, r$sp[k], r$ge[k], r$role[k], .bc_nz(input$pc_col, "species"), ax, isTRUE(input$pc_ell), isTRUE(input$pc_lab),
                      .bc_nz(input$pal_pcoa, "Dark 3"), fig_fs("pcoa"))
    })
    plot_mat <- shiny::reactive({
      r <- res()
      if (identical(input$heat_level, "sp")) {
        M <- .bc_spmatrix(r$D, r$sp, r$role, include_out = TRUE)
        shiny::validate(shiny::need(nrow(M) >= 2, "At least two species are needed."))
        .bc_plot_spheat(M, r$model, mat_ms(), fig_fs("mat", 12))
      } else {
        k <- r$keep
        shiny::validate(shiny::need(length(k) >= 2, "Not enough sequences with complete distances."))
        .bc_plot_heat(r$Dc, r$sp[k], r$role[k], TRUE, .bc_nz(input$heat_order, "group"), r$model, mat_ms(), fig_fs("mat", 12))
      }
    })
    plot_scat <- shiny::reactive({
      r <- res(); d <- r$species
      shiny::validate(shiny::need(any(!is.na(d$Max_intra) & !is.na(d$Min_inter)),
                                  "No species has two or more sequences together with another species to compare with."))
      v <- fig_cols("scat")
      .bc_plot_scatter(d, c("barcode gap" = v[["gap"]], "overlap" = v[["overlap"]]), fig_fs("scat"))
    })
    plot_thr <- shiny::reactive({
      r <- res(); shiny::validate(shiny::need(!is.null(r$thr), "Both intraspecific and interspecific pairs are needed."))
      u <- input$thr_user
      .bc_plot_thr(r$thr, if (!is.null(u) && length(u) == 1L) u else NA_real_, r$model, unname(fig_cols("thr")), fig_fs("thr"))
    })
    tree <- shiny::reactive({ r <- res(); .bc_tree(r, .bc_nz(input$tree_root, "outgroup")) })
    tree_layout <- function() .bc_nz(input$tree_layout, "phylogram")
    draw_tree <- function() .bc_draw_tree(res(), tree(), .bc_nz(input$tree_col, "species"), tree_layout(),
                                          .bc_nz(input$tree_cex, 0.9), .bc_nz(input$tree_lab, "sequence"),
                                          pal = .bc_nz(input$pal_tree, "Dark 3"))
    tree_h <- function() {
      tr <- tryCatch(tree(), error = function(e) NULL)
      if (tree_layout() %in% c("phylogram", "cladogram")) .bc_tree_h(tr, .bc_nz(input$tree_cex, 0.9)) else 640
    }

    output$p_dens <- shiny::renderPlot(plot_dens())
    output$p_box <- shiny::renderPlot(plot_box())
    output$p_pcoa <- shiny::renderPlot(plot_pcoa())
    output$p_mat <- shiny::renderPlot(plot_mat())
    output$p_scat <- shiny::renderPlot(plot_scat())
    output$p_thr <- shiny::renderPlot(plot_thr())
    output$p_tree <- shiny::renderPlot(draw_tree(), height = function() tree_h())

    output$pcoa_note <- shiny::renderUI({
      r <- res(); shiny::req(r$pcoa)
      shiny::helpText(sprintf("PCoA of the %s distance matrix (%d sequences). The percentages are the share of the positive eigenvalues. ",
                              r$model, length(r$keep)),
                      if (r$pcoa$neg_pct > 5) sprintf("Negative eigenvalues make up %.1f%% of the total: the matrix is not fully Euclidean, so read the axes with some care.", r$pcoa$neg_pct) else "")
    })
    output$thr_note <- shiny::renderUI({
      r <- res(); shiny::req(r$thr)
      shiny::helpText(sprintf("Suggested threshold (smallest total error): %s. At this value %.1f%% of the conspecific pairs are split and %.1f%% of the heterospecific pairs are lumped. ",
                              f3(r$thr$best), 100 * r$thr$fn, 100 * r$thr$fp),
                      "The dotted green line is the threshold you typed for the best close match.")
    })

    # ---- tables -------------------------------------------------------------------------------------------
    output$tbl_cat <- DT::renderDT({ r <- res(); shiny::req(r$catstats); dtab(r$catstats, "t") })
    output$tbl_kw <- DT::renderDT({ r <- res(); shiny::req(r$tests$kw); dtab(r$tests$kw, "t") })
    output$tbl_pw <- DT::renderDT({ r <- res(); shiny::req(r$tests$pw); dtab(r$tests$pw, "t") })
    output$tbl_ivi <- DT::renderDT({ r <- res(); shiny::req(r$tests$ivi); dtab(r$tests$ivi, "t") })
    output$tbl_perm <- DT::renderDT({ r <- res(); shiny::req(r$perm); dtab(r$perm, "t") })
    output$tbl_species <- DT::renderDT({ r <- res(); dtab(r$species) })
    output$tbl_pairs <- DT::renderDT({ r <- res(); dtab(r$pairs[, setdiff(names(r$pairs), c("Genus_1", "Genus_2", "Role_1", "Role_2", "Use"))]) })
    output$tbl_mat <- DT::renderDT({
      r <- res()
      M <- if (identical(input$heat_level, "sp")) .bc_spmatrix(r$D, r$sp, r$role, TRUE) else r$D
      dtab(data.frame(Name = rownames(M), M, check.names = FALSE, stringsAsFactors = FALSE), n = 10)
    })
    ident <- shiny::reactive({ r <- res(); .bc_ident(r$D, r$sp, r$role, thr_eff(r)) })
    output$tbl_ident_sum <- DT::renderDT({ id <- ident(); shiny::req(id); dtab(.bc_ident_summary(id), "t") })
    output$tbl_ident <- DT::renderDT({ id <- ident(); shiny::req(id); dtab(id) })
    output$tbl_mono <- DT::renderDT({
      r <- res(); tr <- tree(); shiny::req(tr)
      lab <- if (identical(input$tree_col, "genus")) r$ge[r$keep] else r$sp[r$keep]
      names(lab) <- rownames(r$Dc)
      d <- .bc_mono(tr, lab); shiny::req(d)
      dtab(d, "t")
    })

    # ---- downloads ----------------------------------------------------------------------------------------
    summ_df <- function(r) {
      s <- r$summary
      data.frame(Item = c("Model", "Sequences", "Species", "Genera", "Outgroup sequences", "Outgroup in the gap", "Intraspecific pairs", "Interspecific pairs",
                          "Mean intraspecific", "Largest intraspecific", "Mean interspecific", "Smallest interspecific",
                          "Gap (smallest inter - largest intra)", "Species with a gap", "Permutations", "Sequences left out (missing distances)"),
                 Value = c(s$model, s$n_seq, s$n_species, s$n_genera, s$n_out, if (isTRUE(s$gap_out)) "yes" else "no", s$n_intra, s$n_inter,
                           f3(s$mean_intra), f3(s$max_intra), f3(s$mean_inter), f3(s$min_inter), f3(s$gap),
                           sprintf("%d of %d", s$n_with_gap, s$n_testable), r$nperm, length(r$dropped)), stringsAsFactors = FALSE)
    }
    output$dl_species <- shiny::downloadHandler(function() "barcode_gap_species.csv",
      function(file) { shiny::req(rv$res); utils::write.csv(rv$res$species, file, row.names = FALSE) })
    output$dl_pairs <- shiny::downloadHandler(function() "barcode_gap_pairs.csv",
      function(file) { shiny::req(rv$res); utils::write.csv(rv$res$pairs, file, row.names = FALSE) })
    output$dl_xlsx <- shiny::downloadHandler(function() "barcode_gap.xlsx",
      function(file) {
        r <- shiny::req(rv$res)
        tryCatch({
          sh <- list(Summary = summ_df(r), Species = r$species, Category_stats = r$catstats,
                     Kruskal_Wallis = r$tests$kw, Pairwise_Wilcoxon = r$tests$pw, Intra_vs_inter = r$tests$ivi,
                     PERMANOVA_ANOSIM = r$perm,
                     Species_matrix = { M <- .bc_spmatrix(r$D, r$sp, r$role, TRUE); data.frame(Species = rownames(M), M, check.names = FALSE) },
                     Distance_matrix = data.frame(Sequence = rownames(r$D), r$D, check.names = FALSE),
                     Identification = ident(), Pairs = r$pairs)
          .sg_xlsx(Filter(Negate(is.null), sh), file)
        }, error = function(e) {
          shiny::showNotification(paste("Excel file not written:", conditionMessage(e)), type = "error", duration = 15)
          writeLines(conditionMessage(e), file)
        })
      })
    fig_fun <- function(tab) switch(.bc_nz(tab, ""),
      dens = function() print(plot_dens()), box = function() print(plot_box()), pcoa = function() print(plot_pcoa()),
      mat = function() print(plot_mat()), scat = function() print(plot_scat()), thr = function() print(plot_thr()),
      tree = function() draw_tree(), NULL)
    output$dl_fig <- shiny::downloadHandler(
      filename = function() paste0("barcode_", .bc_nz(input$tabs, "figure"), ".", .bc_nz(input$fig_fmt, "png")),
      content = function(file) {
        fn <- fig_fun(input$tabs)
        if (is.null(fn)) {
          shiny::showNotification("The open tab has no figure: open Distribution, Box plot, PCoA, Distance matrix, Max intra vs nearest inter, Threshold or Tree.", type = "warning", duration = 10)
          writeLines("This tab has no figure.", file); return(invisible())
        }
        w <- .bc_nz(input$fig_w, 8); h <- .bc_nz(input$fig_h, 6)
        switch(.bc_nz(input$fig_fmt, "png"),
          pdf = grDevices::pdf(file, width = w, height = h),
          svg = grDevices::svg(file, width = w, height = h),
          jpg = grDevices::jpeg(file, width = w, height = h, units = "in", res = 300, quality = 95),
          grDevices::png(file, width = w, height = h, units = "in", res = 300))
        on.exit(grDevices::dev.off(), add = TRUE)
        fn()
      })

    # ---- project state ------------------------------------------------------------------------------------
    get_state <- function() shiny::isolate(list(groups = rv$groups, res = rv$res))
    set_state <- function(s) {
      g <- s$groups
      if (!is.null(g) && !("Species" %in% names(g)) && "Group" %in% names(g)) {   # project saved by the previous version
        g <- data.frame(Sequence = g$Sequence, Species = g$Group, Genus = .bc_group(g$Sequence, "genus"),
                        Role = "Ingroup", stringsAsFactors = FALSE)
      }
      rv$groups <- g
      rv$res <- if (!is.null(s$res) && !is.null(s$res$pairs$Category) && !("Outgroup" %in% s$res$pairs$Category)) s$res else NULL   # results of the older 4-category version are dropped
      rv$gver <- rv$gver + 1L
      invisible()
    }
    # Reset button above a figure: default colours, palette / scale and font size
    for (fig in c("dens", "box", "pcoa", "mat", "scat", "thr", "tree")) local({
      f <- fig
      shiny::observeEvent(input[[paste0("reset_", f)]], {
        for (cc in .bc_figcols[[f]])
          session$sendCustomMessage("bc_setcolor", list(id = session$ns(paste0("c_", f, "_", cc$key)), value = cc$default))
        if (!identical(f, "tree")) shiny::updateSliderInput(session, paste0("fs_", f), value = if (identical(f, "mat")) 12 else 13)
        if (identical(f, "pcoa")) shiny::updateSelectInput(session, "pal_pcoa", selected = "Dark 3")
        if (identical(f, "mat")) shiny::updateSelectInput(session, "mat_scale", selected = "viridis")
        if (identical(f, "tree")) { shiny::updateSelectInput(session, "pal_tree", selected = "Dark 3"); shiny::updateSliderInput(session, "tree_cex", value = 0.9) }
      }, ignoreInit = TRUE)
    })
    list(get_state = get_state, set_state = set_state, result = shiny::reactive(rv$res))
  })
}
