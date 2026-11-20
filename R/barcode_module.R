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

# What is compared. "species": groups = species, the next level up = genus (the usual barcode-gap analysis).
# "genus": groups = genera, the next level up = family (for markers that are too conserved to tell species apart).
.bc_levels <- list(
  species = list(unit = "species", parent = "genus",  cats = c("Intraspecific", "Congeneric", "Intergeneric"),
                 ucol = "Species", pcol = "Genus"),
  genus   = list(unit = "genus",   parent = "family", cats = c("Intrageneric", "Same family", "Different family"),
                 ucol = "Genus", pcol = "Family"))
# The three categories of the pairs of one result (the names are kept on the pair table).
.bc_cats <- function(pairs = NULL) {
  a <- if (!is.null(pairs)) attr(pairs, "cats") else NULL
  if (is.null(a)) names(.bc_cat_cols) else a
}
# colours by position (1st = within a group, 2nd = same parent, 3rd = different parent)
.bc_cat_pal <- function(cats, cols = .bc_cat_cols) stats::setNames(unname(cols)[seq_along(cats)], cats)

# Reference markers. level = the level the marker is usually good for; thr = indicative cut-offs (distance, not identity);
# coding = protein-coding (codon positions can be chosen). The numbers are guides taken from the literature, not rules.
.bc_markers <- list(
  coi = list(label = "COI (animals, DNA barcode)", level = "species", coding = TRUE,
    note = paste("COI separates most animal species: the barcode gap is the usual test. A common rule of thumb is that the distance between species is about ten times",
                 "the distance inside a species (about 2-3% in K2P). Use the 3rd codon position to see the most variable sites."),
    thr = c("about 3% (common cut-off)" = 0.03)),
  cytb = list(label = "Cytochrome b (vertebrates)", level = "species", coding = TRUE,
    note = "Cytochrome b is variable enough for species and for populations of vertebrates. Species level; the 3rd codon position is the most informative.", thr = numeric()),
  mt16s = list(label = "16S rRNA, mitochondrial (animals)", level = "species", coding = FALSE,
    note = "Mitochondrial 16S is more conserved than COI and has gaps: expect weaker species separation than COI. Pairwise deletion is the safer choice.", thr = numeric()),
  s12 = list(label = "12S rRNA (fish, vertebrates)", level = "species", coding = FALSE,
    note = "12S is short and conserved: several close species can be identical. Check the 'Identification' tab, and use the genus level for deeper comparisons.", thr = numeric()),
  b16s = list(label = "16S rRNA (bacteria, archaea)", level = "genus", coding = FALSE,
    note = paste("16S is highly conserved in prokaryotes: it is reliable for genus (and family) and often too conserved for species, where closely related species can share almost identical sequences.",
                 "Indicative identity cut-offs: about 98.65% for species, 94.5% for genus and 86.5% for family (Kim et al. 2014; Yarza et al. 2014), shown here as distances."),
    thr = c("species 98.65% identity" = 0.0135, "genus 94.5%" = 0.055, "family 86.5%" = 0.135)),
  s18 = list(label = "18S rRNA (eukaryotes)", level = "genus", coding = FALSE,
    note = "18S is very conserved: it separates families, orders and often genera, rarely species. Use the genus level.", thr = numeric()),
  s28 = list(label = "28S rRNA D1-D3 (fungi, animals)", level = "genus", coding = FALSE,
    note = "28S D1-D3 is more variable than 18S but still conserved: good for genus, sometimes species. Compare both levels.", thr = numeric()),
  rbcl = list(label = "rbcL (plants)", level = "genus", coding = TRUE,
    note = paste("rbcL is conserved in plants: it is good for families and genera and weak for species, because close species share the same sequence.",
                 "Test the genus level, or combine it with matK / ITS for species."), thr = numeric()),
  matk = list(label = "matK (plants)", level = "species", coding = TRUE,
    note = "matK evolves faster than rbcL and discriminates more species, but many close species still overlap. The species level is tried first; look at the genus level if there is no gap.", thr = numeric()),
  trnh = list(label = "trnH-psbA (plants, spacer)", level = "species", coding = FALSE,
    note = "A non-coding spacer: variable and good for species, but it has many indels that complicate the alignment (pairwise deletion).", thr = numeric()),
  itsf = list(label = "ITS (fungi, official barcode)", level = "species", coding = FALSE,
    note = paste("ITS is very variable. It is the fungal barcode, with a typical species cut-off near 3% (97% identity) that varies between groups.",
                 "Within one species the variation can be large (several copies per genome); alignments are reliable only between close species, so align within a genus."),
    thr = c("about 3% (97% identity)" = 0.03)),
  itsp = list(label = "ITS / ITS2 (plants, algae)", level = "species", coding = FALSE,
    note = "ITS is the most variable common marker: strong for species, but with indels, copy variation and hybrids. Align within a genus or family; it can be too variable for distant taxa.", thr = numeric()),
  other = list(label = "Another marker (no advice)", level = "species", coding = FALSE, note = "", thr = numeric()))
.bc_marker_choices <- function() stats::setNames(names(.bc_markers), vapply(.bc_markers, function(m) m$label, ""))

.bc_nz <- function(x, d) {
  if (is.null(x) || length(x) == 0L || (length(x) == 1L && (is.na(x) || identical(x, "")))) d else x
}

# ---- groups ------------------------------------------------------------------------------------------------
# Name of the group of each sequence. Leading accession numbers (MG272989_, NC_012920.1_, ...) are removed first.
# Taxon names are read the way a taxonomist writes them:
#   "Genus species"                  -> Genus_species
#   "Genus aff. species", "cf."      -> Genus_aff_species   (an uncertain identification is NOT the same species as the named one)
#   "Genus sp. X", "Genus sp. 12"    -> Genus_sp_X          ("sp." alone, or followed by isolate / strain / voucher ..., stays a species of its own)
#   "Genus species var. / subsp. x"  -> Genus_species       (the infraspecific rank is ignored)
#   "Genus x species" (hybrid)       -> its own group
# A name that cannot be read as a taxon (no letters, or the epithet is a word such as isolate or voucher) becomes a group of its own.
.bc_qual_unsure <- c("aff", "cf", "nr")
.bc_qual_sp     <- c("sp", "spp", "sp.", "spp.")
.bc_skip_words  <- c("isolate", "strain", "voucher", "clone", "specimen", "haplotype", "genotype", "mitochondrion", "mitochondrial",
                     "chloroplast", "plastid", "small", "large", "partial", "complete", "gene", "genes", "rrna", "16s", "18s", "28s",
                     "its", "internal", "ribosomal", "cytochrome", "oxidase", "maturase", "ribulose", "voucher", "from", "uncultured")
.bc_species_of <- function(t) {
  t <- t[nzchar(t)]
  if (length(t) < 2L || !grepl("^[A-Za-z]", t[1])) return("")
  g <- gsub("[^A-Za-z-]", "", t[1])
  w <- tolower(gsub("[,;()]", "", t[-1])); wd <- sub("\\.$", "", w)
  if (!length(wd)) return("")
  if (wd[1] %in% .bc_qual_unsure) {
    if (length(wd) >= 2L && grepl("^[a-z][a-z-]+$", wd[2]) && !(wd[2] %in% .bc_skip_words)) return(paste0(g, "_", wd[1], "_", wd[2]))
    return("")
  }
  if (wd[1] %in% c("sp", "spp")) {
    if (length(wd) >= 2L && !(wd[2] %in% .bc_skip_words) && grepl("^[A-Za-z0-9][A-Za-z0-9._-]*$", wd[2]))
      return(paste0(g, "_sp_", gsub("[^A-Za-z0-9]", "", t[2 + 1])))
    return("")
  }
  if (wd[1] %in% c("x", "\u00d7")) return("")
  if (grepl("^[a-z][a-z-]+$", wd[1]) && !(wd[1] %in% .bc_skip_words)) return(paste0(g, "_", wd[1]))
  ""
}
.bc_group <- function(nm, rule = "species") {
  x <- sub("^[A-Z]{1,6}_?[0-9]{3,}(\\.[0-9]+)?[_ ]+", "", nm)
  tok <- strsplit(x, "[_ ]+")
  out <- vapply(tok, function(t) {
    t <- t[nzchar(t)]
    if (!length(t) || !grepl("^[A-Za-z]", t[1])) return("")   # no taxon name (e.g. "1392-2958"): the full sequence name is used
    if (rule %in% c("genus", "first")) return(gsub("[^A-Za-z-]", "", t[1]))
    .bc_species_of(t)
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
# the key of a colour box stays "Intraspecific" / "Congeneric" / "Intergeneric"; the label follows the level that is compared
.bc_cat_cc <- function(cats = names(.bc_cat_cols))
  lapply(seq_along(.bc_cat_cols), function(i) .bc_cc(names(.bc_cat_cols)[i], cats[i], unname(.bc_cat_cols[i])))
.bc_figcols <- list(
  dens = .bc_cat_cc(), box = .bc_cat_cc(),
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
.bc_figbar <- function(ns, fig, cats = NULL) {
  lst <- if (fig %in% c("dens", "box") && !is.null(cats)) .bc_cat_cc(cats) else .bc_figcols[[fig]]
  boxes <- lapply(lst, function(cc) .bc_colorbox(ns(paste0("c_", fig, "_", cc$key)), cc$label, cc$default))
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
# Matrix of lower-case bases (rows = sequences) -> distance matrix. model: "K80" (K2P), "raw" (p-distance), "JC69", "TN93", "T92".
# gamma: shape of the gamma rate variation (0 or NA = equal rates). Gaps and ambiguity codes are always ignored pair by pair.
.bc_dist_mat <- function(m, model = "K80", gamma = 0) {
  if (!requireNamespace("ape", quietly = TRUE)) stop("The 'ape' package is required: install.packages(\"ape\")")
  g <- if (is.finite(gamma) && gamma > 0 && model %in% c("K80", "JC69", "TN93")) gamma else FALSE   # ape has no gamma for raw and T92
  D <- ape::dist.dna(ape::as.DNAbin(m), model = model, gamma = g, pairwise.deletion = TRUE, as.matrix = TRUE)
  D[!is.finite(D)] <- NA_real_
  D[!is.na(D) & D <= 0] <- 0        # K2P can give -0 / tiny negative values for identical sequences
  D
}
# The alignment as a matrix with the sites the user chose: complete deletion (a column with a gap or an ambiguity code in any
# sequence is dropped) and the codon positions (pos: any of 1, 2, 3; frame = column where the reading frame starts).
.bc_sites <- function(aln, pairwise = TRUE, pos = 1:3, frame = 1L) {
  m <- do.call(rbind, strsplit(tolower(unname(toupper(aln))), ""))
  rownames(m) <- names(aln)
  keep <- rep(TRUE, ncol(m))
  if (!pairwise) keep <- apply(m, 2, function(x) all(x %in% c("a", "c", "g", "t")))
  if (length(pos) && length(pos) < 3L) {
    ph <- ((seq_len(ncol(m)) - as.integer(frame)) %% 3L) + 1L
    keep <- keep & ph %in% as.integer(pos)
  }
  m[, keep, drop = FALSE]
}
.bc_dist <- function(aln, model = "K80", gamma = 0, pairwise = TRUE, pos = 1:3, frame = 1L) {
  m <- .bc_sites(aln, pairwise, pos, frame)
  if (ncol(m) < 10L) stop("Fewer than 10 sites are left with these settings (complete deletion or codon positions): choose pairwise deletion or more positions.")
  .bc_dist_mat(m, model, gamma)
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
# targets: names of the target groups; when given, only the pairs that involve a target are used and only targets get a row in the species table.
.bc_analyse <- function(D, sp, ge, ro, model = "K80", gap_out = FALSE, cats = names(.bc_cat_cols), targets = NULL) {
  nm <- rownames(D)
  ix <- which(lower.tri(D), arr.ind = TRUE)
  i <- ix[, 2]; j <- ix[, 1]
  pairs <- data.frame(Sequence_1 = nm[i], Sequence_2 = nm[j], Species_1 = sp[i], Species_2 = sp[j],
                      Genus_1 = ge[i], Genus_2 = ge[j], Role_1 = ro[i], Role_2 = ro[j],
                      Distance = D[ix], stringsAsFactors = FALSE)
  pairs$Category <- ifelse(pairs$Species_1 == pairs$Species_2, cats[1],
                      ifelse(pairs$Genus_1 == pairs$Genus_2, cats[2], cats[3]))
  attr(pairs, "cats") <- cats
  use <- if (gap_out) rep(TRUE, length(sp)) else ro != "Outgroup"
  pairs$Use <- use[i] & use[j]
  if (length(targets)) pairs$Use <- pairs$Use & (pairs$Species_1 %in% targets | pairs$Species_2 %in% targets)

  pu <- pairs[pairs$Use & !is.na(pairs$Distance), , drop = FALSE]
  sl <- unique(sp[use])
  if (length(targets)) sl <- intersect(sl, targets)
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
  cat_mean <- vapply(cats, function(k) {
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
  rows <- lapply(.bc_cats(pairs), function(k) {
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
  d$Category <- factor(d$Category, levels = .bc_cats(pairs))
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
    out$ivi <- data.frame(Test = "Mann-Whitney: within groups vs between groups", n_intra = n1, n_inter = n2,
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

# The grouping tests use the ingroup only (an outgroup far away would dominate them) unless with_out is TRUE;
# "Ingroup vs outgroup" always uses every sequence. labs = names of the two levels.
.bc_permtests <- function(Dc, sp, ge, ro, nperm = 999, seed = 1, with_out = FALSE, labs = c("Species", "Genus")) {
  rows <- list()
  add <- function(label, grp, idx = seq_len(nrow(Dc))) {
    Dx <- Dc[idx, idx, drop = FALSE]; grp <- factor(grp[idx])
    if (nlevels(grp) < 2L || nrow(Dx) < 4L) return(invisible())
    pm <- .bc_permanova(Dx, grp, nperm, seed); an <- .bc_anosim(Dx, grp, nperm, seed)
    if (!is.null(pm))
      rows[[length(rows) + 1L]] <<- data.frame(Test = "PERMANOVA", Grouping = label, N = nrow(Dx), Groups = nlevels(grp),
        Statistic = pm$stat, R2 = pm$R2, P = pm$P, Permutations = nperm, stringsAsFactors = FALSE)
    if (!is.null(an))
      rows[[length(rows) + 1L]] <<- data.frame(Test = "ANOSIM", Grouping = label, N = nrow(Dx), Groups = nlevels(grp),
        Statistic = an$stat, R2 = NA_real_, P = an$P, Permutations = nperm, stringsAsFactors = FALSE)
  }
  ii <- if (with_out) seq_len(nrow(Dc)) else which(ro != "Outgroup")
  add(labs[1], sp, ii); add(labs[2], ge, ii)
  if (any(ro == "Outgroup") && any(ro != "Outgroup")) add("Ingroup vs outgroup", ro)
  if (length(rows)) do.call(rbind, rows) else NULL
}

# ---- distances within and between groups (mean, net, bootstrap standard error) ---------------------------------------
# D: distance matrix; grp: group of every row. Returns list(within = data.frame, between = data.frame).
# The net distance between groups X and Y is dXY - (dX + dY) / 2 (Nei and Li 1979), with dX the mean within group X.
.bc_group_means <- function(D, grp) {
  lv <- unique(grp); k <- length(lv)
  W <- stats::setNames(rep(NA_real_, k), lv); N <- stats::setNames(as.integer(table(factor(grp, levels = lv))), lv)
  B <- matrix(NA_real_, k, k, dimnames = list(lv, lv))
  for (a in seq_len(k)) {
    ia <- which(grp == lv[a])
    if (length(ia) > 1L) { s <- D[ia, ia, drop = FALSE]; v <- s[lower.tri(s)]; v <- v[!is.na(v)]; if (length(v)) W[a] <- mean(v) }
    for (b in seq_len(k)) if (b > a) {
      ib <- which(grp == lv[b]); v <- as.vector(D[ia, ib, drop = FALSE]); v <- v[!is.na(v)]
      if (length(v)) B[a, b] <- B[b, a] <- mean(v)
    }
  }
  list(lv = lv, N = N, W = W, B = B)
}
# boot: matrix of the sites (see .bc_sites) and the settings to resample the columns nboot times (standard errors).
.bc_group_dist <- function(D, grp, m = NULL, model = "K80", gamma = 0, nboot = 0L, seed = 1L) {
  g0 <- .bc_group_means(D, grp)
  lv <- g0$lv; k <- length(lv)
  net <- function(g) { n <- g$B; for (a in seq_len(k)) for (b in seq_len(k)) if (a != b) n[a, b] <- g$B[a, b] - mean(c(g$W[a], g$W[b]), na.rm = FALSE); n }
  Net0 <- net(g0)
  seW <- rep(NA_real_, k); seB <- matrix(NA_real_, k, k); seN <- seB
  if (!is.null(m) && nboot >= 10L) {
    set.seed(seed)
    bw <- matrix(NA_real_, nboot, k); bb <- array(NA_real_, c(nboot, k, k)); bn <- bb
    for (i in seq_len(nboot)) {
      Di <- tryCatch(.bc_dist_mat(m[, sample.int(ncol(m), replace = TRUE), drop = FALSE], model, gamma), error = function(e) NULL)
      if (is.null(Di)) next
      gi <- .bc_group_means(Di, grp)
      bw[i, ] <- gi$W; bb[i, , ] <- gi$B; bn[i, , ] <- net(gi)
    }
    sd2 <- function(x) { x <- x[is.finite(x)]; if (length(x) >= 2L) stats::sd(x) else NA_real_ }
    seW <- apply(bw, 2, sd2)
    for (a in seq_len(k)) for (b in seq_len(k)) if (a != b) { seB[a, b] <- sd2(bb[, a, b]); seN[a, b] <- sd2(bn[, a, b]) }
  }
  within <- data.frame(Group = lv, N = unname(g0$N), Mean_within = unname(g0$W), SE = seW, stringsAsFactors = FALSE)
  ix <- which(upper.tri(g0$B), arr.ind = TRUE)
  between <- data.frame(Group_1 = lv[ix[, 1]], Group_2 = lv[ix[, 2]], Mean_between = g0$B[ix], SE = seB[ix],
                        Net_between = Net0[ix], SE_net = seN[ix], stringsAsFactors = FALSE)
  list(within = within, between = between, mat = g0$B, net = Net0)
}

# ---- PCoA, threshold, identification, tree ------------------------------------------------------------------
# use: logical vector (rows of Dc to include); NULL = all
.bc_pcoa <- function(Dc, use = NULL) {
  if (!is.null(use)) Dc <- Dc[use, use, drop = FALSE]
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
.bc_ident <- function(D, sp, ro, thr = NA_real_, targets = NULL) {
  idx <- which(ro != "Outgroup")
  if (length(idx) < 3L) return(NULL)
  D2 <- D[idx, idx, drop = FALSE]; s2 <- sp[idx]; nm <- rownames(D2)
  qs <- if (length(targets)) which(s2 %in% targets) else seq_along(idx)
  if (!length(qs)) return(NULL)
  rows <- lapply(qs, function(q) {
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
.bc_run <- function(D, sp, ge, ro, model = "K80", gap_out = FALSE, nperm = 999, seed = 1, cats = names(.bc_cat_cols), targets = NULL) {
  an <- .bc_analyse(D, sp, ge, ro, model, gap_out, cats, targets)
  kp <- .bc_complete(D)
  Dc <- D[kp, kp, drop = FALSE]
  c(an, list(D = D, Dc = Dc, sp = sp, ge = ge, role = ro, keep = kp,
             dropped = rownames(D)[setdiff(seq_len(nrow(D)), kp)],
             catstats = .bc_catstats(an$pairs),
             tests = .bc_tests(an$pairs),
             perm = if (length(kp) >= 4L) .bc_permtests(Dc, sp[kp], ge[kp], ro[kp], nperm, seed, gap_out,
                                                         c(.bc_levels[[.bc_lvl_of(cats)]]$ucol, .bc_levels[[.bc_lvl_of(cats)]]$pcol)) else NULL,
             pcoa = if (length(kp) >= 3L) .bc_pcoa(Dc) else NULL,
             thr = .bc_threshold(an$pairs),
             nperm = nperm, model = model, cats = cats, level = .bc_lvl_of(cats), targets = targets))
}
.bc_lvl_of <- function(cats) if (identical(cats[1], "Intrageneric")) "genus" else "species"

# ---- figures -----------------------------------------------------------------------------------------------
.bc_theme <- function(fs = 13) list(ggplot2::theme_minimal(base_size = fs), ggplot2::theme(legend.position = "top"))

.bc_plot_dens <- function(pairs, mode = "dens", adjust = 1, model = "K80", thr = NULL, cols = .bc_cat_cols, fs = 13, refs = NULL) {
  d <- pairs[pairs$Use & !is.na(pairs$Distance), , drop = FALSE]
  tab <- table(d$Category)
  d$Category <- factor(d$Category, levels = .bc_cats(pairs))
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
  if (length(refs))      # indicative cut-offs of the marker (distances), dotted
    p <- p + ggplot2::geom_vline(xintercept = unname(refs), linetype = 3, colour = "#6f42c1", linewidth = 0.6) +
      ggplot2::annotate("text", x = unname(refs), y = Inf, label = paste0(" ", names(refs)), hjust = 0, vjust = 3.4, angle = 90,
                        colour = "#6f42c1", size = fs * 0.7 / ggplot2::.pt)
  if (!is.null(thr) && is.finite(thr))
    p <- p + ggplot2::geom_vline(xintercept = thr, linetype = 2, colour = "grey25") +
      ggplot2::annotate("text", x = thr, y = Inf, label = sprintf(" threshold %.4f", thr), hjust = 0, vjust = 1.6, size = fs * 0.8 / ggplot2::.pt)
  p
}

.bc_plot_box <- function(pairs, kw_p = NULL, violin = FALSE, pts = TRUE, model = "K80", cols = .bc_cat_cols, fs = 13) {
  d <- pairs[pairs$Use & !is.na(pairs$Distance), , drop = FALSE]
  lev <- .bc_cats(pairs)[.bc_cats(pairs) %in% d$Category]
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
  # very long sequence names (GenBank titles) would use the whole width and squeeze the matrix to a strip: cut them
  lab <- if (show) { l <- nm[ord]; ifelse(nchar(l) > 42L, paste0(substr(l, 1L, 41L), "\u2026"), l) } else NULL
  p <- ggplot2::ggplot(df, ggplot2::aes(x, y, fill = Distance)) + ggplot2::geom_raster() +
    .bc_fill(ms, "grey85", model) +
    ggplot2::scale_x_continuous(breaks = if (show) seq_len(n) else NULL, labels = lab, expand = c(0, 0)) +
    ggplot2::scale_y_continuous(breaks = if (show) n + 1 - seq_len(n) else NULL, labels = lab, expand = c(0, 0)) +
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

.bc_plot_scatter <- function(sp_tab, cols = c("barcode gap" = "#27ae60", "overlap" = "#c0392b"), fs = 13, unit = "species") {
  d <- sp_tab[!is.na(sp_tab$Max_intra) & !is.na(sp_tab$Min_inter), , drop = FALSE]
  d$Result <- ifelse(d$Gap > 0, "barcode gap", "overlap")
  lim <- c(0, max(d$Max_intra, d$Min_inter) * 1.08)
  ggplot2::ggplot(d, ggplot2::aes(Max_intra, Min_inter, colour = Result, label = Species)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey40") +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_text(vjust = -0.8, size = fs * 0.75 / ggplot2::.pt, check_overlap = TRUE, show.legend = FALSE) +
    ggplot2::scale_colour_manual(values = cols, name = NULL) +
    ggplot2::coord_cartesian(xlim = lim, ylim = lim) +
    ggplot2::labs(x = paste("Largest distance inside a", unit), y = paste("Distance to the nearest other", unit)) + .bc_theme(fs)
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
  f <- pick("family", c(i, j, g)); if (!is.na(f)) out$Family <- as.character(x[[f]])
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


# ---- MEGA-like analysis menu and preferences table ----------------------------------------------------------------
.bc_analyses <- function() list(
  pairwise = list(title = "Compute pairwise distances", icon = "table-cells", desc = "the distance between every two sequences"),
  overall  = list(title = "Compute overall mean distance", icon = "equals", desc = "one mean over all the sequences"),
  within   = list(title = "Compute within group mean distance", icon = "bullseye", desc = "the mean inside each group"),
  between  = list(title = "Compute between group mean distance", icon = "arrows-left-right", desc = "the mean between every two groups"),
  net      = list(title = "Compute net between group mean distances", icon = "scale-balanced", desc = "between minus the variation inside the groups"),
  barcode  = list(title = "Barcode gap and species identification", icon = "barcode", desc = "intraspecific vs interspecific distances, tests, tree"))
.bcp_sec <- function(txt) shiny::div(class = "bcp-sec", txt)
.bcp_row <- function(label, ...) shiny::div(class = "bcp-row", shiny::div(class = "bcp-o", label), shiny::div(class = "bcp-s", ...))
.bc_style <- function() shiny::tags$style(shiny::HTML("
  .bc-an .shiny-options-group{display:flex;flex-direction:column;gap:6px}
  .bc-an .radio,.bc-an .form-check{margin:0;padding:0}
  .bc-an .radio label,.bc-an .form-check label{display:flex;align-items:flex-start;gap:10px;width:100%;margin:0;padding:8px 12px;border:1px solid #dbe2ec;
    border-radius:9px;background:#fff;cursor:pointer;font-weight:600;transition:background .12s,border-color .12s}
  .bc-an .radio label:hover,.bc-an .form-check label:hover{background:#f3f7ff;border-color:#9db8ea}
  .bc-an input[type=radio]{position:static!important;margin:4px 0 0!important;flex:none}
  .bc-an label:has(input:checked){background:#eaf1ff;border-color:#2563eb;box-shadow:0 0 0 2px rgba(37,99,235,.14)}
  .bc-an label>span{display:flex;align-items:flex-start;gap:10px}
  .bc-ai{width:20px;text-align:center;color:#2563eb;margin-top:2px}
  .bc-at{display:flex;flex-direction:column;line-height:1.25}
  .bc-ad{font-weight:400;color:#64748b;font-size:12px}
  .bcp{border:1px solid #d5dbe3;border-radius:6px;overflow:hidden;margin-bottom:10px}
  .bcp-sec{background:#f1f5f9;color:#0f8b8d;font-weight:700;text-transform:uppercase;font-size:12px;letter-spacing:.04em;padding:5px 10px;border-top:1px solid #d5dbe3}
  .bcp>.bcp-sec:first-child{border-top:0}
  .bcp-row{display:grid;grid-template-columns:42% 58%;align-items:center;border-top:1px solid #eef1f5}
  .bcp-o{padding:6px 10px;font-size:13px;color:#334155}
  .bcp-s{padding:4px 8px 4px 4px;border-left:1px solid #eef1f5}
  .bcp-s .form-group,.bcp-s .shiny-input-container{margin:0!important;width:100%!important}
  .bcp-s .form-control,.bcp-s .form-select,.bcp-s select{font-size:13px;padding:3px 8px;min-height:30px}
  .bcp-s .selectize-input{min-height:30px!important;padding:4px 26px 4px 8px!important;font-size:13px!important;line-height:20px!important;box-shadow:none}
  .bcp-s .selectize-dropdown{font-size:13px}
  .bcp-s .radio,.bcp-s .form-check{margin:0;font-size:13px}
  .bcp .help-block,.bcp>.shiny-html-output{font-size:12px;padding:6px 10px}
"))
.bc_display_ui <- function(ns) shiny::tagList(
          shiny::radioButtons(ns("dens_mode"), "Distribution plot", c("Density curves" = "dens", "Histogram" = "hist", "Both" = "both"), inline = TRUE),
          shiny::sliderInput(ns("dens_adj"), "Density smoothing", 0.3, 3, 1, step = 0.1),
          shiny::checkboxInput(ns("dens_thr"), "Mark the suggested threshold", TRUE),
          shiny::radioButtons(ns("box_style"), "Box plot", c("Box" = "box", "Violin + box" = "violin"), inline = TRUE),
          shiny::checkboxInput(ns("box_pts"), "Show the pairs as points", TRUE),
          shiny::checkboxInput(ns("ref_thr"), "Show the reference cut-offs of the marker (distribution plot)", TRUE),
          shiny::selectInput(ns("pc_col"), "PCoA colour", c("Species" = "species", "Genus" = "genus", "Ingroup / outgroup" = "role")),
          shiny::checkboxInput(ns("pc_out"), "Include the outgroup in the PCoA", FALSE),
          shiny::selectInput(ns("pc_axes"), "PCoA axes", c("1 and 2" = "1,2", "1 and 3" = "1,3", "2 and 3" = "2,3")),
          shiny::checkboxInput(ns("pc_ell"), "95% ellipses (groups with 4+ sequences)", TRUE),
          shiny::checkboxInput(ns("pc_lab"), "Sequence labels", FALSE),
          shiny::radioButtons(ns("heat_level"), "Distance matrix", c("Sequences" = "seq", "Species means" = "sp"), inline = TRUE),
          shiny::radioButtons(ns("heat_order"), "Order of the sequences", c("By group" = "group", "By clustering" = "cluster"), inline = TRUE),
          shiny::selectInput(ns("tree_col"), "Tree colour", c("Species" = "species", "Genus" = "genus")),
          shiny::selectInput(ns("tree_layout"), "Tree layout", c("Phylogram" = "phylogram", "Cladogram" = "cladogram", "Unrooted" = "unrooted", "Fan" = "fan", "Radial" = "radial")),
          shiny::selectInput(ns("tree_root"), "Tree root", c("Outgroup (if any)" = "outgroup", "Midpoint" = "midpoint", "None" = "none")),
          shiny::selectInput(ns("tree_lab"), "Tip label", c("Sequence name" = "sequence", "Species" = "species")),
          shiny::numericInput(ns("thr_user"), "Threshold for best close match (empty = suggested)", NA, min = 0, step = 0.005))

# ==============================================================================
# UI
# ==============================================================================
barcodeUI <- function(id) {
  ns <- shiny::NS(id)
  shiny::div(class = "container-fluid",
    .bc_style(),
    # lets the Reset buttons put a colour back into a colour box
    shiny::tags$script(shiny::HTML("(function w(){ if (!(window.Shiny && Shiny.addCustomMessageHandler)) return setTimeout(w, 30);
      Shiny.addCustomMessageHandler('bc_setcolor', function(m){ var e = document.getElementById(m.id);
        if (e) { e.value = m.value; Shiny.setInputValue(m.id, m.value); } }); })();")),
    shiny::fluidRow(class = "sg-siderow",
      shiny::column(4,
        .card("Taxa and groups",
          shiny::selectInput(ns("src"), "Alignment from", c("An aligned FASTA file (open it here)" = "file")),
          shiny::conditionalPanel(sprintf("input['%s'] == 'file'", ns("src")),
            shiny::fileInput(ns("fa"), NULL, buttonLabel = "Open aligned FASTA...", placeholder = "",
                             accept = c(".fa", ".fasta", ".fas", ".fna", ".aln", ".txt"))),
          shiny::uiOutput(ns("src_info")),
          shiny::actionButton(ns("open_groups"), "Select taxa and groups...", icon = shiny::icon("users-rectangle"), class = "btn-primary"),
          shiny::uiOutput(ns("grp_summary"))),
        shiny::conditionalPanel(sprintf("output['%s']", ns("has_data")),
          shiny::conditionalPanel(sprintf("!input['%s'] || !document.documentElement.classList.contains('sg-embed')", ns("analysis")),
          .card("Analysis", class = "sg-nofold",
            shiny::div(class = "bc-an", shiny::radioButtons(ns("analysis"), NULL, selected = character(0),
              choiceValues = names(.bc_analyses()),
              choiceNames = unname(lapply(.bc_analyses(), function(i)
                shiny::tagList(shiny::span(class = "bc-ai", shiny::icon(i$icon)),
                               shiny::span(class = "bc-at", i$title, shiny::tags$small(class = "bc-ad", i$desc))))))))),
          shiny::conditionalPanel(sprintf("input['%s']", ns("analysis")),
            .card("Analysis preferences", class = "sg-nofold",
              shiny::div(class = "bcp",
                .bcp_sec("Analysis"),
                .bcp_row("Scope", shiny::uiOutput(ns("scope_txt"))),
                .bcp_row("Marker (gene region)", shiny::selectInput(ns("marker"), NULL, .bc_marker_choices(), "other", width = "100%")),
                shiny::uiOutput(ns("marker_note")),
                shiny::conditionalPanel(sprintf("['within','between','net','barcode'].indexOf(input['%s']) >= 0", ns("analysis")),
                  .bcp_row("Compare at the level of", shiny::selectInput(ns("level"), NULL,
                    c("Species" = "species", "Genus" = "genus"), width = "100%"))),
                shiny::conditionalPanel(sprintf("input['%s'] == 'barcode'", ns("analysis")),
                  shiny::uiOutput(ns("focus_ui"))),
                .bcp_sec("Estimate variance"),
                .bcp_row("Bootstrap replicates (0 = none)", shiny::numericInput(ns("nboot"), NULL, 0, min = 0, max = 1000, step = 50, width = "100%")),
                .bcp_sec("Substitution model"),
                .bcp_row("Substitutions type", shiny::tags$em("Nucleotide")),
                .bcp_row("Model / method", shiny::selectInput(ns("model"), NULL,
                  c("Kimura 2-parameter (K2P)" = "K80", "p-distance (raw)" = "raw", "Jukes-Cantor (JC69)" = "JC69",
                    "Tamura 3-parameter (T92)" = "T92", "Tamura-Nei (TN93)" = "TN93"), width = "100%")),
                .bcp_sec("Rates and patterns"),
                .bcp_row("Gamma shape (0 = same rate)", shiny::numericInput(ns("gamma"), NULL, 0, min = 0, step = 0.1, width = "100%")),
                .bcp_sec("Data subset to use"),
                .bcp_row("Gaps / missing data", shiny::radioButtons(ns("del"), NULL, c("Pairwise deletion" = "pair", "Complete deletion" = "complete"))),
                shiny::uiOutput(ns("codon_ui")),
                shiny::conditionalPanel(sprintf("input['%s'] == 'barcode'", ns("analysis")),
                  .bcp_sec("Grouping tests"),
                  .bcp_row("Include the outgroup", shiny::checkboxInput(ns("gap_out"), "in the barcode gap and statistics", FALSE)),
                  .bcp_row("Permutations", shiny::numericInput(ns("nperm"), NULL, 999, min = 99, max = 9999, step = 100, width = "100%")),
                  .bcp_row("Random seed", shiny::numericInput(ns("seed"), NULL, 1, min = 1, step = 1, width = "100%")))),
              shiny::actionButton(ns("run"), "Compute", icon = shiny::icon("play"), class = "btn-primary")),
            shiny::conditionalPanel(sprintf("input['%s'] == 'barcode'", ns("analysis")),
              .card("Display", class = "sg-nofold", .bc_display_ui(ns)))))),
      shiny::column(8,
        shiny::uiOutput(ns("msg")),
        shiny::uiOutput(ns("boxes")),
        shiny::tabsetPanel(id = ns("tabs"),
          shiny::tabPanel("Result", value = "mega", shiny::br(), shiny::uiOutput(ns("mega_head")), DT::DTOutput(ns("tbl_mega")),
            shiny::br(), shiny::downloadButton(ns("dl_mega"), "Table (CSV)")),
          shiny::tabPanel("Distribution", value = "dens", shiny::br(), shiny::uiOutput(ns("fb_dens")), shiny::plotOutput(ns("p_dens"), height = "400px"),
            shiny::helpText("Distances between sequences of the same species (intraspecific), of different species of the same genus (congeneric), ",
                            "and of different genera (intergeneric; an outgroup of another genus is counted here when it is included). ",
                            "A barcode gap exists when the intraspecific curve ends before the other curves begin.")),
          shiny::tabPanel("Box plot", value = "box", shiny::br(), shiny::uiOutput(ns("fb_box")), shiny::plotOutput(ns("p_box"), height = "430px")),
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
          shiny::tabPanel("Group distances", value = "gd", shiny::br(),
            shiny::radioButtons(ns("gd_level"), NULL, c("Species (the compared groups)" = "unit", "Genera (the next level)" = "parent"), inline = TRUE),
            shiny::h5("Mean distance within each group"), DT::DTOutput(ns("tbl_gd_in")),
            shiny::h5("Mean distance between groups, and the net distance"), DT::DTOutput(ns("tbl_gd_bt")),
            shiny::helpText("Mean_between = the mean of all distances between a sequence of one group and a sequence of the other. ",
                            "Net_between = Mean_between - (mean within group 1 + mean within group 2) / 2: the divergence that remains after the variation inside the groups is removed. ",
                            "SE = bootstrap standard error over the sites (fill in 'Bootstrap replicates' in the preferences; empty otherwise). A group with one sequence has no within-group mean.")),
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
        shiny::conditionalPanel(sprintf("!input['%1$s'] || input['%1$s'] == 'barcode'", ns("analysis")),
        sg_dlgroup(label = "Save",
          shiny::downloadButton(ns("dl_xlsx"), "Excel (all tables)"),
          shiny::downloadButton(ns("dl_species"), "Species (CSV)"),
          shiny::downloadButton(ns("dl_pairs"), "Pairs (CSV)"),
          shiny::div(style = "width:90px", shiny::selectInput(ns("fig_fmt"), "Figure", c("png", "pdf", "svg", "jpg"), width = "100%")),
          shiny::div(style = "width:70px", shiny::numericInput(ns("fig_w"), "Width", 8, min = 2, max = 30, width = "100%")),
          shiny::div(style = "width:70px", shiny::numericInput(ns("fig_h"), "Height", 6, min = 2, max = 30, width = "100%")),
          shiny::downloadButton(ns("dl_fig"), "Download the figure of the open tab"))))))
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
    rv <- shiny::reactiveValues(groups = NULL, targets = character(), gver = 0L, res = NULL, msg = NULL, file_aln = NULL)
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
    # Qt window, Navigation area: "Open in Genetic distances" sends the path of an aligned FASTA file on this computer
    shiny::observeEvent(input$qt_aln_path, {
      a <- tryCatch(.bc_read_fasta(input$qt_aln_path), error = function(e) NULL)
      if (!ok_aln(a)) { set_msg("The file needs at least 3 sequences in FASTA format.", "warning"); return() }
      if (length(unique(nchar(a))) != 1L) { set_msg("The sequences in this file do not all have the same length, so it is not an alignment. Align it first in the MSA tab.", "danger"); return() }
      rv$file_aln <- a; set_msg(NULL)
      shiny::updateSelectInput(session, "src", selected = "file")
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
      fa <- rep("", length(nm))
      if (keep && !is.null(old) && "Family" %in% names(old)) { m <- match(nm, old$Sequence); h <- !is.na(m); fa[h] <- old$Family[m[h]] }
      rv$groups <- data.frame(Sequence = nm, Species = sp, Genus = ge, Family = fa, Role = ro, stringsAsFactors = FALSE)
      rv$targets <- if (keep) intersect(rv$targets, unique(sp)) else character()
      rv$gver <- rv$gver + 1L
    }
    msg2 <- function(txt, type = "info") {
      set_msg(txt, type)
      if (!is.null(txt)) shiny::showNotification(txt, type = switch(type, danger = "error", warning = "warning", "message"), duration = 12)
    }
    apply_groups <- function() {
      if (is.null(src())) { msg2("Choose or open an alignment first.", "warning"); return(invisible()) }
      mode <- .bc_nz(input$gsrc, "names")
      if (identical(mode, "names")) { detect(FALSE); msg2(NULL); return(invisible()) }
      gt <- NULL; label <- ""
      if (identical(mode, "file")) {
        f <- input$gfile
        if (is.null(f)) { msg2("Open a group table first.", "warning"); return(invisible()) }
        gt <- .bc_as_groups(.bc_read_table(f$datapath)); label <- "the group table"
      } else if (startsWith(mode, "g:")) {
        label <- sub("^g:", "", mode); gt <- .bc_as_groups(group_sources()[[label]])
      } else if (startsWith(mode, "t:")) {
        nm <- sub("^t:", "", mode); thr <- .bc_nz(input$gthr, 0.02)
        gt <- .bc_tree_groups(tree_sources()[[nm]], thr); label <- sprintf("%s, clades at distance <= %s", nm, thr)
      }
      if (is.null(gt) || !nrow(gt)) {
        msg2("No usable groups there: a table needs a column with the sequence names and a column with the group.", "warning")
        return(invisible())
      }
      detect(FALSE)                                      # start from the names, then overwrite what matches
      g <- rv$groups
      m <- match(.bc_key(g$Sequence), .bc_key(gt$Sequence))
      idx <- which(!is.na(m))
      val <- trimws(gt$Species[m[idx]]); good <- !is.na(val) & nzchar(val)
      idx <- idx[good]; val <- val[good]
      if (!length(idx)) { msg2(sprintf("None of the sequence names of the alignment was found in %s.", label), "warning"); return(invisible()) }
      if (!any(grepl("[A-Za-z]", val))) val <- paste0("Group_", val)       # groups numbered 1, 2, 3 ...
      g$Species[idx] <- val
      if ("Genus" %in% names(gt)) { x <- gt$Genus[m[idx]]; g$Genus[idx] <- ifelse(is.na(x) | !nzchar(x), g$Genus[idx], x) }
      if ("Family" %in% names(gt)) { x <- gt$Family[m[idx]]; g$Family[idx] <- ifelse(is.na(x) | !nzchar(x), g$Family[idx], x) }
      if ("Role" %in% names(gt)) g$Role[idx] <- gt$Role[m[idx]]
      rv$groups <- g; rv$gver <- rv$gver + 1L
      msg2(sprintf("Groups taken from %s: %d of %d sequences found, %d groups.%s", label, length(idx), nrow(g), length(unique(val)),
                      if (length(idx) < nrow(g)) " The other sequences keep the group read from their names (edit them in the table)." else ""),
              if (length(idx) < nrow(g)) "warning" else "info")
    }
    shiny::observeEvent(input$detect, apply_groups())
    shiny::observeEvent(src(), {
      if (is.null(src())) { rv$groups <- NULL; rv$res <- NULL; rv$targets <- character(); rv$gver <- rv$gver + 1L }
      else detect(TRUE)
    }, ignoreNULL = FALSE)
    # ---- the "Select taxa and groups" window -------------------------------------------------------------
    push_choices <- function() shiny::isolate({
      ch <- c("The sequence names" = "names", "A group table (CSV / TSV)" = "file")
      gs <- group_sources()
      for (n in names(gs)) if (!is.null(gs[[n]])) ch <- c(ch, stats::setNames(paste0("g:", n), n))
      ts <- tree_sources()
      for (n in names(ts)) if (!is.null(ts[[n]])) ch <- c(ch, stats::setNames(paste0("t:", n), paste0(n, " (clades)")))
      shiny::updateSelectInput(session, "gsrc", choices = ch, selected = "names")
      push_pick()
    })
    push_pick <- function() shiny::isolate({
      g <- rv$groups; a <- src()
      pk <- list()
      if (!is.null(g)) pk <- list(Species = sort(unique(g$Species)), Genus = sort(unique(g$Genus)),
                                  Family = sort(setdiff(unique(g$Family), "")), Role = c("Outgroup", "Ingroup"))
      pk <- Filter(length, pk)
      shiny::updateSelectizeInput(session, "gm_pick", choices = c(list(""), pk), selected = "", server = FALSE)
      ch <- if (is.null(a)) character() else {
        l <- list(Sequences = unname(names(a)))
        if (!is.null(g)) { l$Species <- sort(unique(g$Species)); l$Genera <- sort(unique(g$Genus)) }
        l
      }
      shiny::updateSelectizeInput(session, "outgroup", choices = ch, selected = character(0), server = FALSE)
    })
    groups_modal <- function() {
      ns <- session$ns
      shiny::modalDialog(size = "xl", easyClose = TRUE, title = "Select taxa and groups",
        footer = shiny::tagList(shiny::downloadButton(ns("dl_groups"), "Export the groups (CSV)"), shiny::modalButton("Done")),
        shiny::tagList(shiny::br(),
            shiny::fluidRow(
              shiny::column(8, DT::DTOutput(ns("tbl"))),
              shiny::column(4,
                shiny::div(class = "bcg-step", shiny::span("1"), "Select the sequences in the table"),
                shiny::helpText("Click rows (Ctrl or Shift for several), or pick a whole species, genus or group:"),
                shiny::selectizeInput(ns("gm_pick"), NULL, choices = "", options = list(placeholder = "Select the members of a group...")),
                shiny::div(style = "display:flex;gap:6px;flex-wrap:wrap;margin-bottom:4px",
                  shiny::actionButton(ns("gm_all"), "All", class = "btn-sm"), shiny::actionButton(ns("gm_none"), "None", class = "btn-sm"),
                  shiny::actionButton(ns("gm_inv"), "Invert", class = "btn-sm"), shiny::uiOutput(ns("gm_count"), inline = TRUE)),
                shiny::div(class = "bcg-step", shiny::span("2"), "Press what they are"),
                shiny::div(class = "bcg-btns",
                  shiny::actionButton(ns("gm_target"), "Target species", class = "bcg-b bcg-t",
                    title = "The selected sequences become ONE species and the species you want to test (intraspecific distances are measured inside it)"),
                  shiny::actionButton(ns("gm_cong"), "Congeneric", class = "bcg-b bcg-c",
                    title = "Other species of the same genus as the target species"),
                  shiny::actionButton(ns("gm_other"), "Other genus", class = "bcg-b bcg-i",
                    title = "Species of another genus (intergeneric)"),
                  shiny::actionButton(ns("gm_out"), "Outgroup", class = "bcg-b bcg-o",
                    title = "Left out of the barcode gap and the grouping tests unless you include it in the preferences")),
                shiny::textInput(ns("gm_name"), NULL, placeholder = "Name of the target species (optional; the name in the table is used)"),
                shiny::div(style = "display:flex;gap:6px;flex-wrap:wrap;margin:4px 0 8px",
                  shiny::actionButton(ns("gm_reset"), "Reset the selected (read from the names again)", class = "btn-sm", icon = shiny::icon("rotate-left"))),
                shiny::uiOutput(ns("gm_summary")),
                shiny::tags$details(style = "margin-top:8px",
                  shiny::tags$summary(style = "cursor:pointer;color:#52606d", "More: name a species, genus or family yourself"),
                  shiny::helpText("Type a name above, select the sequences, then press:"),
                  shiny::div(style = "display:flex;gap:6px;flex-wrap:wrap;margin-bottom:6px",
                    shiny::actionButton(ns("gm_sp"), "One species", class = "btn-sm"), shiny::actionButton(ns("gm_ge"), "One genus", class = "btn-sm"),
                    shiny::actionButton(ns("gm_fa"), "One family", class = "btn-sm"), shiny::actionButton(ns("gm_in"), "Ingroup", class = "btn-sm"),
                    shiny::actionButton(ns("gm_ungroup"), "Ungroup", class = "btn-sm")),
                  shiny::helpText("Species, Genus and Family can also be edited by double-clicking a cell. The family is only needed when you compare genera.")))),
            shiny::tags$style(shiny::HTML(paste0(
              ".bcg-step{font-weight:600;margin:6px 0 2px;display:flex;align-items:center;gap:8px}",
              ".bcg-step span{background:#2c7be5;color:#fff;border-radius:50%;width:22px;height:22px;display:inline-flex;align-items:center;justify-content:center;font-size:12px}",
              ".bcg-btns{display:grid;grid-template-columns:1fr 1fr;gap:8px;margin:6px 0 8px}",
              ".bcg-b{font-weight:600;padding:9px 6px;border-width:2px}",
              ".bcg-t{border-color:#2c7be5;color:#2c7be5}.bcg-t:hover{background:#2c7be5;color:#fff}",
              ".bcg-c{border-color:#f5a623;color:#b9750a}.bcg-c:hover{background:#f5a623;color:#fff}",
              ".bcg-i{border-color:#e5484d;color:#e5484d}.bcg-i:hover{background:#e5484d;color:#fff}",
              ".bcg-o{border-color:#7b8794;color:#52606d}.bcg-o:hover{background:#7b8794;color:#fff}")))))
    }
    shiny::observeEvent(input$open_groups, {
      if (is.null(src())) { set_msg("Choose or open an alignment first (card 1).", "warning"); return() }
      if (is.null(rv$groups)) detect(FALSE)
      shiny::showModal(groups_modal())
      session$onFlushed(function() push_choices(), once = TRUE)
    })
    shiny::observeEvent(input$apply_out, {
      g <- rv$groups; shiny::req(g)
      g <- .bc_set_outgroup(g, input$outgroup)
      rv$groups <- g; rv$gver <- rv$gver + 1L
      shiny::showNotification(sprintf("%d sequence(s) are now the outgroup.", sum(g$Role == "Outgroup")), type = "message")
    })
    shiny::observeEvent(input$clear_out, {
      g <- rv$groups; shiny::req(g)
      g$Role <- "Ingroup"; rv$groups <- g; rv$gver <- rv$gver + 1L
      shiny::updateSelectizeInput(session, "outgroup", selected = character(0))
    })

    # what every sequence is, relative to the target species: Target / Congeneric / Other genus / Outgroup
    status_of <- function(g, targets) {
      if (is.null(g)) return(character())
      if (!length(targets)) return(ifelse(g$Role == "Outgroup", "Outgroup", "-"))
      tg <- unique(g$Genus[g$Species %in% targets])
      ifelse(g$Role == "Outgroup", "Outgroup", ifelse(g$Species %in% targets, "Target",
             ifelse(g$Genus %in% tg, "Congeneric", "Other genus")))
    }
    disp <- function(g) data.frame(Sequence = g$Sequence, Status = status_of(g, rv$targets), Species = g$Species, Genus = g$Genus,
                                   Family = g$Family, Role = g$Role, stringsAsFactors = FALSE)
    output$tbl <- DT::renderDT({
      g <- shiny::isolate(rv$groups); shiny::req(g)
      d <- DT::datatable(disp(g), rownames = FALSE, selection = list(mode = "multiple", target = "row"),
                    editable = list(target = "cell", disable = list(columns = c(0, 1))),
                    options = list(pageLength = -1, lengthMenu = list(c(10, 50, -1), c("10", "50", "All")), dom = "lfrtip",
                                   scrollX = TRUE, scrollY = "430px", scrollCollapse = TRUE,
                                   columnDefs = list(list(targets = 0, render = DT::JS(
                                     "function(d, t){ return (t === 'display' && d.length > 55) ? '<span title=\"' + d + '\">' + d.substr(0, 55) + '&hellip;</span>' : d; }")))))
      DT::formatStyle(d, "Status", fontWeight = "bold",
        color = DT::styleEqual(c("Target", "Congeneric", "Other genus", "Outgroup"), c("#2c7be5", "#b9750a", "#e5484d", "#52606d")))
    })
    tbl_proxy <- DT::dataTableProxy("tbl")
    shiny::observeEvent(list(rv$gver, rv$targets), {
      g <- rv$groups; if (!is.null(g)) DT::replaceData(tbl_proxy, disp(g), resetPaging = FALSE, rownames = FALSE, clearSelection = "none")
    }, ignoreInit = TRUE)
    shiny::observeEvent(input$tbl_cell_edit, {
      info <- input$tbl_cell_edit; g <- rv$groups
      nm <- c("Sequence", "Status", "Species", "Genus", "Family", "Role")[as.integer(info$col) + 1L]
      if (!(nm %in% c("Species", "Genus", "Family", "Role"))) return()
      v <- trimws(as.character(info$value)); if (identical(nm, "Role")) v <- .bc_role(v) else v <- gsub("[ ]+", "_", v)
      g[[nm]][as.integer(info$row)] <- v
      rv$groups <- g; rv$gver <- rv$gver + 1L
    })
    sel_rows <- function() {
      i <- input$tbl_rows_selected
      if (!length(i)) { shiny::showNotification("Select one or more sequences in the table first.", type = "warning"); return(NULL) }
      i
    }
    commit_groups <- function(g) {
      rv$targets <- intersect(rv$targets, unique(g$Species))          # a target species that no longer exists is dropped
      rv$groups <- g; rv$gver <- rv$gver + 1L
    }
    top <- function(x) { x <- x[!is.na(x) & nzchar(x)]; if (!length(x)) "" else names(sort(table(x), decreasing = TRUE))[1] }
    set_col <- function(col, val) {
      i <- sel_rows(); if (is.null(i)) return()
      g <- rv$groups; g[[col]][i] <- val; commit_groups(g)
    }
    name_ok <- function() {
      nm <- trimws(.bc_nz(input$gm_name, ""))
      if (!nzchar(nm)) { shiny::showNotification("Type the name first.", type = "warning"); return(NULL) }
      gsub("[ ]+", "_", nm)
    }
    shiny::observeEvent(input$gm_sp, { n <- name_ok(); if (!is.null(n)) set_col("Species", n) })
    shiny::observeEvent(input$gm_ge, { n <- name_ok(); if (!is.null(n)) set_col("Genus", n) })
    shiny::observeEvent(input$gm_fa, { n <- name_ok(); if (!is.null(n)) set_col("Family", n) })
    shiny::observeEvent(input$gm_in, set_col("Role", "Ingroup"))
    shiny::observeEvent(input$gm_out, set_col("Role", "Outgroup"))
    shiny::observeEvent(input$gm_ungroup, {
      i <- sel_rows(); if (is.null(i)) return()
      g <- rv$groups; g$Species[i] <- g$Sequence[i]; commit_groups(g)
    })
    # Target species: the selected sequences become one species (the typed name, else the commonest name among them)
    shiny::observeEvent(input$gm_target, {
      i <- sel_rows(); if (is.null(i)) return()
      g <- rv$groups
      nm <- trimws(.bc_nz(input$gm_name, "")); nm <- if (nzchar(nm)) gsub("[ ]+", "_", nm) else top(g$Species[i])
      g$Species[i] <- nm; g$Genus[i] <- top(g$Genus[i]); g$Role[i] <- "Ingroup"
      rv$targets <- union(rv$targets, nm)
      commit_groups(g)
      shiny::updateTextInput(session, "gm_name", value = "")
      shiny::showNotification(sprintf("Target species '%s': %d sequences (genus %s).", nm, length(i), g$Genus[i[1]]), type = "message")
    })
    shiny::observeEvent(input$gm_cong, {
      g <- rv$groups
      if (!length(rv$targets)) { shiny::showNotification("Choose the target species first (select its sequences, press 'Target species').", type = "warning"); return() }
      i <- sel_rows(); if (is.null(i)) return()
      skip <- g$Species[i] %in% rv$targets
      i <- i[!skip]
      if (!length(i)) { shiny::showNotification("Those sequences belong to a target species.", type = "warning"); return() }
      gen <- top(g$Genus[g$Species %in% rv$targets])
      g$Genus[i] <- gen; g$Role[i] <- "Ingroup"
      commit_groups(g)
      shiny::showNotification(sprintf("%d sequence(s) are now congeneric (genus %s).%s", length(i), gen,
        if (any(skip)) sprintf(" %d target sequence(s) were skipped.", sum(skip)) else ""), type = "message")
    })
    shiny::observeEvent(input$gm_other, {
      g <- rv$groups
      i <- sel_rows(); if (is.null(i)) return()
      tgen <- unique(g$Genus[g$Species %in% rv$targets])
      own <- .bc_group(g$Sequence[i], "genus")
      own <- ifelse(own %in% tgen, paste0(own, "_other"), own)
      g$Genus[i] <- own; g$Role[i] <- "Ingroup"
      i2 <- i[g$Species[i] %in% rv$targets]; if (length(i2)) rv$targets <- setdiff(rv$targets, unique(g$Species[i2]))
      commit_groups(g)
    })
    shiny::observeEvent(input$gm_reset, {
      i <- sel_rows(); if (is.null(i)) return()
      g <- rv$groups; nm <- g$Sequence[i]
      g$Species[i] <- .bc_group(nm, .bc_nz(input$rule, "species")); g$Genus[i] <- .bc_group(nm, "genus"); g$Family[i] <- ""; g$Role[i] <- "Ingroup"
      commit_groups(g)
    })
    shiny::observeEvent(input$gm_all, DT::selectRows(tbl_proxy, seq_len(nrow(rv$groups))))
    shiny::observeEvent(input$gm_none, DT::selectRows(tbl_proxy, NULL))
    shiny::observeEvent(input$gm_inv, {
      n <- nrow(rv$groups); DT::selectRows(tbl_proxy, setdiff(seq_len(n), input$tbl_rows_selected))
    })
    shiny::observeEvent(input$gm_pick, {
      v <- input$gm_pick; g <- rv$groups
      if (is.null(v) || !nzchar(v) || is.null(g)) return()
      hit <- which(g$Species == v | g$Genus == v | g$Family == v | (v %in% c("Outgroup", "Ingroup") & g$Role == v))
      DT::selectRows(tbl_proxy, hit)
    })
    output$gm_count <- shiny::renderUI({
      n <- length(input$tbl_rows_selected)
      shiny::span(style = "color:#6b7785;font-size:12.5px;align-self:center", if (n) sprintf("%d selected", n) else "nothing selected")
    })
    # counts of the groups (shown in the window and under the button of the card)
    grp_counts <- function(g) {
      if (is.null(g)) return(NULL)
      list(n = nrow(g), sp = length(unique(g$Species)), ge = length(unique(g$Genus)), fa = length(setdiff(unique(g$Family), "")),
           out = sum(g$Role == "Outgroup"), single = sum(table(g$Species) == 1L))
    }
    output$gm_summary <- shiny::renderUI({
      rv$gver; g <- rv$groups; shiny::req(g); st <- status_of(g, rv$targets); tg <- rv$targets
      chip <- function(txt, col) shiny::span(style = sprintf("display:inline-block;background:%s;color:#fff;border-radius:12px;padding:2px 10px;margin:2px 4px 2px 0;font-size:12.5px", col), txt)
      if (!length(tg))
        return(shiny::div(style = "background:#f6f8fa;border-radius:8px;padding:8px 10px;font-size:12.5px;color:#52606d",
          "No target species yet. Without one, every species is compared with every other (the whole data set); ",
          "with one, the results are about that species: distances inside it, to its congeners and to other genera.",
          shiny::div(style = "margin-top:4px", sprintf("%d sequences | %d species | %d genera | %d outgroup", nrow(g), length(unique(g$Species)), length(unique(g$Genus)), sum(g$Role == "Outgroup")))))
      shiny::div(style = "background:#f6f8fa;border-radius:8px;padding:8px 10px;font-size:12.5px",
        shiny::div(chip(sprintf("Target: %s (%d)", paste(tg, collapse = ", "), sum(st == "Target")), "#2c7be5")),
        shiny::div(chip(sprintf("Congeneric: %d seq. in %d species", sum(st == "Congeneric"), length(unique(g$Species[st == "Congeneric"]))), "#f5a623"),
                   chip(sprintf("Other genus: %d seq.", sum(st == "Other genus")), "#e5484d"),
                   chip(sprintf("Outgroup: %d seq.", sum(st == "Outgroup")), "#7b8794")))
    })
    output$grp_summary <- shiny::renderUI({
      rv$gver; g <- rv$groups
      if (is.null(g)) return(shiny::helpText("No groups yet: open an alignment, then press the button."))
      k <- grp_counts(g); st <- status_of(g, rv$targets)
      shiny::tagList(
        shiny::div(style = "margin-top:6px;font-weight:600", sprintf("%d taxa | %d species | %d genera | %d outgroup", k$n, k$sp, k$ge, k$out)),
        if (length(rv$targets)) shiny::div(style = "font-size:12.5px;color:#2c7be5",
          sprintf("Target species: %s (%d seq.), congeneric %d seq., other genus %d seq.", paste(rv$targets, collapse = ", "),
                  sum(st == "Target"), sum(st == "Congeneric"), sum(st == "Other genus"))),
        if (k$single == k$sp && k$sp > 1L)
          shiny::div(class = "text-danger", style = "font-size:12.5px", "Every taxon is its own species: there is no intraspecific distance. Check the names, or group them (button above)."))
    })
    output$dl_groups <- shiny::downloadHandler(function() "groups.csv",
      function(file) { shiny::req(rv$groups); utils::write.csv(rv$groups, file, row.names = FALSE) })

    output$focus_ui <- shiny::renderUI({
      tg <- rv$targets
      if (!length(tg)) return(NULL)
      shiny::checkboxInput(session$ns("focus"), sprintf("Only the target species (%s): its own distances, to congeners and to other genera", paste(utils::head(tg, 3), collapse = ", ")), TRUE)
    })
    # ---- marker advice ------------------------------------------------------------------------------------
    shiny::observeEvent(input$marker, {
      m <- .bc_markers[[.bc_nz(input$marker, "other")]]
      if (!identical(input$marker, "other")) shiny::updateSelectInput(session, "level", selected = m$level)
    }, ignoreInit = TRUE)
    output$marker_note <- shiny::renderUI({
      m <- .bc_markers[[.bc_nz(input$marker, "other")]]
      if (!nzchar(m$note)) return(NULL)
      shiny::div(style = "font-size:12.5px;color:#52606d;background:#f6f8fa;border-radius:8px;padding:8px 10px;margin-bottom:10px",
        shiny::strong(if (identical(m$level, "genus")) "Best compared at the genus level. " else "Usually compared at the species level. "), m$note)
    })
    output$codon_ui <- shiny::renderUI({
      m <- .bc_markers[[.bc_nz(input$marker, "other")]]
      if (!isTRUE(m$coding)) return(NULL)
      shiny::tagList(
        shiny::checkboxGroupInput(session$ns("pos"), "Codon positions to use", c("1st" = 1, "2nd" = 2, "3rd" = 3), selected = 1:3, inline = TRUE),
        shiny::numericInput(session$ns("frame"), "The reading frame starts at column", 1, min = 1, max = 3, step = 1))
    })

    # ---- the analysis menu (like MEGA: Analysis > Distance) -------------------------------------------------
    output$has_data <- shiny::reactive(!is.null(src()))
    shiny::outputOptions(output, "has_data", suspendWhenHidden = FALSE)
    output$scope_txt <- shiny::renderUI(shiny::tags$em(.bc_analyses()[[.bc_nz(input$analysis, "barcode")]]$title))
    old_tabs <- c("dens", "box", "tests", "pcoa", "mat", "gd", "sp", "scat", "thr", "ident", "tree", "pairs")
    shiny::observeEvent(input$set_analysis, {
      if (is.null(src())) { shiny::showNotification("Open an alignment first (card 1), then choose the analysis.", type = "warning"); return() }
      shiny::updateRadioButtons(session, "analysis", selected = input$set_analysis)
    })
    shiny::observeEvent(input$analysis, {
      an <- input$analysis; rv$mega <- NULL
      if (is.null(an) || identical(an, "barcode")) {
        for (t in old_tabs) shiny::showTab("tabs", t, session = session)
        shiny::hideTab("tabs", "mega", session = session)
        if (identical(input$tabs, "mega")) shiny::updateTabsetPanel(session, "tabs", selected = "dens")
      } else {
        for (t in old_tabs) shiny::hideTab("tabs", t, session = session)
        shiny::showTab("tabs", "mega", select = TRUE, session = session)
      }
    }, ignoreNULL = FALSE)

    # the analyses that need no barcode statistics: pairwise, overall mean, within / between / net group means
    run_simple <- function(an, a, g) {
      lvl <- .bc_nz(input$level, "species"); nm0 <- names(a); m <- match(nm0, g$Sequence)
      sp <- trimws(g$Species[m]); ge <- trimws(g$Genus[m])
      sp[is.na(sp) | !nzchar(sp)] <- "(unassigned)"; ge[is.na(ge) | !nzchar(ge)] <- "(unassigned)"
      un <- if (identical(lvl, "genus")) ge else sp
      if (an %in% c("between", "net") && length(unique(un)) < 2L) {
        set_msg("At least two groups are needed. Group the sequences with 'Select taxa and groups'.", "warning"); return() }
      model <- .bc_nz(input$model, "K80"); gam <- suppressWarnings(as.numeric(.bc_nz(input$gamma, 0))); if (is.na(gam)) gam <- 0
      pos <- as.integer(.bc_nz(input$pos, 1:3)); fr <- as.integer(.bc_nz(input$frame, 1L))
      nboot <- as.integer(.bc_nz(input$nboot, 0L)); seed <- as.integer(.bc_nz(input$seed, 1))
      pairwise <- !identical(input$del, "complete")
      r <- shiny::withProgress(message = "Calculating distances ...", value = 0.4, tryCatch({
        a2 <- a; names(a2) <- make.unique(nm0)
        ms <- .bc_sites(a2, pairwise, pos, fr)
        if (ncol(ms) < 10L) stop("fewer than 10 sites are left with these settings (complete deletion / codon positions)")
        D <- .bc_dist_mat(ms, model, gam)
        list(D = D, ms = ms)
      }, error = function(e) { set_msg(paste("Analysis failed:", conditionMessage(e)), "danger"); NULL }))
      if (is.null(r)) return()
      D <- r$D; ms <- r$ms; n <- nrow(D)
      mlab <- c(K80 = "Kimura 2-parameter", raw = "p-distance", JC69 = "Jukes-Cantor", T92 = "Tamura 3-parameter", TN93 = "Tamura-Nei")[[model]]
      head <- sprintf("%s | %d sequences | %d sites | %s%s%s", .bc_analyses()[[an]]$title, n, ncol(ms), mlab,
                      if (gam > 0 && model %in% c("K80", "JC69", "TN93")) sprintf(" + G (%.2f)", gam) else "",
                      if (pairwise) " | pairwise deletion" else " | complete deletion")
      se_note <- if (nboot >= 10L) sprintf("SE = bootstrap over the sites (%d replicates).", nboot) else "No standard errors: set 'Bootstrap replicates' (10 or more) in the preferences."
      drop_se <- function(d) { for (k in grep("^SE", names(d), value = TRUE)) if (all(is.na(d[[k]]))) d[[k]] <- NULL; d }
      tab <- NULL; note <- NULL
      if (identical(an, "pairwise")) {
        tab <- data.frame(Sequence = paste0(seq_len(n), ". ", rownames(D)), stats::setNames(as.data.frame(round(D, 4)), seq_len(n)), check.names = FALSE)
        note <- "Each column number is the sequence with that number in the first column. Empty = no distance (no shared sites)."
      } else if (identical(an, "overall")) {
        gd <- shiny::withProgress(message = "Bootstrap ...", value = 0.6, .bc_group_dist(D, rep("All", n), ms, model, gam, nboot, seed))
        tab <- data.frame(Item = c("Sequences", "Sites", "Model", "Overall mean distance", if (nboot >= 10L) "SE (bootstrap)"),
                          Value = as.character(c(n, ncol(ms), mlab, round(gd$within$Mean_within, 4), if (nboot >= 10L) round(gd$within$SE, 4))), stringsAsFactors = FALSE)
        note <- se_note
      } else {
        gd <- shiny::withProgress(message = "Group distances ...", value = 0.6, .bc_group_dist(D, un, ms, model, gam, nboot, seed))
        tab <- switch(an, within = drop_se(gd$within),
                      between = drop_se(gd$between[, c("Group_1", "Group_2", "Mean_between", "SE")]),
                      net = drop_se(gd$between[, c("Group_1", "Group_2", "Net_between", "SE_net")]))
        note <- paste(se_note, switch(an,
          within = "A group with one sequence has no within-group mean.",
          between = "Mean of all the distances between a sequence of one group and a sequence of the other.",
          net = "Net = mean between - (mean within group 1 + mean within group 2) / 2."))
      }
      rv$mega <- list(an = an, head = head, tab = tab, note = note); set_msg(NULL)
    }

    # ---- calculation --------------------------------------------------------------------------------------
    shiny::observeEvent(input$run, {
      a <- src(); g <- rv$groups
      if (is.null(a) || is.null(g)) { set_msg("Provide an alignment with at least 3 sequences first.", "warning"); return() }
      an <- .bc_nz(input$analysis, "barcode")
      if (!identical(an, "barcode")) { run_simple(an, a, g); return() }
      lvl <- .bc_nz(input$level, "species"); L <- .bc_levels[[lvl]]
      nm0 <- names(a); m <- match(nm0, g$Sequence)
      sp <- trimws(g$Species[m]); ge <- trimws(g$Genus[m]); ro <- .bc_role(g$Role[m])
      fa <- if ("Family" %in% names(g)) trimws(g$Family[m]) else rep("", length(m))
      sp[is.na(sp) | !nzchar(sp)] <- "(unassigned)"; ge[is.na(ge) | !nzchar(ge)] <- "(unassigned)"
      fa[is.na(fa) | !nzchar(fa)] <- paste0("(no family) ", ge[is.na(fa) | !nzchar(fa)])   # no family given: every genus stands alone
      un <- if (identical(lvl, "genus")) ge else sp           # the compared groups
      pa <- if (identical(lvl, "genus")) fa else ge           # the level above them
      gap_out <- isTRUE(input$gap_out)
      # target species chosen in 'Select taxa and groups' (at the genus level: their genera)
      tg <- NULL
      if (length(rv$targets) && !identical(input$focus, FALSE)) {
        tg <- unique(un[sp %in% rv$targets]); if (!length(tg)) tg <- NULL
      }
      use <- if (gap_out) rep(TRUE, length(un)) else ro != "Outgroup"
      if (length(unique(un[use])) < 2) {
        set_msg(sprintf("At least two %s are needed in the ingroup (or tick 'Include the outgroup').", L$unit), "warning"); return()
      }
      if (identical(lvl, "species") && all(table(un[use]) == 1L)) {
        set_msg("Every sequence is its own species, so there are no intraspecific distances. Group the sequences with 'Select taxa and groups'.", "warning"); return()
      }
      a2 <- a; names(a2) <- make.unique(nm0)
      model <- .bc_nz(input$model, "K80"); gam <- suppressWarnings(as.numeric(.bc_nz(input$gamma, 0))); if (is.na(gam)) gam <- 0
      pairwise <- !identical(input$del, "complete")
      pos <- as.integer(.bc_nz(input$pos, 1:3)); fr <- as.integer(.bc_nz(input$frame, 1L))
      nboot <- as.integer(.bc_nz(input$nboot, 0L)); seed <- as.integer(.bc_nz(input$seed, 1))
      r <- shiny::withProgress(message = "Calculating distances and tests ...", value = 0.3, {
        tryCatch({
          ms <- .bc_sites(a2, pairwise, pos, fr)
          if (ncol(ms) < 10L) stop("fewer than 10 sites are left with these settings (complete deletion / codon positions)")
          D <- .bc_dist_mat(ms, model, gam)
          out <- .bc_run(D, un, pa, ro, model, gap_out, as.integer(.bc_nz(input$nperm, 999)), seed, L$cats, tg)
          shiny::setProgress(0.7, message = if (nboot >= 10L) sprintf("Bootstrap (%d replicates) ...", nboot) else "Group distances ...")
          ix <- if (gap_out) seq_along(un) else which(ro != "Outgroup")
          sub <- function(v) v[ix]
          Dg <- D[ix, ix, drop = FALSE]; mg <- ms[ix, , drop = FALSE]
          out$gd <- list(unit = .bc_group_dist(Dg, sub(un), mg, model, gam, nboot, seed),
                         parent = .bc_group_dist(Dg, sub(pa), mg, model, gam, 0L, seed))
          out$sites <- ncol(ms); out$gamma <- gam; out$pairwise <- pairwise
          out
        }, error = function(e) { set_msg(paste("Analysis failed:", conditionMessage(e)), "danger"); NULL })
      })
      if (is.null(r)) return()
      rv$res <- r
      shiny::updateSelectInput(session, "pc_col", choices = stats::setNames(c("species", "genus", "role"), c(L$ucol, L$pcol, "Ingroup / outgroup")),
                               selected = shiny::isolate(input$pc_col))
      notes <- character()
      if (gam > 0 && !(model %in% c("K80", "JC69", "TN93")))
        notes <- c(notes, "The gamma shape was ignored: it is available for K2P, Jukes-Cantor and Tamura-Nei only.")
      if (length(r$dropped))
        notes <- c(notes, sprintf("%d sequence(s) with missing distances were left out of the PCoA, PERMANOVA/ANOSIM and the tree: %s.",
                                  length(r$dropped), paste(utils::head(r$dropped, 5), collapse = ", ")))
      if (identical(lvl, "species")) {
        bad <- unique(r$sp[!grepl("^[A-Za-z][A-Za-z.-]*_[A-Za-z]", r$sp)])
        if (length(bad))
          notes <- c(notes, sprintf("These species names do not look like Genus_species: %s. Edit them in 'Select taxa and groups' (a sequence without a species name counts as its own species).",
                                    paste(substr(utils::head(bad, 4), 1L, 50L), collapse = ", ")))
      } else if (all(grepl("^\\(no family\\)", r$ge)))
        notes <- c(notes, "No family was given, so every genus stands alone: the middle category (same family) is empty. Add families in 'Select taxa and groups' if you want it.")
      z <- r$pairs[r$pairs$Use & !is.na(r$pairs$Distance) & r$pairs$Distance == 0 & r$pairs$Species_1 != r$pairs$Species_2, , drop = FALSE]
      if (nrow(z))
        notes <- c(notes, sprintf("%d pair(s) of different %s have distance 0 (identical sequences), which removes the barcode gap, for example %s vs %s (%s / %s). Check the labels, or the marker may be too conserved.",
                                  nrow(z), L$unit, z$Sequence_1[1], z$Sequence_2[1], z$Species_1[1], z$Species_2[1]))
      mk <- .bc_markers[[.bc_nz(input$marker, "other")]]
      s <- r$summary
      if (identical(lvl, "species") && s$n_testable >= 3L && s$n_with_gap / s$n_testable < 0.5)
        notes <- c(notes, sprintf("Only %d of %d species have a barcode gap.%s", s$n_with_gap, s$n_testable,
          if (identical(mk$level, "genus")) sprintf(" This is expected for %s: it is usually too conserved to separate species. Set 'Compare at the level of' to Genus.", mk$label)
          else " Close species may share this marker: check the Identification and Tree tabs, or use a more variable marker (ITS, matK, a multi-locus set)."))
      set_msg(if (length(notes)) paste(notes, collapse = " ") else NULL, "warning")
    })

    output$msg <- shiny::renderUI({
      if (is.null(rv$msg)) return(NULL)
      shiny::div(class = paste0("alert alert-", rv$msg$type), rv$msg$txt)
    })
    res <- shiny::reactive({ r <- rv$res; shiny::req(r); r })
    cats_now <- shiny::reactive(if (is.null(rv$res$cats)) names(.bc_cat_cols) else rv$res$cats)
    # the columns that hold the compared groups are called Species / Genus, or Genus / Family at the genus level
    lvl_cols <- function(r, d) {
      L <- .bc_levels[[.bc_nz(r$level, "species")]]
      if (identical(L$ucol, "Species")) return(d)
      names(d)[names(d) == "Genus"] <- "Family"; names(d)[names(d) == "Species"] <- "Genus"; d
    }
    f3 <- function(x) if (is.null(x) || is.na(x)) "-" else formatC(round(x, 4) + 0, format = "f", digits = 4)

    output$mega_head <- shiny::renderUI({
      m <- rv$mega
      if (is.null(m)) return(shiny::helpText("Choose the analysis, set the preferences and press Compute."))
      shiny::tagList(shiny::h5(m$head), shiny::helpText(m$note))
    })
    output$tbl_mega <- DT::renderDT({ m <- rv$mega; shiny::req(m); dtab(m$tab, "ftp", 15) })
    output$dl_mega <- shiny::downloadHandler(function() paste0("distances_", .bc_nz(rv$mega$an, "result"), ".csv"),
      function(file) { shiny::req(rv$mega); utils::write.csv(rv$mega$tab, file, row.names = FALSE) })

    output$boxes <- shiny::renderUI({
      shiny::req(is.null(input$analysis) || identical(input$analysis, "barcode"))
      r <- res(); s <- r$summary
      verdict <- if (is.na(s$gap)) "-" else if (s$gap > 0) "Yes" else "No"
      cm <- s$cat_mean
      shiny::tagList(
        shiny::fluidRow(
          shiny::column(3, .infobox("layer-group", if (identical(r$level, "genus")) "Genera | families | outgroup" else "Species | genera | outgroup",
                                    sprintf("%d | %d | %d seq.", s$n_species, s$n_genera, s$n_out), "#2c7be5")),
          shiny::column(3, .infobox("arrows-left-right", if (identical(r$level, "genus")) "Mean within | same family | other" else "Mean intra | congen. | intergen.",
                                    paste(f3(cm[1]), "|", f3(cm[2]), "|", f3(cm[3])), "#6f42c1")),
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
      mk <- .bc_markers[[.bc_nz(input$marker, "other")]]
      .bc_plot_dens(r$pairs, mode, .bc_nz(input$dens_adj, 1), r$model,
                    if (isTRUE(input$dens_thr) && !is.null(r$thr)) r$thr$best else NULL, .bc_cat_pal(r$cats, fig_cols("dens")), fig_fs("dens"),
                    if (isTRUE(input$ref_thr) && length(mk$thr)) mk$thr else NULL)
    })
    plot_box <- shiny::reactive({
      r <- res()
      .bc_plot_box(r$pairs, if (!is.null(r$tests$kw)) r$tests$kw$P else NULL, identical(input$box_style, "violin"),
                   isTRUE(input$box_pts), r$model, .bc_cat_pal(r$cats, fig_cols("box")), fig_fs("box"))
    })
    plot_pcoa <- shiny::reactive({
      r <- res()
      shiny::validate(shiny::need(!is.null(r$pcoa), "The PCoA needs at least 3 sequences with a complete distance matrix."))
      ax <- as.integer(strsplit(.bc_nz(input$pc_axes, "1,2"), ",")[[1]])
      k <- r$keep
      use <- if (isTRUE(input$pc_out)) rep(TRUE, length(k)) else r$role[k] != "Outgroup"
      shiny::validate(shiny::need(sum(use) >= 3L, "Fewer than 3 ingroup sequences with a complete matrix: tick 'Include the outgroup in the PCoA'."))
      pc <- .bc_pcoa(r$Dc, use)
      shiny::validate(shiny::need(!is.null(pc), "The PCoA could not be computed."))
      .bc_plot_pcoa(pc, r$sp[k][use], r$ge[k][use], r$role[k][use], .bc_nz(input$pc_col, "species"), ax, isTRUE(input$pc_ell), isTRUE(input$pc_lab),
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
      .bc_plot_scatter(d, c("barcode gap" = v[["gap"]], "overlap" = v[["overlap"]]), fig_fs("scat"), .bc_levels[[.bc_nz(r$level, "species")]]$unit)
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
      k <- r$keep
      use <- if (isTRUE(input$pc_out)) rep(TRUE, length(k)) else r$role[k] != "Outgroup"
      pc <- if (sum(use) >= 3L) .bc_pcoa(r$Dc, use) else NULL; shiny::req(pc)
      shiny::helpText(sprintf("PCoA of the %s distance matrix (%d sequences%s). The percentages are the share of the positive eigenvalues. ",
                              r$model, sum(use), if (isTRUE(input$pc_out) || !any(r$role[k] == "Outgroup")) "" else ", outgroup left out"),
                      if (pc$neg_pct > 5) sprintf("Negative eigenvalues make up %.1f%% of the total: the matrix is not fully Euclidean, so read the axes with some care. ", pc$neg_pct) else "",
                      if (isTRUE(input$pc_out) && any(r$role[k] == "Outgroup")) "A distant outgroup takes most of the first axis and squeezes the ingroup together: untick the option to see the ingroup alone." else "")
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
    output$fb_dens <- shiny::renderUI(.bc_figbar(session$ns, "dens", cats_now()))
    output$fb_box <- shiny::renderUI(.bc_figbar(session$ns, "box", cats_now()))
    gd_now <- function() { r <- res(); shiny::req(r$gd); r$gd[[if (identical(input$gd_level, "parent")) "parent" else "unit"]] }
    output$tbl_gd_in <- DT::renderDT({ g <- gd_now(); dtab(g$within, "t") })
    output$tbl_gd_bt <- DT::renderDT({ g <- gd_now(); dtab(g$between, "ftp", 12) })
    output$tbl_species <- DT::renderDT({ r <- res(); dtab(lvl_cols(r, r$species)) })
    output$tbl_pairs <- DT::renderDT({
      r <- res(); d <- r$pairs[, setdiff(names(r$pairs), c("Genus_1", "Genus_2", "Role_1", "Role_2", "Use"))]
      if (identical(r$level, "genus")) names(d) <- sub("^Species_", "Genus_", names(d))
      dtab(d)
    })
    output$tbl_mat <- DT::renderDT({
      r <- res()
      M <- if (identical(input$heat_level, "sp")) .bc_spmatrix(r$D, r$sp, r$role, TRUE) else r$D
      dtab(data.frame(Name = rownames(M), M, check.names = FALSE, stringsAsFactors = FALSE), n = 10)
    })
    ident <- shiny::reactive({ r <- res(); .bc_ident(r$D, r$sp, r$role, thr_eff(r), r$targets) })
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
                     Group_within = r$gd$unit$within, Group_between = r$gd$unit$between,
                     Parent_within = r$gd$parent$within, Parent_between = r$gd$parent$between,
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
    get_state <- function() shiny::isolate(list(groups = rv$groups, res = rv$res, targets = rv$targets))
    set_state <- function(s) {
      g <- s$groups
      if (!is.null(g) && !("Species" %in% names(g)) && "Group" %in% names(g)) {   # project saved by the previous version
        g <- data.frame(Sequence = g$Sequence, Species = g$Group, Genus = .bc_group(g$Sequence, "genus"),
                        Role = "Ingroup", stringsAsFactors = FALSE)
      }
      if (!is.null(g) && !("Family" %in% names(g))) g$Family <- ""
      if (!is.null(g)) g <- g[, c("Sequence", "Species", "Genus", "Family", "Role")]
      rv$groups <- g
      rv$targets <- if (!is.null(s$targets)) as.character(s$targets) else character()
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
