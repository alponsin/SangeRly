# Sangerly - "Project and report" tab, "About" tab
#   * marker (gene) presets: primers, product, genetic code
#   * save / open a whole project (one .rds file)
#   * synthetic demo data
#   * overall report (HTML) with an automatically written Methods paragraph and the references

SANGERLY_VERSION <- "1.0.0"
# Software citation (GitHub). Change these two lines to the real author names and repository address.
SANGERLY_AUTHORS <- "Alponsin & Putra, A. E."
SANGERLY_URL <- "https://github.com/alponsin/sangerly"
.sangerly_citation <- function()
  sprintf("%s (%s). Sangerly: R Shiny application for Sanger sequence assembly, BLAST identification, alignment and phylogeny (version %s) [Computer software]. %s",
          SANGERLY_AUTHORS, format(Sys.Date(), "%Y"), SANGERLY_VERSION, SANGERLY_URL)

# Typical markers. The primer sequences are the commonly published ones: compare them with the primers you really used.
GENE_PRESETS <- list(
  COI = list(label = "COI / cox1 (animal DNA barcode)", gene = "COI", feature = "CDS", product = "cytochrome c oxidase subunit I",
             location = "mitochondrion", code = "5", size = "about 650 bp (barcode region)",
             primers = c("LCO1490 (forward)" = "GGTCAACAAATCATAAAGATATTGG", "HCO2198 (reverse)" = "TAAACTTCAGGGTGACCAAAAAATCA"),
             note = "Animals. Use genetic code 2 for vertebrates, 5 for most invertebrates."),
  ITS = list(label = "ITS (nuclear rDNA; fungi and plants)", gene = "", feature = "misc_RNA",
             product = "contains internal transcribed spacer 1, 5.8S ribosomal RNA gene, and internal transcribed spacer 2",
             location = "", code = "1", size = "about 500-800 bp",
             primers = c("ITS1 (forward)" = "TCCGTAGGTGAACCTGCGG", "ITS4 (reverse)" = "TCCTCCGCTTATTGATATGC"),
             note = "Fungi (and many plants). Nuclear, so no genome location."),
  `16S` = list(label = "16S rRNA (bacteria)", gene = "", feature = "rRNA", product = "16S ribosomal RNA",
               location = "", code = "11", size = "about 1,400 bp (27F-1492R)",
               primers = c("27F (forward)" = "AGAGTTTGATCMTGGCTCAG", "1492R (reverse)" = "GGTTACCTTGTTACGACTT"),
               note = "Bacteria and archaea. For animal mitochondrial 16S choose the genome 'mitochondrion' in the GenBank tab."),
  rbcL = list(label = "rbcL (plant plastid barcode)", gene = "rbcL", feature = "CDS",
              product = "ribulose-1,5-bisphosphate carboxylase/oxygenase large subunit",
              location = "chloroplast", code = "11", size = "about 550-600 bp (rbcLa)",
              primers = c("rbcLa-F (forward)" = "ATGTCACCACAAACAGAGACTAAAGC", "rbcLa-R (reverse)" = "GTAAAATCAAGTCCACCRCG"),
              note = "Land plants."),
  cytb = list(label = "cytochrome b (vertebrate mitochondrial)", gene = "CYTB", feature = "CDS", product = "cytochrome b",
              location = "mitochondrion", code = "2", size = "about 1,000-1,140 bp",
              primers = c("L14841 (forward)" = "AAAAAGCTTCCATCCAACATCTCAGCATGATGAAA", "H15915 (reverse)" = "AACTGCAGTCATCTCCGGTTTACAAGAC"),
              note = "Mammals, birds, reptiles, fish (genetic code 2).")
)
.marker_choices <- c("(not specified)" = "none", stats::setNames(names(GENE_PRESETS), vapply(GENE_PRESETS, function(x) x$label, "")))

# References that the Methods paragraph can cite
.REFS <- c(
  blast1 = "Altschul SF, Gish W, Miller W, Myers EW, Lipman DJ (1990) Basic local alignment search tool. Journal of Molecular Biology 215: 403-410.",
  blast2 = "Camacho C, Coulouris G, Avagyan V, Ma N, Papadopoulos J, Bealer K, Madden TL (2009) BLAST+: architecture and applications. BMC Bioinformatics 10: 421.",
  decipher = "Wright ES (2016) Using DECIPHER v2.0 to analyze big biological sequence data in R. The R Journal 8(1): 352-359.",
  msa = "Bodenhofer U, Bonatesta E, Horejs-Kainrath C, Hochreiter S (2015) msa: an R package for multiple sequence alignment. Bioinformatics 31: 3997-3999.",
  trimal = "Capella-Gutierrez S, Silla-Martinez JM, Gabaldon T (2009) trimAl: a tool for automated alignment trimming in large-scale phylogenetic analyses. Bioinformatics 25: 1972-1973.",
  ape = "Paradis E, Schliep K (2019) ape 5.0: an environment for modern phylogenetics and evolutionary analyses in R. Bioinformatics 35: 526-528.",
  nj = "Saitou N, Nei M (1987) The neighbor-joining method: a new method for reconstructing phylogenetic trees. Molecular Biology and Evolution 4: 406-425.",
  k2p = "Kimura M (1980) A simple method for estimating evolutionary rates of base substitutions through comparative studies of nucleotide sequences. Journal of Molecular Evolution 16: 111-120.",
  iqtree3 = "Wong TKF, Ly-Trong N, Ren H, Demotte P, Banos H, Roger AJ, Susko E, Bielow C, De Maio N, Goldman N, Hahn MW, dos Reis M, Vinh LS, Huttley G, Lanfear R, Minh BQ (2026) IQ-TREE 3: phylogenomic inference software using complex evolutionary models. Molecular Biology and Evolution msag117. https://doi.org/10.1093/molbev/msag117",
  iqtree = "Minh BQ, Schmidt HA, Chernomor O, Schrempf D, Woodhams MD, von Haeseler A, Lanfear R (2020) IQ-TREE 2: new models and efficient methods for phylogenetic inference in the genomic era. Molecular Biology and Evolution 37: 1530-1534.",
  mf = "Kalyaanamoorthy S, Minh BQ, Wong TKF, von Haeseler A, Jermiin LS (2017) ModelFinder: fast model selection for accurate phylogenetic estimates. Nature Methods 14: 587-589.",
  ufboot = "Hoang DT, Chernomor O, von Haeseler A, Minh BQ, Vinh LS (2018) UFBoot2: improving the ultrafast bootstrap approximation. Molecular Biology and Evolution 35: 518-522.",
  alrt = "Guindon S, Dufayard JF, Lefort V, Anisimova M, Hordijk W, Gascuel O (2010) New algorithms and methods to estimate maximum-likelihood phylogenies: assessing the performance of PhyML 3.0. Systematic Biology 59: 307-321.",
  tvbot = "Xie J, Chen Y, Cai G, Cai R, Hu Z, Wang H (2023) Tree Visualization By One Table (tvBOT): a web application for visualizing, modifying and annotating phylogenetic trees. Nucleic Acids Research 51: W587-W592. https://doi.org/10.1093/nar/gkad359",
  itol = "Letunic I, Bork P (2024) Interactive Tree of Life (iTOL) v6: recent updates to the phylogenetic tree display and annotation tool. Nucleic Acids Research 52: W78-W82.")

# ---- synthetic demo ------------------------------------------------------------------------------------------
# Three simulated "species" (about 8% from a common ancestor) with five sequences each (0.6% variation inside a species).
# NOT real biological data: only to try the program.
.demo_records <- function() {
  old <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
  on.exit(if (!is.null(old)) assign(".Random.seed", old, envir = globalenv()), add = TRUE)
  set.seed(42)
  b <- c("A", "C", "G", "T")
  mutate <- function(s, rate) {
    ch <- strsplit(s, "")[[1]]
    for (i in which(stats::runif(length(ch)) < rate)) ch[i] <- sample(setdiff(b, ch[i]), 1)
    paste(ch, collapse = "")
  }
  root <- paste(sample(b, 600, TRUE), collapse = "")
  sp <- c(speciesA = mutate(root, 0.08), speciesB = mutate(root, 0.08), speciesC = mutate(root, 0.08))
  do.call(rbind, lapply(names(sp), function(n)
    data.frame(id = sprintf("Demo_%s_%02d", n, 1:5), seq = vapply(1:5, function(i) mutate(sp[[n]], 0.006), ""),
               stringsAsFactors = FALSE)))
}

# ---- Methods paragraph --------------------------------------------------------------------------------------------
.nj_model_name <- function(m) {
  switch(as.character(m), K80 = "the Kimura two-parameter model (K2P)", JC69 = "the Jukes-Cantor model", F81 = "the F81 model",
         K81 = "the Kimura three-parameter model", F84 = "the F84 model", TN93 = "the Tamura-Nei model", raw = "uncorrected p-distances",
         logdet = "LogDet distances", as.character(m))
}

# Returns list(text = paragraphs separated by a blank line, refs = keys of .REFS that were cited)
.methods_build <- function(marker, asm_aln, consensus, blast_p, blast_res, msa_meta, aln_all, aln_trim, nj, iq, bc) {
  parts <- character(); refs <- character()
  add <- function(x) parts <<- c(parts, x)
  cite <- function(...) refs <<- c(refs, c(...))
  mk <- if (!identical(marker, "none")) GENE_PRESETS[[marker]] else NULL

  if (!is.null(mk))
    add(sprintf("The target region (%s) was amplified with the primers %s (%s) and %s (%s).", mk$label,
                names(mk$primers)[1], mk$primers[[1]], names(mk$primers)[2], mk$primers[[2]]))

  if (!is.null(asm_aln) && isTRUE(asm_aln$ok)) {
    P <- asm_aln$P; n <- length(asm_aln$labels)
    trim_txt <- if (identical(P$trim_method, "M1")) sprintf("the modified Mott algorithm (error cutoff %g)", P$mott_cutoff)
                else sprintf("a sliding window (%d bases, mean Phred quality of at least %g)", as.integer(P$win), P$q_trim)
    how <- if (identical(asm_aln$mode, "ref")) "mapped to a reference sequence by local alignment"
           else "overlapped (forward read against the reverse complement of the reverse read)"
    add(sprintf(paste0("Sanger chromatograms (%d read%s) were quality-trimmed with %s and %s. The consensus was called column by column in Sangerly ",
                       "(version %s) from the Phred quality (good >= Q%g, poor < Q%g) and the secondary-peak ratio (>= %g, double peak); ",
                       "conflicting positions were inspected on the chromatograms and edited manually where necessary."),
                n, if (n == 1) "" else "s", trim_txt, how, SANGERLY_VERSION, P$q_good, P$q_poor, P$ratio))
  }

  if (!is.null(blast_res) && nrow(blast_res) && !is.null(blast_p)) {
    cite("blast1", "blast2")
    add(sprintf("Sequence identity was checked with %s against the NCBI %s database%s (E-value <= %g, up to %s target sequences; Altschul et al. 1990; Camacho et al. 2009).",
                toupper(.nz(blast_p$program, "blastn")), .nz(blast_p$db, "nucleotide"),
                if (!is.null(blast_p$date) && !is.na(blast_p$date)) paste0(", accessed ", blast_p$date) else "",
                .nz(blast_p$evalue, 0.05), .nz(blast_p$hitlist, "50")))
  }

  if (!is.null(aln_all)) {
    al <- msa_meta$align; tr <- msa_meta$trim
    nc <- nchar(aln_all[[1]]); ns <- length(aln_all)
    if (isTRUE(al$prealigned)) s1 <- sprintf("An alignment of %d sequences supplied by the user was used", ns)
    else {
      lab <- .nz(al$label, "")
      if (grepl("DECIPHER", lab)) { cite("decipher"); prog <- "DECIPHER"; cit <- "Wright 2016" }
      else if (nzchar(lab)) { cite("msa"); prog <- paste(sub(" \\(msa\\)$", "", lab), "(R package msa)"); cit <- "Bodenhofer et al. 2015" }
      else { prog <- "a multiple sequence aligner"; cit <- "" }
      s1 <- sprintf("%d sequences were aligned with %s%s%s", ns, prog, if (nzchar(cit)) paste0("; ", cit) else "",
                    if (isTRUE(al$orient)) ", reverse-complementing the sequences that fitted better in the opposite orientation" else "")
    }
    s2 <- ""
    if (!is.null(aln_trim)) {
      nk <- nchar(aln_trim[[1]])
      eng <- .nz(tr$engine, "native")
      if (identical(eng, "trimal")) {
        cite("trimal")
        fl <- .nz(tr$tflag, "-automated1")
        s2 <- sprintf(" The alignment (%d columns) was trimmed with trimAl (%s; Capella-Gutierrez et al. 2009), leaving %d columns.", nc,
                      if (identical(fl, "custom")) sprintf("-gt %g%s%s", .nz(tr$t_gt, 0.5), if (isTRUE(tr$t_st > 0)) sprintf(" -st %g", tr$t_st) else "",
                                                          if (isTRUE(tr$t_cons > 0)) sprintf(" -cons %g", tr$t_cons) else "") else sub("^-", "", fl), nk)
      } else {
        rule <- switch(.nz(tr$mode, "gt"),
          gt = sprintf("columns in which fewer than %g%% of the sequences had a base were removed (gap threshold)", 100 * .nz(tr$gt, 0.5)),
          st = sprintf("columns with a conservation below %g were removed", .nz(tr$st, 0.5)),
          both = sprintf("columns with fewer than %g%% bases or a conservation below %g were removed", 100 * .nz(tr$gt, 0.5), .nz(tr$st, 0.5)),
          ends = sprintf("the ragged ends were cut where fewer than %g%% of the sequences had a base", 100 * .nz(tr$ends_cov, 0.5)),
          none = "no columns were removed", "poorly aligned columns were removed")
        s2 <- sprintf(" The alignment (%d columns) was then trimmed in Sangerly with trimAl-style rules: %s, leaving %d columns.", nc, rule, nk)
      }
    }
    add(paste0(s1, ".", s2))
  }

  if (!is.null(nj) && !is.null(nj$tree)) {
    cite("ape", "nj"); if (identical(nj$model, "K80")) cite("k2p")
    add(sprintf("A neighbour-joining tree (Saitou and Nei 1987) was built with the R package ape (Paradis and Schliep 2019) from distances calculated with %s%s (pairwise deletion of gaps)%s.",
                .nj_model_name(nj$model), if (isTRUE(nj$alpha > 0)) sprintf(" with gamma-distributed rates (alpha = %.3f)", nj$alpha) else "",
                if (isTRUE(nj$B > 0)) sprintf("; node support was estimated with %d bootstrap replicates", as.integer(nj$B)) else ""))
  }
  if (!is.null(iq) && !is.null(iq$tree)) {
    v3 <- !is.null(iq$version) && !is.na(iq$version) && grepl("^3", iq$version)
    cite(if (v3) "iqtree3" else "iqtree"); o <- .nz(iq$opts, list())
    ver <- if (!is.null(iq$version) && !is.na(iq$version)) paste0(" version ", iq$version) else ""
    best <- iq$report$best
    mtxt <- if (isTRUE(o$model == "fixed")) sprintf("the %s substitution model", .nz(o$fixed, "")) else if (!is.na(best)) {
      cite("mf"); sprintf("the best-fit substitution model (%s) chosen by ModelFinder according to the Bayesian information criterion (Kalyaanamoorthy et al. 2017)", best)
    } else "the substitution model chosen by ModelFinder"
    sup <- if (isTRUE(o$boot == "ufboot")) { cite("ufboot"); sprintf("; branch support was assessed with %d ultrafast bootstrap replicates (Hoang et al. 2018)", as.integer(.nz(o$nboot, 1000))) }
           else if (isTRUE(o$boot == "std")) sprintf("; branch support was assessed with %d standard bootstrap replicates", as.integer(.nz(o$nboot, 100))) else ""
    al <- if (isTRUE(o$alrt)) { cite("alrt"); " and the SH-aLRT test (1000 replicates; Guindon et al. 2010)" } else ""
    add(sprintf("A maximum-likelihood tree was inferred with IQ-TREE%s (%s) under %s%s%s.", ver, if (v3) "Wong et al. 2026" else "Minh et al. 2020", mtxt, sup, al))
  }

  if (!is.null(bc)) {
    s <- bc$summary; cite("k2p")
    add(sprintf(paste0("Genetic distances within and between groups (%d groups, %d sequences) were calculated with %s and pairwise deletion. ",
                       "The barcode gap was assessed as the difference between the smallest interspecific and the largest intraspecific distance ",
                       "(mean intraspecific distance %.2f%%, mean interspecific distance %.2f%%; %d of %d testable groups showed a gap)."),
                s$n_groups, s$n_seq, if (identical(s$model, "K80")) "the Kimura two-parameter model (Kimura 1980)" else paste("the", s$model, "model"),
                100 * s$mean_intra, 100 * s$mean_inter, s$n_with_gap, s$n_testable))
  }
  list(text = if (length(parts)) paste(parts, collapse = "\n\n")
              else "Nothing to describe yet: run at least one step (assembly, BLAST, alignment, tree).",
       refs = unique(refs))
}

# ---- HTML report -----------------------------------------------------------------------------------------------------
.tree_img <- function(tr, title) {
  f <- tempfile(fileext = ".png"); on.exit(unlink(f), add = TRUE)
  n <- ape::Ntip(tr)
  grDevices::png(f, width = 1800, height = round(max(1100, 45 * n + 250)), res = 200)
  tryCatch(.plot_tree(tr, "phylogram", TRUE, 0, 0.9, TRUE, title), finally = grDevices::dev.off())
  sprintf("<img style='max-width:100%%' alt='%s' src='data:image/png;base64,%s'>", htmltools::htmlEscape(title),
          jsonlite::base64_enc(readBin(f, "raw", file.size(f))))
}

.report_html <- function(ctx) {
  esc <- htmltools::htmlEscape
  has <- function(k) k %in% ctx$sections
  safe <- function(expr) tryCatch(expr, error = function(e)
    sprintf("<p class='warn'>This section could not be built: %s</p>", esc(conditionMessage(e))))
  css <- paste("body{font-family:system-ui,Segoe UI,Arial,sans-serif;max-width:1100px;margin:24px auto;padding:0 16px;color:#222;line-height:1.5}",
               "h1{color:#0b2a5b}h2{border-bottom:2px solid #d6e4f0;padding-bottom:4px;margin-top:34px;color:#0b3d91}.meta{color:#6b7785}",
               "table{border-collapse:collapse;font-size:12.5px;margin:8px 0}th,td{border:1px solid #ccd;padding:3px 7px;text-align:left}th{background:#eef3f8}",
               "pre{background:#f6f8fa;padding:10px;overflow:auto;white-space:pre-wrap;word-break:break-all;font-size:12px}",
               ".warn{background:#fff6dd;padding:8px 12px;border-radius:4px}.refs li{margin-bottom:4px}",
               ".msa-box{overflow-x:auto;border:1px solid #dfe3e8;padding:8px}.msa-box pre{margin:0 0 8px;font-family:Consolas,monospace;font-size:11px;line-height:1.3;background:transparent;padding:0}",
               ".mA{background:#c6f0d6}.mC{background:#cfe4fb}.mG{background:#e3e3e3}.mT{background:#ffd0d0}.mN{background:#f4f0c8}.mgap{color:#9aa5b1}.ruler{color:#7b8794}", sep = "\n")
  h <- c(sprintf("<!doctype html><html><head><meta charset='utf-8'><title>%s</title><style>%s</style></head><body>", esc(ctx$title), css),
         sprintf("<h1>%s</h1>", esc(ctx$title)),
         sprintf("<p class='meta'>%sGenerated %s by Sangerly %s (%s).</p>", if (nzchar(ctx$author)) paste0(esc(ctx$author), ". ") else "",
                 format(Sys.time(), "%Y-%m-%d %H:%M"), SANGERLY_VERSION, esc(R.version.string)))
  if (has("methods"))
    h <- c(h, "<h2>Methods</h2>", safe(paste0("<p>", esc(strsplit(ctx$methods$text, "\n\n", fixed = TRUE)[[1]]), "</p>", collapse = "")))
  if (has("assembly") && !is.null(ctx$consensus))
    h <- c(h, "<h2>Assembly</h2>", safe({
      a <- ctx$asm_aln
      c(if (!is.null(a) && length(a$notes)) paste0("<p class='warn'>", paste(esc(a$notes), collapse = "<br>"), "</p>"),
        "<h3>Consensus</h3>", sprintf("<pre>%s</pre>", esc(ctx$consensus)))
    }))
  if (has("blast") && !is.null(ctx$blast_res) && nrow(ctx$blast_res))
    h <- c(h, "<h2>BLAST</h2>", safe({
      r <- ctx$blast_res
      top <- do.call(rbind, lapply(split(r, factor(r$Query, levels = unique(r$Query))), utils::head, ctx$top_n))
      keep <- intersect(c("Query", "Rank", "Species", "Accession", "Identity", "QueryCover", "Evalue", "BitScore"), names(top))
      top <- top[, keep, drop = FALSE]
      for (k in intersect(c("Identity", "QueryCover", "BitScore"), names(top))) top[[k]] <- round(top[[k]], 1)
      c(sprintf("<p>Database: %s. Program: %s. %s</p>", esc(.nz(ctx$blast_p$db, "-")), esc(.nz(ctx$blast_p$program, "-")),
                if (!is.null(ctx$blast_p$rid)) paste("RID:", esc(ctx$blast_p$rid)) else ""),
        sprintf("<p>The best %d hit(s) per query sequence:</p>", ctx$top_n), .html_table(top, 400))
    }))
  if (has("msa") && !is.null(ctx$aln_all))
    h <- c(h, "<h2>Alignment</h2>", safe({
      a <- if (!is.null(ctx$aln_trim)) ctx$aln_trim else ctx$aln_all
      c(sprintf("<p>%d sequences; whole alignment %d columns%s.</p>", length(ctx$aln_all), nchar(ctx$aln_all[[1]]),
                if (!is.null(ctx$aln_trim)) sprintf(", trimmed result %d columns (shown below)", nchar(ctx$aln_trim[[1]])) else ""),
        as.character(msa_html(msa_matrix(a), width = 80L, max_cols = 1600L)))
    }))
  if (has("tree") && (!is.null(ctx$nj$tree) || !is.null(ctx$iq$tree)))
    h <- c(h, "<h2>Phylogeny</h2>", safe({
      out <- character()
      if (!is.null(ctx$nj$tree)) out <- c(out, "<h3>Neighbour-joining</h3>", .tree_img(ctx$nj$tree, "Neighbour-joining"),
                                          sprintf("<p>Newick:</p><pre>%s</pre>", esc(ape::write.tree(ctx$nj$tree))))
      if (!is.null(ctx$iq$tree)) out <- c(out, sprintf("<h3>IQ-TREE%s</h3>", if (!is.na(ctx$iq$report$best)) paste0(" (", esc(ctx$iq$report$best), ")") else ""),
                                          .tree_img(ctx$iq$tree, "IQ-TREE"), sprintf("<p>Newick:</p><pre>%s</pre>", esc(ape::write.tree(ctx$iq$tree))))
      out
    }))
  if (has("barcode") && !is.null(ctx$bc))
    h <- c(h, "<h2>Barcode gap</h2>", safe({
      s <- ctx$bc$summary; g <- ctx$bc$groups
      for (k in c("Max_intra", "Mean_intra", "Min_inter", "Gap")) g[[k]] <- round(g[[k]], 4)
      c(sprintf("<p>Model %s. Mean intraspecific %.4f, mean interspecific %.4f; largest intraspecific %.4f, smallest interspecific %.4f; groups with a gap: %d of %d.</p>",
                esc(s$model), s$mean_intra, s$mean_inter, s$max_intra, s$min_inter, s$n_with_gap, s$n_testable), .html_table(g, 400))
    }))
  if (has("refs") && length(ctx$methods$refs))
    h <- c(h, "<h2>References</h2>", paste0("<ol class='refs'>", paste0("<li>", esc(sort(.REFS[ctx$methods$refs])), "</li>", collapse = ""), "</ol>"),
           "<p class='meta'>Check every reference against the original publication before you cite it.</p>")
  paste(c(h, "</body></html>"), collapse = "\n")
}

# ---- UI -----------------------------------------------------------------------------------------------------------------
projectUI <- function(id) {
  ns <- shiny::NS(id)
  shiny::div(class = "container-fluid",
    shiny::fluidRow(
      shiny::column(5,
        .card("Marker (gene)",
          shiny::selectInput(ns("marker"), "What did you sequence?", .marker_choices, "none"),
          shiny::uiOutput(ns("marker_info")),
          shiny::helpText("The marker fills the primers in the Methods paragraph and the feature in the GenBank tab. ",
                          "Copy the primers into the Assembly tab if you want to crop them from the consensus.")),
        .card("Project file",
          shiny::helpText("One .rds file keeps the BLAST results, the MSA list, alignment and trimming, the trees and the barcode-gap result. ",
                          "Raw chromatograms are not included: use the Assembly tab's own session file for those."),
          shiny::downloadButton(ns("dl_proj"), "Save project (.rds)"),
          shiny::br(), shiny::br(),
          shiny::fileInput(ns("open_proj"), "Open a project (.rds)", accept = ".rds", buttonLabel = "Browse...", placeholder = ""),
          shiny::helpText("Open only project files that you saved yourself: an .rds file from an unknown source is not safe."),
          shiny::uiOutput(ns("proj_msg"))),
        .card("Try the program",
          shiny::helpText("Loads 15 simulated sequences (3 'species' x 5) into the MSA tab. They are made up by the program and are NOT real data; ",
                          "use them to try alignment, trees and the barcode gap."),
          shiny::actionButton(ns("demo"), "Load the synthetic demo", icon = shiny::icon("flask")))),
      shiny::column(7,
        .card("Methods paragraph (written from what you did)",
          shiny::verbatimTextOutput(ns("methods"), placeholder = TRUE),
          shiny::tags$style(shiny::HTML(sprintf("#%s{white-space:pre-wrap;word-break:normal}", ns("methods")))),
          shiny::div(style = "display:flex;gap:6px;flex-wrap:wrap;margin-top:8px",
            shiny::downloadButton(ns("dl_methods"), "Methods + references (.txt)")),
          shiny::helpText("The text follows the settings that were actually used. Read it, correct it and adapt it to your journal; ",
                          "details that the program cannot know (reagents, PCR conditions, vouchers) must be added by you.")),
        .card("Overall report (HTML)",
          shiny::textInput(ns("title"), "Title", "Sangerly report", width = "100%"),
          shiny::textInput(ns("author"), "Author / sample set (optional)", "", width = "100%"),
          shiny::checkboxGroupInput(ns("sections"), "Sections",
            c("Methods paragraph" = "methods", "Assembly (consensus)" = "assembly", "BLAST (best hits)" = "blast", "Alignment" = "msa",
              "Trees" = "tree", "Barcode gap" = "barcode", "References" = "refs"),
            selected = c("methods", "assembly", "blast", "msa", "tree", "barcode", "refs")),
          shiny::numericInput(ns("top_n"), "Hits per query shown in the BLAST table", 5, min = 1, max = 50, step = 1, width = "200px"),
          shiny::downloadButton(ns("dl_report"), "Download the report (HTML)", class = "btn-primary"),
          shiny::helpText("Sections without results are left out. The file is self-contained: open it in a browser, or print it to PDF.")))))
}

aboutUI <- function(id) {
  ns <- shiny::NS(id)
  i <- function(en, idn) shiny::tags$span(class = "i18n", `data-en` = en, `data-id` = idn, en)
  shiny::div(class = "container", style = "max-width:900px",
    shiny::div(class = "sg-hero", shiny::h1(shiny::icon("dna"), " Sangerly"), shiny::p(paste("Version", SANGERLY_VERSION))),
    .card(i("About", "Tentang"),
      shiny::p(i("Sangerly is an R Shiny application for Sanger sequencing: assembly of forward and reverse reads, species identification with BLAST, multiple sequence alignment and trimming, phylogenetic trees (neighbour-joining and IQ-TREE), barcode-gap analysis and preparation of files for GenBank.",
                 "Sangerly adalah aplikasi R Shiny untuk sekuensing Sanger: perakitan read forward dan reverse, identifikasi spesies dengan BLAST, multiple sequence alignment dan trimming, pohon filogenetik (neighbour-joining dan IQ-TREE), analisis barcode gap, dan penyiapan berkas untuk GenBank.")),
      shiny::p(i("The program runs on your own computer; only BLAST, NCBI Taxonomy and the tree websites that you open yourself use the internet.",
                 "Program berjalan di komputer Anda sendiri; hanya BLAST, NCBI Taxonomy, dan situs pohon yang Anda buka sendiri yang memakai internet."))),
    .card(i("How to cite", "Cara mengutip"),
      shiny::p(i("Until a paper is published, cite Sangerly as software from its GitHub repository:",
                 "Sebelum artikelnya terbit, kutip Sangerly sebagai perangkat lunak dari repositori GitHub-nya:")),
      shiny::tags$blockquote(style = "font-size:14px", .sangerly_citation()),
      shiny::helpText(i("Edit SANGERLY_AUTHORS and SANGERLY_URL at the top of R/project_module.R (and CITATION.cff in the repository). Also cite the programs that did the work (list below, and the References of your report).",
                        "Ubah SANGERLY_AUTHORS dan SANGERLY_URL di bagian atas R/project_module.R (dan CITATION.cff di repositori). Kutip juga program yang mengerjakan analisisnya (daftar di bawah dan References pada laporan Anda)."))),
    .card(i("Programs and packages", "Program dan paket"), shiny::tableOutput(ns("versions"))),
    .card(i("References of the methods", "Pustaka metode"),
      shiny::tags$ol(lapply(.REFS[c("blast1", "blast2", "decipher", "msa", "trimal", "ape", "nj", "k2p", "iqtree3", "iqtree", "mf", "ufboot", "alrt", "itol", "tvbot")],
                            function(x) shiny::tags$li(x))),
      shiny::helpText("Check each reference against the original publication before you cite it.")))
}

# ---- server -------------------------------------------------------------------------------------------------------------------
aboutServer <- function(id) {
  shiny::moduleServer(id, function(input, output, session) {
    output$versions <- shiny::renderTable({
      pk <- c("shiny", "ape", "Biostrings", "sangerseqR", "pwalign", "DECIPHER", "msa", "phangorn", "processx", "DT", "ggplot2")
      v <- vapply(pk, function(p) if (nzchar(system.file(package = p))) as.character(utils::packageVersion(p)) else "not installed", "")
      prog <- function(label, path, arg) {
        if (!nzchar(path)) return(c(label, "not found"))
        r <- tryCatch(suppressWarnings(system2(path, arg, stdout = TRUE, stderr = TRUE)), error = function(e) "could not be started")
        c(label, paste(utils::head(r[nzchar(r)], 1), collapse = " "))
      }
      pr <- rbind(prog("IQ-TREE", .iq_find(""), "--version"), prog("trimAl", .msa_find_trimal(""), "--version"))
      data.frame(Component = c("R", pk, pr[, 1]), Version = c(as.character(getRversion()), v, pr[, 2]), stringsAsFactors = FALSE)
    })
  })
}

projectServer <- function(id, asm, blast, msa, trees, barcode, parent_session) {
  shiny::moduleServer(id, function(input, output, session) {
    rv <- shiny::reactiveValues(msg = NULL)
    marker <- shiny::reactive(.nz(input$marker, "none"))

    output$marker_info <- shiny::renderUI({
      m <- GENE_PRESETS[[marker()]]
      if (is.null(m)) return(NULL)
      shiny::tagList(
        shiny::tags$table(class = "table table-sm",
          shiny::tags$tr(shiny::tags$th("Size"), shiny::tags$td(m$size)),
          shiny::tags$tr(shiny::tags$th("Feature"), shiny::tags$td(paste(m$feature, "-", m$product))),
          shiny::tags$tr(shiny::tags$th("Genetic code"), shiny::tags$td(m$code)),
          lapply(seq_along(m$primers), function(k)
            shiny::tags$tr(shiny::tags$th(names(m$primers)[k]), shiny::tags$td(shiny::tags$code(m$primers[[k]]))))),
        shiny::helpText(m$note, " Check the primers against the ones you really used."))
    })

    # ---- save / open ------------------------------------------------------------------------------------------
    collect <- function() list(
      app = "Sangerly", version = SANGERLY_VERSION, saved = Sys.time(), r = R.version.string, marker = marker(),
      consensus = shiny::isolate(asm$consensus()), blast = blast$get_state(), msa = msa$get_state(),
      tree = trees$get_state(), barcode = barcode$get_state())
    output$dl_proj <- shiny::downloadHandler(
      filename = function() paste0("sangerly_project_", format(Sys.time(), "%Y%m%d_%H%M"), ".rds"),
      content = function(file) saveRDS(collect(), file))

    shiny::observeEvent(input$open_proj, {
      s <- tryCatch(readRDS(input$open_proj$datapath), error = function(e) NULL)
      if (!is.list(s) || !identical(s$app, "Sangerly")) {
        rv$msg <- list(txt = "That file is not a Sangerly project.", type = "danger"); return()
      }
      res <- tryCatch({
        if (!is.null(s$blast)) blast$set_state(s$blast)
        if (!is.null(s$msa)) msa$set_state(s$msa)
        if (!is.null(s$tree)) trees$set_state(s$tree)
        if (!is.null(s$barcode)) barcode$set_state(s$barcode)
        if (!is.null(s$marker)) shiny::updateSelectInput(session, "marker", selected = s$marker)
        NULL
      }, error = function(e) conditionMessage(e))
      if (!is.null(res)) { rv$msg <- list(txt = paste("The project could not be restored completely:", res), type = "danger"); return() }
      rv$msg <- list(txt = sprintf("Project opened (saved %s with Sangerly %s).", format(s$saved, "%Y-%m-%d %H:%M"), .nz(s$version, "?")), type = "success")
      shiny::updateNavbarPage(parent_session, "nav",
        selected = if (!is.null(s$tree) && (!is.null(s$tree$nj) || !is.null(s$tree$iq))) "tree" else if (!is.null(s$msa$aln)) "msa" else "blast")
    })
    output$proj_msg <- shiny::renderUI({
      if (is.null(rv$msg)) NULL else shiny::div(class = paste0("alert alert-", rv$msg$type), rv$msg$txt)
    })

    # ---- demo ---------------------------------------------------------------------------------------------------
    shiny::observeEvent(input$demo, {
      msa$clear()
      n <- msa$add(.demo_records(), "Demo")
      shiny::showNotification(sprintf("%d simulated sequences were added to the MSA tab. Press Align there.", nrow(msa$pool())),
                              type = "message", duration = 10)
      shiny::updateNavbarPage(parent_session, "nav", selected = "msa")
    })

    # ---- Methods paragraph and report -------------------------------------------------------------------------
    methods <- shiny::reactive({
      .methods_build(marker(), asm$aln(), NULL, blast$params(), blast$results(), msa$meta(), msa$aln(), msa$trimmed(),
                     trees$nj(), trees$iq(), barcode$result())
    })
    output$methods <- shiny::renderText(methods()$text)
    output$dl_methods <- shiny::downloadHandler(function() "methods.txt", function(file) {
      m <- methods()
      writeLines(c(m$text, "", if (length(m$refs)) c("References", "", .REFS[sort(m$refs)])), file, useBytes = TRUE)
    })
    output$dl_report <- shiny::downloadHandler(
      filename = function() paste0("sangerly_report_", format(Sys.time(), "%Y%m%d_%H%M"), ".html"),
      content = function(file) {
        cons <- asm$consensus()
        cons_txt <- if (is.null(cons)) NULL else if (length(cons) == 1 && grepl("^>", cons[[1]])) cons[[1]]
                    else paste0(">", .nz(names(cons), paste0("seq", seq_along(cons))), "\n", cons, collapse = "\n")
        ctx <- list(title = .nz(input$title, "Sangerly report"), author = trimws(.nz(input$author, "")),
                    sections = input$sections, top_n = as.integer(.nz(input$top_n, 5)),
                    methods = methods(), asm_aln = asm$aln(), consensus = cons_txt,
                    blast_p = blast$params(), blast_res = blast$results(),
                    aln_all = msa$aln(), aln_trim = msa$trimmed(), nj = trees$nj(), iq = trees$iq(), bc = barcode$result())
        writeLines(.report_html(ctx), file, useBytes = TRUE)
      })

    list(marker = marker)
  })
}
