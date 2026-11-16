# SangeRly - "Phylogeny" tab
#   * Neighbour-joining (package ape): distance model chosen by you, or "auto" = the best-fit model found
#     by IQ-TREE ModelFinder (mapped to the closest distance model that ape offers), optional bootstrap.
#   * IQ-TREE (external program): ModelFinder (-m MFP / MF) or a fixed model, ultrafast / standard
#     bootstrap, SH-aLRT, outgroup, threads. It runs in the background with a live log and a Cancel button
#     (needs the 'processx' package; without it the run blocks the app until it ends).
#
# Packages:  install.packages(c("ape", "processx"))      optional: phangorn (midpoint rooting)
# Program :  iqtree2 / iqtree3 (https://iqtree.github.io): put it in the app folder (IQTREE/...), on the PATH,
#            or give its full path in the tab.

# ---------------------------------------------------------------------------------------------------
# Core functions (no Shiny)
# ---------------------------------------------------------------------------------------------------

# IQ-TREE: path typed by the user, else a copy inside the app folder (IQTREE/, iqtree/, any folder named iqtree*,
# searched recursively, newest version first), else the PATH.
.iq_find <- function(user = "") {
  user <- trimws(user)
  if (nzchar(user)) return(if (file.exists(user)) normalizePath(user, winslash = "/") else "")
  pat <- "^iqtree[0-9]*(\\.exe)?$"
  # also any folder of the app whose name starts with "iqtree" (for example the unzipped iqtree-3.0.1-Windows)
  extra <- list.dirs(".", recursive = FALSE, full.names = FALSE)
  extra <- extra[grepl("^iqtree", extra, ignore.case = TRUE)]
  for (d in unique(c("IQTREE", "iqtree", "IQ-TREE", "bin", extra))) {
    if (!dir.exists(d)) next
    f <- list.files(d, pattern = pat, recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
    f <- f[!dir.exists(f)]
    if (length(f)) return(normalizePath(sort(f, decreasing = TRUE)[1], winslash = "/"))
  }
  for (nm in c("iqtree3", "iqtree2", "iqtree")) {
    w <- unname(Sys.which(nm))
    if (nzchar(w)) return(w)
  }
  ""
}

# Version number of IQ-TREE from the first lines of its log ("IQ-TREE multicore version 3.0.1 ..."); NA if absent.
.iq_version <- function(log) {
  h <- grep("IQ-TREE.*version", log, value = TRUE)
  if (!length(h)) return(NA_character_)
  v <- regmatches(h[1], regexpr("version [0-9][0-9.]*", h[1]))
  if (!length(v)) NA_character_ else sub("version ", "", v)
}

# Clean sequence names (the tree programs dislike spaces and symbols) and drop sequences without any base.
.tree_prepare <- function(aln) {
  if (is.null(aln) || length(aln) < 3) return(NULL)
  aln <- gsub("U", "T", toupper(aln))
  nm <- names(aln)
  safe <- make.unique(gsub("[^A-Za-z0-9._-]", "_", nm), sep = "_")
  names(aln) <- safe
  nobase <- nchar(gsub("[-?.NX]", "", aln)) == 0
  dropped <- safe[nobase]
  aln <- aln[!nobase]
  if (length(aln) < 3) return(NULL)
  list(aln = aln, dropped = dropped, renamed = sum(safe != nm))
}

.tree_matrix <- function(aln) {
  m <- do.call(rbind, strsplit(tolower(unname(aln)), ""))
  rownames(m) <- names(aln)
  m
}

# Map an IQ-TREE model name (e.g. "GTR+F+I+G4") to the closest distance model of ape::dist.dna.
.iq_to_ape <- function(best) {
  base <- toupper(sub("\\+.*$", "", trimws(best)))
  model <- switch(base,
    "JC" = , "JC69" = "JC69",
    "F81" = "F81",
    "K80" = , "K2P" = "K80",
    "K81" = , "K3P" = "K81",
    "TN" = , "TN93" = , "HKY" = "TN93",
    "TN93")                                           # TPM*, TIM*, TVM*, SYM, GTR...: nearest model available in ape
  list(model = model, exact = base %in% c("JC", "JC69", "F81", "K80", "K2P", "K81", "K3P", "TN", "TN93"))
}

# Read the IQ-TREE report (.iqtree): best models, gamma alpha, log-likelihood, ModelFinder table.
.iq_parse_report <- function(path) {
  out <- list(best = NA_character_, aic = NA_character_, aicc = NA_character_, bic = NA_character_,
              alpha = NA_real_, loglik = NA_real_, table = NULL)
  if (!file.exists(path)) return(out)
  ln <- readLines(path, warn = FALSE)
  pick <- function(rx) {
    h <- grep(rx, ln, value = TRUE)
    if (length(h)) trimws(sub(rx, "", h[1])) else NA_character_
  }
  out$bic  <- pick("^Bayesian Information Criterion: *")
  out$aic  <- pick("^Akaike Information Criterion: *")
  out$aicc <- pick("^Corrected Akaike Information Criterion: *")
  b <- pick(".*Best-fit model according to BIC: *")
  if (is.na(b)) {
    h <- grep("^Best-fit model: ", ln, value = TRUE)
    if (length(h)) b <- trimws(sub(" chosen according.*$", "", sub("^Best-fit model: *", "", h[1])))
  }
  if (is.na(b)) b <- pick("^Model of substitution: *")
  out$best <- b
  a <- pick("^Gamma shape alpha: *")
  if (!is.na(a)) out$alpha <- suppressWarnings(as.numeric(strsplit(a, "\\s+")[[1]][1]))
  l <- pick("^Log-likelihood of the tree: *")
  if (!is.na(l)) out$loglik <- suppressWarnings(as.numeric(strsplit(l, "\\s+")[[1]][1]))
  i <- grep("^List of models sorted by BIC scores", ln)
  if (length(i)) {
    j <- i[1] + 1L
    while (j <= length(ln) && !grepl("^Model\\s", ln[j])) j <- j + 1L
    rows <- list(); j <- j + 1L
    while (j <= length(ln) && nzchar(trimws(ln[j]))) {
      tok <- strsplit(trimws(ln[j]), "\\s+")[[1]]
      num <- suppressWarnings(as.numeric(tok[-1]))
      num <- num[!is.na(num)]
      if (length(num) >= 7)
        rows[[length(rows) + 1L]] <- data.frame(Model = tok[1], LogL = num[1], AIC = num[2], wAIC = num[3],
                                                AICc = num[4], wAICc = num[5], BIC = num[6], wBIC = num[7],
                                                stringsAsFactors = FALSE)
      j <- j + 1L
    }
    if (length(rows)) out$table <- do.call(rbind, rows)
  }
  out
}

# Rough 0-100 progress from the IQ-TREE log (patterns of the IQ-TREE 2/3 log; an estimate, not exact).
.iq_progress <- function(log, nboot = 0) {
  if (!length(log)) return(0L)
  tx <- paste(log, collapse = "\n")
  has <- function(rx) grepl(rx, tx)
  if (has("Total wall-clock time used")) return(100L)
  p <- 2
  if (has("Alignment most likely contains|Alignment has")) p <- 5
  if (has("ModelFinder will test")) {
    p <- 10
    n <- suppressWarnings(as.numeric(sub(".*ModelFinder will test up to ([0-9]+).*", "\\1",
                                         regmatches(tx, regexpr("ModelFinder will test up to [0-9]+", tx)))))
    done <- sum(grepl("^\\s*[0-9]+\\s+[A-Za-z0-9+]+\\s+-?[0-9]+\\.[0-9]+", log))
    if (length(n) == 1 && !is.na(n) && n > 0) p <- 10 + 35 * min(1, done / n) else p <- 25
  }
  if (has("Best-fit model: ") || has("Best-fit model according")) p <- max(p, 45)
  if (has("Perform fast likelihood tree search|Create initial parsimony tree")) p <- max(p, 15)
  if (has("Optimizing NNI|nearest neighbor interchange")) p <- max(p, 55)
  if (has("Computing ML distances|Estimate model parameters")) p <- max(p, 50)
  if (has("ultrafast bootstrap|UFBoot|Generating [0-9]+ samples")) {
    p <- max(p, 70)
    b <- sum(grepl("^Bootstrap tree #?[0-9]+", log))
    if (nboot > 0 && b > 0) p <- max(p, 70 + 25 * min(1, b / nboot))
  }
  if (has("Consensus tree written|Analysis results written")) p <- max(p, 97)
  as.integer(min(99, round(p)))
}

.iq_args <- function(fasta, prefix, model_mode = "MFP", fixed_model = "", boot_mode = "none", nboot = 1000,
                     alrt = FALSE, threads = "AUTO", outgroup = NULL, seed = NULL, extra = "", quote_paths = FALSE) {
  q <- function(x) if (quote_paths) shQuote(x) else x
  a <- c("-s", q(fasta), "--prefix", q(prefix))
  a <- c(a, "-m", switch(model_mode, MFP = "MFP", MF = "MF", trimws(fixed_model)))
  if (!identical(model_mode, "MF")) {
    if (identical(boot_mode, "ufboot")) a <- c(a, "-B", as.character(max(1000, as.integer(nboot))), "-bnni")
    if (identical(boot_mode, "std"))    a <- c(a, "-b", as.character(as.integer(nboot)))
    if (isTRUE(alrt)) a <- c(a, "-alrt", "1000")
  }
  if (nzchar(trimws(threads))) a <- c(a, "-T", trimws(threads))
  if (length(outgroup)) a <- c(a, "-o", paste(outgroup, collapse = ","))
  if (!is.null(seed) && !is.na(seed)) a <- c(a, "-seed", as.character(as.integer(seed)))
  if (nzchar(trimws(extra))) a <- c(a, strsplit(trimws(extra), "\\s+")[[1]])
  c(a, "-redo")
}

# Distances + neighbour-joining. Returns list(tree, D).
.nj_build <- function(x, model, alpha) {
  d <- function(mm) {
    g <- if (!is.null(alpha) && is.finite(alpha) && alpha > 0 && !(model %in% c("raw", "logdet"))) alpha else FALSE
    ape::dist.dna(ape::as.DNAbin(mm), model = model, gamma = g, pairwise.deletion = TRUE, as.matrix = TRUE)
  }
  D <- d(x)
  if (any(!is.finite(D)))
    stop("Some distances cannot be computed with this model (saturated or non-overlapping sequences). Try 'raw' (p-distance) or JC69/K80.")
  tr <- ape::nj(stats::as.dist(D))
  tr$edge.length[tr$edge.length < 0] <- 0
  list(tree = tr, D = D, dist_fun = d)
}

.root_tree <- function(tr, mode = "none", outgroup = NULL) {
  if (identical(mode, "outgroup") && length(outgroup)) {
    og <- intersect(outgroup, tr$tip.label)
    if (length(og) && length(og) < length(tr$tip.label)) {
      r <- tryCatch(ape::root(tr, outgroup = og, resolve.root = TRUE, edgelabel = TRUE),
                    error = function(e) tryCatch(ape::root(tr, outgroup = og, resolve.root = TRUE), error = function(e2) tr))
      return(r)
    }
  }
  if (identical(mode, "midpoint") && requireNamespace("phangorn", quietly = TRUE))
    return(tryCatch(phangorn::midpoint(tr), error = function(e) tr))
  tr
}

# summary of the branch support of a tree (the last number of "SH-aLRT/UFBoot" counts)
# Tanglegram of two trees: tree A (left) and tree B (right); the names stand between the trees, one column per tree, so nothing overlaps.
# Each line joins the same sequence in both trees. Tips are spread evenly (rows are the same size in both trees).
.sg_tanglegram <- function(a, b, line_col = "#2563eb", titles = c("NJ tree", "IQ-TREE tree")) {
  n <- ape::Ntip(a); lab <- a$tip.label
  graphics::par(mar = c(0, 0, 0, 0)); graphics::plot.new()          # a plot must exist before text widths can be measured
  din <- graphics::par("din")[1]
  cex <- 0.7; gap <- max(0.9, 0.12 * din); tmin <- 1.3; disp <- lab
  wl <- max(graphics::strwidth(disp, units = "inches", cex = cex))
  if (2 * wl + gap + 2 * tmin > din) {                               # long names: smaller letters first, then shorter names
    cex <- max(0.42, cex * (din - gap - 2 * tmin) / (2 * wl))
    wl <- max(graphics::strwidth(disp, units = "inches", cex = cex))
    if (2 * wl + gap + 2 * tmin > din) {
      keep <- max(12, floor(max(nchar(lab)) * (din - gap - 2 * tmin) / (2 * wl)))
      disp <- ifelse(nchar(lab) > keep, paste0(substr(lab, 1, keep - 1), "\u2026"), lab)
      wl <- max(graphics::strwidth(disp, units = "inches", cex = cex))
    }
  }
  names(disp) <- lab
  tw <- max(0.8, (din - 2 * wl - gap - 0.2) / 2)                     # width of each tree (inches)
  ylim <- c(0.5, n + 1.8)
  graphics::par(xaxs = "i", yaxs = "i", mar = c(0, 0, 0, 0))
  env <- utils::getFromNamespace(".PlotPhyloEnv", "ape")
  tip_y <- function(tr) get("last_plot.phylo", envir = env)$yy[seq_len(ape::Ntip(tr))]
  graphics::par(fig = c(0, tw / din, 0, 1), new = TRUE)
  ape::plot.phylo(a, show.tip.label = FALSE, no.margin = TRUE, y.lim = ylim, edge.color = "#334155", edge.width = 0.9)
  ya <- tip_y(a); names(ya) <- a$tip.label
  graphics::par(fig = c((din - tw) / din, 1, 0, 1), new = TRUE)
  ape::plot.phylo(b, direction = "leftwards", show.tip.label = FALSE, no.margin = TRUE, y.lim = ylim, edge.color = "#334155", edge.width = 0.9)
  yb <- tip_y(b); names(yb) <- b$tip.label
  graphics::par(fig = c(0, 1, 0, 1), new = TRUE, mar = c(0, 0, 0, 0))
  graphics::plot.new(); graphics::plot.window(xlim = c(0, din), ylim = ylim, xaxs = "i", yaxs = "i")
  nm <- lab[lab %in% names(yb)]
  x1 <- tw + 0.05 + wl + 0.06; x2 <- din - tw - 0.05 - wl - 0.06
  graphics::segments(x1, ya[nm], x2, yb[nm], col = grDevices::adjustcolor(line_col, 0.55), lwd = 0.9)
  graphics::text(tw + 0.05, ya[nm], disp[nm], adj = c(0, 0.5), cex = cex, col = "#0f172a")
  graphics::text(din - tw - 0.05, yb[nm], disp[nm], adj = c(1, 0.5), cex = cex, col = "#0f172a")
  graphics::text(tw / 2, n + 1.2, titles[1], font = 2, cex = 0.9, col = "#0f172a")
  graphics::text(din - tw / 2, n + 1.2, titles[2], font = 2, cex = 0.9, col = "#0f172a")
  invisible(NULL)
}

.sup_summary <- function(tr, name) {
  v <- if (is.null(tr$node.label)) numeric() else suppressWarnings(as.numeric(sub(".*/", "", as.character(tr$node.label))))
  v <- v[!is.na(v)]
  data.frame(Tree = name, Tips = ape::Ntip(tr), `Nodes with support` = length(v),
             `Mean support` = if (length(v)) round(mean(v), 1) else NA,
             `>= 95` = sum(v >= 95), `70 - 94` = sum(v >= 70 & v < 95), `< 70` = sum(v < 70), check.names = FALSE)
}

.plot_tree <- function(tr, type = "phylogram", support = TRUE, minsup = 0, cex = 0.9, scale = TRUE, title = NULL, dots = FALSE) {
  op <- graphics::par(mar = c(2, 1, 2.2, 1)); on.exit(graphics::par(op))
  ape::plot.phylo(tr, type = type, use.edge.length = !identical(type, "cladogram"), cex = cex,
                  edge.width = 1.6, label.offset = if (type %in% c("phylogram", "cladogram")) 0.002 else 0, main = title)
  if (isTRUE(support) && !is.null(tr$node.label)) {
    lab <- as.character(tr$node.label)
    num <- suppressWarnings(as.numeric(sub(".*/", "", lab)))        # "SH-aLRT/UFBoot": the threshold uses the last value
    lab[is.na(num) | num < minsup] <- ""
    if (isTRUE(dots)) {            # coloured dots: green >= 95, orange 70-94, red < 70
      okn <- which(!is.na(num) & num >= minsup)
      if (length(okn)) ape::nodelabels(node = okn + ape::Ntip(tr), pch = 21, col = "white", cex = cex * 1.25,
                                       bg = ifelse(num[okn] >= 95, "#16a34a", ifelse(num[okn] >= 70, "#d97706", "#dc2626")))
    } else if (any(nzchar(lab))) ape::nodelabels(lab, frame = "none", adj = c(1.15, -0.35), cex = cex * 0.75, col = "#c0392b")
  }
  if (isTRUE(scale) && identical(type, "phylogram")) ape::add.scale.bar(cex = cex * 0.8)
  invisible()
}

# ---------------------------------------------------------------------------------------------------
# Shiny module
# ---------------------------------------------------------------------------------------------------

# A colour box (browser colour picker) that reports its value to Shiny when the user picks a colour.
# Until the user picks one, input$<id> stays NULL and the server uses its default.
.color_pick <- function(id, label, value) {
  shiny::div(style = "display:flex;flex-direction:column;align-items:center;font-size:12px",
    shiny::tags$input(type = "color", id = id, value = value, title = label,
                      style = "width:46px;height:32px;padding:1px;border:1px solid #b8c2cc;border-radius:6px;cursor:pointer;background:#fff",
                      onchange = sprintf("Shiny.setInputValue('%s', this.value)", id)),
    shiny::span(style = "margin-top:2px;color:#52606d", label))
}

# lets the server change a colour box (palette choice, reset): message 'sg_setcolor' with the full id and the value
.sg_color_js <- function() {
  shiny::tags$script(shiny::HTML("
    (function wait(){
      if (!(window.Shiny && Shiny.addCustomMessageHandler)) return setTimeout(wait, 30);
      if (window.sgColorReady) return; window.sgColorReady = true;
      Shiny.addCustomMessageHandler('sg_setcolor', function(m){
        var el = document.getElementById(m.id); if (el) el.value = m.value;
        Shiny.setInputValue(m.id, m.value);
      });
    })();"))
}

treeUI <- function(id) {
  ns <- shiny::NS(id)
  have_ape <- requireNamespace("ape", quietly = TRUE)
  dl <- function(x, lab) shiny::downloadButton(ns(x), lab)
  disp_btns <- function(prefix)
    shiny::div(style = "display:flex;gap:6px;flex-wrap:wrap;margin-top:8px",
               dl(paste0(prefix, "_nwk"), "Newick"), dl(paste0(prefix, "_png"), "PNG"), dl(paste0(prefix, "_pdf"), "PDF"))
  shiny::div(class = "container-fluid",
    if (!have_ape) shiny::div(class = "alert alert-danger", "The 'ape' package is required: install.packages(\"ape\")"),
    shiny::fluidRow(
      shiny::column(4,
        .card("1. Alignment",
          shiny::radioButtons(ns("src"), NULL,
            c("MSA tab: kept columns (after trimming / editing)" = "trim",
              "MSA tab: whole alignment" = "all",
              "An aligned FASTA file" = "file")),
          shiny::conditionalPanel(sprintf("input['%s'] == 'file'", ns("src")),
            shiny::fileInput(ns("fa"), NULL, buttonLabel = "Open aligned FASTA...", placeholder = "",
                             accept = c(".fa", ".fasta", ".fas", ".fna", ".txt"))),
          shiny::uiOutput(ns("aln_info"))),
        .card("2. Neighbour-joining",
          shiny::selectInput(ns("nj_model"), "Distance model",
            c("Auto: best-fit model from IQ-TREE ModelFinder" = "auto", "JC69" = "JC69", "K80 (Kimura 2-parameter)" = "K80",
              "F81" = "F81", "K81" = "K81", "F84" = "F84", "TN93" = "TN93", "p-distance (raw)" = "raw", "LogDet" = "logdet"),
            "K80"),
          shiny::conditionalPanel(sprintf("input['%s'] != 'auto' && input['%s'] != 'raw' && input['%s'] != 'logdet'", ns("nj_model"), ns("nj_model"), ns("nj_model")),
            shiny::numericInput(ns("nj_alpha"), "Gamma shape alpha (0 = equal rates)", 0, min = 0, step = 0.1)),
          shiny::numericInput(ns("nj_boot"), "Bootstrap replicates (0 = none)", 100, min = 0, step = 100),
          shiny::helpText("Auto runs ModelFinder once (IQ-TREE needed) and uses the closest distance model that ape has: ",
                          "HKY, TN, TPM, TIM, TVM, SYM and GTR all become TN93. Use IQ-TREE below for the exact model."),
          shiny::actionButton(ns("run_nj"), "Build NJ tree", icon = shiny::icon("play"), class = "btn-primary")),
        .card("3. IQ-TREE (maximum likelihood)",
          shiny::radioButtons(ns("iq_model"), "Substitution model",
            c("ModelFinder + tree (-m MFP)" = "MFP", "ModelFinder only (-m MF)" = "MF", "Fixed model" = "fixed")),
          shiny::conditionalPanel(sprintf("input['%s'] == 'fixed'", ns("iq_model")),
            shiny::textInput(ns("iq_fixed"), "Model", "GTR+G", placeholder = "HKY+F+G4")),
          shiny::selectInput(ns("iq_boot"), "Branch support",
            c("Ultrafast bootstrap (-B)" = "ufboot", "Standard bootstrap (-b, slow)" = "std", "None" = "none")),
          shiny::conditionalPanel(sprintf("input['%s'] != 'none'", ns("iq_boot")),
            shiny::numericInput(ns("iq_nboot"), "Replicates (ultrafast: at least 1000)", 1000, min = 100, step = 100)),
          shiny::checkboxInput(ns("iq_alrt"), "Also SH-aLRT (-alrt 1000)", FALSE),
          shiny::textInput(ns("iq_threads"), "Threads (AUTO or a number)", "AUTO"),
          shiny::numericInput(ns("iq_seed"), "Random seed (empty = random)", NA, step = 1),
          shiny::textInput(ns("iq_extra"), "Extra IQ-TREE arguments (optional)", "", placeholder = "-mset JC,K80,HKY,GTR"),
          shiny::uiOutput(ns("iq_out_ui")),
          shiny::uiOutput(ns("iq_path_ui")),     # only shown when IQ-TREE is not found automatically
          shiny::uiOutput(ns("iq_found")),
          shiny::div(style = "display:flex;gap:6px;flex-wrap:wrap",
            shiny::actionButton(ns("run_iq"), "Run IQ-TREE", icon = shiny::icon("play"), class = "btn-primary"),
            shiny::actionButton(ns("cancel_iq"), "Cancel", icon = shiny::icon("stop")))),
        .card("4. Display",
          shiny::selectizeInput(ns("outgroup"), "Outgroup (optional; also used by -o in IQ-TREE)", choices = NULL, multiple = TRUE),
          shiny::radioButtons(ns("root"), "Root the drawn tree",
            c("As computed" = "none", "On the outgroup" = "outgroup", "Midpoint (needs phangorn)" = "midpoint")),
          shiny::selectInput(ns("layout"), "Layout", c("Phylogram" = "phylogram", "Cladogram" = "cladogram",
                                                       "Unrooted" = "unrooted", "Fan" = "fan", "Radial" = "radial")),
          shiny::checkboxInput(ns("show_sup"), "Show branch support", TRUE),
          shiny::checkboxInput(ns("sup_dots"), "Support as coloured dots (green >= 95, orange 70-94, red < 70)", FALSE),
          shiny::sliderInput(ns("min_sup"), "Hide support below", 0, 100, 0, 5),
          shiny::sliderInput(ns("cex"), "Label size", 0.4, 2, 0.9, 0.1),
          shiny::checkboxInput(ns("scale"), "Scale bar (phylogram)", TRUE))),
      shiny::column(8,
        shiny::uiOutput(ns("msg")),
        shiny::tabsetPanel(id = ns("tabs"),
          shiny::tabPanel("NJ tree", shiny::br(), shiny::uiOutput(ns("nj_info")),
                          shiny::plotOutput(ns("nj_plot"), height = "auto"), disp_btns("nj")),
          shiny::tabPanel("IQ-TREE tree", shiny::br(), shiny::uiOutput(ns("iq_run_ui")),
                          shiny::plotOutput(ns("iq_plot"), height = "auto"), disp_btns("iq")),
          shiny::tabPanel("Interactive view", shiny::br(),
                          shiny::tags$script(src = "sangerly_treeview.js"),
                          shiny::div(id = ns("itree"), class = "sg-itree"),
                          shiny::helpText("The NJ and IQ-TREE trees as you rooted them in the Display box. Rectangular, circular or radial layout, ",
                                          "collapse a clade by clicking its node, flip it with a double-click, search a name, colour the tips by the first word ",
                                          "of the name (genus), and save the picture as SVG (editable in Inkscape / Illustrator) or PNG.")),
          shiny::tabPanel("Compare trees", shiny::br(), shiny::uiOutput(ns("cmp_info")),
                          shiny::plotOutput(ns("cmp_plot"), height = "auto"),
                          shiny::helpText("Each line joins the same sequence in the NJ tree (left) and the IQ-TREE tree (right). Lines that cross show where ",
                                          "the two methods disagree. The Robinson-Foulds distance counts the clades found in only one of the trees.")),
          shiny::tabPanel("ModelFinder", shiny::br(), shiny::uiOutput(ns("mf_best")),
                          DT::dataTableOutput(ns("mf_tbl")),
                          shiny::helpText("Models sorted by BIC (lower is better). w-BIC is the model weight.")),
          shiny::tabPanel("IQ-TREE log and report", shiny::br(),
                          shiny::verbatimTextOutput(ns("iq_log"), placeholder = TRUE),
                          shiny::div(style = "display:flex;gap:6px;flex-wrap:wrap;margin-top:8px",
                                     dl("dl_report", "Report (.iqtree)"), dl("dl_log", "Log (.log)"), dl("dl_treefile", "Treefile"))),
          shiny::tabPanel("Save results", shiny::br(),
            shiny::textInput(ns("sv_name"), "Name for the files (empty = the default names)", "",
                             placeholder = "e.g. COI_Odonata", width = "100%"),
            shiny::uiOutput(ns("sv_local_ui")),
            shiny::helpText("The Download buttons under the trees and tables use this name as well; the browser decides where ",
                            "they are saved. To choose the folder yourself, use 'Save to this folder' (desktop version).")),
          shiny::tabPanel("Distances", shiny::br(), shiny::uiOutput(ns("dist_info")), .sg_color_js(),
            .card("Heat map appearance",
              shiny::tags$details(open = NA,
                shiny::tags$summary(style = "cursor:pointer;font-weight:600;margin-bottom:10px", "Settings (click to hide or show)"),
                shiny::fluidRow(
                  shiny::column(4,
                    shiny::radioButtons(ns("hm_layout"), "Shape",
                      c("Lower triangle (default)" = "lower", "Full square (every pair shown twice)" = "full")),
                    shiny::radioButtons(ns("hm_scale"), "Colour scale",
                      c("Whole matrix (legend = real min / max)" = "matrix", "Each row (own min to max)" = "row"))),
                  shiny::column(4,
                    shiny::selectInput(ns("hm_palette"), "Colour palette", .HM_PALETTE_CHOICES),
                    shiny::div(style = "display:flex;gap:14px;flex-wrap:wrap",
                      .color_pick(ns("hm_low"), "Low", "#0000FF"), .color_pick(ns("hm_mid"), "Middle", "#FFFFFF"),
                      .color_pick(ns("hm_high"), "High", "#FF0000"), .color_pick(ns("hm_grey"), "Empty cells", "#BFBFBF")),
                    shiny::helpText("Choose a palette, then click a colour box to change that colour. 'Empty cells' is the unused ",
                                    "upper triangle. The numbers turn white or black by themselves so they stay readable.")),
                  shiny::column(4,
                    shiny::checkboxInput(ns("hm_num"), "Numbers in the cells", TRUE),
                    shiny::numericInput(ns("hm_dig"), "Decimals", 2, min = 0, max = 5, step = 1),
                    shiny::sliderInput(ns("hm_numsize"), "Number size", 5, 14, 8, 1),
                    shiny::sliderInput(ns("hm_labsize"), "Name size", 6, 14, 9, 1))),
                shiny::fluidRow(
                  shiny::column(4, shiny::checkboxInput(ns("hm_border"), "White lines between the cells", TRUE)),
                  shiny::column(5, shiny::textInput(ns("hm_title"), "Title (optional)", "", width = "100%")),
                  shiny::column(3, shiny::div(style = "margin-top:25px",
                    shiny::actionButton(ns("hm_reset"), "Reset appearance", icon = shiny::icon("rotate-left"))))))),
            shiny::helpText("The picture grows with the number of sequences and decimals so every number stays readable; scroll inside the frame."),
            shiny::div(style = "overflow:auto;max-height:82vh;border:1px solid #e3e8ee;border-radius:6px;background:#fff",
              shiny::plotOutput(ns("dist_heat"), width = "auto", height = "auto")),
            shiny::div(style = "display:flex;gap:6px;flex-wrap:wrap;margin:8px 0",
              dl("dl_dist_xlsx", "Excel (matrix + pairs)"), dl("dl_dist", "Matrix (CSV)"), dl("dl_pairs", "Pairs (CSV)"),
              dl("dl_hm_svg", "Heat map (SVG)"), dl("dl_hm_jpg", "Heat map (JPG)"),
              dl("dl_hm_png", "Heat map (PNG)"), dl("dl_hm_pdf", "Heat map (PDF)")),
            DT::dataTableOutput(ns("dist_tbl")))))))
}

treeServer <- function(id, aln_all = shiny::reactive(NULL), aln_trim = shiny::reactive(NULL)) {
  shiny::moduleServer(id, function(input, output, session) {
    rv <- shiny::reactiveValues(file_aln = NULL, nj = NULL, iq = NULL, msg = NULL, running = FALSE, proc = NULL,
                                iq_log = character(), iq_pct = 0L, t0 = NULL, tick = 0L, dir = NULL, prefix = NULL,
                                nboot = 0, opts = NULL)
    set_msg <- function(txt, type = "info") rv$msg <- if (is.null(txt)) NULL else list(txt = txt, type = type)

    # ---- alignment --------------------------------------------------------------------------------
    shiny::observeEvent(input$fa, {
      recs <- tryCatch(blast_parse_records(paste(readLines(input$fa$datapath, warn = FALSE), collapse = "\n")), error = function(e) NULL)
      if (is.null(recs) || nrow(recs) < 3) { set_msg("The file needs at least 3 aligned sequences.", "warning"); rv$file_aln <- NULL; return() }
      s <- stats::setNames(gsub("[.?]", "-", recs$seq), recs$id)
      if (length(unique(nchar(s))) != 1) { set_msg("The sequences in this file are not all the same length, so it is not an alignment.", "danger"); rv$file_aln <- NULL; return() }
      rv$file_aln <- s; set_msg(NULL)
    })
    # Qt dialog "Build a tree": an aligned FASTA file chosen there (a path on this computer)
    shiny::observeEvent(input$qt_fasta, {
      p <- as.character(input$qt_fasta)
      recs <- tryCatch(blast_parse_records(paste(readLines(p, warn = FALSE), collapse = "\n")), error = function(e) NULL)
      if (is.null(recs) || nrow(recs) < 3) { set_msg("The file needs at least 3 aligned sequences.", "warning"); rv$file_aln <- NULL; return() }
      s <- stats::setNames(gsub("[.?]", "-", recs$seq), recs$id)
      if (length(unique(nchar(s))) != 1) { set_msg("The sequences in this file are not all the same length, so it is not an alignment.", "danger"); rv$file_aln <- NULL; return() }
      rv$file_aln <- s; set_msg(NULL)
    })
    # the alignment was cleared in the MSA tab: the old trees and distances must not stay on screen
    shiny::observeEvent(list(aln_all(), aln_trim()), {
      if (is.null(aln_all()) && is.null(aln_trim()) && !isTRUE(rv$running)) {
        rv$nj <- NULL; rv$iq <- NULL; rv$iq_log <- character(); rv$iq_pct <- 0L
      }
    }, ignoreInit = TRUE)
    prep <- shiny::reactive({
      x <- switch(.nz(input$src, "trim"), trim = aln_trim(), all = aln_all(), file = rv$file_aln)
      .tree_prepare(x)
    })
    aln <- shiny::reactive({ p <- prep(); if (is.null(p)) NULL else p$aln })
    shiny::observe({
      a <- aln()
      shiny::updateSelectizeInput(session, "outgroup", choices = if (is.null(a)) character() else names(a), server = FALSE)
    })
    output$aln_info <- shiny::renderUI({
      p <- prep()
      if (is.null(p)) return(shiny::helpText("No alignment yet (at least 3 sequences). Align in the MSA tab, or open an aligned FASTA."))
      a <- p$aln; m <- .tree_matrix(a)
      nvar <- sum(apply(m, 2, function(cc) length(unique(cc[!cc %in% c("-", "n", "?", "x")])) > 1))
      shiny::tagList(
        shiny::helpText(sprintf("%d sequences, %d columns, %d variable sites, %.1f%% gaps.", length(a), ncol(m), nvar,
                                100 * mean(m == "-"))),
        if (p$renamed > 0) shiny::helpText(sprintf("%d name(s) were cleaned (spaces and symbols become '_') for the tree programs.", p$renamed)),
        if (length(p$dropped)) shiny::div(class = "alert alert-warning", "Left out (no bases at all): ", paste(p$dropped, collapse = ", ")))
    })

    # ---- IQ-TREE binary -----------------------------------------------------------------------------
    output$iq_path_ui <- shiny::renderUI({
      if (nzchar(.iq_find(""))) NULL      # bundled or on the PATH: no need to ask for a path
      else shiny::textInput(session$ns("iq_path"), "Path to IQ-TREE (not found automatically)", "",
                            placeholder = "D:/alponsin/sangerly/IQTREE/bin/iqtree3.exe")
    })
    iq_bin <- shiny::reactive(.iq_find(.nz(input$iq_path, "")))
    output$iq_out_ui <- shiny::renderUI({
      if (!(nzchar(Sys.getenv("SANGERLY_DESKTOP")) || interactive())) return(NULL)
      shiny::tagList(
        shiny::textInput(session$ns("iq_outdir"), "Folder for the IQ-TREE files (empty = temporary folder)", "",
                         placeholder = "C:/Users/you/Documents/Sangerly/trees", width = "100%"),
        shiny::actionButton(session$ns("iq_outbrowse"), "Browse folder..."),
        shiny::helpText("Each run gets its own sub-folder here (named after the name in 'Save results', or 'iqtree', plus the date and time) ",
                        "with the log, report, treefile and all the other IQ-TREE files."))
    })
    shiny::observeEvent(input$iq_outbrowse, {
      sg_dirpick(session, function(d) shiny::updateTextInput(session, "iq_outdir", value = gsub("\\\\", "/", d)),
                 "Choose the folder for the IQ-TREE files", if (is.null(input$iq_outdir)) "" else input$iq_outdir)
    })
    output$iq_found <- shiny::renderUI({
      last <- if (!is.null(rv$dir) && nzchar(trimws(if (is.null(input$iq_outdir)) "" else input$iq_outdir)))
        shiny::helpText(paste("Last run folder:", rv$dir)) else NULL
      if (nzchar(iq_bin())) shiny::tagList(shiny::helpText(paste("IQ-TREE found:", iq_bin())), last)
      else shiny::div(class = "alert alert-warning",
                      "IQ-TREE was not found. Download iqtree2/iqtree3 from iqtree.github.io, put it in the app folder ",
                      "(IQTREE/...), add it to the PATH, or type its full path above. The NJ tree works without it ",
                      "(except the 'Auto' model).")
    })

    # ---- IQ-TREE: result parsing -----------------------------------------------------------------------
    finish_iq <- function(exit = NA) {
      pre <- rv$prefix; tf <- paste0(pre, ".treefile"); lg <- paste0(pre, ".log")
      rv$iq_log <- if (file.exists(lg)) readLines(lg, warn = FALSE) else rv$iq_log
      if (!file.exists(tf)) {
        set_msg(paste("IQ-TREE ended without a tree.", if (!is.na(exit)) paste("Exit status", exit, "."),
                      "Last lines of the log:", paste(utils::tail(rv$iq_log, 6), collapse = " | ")), "danger")
        rv$iq_pct <- 0L; return(invisible())
      }
      rep <- .iq_parse_report(paste0(pre, ".iqtree"))
      use <- if (identical(rv$opts$boot, "std") && file.exists(paste0(pre, ".contree"))) paste0(pre, ".contree") else tf
      tr <- tryCatch(ape::read.tree(use), error = function(e) NULL)
      rv$iq <- list(tree = tr, report = rep, prefix = pre, aln = rv$opts$aln, treefile = tf, mf_only = identical(rv$opts$model, "MF"),
                    opts = rv$opts[c("model", "boot", "nboot", "alrt", "fixed")], version = .iq_version(rv$iq_log))
      rv$iq_pct <- 100L
      set_msg(NULL)
      shiny::showNotification(sprintf("IQ-TREE finished%s.", if (!is.na(rep$best)) paste0(" (model ", rep$best, ")") else ""), type = "message")
      invisible()
    }

    # folder of one IQ-TREE run: the folder chosen by the researcher (a sub-folder with date and time), else a temporary one
    run_dir <- function(kind) {
      od <- trimws(if (is.null(input$iq_outdir)) "" else input$iq_outdir)
      nm <- sname()
      if (nzchar(od)) {
        base <- file.path(od, paste0(if (nzchar(nm)) nm else kind, "_", format(Sys.time(), "%Y%m%d_%H%M%S")))
        ok <- tryCatch({ dir.create(base, showWarnings = FALSE, recursive = TRUE); dir.exists(base) }, error = function(e) FALSE)
        if (isTRUE(ok)) return(base)
        set_msg("The folder for the IQ-TREE files cannot be created; a temporary folder is used instead.", "warning")
      }
      d <- file.path(tempdir(), format(Sys.time(), paste0(kind, "_%H%M%S")))
      dir.create(d, showWarnings = FALSE, recursive = TRUE)
      d
    }
    run_base <- function() { nm <- sname(); if (nzchar(nm)) nm else "run" }

    launch_iq <- function() {
      a <- aln(); bin <- iq_bin()
      if (is.null(a)) { set_msg("Provide an alignment with at least 3 sequences first.", "warning"); return(invisible()) }
      if (!nzchar(bin)) { set_msg("IQ-TREE was not found (see the message under the IQ-TREE settings).", "warning"); return(invisible()) }
      if (isTRUE(rv$running)) { set_msg("IQ-TREE is already running.", "info"); return(invisible()) }
      dir <- run_dir("iqtree")
      fasta <- file.path(dir, "alignment.fasta"); prefix <- file.path(dir, run_base())
      writeLines(msa_fasta_text(a), fasta)
      og <- input$outgroup
      opts <- list(model = input$iq_model, boot = input$iq_boot, nboot = .nz(input$iq_nboot, 1000), aln = a,
                   alrt = isTRUE(input$iq_alrt), fixed = .nz(input$iq_fixed, ""))
      use_px <- requireNamespace("processx", quietly = TRUE)
      args <- .iq_args(fasta, prefix, input$iq_model, .nz(input$iq_fixed, ""), input$iq_boot, opts$nboot,
                       isTRUE(input$iq_alrt), .nz(input$iq_threads, "AUTO"), og,
                       if (is.null(input$iq_seed) || is.na(input$iq_seed)) NULL else input$iq_seed,
                       .nz(input$iq_extra, ""), quote_paths = !use_px)
      rv$dir <- dir; rv$prefix <- prefix; rv$opts <- opts; rv$t0 <- Sys.time(); rv$iq <- NULL
      rv$iq_log <- character(); rv$iq_pct <- 1L; set_msg(NULL)
      if (use_px) {
        rv$proc <- tryCatch(processx::process$new(bin, args, wd = dir, stdout = file.path(dir, "stdout.txt"),
                                                  stderr = "2>&1", cleanup = TRUE),
                            error = function(e) { set_msg(paste("Could not start IQ-TREE:", conditionMessage(e)), "danger"); NULL })
        if (!is.null(rv$proc)) rv$running <- TRUE
      } else {
        shiny::showNotification("The 'processx' package is not installed, so IQ-TREE runs in the foreground and the app waits. install.packages(\"processx\") for a background run with a live log.",
                                type = "warning", duration = 10)
        out <- shiny::withProgress(message = "IQ-TREE is running (the app waits)...", value = 0.5,
          tryCatch(suppressWarnings(system2(bin, args, stdout = TRUE, stderr = TRUE, wait = TRUE)),
                   error = function(e) { set_msg(paste("Could not run IQ-TREE:", conditionMessage(e)), "danger"); NULL }))
        if (!is.null(out)) finish_iq()
      }
      invisible()
    }
    shiny::observeEvent(input$run_iq, launch_iq())

    # background run: poll the log once a second
    shiny::observe({
      shiny::req(rv$running)
      shiny::invalidateLater(1000, session)
      p <- shiny::isolate(rv$proc); pre <- shiny::isolate(rv$prefix)
      lg <- paste0(pre, ".log")
      if (file.exists(lg)) rv$iq_log <- readLines(lg, warn = FALSE)
      rv$iq_pct <- .iq_progress(shiny::isolate(rv$iq_log), if (identical(shiny::isolate(rv$opts$boot), "none")) 0 else shiny::isolate(rv$opts$nboot))
      rv$tick <- shiny::isolate(rv$tick) + 1L
      if (!is.null(p) && !p$is_alive()) { rv$running <- FALSE; shiny::isolate(finish_iq(p$get_exit_status())) }
    })
    shiny::observeEvent(input$cancel_iq, {
      if (isTRUE(rv$running) && !is.null(rv$proc)) {
        try(rv$proc$kill(), silent = TRUE); rv$running <- FALSE; rv$iq_pct <- 0L; set_msg("IQ-TREE was cancelled.", "info")
      } else set_msg("No IQ-TREE run is active.", "info")
    })
    output$iq_run_ui <- shiny::renderUI({
      rv$tick
      if (isTRUE(rv$running)) {
        el <- as.integer(difftime(Sys.time(), rv$t0, units = "secs"))
        shiny::tagList(
          shiny::div(class = "progress", style = "height:22px",
            shiny::div(class = "progress-bar progress-bar-striped progress-bar-animated", role = "progressbar",
                       style = sprintf("width:%d%%", max(rv$iq_pct, 3L)), sprintf("%d %%", rv$iq_pct))),
          shiny::helpText(sprintf("Running for %d min %02d s. The percentage is estimated from the IQ-TREE log (see the log tab).", el %/% 60L, el %% 60L)))
      } else if (!is.null(rv$iq) && isTRUE(rv$iq$mf_only))
        shiny::div(class = "alert alert-info", "ModelFinder only: see the ModelFinder tab (no tree is built with -m MF).")
    })

    # ---- NJ ---------------------------------------------------------------------------------------------
    shiny::observeEvent(input$run_nj, {
      a <- aln()
      if (is.null(a)) { set_msg("Provide an alignment with at least 3 sequences first.", "warning"); return() }
      if (!requireNamespace("ape", quietly = TRUE)) { set_msg("The 'ape' package is not installed: install.packages(\"ape\")", "danger"); return() }
      prog <- .msa_progress(session, "Neighbour-joining"); on.exit(prog$close(), add = TRUE)
      model <- .nz(input$nj_model, "K80"); alpha <- .nz(input$nj_alpha, 0); note <- NULL
      if (identical(model, "auto")) {
        rep <- NULL
        if (!is.null(rv$iq) && identical(rv$iq$aln, a)) rep <- rv$iq$report
        if (is.null(rep) || is.na(rep$best)) {
          bin <- iq_bin()
          if (!nzchar(bin)) { set_msg("'Auto' needs IQ-TREE for ModelFinder, and it was not found. Pick a model by hand or install IQ-TREE.", "warning"); return() }
          prog$set(5, "running ModelFinder (IQ-TREE); this can take a few minutes")
          dir <- run_dir("mf")
          fasta <- file.path(dir, "alignment.fasta"); writeLines(msa_fasta_text(a), fasta)
          args <- .iq_args(fasta, file.path(dir, "mf"), "MF", threads = .nz(input$iq_threads, "AUTO"), quote_paths = TRUE)
          out <- tryCatch(suppressWarnings(system2(bin, args, stdout = TRUE, stderr = TRUE)), error = function(e) NULL)
          rep <- .iq_parse_report(file.path(dir, "mf.iqtree"))
          if (is.na(rep$best)) { set_msg("ModelFinder gave no result. Pick a model by hand or look at the IQ-TREE tab.", "danger"); return() }
          rv$iq <- list(tree = NULL, report = rep, prefix = file.path(dir, "mf"), aln = a, treefile = NULL, mf_only = TRUE)
        }
        mp <- .iq_to_ape(rep$best); model <- mp$model
        alpha <- if (grepl("\\+G", rep$best) && is.finite(rep$alpha)) rep$alpha else 0
        note <- sprintf("ModelFinder: %s -> distance model %s%s%s.", rep$best, model,
                        if (alpha > 0) sprintf(" with gamma alpha = %.3f", alpha) else "",
                        if (!mp$exact) " (closest model available in ape)" else "")
        if (grepl("\\+R", rep$best)) note <- paste(note, "FreeRate (+R) cannot be used in ape distances; rates were left equal.")
      }
      prog$set(15, "distances and tree")
      res <- tryCatch(.nj_build(.tree_matrix(a), model, alpha), error = function(e) { set_msg(conditionMessage(e), "danger"); NULL })
      if (is.null(res)) return()
      tr <- res$tree; B <- as.integer(.nz(input$nj_boot, 0))
      if (!is.na(B) && B > 0) {
        m <- .tree_matrix(a); nc <- ncol(m); trees <- list()
        for (b in seq_len(B)) {
          prog$set(15 + 80 * b / B, sprintf("bootstrap %d/%d", b, B))
          rep_tr <- tryCatch(.nj_build(m[, sample.int(nc, nc, replace = TRUE), drop = FALSE], model, alpha)$tree, error = function(e) NULL)
          if (!is.null(rep_tr)) trees[[length(trees) + 1L]] <- rep_tr
        }
        if (length(trees) >= 2) {
          class(trees) <- "multiPhylo"
          bp <- tryCatch(ape::prop.clades(tr, trees), error = function(e) NULL)
          if (!is.null(bp)) tr$node.label <- ifelse(is.na(bp), "", as.character(round(100 * bp / length(trees))))
          if (length(trees) < B) note <- paste(note, sprintf("%d of %d bootstrap replicates could not be computed and were skipped.", B - length(trees), B))
        } else note <- paste(note, "Too many bootstrap replicates failed; no support values.")
      }
      prog$set(100, "done")
      rv$nj <- list(tree = tr, D = res$D, model = model, alpha = alpha, B = B, note = note, n = length(a))
      set_msg(NULL)
    })

    # ---- drawing ---------------------------------------------------------------------------------------------
    shown <- function(tr) .root_tree(tr, .nz(input$root, "none"), input$outgroup)
    draw <- function(tr, title) {
      .plot_tree(shown(tr), .nz(input$layout, "phylogram"), isTRUE(input$show_sup), .nz(input$min_sup, 0),
                 .nz(input$cex, 0.9), isTRUE(input$scale), title, isTRUE(input$sup_dots))
    }
    plot_h <- function(tr) {
      if (is.null(tr)) return(60)
      n <- ape::Ntip(tr)
      if (.nz(input$layout, "phylogram") %in% c("phylogram", "cladogram")) max(380, 22 * n * max(0.6, .nz(input$cex, 0.9)) + 80) else 640
    }
    nj_tree <- shiny::reactive(if (is.null(rv$nj)) NULL else rv$nj$tree)
    iq_tree <- shiny::reactive(if (is.null(rv$iq)) NULL else rv$iq$tree)
    # interactive viewer: Newick of the drawn (rooted) trees
    shiny::observe({
      trs <- list()
      nj <- nj_tree(); iq <- iq_tree()
      if (!is.null(nj)) trs[["NJ"]] <- tryCatch(ape::write.tree(shown(nj)), error = function(e) NULL)
      if (!is.null(iq)) trs[["IQ-TREE"]] <- tryCatch(ape::write.tree(shown(iq)), error = function(e) NULL)
      trs <- trs[!vapply(trs, is.null, TRUE)]
      session$sendCustomMessage("sg_tree_data", list(id = session$ns("itree"), trees = trs))
    })
    # compare the NJ and the IQ-TREE tree
    cmp_pair <- shiny::reactive({
      a <- nj_tree(); b <- iq_tree(); shiny::req(a, b)
      a <- shown(a); b <- shown(b)
      common <- intersect(a$tip.label, b$tip.label)
      list(a = a, b = b, common = common)
    })
    output$cmp_info <- shiny::renderUI({
      if (is.null(nj_tree()) || is.null(iq_tree()))
        return(shiny::div(class = "alert alert-info", "Build the NJ tree and the IQ-TREE tree first; then they are compared here."))
      x <- cmp_pair()
      if (length(x$common) < 4) return(shiny::div(class = "alert alert-warn", "The two trees share fewer than 4 sequence names."))
      rf <- tryCatch(ape::dist.topo(ape::unroot(ape::keep.tip(x$a, x$common)), ape::unroot(ape::keep.tip(x$b, x$common))), error = function(e) NA)
      mx <- 2 * (length(x$common) - 3)
      sm <- rbind(.sup_summary(x$a, "NJ"), .sup_summary(x$b, "IQ-TREE"))
      shiny::tagList(
        shiny::div(class = "alert alert-info",
          shiny::strong(sprintf("%d shared sequences. ", length(x$common))),
          if (is.finite(rf)) sprintf("Robinson-Foulds distance %d of at most %d (%.0f%% of the clades differ).", as.integer(rf), as.integer(mx), 100 * rf / max(mx, 1))
          else "The Robinson-Foulds distance could not be computed."),
        shiny::tags$table(class = "table table-sm",
          shiny::tags$thead(shiny::tags$tr(lapply(names(sm), shiny::tags$th))),
          shiny::tags$tbody(lapply(seq_len(nrow(sm)), function(i) shiny::tags$tr(lapply(sm[i, ], function(v) shiny::tags$td(as.character(v))))))))
    })
    output$cmp_plot <- shiny::renderPlot({
      x <- cmp_pair(); shiny::req(length(x$common) >= 4)
      a <- ape::ladderize(ape::keep.tip(x$a, x$common)); b <- ape::ladderize(ape::keep.tip(x$b, x$common))
      ok <- tryCatch({ .sg_tanglegram(a, b); TRUE }, error = function(e) FALSE)
      if (!ok) {                                                      # fallback: the plain ape tanglegram
        graphics::par(new = FALSE, fig = c(0, 1, 0, 1), mar = c(1, 1, 1, 1)); graphics::plot.new()
        ape::cophyloplot(a, b, assoc = cbind(x$common, x$common), length.line = 3, space = 90, gap = 2,
                         col = "#2563eb", lwd = 1.2, show.tip.label = TRUE, font = 1, cex = 0.6)
      }
    }, height = function() { a <- tryCatch(nj_tree(), error = function(e) NULL); if (is.null(a)) 60 else max(460, 15 * ape::Ntip(a) + 120) })
    output$nj_plot <- shiny::renderPlot({
      tr <- nj_tree(); shiny::req(tr)
      draw(tr, sprintf("Neighbour-joining (%s)", rv$nj$model))
    }, height = function() plot_h(nj_tree()))
    output$iq_plot <- shiny::renderPlot({
      tr <- iq_tree(); shiny::req(tr)
      draw(tr, sprintf("IQ-TREE%s", if (!is.na(rv$iq$report$best)) paste0(" (", rv$iq$report$best, ")") else ""))
    }, height = function() plot_h(iq_tree()))
    output$nj_info <- shiny::renderUI({
      r <- rv$nj
      if (is.null(r)) return(shiny::helpText("Press 'Build NJ tree'."))
      shiny::helpText(paste0(r$n, " sequences, model ", r$model, if (r$alpha > 0) sprintf(", gamma alpha %.3f", r$alpha) else "",
                             if (r$B > 0) paste0(", ", r$B, " bootstrap replicates (support in %)") else ", no bootstrap", ". ", .nz(r$note, "")))
    })

    # ---- ModelFinder, log, distances ---------------------------------------------------------------------
    output$mf_best <- shiny::renderUI({
      r <- rv$iq$report
      if (is.null(r) || is.na(r$best)) return(shiny::helpText("Run IQ-TREE with ModelFinder (-m MFP or -m MF) to see the result here."))
      shiny::tagList(
        shiny::div(class = "alert alert-success", shiny::tags$b("Best-fit model: "), r$best),
        shiny::helpText(paste0("Best by AIC: ", .nz(r$aic, "-"), "  |  AICc: ", .nz(r$aicc, "-"), "  |  BIC: ", .nz(r$bic, "-"),
                               if (is.finite(r$loglik)) sprintf("  |  tree log-likelihood: %.2f", r$loglik) else "")))
    })
    output$mf_tbl <- DT::renderDataTable({
      tb <- rv$iq$report$table; shiny::req(tb)
      dt <- DT::datatable(tb, rownames = FALSE, options = list(pageLength = 15, scrollX = TRUE, dom = "tp"))
      dt <- DT::formatRound(dt, c("LogL", "AIC", "AICc", "BIC"), 2)
      DT::formatRound(dt, c("wAIC", "wAICc", "wBIC"), 4)
    })
    output$iq_log <- shiny::renderText({
      rv$tick
      lg <- rv$iq_log
      if (!length(lg)) "No IQ-TREE run yet." else paste(utils::tail(lg, 60), collapse = "\n")
    })
    output$dist_tbl <- DT::renderDataTable({
      D <- rv$nj$D; shiny::req(D)
      DT::datatable(round(as.data.frame(D), 5), options = list(pageLength = 15, scrollX = TRUE, dom = "tp"))
    })

    # ---- downloads -----------------------------------------------------------------------------------------------
    # ---- distance heat map and Excel export ------------------------------------------------------------------------
    output$dist_info <- shiny::renderUI({
      shiny::req(rv$nj)
      shiny::helpText(sprintf("Distances from the NJ run: model %s%s, %d sequences. Choose the model in the settings on the left and run NJ again to change them.",
                              rv$nj$model, if (!is.null(rv$nj$alpha) && is.finite(rv$nj$alpha) && rv$nj$alpha > 0)
                                sprintf(" (gamma alpha = %s)", format(rv$nj$alpha)) else "", nrow(rv$nj$D)))
    })
    hm_dig <- shiny::reactive({
      d <- input$hm_dig
      if (is.null(d) || is.na(d) || d < 0) 2L else as.integer(min(5, d))
    })
    # appearance chosen by the user (colour boxes stay NULL until touched, so the defaults apply)
    hm_args <- shiny::reactive(list(
      low = .nz(input$hm_low, "#0000FF"), mid = .nz(input$hm_mid, "#FFFFFF"), high = .nz(input$hm_high, "#FF0000"),
      na_col = .nz(input$hm_grey, "#BFBFBF"),
      num_pt = .nz(input$hm_numsize, 8), ax_pt = .nz(input$hm_labsize, 9),
      layout = .nz(input$hm_layout, "lower"), border = !identical(input$hm_border, FALSE)))
    heat_plot <- shiny::reactive({
      D <- rv$nj$D; shiny::req(D)
      sc <- input$hm_scale
      if (is.null(sc)) sc <- "matrix"
      ttl <- trimws(.nz(input$hm_title, ""))
      do.call(.dist_heatmap, c(list(D, digits = hm_dig(), show_num = isTRUE(input$hm_num), scale_by = sc,
                                    title = if (nzchar(ttl)) ttl else NULL), hm_args()))
    })
    hm_dims <- function() {
      D <- rv$nj$D
      if (is.null(D)) return(list(w = 1, h = 1))
      a <- hm_args()
      .dist_dims(D, hm_dig(), isTRUE(input$hm_num), a$num_pt, a$ax_pt)
    }
    output$dist_heat <- shiny::renderPlot({ heat_plot() }, res = 96,
      width = function() as.integer(min(20000, max(300, 96 * hm_dims()$w))),
      height = function() as.integer(min(20000, max(100, 96 * hm_dims()$h))))
    save_hm_file <- function(type, file) {
      shiny::req(rv$nj); a <- hm_args()
      .dist_save(heat_plot(), file, type, rv$nj$D, hm_dig(), isTRUE(input$hm_num), a$num_pt, a$ax_pt)
    }
    save_hm <- function(type) function(file) save_hm_file(type, file)

    # choosing a palette fills the colour boxes; "Reset appearance" puts every heat map setting back to its default
    set_col <- function(id, value) session$sendCustomMessage("sg_setcolor", list(id = session$ns(id), value = value))
    set_palette <- function(key) {
      p <- .hm_palette(key)
      set_col("hm_low", p[1]); set_col("hm_mid", p[2]); set_col("hm_high", p[3])
    }
    shiny::observeEvent(input$hm_palette, set_palette(input$hm_palette), ignoreInit = TRUE)
    shiny::observeEvent(input$hm_reset, {
      shiny::updateSelectInput(session, "hm_palette", selected = "bwr")
      set_palette("bwr"); set_col("hm_grey", "#BFBFBF")
      shiny::updateRadioButtons(session, "hm_layout", selected = "lower")
      shiny::updateRadioButtons(session, "hm_scale", selected = "matrix")
      shiny::updateCheckboxInput(session, "hm_num", value = TRUE)
      shiny::updateNumericInput(session, "hm_dig", value = 2)
      shiny::updateSliderInput(session, "hm_numsize", value = 8)
      shiny::updateSliderInput(session, "hm_labsize", value = 9)
      shiny::updateCheckboxInput(session, "hm_border", value = TRUE)
      shiny::updateTextInput(session, "hm_title", value = "")
    })
    output$dl_hm_svg <- shiny::downloadHandler(function() fn("distance_heatmap.svg", "_distance_heatmap.svg"), save_hm("svg"))
    output$dl_hm_jpg <- shiny::downloadHandler(function() fn("distance_heatmap.jpg", "_distance_heatmap.jpg"), save_hm("jpg"))
    output$dl_hm_png <- shiny::downloadHandler(function() fn("distance_heatmap.png", "_distance_heatmap.png"), save_hm("png"))
    output$dl_hm_pdf <- shiny::downloadHandler(function() fn("distance_heatmap.pdf", "_distance_heatmap.pdf"), save_hm("pdf"))
    output$dl_pairs <- shiny::downloadHandler(function() fn("distances_pairs.csv", "_distance_pairs.csv"),
      function(file) { shiny::req(rv$nj); utils::write.csv(.dist_pairs(rv$nj$D), file, row.names = FALSE) })
    output$dl_dist_xlsx <- shiny::downloadHandler(function() fn("genetic_distances.xlsx", "_distances.xlsx"),
      function(file) {
        shiny::req(rv$nj)
        D <- as.matrix(rv$nj$D)
        info <- data.frame(Item = c("Distance model", "Gamma alpha", "Number of sequences", "Created"),
                           Value = c(as.character(rv$nj$model),
                                     if (is.null(rv$nj$alpha) || !is.finite(rv$nj$alpha)) "none" else format(rv$nj$alpha),
                                     as.character(nrow(D)), format(Sys.time(), "%Y-%m-%d %H:%M")),
                           stringsAsFactors = FALSE)
        tryCatch(.sg_xlsx(list(Matrix = .dist_matrix_df(D), Pairs = .dist_pairs(D), Info = info), file),
          error = function(e) {
            shiny::showNotification(paste("Excel file not written:", conditionMessage(e)), type = "error", duration = 15)
            writeLines(conditionMessage(e), file)
          })
      })

    save_nwk <- function(tr) function(file) { shiny::req(tr()); ape::write.tree(shown(tr()), file) }
    save_png <- function(tr, ttl) function(file) {
      shiny::req(tr()); n <- ape::Ntip(tr())
      grDevices::png(file, width = 2400, height = round(max(1800, 60 * n + 300)), res = 220); on.exit(grDevices::dev.off()); draw(tr(), ttl)
    }
    save_pdf <- function(tr, ttl) function(file) {
      shiny::req(tr()); n <- ape::Ntip(tr())
      grDevices::pdf(file, width = 10, height = max(8, 0.3 * n + 2)); on.exit(grDevices::dev.off()); draw(tr(), ttl)
    }
    output$nj_nwk <- shiny::downloadHandler(function() fn("tree_nj.nwk", "_nj.nwk"), save_nwk(nj_tree))
    output$nj_png <- shiny::downloadHandler(function() fn("tree_nj.png", "_nj.png"), save_png(nj_tree, "Neighbour-joining"))
    output$nj_pdf <- shiny::downloadHandler(function() fn("tree_nj.pdf", "_nj.pdf"), save_pdf(nj_tree, "Neighbour-joining"))
    output$iq_nwk <- shiny::downloadHandler(function() fn("tree_iqtree.nwk", "_iqtree.nwk"), save_nwk(iq_tree))
    output$iq_png <- shiny::downloadHandler(function() fn("tree_iqtree.png", "_iqtree.png"), save_png(iq_tree, "IQ-TREE"))
    output$iq_pdf <- shiny::downloadHandler(function() fn("tree_iqtree.pdf", "_iqtree.pdf"), save_pdf(iq_tree, "IQ-TREE"))
    copy_if <- function(ext) function(file) { shiny::req(rv$iq); f <- paste0(rv$iq$prefix, ext); shiny::req(file.exists(f)); file.copy(f, file, overwrite = TRUE) }
    output$dl_report <- shiny::downloadHandler(function() fn("iqtree_report.iqtree", "_iqtree.iqtree"), copy_if(".iqtree"))
    output$dl_log <- shiny::downloadHandler(function() fn("iqtree_run.log", "_iqtree.log"), copy_if(".log"))
    output$dl_treefile <- shiny::downloadHandler(function() fn("iqtree_run.treefile", "_iqtree.treefile"), copy_if(".treefile"))
    output$dl_dist <- shiny::downloadHandler(function() fn("nj_distances.csv", "_nj_distances.csv"),
      function(file) { shiny::req(rv$nj); utils::write.csv(as.data.frame(rv$nj$D), file) })

    # ---- names and folder chosen by the researcher ----------------------------------------------------------------
    ns <- session$ns
    sname <- function() { x <- input$sv_name; if (is.null(x)) "" else gsub("[^A-Za-z0-9._-]+", "_", trimws(x)) }
    fn <- function(default, suffix) { nm <- sname(); if (nzchar(nm)) paste0(nm, suffix) else default }
    local_mode <- nzchar(Sys.getenv("SANGERLY_DESKTOP")) || interactive()
    output$sv_local_ui <- shiny::renderUI({
      if (!local_mode) return(shiny::helpText("Saving into a folder is available in the desktop version; here use the Download buttons."))
      shiny::tagList(
        shiny::textInput(ns("sv_dir"), "Save into this folder (path on this computer)", "",
                         placeholder = "C:/Users/you/Documents/Sangerly", width = "100%"),
        shiny::actionButton(ns("sv_browse"), "Browse folder..."),
        shiny::checkboxGroupInput(ns("sv_what"), "What to save",
          c("NJ tree: Newick" = "nj_nwk", "NJ tree: PNG" = "nj_png", "NJ tree: PDF" = "nj_pdf",
            "IQ-TREE tree: Newick" = "iq_nwk", "IQ-TREE tree: PNG" = "iq_png", "IQ-TREE tree: PDF" = "iq_pdf",
            "IQ-TREE report (.iqtree)" = "iq_rep", "IQ-TREE log (.log)" = "iq_log", "IQ-TREE treefile (with support values)" = "iq_tf",
            "ALL files of the IQ-TREE run (own sub-folder)" = "iq_all", "ModelFinder table (CSV)" = "mf_csv",
            "Distance matrix (CSV)" = "d_csv", "Distances: Excel (matrix + pairs)" = "d_xlsx",
            "Distance heat map: SVG" = "hm_svg", "Distance heat map: JPG" = "hm_jpg"),
          selected = c("nj_nwk", "nj_png", "iq_nwk", "iq_png", "iq_rep", "iq_tf", "d_xlsx")),
        shiny::actionButton(ns("sv_save"), "Save to this folder", icon = shiny::icon("floppy-disk"), class = "btn-primary"),
        shiny::uiOutput(ns("sv_msg")))
    })
    shiny::observeEvent(input$sv_browse, {
      sg_dirpick(session, function(d) shiny::updateTextInput(session, "sv_dir", value = gsub("\\\\", "/", d)),
                 "Choose the folder where the results are saved", if (is.null(input$sv_dir)) "" else input$sv_dir)
    })
    sv_res <- shiny::reactiveVal(NULL)
    output$sv_msg <- shiny::renderUI({
      x <- sv_res(); if (is.null(x)) return(NULL)
      shiny::div(class = "alert alert-info", style = "margin-top:8px",
        if (length(x$saved)) shiny::tagList(shiny::strong(sprintf("%d item(s) saved:", length(x$saved))), shiny::tags$br(),
                                           paste(x$saved, collapse = "\n")),
        if (length(x$skipped)) shiny::tagList(shiny::tags$br(), shiny::strong("Not saved:"), shiny::tags$br(),
                                              paste(x$skipped, collapse = "\n")))
    })
    shiny::observeEvent(input$sv_save, {
      d <- trimws(if (is.null(input$sv_dir)) "" else input$sv_dir)
      if (!nzchar(d)) { shiny::showNotification("Choose or type a folder first.", type = "warning"); return() }
      if (!dir.exists(d) && !isTRUE(tryCatch(dir.create(d, recursive = TRUE), warning = function(w) FALSE, error = function(e) FALSE))) {
        shiny::showNotification("That folder does not exist and cannot be created.", type = "error"); return()
      }
      what <- input$sv_what
      if (!length(what)) { shiny::showNotification("Tick at least one item to save.", type = "warning"); return() }
      pre <- if (nzchar(sname())) sname() else "tree"
      saved <- character(); skipped <- character()
      path_for <- function(suffix) {        # never overwrite silently: the time is added to the name of an existing file
        f <- file.path(d, paste0(pre, suffix))
        if (file.exists(f)) f <- sub("(\\.[^.]+)$", paste0("_", format(Sys.time(), "%H%M%S"), "\\1"), f)
        f
      }
      add <- function(label, suffix, fun) {
        f <- path_for(suffix)
        r <- tryCatch({ fun(f); NULL }, error = function(e) if (inherits(e, "shiny.silent.error")) "nothing to save yet" else conditionMessage(e))
        if (is.null(r) && file.exists(f)) saved <<- c(saved, normalizePath(f, winslash = "/"))
        else skipped <<- c(skipped, sprintf("%s (%s)", label, if (is.null(r)) "no file produced" else r))
      }
      has_iq <- !is.null(rv$iq) && !is.null(rv$iq$prefix)
      cp <- function(ext) function(f) { shiny::req(has_iq); g <- paste0(rv$iq$prefix, ext); shiny::req(file.exists(g)); file.copy(g, f, overwrite = TRUE) }
      if ("nj_nwk" %in% what) add("NJ tree (Newick)", "_nj.nwk", save_nwk(nj_tree))
      if ("nj_png" %in% what) add("NJ tree (PNG)", "_nj.png", save_png(nj_tree, "Neighbour-joining"))
      if ("nj_pdf" %in% what) add("NJ tree (PDF)", "_nj.pdf", save_pdf(nj_tree, "Neighbour-joining"))
      if ("iq_nwk" %in% what) add("IQ-TREE tree (Newick)", "_iqtree.nwk", save_nwk(iq_tree))
      if ("iq_png" %in% what) add("IQ-TREE tree (PNG)", "_iqtree.png", save_png(iq_tree, "IQ-TREE"))
      if ("iq_pdf" %in% what) add("IQ-TREE tree (PDF)", "_iqtree.pdf", save_pdf(iq_tree, "IQ-TREE"))
      if ("iq_rep" %in% what) add("IQ-TREE report", "_iqtree.iqtree", cp(".iqtree"))
      if ("iq_log" %in% what) add("IQ-TREE log", "_iqtree.log", cp(".log"))
      if ("iq_tf" %in% what)  add("IQ-TREE treefile", "_iqtree.treefile", cp(".treefile"))
      if ("mf_csv" %in% what) add("ModelFinder table", "_modelfinder.csv", function(f) {
        shiny::req(has_iq); tb <- rv$iq$report$table; shiny::req(tb); utils::write.csv(tb, f, row.names = FALSE) })
      if ("d_csv" %in% what)  add("Distance matrix", "_nj_distances.csv", function(f) { shiny::req(rv$nj); utils::write.csv(as.data.frame(rv$nj$D), f) })
      if ("d_xlsx" %in% what) add("Distances (Excel)", "_distances.xlsx", function(f) {
        shiny::req(rv$nj); D <- as.matrix(rv$nj$D)
        .sg_xlsx(list(Matrix = .dist_matrix_df(D), Pairs = .dist_pairs(D)), f) })
      if ("hm_svg" %in% what) add("Heat map (SVG)", "_distance_heatmap.svg", function(f) save_hm_file("svg", f))
      if ("hm_jpg" %in% what) add("Heat map (JPG)", "_distance_heatmap.jpg", function(f) save_hm_file("jpg", f))
      if ("iq_all" %in% what) {
        run_folder <- if (has_iq) dirname(rv$iq$prefix) else NA_character_
        in_tmp <- !is.na(run_folder) && startsWith(normalizePath(run_folder, winslash = "/"), normalizePath(tempdir(), winslash = "/"))
        if (!has_iq) skipped <- c(skipped, "All IQ-TREE files (no IQ-TREE run yet)")
        else if (!in_tmp) saved <- c(saved, paste0(normalizePath(run_folder, winslash = "/"), "  (the IQ-TREE run was already written here)"))
        else {
          src <- list.files(dirname(rv$iq$prefix), full.names = TRUE)
          sub <- file.path(d, paste0(pre, "_iqtree_files"))
          if (dir.exists(sub)) sub <- paste0(sub, "_", format(Sys.time(), "%H%M%S"))
          dir.create(sub, showWarnings = FALSE, recursive = TRUE)
          ok <- file.copy(src, sub, overwrite = TRUE)
          saved <- c(saved, sprintf("%s  (%d files)", normalizePath(sub, winslash = "/"), sum(ok)))
        }
      }
      sv_res(list(saved = saved, skipped = skipped))
    })

    output$msg <- shiny::renderUI({
      if (is.null(rv$msg)) return(NULL)
      shiny::div(class = paste0("alert alert-", rv$msg$type), rv$msg$txt)
    })

    # ---- project file: what is saved / restored (used by the Project tab) -----------------------------
    get_state <- function() shiny::isolate(list(nj = rv$nj, iq = rv$iq))
    set_state <- function(s) { rv$nj <- s$nj; rv$iq <- s$iq; rv$msg <- NULL; invisible() }
    list(get_state = get_state, set_state = set_state, nj = shiny::reactive(rv$nj), iq = shiny::reactive(rv$iq))
  })
}
