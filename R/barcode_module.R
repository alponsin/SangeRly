# Sangerly - "Barcode gap" tab
#   Distances within and between groups (normally species), as used in DNA barcoding:
#   intraspecific distance (inside a group) versus the distance to the nearest other group.
#   A barcode gap exists when the smallest interspecific distance is larger than the largest intraspecific distance.
# The groups are taken from the sequence names (Genus_species after the accession) and can be edited in the table.

# Name of the group of each sequence. Leading accession numbers (MG272989_, NC_012920.1_, ...) are removed first.
.bc_group <- function(nm, rule = "species") {
  x <- sub("^[A-Z]{1,6}_?[0-9]{3,}(\\.[0-9]+)?[_ ]+", "", nm)
  tok <- strsplit(x, "[_ ]+")
  out <- vapply(tok, function(t) {
    t <- t[nzchar(t)]
    if (!length(t)) return("")
    switch(rule, genus = t[1], first = t[1], paste(utils::head(t, 2), collapse = "_"))
  }, "")
  out[!nzchar(out)] <- nm[!nzchar(out)]
  out
}

# Distance matrix of an aligned named character vector. model: "K80" (K2P), "raw" (p-distance), "JC69", "TN93".
.bc_dist <- function(aln, model = "K80") {
  if (!requireNamespace("ape", quietly = TRUE)) stop("The 'ape' package is required: install.packages(\"ape\")")
  m <- do.call(rbind, strsplit(tolower(unname(toupper(aln))), ""))
  rownames(m) <- names(aln)
  D <- ape::dist.dna(ape::as.DNAbin(m), model = model, pairwise.deletion = TRUE, as.matrix = TRUE)
  D[!is.finite(D)] <- NA_real_
  D
}

# D: distance matrix; grp: group of each row of D.
.bc_analyse <- function(D, grp, model = "K80") {
  nm <- rownames(D)
  ix <- which(lower.tri(D), arr.ind = TRUE)
  pairs <- data.frame(Sequence_1 = nm[ix[, 2]], Sequence_2 = nm[ix[, 1]], Group_1 = grp[ix[, 2]], Group_2 = grp[ix[, 1]],
                      Distance = D[ix], stringsAsFactors = FALSE)
  pairs$Type <- ifelse(pairs$Group_1 == pairs$Group_2, "Intraspecific", "Interspecific")
  gl <- unique(grp)
  rows <- lapply(gl, function(g) {
    n <- sum(grp == g)
    intra <- pairs$Distance[pairs$Group_1 == g & pairs$Group_2 == g]
    inter <- pairs[xor(pairs$Group_1 == g, pairs$Group_2 == g) & !is.na(pairs$Distance), , drop = FALSE]
    other <- ifelse(inter$Group_1 == g, inter$Group_2, inter$Group_1)
    k <- if (nrow(inter)) which.min(inter$Distance) else NA_integer_
    mx <- if (any(!is.na(intra))) max(intra, na.rm = TRUE) else NA_real_
    mn <- if (!is.na(k)) inter$Distance[k] else NA_real_
    data.frame(Group = g, N = n,
               Max_intra = mx,
               Mean_intra = if (any(!is.na(intra))) mean(intra, na.rm = TRUE) else NA_real_,
               Nearest_group = if (!is.na(k)) other[k] else NA_character_,
               Min_inter = mn,
               Gap = if (is.na(mx) || is.na(mn)) NA_real_ else mn - mx,
               Verdict = if (is.na(mn)) "no other group"
                         else if (is.na(mx)) "single sequence (no intraspecific distance)"
                         else if (mn > mx) "barcode gap" else "overlap (no gap)",
               stringsAsFactors = FALSE)
  })
  groups <- do.call(rbind, rows)
  intra_all <- pairs$Distance[pairs$Type == "Intraspecific"]; inter_all <- pairs$Distance[pairs$Type == "Interspecific"]
  rng <- function(f, x) if (any(!is.na(x))) f(x, na.rm = TRUE) else NA_real_
  summ <- list(model = model, n_seq = length(nm), n_groups = length(gl),
               n_intra = sum(!is.na(intra_all)), n_inter = sum(!is.na(inter_all)),
               mean_intra = rng(mean, intra_all), max_intra = rng(max, intra_all),
               mean_inter = rng(mean, inter_all), min_inter = rng(min, inter_all))
  summ$gap <- if (is.na(summ$max_intra) || is.na(summ$min_inter)) NA_real_ else summ$min_inter - summ$max_intra
  summ$n_with_gap <- sum(groups$Verdict == "barcode gap")
  summ$n_testable <- sum(groups$Verdict %in% c("barcode gap", "overlap (no gap)"))
  list(pairs = pairs, groups = groups, summary = summ)
}

barcodeUI <- function(id) {
  ns <- shiny::NS(id)
  shiny::div(class = "container-fluid",
    shiny::fluidRow(
      shiny::column(4,
        .card("1. Alignment and groups",
          shiny::radioButtons(ns("src"), NULL, c("MSA tab: kept columns (trimmed)" = "trim", "MSA tab: whole alignment" = "all")),
          shiny::uiOutput(ns("src_info")),
          shiny::selectInput(ns("rule"), "Group = ",
            c("Genus_species (the name after the accession)" = "species", "Genus only" = "genus", "First word of the name" = "first")),
          shiny::actionButton(ns("detect"), "Detect groups from the names", icon = shiny::icon("wand-magic-sparkles")),
          shiny::helpText("Edit the Group column by double-clicking a cell. Sequences with the same group name are treated as one species."),
          DT::DTOutput(ns("tbl"))),
        .card("2. Distance",
          shiny::selectInput(ns("model"), "Model",
            c("Kimura 2-parameter (K2P, standard in barcoding)" = "K80", "p-distance (raw)" = "raw", "Jukes-Cantor (JC69)" = "JC69", "Tamura-Nei (TN93)" = "TN93")),
          shiny::helpText("Pairwise deletion of gaps and missing data."),
          shiny::actionButton(ns("run"), "Calculate", icon = shiny::icon("play"), class = "btn-primary"))),
      shiny::column(8,
        shiny::uiOutput(ns("msg")),
        shiny::uiOutput(ns("boxes")),
        shiny::tabsetPanel(
          shiny::tabPanel("Histogram", shiny::br(), shiny::plotOutput(ns("p_hist"), height = "380px"),
            shiny::helpText("Blue = distances between sequences of the same group, red = between different groups. ",
                            "Dashed lines: the largest intraspecific and the smallest interspecific distance. A gap exists when the red part starts to the right of the blue part.")),
          shiny::tabPanel("Per group", shiny::br(), DT::DTOutput(ns("tbl_groups")),
            shiny::helpText("Max_intra = largest distance inside the group. Min_inter = distance to the nearest sequence of another group. ",
                            "Gap = Min_inter - Max_intra (positive: the group is separated by a barcode gap).")),
          shiny::tabPanel("Max intra vs nearest inter", shiny::br(), shiny::plotOutput(ns("p_scatter"), height = "460px"),
            shiny::helpText("Each point is a group. Points above the diagonal line have a barcode gap; points on or below it overlap with another group. ",
                            "Groups with one sequence have no intraspecific distance and are not drawn.")),
          shiny::tabPanel("Pairs", shiny::br(), DT::DTOutput(ns("tbl_pairs")))),
        shiny::br(),
        shiny::div(style = "display:flex;gap:6px;flex-wrap:wrap",
          shiny::downloadButton(ns("dl_xlsx"), "Excel (summary, groups, pairs)"),
          shiny::downloadButton(ns("dl_groups"), "Groups (CSV)"),
          shiny::downloadButton(ns("dl_pairs"), "Pairs (CSV)")))))
}

barcodeServer <- function(id, aln_trim = shiny::reactive(NULL), aln_all = shiny::reactive(NULL)) {
  shiny::moduleServer(id, function(input, output, session) {
    rv <- shiny::reactiveValues(groups = NULL, gver = 0L, res = NULL, msg = NULL)
    set_msg <- function(txt, type = "info") rv$msg <- if (is.null(txt)) NULL else list(txt = txt, type = type)

    src <- shiny::reactive({
      a <- if (identical(input$src, "all")) aln_all() else aln_trim()
      if (is.null(a) || length(a) < 3) NULL else a
    })
    output$src_info <- shiny::renderUI({
      a <- src()
      if (is.null(a)) shiny::helpText("No alignment yet (at least 3 sequences). Align in the MSA tab first.")
      else shiny::helpText(sprintf("%d sequences, %d columns.", length(a), nchar(a[[1]])))
    })

    detect <- function(keep = FALSE) {
      a <- src(); if (is.null(a)) return(invisible())
      g <- .bc_group(names(a), .nz(input$rule, "species"))
      old <- rv$groups
      if (keep && !is.null(old)) {                       # keep what the user already edited
        m <- match(names(a), old$Sequence); g[!is.na(m)] <- old$Group[m[!is.na(m)]]
      }
      rv$groups <- data.frame(Sequence = names(a), Group = g, stringsAsFactors = FALSE)
      rv$gver <- rv$gver + 1L
    }
    shiny::observeEvent(input$detect, detect(FALSE))
    shiny::observeEvent(src(), {
      if (is.null(src())) { rv$groups <- NULL; rv$res <- NULL; rv$gver <- rv$gver + 1L }
      else detect(TRUE)
    }, ignoreNULL = FALSE)

    output$tbl <- DT::renderDT({
      rv$gver
      g <- shiny::isolate(rv$groups); shiny::req(g)
      DT::datatable(g, rownames = FALSE, selection = "none",
                    editable = list(target = "cell", disable = list(columns = 0)),
                    options = list(pageLength = 10, scrollX = TRUE, dom = "ftp"))
    })
    shiny::observeEvent(input$tbl_cell_edit, {
      rv$groups <- DT::editData(rv$groups, input$tbl_cell_edit, rownames = FALSE)
    })

    shiny::observeEvent(input$run, {
      a <- src(); g <- rv$groups
      if (is.null(a) || is.null(g)) { set_msg("Provide an alignment with at least 3 sequences first.", "warning"); return() }
      grp <- trimws(g$Group[match(names(a), g$Sequence)])
      grp[is.na(grp) | !nzchar(grp)] <- "(unassigned)"
      if (length(unique(grp)) < 2) { set_msg("There is only one group: the groups must differ to compare within and between them.", "warning"); return() }
      D <- tryCatch(.bc_dist(a, .nz(input$model, "K80")), error = function(e) { set_msg(conditionMessage(e), "danger"); NULL })
      if (is.null(D)) return()
      rv$res <- .bc_analyse(D, grp, .nz(input$model, "K80"))
      set_msg(NULL)
    })

    output$msg <- shiny::renderUI({
      if (is.null(rv$msg)) return(NULL)
      shiny::div(class = paste0("alert alert-", rv$msg$type), rv$msg$txt)
    })
    f3 <- function(x) if (is.null(x) || is.na(x)) "-" else formatC(x, format = "f", digits = 4)
    output$boxes <- shiny::renderUI({
      r <- rv$res; shiny::req(r); s <- r$summary
      verdict <- if (is.na(s$gap)) "-" else if (s$gap > 0) "Yes" else "No"
      shiny::tagList(
        shiny::fluidRow(
          shiny::column(3, .infobox("layer-group", "Groups", sprintf("%d (%d sequences)", s$n_groups, s$n_seq), "#2c7be5")),
          shiny::column(3, .infobox("arrows-left-right", "Mean intra / inter", paste(f3(s$mean_intra), "/", f3(s$mean_inter)), "#6f42c1")),
          shiny::column(3, .infobox("ruler", "Max intra | Min inter", paste(f3(s$max_intra), "|", f3(s$min_inter)), "#e67e22")),
          shiny::column(3, .infobox("circle-check", "Barcode gap overall", verdict, if (identical(verdict, "Yes")) "#27ae60" else "#c0392b"))),
        shiny::helpText(sprintf("Model %s. Groups with a gap: %d of %d groups that can be tested.", s$model, s$n_with_gap, s$n_testable)))
    })

    output$p_hist <- shiny::renderPlot({
      r <- rv$res; shiny::req(r)
      d <- r$pairs[!is.na(r$pairs$Distance), , drop = FALSE]; shiny::req(nrow(d) > 0)
      d$Type <- factor(d$Type, levels = c("Intraspecific", "Interspecific"))
      s <- r$summary
      p <- ggplot2::ggplot(d, ggplot2::aes(Distance, fill = Type)) +
        ggplot2::geom_histogram(bins = 40, alpha = 0.7, position = "identity", colour = "white", linewidth = 0.2) +
        ggplot2::scale_fill_manual(values = c(Intraspecific = "#2c7be5", Interspecific = "#e5484d"), name = NULL, drop = FALSE) +
        ggplot2::labs(x = paste("Distance (", s$model, ")", sep = ""), y = "Number of pairs") +
        ggplot2::theme_minimal(base_size = 13) + ggplot2::theme(legend.position = "top")
      v <- c(s$max_intra, s$min_inter); v <- v[!is.na(v)]
      if (length(v)) p <- p + ggplot2::geom_vline(xintercept = v, linetype = 2, colour = "grey30")
      p
    })
    output$tbl_groups <- DT::renderDT({
      r <- rv$res; shiny::req(r)
      d <- r$groups
      for (k in c("Max_intra", "Mean_intra", "Min_inter", "Gap")) d[[k]] <- round(d[[k]], 4)
      DT::datatable(d, rownames = FALSE, options = list(pageLength = 15, scrollX = TRUE, dom = "ftp"))
    })
    output$p_scatter <- shiny::renderPlot({
      r <- rv$res; shiny::req(r)
      d <- r$groups[!is.na(r$groups$Max_intra) & !is.na(r$groups$Min_inter), , drop = FALSE]
      shiny::validate(shiny::need(nrow(d) > 0, "No group has two or more sequences with another group to compare with."))
      d$Result <- ifelse(d$Gap > 0, "barcode gap", "overlap")
      lim <- c(0, max(d$Max_intra, d$Min_inter) * 1.08)
      ggplot2::ggplot(d, ggplot2::aes(Max_intra, Min_inter, colour = Result, label = Group)) +
        ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey40") +
        ggplot2::geom_point(size = 3) +
        ggplot2::geom_text(vjust = -0.8, size = 3.3, check_overlap = TRUE, show.legend = FALSE) +
        ggplot2::scale_colour_manual(values = c("barcode gap" = "#27ae60", "overlap" = "#c0392b"), name = NULL) +
        ggplot2::coord_cartesian(xlim = lim, ylim = lim) +
        ggplot2::labs(x = "Largest intraspecific distance", y = "Distance to the nearest other group") +
        ggplot2::theme_minimal(base_size = 13) + ggplot2::theme(legend.position = "top")
    })
    output$tbl_pairs <- DT::renderDT({
      r <- rv$res; shiny::req(r)
      d <- r$pairs; d$Distance <- round(d$Distance, 5)
      DT::datatable(d, rownames = FALSE, options = list(pageLength = 15, scrollX = TRUE, dom = "ftp"))
    })

    summ_df <- function(r) {
      s <- r$summary
      data.frame(Item = c("Model", "Sequences", "Groups", "Intraspecific pairs", "Interspecific pairs", "Mean intraspecific", "Largest intraspecific",
                          "Mean interspecific", "Smallest interspecific", "Gap (smallest inter - largest intra)", "Groups with a gap"),
                 Value = c(s$model, s$n_seq, s$n_groups, s$n_intra, s$n_inter, f3(s$mean_intra), f3(s$max_intra), f3(s$mean_inter),
                           f3(s$min_inter), f3(s$gap), sprintf("%d of %d", s$n_with_gap, s$n_testable)), stringsAsFactors = FALSE)
    }
    output$dl_groups <- shiny::downloadHandler(function() "barcode_gap_groups.csv",
      function(file) { shiny::req(rv$res); utils::write.csv(rv$res$groups, file, row.names = FALSE) })
    output$dl_pairs <- shiny::downloadHandler(function() "barcode_gap_pairs.csv",
      function(file) { shiny::req(rv$res); utils::write.csv(rv$res$pairs, file, row.names = FALSE) })
    output$dl_xlsx <- shiny::downloadHandler(function() "barcode_gap.xlsx",
      function(file) {
        shiny::req(rv$res)
        tryCatch(.sg_xlsx(list(Summary = summ_df(rv$res), Groups = rv$res$groups, Pairs = rv$res$pairs), file),
          error = function(e) {
            shiny::showNotification(paste("Excel file not written:", conditionMessage(e)), type = "error", duration = 15)
            writeLines(conditionMessage(e), file)
          })
      })

    get_state <- function() shiny::isolate(list(groups = rv$groups, res = rv$res))
    set_state <- function(s) { rv$groups <- s$groups; rv$res <- s$res; rv$gver <- rv$gver + 1L; invisible() }
    list(get_state = get_state, set_state = set_state, result = shiny::reactive(rv$res))
  })
}
