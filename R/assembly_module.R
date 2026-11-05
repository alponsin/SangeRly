# ==============================================================================
# Sangerly - Viewer + Shiny assembly module (v5)
# Needs (sourced first): assembly_core.R, pileup_core.R, extras_core.R
# Packages: shiny, ggplot2, DT, httr, jsonlite, Biostrings
# ==============================================================================

suppressPackageStartupMessages({ library(shiny); library(ggplot2) })

BASE_COLS  <- c(A = "#5fd39a", C = "#58a9e0", G = "#dcdcdc", T = "#ff7b7b", N = "#9aa0a6", `-` = "#ffffff")
TRACE_COLS <- c(A = "#1a9850", C = "#2166ac", G = "#222222", T = "#d73027")


sangerly_css <- function() tags$style(HTML("
  body{background:#f4f6f9}
  .navbar{box-shadow:0 2px 6px rgba(0,0,0,.15)}
  .alert{padding:10px 14px;margin-bottom:10px;border-radius:6px}
  .alert-info{background:#e8f1fb;border:1px solid #c9ddf3}
  .alert-danger{background:#fde8e8;border:1px solid #f3c9c9}
  .alert-warn{background:#fff6dd;border:1px solid #f1e0a6}
  .legend-chip{display:inline-block;padding:1px 8px;margin:0 6px 4px 0;border-radius:3px;font-size:12px}
  .card{border:0;box-shadow:0 1px 3px rgba(0,0,0,.12);border-radius:8px}
  .card-header{background:#fff;border-bottom:2px solid #2c7be5;color:#1f3c5f}
  .ibox{display:flex;align-items:stretch;background:#fff;border-radius:8px;box-shadow:0 1px 3px rgba(0,0,0,.12);margin-bottom:12px;overflow:hidden}
  .ibox-ico{width:58px;display:flex;align-items:center;justify-content:center;color:#fff;font-size:22px}
  .ibox-txt{padding:8px 12px}
  .ibox-lab{font-size:12px;color:#6b7785;text-transform:uppercase;letter-spacing:.4px}
  .ibox-val{font-size:22px;font-weight:600;color:#1f3c5f;line-height:1.2}
  .sidebar-sec{background:#fff;border-radius:8px;box-shadow:0 1px 3px rgba(0,0,0,.12);padding:10px 14px;margin-bottom:12px}
  .sidebar-sec summary{font-weight:600;color:#1f3c5f;cursor:pointer;margin-bottom:6px}
  .step{display:flex;align-items:flex-start;margin-bottom:10px}
  .step-n{background:#2c7be5;color:#fff;border-radius:50%;width:26px;height:26px;display:flex;align-items:center;justify-content:center;font-weight:600;margin-right:10px;flex:none}
  .dataTables_wrapper{font-size:13px}
  pre.sx{font-family:Consolas,Menlo,monospace;font-size:12px;line-height:1.35;background:#fff;border:1px solid #e3e8ee;
         border-radius:6px;padding:10px;overflow-x:auto;white-space:pre}
  .sx1{background:#fbd5d0}
  .sx2{background:#74c69d;font-weight:700}
  .sx3{color:#1f5fa8;font-weight:700}
"))

.card <- function(title, ..., class = "") {
  div(class = paste("card mb-3", class),
      div(class = "card-header fw-semibold", title),
      div(class = "card-body", ...))
}
.infobox <- function(ico, label, value, color = "#2c7be5") {
  div(class = "ibox",
      div(class = "ibox-ico", style = paste0("background:", color), icon(ico)),
      div(class = "ibox-txt", div(class = "ibox-lab", label), div(class = "ibox-val", value)))
}


# Trimming quality plot (green = raw read, orange = region kept), as in sangeranalyseR
plot_phred <- function(rd) {
  q <- rd$raw$qual; n <- length(q); f <- rd$qc$from; t <- rd$qc$to
  df <- data.frame(i = seq_len(n), q = q, used = seq_len(n) >= f & seq_len(n) <= t)
  ggplot(df, aes(i, q, fill = used)) + geom_col(width = 1) +
    scale_fill_manual(values = c(`TRUE` = "#f4a261", `FALSE` = "#cfd8dc"), guide = "none") +
    geom_hline(yintercept = c(20, 30), linetype = 3, colour = "grey40") +
    annotate("rect", xmin = .5, xmax = n + .5, ymin = -6, ymax = -3, fill = "#2a9d8f") +
    annotate("rect", xmin = f - .5, xmax = t + .5, ymin = -11, ymax = -8, fill = "#f4a261") +
    labs(x = "base index (raw read)", y = "Phred",
         subtitle = "green = raw read | orange = region kept after trimming | dotted lines = Q20 and Q30") +
    theme_minimal(base_size = 11) + theme(panel.grid.minor = element_blank())
}

# Chromatogram of one read (original sequencer orientation), several rows, secondary calls
plot_read_chromatogram <- function(rd, start = 1, nper = 25, nrows = 3, ratio = 0.33, show_trim = TRUE) {
  rs <- rd$raw$seq; rq <- rd$raw$qual; rp <- rd$raw$ploc; n <- length(rs)
  tr <- rd$traces
  if (isTRUE(rd$rc)) colnames(tr) <- unname(COMP_CHAR[colnames(tr)])
  start <- max(1L, min(as.integer(start), n))
  L <- list(); LB <- list(); TR <- list(); strips <- character()
  for (r in seq_len(nrows)) {
    a <- start + (r - 1L) * nper
    if (a > n) break
    b <- min(n, a + nper - 1L); ii <- a:b
    sp <- if (length(ii) > 1) stats::median(diff(rp[ii])) else 10
    if (!is.finite(sp) || sp <= 0) sp <- 10
    sc <- floor(rp[a] - sp / 2):ceiling(rp[b] + sp / 2)
    sc <- sc[sc >= 1 & sc <= nrow(tr)]
    mx <- max(tr[sc, , drop = FALSE]); if (!is.finite(mx) || mx <= 0) mx <- 1
    strip <- sprintf("bases %d-%d", a, b)
    strips <- c(strips, strip)
    for (ch in colnames(tr))
      L[[length(L) + 1]] <- data.frame(x = sc, y = tr[sc, ch] / mx, ch = ch, strip = strip, stringsAsFactors = FALSE)
    sec <- rep("", length(ii))
    for (j in seq_along(ii)) {
      p0 <- min(max(rp[ii[j]], 1L), nrow(tr))
      v <- apply(tr[max(1, p0 - 1):min(nrow(tr), p0 + 1), , drop = FALSE], 2, max)
      o <- order(v, decreasing = TRUE)
      if (v[o[1]] > 0 && v[o[2]] / v[o[1]] >= ratio) sec[j] <- names(v)[o[2]]
    }
    LB[[length(LB) + 1]] <- data.frame(x = rp[ii], b = rs[ii], q = rq[ii], sec = sec, strip = strip, stringsAsFactors = FALSE)
    if (show_trim) {
      out <- ii < rd$qc$from | ii > rd$qc$to
      if (any(out))
        TR[[length(TR) + 1]] <- data.frame(xmin = rp[ii[out]] - sp / 2, xmax = rp[ii[out]] + sp / 2, strip = strip, stringsAsFactors = FALSE)
    }
  }
  if (!length(L)) return(NULL)
  L <- do.call(rbind, L); LB <- do.call(rbind, LB)
  L$strip <- factor(L$strip, levels = strips); LB$strip <- factor(LB$strip, levels = strips)
  p <- ggplot()
  if (length(TR)) {
    TR <- do.call(rbind, TR); TR$strip <- factor(TR$strip, levels = strips)
    p <- p + geom_rect(data = TR, aes(xmin = xmin, xmax = xmax), ymin = -0.15, ymax = 1.3, fill = "#e74c3c", alpha = .18)
  }
  p + geom_line(data = L, aes(x = x, y = y, colour = ch), linewidth = .5) +
    scale_colour_manual(values = TRACE_COLS, name = "channel") +
    geom_text(data = LB, aes(x = x, y = 1.2, label = b, colour = b), size = 3.8, fontface = "bold", show.legend = FALSE) +
    geom_text(data = LB[LB$sec != "", , drop = FALSE], aes(x = x, y = 1.08, label = sec), size = 2.9, colour = "#7f8c8d") +
    geom_text(data = LB, aes(x = x, y = -0.08, label = q), size = 2.2, colour = "#999999") +
    facet_wrap(~strip, ncol = 1, scales = "free_x") +
    coord_cartesian(ylim = c(-0.15, 1.3)) +
    theme_minimal(base_size = 10) +
    theme(axis.title = element_blank(), axis.text = element_blank(), panel.grid = element_blank(),
          legend.position = "bottom", strip.text = element_text(hjust = 0, face = "bold"))
}

plot_dist_heat <- function(D) {
  df <- as.data.frame(as.table(D)); names(df) <- c("a", "b", "d")
  ggplot(df, aes(a, b, fill = d)) + geom_tile(colour = "white") +
    geom_text(aes(label = ifelse(is.na(d), "", sprintf("%.3f", d))), size = 3) +
    scale_fill_gradient(low = "#f7fbff", high = "#2171b5", na.value = "#eeeeee", name = "p-distance") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 40, hjust = 1), axis.title = element_blank(), panel.grid = element_blank())
}

# Consensus as one or several contigs (split at coverage gaps).
# Positions deleted by the user (Type "ColDel") ALWAYS vanish from the consensus.
# Real gaps ("-" from the reads / reference alignment or set manually) stay as "-" only if keep_gaps = TRUE.
consensus_contigs <- function(cf, split = TRUE, keep_gaps = FALSE) {
  b <- cf$Final; nocov <- cf$Type == "NoCov"
  b[cf$Type == "ColDel"] <- ""
  rm_gap <- function(x) if (keep_gaps) x else gsub("-", "", x)
  if (!split || !any(nocov)) return(list(consensus = rm_gap(paste(b, collapse = ""))))
  r <- rle(nocov); ends <- cumsum(r$lengths); starts <- ends - r$lengths + 1
  out <- list()
  for (j in which(!r$values)) {
    nm <- sprintf("contig%d_col%d-%d", length(out) + 1, cf$Col[starts[j]], cf$Col[ends[j]])
    if (!all(is.na(cf$RefPos[starts[j]:ends[j]])))
      nm <- sprintf("%s_ref%d-%d", nm, cf$RefPos[starts[j]], cf$RefPos[ends[j]])
    seq1 <- rm_gap(paste(b[starts[j]:ends[j]], collapse = ""))
    if (nzchar(gsub("-", "", seq1))) out[[nm]] <- seq1
  }
  if (!length(out)) out[["consensus"]] <- ""
  out
}

# Overview of the six frames: ORFs (red), start Met (green dot), stop codons (blue ticks)
plot_frames <- function(six) {
  fn <- vapply(six$frames, function(f) f$name, "")
  yv <- setNames(seq(6, 1), fn)
  o <- six$orfs
  stp <- do.call(rbind, lapply(six$frames, function(f) {
    x <- f$mid[f$aa == "*"]
    if (length(x)) data.frame(frame = f$name, x = x, stringsAsFactors = FALSE)
  }))
  p <- ggplot() + geom_hline(yintercept = unname(yv), colour = "#e3e8ee")
  if (nrow(o)) {
    y <- unname(yv[o$Frame])
    rc <- data.frame(xmin = o$Start_nt - .5, xmax = o$End_nt + .5, y = y)
    sx <- ifelse(startsWith(o$Frame, "5"), o$Start_nt + 1, o$End_nt - 1)
    p <- p + geom_rect(data = rc, aes(xmin = xmin, xmax = xmax, ymin = y - .3, ymax = y + .3), fill = "#e74c3c", alpha = .4) +
      geom_point(data = data.frame(x = sx, y = y), aes(x = x, y = y), colour = "#27ae60", size = 2)
  }
  if (!is.null(stp)) {
    stp$y <- unname(yv[stp$frame])
    p <- p + geom_segment(data = stp, aes(x = x, xend = x, y = y - .3, yend = y + .3), colour = "#1f5fa8", linewidth = .3)
  }
  p + scale_y_continuous(breaks = unname(yv), labels = names(yv), limits = c(0.5, 6.5)) +
    scale_x_continuous(limits = c(1, max(2, six$n)), expand = c(0.01, 0)) +
    labs(x = "position in the translated sequence (nt)", y = NULL,
         subtitle = "red = ORF (Met to stop) | green dot = start Met | blue ticks = stop codons") +
    theme_minimal(base_size = 11) + theme(panel.grid = element_blank())
}

.sv_local <- function() nzchar(Sys.getenv("SANGERLY_DESKTOP")) || interactive()

fasta_of <- function(lst) paste0(">", names(lst), "\n", unlist(lst), collapse = "\n")

# ------------------------------------------------------------------------------
# UI
# ------------------------------------------------------------------------------
assemblyUI <- function(id) {
  ns <- NS(id)
  sidebarLayout(
    sidebarPanel(
      width = 3,
      div(class = "sidebar-sec",
        tags$details(open = NA, tags$summary("1. Input"),
          radioButtons(ns("mode"), NULL,
            c("De novo: Forward + Reverse overlap (no reference)" = "denovo",
              "Reference-guided: many reads / tiled regions" = "ref")),
          fileInput(ns("reads"), "Add .ab1 files (pick as many times as you like, from any folder)",
                    accept = ".ab1", multiple = TRUE, buttonLabel = "Add files...", placeholder = "no file chosen"),
          if (.sv_local()) tagList(
            textInput(ns("dir_path"), "Or add a whole folder (path on this computer)", "", placeholder = "C:/data/run1"),
            div(style = "display:flex;gap:10px;align-items:center;margin-bottom:8px",
                actionButton(ns("dir_add"), "Add folder"),
                checkboxInput(ns("dir_rec"), "with subfolders", TRUE)),
            if (.Platform$OS.type == "windows")
              div(style = "display:flex;gap:8px;flex-wrap:wrap", class = "mb-2",
                  actionButton(ns("browse_add"), "Browse files..."),
                  actionButton(ns("browse_dir"), "Browse folder..."))),
          uiOutput(ns("files_count")),
          div(style = "display:flex;gap:6px;margin-bottom:8px",
              actionButton(ns("files_rm"), "Remove selected", icon = icon("minus")),
              actionButton(ns("files_clear"), "Clear all")),
          conditionalPanel(sprintf("input['%s'] == 'ref'", ns("mode")),
            fileInput(ns("ref_file"), "Reference region(s) (FASTA; several records allowed)", accept = c(".fa", ".fasta", ".fas", ".txt")),
            fluidRow(column(8, textInput(ns("ref_acc"), "or NCBI accession", "")),
                     column(4, br(), actionButton(ns("ref_fetch"), "Fetch"))),
            textOutput(ns("ref_info"))))),
      div(class = "sidebar-sec",
        tags$details(open = NA, tags$summary("2. Trimming and thresholds"),
          radioButtons(ns("trim_method"), "Trimming method",
            c("M2: sliding window (Trimmomatic-like)" = "M2", "M1: modified Mott (Phred/Phrap-like)" = "M1")),
          conditionalPanel(sprintf("input['%s'] == 'M2'", ns("trim_method")),
            numericInput(ns("q_trim"), "Window mean Phred >=", 20, 0, 60),
            numericInput(ns("win"), "Window size", 10, 2, 40)),
          conditionalPanel(sprintf("input['%s'] == 'M1'", ns("trim_method")),
            numericInput(ns("mott_cutoff"), "Mott cutoff (0-1)", 0.0001, 0, 1, step = 0.0001),
            helpText("The sangeranalyseR default (0.0001) is very strict: only bases with Phred > 40 score positive. ",
                     "Raise it (e.g. 0.001 - 0.01) if reads are trimmed too short.")),
          numericInput(ns("q_good"), "'Good' Phred >=", 30, 10, 60),
          numericInput(ns("q_poor"), "'Poor' Phred <", 20, 5, 50),
          numericInput(ns("ratio"), "Signal ratio cutoff (double peak)", 0.33, 0.05, 1, step = 0.01),
          selectInput(ns("type4"), "Good reads disagree (type 4)",
                      c("Mark N" = "N", "IUPAC code" = "iupac", "Follow best-supported base" = "F")),
          checkboxInput(ns("auto_del"), "Delete bases below a Phred threshold while assembling", FALSE),
          numericInput(ns("auto_q"), "Auto-delete threshold (Phred <)", 15, 0, 60))),
      actionButton(ns("go"), "ASSEMBLE / RECALCULATE", icon = icon("play"), class = "btn-danger w-100 mb-3"),
      div(class = "sidebar-sec",
        tags$details(tags$summary("3. Session"),
          downloadButton(ns("dl_sess"), "Save session (.rds)"),
          fileInput(ns("sess_file"), "Load a Sangerly session (.rds)", accept = ".rds"),
          helpText("Only load session files that you saved yourself from Sangerly.")))
    ),
    mainPanel(
      width = 9,
      uiOutput(ns("summary")),
      uiOutput(ns("sample_bar")),
      tabsetPanel(
        id = ns("tabs"),
        tabPanel("Input files", value = "files", br(),
          .card("Files to process",
            fluidRow(
              column(5, radioButtons(ns("group_mode"), "How files are grouped into samples",
                c("By file name (the F / R tag is removed: S1_F + S1_R = sample S1)" = "name",
                  "By file name, ignoring long run numbers (1st BASE / Macrogen IDs)" = "name_norun",
                  "All files are one sample" = "all"), selected = "name_norun")),
              column(7, helpText("Files accumulate: use the Add files button as often as you like, also with files from different folders. ",
                                 "The Sample column can be edited by double-clicking a cell; files with the same sample name are assembled together, ",
                                 "and all samples are processed in one run with ASSEMBLE. Changing the grouping rule overwrites manual edits."))),
            DT::DTOutput(ns("files_table")))),
        tabPanel("Samples", value = "smp", br(),
          uiOutput(ns("smp_head")),
          DT::DTOutput(ns("smp_table")),
          br(),
          downloadButton(ns("dl_all_fasta"), "All consensus sequences (FASTA)"),
          downloadButton(ns("dl_all_csv"), "Summary (CSV)"),
          helpText("Click a row to open that sample in the Alignment tab. The crop settings of the Consensus and export tab apply to every sample.")),
        tabPanel("Overview", value = "ring", br(), uiOutput(ns("boxes")), uiOutput(ns("welcome"))),

        tabPanel("Alignment", value = "aln", br(),
          fluidRow(
            column(4, tags$label("Jump to base number"),
                   div(style = "display:flex;gap:6px",
                       numericInput(ns("jump_pos"), NULL, value = NA, min = 1, step = 1, width = "120px"),
                       actionButton(ns("jump_go"), "Go", icon = icon("magnifying-glass"))),
                   uiOutput(ns("jump_hint"))),
            column(7, br(), actionButton(ns("prev"), "< Needs review"), actionButton(ns("nxt"), "Needs review >"),
                   actionButton(ns("snp_nxt"), "SNP >"))),
          div(tags$span(class = "legend-chip", style = "background:#f7e08a", "needs review"),
              tags$span(class = "legend-chip", style = "background:#bfe6cf", "decided manually"),
              tags$span(class = "legend-chip", style = "background:#9aa5b1", "position deleted"),
              tags$span(class = "legend-chip", style = "border:2px solid #c0392b", "read differs from consensus"),
              tags$span(class = "legend-chip", style = "background:#eee", "x = base deleted"),
              tags$span(class = "legend-chip", style = "background:#eee", "base transparency = Phred")),
          div(class = "text-muted small mb-2",
              "Drag the alignment to move along it. Click a base to select it; click the Consensus row to select the whole position. ",
              "Shift+click or Shift+arrows extend the selection. ", tags$kbd("Delete"), " / ", tags$kbd("Backspace"),
              " removes the selection: on the Consensus row the whole position (all reads and the consensus), on a read row only that read's bases. ",
              tags$kbd("Ctrl+Z"), " undoes. Ctrl+wheel zooms. ",
              tags$b("Edit the consensus from the keyboard: "), "select a position (or a range) and type ",
              tags$kbd("A"), tags$kbd("C"), tags$kbd("G"), tags$kbd("T"), tags$kbd("N"), " (also IUPAC letters, or ", tags$kbd("-"),
              ") to set the consensus base; ", tags$kbd("="), " (or ", tags$kbd("."), " / ", tags$kbd("Enter"), ") uses the reference base; ",
              tags$kbd("0"), " (or ", tags$kbd("Esc"), ") returns to the rule decision. The selection moves on to the next position, so you can keep typing. ",
              "Typing always edits the consensus, whichever row is selected. Click the viewer first so that it receives the keys."),
          svViewerUI(ns("viewer")),
          .card("Delete low-Phred bases",
            fluidRow(
              column(6,
                numericInput(ns("del_q"), "Delete bases with Phred <", 20, 0, 60),
                radioButtons(ns("del_scope"), "What happens to the consensus",
                  c("Reads only: recompute the consensus from the remaining reads" = "reads",
                    "Reads + consensus position when no usable base remains" = "orphans",
                    "Reads + consensus position wherever a low base occurred" = "columns"),
                  selected = "orphans"),
                actionButton(ns("del_apply"), "Delete", icon = icon("scissors")),
                actionButton(ns("del_restore_all"), "Restore everything")),
              column(6,
                strong("Delete a stretch of positions (reads and consensus)"),
                fluidRow(column(6, numericInput(ns("col_from"), "From", NA, 1)),
                         column(6, numericInput(ns("col_to"), "To", NA, 1))),
                uiOutput(ns("col_hint")),
                actionButton(ns("col_del_go"), "Delete these positions", icon = icon("scissors")),
                hr(),
                strong("Delete a stretch in ONE read only"),
                fluidRow(column(5, selectInput(ns("rng_read"), "Read", choices = NULL)),
                         column(3, numericInput(ns("rng_from"), "From col", 1, 1)),
                         column(3, numericInput(ns("rng_to"), "To col", 1, 1))),
                actionButton(ns("rng_go"), "Delete in this read"))),
            uiOutput(ns("del_count")),
            helpText("Deleted bases are ignored by decisions, the consensus, SNPs and genetic distances. ",
                     "A deleted position is removed from the consensus as well (the consensus gets shorter), ",
                     "not replaced by N. Everything can be restored. For a single base, click its column and use the panel below.")),
          uiOutput(ns("detail"))),

        tabPanel("Read and QC", value = "read", br(),
          fluidRow(
            column(4, selectInput(ns("qc_read"), "Select read", choices = NULL)),
            column(3, br(), uiOutput(ns("qc_toggle"))),
            column(5, uiOutput(ns("qc_msg")))),
          uiOutput(ns("qc_boxes")),
          .card("Quality and trimming", plotOutput(ns("phred_plot"), height = 260)),
          .card("Read chromatogram (original orientation)",
            fluidRow(
              column(2, numericInput(ns("chr_nper"), "Bases per row", 25, 10, 60)),
              column(2, numericInput(ns("chr_nrows"), "Rows", 3, 1, 6)),
              column(2, numericInput(ns("chr_ratio"), "Signal ratio cutoff", 0.33, 0.05, 1, step = 0.01)),
              column(3, checkboxInput(ns("chr_trim"), "Mark trimmed region (red)", TRUE)),
              column(3, sliderInput(ns("chr_start"), "Start at base", 1, 100, 1, step = 1))),
            plotOutput(ns("chr_plot"), height = "auto")),
          .card("QC of all reads", tableOutput(ns("qc")))),

        tabPanel("Needs review", value = "todo", br(),
          p("Click a row to jump to that position in the Alignment tab."),
          DT::DTOutput(ns("todo"))),

        tabPanel("SNPs and variants", value = "snp", br(),
          uiOutput(ns("snp_head")),
          fluidRow(column(5, div(style = "display:flex;gap:6px",
                    numericInput(ns("snp_pos"), "Jump to reference base number", NA, 1, step = 1, width = "240px"),
                    div(style = "margin-top:32px", actionButton(ns("snp_go"), "Go", icon = icon("magnifying-glass"))))),
                   column(7, br(), downloadButton(ns("dl_snp"), "SNPs (CSV)"))),
          DT::DTOutput(ns("snps")),
          h5("Indels"), DT::DTOutput(ns("indels")),
          h5("Read differences from the consensus"), DT::DTOutput(ns("diffs"))),

        tabPanel("Translation", value = "trl", br(),
          fluidRow(column(5, selectInput(ns("trl_code"), "Genetic code", GENETIC_CODES, selected = "1")),
                   column(3, numericInput(ns("trl_minaa"), "Minimum ORF length (aa)", 30, 1, 5000)),
                   column(4, br(), downloadButton(ns("dl_orf"), "ORF proteins (FASTA)"),
                          downloadButton(ns("dl_trl"), "Six frames (text)"))),
          uiOutput(ns("trl_info")),
          .card("ORF map (six frames)", plotOutput(ns("trl_plot"), height = 280)),
          .card("Open reading frames", DT::DTOutput(ns("trl_orfs"))),
          .card("Six-frame translation", uiOutput(ns("trl_html")))),

        tabPanel("Sequence analysis", value = "sa", br(),
          uiOutput(ns("sa_info")),
          .card("Composition",
            DT::DTOutput(ns("sa_comp")),
            downloadButton(ns("dl_sa_comp"), "CSV"),
            helpText("GC and AT % are calculated on A, C, G and T only (ambiguous bases are counted in Other). ",
                     "GC skew = (G - C) / (G + C).")),
          .card("Restriction sites",
            fluidRow(column(4, uiOutput(ns("sa_contig_ui"))),
                     column(4, selectInput(ns("sa_filter"), "Show",
                       c("Enzymes that cut" = "cut", "Single cutters" = "single", "Enzymes that do not cut" = "none")))),
            selectizeInput(ns("sa_enz"), "Enzymes (leave empty = all in the list)",
                           choices = .sa_enzymes$Enzyme, multiple = TRUE, width = "100%"),
            DT::DTOutput(ns("sa_rst")),
            downloadButton(ns("dl_sa_rst"), "CSV"),
            h5("Fragments of a linear digest with the chosen enzymes"),
            DT::DTOutput(ns("sa_frag")),
            helpText("Common enzymes with palindromic sites; the position is the base after which the top strand is cut. ",
                     "Positions refer to the consensus without gaps."))),

        tabPanel("Distances and dendrogram", value = "phy", br(),
          fluidRow(column(6, plotOutput(ns("dist_heat"), height = 380)),
                   column(6, plotOutput(ns("dendro"), height = 380))),
          helpText("Distance = proportion of differing bases (p-distance) over columns covered by both reads (at least 30). ",
                   "Dendrogram: UPGMA on those distances; needs at least 3 reads.")),

        tabPanel("Consensus and export", value = "exp", br(),
          .card("Crop the consensus (optional)",
            radioButtons(ns("crop_mode"), NULL, inline = TRUE,
              c("No cropping" = "none", "Follow the primers" = "primers", "Follow the reference region" = "ref")),
            conditionalPanel(sprintf("input['%s'] == 'primers'", ns("crop_mode")),
              fluidRow(column(5, textInput(ns("pr_f"), "Forward primer (5'->3')", "")),
                       column(5, textInput(ns("pr_r"), "Reverse primer (5'->3')", "")),
                       column(2, numericInput(ns("pr_mm"), "Max mismatches", 2, 0, 8))),
              checkboxInput(ns("pr_remove"), "Remove the primer sequences (keep only the insert between them)", TRUE),
              helpText("IUPAC codes are allowed in the primers. The reverse primer is searched as its reverse complement. ",
                       "If the primers match the opposite strand, the consensus is reverse-complemented automatically.")),
            conditionalPanel(sprintf("input['%s'] == 'ref'", ns("crop_mode")),
              fluidRow(column(4, numericInput(ns("cr_from"), "Reference position from", NA, 1)),
                       column(4, numericInput(ns("cr_to"), "to", NA, 1))),
              helpText("Reference base numbers (reference-guided mode). Leave a box empty to not cut on that side.")),
            conditionalPanel(sprintf("input['%s'] != 'none'", ns("crop_mode")),
              numericInput(ns("cr_cov"), "Also trim the ends covered by fewer than N reads", 1, 1, 100, width = "320px"),
              helpText("For a Forward + Reverse pair, N = 2 keeps only the overlap.")),
            uiOutput(ns("crop_msg"))),
          checkboxInput(ns("keep_gaps"), "Keep real gaps as '-' in the consensus (positions you deleted are always removed, never shown as gaps; translation and BLAST ignore the gaps)", TRUE),
          checkboxInput(ns("split_gaps"), "Split into several contigs at coverage gaps (otherwise gaps are filled with N)", TRUE),
          verbatimTextOutput(ns("cons")),
          downloadButton(ns("dl_fasta"), "Consensus FASTA"),
          downloadButton(ns("dl_aln"), "Alignment FASTA"),
          downloadButton(ns("dl_csv"), "Per-position decisions (CSV)"),
          downloadButton(ns("dl_var"), "Variants (CSV)"),
          downloadButton(ns("dl_report"), "HTML report"))
      )
    )
  )
}

# ------------------------------------------------------------------------------
# Server
# ------------------------------------------------------------------------------
assemblyServer <- function(id, email = "") {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    get_email <- function() if (is.function(email)) email() else email
    rv <- reactiveValues(aln = NULL, ov = character(), sel = NA_integer_,
                         mask = NULL, delcol = NULL, excl = character(), undo = list(), view_first = 1L,
                         files = data.frame(name = character(), path = character(), size = numeric(), key = character(),
                                            where = character(), sample = character(), stringsAsFactors = FALSE),
                         refs = NULL, samples = list(), cur = NULL)

    params <- reactive({
      P <- default_params()
      P$trim_method <- input$trim_method; P$q_trim <- input$q_trim; P$win <- input$win
      P$mott_cutoff <- input$mott_cutoff; P$q_good <- input$q_good; P$q_poor <- input$q_poor
      P$ratio <- input$ratio; P$type4 <- input$type4
      P
    })

    # --- input files: they ACCUMULATE (a fileInput alone would keep only the last selection) ----------------
    up_dir <- file.path(tempdir(), paste0("sangerly_up_", paste(sample(c(letters, 0:9), 8, TRUE), collapse = "")))
    dir.create(up_dir, showWarnings = FALSE, recursive = TRUE)
    up_n <- 0L
    add_files <- function(name, path, size, key, where) {
      if (!length(path)) return(invisible())
      f <- rv$files
      keep <- !(key %in% f$key) & !duplicated(key)
      nskip <- sum(!keep)
      if (any(keep)) {
        nd <- data.frame(name = name[keep], path = path[keep], size = size[keep], key = key[keep],
                         where = where[keep], stringsAsFactors = FALSE)
        nd$sample <- group_label(nd$name, .nz(input$group_mode, "name_norun"))
        rv$files <- rbind(f, nd)
      }
      showNotification(sprintf("%d file(s) added%s. %d in the list.", sum(keep),
                               if (nskip) sprintf(", %d already in the list (skipped)", nskip) else "",
                               nrow(rv$files)), type = "message", duration = 4)
    }
    observeEvent(input$reads, {
      u <- input$reads; req(nrow(u) > 0)
      dest <- character(nrow(u))
      for (i in seq_len(nrow(u))) {
        up_n <<- up_n + 1L
        dest[i] <- file.path(up_dir, sprintf("%05d_%s", up_n, u$name[i]))
        file.copy(u$datapath[i], dest[i], overwrite = TRUE)
      }
      add_files(u$name, dest, u$size, paste(u$name, u$size, sep = "|"), rep("(uploaded)", nrow(u)))
    })
    add_dir <- function(d) {
      d <- trimws(gsub("^\"|\"$", "", d))
      if (!nzchar(d) || !dir.exists(d)) { showNotification("That folder does not exist.", type = "warning"); return() }
      p <- list.files(d, pattern = "\\.ab1$", ignore.case = TRUE, recursive = isTRUE(input$dir_rec), full.names = TRUE)
      if (!length(p)) { showNotification("No .ab1 files in that folder.", type = "warning"); return() }
      p <- normalizePath(p, winslash = "/")
      add_files(basename(p), p, file.size(p), p, dirname(p))
    }
    observeEvent(input$dir_add, add_dir(.nz(input$dir_path, "")))
    observeEvent(input$browse_dir, {
      d <- tryCatch(utils::choose.dir(caption = "Choose the folder with the .ab1 files"), error = function(e) NA_character_)
      if (is.na(d) || !nzchar(d)) return()
      add_dir(d)
    })
    observeEvent(input$browse_add, {
      p <- tryCatch(utils::choose.files(filters = matrix(c("AB1 files (*.ab1)", "*.ab1"), 1), multi = TRUE),
                    error = function(e) character())
      if (!length(p)) return()
      p <- normalizePath(p, winslash = "/")
      add_files(basename(p), p, file.size(p), p, dirname(p))
    })
    observeEvent(input$files_clear, { rv$files <- rv$files[0, , drop = FALSE] })
    observeEvent(input$files_rm, {
      s <- input$files_table_rows_selected
      if (!length(s)) { showNotification("Select rows in the list first (click them).", type = "message"); return() }
      rv$files <- rv$files[-s, , drop = FALSE]
    })
    observeEvent(input$group_mode, {
      if (nrow(rv$files)) rv$files$sample <- group_label(rv$files$name, input$group_mode)
    }, ignoreInit = TRUE)
    output$files_count <- renderUI({
      n <- nrow(rv$files)
      if (!n) return(helpText("No files yet. Use Add files... (repeat for other folders)."))
      helpText(sprintf("%d file(s) in %d sample(s). See the Input files tab.", n, length(unique(rv$files$sample))))
    })
    output$files_table <- DT::renderDT({
      d <- rv$files
      show <- data.frame(File = d$name, Sample = d$sample,
                         Direction = vapply(d$name, function(x) .nz(ab1_direction(x), ""), "", USE.NAMES = FALSE),
                         Folder = d$where, stringsAsFactors = FALSE)
      DT::datatable(show, rownames = FALSE, selection = "multiple", editable = list(target = "cell", disable = list(columns = c(0, 2, 3))),
                    options = list(pageLength = 25, scrollX = TRUE, order = list()))
    })
    observeEvent(input$files_table_cell_edit, {
      i <- input$files_table_cell_edit
      r <- as.integer(i$row); v <- trimws(as.character(i$value))
      if (is.na(r) || r < 1 || r > nrow(rv$files)) return()
      if (!nzchar(v)) { showNotification("The sample name cannot be empty.", type = "warning"); rv$files$sample <- rv$files$sample; return() }
      rv$files$sample[r] <- v
    })

    # --- reference(s) -----------------------------------------------------------------------------------------
    observeEvent(input$ref_file, {
      rv$refs <- tryCatch(read_fasta_all(input$ref_file$datapath),
        error = function(e) { showNotification(conditionMessage(e), type = "error"); NULL })
    })
    observeEvent(input$ref_fetch, {
      req(nzchar(trimws(input$ref_acc)))
      rv$refs <- tryCatch(list(fetch_ncbi_fasta(trimws(input$ref_acc), get_email())),
        error = function(e) { showNotification(conditionMessage(e), type = "error"); NULL })
    })
    output$ref_info <- renderText({
      if (is.null(rv$refs)) "No reference loaded."
      else if (length(rv$refs) == 1) sprintf("%s (%d bp)", substr(rv$refs[[1]]$name, 1, 50), nchar(rv$refs[[1]]$seq))
      else sprintf("%d reference sequences (the best match is chosen for each sample).", length(rv$refs))
    })

    # --- one alignment per sample; switching stores / restores the editing state -----------------------------
    cur_state <- function() list(aln = rv$aln, ov = rv$ov, mask = rv$mask, delcol = rv$delcol, excl = rv$excl, undo = rv$undo)
    load_state <- function(st) {
      rv$aln <- st$aln; rv$ov <- .nz(st$ov, character()); rv$mask <- st$mask; rv$delcol <- st$delcol
      rv$excl <- .nz(st$excl, character()); rv$undo <- .nz(st$undo, list()); rv$sel <- NA_integer_
    }
    store_cur <- function() if (!is.null(rv$cur) && !is.null(rv$aln)) rv$samples[[rv$cur]] <- cur_state()
    snapshot <- function() {            # all samples, the current one live
      s <- rv$samples
      if (!is.null(rv$cur) && !is.null(rv$aln) && rv$cur %in% names(s)) s[[rv$cur]] <- cur_state()
      s
    }
    switch_sample <- function(nm) {
      if (is.null(nm) || !nzchar(nm) || identical(nm, rv$cur) || !(nm %in% names(rv$samples))) return(invisible())
      store_cur(); rv$cur <- nm; load_state(rv$samples[[nm]])
    }

    # --- assemble all samples ---------------------------------------------------------------------------------
    observeEvent(input$go, {
      fl <- rv$files
      if (!nrow(fl)) {
        showNotification("Add the .ab1 files first (Add files... button).", type = "warning")
        updateTabsetPanel(session, "tabs", selected = "files"); return()
      }
      if (input$mode == "ref" && is.null(rv$refs)) { showNotification("Load a reference first.", type = "warning"); return() }
      P <- params()
      res <- withProgress(message = "Assembling samples", value = 0, {
        build_batch(fl, input$mode, P, rv$refs, progress = function(i, n, nm)
          setProgress(value = (i - 1) / n, detail = sprintf("%d / %d: %s", i, n, nm)))
      })
      states <- lapply(res, function(aln) {
        if (!isTRUE(aln$ok))
          return(list(aln = aln, ov = character(), mask = NULL, delcol = NULL, excl = character(), undo = list()))
        m <- matrix(FALSE, nrow(aln$base), ncol(aln$base)); dc <- rep(FALSE, ncol(aln$base))
        if (isTRUE(input$auto_del)) {
          r <- aln_delete_low(aln, m, dc, input$auto_q, "orphans"); m <- r$mask; dc <- r$delcol
        }
        list(aln = aln, ov = character(), mask = m, delcol = dc, excl = character(), undo = list())
      })
      isok <- vapply(states, function(s) isTRUE(s$aln$ok), TRUE)
      rv$samples <- states
      first <- if (any(isok)) names(states)[which(isok)[1]] else names(states)[1]
      rv$cur <- first; load_state(states[[first]])
      if (any(!isok))
        showNotification(sprintf("%d of %d sample(s) could not be assembled (%s). Details in the Samples tab.",
                                 sum(!isok), length(isok), paste(names(states)[!isok], collapse = ", ")),
                         type = "warning", duration = NULL)
      updateTabsetPanel(session, "tabs", selected = if (length(states) > 1) "smp" else "aln")
    })

    ok <- reactive(!is.null(rv$aln) && isTRUE(rv$aln$ok) && !is.null(rv$mask) && !is.null(rv$delcol))
    emask <- reactive({
      req(ok()); m <- rv$mask
      if (length(rv$excl)) m[rv$aln$labels %in% rv$excl, ] <- TRUE
      aln_mask_all(m, rv$delcol)
    })
    cols <- reactive({ req(ok()); aln_decide(rv$aln, emask(), rv$delcol) })
    fin <- reactive({ req(ok()); aln_final(rv$aln, cols(), rv$ov, emask()) })
    todo <- reactive({ cf <- fin(); cf[cf$Review & !cf$Manual, , drop = FALSE] })
    snps <- reactive({ req(ok()); aln_snps(rv$aln, fin(), emask()) })
    is_ref <- reactive(ok() && !is.null(rv$aln$ref))

    # --- sync controls after an alignment is built ---------------------------------------
    observeEvent(rv$aln, {
      req(ok())
      labs <- rv$aln$labels
      updateSelectInput(session, "rng_read", choices = labs)
      updateNumericInput(session, "rng_to", value = ncol(rv$aln$base), max = ncol(rv$aln$base))
      updateNumericInput(session, "rng_from", max = ncol(rv$aln$base))
      grp <- split(labs, rv$aln$strand)
      names(grp) <- ifelse(names(grp) == "Forward", "Forward (as is)", "Reverse (reverse-complemented)")
      updateSelectInput(session, "qc_read", choices = grp)
    })
    output$jump_hint <- renderUI({
      if (!ok()) return(NULL)
      cf <- fin()
      if (is_ref()) helpText(sprintf("Reference base number (covered by reads: %d-%d).", min(cf$RefPos), max(cf$RefPos)))
      else helpText(sprintf("Column number (1-%d).", nrow(cf)))
    })
    output$col_hint <- renderUI({
      if (!ok()) return(NULL)
      helpText(if (is_ref()) "Reference base numbers." else "Column numbers.")
    })

    # --- navigation -------------------------------------------------------------------------
    goto <- function(k, switch_tab = FALSE) {
      rv$sel <- as.integer(k)
      session$sendCustomMessage("sv_goto", list(id = ns("viewer"), col = as.integer(k)))
      if (switch_tab) updateTabsetPanel(session, "tabs", selected = "aln")
    }
    do_jump <- function(v, switch_tab = FALSE) {
      if (is.null(v) || is.na(v)) return()
      cf <- fin()
      if (is_ref()) {
        k <- which(cf$RefPos == v)
        if (!length(k)) {
          showNotification(sprintf("Reference position %s is not covered by any read (covered: %d-%d).",
                                   format(v), min(cf$RefPos), max(cf$RefPos)), type = "warning"); return()
        }
      } else {
        k <- as.integer(v)
        if (k < 1 || k > nrow(cf)) { showNotification(sprintf("Column must be 1-%d.", nrow(cf)), type = "warning"); return() }
      }
      goto(k[1], switch_tab)
    }
    observeEvent(input$jump_go, { req(ok()); do_jump(input$jump_pos) })
    observeEvent(input$snp_go, { req(ok()); do_jump(input$snp_pos, TRUE) })

    # --- viewer (JavaScript canvas) -> server: selection, deletion, undo, visible window ----------
    observeEvent(input$viewer_sel, {
      req(ok()); k <- suppressWarnings(as.integer(input$viewer_sel$col))
      if (!is.na(k) && k >= 1 && k <= ncol(rv$aln$base)) rv$sel <- k
    })
    observeEvent(input$viewer_view, {
      f <- suppressWarnings(as.integer(input$viewer_view$first))
      if (!is.na(f)) rv$view_first <- f
    })
    push_undo <- function() {
      rv$undo <- utils::head(c(list(list(mask = rv$mask, delcol = rv$delcol, ov = rv$ov)), rv$undo), 50)
    }
    observeEvent(input$viewer_del, {
      req(ok()); x <- input$viewer_del; nc <- ncol(rv$aln$base)
      a <- max(1L, as.integer(x$from)); b <- min(nc, as.integer(x$to)); i <- as.integer(x$read)
      if (is.na(a) || is.na(b) || a > b) return()
      ks <- a:b
      if (is.na(i) || i < 1L) {                      # Consensus row: whole positions
        ks <- ks[!rv$delcol[ks]]
        if (!length(ks)) { showNotification("Those positions are already deleted.", type = "message"); return() }
        push_undo(); dc <- rv$delcol; dc[ks] <- TRUE; rv$delcol <- dc
        showNotification(sprintf("Deleted %d position(s) from all reads and the consensus. Ctrl+Z undoes.", length(ks)),
                         type = "message")
      } else {                                       # a read row: only that read's bases
        if (i > nrow(rv$aln$base)) return()
        ks <- ks[!is.na(rv$aln$base[i, ks]) & !rv$mask[i, ks]]
        if (!length(ks)) { showNotification("Nothing to delete there (no base, or already deleted).", type = "message"); return() }
        push_undo(); m <- rv$mask; m[i, ks] <- TRUE; rv$mask <- m
        showNotification(sprintf("Deleted %d base(s) of read %s. Ctrl+Z undoes.", length(ks), rv$aln$labels[i]),
                         type = "message")
      }
    })
    observeEvent(input$viewer_undo, {
      req(ok())
      if (!length(rv$undo)) { showNotification("Nothing to undo.", type = "message"); return() }
      s <- rv$undo[[1]]; rv$undo <- rv$undo[-1]
      rv$mask <- s$mask; rv$delcol <- s$delcol
      if (!is.null(s$ov)) rv$ov <- s$ov
    })
    # Keyboard editing of the consensus: A/C/G/T/N/IUPAC/"-" set the base, "=" the reference base, "0" resets
    observeEvent(input$viewer_set, {
      req(ok()); x <- input$viewer_set; nc <- ncol(rv$aln$base)
      a <- max(1L, as.integer(x$from)); b <- min(nc, as.integer(x$to)); tok <- as.character(x$base)
      if (is.na(a) || is.na(b) || a > b) return()
      ks <- (a:b)[!rv$delcol[a:b]]
      if (!length(ks)) { showNotification("Those positions are deleted; restore them first (Ctrl+Z).", type = "message"); return() }
      cf <- fin()
      if (tok == "REF" && !is_ref()) { showNotification("The reference base is only available in reference-guided mode.", type = "message"); return() }
      push_undo()
      if (tok == "RULE") {
        rv$ov <- rv$ov[!(names(rv$ov) %in% as.character(ks))]
      } else {
        bases <- if (tok == "REF") cf$RefBase[ks] else rep(tok, length(ks))
        okb <- !is.na(bases) & nzchar(bases)
        ov <- rv$ov; ov[as.character(ks[okb])] <- bases[okb]; rv$ov <- ov
      }
    })
    observeEvent(input$nxt, {
      req(ok()); t <- todo()$Col
      if (!length(t)) { showNotification("Nothing left to review.", type = "message"); return() }
      cur <- if (is.na(rv$sel)) rv$view_first - 1L else rv$sel
      nx <- t[t > cur]; goto(if (length(nx)) nx[1] else t[1])
    })
    observeEvent(input$prev, {
      req(ok()); t <- todo()$Col
      if (!length(t)) return()
      cur <- if (is.na(rv$sel)) rv$view_first + 1L else rv$sel
      pv <- t[t < cur]; goto(if (length(pv)) pv[length(pv)] else t[length(t)])
    })
    observeEvent(input$snp_nxt, {
      req(ok())
      if (!is_ref()) { showNotification("SNP navigation needs reference-guided mode.", type = "message"); return() }
      t <- snps()$Col
      if (!length(t)) { showNotification("No SNPs relative to the reference.", type = "message"); return() }
      cur <- if (is.na(rv$sel)) rv$view_first - 1L else rv$sel
      nx <- t[t > cur]; goto(if (length(nx)) nx[1] else t[1])
    })

    # --- overview ---------------------------------------------------------------------------------
    do_crop <- function(cf, isref) {
      notes <- character()
      mode <- if (is.null(input$crop_mode)) "none" else input$crop_mode
      keep <- rep(TRUE, nrow(cf))
      if (mode != "none") {
        if (mode == "ref") {
          if (!isref) notes <- c(notes, "Reference cropping needs reference-guided mode; ignored.")
          else {
            a <- input$cr_from; b <- input$cr_to
            if (!is.null(a) && !is.na(a)) keep <- keep & !is.na(cf$RefPos) & cf$RefPos >= a
            if (!is.null(b) && !is.na(b)) keep <- keep & !is.na(cf$RefPos) & cf$RefPos <= b
          }
        }
        mc <- input$cr_cov
        if (!is.null(mc) && !is.na(mc) && mc > 1) {
          w <- which(cf$Cov >= mc)
          idx <- seq_len(nrow(cf))
          keep <- if (length(w)) keep & idx >= min(w) & idx <= max(w) else rep(FALSE, nrow(cf))
        }
        if (!all(keep)) {
          kk <- which(keep)
          notes <- c(notes, if (length(kk)) sprintf("Cropped to alignment columns %d-%d%s (%d of %d columns).",
                       min(kk), max(kk), if (isref && !all(is.na(cf$RefPos[kk]))) sprintf(", reference %d-%d", min(cf$RefPos[kk]), max(cf$RefPos[kk])) else "",
                       length(kk), nrow(cf))
                     else "Nothing is left after cropping: relax the crop settings.")
        }
      }
      ct <- consensus_contigs(cf[keep, , drop = FALSE], isTRUE(input$split_gaps), keep_gaps = isTRUE(input$keep_gaps) && mode != "primers")
      if (mode == "primers") {
        if (!nzchar(.clean_primer(input$pr_f)) && !nzchar(.clean_primer(input$pr_r))) {
          notes <- c(notes, "Enter at least one primer.")
        } else {
          mm <- if (is.null(input$pr_mm) || is.na(input$pr_mm)) 2L else as.integer(input$pr_mm)
          for (nm in names(ct)) {
            r <- crop_by_primers(ct[[nm]], input$pr_f, input$pr_r, isTRUE(input$pr_remove), mm)
            ct[[nm]] <- r$seq
            notes <- c(notes, paste0(nm, ": ", r$notes))
          }
        }
      }
      if (isTRUE(input$keep_gaps) && mode == "primers")
        notes <- c(notes, "Primer cropping works on the sequence without gaps, so real gaps are not shown as '-' here.")
      list(contigs = ct, notes = notes)
    }
    crop_res <- reactive({ req(ok()); do_crop(fin(), is_ref()) })
    contigs <- reactive(crop_res()$contigs)
    output$crop_msg <- renderUI({
      n <- crop_res()$notes
      if (length(n)) div(class = "alert alert-info", tags$ul(lapply(n, tags$li)))
    })
    output$summary <- renderUI({
      if (is.null(rv$aln) || isTRUE(rv$aln$ok)) return(NULL)
      div(class = "alert alert-danger", rv$aln$msg)
    })
    output$welcome <- renderUI({
      if (ok()) {
        return(tagList(
          if (length(rv$aln$notes)) div(class = "alert alert-warn", tags$ul(lapply(rv$aln$notes, tags$li))),
          .card("Consensus", verbatimTextOutput(ns("cons_prev")))))
      }
      step <- function(n, txt) div(class = "step", div(class = "step-n", n), div(txt))
      .card("Start here",
        step(1, "Choose a mode and add the .ab1 files in the left panel (as many times as you like, from any folder; plus a reference for reference-guided mode). Files with the same sample name are assembled together, and all samples run in one go."),
        step(2, "Check the Input files tab (edit the Sample column if the grouping is wrong), adjust trimming if needed, then press ASSEMBLE."),
        step(3, "Review the yellow positions in the Alignment tab: delete low-Phred bases or whole positions, or decide manually."),
        step(4, "Check SNPs, translation and read distances; download the consensus or the HTML report."),
        step(5, "Continue to the BLAST tab for species identification."))
    })
    output$cons_prev <- renderText({
      req(ok()); ct <- contigs()
      paste(vapply(names(ct), function(n) sprintf(">%s (%d bp)\n%s", n, nchar(gsub("-", "", ct[[n]])), ct[[n]]), ""), collapse = "\n\n")
    })
    output$boxes <- renderUI({
      req(ok()); aln <- rv$aln; B <- aln_eff_base(aln, emask())
      qm <- mean(aln$q[!is.na(B) & B != "-"], na.rm = TRUE)
      fluidRow(
        column(4, .infobox("dna", "Reads", nrow(aln$base), "#2c7be5")),
        column(4, .infobox("ruler-horizontal", "Alignment columns", ncol(aln$base), "#6f42c1")),
        column(4, .infobox("check", "Consensus (bp)", paste(nchar(gsub("-", "", unlist(contigs()))), collapse = " + "), "#27ae60")),
        column(4, .infobox("flag", "Needs review", nrow(todo()), "#e67e22")),
        column(4, .infobox("scissors", "Deleted bases / positions",
                           sprintf("%d / %d", sum(emask() & !is.na(aln$base) & !matrix(rv$delcol, nrow(aln$base), ncol(aln$base), byrow = TRUE)),
                                   sum(rv$delcol)), "#c0392b")),
        column(4, .infobox("chart-bar", "Mean Phred (used bases)", sprintf("%.1f", qm), "#17a2b8")))
    })

    # --- samples: selector, summary table, batch exports -------------------------------------------------------
    st_cf <- function(st) {
      m <- st$mask
      if (length(st$excl)) m[st$aln$labels %in% st$excl, ] <- TRUE
      em <- aln_mask_all(m, st$delcol)
      cols <- aln_decide(st$aln, em, st$delcol)
      list(cf = aln_final(st$aln, cols, st$ov, em), em = em)
    }
    sample_contigs <- function(st) do_crop(st_cf(st)$cf, !is.null(st$aln$ref))$contigs
    output$sample_bar <- renderUI({
      nm <- names(rv$samples); if (length(nm) < 2) return(NULL)
      okv <- vapply(rv$samples, function(s) isTRUE(s$aln$ok), TRUE)
      ch <- setNames(nm, paste0(nm, ifelse(okv, "", " (failed)")))
      div(class = "sidebar-sec", style = "display:flex;gap:10px;align-items:flex-end;flex-wrap:wrap",
          div(style = "flex:1;min-width:220px;margin-bottom:-12px",
              selectInput(ns("sample_sel"), sprintf("Sample (%d in this run)", length(nm)), ch,
                          selected = isolate(rv$cur), width = "100%")),
          div(actionButton(ns("smp_prev"), "< Previous"), actionButton(ns("smp_next"), "Next >")))
    })
    observeEvent(input$sample_sel, switch_sample(input$sample_sel), ignoreInit = TRUE)
    step_sample <- function(d) {
      nm <- names(rv$samples); if (length(nm) < 2) return()
      k <- match(rv$cur, nm); if (is.na(k)) k <- 1L
      updateSelectInput(session, "sample_sel", selected = nm[((k - 1L + d) %% length(nm)) + 1L])
    }
    observeEvent(input$smp_prev, step_sample(-1L))
    observeEvent(input$smp_next, step_sample(1L))

    smp_df <- reactive({
      s <- snapshot(); if (!length(s)) return(NULL)
      do.call(rbind, lapply(names(s), function(nm) {
        st <- s[[nm]]
        if (!isTRUE(st$aln$ok))
          return(data.frame(Sample = nm, Status = "failed", Reads = NA_integer_, Columns = NA_integer_, Consensus_bp = "",
                            Needs_review = NA_integer_, Mean_Phred = NA_real_, Note = .nz(st$aln$msg, ""), stringsAsFactors = FALSE))
        x <- st_cf(st); ct <- do_crop(x$cf, !is.null(st$aln$ref))$contigs
        data.frame(Sample = nm, Status = "ok", Reads = nrow(st$aln$base), Columns = ncol(st$aln$base),
                   Consensus_bp = paste(nchar(gsub("-", "", unlist(ct))), collapse = " + "),
                   Needs_review = sum(x$cf$Review & !x$cf$Manual),
                   Mean_Phred = round(mean(unlist(lapply(st$aln$reads, function(r) r$qual))), 1),
                   Note = paste(utils::head(.nz(st$aln$notes, ""), 2), collapse = " "), stringsAsFactors = FALSE)
      }))
    })
    output$smp_head <- renderUI({
      d <- smp_df()
      if (is.null(d)) return(div(class = "alert alert-info", "No samples yet. Add files, then press ASSEMBLE."))
      div(class = if (any(d$Status != "ok")) "alert alert-warn" else "alert alert-info",
          sprintf("%d sample(s): %d assembled, %d failed.", nrow(d), sum(d$Status == "ok"), sum(d$Status != "ok")))
    })
    output$smp_table <- DT::renderDT({
      d <- smp_df(); req(d)
      DT::datatable(d, rownames = FALSE, selection = "single", options = list(pageLength = 25, scrollX = TRUE)) |>
        DT::formatStyle("Status", color = DT::styleEqual(c("ok", "failed"), c("#1e7e34", "#c0392b")), fontWeight = "bold")
    })
    observeEvent(input$smp_table_rows_selected, {
      s <- input$smp_table_rows_selected; req(length(s) == 1)
      d <- smp_df(); nm <- d$Sample[s]
      switch_sample(nm); updateSelectInput(session, "sample_sel", selected = nm)
      if (d$Status[s] == "ok") updateTabsetPanel(session, "tabs", selected = "aln")
      DT::selectRows(DT::dataTableProxy("smp_table"), NULL)
    })
    longest_per_sample <- function() {          # one sequence per sample, for BLAST
      s <- snapshot(); out <- character()
      for (nm in names(s)) {
        if (!isTRUE(s[[nm]]$aln$ok)) next
        ct <- sample_contigs(s[[nm]]); ct <- ct[nzchar(unlist(ct))]
        if (!length(ct)) next
        w <- which.max(nchar(gsub("-", "", unlist(ct))))
        out[if (length(ct) > 1) paste0(nm, "_", names(ct)[w]) else nm] <- gsub("-", "", ct[[w]])
      }
      out
    }
    output$dl_all_fasta <- downloadHandler(
      filename = function() paste0("sangerly_consensus_all_", format(Sys.time(), "%Y%m%d_%H%M"), ".fasta"),
      content = function(file) {
        s <- isolate(snapshot()); out <- list()
        for (nm in names(s)) {
          if (!isTRUE(s[[nm]]$aln$ok)) next
          ct <- sample_contigs(s[[nm]]); ct <- ct[nzchar(unlist(ct))]
          for (j in seq_along(ct)) out[[if (length(ct) == 1) nm else paste0(nm, "_", names(ct)[j])]] <- ct[[j]]
        }
        writeLines(if (length(out)) fasta_of(out) else "", file)
      })
    output$dl_all_csv <- downloadHandler(
      filename = function() paste0("sangerly_samples_", format(Sys.time(), "%Y%m%d_%H%M"), ".csv"),
      content = function(file) utils::write.csv(isolate(smp_df()), file, row.names = FALSE))

    # --- alignment viewer: send the data; the browser draws, pans and selects ------------------------
    vid <- ns("viewer")
    observeEvent(rv$aln, {
      if (is.null(rv$aln) || !isTRUE(rv$aln$ok)) session$sendCustomMessage("sv_clear", list(id = vid))
      else session$sendCustomMessage("sv_static", sv_static(rv$aln, vid))
    }, ignoreNULL = FALSE)
    observe({
      req(ok())
      session$sendCustomMessage("sv_dyn", sv_dyn(fin(), emask(), rv$delcol, vid))
    })

    # --- deleting bases and positions ----------------------------------------------------------------------
    output$del_count <- renderUI({
      req(ok()); aln <- rv$aln
      dcm <- matrix(rv$delcol, nrow(aln$base), ncol(aln$base), byrow = TRUE)
      nb <- sum(rv$mask & !is.na(aln$base) & !dcm)
      p(strong(sprintf("%d base(s) and %d whole position(s) currently deleted", nb, sum(rv$delcol))),
        if (length(rv$excl)) sprintf(" (excluded reads: %s)", paste(rv$excl, collapse = ", ")))
    })
    observeEvent(input$del_apply, {
      req(ok())
      r <- aln_delete_low(rv$aln, rv$mask, rv$delcol, input$del_q, input$del_scope)
      rv$mask <- r$mask; rv$delcol <- r$delcol
      showNotification(sprintf("Deleted %d base(s) and %d consensus position(s) (Phred < %g).",
                               r$n_bases, r$n_cols, input$del_q), type = "message")
    })
    observeEvent(input$del_restore_all, {
      req(ok()); m <- rv$mask; m[] <- FALSE; rv$mask <- m
      rv$delcol <- rep(FALSE, ncol(rv$aln$base)); rv$excl <- character()
      showNotification("All bases, positions and reads restored.", type = "message")
    })
    observeEvent(input$col_del_go, {
      req(ok()); a <- input$col_from; b <- input$col_to
      if (is.na(a) || is.na(b) || a > b) { showNotification("Invalid range.", type = "warning"); return() }
      cf <- fin()
      ks <- if (is_ref()) which(cf$RefPos >= a & cf$RefPos <= b) else which(cf$Col >= a & cf$Col <= b)
      if (!length(ks)) { showNotification("No alignment position in that range.", type = "warning"); return() }
      dc <- rv$delcol; dc[ks] <- TRUE; rv$delcol <- dc
      showNotification(sprintf("%d position(s) deleted from reads and consensus.", length(ks)), type = "message")
    })
    observeEvent(input$rng_go, {
      req(ok()); i <- match(input$rng_read, rv$aln$labels)
      a <- input$rng_from; b <- input$rng_to
      if (is.na(i) || is.na(a) || is.na(b) || a > b) { showNotification("Invalid range.", type = "warning"); return() }
      b <- min(b, ncol(rv$aln$base)); m <- rv$mask; m[i, a:b] <- TRUE; rv$mask <- m
    })

    # --- details of the selected position --------------------------------------------------------------------------
    output$detail <- renderUI({
      req(ok()); k <- rv$sel
      if (is.na(k)) return(.card("Selected position", p(tags$em("Click a column in the alignment, use the navigation buttons, or jump to a base number to see the evidence and decide."))))
      r <- fin()[k, ]
      head_txt <- sprintf("Column %d%s | type %s | decision: %s (%s)", k,
                          if (!is.na(r$RefPos)) sprintf(" | reference position %d", r$RefPos) else "", r$Type, r$Final, r$Source)
      if (r$Type == "ColDel")
        return(.card(head_txt, p("This position was deleted: it is removed from every read and from the consensus."),
                     actionButton(ns("restore_pos"), "Restore this position")))
      B <- aln_eff_base(rv$aln, emask())
      obs <- unique(B[, k]); obs <- obs[!is.na(obs)]
      nuc <- obs[obs %in% c("A", "C", "G", "T")]
      ref_b <- if (!is.null(r$RefBase) && !is.na(r$RefBase) && nzchar(r$RefBase)) r$RefBase else NULL
      ch <- unique(c(r$Auto, ref_b, nuc, "A", "C", "G", "T", "N", "-"))
      if (length(nuc) >= 2) ch <- unique(c(ch, .iupac(nuc[1:2])))
      cov_reads <- rv$aln$labels[!is.na(rv$aln$base[, k])]
      .card(head_txt,
        p(strong("Rule recommendation: "), r$Recommendation, " -> ", code(r$Auto)),
        p(tags$em(r$Reason)),
        tableOutput(ns("evid")),
        fluidRow(
          column(3, selectInput(ns("pick"), "Final base", ch, selected = r$Final)),
          column(9, br(), actionButton(ns("apply_pick"), "Apply this base"),
                 actionButton(ns("accept_rule"), "Accept rule recommendation"),
                 if (!is.null(ref_b)) actionButton(ns("use_ref"), sprintf("Use reference base (%s)", ref_b),
                                                   class = "btn-outline-primary"),
                 actionButton(ns("reset_pos"), "Reset"))),
        fluidRow(
          column(4, if (length(cov_reads)) selectInput(ns("mask_read"), "Delete / restore the base of read", cov_reads)),
          column(8, br(),
                 if (length(cov_reads)) actionButton(ns("mask_cell"), "Delete this read's base here", icon = icon("scissors")),
                 if (length(cov_reads)) actionButton(ns("unmask_cell"), "Restore"),
                 actionButton(ns("del_pos"), "Delete whole position (all reads + consensus)", icon = icon("scissors"),
                              class = "btn-outline-danger"))))
    })
    output$evid <- renderTable({ req(ok(), !is.na(rv$sel)); aln_evidence_table(rv$aln, rv$sel, emask()) })

    set_final <- function(k, base) rv$ov[as.character(k)] <- base
    observeEvent(input$apply_pick, { req(!is.na(rv$sel)); set_final(rv$sel, input$pick) })
    observeEvent(input$use_ref, {
      req(ok(), !is.na(rv$sel)); b <- fin()$RefBase[rv$sel]
      if (!is.na(b) && nzchar(b)) set_final(rv$sel, b)
    })
    observeEvent(input$accept_rule, { req(!is.na(rv$sel)); set_final(rv$sel, fin()$Auto[rv$sel]) })
    observeEvent(input$reset_pos, {
      req(!is.na(rv$sel)); rv$ov <- rv$ov[names(rv$ov) != as.character(rv$sel)]
    })
    observeEvent(input$mask_cell, {
      req(ok(), !is.na(rv$sel)); i <- match(input$mask_read, rv$aln$labels); req(!is.na(i))
      m <- rv$mask; m[i, rv$sel] <- TRUE; rv$mask <- m
    })
    observeEvent(input$unmask_cell, {
      req(ok(), !is.na(rv$sel)); i <- match(input$mask_read, rv$aln$labels); req(!is.na(i))
      m <- rv$mask; m[i, rv$sel] <- FALSE; rv$mask <- m
    })
    observeEvent(input$del_pos, {
      req(ok(), !is.na(rv$sel)); dc <- rv$delcol; dc[rv$sel] <- TRUE; rv$delcol <- dc
    })
    observeEvent(input$restore_pos, {
      req(ok(), !is.na(rv$sel)); dc <- rv$delcol; dc[rv$sel] <- FALSE; rv$delcol <- dc
    })

    # --- tab: needs review ------------------------------------------------------------------------------------------------
    output$todo <- DT::renderDT({
      req(ok()); d <- todo()
      show <- data.frame(Col = d$Col, RefPos = d$RefPos, Cov = d$Cov, Type = d$Type, Rule = d$Auto,
                         Recommendation = d$Recommendation, stringsAsFactors = FALSE)
      DT::datatable(show, rownames = FALSE, selection = "single", options = list(pageLength = 15, scrollX = TRUE))
    })
    observeEvent(input$todo_rows_selected, {
      s <- input$todo_rows_selected; req(length(s) == 1)
      goto(todo()$Col[s], TRUE)
      DT::selectRows(DT::dataTableProxy("todo"), NULL)
    })

    # --- tab: SNPs and variants ---------------------------------------------------------------------------------------------
    output$snp_head <- renderUI({
      if (!ok()) return(p("Assemble first."))
      if (!is_ref()) return(div(class = "alert alert-warn",
        "SNPs relative to a reference need reference-guided mode. In de novo mode, see the read-difference table below."))
      s <- snps()
      div(class = "alert alert-info", sprintf("%d position(s) differ from the reference (%d single-base SNPs, %d ambiguous/IUPAC). Click a row to jump to it in the Alignment tab.",
          nrow(s), sum(s$Alt %in% c("A", "C", "G", "T")), sum(!s$Alt %in% c("A", "C", "G", "T"))))
    })
    output$snps <- DT::renderDT({
      req(ok()); DT::datatable(snps(), rownames = FALSE, selection = "single", options = list(pageLength = 10, scrollX = TRUE))
    })
    observeEvent(input$snps_rows_selected, {
      s <- input$snps_rows_selected; req(length(s) == 1)
      goto(snps()$Col[s], TRUE)
      DT::selectRows(DT::dataTableProxy("snps"), NULL)
    })
    output$indels <- DT::renderDT({
      req(ok()); DT::datatable(aln_indels(rv$aln, fin()), rownames = FALSE, options = list(pageLength = 8, scrollX = TRUE))
    })
    output$diffs <- DT::renderDT({
      req(ok()); DT::datatable(aln_read_diffs(rv$aln, fin(), emask()), rownames = FALSE, options = list(pageLength = 10, scrollX = TRUE))
    })

    # --- tab: read and QC ------------------------------------------------------------------------------------------------------
    qc_i <- reactive({ req(ok(), input$qc_read); i <- match(input$qc_read, rv$aln$labels); req(!is.na(i)); i })
    observeEvent(input$qc_read, {
      req(ok()); i <- match(input$qc_read, rv$aln$labels); req(!is.na(i))
      n <- length(rv$aln$reads[[i]]$raw$seq)
      updateSliderInput(session, "chr_start", max = max(1L, n), value = max(1L, rv$aln$reads[[i]]$qc$from))
    })
    output$qc_toggle <- renderUI({
      req(ok(), input$qc_read)
      out <- input$qc_read %in% rv$excl
      actionButton(ns("qc_toggle_btn"), if (out) "Include this read again" else "Exclude this read from the alignment",
                   class = if (out) "btn-success" else "btn-outline-danger")
    })
    observeEvent(input$qc_toggle_btn, {
      req(ok(), input$qc_read)
      rv$excl <- if (input$qc_read %in% rv$excl) setdiff(rv$excl, input$qc_read) else c(rv$excl, input$qc_read)
    })
    output$qc_msg <- renderUI({
      req(ok(), input$qc_read)
      if (input$qc_read %in% rv$excl) div(class = "alert alert-warn", "This read is excluded (all its bases are deleted). Results are recomputed without it.")
    })
    output$qc_boxes <- renderUI({
      req(ok()); rd <- rv$aln$reads[[qc_i()]]; q <- rd$raw$qual; u <- rd$qual
      fluidRow(
        column(3, .infobox("dna", "Raw length", length(q), "#6c757d")),
        column(3, .infobox("scissors", "After trimming", sprintf("%d (bases %d-%d)", length(u), rd$qc$from, rd$qc$to), "#f4a261")),
        column(3, .infobox("chart-bar", "Mean Phred raw -> used", sprintf("%.1f -> %.1f", mean(q), mean(u)), "#17a2b8")),
        column(3, .infobox("check", "Q20 / Q30 (used)", sprintf("%.0f%% / %.0f%%", 100 * mean(u >= 20), 100 * mean(u >= 30)), "#27ae60")))
    })
    output$phred_plot <- renderPlot({ req(ok()); plot_phred(rv$aln$reads[[qc_i()]]) }, res = 96)
    chr_h <- reactive(max(200, 175 * max(1, input$chr_nrows) + 50))
    output$chr_plot <- renderPlot({
      req(ok())
      p <- plot_read_chromatogram(rv$aln$reads[[qc_i()]], input$chr_start, input$chr_nper, input$chr_nrows,
                                  input$chr_ratio, isTRUE(input$chr_trim))
      req(p); p
    }, height = function() chr_h(), res = 96)
    qc_df <- reactive({
      req(ok()); aln <- rv$aln
      q <- do.call(rbind, lapply(seq_along(aln$reads), function(i) {
        rd <- aln$reads[[i]]
        data.frame(Read = aln$labels[i], Strand = aln$strand[i], Raw_length = rd$qc$n_full, Used = length(rd$qual),
                   Mean_Phred = round(mean(rd$qual), 1), Q20_pct = round(100 * mean(rd$qual >= 20), 1),
                   Q30_pct = round(100 * mean(rd$qual >= 30), 1),
                   Status = if (aln$labels[i] %in% rv$excl) "excluded" else "used", stringsAsFactors = FALSE)
      }))
      if (!is.null(aln$map_info)) q <- merge(q, aln$map_info[, c("Read", "Identity_to_ref")], by = "Read", sort = FALSE)
      q
    })
    output$qc <- renderTable(qc_df())

    # --- tab: translation (six frames, like ExPASy Translate) -------------------------------------------------------
    trl_src <- reactive({
      req(ok()); ct <- contigs(); w <- which.max(nchar(unlist(ct)))
      list(name = names(ct)[w], seq = gsub("-", "", ct[[w]]))
    })
    six <- reactive({
      s <- trl_src(); req(nzchar(s$seq), !is.na(input$trl_minaa))
      translate_six(s$seq, input$trl_code, input$trl_minaa)
    })
    output$trl_info <- renderUI({
      s <- trl_src(); t <- six()
      div(class = "alert alert-info",
          sprintf("Translating %s (%d nt) in all six frames, genetic code %s. %d ORF(s) of at least %d aa (Met to stop) are highlighted. The sequence is the (cropped) consensus from the Consensus and export tab.",
                  s$name, t$n, input$trl_code, nrow(t$orfs), t$min_aa))
    })
    output$trl_plot <- renderPlot({ plot_frames(six()) }, res = 96)
    output$trl_orfs <- DT::renderDT({
      DT::datatable(six()$orfs, rownames = FALSE, options = list(pageLength = 8, scrollX = TRUE))
    })
    output$trl_html <- renderUI({ six_frame_html(six(), 60L) })
    output$dl_trl <- downloadHandler(
      filename = function() paste0(base_name(), "_six_frames.txt"),
      content = function(file) writeLines(six_frame_text(six()), file))
    output$dl_orf <- downloadHandler(
      filename = function() paste0(base_name(), "_ORFs.faa"),
      content = function(file) {
        o <- six()$orfs
        if (!nrow(o)) { writeLines("", file); return() }
        wrap <- function(x) paste(substring(x, seq(1, nchar(x), by = 60), pmin(nchar(x), seq(60, nchar(x) + 59, by = 60))), collapse = "\n")
        hdr <- sprintf(">%s_%s_%s_frame%s_nt%d-%d_%daa", base_name(), o$ORF, ifelse(startsWith(o$Frame, "5"), "fwd", "rev"),
                       sub(".*Frame ", "", o$Frame), o$Start_nt, o$End_nt, o$Length_aa)
        writeLines(paste0(hdr, "\n", vapply(o$Protein, wrap, "", USE.NAMES = FALSE)), file)
      })

    # --- tab: sequence analysis (composition and restriction sites of the consensus) ---------------------------------------
    sa_df <- reactive({
      req(ok()); ct <- contigs()
      ids <- if (is.null(names(ct))) paste0("contig", seq_along(ct)) else names(ct)
      seqs <- vapply(ct, function(x) toupper(gsub("[^A-Za-z]", "", x)), "")
      keep <- nzchar(seqs); req(any(keep))
      data.frame(id = ids[keep], seq = unname(seqs[keep]), stringsAsFactors = FALSE)
    })
    output$sa_info <- renderUI({
      d <- sa_df()
      div(class = "alert alert-info",
          sprintf("Sequence: the (cropped) consensus of %s, %d contig(s). Gaps are ignored.",
                  if (!is.null(rv$cur)) rv$cur else "this sample", nrow(d)))
    })
    sa_comp <- reactive({
      d <- sa_df()
      do.call(rbind, lapply(seq_len(nrow(d)), function(i) .sa_composition(d$id[i], d$seq[i])))
    })
    output$sa_comp <- DT::renderDT(
      DT::datatable(sa_comp(), rownames = FALSE, options = list(pageLength = 8, scrollX = TRUE)))
    output$dl_sa_comp <- downloadHandler(
      filename = function() paste0(base_name(), "_composition.csv"),
      content = function(file) utils::write.csv(sa_comp(), file, row.names = FALSE))

    output$sa_contig_ui <- renderUI({
      d <- sa_df(); prev <- isolate(input$sa_contig)
      selectInput(ns("sa_contig"), "Contig", d$id,
                  selected = if (!is.null(prev) && prev %in% d$id) prev else d$id[which.max(nchar(d$seq))])
    })
    sa_enz <- reactive({
      e <- input$sa_enz
      if (is.null(e) || !length(e)) .sa_enzymes else .sa_enzymes[.sa_enzymes$Enzyme %in% e, , drop = FALSE]
    })
    sa_sites <- reactive({
      d <- sa_df(); k <- match(input$sa_contig, d$id); req(!is.na(k))
      list(n = nchar(d$seq[k]), sites = .sa_find_sites(d$seq[k], sa_enz()))
    })
    sa_rst_df <- reactive({
      r <- sa_sites(); enz <- sa_enz()
      d <- do.call(rbind, lapply(seq_len(nrow(enz)), function(i) {
        p <- if (is.null(r$sites)) integer(0) else sort(r$sites$CutAfter[r$sites$Enzyme == enz$Enzyme[i]])
        data.frame(Enzyme = enz$Enzyme[i], Site = enz$Site[i], Cuts = length(p),
                   Cut_after_position = paste(p, collapse = ", "), stringsAsFactors = FALSE)
      }))
      flt <- input$sa_filter
      if (is.null(flt)) flt <- "cut"
      keep <- switch(flt, cut = d$Cuts > 0, single = d$Cuts == 1, none = d$Cuts == 0)
      d <- d[keep, , drop = FALSE]
      d[order(d$Cuts, d$Enzyme), , drop = FALSE]
    })
    output$sa_rst <- DT::renderDT({
      d <- sa_rst_df()
      validate(need(nrow(d) > 0, "Nothing to show with this filter."))
      DT::datatable(d, rownames = FALSE, options = list(pageLength = 10, scrollX = TRUE))
    })
    output$dl_sa_rst <- downloadHandler(
      filename = function() paste0(base_name(), "_restriction_sites.csv"),
      content = function(file) utils::write.csv(sa_rst_df(), file, row.names = FALSE))
    output$sa_frag <- DT::renderDT({
      validate(need(length(input$sa_enz) > 0, "Choose one or more enzymes above to see the fragments."))
      r <- sa_sites()
      cuts <- if (is.null(r$sites)) integer(0) else
        sort(unique(r$sites$CutAfter[r$sites$CutAfter >= 1 & r$sites$CutAfter < r$n]))
      validate(need(length(cuts) > 0, "The chosen enzyme(s) do not cut this sequence."))
      st <- c(1L, cuts + 1L); en <- c(cuts, r$n)
      DT::datatable(data.frame(Fragment = seq_along(st), Start = st, End = en, Length = en - st + 1L),
                    rownames = FALSE, options = list(pageLength = 10, scrollX = TRUE))
    })

    # --- tab: distances and dendrogram --------------------------------------------------------------------------------------------------
    dmat <- reactive({ req(ok()); aln_dist(rv$aln, emask()) })
    output$dist_heat <- renderPlot({ req(ok()); plot_dist_heat(dmat()) }, res = 96)
    output$dendro <- renderPlot({
      req(ok()); hc <- dist_to_hclust(dmat())
      if (is.null(hc)) { plot.new(); text(0.5, 0.5, "A dendrogram needs at least 3 reads.") }
      else plot(hc, main = "Read dendrogram (UPGMA, p-distance)", xlab = "", sub = "", hang = -1)
    }, res = 96)

    # --- tab: consensus and export ----------------------------------------------------------------------------------------------------------
    output$cons <- renderText({
      if (!ok()) return("(no consensus yet)")
      ct <- contigs()
      paste(vapply(names(ct), function(n) sprintf(">%s (%d bp)\n%s", n, nchar(gsub("-", "", ct[[n]])), ct[[n]]), ""), collapse = "\n\n")
    })
    base_name <- reactive({
      gsub("[^A-Za-z0-9._-]+", "_", if (!is.null(rv$cur)) rv$cur else "sangerly")
    })
    output$dl_fasta <- downloadHandler(
      filename = function() paste0(base_name(), "_consensus.fasta"),
      content = function(file) { ct <- contigs(); names(ct) <- paste0(base_name(), "_", names(ct)); writeLines(fasta_of(ct), file) })
    output$dl_aln <- downloadHandler(
      filename = function() paste0(base_name(), "_alignment.fasta"),
      content = function(file) writeLines(aln_alignment_fasta(rv$aln, fin(), emask()), file))
    output$dl_csv <- downloadHandler(
      filename = function() paste0(base_name(), "_decisions.csv"),
      content = function(file) utils::write.csv(fin(), file, row.names = FALSE))
    output$dl_var <- downloadHandler(
      filename = function() paste0(base_name(), "_variants.csv"),
      content = function(file) utils::write.csv(aln_variants(rv$aln, fin()), file, row.names = FALSE))
    output$dl_snp <- downloadHandler(
      filename = function() paste0(base_name(), "_snps.csv"),
      content = function(file) utils::write.csv(snps(), file, row.names = FALSE))
    output$dl_report <- downloadHandler(
      filename = function() paste0(base_name(), "_report.html"),
      content = function(file) {
        cf <- fin(); t <- tryCatch(six(), error = function(e) NULL)
        stops <- if (is.null(t)) NULL else {
          o <- t$orfs
          list(desc = sprintf("Genetic code %s: %d ORF(s) of at least %d aa in six frames%s.", input$trl_code, nrow(o), t$min_aa,
                              if (nrow(o)) sprintf(" (longest: %s, %d aa)", o$Frame[1], o$Length_aa[1]) else ""),
               aa = if (nrow(o)) o$Protein[1] else "(no ORF)",
               table = o[, setdiff(names(o), "Protein"), drop = FALSE])
        }
        aln <- rv$aln
        dcm <- matrix(rv$delcol, nrow(aln$base), ncol(aln$base), byrow = TRUE)
        html <- build_report_html(aln, cf, contigs(), qc_df(), snps(), aln_indels(aln, cf),
                                  aln_read_diffs(aln, cf, emask()), stops,
                                  n_deleted_bases = sum(rv$mask & !is.na(aln$base) & !dcm),
                                  n_deleted_cols = sum(rv$delcol), excluded = rv$excl)
        writeLines(html, file, useBytes = TRUE)
      })

    # --- save / load session ------------------------------------------------------------------------------------------------------------------------
    output$dl_sess <- downloadHandler(
      filename = function() paste0("sangerly_session_", format(Sys.time(), "%Y%m%d_%H%M"), ".rds"),
      content = function(file) {
        saveRDS(list(version = 5L, samples = isolate(snapshot()), cur = isolate(rv$cur), refs = isolate(rv$refs)), file)
      })
    observeEvent(input$sess_file, {
      x <- tryCatch(readRDS(input$sess_file$datapath), error = function(e) NULL)
      if (is.list(x) && identical(x$version, 5L) && is.list(x$samples) && length(x$samples)) {
        states <- x$samples; cur <- if (!is.null(x$cur) && x$cur %in% names(states)) x$cur else names(states)[1]
      } else if (is.list(x) && identical(x$version, 4L) && !is.null(x$aln$base) && !is.null(x$mask) && !is.null(x$delcol)) {
        states <- list(sample = list(aln = x$aln, ov = x$ov, mask = x$mask, delcol = x$delcol, excl = x$excl, undo = list()))
        cur <- "sample"; x$refs <- if (!is.null(x$ref)) list(x$ref) else NULL
      } else {
        showNotification("This is not a valid Sangerly session file (v4 or v5).", type = "error"); return()
      }
      rv$samples <- states; rv$cur <- cur; rv$refs <- x$refs
      load_state(states[[cur]])
      updateTabsetPanel(session, "tabs", selected = if (length(states) > 1) "smp" else "aln")
      showNotification(sprintf("Session loaded (%d sample(s)).", length(states)), type = "message")
    })

    list(consensus = reactive({
           x <- longest_per_sample(); if (!length(x)) return(NULL)
           if (length(x) == 1) paste0(">", names(x), "\n", x) else x
         }),
         aln = reactive(rv$aln), decisions = fin)
  })
}
