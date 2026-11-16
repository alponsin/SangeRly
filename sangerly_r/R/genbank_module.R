# SangeRly - "GenBank export" tab
#   Prepares your sequences for submission to GenBank:
#     * FASTA with source modifiers in the header ([organism=...] [isolate=...] ...), ready for BankIt or table2asn
#     * feature table (.tbl) for table2asn / tbl2asn: gene + CDS / rRNA / misc_RNA of the marker
#     * source table (tab separated) for table2asn (-src-file)
#   Always read the files and check them in BankIt or with table2asn before submitting: this tab prepares, NCBI validates.

.gb_cols <- c("SeqID", "Organism", "Isolate", "Country", "Collection_date", "Specimen_voucher", "Host", "Note")

.gb_clean_id <- function(x) {
  x <- gsub("[^A-Za-z0-9_.-]", "_", x)
  substr(x, 1, 25)                                    # GenBank sequence names: at most 25 characters
}

# Organism guessed from the name: "MG272989_Homo_sapiens" -> "Homo sapiens"
.gb_organism <- function(nm) {
  x <- sub("^[A-Z]{1,6}_?[0-9]{3,}(\\.[0-9]+)?[_ ]+", "", nm)
  tok <- strsplit(x, "[_ ]+")
  vapply(tok, function(t) { t <- t[nzchar(t)]; if (length(t) >= 2) paste(t[1:2], collapse = " ") else "" }, "")
}

genbankUI <- function(id) {
  ns <- shiny::NS(id)
  shiny::div(class = "container-fluid",
    shiny::fluidRow(
      shiny::column(4,
        .card("1. Sequences",
          shiny::helpText("Takes the sequence list of the MSA tab (gaps removed). BLAST hits are other people's sequences and are left out unless you tick the box."),
          shiny::checkboxInput(ns("with_blast"), "Also include sequences that came from BLAST", FALSE),
          shiny::actionButton(ns("load"), "Load / refresh from the MSA list", icon = shiny::icon("arrows-rotate"), class = "btn-primary"),
          shiny::uiOutput(ns("load_info"))),
        .card("2. Defaults for the empty cells",
          shiny::textInput(ns("d_org"), "Organism (scientific name)", "", placeholder = "Homo sapiens"),
          shiny::textInput(ns("d_country"), "Country", "", placeholder = "Indonesia: West Sumatra"),
          shiny::textInput(ns("d_date"), "Collection date", "", placeholder = "2025-08-14 or Aug-2025"),
          shiny::textInput(ns("d_host"), "Host (if parasite / symbiont)", ""),
          shiny::actionButton(ns("fill"), "Fill the empty cells with these defaults")),
        .card("3. Feature (from the marker chosen in the Project tab)",
          shiny::textInput(ns("f_gene"), "Gene symbol", "", placeholder = "COI"),
          shiny::selectInput(ns("f_type"), "Feature", c("CDS (protein coding)" = "CDS", "rRNA" = "rRNA", "misc_RNA (ITS region)" = "misc_RNA", "none (sequence only)" = "none")),
          shiny::textInput(ns("f_product"), "Product / note", "", width = "100%"),
          shiny::selectInput(ns("f_loc"), "Genome", c("nuclear / not specified" = "", "mitochondrion" = "mitochondrion", "chloroplast" = "chloroplast")),
          shiny::selectInput(ns("f_code"), "Genetic code (CDS)", c("1 Standard" = "1", "2 Vertebrate mitochondrial" = "2", "4 Mold / protozoan mitochondrial" = "4",
                                                                  "5 Invertebrate mitochondrial" = "5", "9 Echinoderm / flatworm mitochondrial" = "9", "11 Bacterial / plastid" = "11"), "1"),
          shiny::numericInput(ns("f_frame"), "codon_start (1, 2 or 3)", 1, min = 1, max = 3, step = 1),
          shiny::checkboxInput(ns("f_partial"), "Partial at both ends (usual for an amplicon)", TRUE))),
      shiny::column(8,
        .card("Source information of each sequence (double-click a cell to edit)",
          DT::DTOutput(ns("tbl")),
          shiny::helpText("Country: 'Country: region'. Collection date: YYYY, YYYY-MM, YYYY-MM-DD or DD-Mon-YYYY. Isolate: your sample code. ",
                          "Specimen_voucher: museum / herbarium code if any.")),
        shiny::uiOutput(ns("checks")),
        shiny::div(style = "display:flex;gap:6px;flex-wrap:wrap",
          shiny::downloadButton(ns("dl_fasta"), "FASTA with source modifiers"),
          shiny::downloadButton(ns("dl_tbl"), "Feature table (.tbl)"),
          shiny::downloadButton(ns("dl_src"), "Source table (.src, tab separated)")),
        shiny::helpText("BankIt (web): upload the FASTA, then add the feature in the form (or use the same information). ",
                        "table2asn (command line): table2asn -t template.sbt -i file.fsa -f file.tbl -src-file file.src -V v"))))
}

genbankServer <- function(id, pool = shiny::reactive(NULL), marker = shiny::reactive("none")) {
  shiny::moduleServer(id, function(input, output, session) {
    rv <- shiny::reactiveValues(meta = NULL, seq = NULL, ver = 0L)

    # the marker chosen in the Project tab fills the feature settings
    shiny::observeEvent(marker(), {
      m <- GENE_PRESETS[[marker()]]
      if (is.null(m)) return()
      shiny::updateTextInput(session, "f_gene", value = m$gene)
      shiny::updateSelectInput(session, "f_type", selected = m$feature)
      shiny::updateTextInput(session, "f_product", value = m$product)
      shiny::updateSelectInput(session, "f_loc", selected = m$location)
      shiny::updateSelectInput(session, "f_code", selected = m$code)
    }, ignoreInit = FALSE)

    shiny::observeEvent(input$load, {
      p <- pool()
      if (is.null(p) || !nrow(p)) { shiny::showNotification("The MSA tab has no sequences yet.", type = "warning"); return() }
      if (!isTRUE(input$with_blast)) p <- p[!p$source %in% c("BLAST", "Demo"), , drop = FALSE]
      if (!nrow(p)) {
        shiny::showNotification("Only BLAST / demo sequences are in the MSA list. Add your own sequences there (Assembly or FASTA), or tick the box.", type = "warning")
        return()
      }
      seqs <- toupper(gsub("[^A-Za-z]", "", p$seq))
      ids <- make.unique(.gb_clean_id(p$id), sep = "_")
      rv$seq <- stats::setNames(seqs, ids)
      m <- data.frame(SeqID = ids, Organism = .gb_organism(p$id), Isolate = p$id, Country = "", Collection_date = "",
                      Specimen_voucher = "", Host = "", Note = "", stringsAsFactors = FALSE)
      rv$meta <- m; rv$ver <- rv$ver + 1L
    })
    output$load_info <- shiny::renderUI({
      if (is.null(rv$meta)) shiny::helpText("Nothing loaded yet.")
      else shiny::helpText(sprintf("%d sequence(s) loaded.", nrow(rv$meta)))
    })

    shiny::observeEvent(input$fill, {
      m <- rv$meta; shiny::req(m)
      fl <- function(col, val) if (nzchar(trimws(val))) m[[col]][!nzchar(trimws(m[[col]]))] <<- trimws(val)
      fl("Organism", input$d_org); fl("Country", input$d_country); fl("Collection_date", input$d_date); fl("Host", input$d_host)
      rv$meta <- m; rv$ver <- rv$ver + 1L
    })

    output$tbl <- DT::renderDT({
      rv$ver
      m <- shiny::isolate(rv$meta); shiny::req(m)
      DT::datatable(m, rownames = FALSE, selection = "none",
                    editable = list(target = "cell", disable = list(columns = 0)),
                    options = list(pageLength = 15, scrollX = TRUE, dom = "ftp"))
    })
    shiny::observeEvent(input$tbl_cell_edit, {
      rv$meta <- DT::editData(rv$meta, input$tbl_cell_edit, rownames = FALSE)
    })

    mod <- function(name, v) if (!is.na(v) && nzchar(trimws(v))) sprintf(" [%s=%s]", name, trimws(v)) else ""
    feat_type <- function() .nz(input$f_type, "none")
    header <- function(i) {
      m <- rv$meta[i, ]
      loc <- .nz(input$f_loc, "")
      code_mod <- if (identical(feat_type(), "CDS")) {
        if (identical(loc, "mitochondrion")) mod("mgcode", input$f_code) else mod("gcode", input$f_code)
      } else ""
      paste0(">", m$SeqID, mod("organism", m$Organism), mod("isolate", m$Isolate), mod("country", m$Country),
             mod("collection_date", m$Collection_date), mod("specimen_voucher", m$Specimen_voucher), mod("host", m$Host),
             mod("note", m$Note), mod("location", loc), code_mod, " [moltype=genomic DNA]")
    }
    wrap60 <- function(s) substring(s, seq(1, max(1, nchar(s)), by = 60), pmin(seq(60, max(60, nchar(s) + 59), by = 60), nchar(s)))

    output$dl_fasta <- shiny::downloadHandler(function() "genbank_submission.fsa", function(file) {
      shiny::req(rv$meta, rv$seq)
      out <- unlist(lapply(seq_len(nrow(rv$meta)), function(i) c(header(i), wrap60(rv$seq[[rv$meta$SeqID[i]]]))))
      writeLines(out, file)
    })

    feature_lines <- function(id, len) {
      ty <- feat_type()
      if (identical(ty, "none")) return(NULL)
      lo <- if (isTRUE(input$f_partial)) "<1" else "1"; hi <- if (isTRUE(input$f_partial)) paste0(">", len) else as.character(len)
      gene <- trimws(.nz(input$f_gene, "")); prod <- trimws(.nz(input$f_product, ""))
      out <- c(sprintf(">Feature %s", id))
      if (nzchar(gene)) out <- c(out, sprintf("%s\t%s\tgene", lo, hi), sprintf("\t\t\tgene\t%s", gene))
      out <- c(out, sprintf("%s\t%s\t%s", lo, hi, ty))
      if (identical(ty, "misc_RNA")) { if (nzchar(prod)) out <- c(out, sprintf("\t\t\tnote\t%s", prod)) }
      else if (nzchar(prod)) out <- c(out, sprintf("\t\t\tproduct\t%s", prod))
      if (identical(ty, "CDS")) out <- c(out, sprintf("\t\t\tcodon_start\t%d", as.integer(.nz(input$f_frame, 1))),
                                         sprintf("\t\t\ttransl_table\t%s", .nz(input$f_code, "1")))
      out
    }
    output$dl_tbl <- shiny::downloadHandler(function() "genbank_submission.tbl", function(file) {
      shiny::req(rv$meta, rv$seq)
      out <- unlist(lapply(rv$meta$SeqID, function(id) feature_lines(id, nchar(rv$seq[[id]]))))
      if (is.null(out)) out <- "# no feature selected"
      writeLines(out, file)
    })
    output$dl_src <- shiny::downloadHandler(function() "genbank_submission.src", function(file) {
      shiny::req(rv$meta)
      m <- rv$meta
      d <- data.frame(Sequence_ID = m$SeqID, organism = m$Organism, isolate = m$Isolate, country = m$Country,
                      collection_date = m$Collection_date, specimen_voucher = m$Specimen_voucher, host = m$Host, note = m$Note,
                      stringsAsFactors = FALSE)
      d <- d[, c(TRUE, vapply(d[-1], function(x) any(nzchar(trimws(x))), TRUE)), drop = FALSE]    # leave out empty columns
      utils::write.table(d, file, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
    })

    output$checks <- shiny::renderUI({
      m <- rv$meta                       # refreshes after every edit of the table
      if (is.null(m) || is.null(rv$seq)) return(NULL)
      msgs <- character()
      lens <- nchar(rv$seq[m$SeqID])
      if (any(lens < 200)) msgs <- c(msgs, sprintf("%d sequence(s) are shorter than 200 bp: GenBank does not accept them (%s).", sum(lens < 200),
                                                    paste(utils::head(m$SeqID[lens < 200], 4), collapse = ", ")))
      if (any(!nzchar(trimws(m$Organism)))) msgs <- c(msgs, sprintf("%d sequence(s) have no organism name.", sum(!nzchar(trimws(m$Organism)))))
      dt <- trimws(m$Collection_date); dt <- dt[nzchar(dt)]
      okd <- grepl("^[0-9]{4}(-[0-9]{2}(-[0-9]{2})?)?$|^[0-9]{2}-[A-Za-z]{3}-[0-9]{4}$|^[A-Za-z]{3}-[0-9]{4}$", dt)
      if (any(!okd)) msgs <- c(msgs, "Some collection dates are not in a format that GenBank accepts (YYYY, YYYY-MM, YYYY-MM-DD, Mon-YYYY or DD-Mon-YYYY).")
      if (any(grepl("N{10,}", rv$seq[m$SeqID]))) msgs <- c(msgs, "Some sequences contain long runs of N: trim them before submitting.")
      if (!length(msgs)) return(shiny::div(class = "alert alert-success", "No problems found in this quick check. NCBI still validates the files."))
      shiny::div(class = "alert alert-warning", shiny::tags$b("Check before submitting:"), shiny::tags$ul(lapply(msgs, shiny::tags$li)))
    })
    invisible(NULL)
  })
}
