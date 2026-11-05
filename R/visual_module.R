# Sangerly - "Visualisation" tab
# Links to external tree viewers (iTOL, ChiPlot) and Newick files prepared for them.
# The viewers are websites: Sangerly itself stays offline; only the file you upload there leaves the computer.

visualUI <- function(id) {
  ns <- shiny::NS(id)
  site <- function(title, url, text) {
    .card(title,
      shiny::p(text),
      shiny::tags$a(href = url, target = "_blank", rel = "noopener noreferrer", class = "btn btn-primary",
                    shiny::icon("arrow-up-right-from-square"), paste(" Open", url)))
  }
  shiny::div(class = "container-fluid",
    shiny::fluidRow(
      shiny::column(6, site("iTOL - Interactive Tree Of Life", "https://itol.embl.de/",
        "Display, annotate and export phylogenetic trees (Newick files).")),
      shiny::column(6, site("ChiPlot", "https://chiplot.online/",
        "Online tool for drawing and styling trees and other charts."))),
    shiny::fluidRow(shiny::column(12,
      .card("Your trees from the Phylogeny tab",
        shiny::uiOutput(ns("status")),
        shiny::checkboxInput(ns("last_only"),
          "Keep only the last support value on the branches (IQ-TREE 'SH-aLRT/UFBoot' becomes 'UFBoot'; plain numbers are easier for the web viewers)", TRUE),
        shiny::div(style = "display:flex;gap:8px;flex-wrap:wrap",
          shiny::downloadButton(ns("dl_nj"), "NJ tree (Newick)"),
          shiny::downloadButton(ns("dl_iq"), "IQ-TREE tree (Newick)")),
        shiny::tags$ol(style = "margin-top:14px",
          shiny::tags$li("Build a tree in the Phylogeny tab, then download its Newick file here."),
          shiny::tags$li("Open iTOL or ChiPlot with the buttons above."),
          shiny::tags$li("Upload the Newick file on that website and annotate or export the figure there.")),
        shiny::helpText("The trees are exactly as computed (not rerooted here); root them in the viewer, ",
                        "or choose the root in the Phylogeny tab first. Sequence names were cleaned there (spaces and symbols became '_'). ",
                        "Please read the terms of use of each website before uploading data.")))))
}

visualServer <- function(id, nj = shiny::reactive(NULL), iq = shiny::reactive(NULL)) {
  shiny::moduleServer(id, function(input, output, session) {
    get_tree <- function(r) if (is.null(r)) NULL else r$tree
    prep <- function(tr) {
      if (isTRUE(input$last_only) && !is.null(tr$node.label))
        tr$node.label <- sub(".*/", "", as.character(tr$node.label))
      tr
    }
    output$status <- shiny::renderUI({
      a <- !is.null(get_tree(nj())); b <- !is.null(get_tree(iq()))
      shiny::helpText(paste0("NJ tree: ", if (a) "ready" else "not built yet", "   |   IQ-TREE tree: ",
                             if (b) "ready" else "not built yet"))
    })
    output$dl_nj <- shiny::downloadHandler(function() "tree_nj_for_web.nwk",
      function(file) { tr <- get_tree(nj()); shiny::req(tr); ape::write.tree(prep(tr), file) })
    output$dl_iq <- shiny::downloadHandler(function() "tree_iqtree_for_web.nwk",
      function(file) { tr <- get_tree(iq()); shiny::req(tr); ape::write.tree(prep(tr), file) })
  })
}
