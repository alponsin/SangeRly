# Sangerly - Sanger assembly (F+R, or many reads / tiled regions against a reference)
# -> alignment viewer -> delete low-Phred bases -> SNPs -> consensus / report
# -> BLAST (with FASTA downloads)
#
# Run from the sangerly/ folder:   shiny::runApp()
#
# Packages:
#   install.packages(c("shiny", "bslib", "ggplot2", "DT", "httr", "xml2", "jsonlite", "BiocManager"))
#   BiocManager::install(c("sangerseqR", "Biostrings", "pwalign", "BiocGenerics"))
#   MSA tab (at least one):  BiocManager::install(c("DECIPHER", "msa"))   [optional: the trimal program]
#   Phylogeny tab:           install.packages(c("ape", "processx"))   [optional: phangorn]   + the IQ-TREE program


library(shiny)
# timing lines in the desktop log (%LOCALAPPDATA%/Sangerly/session.log); does nothing outside the desktop launcher
.tlog <- function(m) if (exists("sglog", mode = "function", envir = globalenv())) get("sglog", envir = globalenv())(m)
.tlog("app.R: shiny loaded")
options(shiny.maxRequestSize = Inf)   # no limit on the size of uploaded files (FASTA, .ab1, ...)
library(ggplot2)
if (!file.exists(file.path("www", "sangerly_viewer.js")))
  stop("www/sangerly_viewer.js was not found in ", getwd(), ". ",
       "Put the www folder next to app.R (it holds the alignment viewer), then run the app again.")
source("R/assembly_core.R")
source("R/seqtools_core.R")
source("R/stats_core.R")
source("R/pileup_core.R")
source("R/batch_core.R")
source("R/extras_core.R")
source("R/viewer_module.R")
source("R/assembly_module.R")
source("R/blast_module.R")
source("R/msa_module.R")
source("R/tree_module.R")
source("R/visual_module.R")
source("R/barcode_module.R")
source("R/genbank_module.R")
source("R/project_module.R")
source("R/help_module.R")
.tlog("app.R: all files sourced")

theme <- bslib::bs_theme(version = 5, bootswatch = "flatly", primary = "#2c7be5",
                         base_font = bslib::font_collection("system-ui", "Segoe UI", "Arial", "sans-serif"))

# ---- Home page + logos ---------------------------------------------------------------------------------
# Put the five logo files in www/img/ with these names (png). If a file is missing, the card and the tab
# still work; only the picture is left out.
# Serve the logo folder explicitly (works wherever the app is started from; checked when the app starts)
logo_dir <- file.path("www", "img")
if (dir.exists(logo_dir)) {
  addResourcePath("img", normalizePath(logo_dir))
  message("Sangerly logos: ", paste(list.files(logo_dir, pattern = "\\.png$"), collapse = ", "))
} else {
  message("Sangerly: folder www/img was NOT found in ", getwd(), " - the logos will not show. ",
          "Put the logo_*.png files in ", file.path(getwd(), "www", "img"), ".")
}

# small helper: a text in two languages (the EN | ID switch in the menu bar changes it without reloading)
i18n <- function(en, idn) tags$span(class = "i18n", `data-en` = en, `data-id` = idn, en)

steps <- list(
  list(id = "assembly", n = 1, img = "img/logo_assembly.png", name = "Assembly", name_id = "Perakitan",
       desc = "Assemble Forward + Reverse .ab1 reads, check the chromatograms, edit the consensus.",
       desc_id = "Rakit read Forward + Reverse (.ab1), periksa kromatogram, dan edit konsensus."),
  list(id = "blast",    n = 2, img = "img/logo_blast.png",    name = "BLAST", name_id = "BLAST",
       desc = "Identify the species: search the consensus against NCBI databases.",
       desc_id = "Identifikasi spesies: cari konsensus pada basis data NCBI."),
  list(id = "msa",      n = 3, img = "img/logo_msa.png",      name = "MSA and trimming", name_id = "MSA dan trimming",
       desc = "Align many sequences and trim the poorly aligned columns.",
       desc_id = "Sejajarkan banyak sekuens dan pangkas kolom yang sejajarnya buruk."),
  list(id = "tree",     n = 4, img = "img/logo_phylo.png",    name = "Phylogeny", name_id = "Filogeni",
       desc = "Build NJ and maximum-likelihood (IQ-TREE) trees with ModelFinder.",
       desc_id = "Bangun pohon NJ dan maximum-likelihood (IQ-TREE) dengan ModelFinder."),
  list(id = "vis",      n = 5, img = "img/logo_visual.png",   name = "Visualisation", name_id = "Visualisasi",
       desc = "Prepare Newick files for iTOL and ChiPlot and draw publication figures.",
       desc_id = "Siapkan berkas Newick untuk iTOL dan ChiPlot dan gambar figur publikasi."))

# extra tools (menu "More")
extras <- list(
  list(id = "project", icon = "folder-open", name = "Project and report", name_id = "Proyek dan laporan",
       desc = "Save or open the whole project, marker presets, Methods paragraph and HTML report.",
       desc_id = "Simpan atau buka seluruh proyek, preset gen, paragraf Methods, dan laporan HTML."),
  list(id = "barcode", icon = "chart-simple", name = "Barcode gap", name_id = "Barcode gap",
       desc = "Distances within and between species (K2P) and the gap between them.",
       desc_id = "Jarak dalam dan antar spesies (K2P) serta celah di antaranya."),
  list(id = "genbank", icon = "file-export", name = "GenBank export", name_id = "Ekspor GenBank",
       desc = "FASTA with source modifiers and a feature table for BankIt / table2asn.",
       desc_id = "FASTA dengan source modifier dan feature table untuk BankIt / table2asn."),
  list(id = "help", icon = "circle-question", name = "Help", name_id = "Bantuan",
       desc = "Installation, use, and what to do about antivirus, slow start and other problems.",
       desc_id = "Instalasi, cara pakai, serta antivirus, buka lambat, dan masalah lain."),
  list(id = "about", icon = "circle-info", name = "About and citation", name_id = "Tentang dan sitasi",
       desc = "Versions of the programs, how to cite, and the references of the methods.",
       desc_id = "Versi program, cara mengutip, dan pustaka metode."))

home_css <- tags$style(HTML("
  .sg-hero{text-align:center;padding:26px 10px 8px}
  .sg-hero h1{font-weight:700;color:#1f3c5f;margin:0}
  .sg-hero p{color:#6b7785;margin:6px 0 0;font-size:16px}
  .sg-grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(250px,1fr));gap:18px;
           max-width:1250px;margin:22px auto 30px;padding:0 14px}
  .sg-card{position:relative;background:#fff;border-radius:14px;box-shadow:0 1px 4px rgba(0,0,0,.13);
           padding:16px 16px 18px;cursor:pointer;text-align:center;border:2px solid transparent;
           transition:transform .15s ease, box-shadow .15s ease, border-color .15s ease}
  .sg-card:hover{transform:translateY(-4px);box-shadow:0 8px 20px rgba(44,123,229,.22);border-color:#2c7be5}
  .sg-card:focus{outline:none;border-color:#2c7be5}
  .sg-step{position:absolute;top:10px;left:10px;width:28px;height:28px;border-radius:50%;background:#2c7be5;
           color:#fff;font-weight:700;display:flex;align-items:center;justify-content:center;font-size:14px}
  .sg-card img{width:100%;height:170px;object-fit:contain;margin-bottom:10px}
  .sg-card h4{margin:0 0 6px;color:#1f3c5f;font-weight:700}
  .sg-card p{margin:0 0 12px;color:#52606d;font-size:13.5px;min-height:54px}
  .sg-flow{text-align:center;color:#6b7785;font-size:13px;margin-bottom:24px}
  .navbar .nav-link img{height:26px;width:auto;margin-right:6px;vertical-align:middle}
  .navbar .navbar-brand{cursor:pointer}
  /* light navbar in the colours of the logos (navy, blue, teal, green) */
  .navbar{background:#fff !important;border-bottom:4px solid transparent;
          border-image:linear-gradient(90deg,#0b3d91,#1e88e5,#14a3b8,#22b07d) 1;
          box-shadow:0 2px 8px rgba(11,61,145,.12)}
  .navbar .navbar-brand{color:#0b2a5b !important;font-weight:700;margin-right:26px}
  .navbar .navbar-brand .fa,.navbar .navbar-brand .fas{color:#14a3b8}
  .navbar .navbar-nav .nav-link{color:#1f3c5f !important;font-weight:600;display:flex;align-items:center;
          padding:8px 14px;border-bottom:3px solid transparent}
  .navbar .navbar-nav .nav-link:hover{color:#1e88e5 !important}
  .navbar .navbar-nav .nav-link.active{color:#0b7f6d !important;border-bottom-color:#22b07d}
  .navbar .nav-link img{height:36px;max-width:70px;object-fit:contain;margin-right:8px}
  .navbar .nav-link .fa,.navbar .nav-link .fas{margin-right:6px}
  .navbar-toggler{border-color:#1f3c5f}
  /* Home page */
  .sg-grid{grid-template-columns:repeat(auto-fit,minmax(230px,1fr));max-width:1500px}
  .sg-step{background:linear-gradient(135deg,#0b3d91,#14a3b8)}
  .sg-card .btn-primary{background:linear-gradient(90deg,#0b3d91,#14a3b8);border:0}
  .sg-card:hover{border-color:#14a3b8;box-shadow:0 8px 20px rgba(20,163,184,.25)}
  .sg-hero h1{color:#0b2a5b}
  /* finished steps get a green tick */
  .sg-card.done{border-color:#22b07d}
  .sg-card.done::after{content:'\\2713';position:absolute;top:10px;right:10px;width:26px;height:26px;border-radius:50%;
          background:#22b07d;color:#fff;display:flex;align-items:center;justify-content:center;font-weight:700;font-size:15px}
  .navbar a.sg-done::after{content:'\\2713';color:#22b07d;margin-left:6px;font-weight:700}
  /* extra tools */
  .sg-extra-title{text-align:center;color:#52606d;font-weight:600;margin:6px 0 -6px}
  .sg-grid.sg-extra{max-width:1100px;margin-top:14px}
  .sg-card.sm{padding:12px 12px 14px}
  .sg-card.sm .sg-ico{font-size:30px;color:#14a3b8;margin:4px 0 6px}
  .sg-card.sm h5{margin:0 0 4px;color:#1f3c5f;font-weight:700}
  .sg-card.sm p{margin:0;color:#52606d;font-size:13px}
  /* language switch */
  #sg_lang{margin-left:auto;display:flex;gap:4px;align-items:center}
  #sg_lang button{border:1px solid #b8c2cc;background:#fff;color:#1f3c5f;border-radius:6px;padding:2px 9px;font-size:12px;font-weight:600}
  #sg_lang button.on{background:#0b3d91;color:#fff;border-color:#0b3d91}
"))

# adds the file's modification time to the address, so the browser never shows an old (cached) picture
img_src <- function(p) {
  f <- file.path("www", p)
  paste0(p, "?v=", if (file.exists(f)) as.integer(file.mtime(f)) else 0L)
}

home_ui <- div(
  div(class = "sg-hero",
      h1(icon("dna"), " Sangerly"),
      p(i18n("Sanger sequencing: from raw reads to species identification and phylogenetic trees. Choose a step to start.",
             "Sekuensing Sanger: dari read mentah sampai identifikasi spesies dan pohon filogenetik. Pilih satu tahap untuk memulai."))),
  div(class = "sg-grid",
      lapply(steps, function(s)
        div(class = "sg-card", tabindex = 0, role = "button", `data-go` = s$id,
            div(class = "sg-step", s$n),
            tags$img(src = img_src(s$img), alt = s$name, onerror = "this.style.visibility='hidden'"),
            h4(i18n(s$name, s$name_id)),
            p(i18n(s$desc, s$desc_id)),
            tags$span(class = "btn btn-primary btn-sm", i18n("Open", "Buka"), " ", icon("arrow-right"))))),
  div(class = "sg-flow", i18n("Assembly  >  BLAST  >  MSA and trimming  >  Phylogeny  >  Visualisation   (you can open any step at any time; a green tick marks a step that has results)",
                              "Perakitan  >  BLAST  >  MSA dan trimming  >  Filogeni  >  Visualisasi   (setiap tahap bisa dibuka kapan saja; tanda centang hijau = tahap sudah punya hasil)")),
  div(class = "sg-extra-title", i18n("More tools", "Alat tambahan")),
  div(class = "sg-grid sg-extra",
      lapply(extras, function(x)
        div(class = "sg-card sm", tabindex = 0, role = "button", `data-go` = x$id,
            div(class = "sg-ico", icon(x$icon)),
            h5(i18n(x$name, x$name_id)),
            p(i18n(x$desc, x$desc_id))))))

# a tab title with its logo next to the name
nav_title <- function(s) tagList(tags$img(src = img_src(s$img), alt = "", onerror = "this.style.display='none'"), i18n(s$name, s$name_id))
step <- function(id) steps[[which(vapply(steps, function(s) s$id, "") == id)]]
extra <- function(id) extras[[which(vapply(extras, function(x) x$id, "") == id)]]
extra_title <- function(id) { x <- extra(id); tagList(icon(x$icon), i18n(x$name, x$name_id)) }

ui <- navbarPage(
  id = "nav",
  title = div(id = "sg_brand", icon("dna"), " Sangerly"),
  windowTitle = "Sangerly",
  theme = theme,
  header = tagList(sangerly_css(), svViewerHead(), home_css,
                   tags$head(tags$link(rel = "icon", type = "image/x-icon", href = "favicon.ico")),
                   # desktop version: if the browser drops the connection (sleeping window), start a new session
                   tags$script(HTML("$(document).on('shiny:disconnected', function(){ setTimeout(function(){ location.reload(); }, 1500); });")),
                   # clicking a card (or the Sangerly name) opens that tab
                   tags$script(HTML("
                     $(document).on('click keydown', '.sg-card', function(e){
                       if (e.type === 'keydown' && e.key !== 'Enter' && e.key !== ' ') return;
                       e.preventDefault();
                       Shiny.setInputValue('go_tab', $(this).data('go'), {priority: 'event'});
                     });
                     $(document).on('click', '#sg_brand', function(){
                       Shiny.setInputValue('go_tab', 'home', {priority: 'event'});
                     });")),
                   # green ticks on the steps that have results (the server sends the list)
                   tags$script(HTML("
                     (function wait(){
                       if (!(window.Shiny && Shiny.addCustomMessageHandler)) return setTimeout(wait, 30);
                       Shiny.addCustomMessageHandler('sg_done', function(m){
                         var ids = m.ids || []; if (!Array.isArray(ids)) ids = [ids];
                         document.querySelectorAll('.sg-card[data-go]').forEach(function(c){
                           c.classList.toggle('done', ids.indexOf(c.getAttribute('data-go')) >= 0); });
                         document.querySelectorAll('.navbar a[data-value]').forEach(function(a){
                           a.classList.toggle('sg-done', ids.indexOf(a.getAttribute('data-value')) >= 0); });
                       });
                     })();")),
                   # language switch EN | ID (Home, menu and About page; the other tabs stay in English)
                   tags$script(HTML("
                     (function(){
                       function setLang(l){
                         try { localStorage.setItem('sg_lang', l); } catch(e) {}
                         document.querySelectorAll('.i18n').forEach(function(el){
                           var t = el.getAttribute('data-' + l); if (t !== null) el.textContent = t; });
                         document.querySelectorAll('#sg_lang button').forEach(function(b){
                           b.classList.toggle('on', b.getAttribute('data-l') === l); });
                       }
                       $(function(){
                         var nav = document.querySelector('.navbar > .container-fluid') || document.querySelector('.navbar');
                         if (nav && !document.getElementById('sg_lang')) {
                           var d = document.createElement('div'); d.id = 'sg_lang';
                           d.innerHTML = '<button type=\"button\" data-l=\"en\">EN</button><button type=\"button\" data-l=\"id\">ID</button>';
                           nav.appendChild(d);
                         }
                         var l = 'en'; try { l = localStorage.getItem('sg_lang') || 'en'; } catch(e) {}
                         setLang(l);
                       });
                       $(document).on('click', '#sg_lang button', function(){ setLang(this.getAttribute('data-l')); });
                     })();"))),
  tabPanel(tagList(icon("house"), "Home"), value = "home", home_ui),
  tabPanel(nav_title(step("assembly")), value = "assembly", assemblyUI("asm")),
  tabPanel(nav_title(step("blast")), value = "blast",
    div(class = "container-fluid",
      fluidRow(column(12, textInput("email", "Contact email for NCBI (recommended; not stored)", "",
                                    placeholder = "name@institute.org", width = "320px"))),
      blastUI("blast"))),
  tabPanel(nav_title(step("msa")), value = "msa", msaUI("msa")),
  tabPanel(nav_title(step("tree")), value = "tree", treeUI("tree")),
  tabPanel(nav_title(step("vis")), value = "vis", visualUI("vis")),
  navbarMenu("More", icon = icon("ellipsis"),
    tabPanel(extra_title("project"), value = "project", projectUI("project")),
    tabPanel(extra_title("barcode"), value = "barcode", barcodeUI("barcode")),
    tabPanel(extra_title("genbank"), value = "genbank", genbankUI("genbank")),
    tabPanel(extra_title("help"), value = "help", helpUI("help")),
    tabPanel(extra_title("about"), value = "about", aboutUI("about")))
)

server <- function(input, output, session) {
  if (exists("sangerly_register_session")) sangerly_register_session(session)
  # open the tab chosen on the Home page
  observeEvent(input$go_tab, updateNavbarPage(session, "nav", selected = input$go_tab))
  email <- reactive(input$email)
  asm <- assemblyServer("asm", email = email)
  # the "Use consensus from the assembly" button in the BLAST tab uses the longest contig
  blast <- blastServer("blast", consensus = asm$consensus, email = email)
  msa <- msaServer("msa", consensus = asm$consensus)
  # "Send to MSA" button of the BLAST tab: the hits (and the query) join the sequence list of the MSA tab
  observeEvent(blast$to_msa(), {
    x <- blast$to_msa(); req(x)
    n <- msa$add(x$recs, "BLAST")
    showNotification(sprintf("%d sequence(s) sent to the MSA tab (%d hit(s)%s). Press Align there.", nrow(x$recs), x$nh,
                             if (x$nq > 0) paste0(" + ", x$nq, " query") else ""), type = "message", duration = 10)
    updateNavbarPage(session, "nav", selected = "msa")
  }, ignoreInit = TRUE)
  trees <- treeServer("tree", aln_all = msa$aln, aln_trim = msa$trimmed)
  visualServer("vis", nj = trees$nj, iq = trees$iq)
  barcode <- barcodeServer("barcode", aln_trim = msa$trimmed, aln_all = msa$aln)
  proj <- projectServer("project", asm = asm, blast = blast, msa = msa, trees = trees, barcode = barcode, parent_session = session)
  genbankServer("genbank", pool = msa$pool, marker = proj$marker)
  aboutServer("about")

  # green ticks on the Home page and in the menu: the steps that already have results
  observe({
    res <- tryCatch(blast$results(), error = function(e) NULL)
    cons <- tryCatch(asm$consensus(), error = function(e) NULL)
    ids <- c(if (!is.null(cons)) "assembly",
             if (!is.null(res) && nrow(res) > 0) "blast",
             if (!is.null(msa$aln())) "msa",
             if (!is.null(trees$nj()) || !is.null(trees$iq())) c("tree", "vis"),
             if (!is.null(barcode$result())) "barcode")
    session$sendCustomMessage("sg_done", list(ids = as.list(ids)))
  })
}

shinyApp(ui, server)
