# SangeRly - Sanger assembly (F+R, or many reads / tiled regions against a reference)
# -> alignment viewer -> delete low-Phred bases -> SNPs -> consensus / report
# -> BLAST (with FASTA downloads)
#
# Run from the sangerly/ folder:   shiny::runApp()
#
# Packages:
#   install.packages(c("shiny", "bslib", "ggplot2", "DT", "httr", "xml2", "jsonlite", "readxl", "writexl", "processx", "BiocManager"))
#   BiocManager::install(c("sangerseqR", "Biostrings", "pwalign", "BiocGenerics"))
#   MSA tab (at least one):  BiocManager::install(c("DECIPHER", "msa"))   [optional: the trimal program]
#   MSA tab, "Reference-guided (like ViralMSA)":  [optional: the minimap2 program. It has no Windows build: on Windows install it inside WSL
#                            (wsl --install, then  sudo apt install minimap2 ); without it a built-in R aligner is used]
#   Phylogeny tab:           install.packages(c("ape", "processx")) install.packages("pegas")  [optional: phangorn]   + the IQ-TREE program
#   Population structure tab: install.packages("pegas")  [optional; haploNet network]   install.packages(c("maps", "mapdata"))  [optional; base maps]


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
source("R/dialog_core.R")
source("R/pileup_core.R")
source("R/batch_core.R")
source("R/extras_core.R")
source("R/msa_ref_core.R")
source("R/ui_helpers.R")             # Save menu, folding cards, small hints
source("R/viewer_module.R")
source("R/assembly_module.R")
source("R/open_module.R")
source("R/blast_module.R")
source("R/msa_module.R")
source("R/tree_module.R")
source("R/visual_module.R")
source("R/barcode_module.R")
source("R/primer_core.R")
source("R/primer_module.R")
source("R/popgen_core.R")
source("R/popgen_module.R")
source("R/genbank_module.R")
source("R/project_module.R")
source("R/help_module.R")
source("R/tutorial_module.R")

.tlog("app.R: all files sourced")

theme <- bslib::bs_theme(version = 5, bootswatch = "flatly", primary = "#2c7be5",
                         base_font = bslib::font_collection("system-ui", "Segoe UI", "Arial", "sans-serif"))

# ---- Home page + logos ---------------------------------------------------------------------------------
# Put the seven logo files in www/img/ with these names (png). If a file is missing, the card and the tab
# still work; only the picture is left out.
# www/img is served by Shiny itself as "img/..." (no addResourcePath: it made Shiny warn about conflicting resource prefixes)
logo_dir <- file.path("www", "img")
if (dir.exists(logo_dir)) {
  message("SangeRly logos: ", paste(list.files(logo_dir, pattern = "\\.png$"), collapse = ", "))
} else {
  message("SangeRly: folder www/img was NOT found in ", getwd(), " - the logos will not show. ",
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
       desc_id = "Siapkan berkas Newick untuk iTOL dan ChiPlot dan gambar figur publikasi."),
  list(id = "popgen",   n = 6, img = "img/logo_popup.png",    name = "Population structure", name_id = "Struktur populasi",
       desc = "Haplotypes, diversity, Tajima's D, Fu's Fs, networks, AMOVA / PhiST, and export for DnaSP, Arlequin and PopART.",
       desc_id = "Haplotipe, keragaman, Tajima's D, Fu's Fs, jaringan, AMOVA / PhiST, serta ekspor untuk DnaSP, Arlequin, dan PopART."),
  list(id = "barcode",  n = 7, img = "img/logo_barcode.png",  name = "Genetic distances", name_id = "Jarak genetik",
       desc = "Distance matrices, barcode gap, PCoA, statistical tests and species identification for groups you define (species, genus, outgroup).",
       desc_id = "Matriks jarak, barcode gap, PCoA, uji statistik, dan identifikasi spesies untuk grup yang Anda tentukan (spesies, genus, outgroup)."),
  list(id = "primer",   n = 8, img = "img/logo_primer.png",   name = "Primer mining and design", name_id = "Penambangan dan desain primer",
       desc = "Find sequences of your target in BOLD and NCBI, design primers for PCR, sequencing, exon junctions or degenerate primers, check their quality and BLAST them.",
       desc_id = "Cari sekuens target di BOLD dan NCBI, rancang primer untuk PCR, sekuensing, sambungan antar-ekson atau primer degenerate, periksa kualitasnya, lalu BLAST."))

# extra tools (menu "More")
extras <- list(
  list(id = "open", icon = "eye", name = "Open and view files", name_id = "Buka dan lihat berkas",
       desc = "Just look at files: AB1 chromatograms, FASTA, FASTQ, GenBank, Newick trees and tables. No assembly needed.",
       desc_id = "Cukup membuka berkas: kromatogram AB1, FASTA, FASTQ, GenBank, pohon Newick, dan tabel. Tanpa assembly."),
  list(id = "project", icon = "folder-open", name = "Project and report", name_id = "Proyek dan laporan",
       desc = "Save or open the whole project, marker presets, Methods paragraph and HTML report.",
       desc_id = "Simpan atau buka seluruh proyek, preset gen, paragraf Methods, dan laporan HTML."),
  list(id = "genbank", icon = "file-export", name = "GenBank export", name_id = "Ekspor GenBank",
       desc = "FASTA with source modifiers and a feature table for BankIt / table2asn.",
       desc_id = "FASTA dengan source modifier dan feature table untuk BankIt / table2asn."),
  list(id = "tutorial", icon = "graduation-cap", name = "Tutorial", name_id = "Tutorial",
       desc = "Step-by-step lessons for every tool, with what you need, what to click and what you get.",
       desc_id = "Pelajaran langkah demi langkah untuk setiap alat: apa yang diperlukan, apa yang diklik, dan hasilnya."),
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
  /* 8 steps: 4 + 4 on wide windows, 3 per row on medium windows, 2 per row on narrow ones */
  .sg-grid:not(.sg-extra){grid-template-columns:repeat(4,minmax(0,1fr));gap:14px;max-width:1500px}
  .sg-grid:not(.sg-extra) .sg-card{padding:12px 10px 14px}
  .sg-grid:not(.sg-extra) .sg-card img{height:130px;margin-bottom:6px}
  .sg-grid:not(.sg-extra) .sg-card h4{font-size:17px}
  .sg-grid:not(.sg-extra) .sg-card p{font-size:12.5px;min-height:70px}
  @media (max-width:1150px){
    .sg-grid:not(.sg-extra){grid-template-columns:repeat(3,minmax(0,1fr));gap:16px}
    .sg-grid:not(.sg-extra) .sg-card img{height:140px}}
  @media (max-width:700px){
    .sg-grid:not(.sg-extra){grid-template-columns:repeat(2,minmax(0,1fr))}
    .sg-grid.sg-extra{grid-template-columns:repeat(2,minmax(0,1fr))}}
  /* the top menu shows logos / icons only (the name is the tooltip) */
  .sg-sr{position:absolute;width:1px;height:1px;overflow:hidden;clip:rect(0,0,0,0);white-space:nowrap}
  .navbar .nav-link .sg-navname{margin-left:4px}
  .navbar .nav-link img{height:44px;max-width:96px;margin-right:0}
  .navbar .nav-link .fa,.navbar .nav-link .fas{margin-right:0;font-size:22px}
  .navbar .nav-link.dropdown-toggle{font-size:0}
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
  .sg-grid.sg-extra{grid-template-columns:repeat(3,minmax(0,1fr));max-width:1100px;margin-top:14px}
  .sg-card.sm{padding:12px 12px 14px}
  .sg-card.sm .sg-ico{font-size:30px;color:#14a3b8;margin:4px 0 6px}
  .sg-card.sm h5{margin:0 0 4px;color:#1f3c5f;font-weight:700}
  .sg-card.sm p{margin:0;color:#52606d;font-size:13px}
  /* New-here link under the Home title */
  .sg-newhere{display:inline-flex;align-items:center;gap:8px;margin-top:12px;padding:6px 16px;border-radius:999px;cursor:pointer;
              background:#eef7ff;border:1px solid #b9d9f7;color:#0b3d91;font-weight:600;font-size:14px;text-decoration:none}
  .sg-newhere:hover{background:#0b3d91;color:#fff;border-color:#0b3d91}
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
      h1(icon("dna"), " SangeRly"),
      p(i18n("Sanger sequencing: from raw reads to species identification and phylogenetic trees. Choose a step to start.",
             "Sekuensing Sanger: dari read mentah sampai identifikasi spesies dan pohon filogenetik. Pilih satu tahap untuk memulai.")),
      tags$a(class = "sg-newhere tt-go", `data-go` = "tutorial", tabindex = 0, role = "button",
             icon("graduation-cap"), i18n("New here? Take the tutorial", "Baru memakai? Ikuti tutorialnya"))),
  div(class = "sg-grid",
      lapply(steps, function(s)
        div(class = "sg-card", tabindex = 0, role = "button", `data-go` = s$id,
            div(class = "sg-step", s$n),
            tags$img(src = img_src(s$img), alt = s$name, onerror = "this.style.visibility='hidden'"),
            h4(i18n(s$name, s$name_id)),
            p(i18n(s$desc, s$desc_id)),
            tags$span(class = "btn btn-primary btn-sm", i18n("Open", "Buka"), " ", icon("arrow-right"))))),
  div(class = "sg-flow", i18n("Assembly  >  BLAST  >  MSA and trimming  >  Phylogeny  >  Visualisation  >  Population structure  >  Genetic distances  >  Primer mining and design   (you can open any step at any time; a green tick marks a step that has results)",
                              "Perakitan  >  BLAST  >  MSA dan trimming  >  Filogeni  >  Visualisasi  >  Struktur populasi  >  Jarak genetik  >  Penambangan dan desain primer   (setiap tahap bisa dibuka kapan saja; tanda centang hijau = tahap sudah punya hasil)")),
  div(class = "sg-extra-title", i18n("More tools", "Alat tambahan")),
  div(class = "sg-grid sg-extra",
      lapply(extras, function(x)
        div(class = "sg-card sm", tabindex = 0, role = "button", `data-go` = x$id,
            div(class = "sg-ico", icon(x$icon)),
            h5(i18n(x$name, x$name_id)),
            p(i18n(x$desc, x$desc_id))))))

# a tab title with its logo next to the name
nav_title <- function(s) tagList(tags$img(src = img_src(s$img), alt = "", onerror = "this.style.display='none'"), i18n(s$name, s$name_id))
# top menu: the logo only; the name shows as a tooltip, and in place of the logo if the picture is missing
nav_logo <- function(s) tags$span(title = s$name,
  tags$img(src = img_src(s$img), alt = s$name,
           onerror = "this.style.display='none'; var n = this.nextElementSibling; if (n) n.style.display = 'inline';"),
  tags$span(class = "sg-navname", style = "display:none", s$name))
# top menu: an icon only, with a tooltip and a hidden text for screen readers
nav_icon <- function(ic, name) tags$span(title = name, icon(ic), tags$span(class = "sg-sr", name))
step <- function(id) steps[[which(vapply(steps, function(s) s$id, "") == id)]]
extra <- function(id) extras[[which(vapply(extras, function(x) x$id, "") == id)]]
extra_title <- function(id) { x <- extra(id); tagList(icon(x$icon), i18n(x$name, x$name_id)) }

ui <- navbarPage(
  id = "nav",
  title = div(id = "sg_brand", icon("dna"), " SangeRly"),
  windowTitle = "SangeRly",
  theme = theme,
  header = tagList(sangerly_css(), svViewerHead(), home_css,
                   tags$link(rel = "stylesheet", href = img_src("sangerly_embed.css")),   # look used inside the Qt window (active only with ?embed=1)
                   tags$script(src = img_src("sangerly_embed.js")),
                   .sg_ui_head(),                                                        # density, side tabs, drag-to-move (R/ui_helpers.R)
                   tags$head(tags$link(rel = "icon", type = "image/x-icon", href = "favicon.ico")),
                   # desktop version: if the browser drops the connection (sleeping window), start a new session
                   tags$script(HTML("$(document).on('shiny:disconnected', function(){ setTimeout(function(){ location.reload(); }, 1500); });")),
                   # clicking a card (or the SangeRly name) opens that tab
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
                         document.querySelectorAll('details.tt-lesson[data-go]').forEach(function(l){
                           l.classList.toggle('tt-hasres', ids.indexOf(l.getAttribute('data-go')) >= 0); });
                       });
                     })();")),
                   # Indonesian for the whole program: www/sangerly_i18n.js translates every text found in www/i18n_id.json
                   tags$script(src = "sangerly_i18n.js"),
                   # language switch EN | ID
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
  tabPanel(nav_icon("house", "Home"), value = "home", home_ui),
  tabPanel(nav_icon(extra("open")$icon, extra("open")$name), value = "open", openUI("open")),
  tabPanel(nav_logo(step("assembly")), value = "assembly", assemblyUI("asm")),
  tabPanel(nav_logo(step("blast")), value = "blast",
    div(class = "container-fluid",
      fluidRow(column(12, textInput("email", "Contact email for NCBI (recommended; not stored)", "",
                                    placeholder = "name@institute.org", width = "320px"))),
      blastUI("blast"))),
  tabPanel(nav_logo(step("msa")), value = "msa", msaUI("msa")),
  tabPanel(nav_logo(step("tree")), value = "tree", treeUI("tree")),
  tabPanel(nav_logo(step("vis")), value = "vis", visualUI("vis")),
  tabPanel(nav_logo(step("popgen")), value = "popgen", popgenUI("popgen")),
  tabPanel(nav_logo(step("barcode")), value = "barcode", barcodeUI("barcode")),
  tabPanel(nav_logo(step("primer")), value = "primer", primerUI("primer")),
  navbarMenu("More", icon = icon("ellipsis"),
    tabPanel(extra_title("project"), value = "project", projectUI("project")),
    tabPanel(extra_title("genbank"), value = "genbank", genbankUI("genbank")),
    tabPanel(extra_title("tutorial"), value = "tutorial", tutorialUI("tutorial")),
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
  open <- openServer("open", msa_add = msa$add, go = function(tab) updateNavbarPage(session, "nav", selected = tab))
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
  popgen <- popgenServer("popgen", aln_trim = msa$trimmed, aln_all = msa$aln)
  barcode <- barcodeServer("barcode", aln_trim = msa$trimmed, aln_all = msa$aln,
    # aligned sequence files opened in "Open and view files"
    aln_extra = open$aln,
    # groups of the Population structure tab: its populations, its groups (sequence sets) and the ancestry clusters
    group_sources = reactive({
      out <- list()
      pp <- tryCatch(popgen$pop(), error = function(e) NULL)
      if (is.data.frame(pp) && nrow(pp) > 0) {
        out[["Population structure: populations"]] <- data.frame(Sequence = pp$Sequence, Group = pp$Population, stringsAsFactors = FALSE)
        if (length(unique(pp$Group)) > 1)
          out[["Population structure: groups (sets)"]] <- data.frame(Sequence = pp$Sequence, Group = pp$Group, stringsAsFactors = FALSE)
      }
      sg <- tryCatch(popgen$struct_groups(), error = function(e) NULL)
      if (length(sg)) out[["Population structure: ancestry clusters"]] <- sg
      out
    }),
    # clades of the NJ and IQ-TREE trees of the Phylogeny tab
    tree_sources = reactive({
      nj <- tryCatch(trees$nj(), error = function(e) NULL); iq <- tryCatch(trees$iq(), error = function(e) NULL)
      out <- list()
      if (!is.null(nj$tree)) out[["Phylogeny: NJ tree"]] <- nj$tree
      if (!is.null(iq$tree)) out[["Phylogeny: IQ-TREE"]] <- iq$tree
      out
    }))
  # step 8: primer mining and design (sequences of BOLD / NCBI / the MSA tab, the assembly consensus, NCBI BLAST)
  primer <- primerServer("primer", consensus = asm$consensus, aln_trim = msa$trimmed, aln_all = msa$aln, email = email)
  observeEvent(msa$to_primer(), {            # MSA tab: 'Send to Primer design'
    rq <- msa$to_primer(); req(!is.null(rq))
    primer$use(rq$kind, rq$src)
    updateNavbarPage(session, "nav", selected = "primer")
  }, ignoreInit = TRUE)
  proj <- projectServer("project", asm = asm, blast = blast, msa = msa, trees = trees, barcode = barcode, parent_session = session,
                        popgen = popgen)   # popgen: the Methods paragraph describes the Structure analysis
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
             if (!is.null(barcode$result())) "barcode",
             if (!is.null(primer$designed())) "primer")
    session$sendCustomMessage("sg_done", list(ids = as.list(ids)))
  })

  # the data of the session, for the "Navigation area" of the Qt window (folders and what they hold); read-only, it changes nothing
  observe({
    it <- function(name, info = "", go, sub = "", kind = "") list(name = as.character(name), info = as.character(info), go = go, sub = sub, kind = kind)
    grp <- function(key, title, items, go, sub = "") list(key = key, title = title, go = go, sub = sub, items = unname(items))
    take <- function(x, n = 300L) utils::head(x, n)
    cons <- tryCatch(asm$consensus(), error = function(e) NULL)
    cn <- if (length(cons)) {
      n <- names(cons)
      if (is.null(n)) n <- if (length(cons) == 1L && startsWith(cons[[1]], ">")) sub("^>([^\n]*).*", "\\1", cons[[1]]) else paste("Consensus", seq_along(cons))
      n } else character()
    res <- tryCatch(blast$results(), error = function(e) NULL)
    pool <- tryCatch(msa$pool(), error = function(e) NULL)
    pid <- if (is.data.frame(pool) && nrow(pool) > 0) as.character(pool$id) else character()
    aln <- tryCatch(msa$aln(), error = function(e) NULL)
    trm <- tryCatch(msa$trimmed(), error = function(e) NULL)
    nj <- tryCatch(trees$nj(), error = function(e) NULL); iq <- tryCatch(trees$iq(), error = function(e) NULL)
    bc <- tryCatch(barcode$result(), error = function(e) NULL)
    pr <- tryCatch(primer$designed(), error = function(e) NULL)
    dims <- function(a) sprintf("%d sequences x %d columns", length(a), nchar(a[[1]]))
    nav <- list(
      grp("assembly", "Assembly (consensus)", lapply(take(cn), function(n) it(n, "", "assembly", kind = "sample")), "assembly"),
      grp("blast", "BLAST", if (is.data.frame(res) && nrow(res) > 0) list(it("BLAST hits", sprintf("%d hits", nrow(res)), "blast")) else list(), "blast"),
      grp("seqs", sprintf("Sequences (%d)", length(pid)), lapply(take(pid), function(n) it(n, "", "msa", "Sequences")), "msa", "Sequences"),
      grp("aln", "Alignments",
          c(if (length(aln)) list(it("Alignment", dims(aln), "msa", "Editor")),
            if (length(trm)) list(it("Trimmed alignment", dims(trm), "msa", "Trimmed result"))), "msa", "Editor"),
      grp("trees", "Trees",
          c(if (!is.null(nj$tree)) list(it("NJ tree", sprintf("%d tips", ape::Ntip(nj$tree)), "tree", "NJ tree")),
            if (!is.null(iq$tree)) list(it("IQ-TREE tree", sprintf("%d tips", ape::Ntip(iq$tree)), "tree", "IQ-TREE tree"))), "tree", "NJ tree"),
      grp("dist", "Genetic distances", if (!is.null(bc)) list(it("Barcode gap result", "", "barcode", "Distribution")) else list(), "barcode"),
      grp("primer", "Primers", if (!is.null(pr)) list(it("Designed primers", "", "primer")) else list(), "primer"))
    session$sendCustomMessage("sg_nav", list(groups = nav))
  })
}

shinyApp(ui, server)
