# SangeRly - "Population genetics" (step 6): Haplotype | Network | AMOVA | Structure | Result
# The calculations are in popgen_core.R; the network viewer is www/sangerly_network.js.
# Packages: ape (needed), DT, pegas (optional: median-joining / haploNet networks).   install.packages(c("ape", "pegas"))

.pg_cols <- function(n) {
  base <- c("#2c7be5", "#e67e22", "#27ae60", "#c0392b", "#8e44ad", "#16a085", "#f1c40f", "#7f8c8d", "#d35400", "#2980b9")
  if (n <= length(base)) base[seq_len(n)] else grDevices::hcl.colors(n, "Dark 3")
}
.pg_progress <- function(session, msg) {
  p <- shiny::Progress$new(session, min = 0, max = 100)
  p$set(0, message = msg, detail = "0 %")
  list(set = function(pct, detail = NULL)
         p$set(value = pct, message = msg, detail = paste0(round(pct), " %", if (!is.null(detail)) paste0(" - ", detail) else "")),
       close = function() p$close())
}
# clean sequence names, drop sequences without any base; at least 3 sequences
.pg_prep_names <- function(aln) {
  if (is.null(aln) || length(aln) < 3) return(NULL)
  nm <- names(aln)
  safe <- make.unique(gsub("[^A-Za-z0-9._-]", "_", nm), sep = "_")
  x <- stats::setNames(toupper(unname(aln)), safe)
  nobase <- nchar(gsub("[-?.NX]", "", x)) == 0
  x <- x[!nobase]
  if (length(x) < 3) return(NULL)
  list(aln = x, dropped = safe[nobase], renamed = sum(safe != nm))
}

.pg_gloss <- c(
  "n" = "Number of sequences", "sites" = "Number of sites used (complete in every sequence)", "s" = "Segregating (variable) sites", "eta" = "Total number of mutations (Eta)",
  "h" = "Number of different haplotypes", "hd" = "Haplotype diversity: the chance that two random sequences are different haplotypes", "sd hd" = "Standard deviation of Hd",
  "k" = "Mean number of pairwise differences between two sequences", "pi" = "Nucleotide diversity: pairwise differences per site",
  "theta hom" = "Theta (4Nu) from the homozygosity of the haplotype frequencies", "theta s" = "Theta from the number of segregating sites (Watterson)",
  "theta k" = "Theta from the expected number of alleles", "theta pi" = "Theta from the mean pairwise differences",
  "tajima d" = "Tajima's D: negative = excess of rare variants (expansion, selection); positive = deficit", "p d" = "P-value of Tajima's D (beta distribution)", "p d sim" = "P-value of Tajima's D from coalescent simulations (as Arlequin)",
  "fs" = "Fu's Fs: strongly negative = excess of haplotypes (expansion); significant when P < 0.02", "p fs" = "P-value of Fu's Fs from simulations",
  "sd k" = "Standard deviation of k", "sd pi" = "Standard deviation of Pi", "sd theta pi" = "Standard deviation of theta(Pi)",
  "f obs" = "Observed homozygosity F of the haplotype frequencies (Ewens-Watterson)", "f mean simulated" = "Mean F of the simulated samples", "p f sim le f obs" = "P(F simulated <= F observed): small = frequencies too even", "p f sim ge f obs" = "P(F simulated >= F observed): small = frequencies too uneven",
  "d.f." = "Degrees of freedom", "sum of squares" = "Sum of squared deviations", "variance components" = "Estimated variance component", "percentage of variation" = "Share of the total variance",
  "fixation index" = "Fixation index (PhiST, PhiSC, PhiCT): how much of the variation lies between populations / groups", "p" = "P-value from permutations", "value" = "Value",
  "haplotype" = "Haplotype name", "count" = "Number of sequences with this haplotype", "frequency" = "Frequency of the haplotype", "sequences" = "Sequences that have this haplotype",
  "d*" = "Fu and Li's D*: negative = excess of singletons (expansion, sweep, purifying selection)", "f*" = "Fu and Li's F*: like D*, also uses the mean pairwise differences", "r2" = "Ramos-Onsins and Rozas' R2: small = excess of singletons (expansion)", "singletons" = "Sites where one sequence differs from all the others", "s biallelic" = "Sites with two nucleotides",
  "population" = "Population (set of sequences)", "group" = "Group of populations (for the hierarchical AMOVA)", "sequence" = "Sequence name")
.pg_tip_js <- function() {
  g <- paste(sprintf("%s:%s", jsonlite_str(names(.pg_gloss)), jsonlite_str(unname(.pg_gloss))), collapse = ",")
  DT::JS(sprintf("function(thead){var g={%s};$(thead).find('th').each(function(){var t=$(this).text().trim().toLowerCase().replace(/[_]+/g,' ');if(g[t])$(this).attr('title',g[t]);});}", g))
}
jsonlite_str <- function(x) paste0("\"", gsub("\"", "\\\\\"", x), "\"")
dtt <- function(...) { w <- DT::datatable(...); w$x$options$headerCallback <- .pg_tip_js(); w }
LENGTHS <- list(c(10, 20, 50, 100, -1), c("10", "20", "50", "100", "All"))      # choices of "show n rows"
.pg_netjs_tag <- function() {
  p <- c("www/sangerly_network.js", "sangerly_network.js", "R/sangerly_network.js")
  p <- p[file.exists(p)]
  if (length(p)) {
    js <- paste(readLines(p[1], warn = FALSE, encoding = "UTF-8"), collapse = "\n")
    return(shiny::tags$script(shiny::HTML(gsub("</script", "<\\/script", js, fixed = TRUE))))
  }
  shiny::tags$script(src = "sangerly_network.js")
}

# ---------------------------------------------------------------------------------------------------
popgenUI <- function(id) {
  ns <- shiny::NS(id)
  dl <- function(x, lab) shiny::downloadButton(ns(x), lab)
  cnd <- function(expr, ...) shiny::conditionalPanel(sprintf(expr, ns("pop_method")), ...)
  chk <- function(cls, lab, val = TRUE) shiny::tags$label(style = "display:flex;align-items:center;gap:6px;margin:0 12px 0 0;font-weight:400",
                                                          shiny::tags$input(type = "checkbox", class = cls, checked = if (val) NA else NULL), lab)
  num <- function(cls, lab, val, w = "62px", step = NULL) shiny::tags$label(style = "display:flex;align-items:center;gap:6px;margin:0 12px 0 0;font-weight:400", lab,
                                                               shiny::tags$input(type = "number", class = paste("form-control form-control-sm", cls), value = val, step = step, style = paste0("width:", w)))
  sel_nv <- function(cls, lab, choices) shiny::tags$label(style = "display:flex;align-items:center;gap:6px;margin:0 12px 0 0;font-weight:400", lab,
    shiny::tags$select(class = paste("form-select form-select-sm", cls), style = "width:auto",
      lapply(seq_along(choices), function(i) shiny::tags$option(value = unname(choices)[i], names(choices)[i]))))
  btn <- function(cls, lab, primary = FALSE) shiny::tags$button(type = "button", class = paste("btn btn-sm", if (primary) "btn-primary" else "btn-light border", cls), lab)
  shiny::div(class = "container-fluid",
    .pg_netjs_tag(),
    shiny::tags$style(shiny::HTML(".pg-steps{display:flex;gap:6px;flex-wrap:wrap;margin-bottom:8px}.pg-chip{flex:1;min-width:150px;border-radius:8px;padding:6px 10px;border:1px solid;font-size:12.5px;line-height:1.35}.pg-ok{background:#e8f7ef;border-color:#9bd8b5}.pg-todo{background:#f3f5f7;border-color:#cfd6dd;color:#52606d}.pg-stale{background:#fff4e0;border-color:#f0c070}.pg-warn{background:#fffbe6;border-color:#ecd56a}details.pg-help{border:1px solid #cfe3f5;background:#f4f9fe;border-radius:8px;padding:6px 10px;margin-bottom:8px}details.pg-help summary{cursor:pointer}.nv-box{background:#fff;border:1px solid #dfe3e8;border-radius:8px;padding:8px;margin-bottom:8px}.nv-bar{display:flex;flex-wrap:wrap;gap:6px 4px;align-items:center;margin-bottom:6px;font-size:12.5px}.nv-canvas{display:none;border:1px solid #dfe3e8;border-radius:6px;cursor:default}.nv-empty{padding:30px;text-align:center;color:#6b7785}")),
    if (!requireNamespace("ape", quietly = TRUE)) shiny::div(class = "alert alert-danger", "The 'ape' package is required: install.packages(\"ape\")"),
    shiny::fluidRow(
      shiny::column(4,
        .card("1. Alignment",
          shiny::radioButtons(ns("src"), NULL,
            c("MSA tab: kept columns (after trimming / editing)" = "trim", "MSA tab: whole alignment" = "all", "An aligned FASTA file" = "file")),
          shiny::conditionalPanel(sprintf("input['%s'] == 'file'", ns("src")),
            shiny::fileInput(ns("fa"), NULL, buttonLabel = "Open aligned FASTA...", placeholder = "", accept = c(".fa", ".fasta", ".fas", ".fna", ".txt"))),
          shiny::uiOutput(ns("aln_info")),
          shiny::uiOutput(ns("msg_aln")),
          shiny::actionButton(ns("load_example"), "Load example data", icon = shiny::icon("flask"), class = "btn-sm btn-outline-secondary"),
          shiny::helpText("60 ITS sequences of Durio zibethinus with a table of locations, to try the program.")),
        .card("2. Populations and groups (sequence sets)",
          shiny::uiOutput(ns("msg_pop")),
          shiny::radioButtons(ns("pop_src"), "Where do the populations come from? (choose one)",
            c("From the sequence names (a rule)" = "names", "From a CSV / Excel file" = "file", "By hand (select sequences, name the set)" = "manual")),
          shiny::conditionalPanel(sprintf("input['%s'] == 'names'", ns("pop_src")),
            shiny::selectInput(ns("pop_method"), "Population taken from the sequence name",
              c("Text before a separator" = "prefix", "First characters" = "chars", "Regular expression (first group)" = "regex", "One population" = "all")),
            cnd("input['%s'] == 'prefix'", shiny::textInput(ns("pop_sep"), "Separator", "_"), shiny::numericInput(ns("pop_nfield"), "Number of fields", 1, min = 1, step = 1)),
            cnd("input['%s'] == 'chars'", shiny::numericInput(ns("pop_nchars"), "Number of characters", 3, min = 1, step = 1)),
            cnd("input['%s'] == 'regex'", shiny::textInput(ns("pop_regex"), "Expression", "^([A-Za-z]+)")),
            shiny::actionButton(ns("apply_rule"), "Apply this rule")),
          shiny::conditionalPanel(sprintf("input['%s'] == 'file'", ns("pop_src")),
            shiny::div(style = "display:flex;gap:8px;flex-wrap:wrap;align-items:center",
              shiny::fileInput(ns("pop_csv"), NULL, buttonLabel = "Open CSV / Excel...", placeholder = "", accept = c(".csv", ".tsv", ".txt", ".xlsx", ".xls"), width = "260px"),
              shiny::downloadButton(ns("pop_tpl"), "Table with my sequence names (CSV)")),
            shiny::uiOutput(ns("csv_state")),
            shiny::helpText("Columns: Sequence (the full name, or just an accession such as MT298407.1), Population and, optionally, Group. Other columns are ignored. ",
                            "Tip: download the table with your sequence names, fill the Population column in Excel, and open it here.")),
          shiny::conditionalPanel(sprintf("input['%s'] == 'manual'", ns("pop_src")),
            shiny::helpText("Select rows in the table below (Ctrl / Shift + click), then give the population and press 'Update population'. The Population and Group cells can also be edited directly.")),
          shiny::uiOutput(ns("pop_srclabel")),
          shiny::h6("Sequences and their sets"),
          DT::dataTableOutput(ns("seq_tbl")),
          shiny::conditionalPanel(sprintf("input['%s'] == 'manual'", ns("pop_src")),
            shiny::div(style = "display:flex;gap:6px;flex-wrap:wrap;margin-top:6px",
              shiny::actionButton(ns("sel_all"), "Select all shown"), shiny::actionButton(ns("sel_none"), "Unselect all")),
            shiny::div(style = "display:flex;gap:6px;align-items:end;flex-wrap:wrap;margin-top:6px",
              shiny::selectizeInput(ns("set_pop"), "Population of the selected", choices = NULL, width = "210px",
                                    options = list(create = TRUE, persist = FALSE, placeholder = "pick one or type a new name")),
              shiny::actionButton(ns("assign_pop"), "Update population"),
              shiny::selectizeInput(ns("set_grp"), "Group of the selected", choices = NULL, width = "210px",
                                    options = list(create = TRUE, persist = FALSE, placeholder = "pick one or type a new name")),
              shiny::actionButton(ns("assign_grp"), "Update group"))),
          shiny::h6(style = "margin-top:12px", "Populations (edit the Group, or rename a population)"),
          DT::dataTableOutput(ns("grp_tbl")),
          shiny::helpText("Group = region, for the hierarchical AMOVA (every population in one group). Leave 'all' for a simple AMOVA.")),
        .card("3. Settings",
          shiny::h6("Arlequin settings file"),
          shiny::fileInput(ns("ars"), NULL, buttonLabel = "Open Arlequin settings (.ars)...", placeholder = "", accept = ".ars", width = "100%"),
          shiny::uiOutput(ns("ars_msg")),
          shiny::uiOutput(ns("msg_set")),
          shiny::uiOutput(ns("ars_save")),
          shiny::helpText("Arlequin: Settings > Save settings. The values it holds (permutations, distance method, gamma, significance level, neutrality tests, exact test) fill the form below, ",
                          "so the AMOVA runs with the same settings as in Arlequin. After changing the form you can save the file again for Arlequin."),
          shiny::h6("Diversity and neutrality (Arlequin: molecular diversity, neutrality tests)"),
          shiny::radioButtons(ns("div_gap"), "Sites used for diversity and neutrality (gaps, N, IUPAC codes)",
            c("Sites complete in all the sequences (one common set of sites)" = "all",
              "Sites complete in each population separately (DnaSP default; Arlequin: 'allowed level of missing data' below, per sample)" = "subset",
              "All sites; gaps ignored pair by pair for k and pi (DnaSP: pairwise deletion)" = "pairwise")),
          shiny::checkboxInput(ns("div_inv"), "Per-site values (pi) divided by the number of variable sites of the whole alignment (as Arlequin does with a DnaSP file in which the invariable sites were removed)", FALSE),
          shiny::helpText("The two rules above only change the per-site values and the sites counted; with an alignment without gaps the choice 'complete in all' and 'complete in each population' give the same numbers. ",
                          "To reproduce a DnaSP + Arlequin analysis: use the DnaSP rule you chose there, and tick the box above if you exported the haplotype file with 'invariable sites: removed'."),
          shiny::numericInput(ns("nsim"), "Simulated samples for the p-value of Fu's Fs (0 = none)", 1000, min = 0, step = 500),
          shiny::numericInput(ns("ew_nsim"), "Ewens-Watterson test: simulated samples (0 = none)", 1000, min = 0, step = 500),
          shiny::h6("Genetic differentiation and gene flow (DnaSP)"),
          shiny::radioButtons(ns("gf_gap"), "Sites with alignment gaps are...",
            c("excluded" = "excluded", "considered (gap as the fifth state)" = "fifth", "excluded only in pairwise comparisons" = "pairwise")),
          shiny::selectInput(ns("gf_genome"), "Genome", c("Diploid-Autosome" = "diploid", "Haploid (mtDNA, cpDNA)" = "haploid", "X-linked" = "X", "Y-linked" = "Y")),
          shiny::numericInput(ns("gf_nperm"), "Permutation test: replicates (0 = none)", 1000, min = 0, step = 500),
          shiny::numericInput(ns("gf_seed"), "Pseudorandom number seed (empty = random)", NA, step = 1),
          shiny::uiOutput(ns("gf_pops_ui")),
          shiny::numericInput(ns("ex_steps"), "Exact test of differentiation: Markov chain steps (0 = none)", 100000, min = 0, step = 50000),
          shiny::numericInput(ns("ex_demem"), "Dememorization steps", 10000, min = 0, step = 5000),
          shiny::checkboxInput(ns("ex_pair"), "Also for every pair of populations (slower)", FALSE),
          shiny::h6("AMOVA and pairwise FST (Arlequin)"),
          shiny::selectInput(ns("dmodel"), "Distance method", c("Pairwise difference" = "N", "Kimura 2-parameter" = "K80", "Tamura-Nei" = "TN93", "Jukes-Cantor" = "JC69")),
          shiny::numericInput(ns("gamma"), "Gamma a value (0 = none)", 0, min = 0, step = 0.1),
          shiny::numericInput(ns("nperm"), "No. of permutations (AMOVA)", 1000, min = 100, step = 500),
          shiny::numericInput(ns("nperm_pair"), "No. of permutations (pairwise FST)", 200, min = 50, step = 100),
          shiny::numericInput(ns("alpha"), "Significance level", 0.05, min = 0.001, max = 0.5, step = 0.01),
          shiny::checkboxInput(ns("slatkin"), "Show Slatkin's linearized FST", FALSE),
          shiny::checkboxInput(ns("reyn"), "Show Reynolds' distance", FALSE),
          shiny::checkboxInput(ns("pipop"), "Show pairwise differences (pi) between populations", FALSE),
          shiny::numericInput(ns("miss"), "Allowed level of missing data per site (0 = none; AMOVA, FST, pi)", 0, min = 0, max = 0.5, step = 0.01),
          shiny::h6("Structure (K clusters, as STRUCTURE / ADMIXTURE)"),
          shiny::selectInput(ns("st_method"), "Method",
            c("Maximum likelihood, EM (fast; as ADMIXTURE / LEA)" = "em", "Bayesian MCMC, Gibbs sampler (as STRUCTURE; slower)" = "mcmc")),
          shiny::selectInput(ns("st_model"), "Model",
            c("Admixture: every sequence can have several ancestries" = "admix", "No admixture: every sequence in one cluster (one non-recombining locus)" = "mix")),
          shiny::fluidRow(shiny::column(6, shiny::numericInput(ns("st_kmin"), "K from", 1, min = 1, step = 1)),
                          shiny::column(6, shiny::numericInput(ns("st_kmax"), "K to", 5, min = 2, step = 1))),
          shiny::numericInput(ns("st_rep"), "Runs per K (the best likelihood is shown)", 5, min = 1, step = 1),
          shiny::conditionalPanel(sprintf("input['%s'] != 'mcmc'", ns("st_method")),
            shiny::numericInput(ns("st_cv"), "Cross-entropy: % of the genotypes hidden (0 = none)", 5, min = 0, max = 30, step = 1),
            shiny::numericInput(ns("st_iter"), "Maximum iterations per run", 300, min = 50, step = 100)),
          shiny::conditionalPanel(sprintf("input['%s'] == 'mcmc'", ns("st_method")),
            shiny::numericInput(ns("st_burn"), "Length of the burn-in (sweeps)", 2000, min = 100, step = 500),
            shiny::numericInput(ns("st_nrep"), "Number of MCMC sweeps after the burn-in", 3000, min = 100, step = 1000),
            shiny::numericInput(ns("st_alpha"), "Alpha (the starting value when it is inferred)", 1, min = 0.01, max = 10, step = 0.1),
            shiny::checkboxInput(ns("st_alpha_inf"), "Infer alpha from the data (as STRUCTURE)", TRUE),
            shiny::helpText("STRUCTURE itself uses 10 000 burn-in and 20 000 sweeps; here fewer sweeps are the default to keep the program fast. Check that the ln likelihood is stable between runs of the same K.")),
          shiny::numericInput(ns("st_seed"), "Seed (empty = random)", NA, step = 1)),
        .card("4. Run",
          shiny::div(style = "display:flex;gap:8px;flex-wrap:wrap",
            shiny::actionButton(ns("run_div"), "Haplotypes, diversity, differentiation", icon = shiny::icon("play"), class = "btn-primary"),
            shiny::actionButton(ns("run_amova"), "AMOVA and pairwise FST", icon = shiny::icon("play"), class = "btn-primary"),
            shiny::actionButton(ns("run_struct"), "Structure (K clusters)", icon = shiny::icon("play"), class = "btn-primary")),
          shiny::uiOutput(ns("msg_run")),
          shiny::hr(),
          shiny::div(style = "display:flex;gap:8px;flex-wrap:wrap;align-items:center",
            dl("proj_save", "Save settings and sets (.rds)"),
            shiny::fileInput(ns("proj_load"), NULL, buttonLabel = "Open saved...", placeholder = "", accept = ".rds", width = "230px")),
          shiny::helpText("This file keeps only the sets and the settings of this step."))),
      shiny::column(8,
        shiny::helpText(shiny::icon("graduation-cap"), " Step-by-step tutorial and how the results compare with Arlequin / DnaSP: menu More > Tutorial > Population structure."),
        shiny::uiOutput(ns("steps")),
        shiny::uiOutput(ns("msg")),
        shiny::tabsetPanel(id = ns("tabs"),
          shiny::tabPanel("Haplotype", shiny::br(),
            shiny::h5("Haplotype table"), DT::dataTableOutput(ns("hap_tbl")),
            shiny::h5(style = "margin-top:14px", "Haplotype distribution (as DnaSP)"), shiny::verbatimTextOutput(ns("hap_text")),
            shiny::h5(style = "margin-top:14px", "Haplotype / DNA sequences data file (as DnaSP: Generate)"),
            shiny::fluidRow(
              shiny::column(4, shiny::radioButtons(ns("hf_format"), "Generate...",
                c("NEXUS haplotype data file" = "nexus", "Arlequin haplotype list (.arp)" = "arlequin", "FASTA haplotypes" = "fasta"))),
              shiny::column(4, shiny::radioButtons(ns("hf_gaps"), "Sites with gaps/missing...", c("Not considered" = "no", "Considered" = "yes"))),
              shiny::column(4, shiny::radioButtons(ns("hf_invar"), "Invariable sites...", c("Removed" = "removed", "Included" = "included")))),
            shiny::div(style = "display:flex;gap:6px;flex-wrap:wrap", dl("dl_hapfile", "Generate haplotype file"), dl("dl_hap", "Haplotype table (CSV)")),
            shiny::h5(style = "margin-top:14px", "Sequence files for other programs"),
            shiny::div(style = "display:flex;gap:6px;flex-wrap:wrap",
              dl("ex_fasta", "FASTA"), dl("ex_nexus", "NEXUS (with TaxSets)"), dl("ex_phylip", "PHYLIP"), dl("ex_phymap", "PHYLIP names (CSV)"), dl("ex_mega", "MEGA"),
              dl("ex_popart", "PopART (haplotypes + traits)"), dl("ex_pop", "Population table (CSV)")),
            shiny::helpText("DnaSP reads FASTA, NEXUS (with the sets defined by TaxSet), PHYLIP, MEGA and Arlequin haplotype files. ",
                            "The Roehl data file (.rdf) of the NETWORK program is not produced yet: DnaSP can make it from the NEXUS or Arlequin file above.")),
          shiny::tabPanel("Network", shiny::br(),
            shiny::div(style = "display:flex;gap:8px;flex-wrap:wrap;align-items:end",
              shiny::radioButtons(ns("net_method"), "Network",
                c("Minimum spanning network" = "msn", "Median-joining (as NETWORK / DnaSP)" = "mjn", "Haplotype network (pegas::haploNet)" = "haplo"), inline = TRUE),
              shiny::conditionalPanel(sprintf("input['%s'] == 'mjn'", ns("net_method")), shiny::numericInput(ns("net_eps"), "epsilon (0 = classical)", 0, min = 0, step = 1, width = "150px")),
              shiny::actionButton(ns("run_net"), "Build network", icon = shiny::icon("play"), class = "btn-primary")),
            shiny::uiOutput(ns("msg_net")),
            shiny::div(style = "margin-bottom:6px", shiny::actionButton(ns("net_save"), "Save SVG + PNG in the output folder", icon = shiny::icon("folder-open"), class = "btn-sm btn-outline-secondary"),
                       shiny::span(class = "text-muted", style = "margin-left:8px;font-size:12.5px", "(the folder is chosen in the Result tab; the SVG / PNG buttons below download through the browser)")),
            shiny::div(id = ns("net"), class = "nv-box",
              shiny::div(class = "nv-bar",
                chk("nv-prop", "Nodes proportional to frequencies"), chk("nv-mut", "Display mutated position"), chk("nv-ticks", "Mutated position lines", FALSE),
                chk("nv-names", "Display node names"), chk("nv-med", "Show median vectors"), chk("nv-mednames", "Median vector names")),
              shiny::div(class = "nv-bar",
                num("nv-min", "Node size min.", 10), num("nv-max", "max.", 50), num("nv-fixed", "Fixed size", 14), num("nv-medsize", "Median size", 6),
                num("nv-minfreq", "Show haplotypes with n >=", 1)),
              shiny::div(class = "nv-bar",
                num("nv-font", "Name font", 11), num("nv-fontmed", "Median vector font", 10), num("nv-fontmut", "Mutated position font", 10),
                chk("nv-fontscale", "Text scales with the zoom", FALSE)),
              shiny::div(class = "nv-bar",
                sel_nv("nv-lblpos", "Names placed", c("left of the node" = "left", "right of the node" = "right", "above the node" = "top", "below the node" = "bottom")),
                sel_nv("nv-mutpos", "Mutated positions placed", c("along the link, one number per position" = "along", "beside the link (one text)" = "beside", "on the link (one text)" = "center", "near the first node" = "start", "near the second node" = "end")),
                num("nv-mutoff", "Distance from the link (px)", 9), num("nv-maxpos", "Max. positions shown (0 = all)", 0), chk("nv-owon", "Pie outline"), num("nv-ow", "thickness (px)", 1.8, step = 0.2), chk("nv-sizeleg", "Size legend (circles)")),
              shiny::div(class = "nv-bar",
                num("nv-h", "Height (px)", 560, "70px"), btn("nv-zin", "Zoom +"), btn("nv-zout", "Zoom -"), btn("nv-fit", "Fit"), btn("nv-relayout", "Re-arrange"), btn("nv-resetlbl", "Reset texts"), btn("nv-undo", "Undo"),
                btn("nv-png", "PNG", TRUE), btn("nv-svg", "SVG", TRUE), shiny::span(class = "nv-info", style = "color:#52606d;margin-left:8px")),
              shiny::div(class = "nv-bar", shiny::span(style = "font-weight:600", "Population colours:"), shiny::span(class = "nv-colors")),
              shiny::tags$canvas(class = "nv-canvas"), shiny::div(class = "nv-empty", "Run 'Haplotypes, diversity, differentiation', then press 'Build network'.")),
            shiny::helpText("Drag a node to move it; drag a haplotype / median-vector name or the mutated positions to move that text; drag the background to pan; wheel = zoom. ",
                            "Double-click a node, a name or the mutated positions of a link to edit the text (empty text on a link = automatic again). ",
                            "Numbers on a link = alignment columns that differ, each one on its own place on the link (with a short mark when 'Mutated position lines' is ticked); in a crowded link the numbers alternate on both sides. ",
                            "Median vectors (median-joining) are unobserved haplotypes inferred by the algorithm; the program removes the ones that are left with two links or fewer, as NETWORK does. The circles at the bottom left show the size of 1 sample and of the round number of samples.")),
          shiny::tabPanel("Map", shiny::br(),
            shiny::p("A pie chart of the haplotypes of every population at its place on the map (as PopART): one slice per haplotype, the area of the pie grows with the number of sequences. ",
                     "Run 'Haplotypes, diversity, differentiation' first, then load the coordinates."),
            shiny::fluidRow(
              shiny::column(4, shiny::fileInput(ns("map_coord"), NULL, buttonLabel = "Coordinates (CSV / Excel)...", placeholder = "Population, Latitude, Longitude", accept = c(".csv", ".txt", ".tsv", ".xlsx", ".xls")),
                             dl("map_tpl", "Coordinates template (CSV)")),
              shiny::column(3, shiny::selectInput(ns("map_base"), "Base map", c("Satellite (Esri World Imagery)" = "sat", "Topographic (Esri)" = "topo", "Ocean (Esri)" = "ocean", "Streets (OpenStreetMap)" = "street",
                                                                               "Light grey (Carto)" = "light", "Dark (Carto)" = "dark", "Terrain (OpenTopoMap)" = "terrain", "Coastlines (offline, no internet needed)" = "offline", "Plain colour" = "none"), "sat"),
                             shiny::textInput(ns("map_title"), "Title", "")),
              shiny::column(2, shiny::sliderInput(ns("map_pie"), "Largest pie (% of width)", 1, 15, 5, 0.5), shiny::sliderInput(ns("map_h"), "Height (px)", 300, 900, 560, 20)),
              shiny::column(3, shiny::checkboxInput(ns("map_lab"), "Write the population names", TRUE), shiny::checkboxInput(ns("map_hleg"), "Legend of the haplotypes", TRUE),
                             shiny::checkboxInput(ns("map_sleg"), "Legend of the size (circles)", TRUE), shiny::checkboxInput(ns("map_scale"), "Scale bar", TRUE))),
            shiny::fluidRow(
              shiny::column(2, shiny::textInput(ns("map_land"), "Land colour (offline / plain)", "#e8e4d8")),
              shiny::column(2, shiny::textInput(ns("map_sea"), "Sea colour (offline / plain)", "#dbe9f4")),
              shiny::column(8, shiny::textInput(ns("map_cols"), "Haplotype colours (optional: #hex, comma-separated, in the order Hap_1, Hap_2...)", ""))),
            shiny::uiOutput(ns("msg_map")),
            shiny::div(id = ns("hmap"), style = "min-height:60px"),
            shiny::helpText("The CSV (or Excel) has one row per population with the columns Population, Latitude and Longitude in decimal degrees (south and west are negative); the population names must be the same as in section 2. ",
                            "The template downloads the populations of the last run. The same coordinates are used by the map of the Structure tab. ",
                            "Satellite, streets and the other base maps are pictures (tiles) loaded from public servers while you look at the map, so they need internet; without internet choose 'Coastlines (offline)'. ",
                            "Drag the map to move it, use the wheel to zoom, drag a pie to a better place. The buttons PNG and SVG above the map save what you see (the SVG has the base map as one picture and the pies, names and legends as editable objects). ",
                            "Please keep the credit of the base map that is written in the corner of the figure when you publish it.")),
          shiny::tabPanel("AMOVA", shiny::br(), shiny::uiOutput(ns("amova_head")), DT::dataTableOutput(ns("amova_tbl")), DT::dataTableOutput(ns("phi_tbl")),
            shiny::h5(style = "margin-top:14px", "Pairwise FST (below the diagonal) and p-value (above); + = significant"), DT::dataTableOutput(ns("pw_tbl")), shiny::uiOutput(ns("amova_extra")),
            shiny::div(style = "margin-top:8px;display:flex;gap:6px;flex-wrap:wrap", dl("dl_amova", "AMOVA (CSV)"), dl("dl_pw", "Pairwise FST (CSV)"))),
          shiny::tabPanel("Structure", shiny::br(),
            shiny::uiOutput(ns("st_head")),
            shiny::uiOutput(ns("msg_struct")),
            shiny::fluidRow(
              shiny::column(4, shiny::sliderInput(ns("st_llh"), "Height of the likelihood figure", 150, 600, 270, 10)),
              shiny::column(4, shiny::sliderInput(ns("st_mh"), "Height of the map", 250, 900, 460, 10))),
            shiny::fluidRow(
              shiny::column(3, shiny::uiOutput(ns("st_k_ui"))),
              shiny::column(3, shiny::selectInput(ns("st_mode"), "Show",
                c("K by population and K by individual" = "both", "K by population" = "pop", "K by individual" = "ind"))),
              shiny::column(3, shiny::selectInput(ns("st_palette"), "Colours", c("As in the paper" = "figure", "Bright" = "bright", "Soft" = "soft", "Automatic" = "auto"))),
              shiny::column(3, shiny::sliderInput(ns("st_h"), "Height of the figure", 240, 900, 460, 10))),
            shiny::fluidRow(
              shiny::column(4, shiny::textInput(ns("st_title"), "Title (for example the region)", "", placeholder = "Region EastAsia--Oceania")),
              shiny::column(4, shiny::textInput(ns("st_order"), "Order of the populations (comma separated)", "", placeholder = "empty = as in the population table")),
              shiny::column(4, shiny::textInput(ns("st_colors"), "Own colours, first cluster first (hex)", "", placeholder = "#EE1D52, #19F0F0, ..."))),
            shiny::fluidRow(
              shiny::column(4, shiny::sliderInput(ns("st_cex"), "Text size", 0.4, 1.6, 0.85, 0.05)),
              shiny::column(4, shiny::checkboxInput(ns("st_sort"), "Sort the individuals of a population by ancestry", FALSE)),
              shiny::column(4, shiny::checkboxInput(ns("st_names"), "Write the name of every sequence", FALSE))),
            shiny::plotOutput(ns("st_plot"), height = "auto"),
            shiny::div(style = "margin:8px 0 14px;display:flex;gap:6px;flex-wrap:wrap", dl("st_png", "Figure (PNG)"), dl("st_pdf", "Figure (PDF)"), dl("st_svg", "Figure (SVG)"),
                       dl("st_dl_q", "Q of every sequence (CSV)"), dl("st_dl_qpop", "Q of every population (CSV)"), dl("st_dl_ll", "Likelihood by K (CSV)")),
            shiny::h5("Which K? Cross-entropy (as LEA / ADMIXTURE), likelihood, BIC and Evanno's delta K"), shiny::plotOutput(ns("st_ll_plot"), height = "auto"),
            shiny::helpText("Cross-entropy: look for the K where the curve reaches its minimum or levels off (the circle marks the minimum). Evanno's delta K is the peak of the curve (best for 2 or more runs and K starting at 1 or 2)."), DT::dataTableOutput(ns("st_ll_tbl")),
            shiny::h5(style = "margin-top:14px", "Ancestry of every population (mean of its sequences)"), DT::dataTableOutput(ns("st_qpop_tbl")),
            shiny::h5(style = "margin-top:14px", "Map of the populations (pie chart of the mean ancestry)"),
            shiny::fluidRow(
              shiny::column(5, shiny::fileInput(ns("st_coord"), NULL, buttonLabel = "Coordinates (CSV)...", placeholder = "Population, Latitude, Longitude", accept = c(".csv", ".txt", ".tsv"))),
              shiny::column(3, shiny::sliderInput(ns("st_pie"), "Pie size (% of the map width)", 1, 12, 3, 0.5)),
              shiny::column(4, shiny::checkboxInput(ns("st_pie_n"), "Pie size grows with the number of sequences", FALSE), shiny::checkboxInput(ns("st_pie_lab"), "Write the population names", TRUE))),
            shiny::plotOutput(ns("st_map"), height = "auto"),
            shiny::div(style = "margin:6px 0 10px", dl("st_map_png", "Map (PNG)")),
            shiny::helpText("CSV with the columns Population, Latitude and Longitude (decimal degrees; one row per population; names as in section 2). ",
                            "Coastlines are added when the 'maps' package is installed: install.packages(\"maps\")."),
            shiny::h5(style = "margin-top:14px", "Ancestry of every sequence (Q)"), DT::dataTableOutput(ns("st_q_tbl")),
            shiny::div(style = "margin-top:8px", shiny::fileInput(ns("st_q_csv"), NULL, buttonLabel = "Draw a Q file from STRUCTURE / ADMIXTURE...", placeholder = "", accept = c(".csv", ".txt", ".tsv", ".Q"), width = "330px")),
            shiny::helpText("The clusters are found from the variable sites of the alignment, either with a maximum-likelihood EM algorithm (fast; the model of ADMIXTURE / LEA) or with the Bayesian Gibbs sampler of STRUCTURE (Pritchard et al. 2000; Q = posterior mean). Several random starts; the run with the best likelihood is drawn. ",
                            "K is chosen with the cross-entropy of hidden genotypes (as LEA and ADMIXTURE), BIC and Evanno's delta K. ",
                            "The sites are treated as independent, as in STRUCTURE's no-linkage model. For one gene or mtDNA the sites are linked, so the likelihood is over-confident: ",
                            "use the figure to describe the groups, choose K with delta K and BIC together with the biology, and prefer the 'no admixture' model for a single locus. ",
                            "For a formal analysis run STRUCTURE / ADMIXTURE on your SNP data and load their Q table here (first column: the sequence name; then one column per cluster): the same figure is drawn. ",
                            "The population of every sequence is taken from section 2, so changing a set re-draws the figure without running again.")),
          shiny::tabPanel("Sliding window", shiny::br(),
            shiny::p("Diversity along the alignment: a window of a fixed number of alignment columns moves along the sequences. Only the sites without gap / N / IUPAC code (the ones used in the Result tab) count; the per-site values are divided by the number of such sites in the window. Run 'Haplotypes, diversity, differentiation' first."),
            shiny::div(style = "display:flex;gap:12px;flex-wrap:wrap;align-items:end",
              shiny::numericInput(ns("sl_w"), "Window (columns)", 100, min = 10, step = 10, width = "140px"),
              shiny::numericInput(ns("sl_s"), "Step (columns)", 25, min = 1, step = 5, width = "130px"),
              shiny::numericInput(ns("sl_min"), "Minimum sites in a window", 10, min = 2, step = 1, width = "190px"),
              shiny::selectInput(ns("sl_stat"), "Statistic to draw", c("Nucleotide diversity (Pi, per site)" = "Pi", "Watterson's theta (per site)" = "ThetaW", "Tajima's D" = "TajimaD", "Segregating sites (S)" = "S", "Haplotype diversity (Hd)" = "Hd", "Haplotypes (h)" = "h"), width = "260px"),
              shiny::checkboxInput(ns("sl_pop"), "Also every population", TRUE),
              shiny::actionButton(ns("run_sl"), "Run the sliding window", icon = shiny::icon("play"), class = "btn-primary", style = "margin-bottom:15px")),
            shiny::uiOutput(ns("msg_slide")),
            shiny::sliderInput(ns("sl_h"), "Height of the figure", 200, 700, 340, 10, width = "300px"),
            shiny::plotOutput(ns("sl_plot"), height = "auto"),
            shiny::div(style = "margin:8px 0;display:flex;gap:6px;flex-wrap:wrap", dl("sl_png", "Figure (PNG)"), dl("sl_csv", "Table (CSV)")),
            DT::dataTableOutput(ns("sl_tbl"))),
          shiny::tabPanel("More tests", shiny::br(),
            shiny::p("Further neutrality and structure tests of DnaSP / Arlequin. Run 'Haplotypes, diversity, differentiation' first (Strobeck's S, Chakraborty's test and Achaz's Y* are in the Result tab). Sub-tab 'Phasing' and 'Hardy-Weinberg' work on the sequences with IUPAC codes (R, Y, S, W, K, M) of diploid individuals."),
            shiny::uiOutput(ns("msg_more")),
            shiny::tabsetPanel(
              shiny::tabPanel("SFS and linkage (LD)", shiny::br(),
                shiny::div(style = "display:flex;gap:12px;flex-wrap:wrap;align-items:end",
                  shiny::uiOutput(ns("mt_unit_ui")),
                  shiny::numericInput(ns("ld_max"), "Maximum number of sites for LD", 300, min = 20, step = 50, width = "200px"),
                  shiny::checkboxInput(ns("ld_fisher"), "Fisher's exact test for every pair", TRUE),
                  shiny::actionButton(ns("run_ld"), "Run SFS and LD", icon = shiny::icon("play"), class = "btn-primary", style = "margin-bottom:15px")),
                shiny::h6("Site frequency spectrum (biallelic sites, folded)"), shiny::plotOutput(ns("sfs_plot"), height = "260px"), DT::dataTableOutput(ns("sfs_tbl")),
                shiny::h6("Linkage disequilibrium"), DT::dataTableOutput(ns("ld_sum")), shiny::plotOutput(ns("ld_plot"), height = "320px"),
                shiny::div(style = "margin:8px 0;display:flex;gap:6px;flex-wrap:wrap", dl("ld_csv", "Pairs (CSV)"), dl("sfs_csv", "SFS (CSV)"), dl("ld_png", "LD figure (PNG)")),
                DT::dataTableOutput(ns("ld_tbl")),
                shiny::helpText("D, D' (Lewontin) and r^2 between all pairs of biallelic sites (sites with gap / N / IUPAC left out). ZnS = mean r^2 of all pairs (Kelly 1997), Za = mean r^2 of adjacent sites, ZZ = Za - ZnS (Rozas et al. 2001). ",
                                "Rm = minimum number of recombination events (Hudson and Kaplan 1985, four-gamete test). The SFS is compared with the neutral expectation for a constant population (theta from S).")),
              shiny::tabPanel("Mantel test", shiny::br(),
                shiny::div(style = "display:flex;gap:12px;flex-wrap:wrap;align-items:end",
                  shiny::selectInput(ns("mn_gen"), "Genetic distance", c("Pairwise FST" = "fst", "FST / (1 - FST) (Rousset 1997)" = "rousset"), width = "260px"),
                  shiny::selectInput(ns("mn_geo"), "Geographic distance", c("km" = "km", "ln(km)" = "ln"), width = "160px"),
                  shiny::numericInput(ns("mn_perm"), "Permutations", 9999, min = 99, step = 1000, width = "130px"),
                  shiny::actionButton(ns("run_mn"), "Run the Mantel test", icon = shiny::icon("play"), class = "btn-primary", style = "margin-bottom:15px")),
                shiny::verbatimTextOutput(ns("mn_text")), shiny::plotOutput(ns("mn_plot"), height = "340px"),
                shiny::div(style = "margin:8px 0", dl("mn_png", "Figure (PNG)")),
                shiny::helpText("Isolation by distance: correlation between the pairwise FST (tab AMOVA, run it first) and the great-circle distance between the populations (coordinates loaded in the Map tab). The p-value comes from permuting the populations of one matrix. At least 4 populations are needed.")),
              shiny::tabPanel("Hardy-Weinberg", shiny::br(),
                shiny::div(style = "display:flex;gap:12px;flex-wrap:wrap;align-items:end",
                  shiny::numericInput(ns("hw_min"), "Minimum individuals per site", 5, min = 3, step = 1, width = "190px"),
                  shiny::actionButton(ns("run_hw"), "Run the Hardy-Weinberg test", icon = shiny::icon("play"), class = "btn-primary", style = "margin-bottom:15px")),
                DT::dataTableOutput(ns("hw_tbl")), shiny::div(style = "margin:8px 0", dl("hw_csv", "Table (CSV)")),
                shiny::helpText("For diploid individuals (Sanger sequences with IUPAC codes at the heterozygous sites): every site with two bases and its heterozygous code is tested with the exact test of Wigginton, Cutler and Abecasis (2005); Fis = 1 - Ho / He. Sites are not independent when they are on one gene. Sequences without any heterozygous code (haploid data, mtDNA, cpDNA) have no genotypes and nothing is tested.")),
              shiny::tabPanel("Phasing (EM)", shiny::br(),
                shiny::div(style = "display:flex;gap:12px;flex-wrap:wrap;align-items:end",
                  shiny::numericInput(ns("em_max"), "Maximum heterozygous sites per individual", 12, min = 2, max = 16, step = 1, width = "250px"),
                  shiny::actionButton(ns("run_em"), "Infer the haplotypes", icon = shiny::icon("play"), class = "btn-primary", style = "margin-bottom:15px")),
                shiny::verbatimTextOutput(ns("em_text")), DT::dataTableOutput(ns("em_freq")), shiny::br(), DT::dataTableOutput(ns("em_ind")),
                shiny::div(style = "margin:8px 0;display:flex;gap:6px;flex-wrap:wrap", dl("em_fa", "Phased sequences (FASTA, 2 per individual)"), dl("em_csv", "Haplotypes (CSV)"), dl("em_ind_csv", "Individuals (CSV)")),
                shiny::helpText("Expectation-maximisation algorithm of Excoffier and Slatkin (1995): the heterozygous sites (IUPAC codes) of each individual are resolved into the two haplotypes that explain the data best, given the haplotype frequencies of the whole sample (as the 'haplotype estimation' of Arlequin and PHASE-like EM in DnaSP). ",
                                "Only columns made of A, C, G, T and two-base IUPAC codes are used. Save the phased FASTA and load it again as an alignment (population by name: ID_a, ID_b) for haplotype networks and diversity.")),
              shiny::tabPanel("Fst outliers (FDIST)", shiny::br(),
                shiny::div(style = "display:flex;gap:12px;flex-wrap:wrap;align-items:end",
                  shiny::numericInput(ns("fo_sim"), "Simulated loci", 10000, min = 1000, step = 1000, width = "150px"),
                  shiny::numericInput(ns("fo_mac"), "Minimum minor-allele count", 2, min = 1, step = 1, width = "190px"),
                  shiny::numericInput(ns("fo_seed"), "Seed (optional)", NA, width = "120px"),
                  shiny::actionButton(ns("run_fo"), "Detect outlier sites", icon = shiny::icon("play"), class = "btn-primary", style = "margin-bottom:15px")),
                shiny::verbatimTextOutput(ns("fo_text")), shiny::plotOutput(ns("fo_plot"), height = "380px"),
                shiny::div(style = "margin:8px 0;display:flex;gap:6px;flex-wrap:wrap", dl("fo_csv", "Table (CSV)"), dl("fo_png", "Figure (PNG)")), DT::dataTableOutput(ns("fo_tbl")),
                shiny::helpText("Beaumont and Nichols (1996) / FDIST2 idea for biallelic sites (SNPs): the Fst and the heterozygosity of every site are compared with a neutral island model simulated by the coalescent, with the same sample sizes, and with the migration rate chosen so that the mean simulated Fst equals the observed one. ",
                                "P_high small = more differentiated than neutral (divergent selection); P_low small = less differentiated (balancing selection). This is for many independent SNPs: sites of one gene are linked, and the few sites then give only an exploratory picture.")))),
          shiny::tabPanel("Result", shiny::br(),
            shiny::div(style = "display:flex;gap:6px;flex-wrap:wrap;margin-bottom:10px", dl("rep_html", "Report (HTML, one file, Arlequin-like)"), dl("rep_tbl", "Report (HTML tables, with Structure)"), dl("rep_txt", "Report (text, DnaSP-like)"), dl("dl_sum", "Summary (CSV)")),
            shiny::div(style = "border:1px solid #dfe3e8;border-radius:8px;padding:8px 10px;margin-bottom:10px;background:#fafbfc",
              shiny::strong("Save the results to a folder on this computer"),
              shiny::div(style = "display:flex;gap:8px;align-items:end;flex-wrap:wrap;margin-top:4px",
                shiny::textInput(ns("rep_name"), "Name", "sangerly", width = "170px"),
                shiny::textInput(ns("rep_dir"), "Output folder", .pg_default_dir(), width = "380px"),
                if (.Platform$OS.type == "windows") shiny::actionButton(ns("rep_browse"), "Browse...", style = "margin-bottom:15px"),
                shiny::actionButton(ns("rep_all"), "Save all results", icon = shiny::icon("folder-open"), class = "btn-primary", style = "margin-bottom:15px"),
                shiny::actionButton(ns("rep_folder"), "Only the Arlequin-style folder", style = "margin-bottom:15px")),
              shiny::uiOutput(ns("rep_saved")),
              shiny::helpText("'Save all results' writes the tables (CSV), the HTML report, the Arlequin-style folder, the Structure files and the network (SVG + PNG, when it is built) into this folder, then opens the folder. The download buttons (above and in the other tabs) use the browser instead: the files go to its download folder (usually Downloads).")),
            shiny::h5("Molecular diversity and neutrality tests"), DT::dataTableOutput(ns("sum_tbl")),
            shiny::helpText("S segregating sites; h haplotypes; Hd haplotype diversity (SD); k mean pairwise differences; Pi per site; theta: Hom (homozygosity), S (Watterson), k (alleles), Pi. ",
                            "Tajima's D: p from the beta distribution. Fu's Fs: p = fraction of coalescent simulations with Fs <= observed (significant when p < 0.02, Fu 1997). ",
                            "Which sites are counted (and what pi is divided by) depends on the choice in section 3 ('Sites used for diversity and neutrality'). Theta(Hom) as Arlequin (Zouros 1979; Chakraborty & Weiss 1991)."),
            shiny::h5(style = "margin-top:14px", "Ewens-Watterson neutrality test (infinite alleles model)"), DT::dataTableOutput(ns("ew_tbl")),
            shiny::h5(style = "margin-top:14px", "Strobeck's S and Chakraborty's test of population amalgamation"), DT::dataTableOutput(ns("ewens_tbl")),
            shiny::helpText("Both use Ewens' (1972) sampling formula. Strobeck's S (DnaSP; Strobeck 1987, Fu 1997): probability of a sample with k or fewer haplotypes (and exactly k) given theta = Theta(Pi) (the same theta as Fu's Fs). ",
                            "Chakraborty's test (Arlequin; Chakraborty 1990): probability of k or more haplotypes given theta = Theta(Hom); a small value means more haplotypes than expected, as in a sample that mixes populations. Not computed with a single haplotype or when theta is not defined."),
            shiny::helpText("F = homozygosity of the haplotype frequencies; the null distribution is simulated given the number of haplotypes. ",
                            "p(F_sim <= F_obs) small = haplotypes too evenly distributed (balancing selection); p(F_sim >= F_obs) small = too uneven (directional selection, expansion). Significant at the 5% level when p < 0.025."),
            shiny::h5(style = "margin-top:14px", "Fu and Li's D* and F*, Ramos-Onsins and Rozas' R2"), DT::dataTableOutput(ns("fuli_tbl")),
            shiny::helpText("Only sites with two nucleotides are used (as DnaSP 6); 'singletons' = sites where one sequence differs from all the others. D* and F* (Fu & Li 1993, no outgroup) are large and negative with an excess of singletons ",
                            "(population growth, purifying selection, a selective sweep) and positive with a deficit. R2 (Ramos-Onsins & Rozas 2002) is small after growth, and the most powerful of these tests for small samples. ",
                            "p(<=) and p(>=): fraction of neutral coalescent samples with the same n and the same S that give a value <= / >= the observed one; for D* and F* a result is significant at 5% when p(<=) < 0.05 (growth) or p(>=) < 0.05 (balancing selection, structure); ",
                            "for R2 only p(<=) is used. The number of samples is the one set for Fu's Fs. Not defined for n < 4. The values of D* and F* can differ slightly from DnaSP, which standardises them with the general variance of Achaz (2009); the p-values are not affected by this."),
            shiny::h5(style = "margin-top:14px", "Genetic differentiation among populations and gene flow"), shiny::verbatimTextOutput(ns("gf_text")),
            shiny::h6("Exact test of population differentiation (haplotype frequencies; Raymond & Rousset 1995)"), shiny::verbatimTextOutput(ns("exact_text")), DT::dataTableOutput(ns("exact_pw")),
            DT::dataTableOutput(ns("gf_tests")), shiny::h6("Pairwise"), DT::dataTableOutput(ns("gf_pair")),
            shiny::helpText("Hs, Hst (Hudson, Boos & Kaplan 1992; weights n - 2); Ks, Kst (weights n); Fst (Hudson, Slatkin & Maddison 1992); Nst (Jukes-Cantor, Lynch & Crease 1990); ",
                            "Snn (Hudson 2000); chi-square of the haplotype table. Not yet included: Ks*, Kst*, Z, Z*, Nei's Gst/DeltaSt/GammaSt (formulas to be checked against DnaSP)."),
            shiny::h5(style = "margin-top:14px", "Mismatch distribution"), shiny::uiOutput(ns("mm_unit_ui")), shiny::plotOutput(ns("mm_plot"), height = "360px")))))
  )
}

# ---------------------------------------------------------------------------------------------------
popgenServer <- function(id, aln_all = shiny::reactive(NULL), aln_trim = shiny::reactive(NULL)) {
  shiny::moduleServer(id, function(input, output, session) {
    rv <- shiny::reactiveValues(file_aln = NULL, sets = NULL, pop_how = "names", st_div = NULL, st_amova = NULL, st_struct = NULL, last_out = NULL, net_built = FALSE, pop = NULL, pop_ver = 0L, div = NULL, amova = NULL, struct = NULL, msg = NULL, slide = NULL)
    set_msg <- function(txt, type = "info") rv$msg <- if (is.null(txt)) NULL else list(txt = txt, type = type)
    # messages are shown next to the part of the page where the action happened (aln, pop, set, run, net, struct); warnings also pop up
    mg <- shiny::reactiveValues()
    set_msg_at <- function(area, txt, type = "info") {
      mg[[area]] <- if (is.null(txt)) NULL else list(txt = txt, type = type)
      if (!is.null(txt) && type %in% c("warning", "danger")) shiny::showNotification(txt, type = if (type == "danger") "error" else "warning", duration = 8)
    }
    for (a in c("aln", "pop", "set", "run", "net", "struct", "slide", "map", "more")) local({
      a2 <- a
      output[[paste0("msg_", a2)]] <- shiny::renderUI({ m <- mg[[a2]]; if (is.null(m)) NULL else shiny::div(class = paste0("alert alert-", m$type), style = "padding:6px 10px;margin:6px 0;font-size:13px;word-break:break-word", m$txt) })
    })
    output$msg <- shiny::renderUI({ if (is.null(rv$msg)) return(NULL); shiny::div(class = paste0("alert alert-", rv$msg$type), rv$msg$txt) })
    nz <- function(x, d) if (is.null(x) || length(x) == 0 || (length(x) == 1 && is.na(x))) d else x

    # ---- alignment ---------------------------------------------------------------------------------
    load_aln_file <- function(p) {
      recs <- tryCatch(blast_parse_records(paste(readLines(p, warn = FALSE), collapse = "\n")), error = function(e) NULL)
      if (is.null(recs) || nrow(recs) < 3) { set_msg_at("aln", "The file needs at least 3 aligned sequences.", "warning"); rv$file_aln <- NULL; return() }
      s <- stats::setNames(gsub("[.?]", "-", recs$seq), recs$id)
      if (length(unique(nchar(s))) != 1) { set_msg_at("aln", "The sequences in this file are not all the same length, so it is not an alignment.", "danger"); rv$file_aln <- NULL; return() }
      rv$file_aln <- s; set_msg_at("aln", NULL)
    }
    shiny::observeEvent(input$fa, load_aln_file(input$fa$datapath))
    shiny::observeEvent(input$qt_fasta, load_aln_file(as.character(input$qt_fasta)))

    prep <- shiny::reactive(.pg_prep_names(switch(nz(input$src, "trim"), trim = aln_trim(), all = aln_all(), file = rv$file_aln)))
    clean <- shiny::reactive({ p <- prep(); shiny::req(p); pg_prepare(p$aln, "excluded") })
    output$aln_info <- shiny::renderUI({
      p <- prep()
      if (is.null(p)) return(shiny::helpText("No alignment yet (at least 3 sequences). Align in the MSA tab, or open an aligned FASTA."))
      cl <- pg_prepare(p$aln, "excluded")
      shiny::tagList(
        shiny::helpText(sprintf("%d sequences, %d columns; %d sites without any gap, N or IUPAC code.", length(p$aln), cl$total, cl$kept)),
        if (cl$kept == 0) shiny::div(class = "alert alert-danger", "No site is free of gaps / N / IUPAC codes in every sequence. Trim the alignment further first."),
        if (cl$kept > 0 && cl$kept < 0.5 * cl$total) shiny::div(class = "alert alert-warning", "More than half of the columns are left out. Trim the ragged ends in the MSA tab to keep more sites."),
        if (length(p$dropped)) shiny::div(class = "alert alert-warning", "Left out (no bases at all): ", paste(p$dropped, collapse = ", ")))
    })

    # ---- populations and groups (sequence sets) ---------------------------------------------------------------
    make_pop <- function(nm) pg_pop_from_names(nm, nz(input$pop_method, "prefix"), nz(input$pop_sep, "_"), nz(input$pop_nfield, 1L), nz(input$pop_nchars, 3L), nz(input$pop_regex, ""))
    grp_ver <- shiny::reactiveVal(0L)
    bump_groups <- function() grp_ver(grp_ver() + 1L)
    bump <- function() { rv$pop_ver <- rv$pop_ver + 1L; bump_groups() }
    # the table loaded from CSV / Excel is remembered, so applying the rule again or re-reading the alignment never wipes it
    apply_sets <- function(d0, d, cols = NULL) {
      if (is.null(cols)) cols <- pg_pick_set_cols(names(d))
      kcol <- cols$key; pcol <- cols$pop; gcol <- if (is.null(cols$grp) || is.na(cols$grp) || !cols$grp %in% names(d)) NULL else cols$grp
      d <- d[!is.na(d[[kcol]]) & nzchar(trimws(d[[kcol]])), , drop = FALSE]
      mt <- pg_match_keys(d0$Sequence, d[[kcol]]); ok <- !is.na(mt$index)
      pv <- trimws(as.character(d[[pcol]][mt$index[ok]])); good <- !is.na(pv) & nzchar(pv)
      okp <- which(ok)[good]
      d0$Population[okp] <- pv[good]
      if (!is.null(gcol)) { gv <- trimws(as.character(d[[gcol]][mt$index[ok]])); gg <- !is.na(gv) & nzchar(gv); d0$Group[which(ok)[gg]] <- gv[gg] }
      list(d0 = d0, ok = ok, mt = mt, kcol = kcol, empty = sum(ok) - sum(good))
    }
    regen <- function() {
      p <- prep()
      if (is.null(p)) { rv$pop <- NULL; bump(); return() }
      src <- nz(input$pop_src, "names")
      d0 <- data.frame(Sequence = names(p$aln), Population = make_pop(names(p$aln)), Group = "all", stringsAsFactors = FALSE); how <- "names"
      if (identical(src, "file") && !is.null(rv$sets)) { r <- apply_sets(d0, rv$sets$table, rv$sets$cols); if (any(r$ok)) { d0 <- r$d0; how <- "file"; rv$sets$matched <- sum(r$ok); rv$sets$n <- nrow(d0) } }
      if (identical(src, "manual") && !is.null(shiny::isolate(rv$pop))) {      # keep what was assigned by hand
        old <- shiny::isolate(rv$pop); i <- match(d0$Sequence, old$Sequence); k <- !is.na(i)
        d0$Population[k] <- old$Population[i[k]]; d0$Group[k] <- old$Group[i[k]]; how <- "manual"
      }
      rv$pop_how <- how
      if (!identical(shiny::isolate(rv$pop), d0)) { rv$pop <- d0; bump() }        # nothing changed = the results stay valid
    }
    to_manual <- function() { if (!identical(nz(input$pop_src, "names"), "manual")) shiny::updateRadioButtons(session, "pop_src", selected = "manual"); rv$pop_how <- "manual" }
    shiny::observeEvent(input$pop_src, { if (identical(input$pop_src, "manual")) rv$pop_how <- "manual" else regen() }, ignoreInit = TRUE)
    output$pop_srclabel <- shiny::renderUI({
      d <- rv$pop; if (is.null(d)) return(NULL); how <- rv$pop_how; src <- nz(input$pop_src, "names")
      rule <- switch(nz(input$pop_method, "prefix"), prefix = sprintf("text before '%s' (field %d)", nz(input$pop_sep, "_"), as.integer(nz(input$pop_nfield, 1L))),
                     chars = sprintf("first %d characters", as.integer(nz(input$pop_nchars, 3L))), regex = "a regular expression", all = "one population for all")
      if (identical(src, "file") && !identical(how, "file")) return(shiny::div(class = "alert alert-warning", style = "padding:6px 10px;font-size:13px", "No file loaded yet: the populations shown now come from the sequence names. Open a CSV / Excel file above."))
      txt <- switch(how, names = paste0("sequence names (", rule, ")"), file = sprintf("%s (%d of %d sequences matched)", rv$sets$name, as.integer(rv$sets$matched), as.integer(rv$sets$n)), manual = "sets defined by hand")
      shiny::div(style = "margin:4px 0 6px;font-size:13px", shiny::strong("Source in use: "), txt, sprintf(" - %d populations, %d sequences.", length(unique(d$Population)), nrow(d)))
    })
    output$pop_tpl <- shiny::downloadHandler("sequence_sets_template.csv", function(file) {
      p <- shiny::req(prep()); nm <- names(p$aln)
      utils::write.csv(data.frame(Sequence = nm, Population = make_pop(nm), Group = "all", stringsAsFactors = FALSE), file, row.names = FALSE)
    })
    shiny::observeEvent(input$load_example, {
      aln <- pg_example(); tb <- pg_example_sets(); cols <- pg_pick_set_cols(names(tb))
      rv$file_aln <- aln; rv$sets <- list(table = tb, name = "Durio_zibethinus_ID (example table)", cols = cols, matched = nrow(tb), n = length(aln))
      shiny::updateRadioButtons(session, "src", selected = "file"); shiny::updateRadioButtons(session, "pop_src", selected = "file")
      set_msg_at("aln", "Example loaded: 60 ITS sequences of Durio zibethinus (654 sites) and the table of their locations (columns ID, Species, Location). The populations are the locations. Now go to section 4 and press the first button.", "success")
    })
    output$csv_state <- shiny::renderUI({
      s <- rv$sets; if (is.null(s)) return(NULL)
      nmx <- names(s$table); sel <- function(x, none = FALSE) if (is.null(x) || is.na(x)) (if (none) "(none)" else nmx[1]) else x
      shiny::div(style = "margin-bottom:6px",
        shiny::tags$span(class = "text-success", sprintf("Table in use: %s (%d rows).", s$name, nrow(s$table))), " ", shiny::actionLink(session$ns("clear_csv"), "Forget the table"),
        shiny::div(style = "display:flex;gap:8px;flex-wrap:wrap;margin-top:4px",
          shiny::selectInput(session$ns("col_key"), "Sequence name in column", nmx, sel(s$cols$key), width = "170px"),
          shiny::selectInput(session$ns("col_pop"), "Population in column", nmx, sel(s$cols$pop), width = "170px"),
          shiny::selectInput(session$ns("col_grp"), "Group in column", c("(none)", nmx), sel(s$cols$grp, TRUE), width = "170px")))
    })
    shiny::observeEvent(list(input$col_key, input$col_pop, input$col_grp), {
      s <- rv$sets; if (is.null(s) || is.null(input$col_key) || is.null(input$col_pop)) return()
      cols <- list(key = input$col_key, pop = input$col_pop, grp = if (is.null(input$col_grp) || identical(input$col_grp, "(none)")) NA_character_ else input$col_grp)
      if (identical(cols, s$cols)) return()
      rv$sets$cols <- cols; regen()
      m <- rv$sets$matched; if (!is.null(m)) set_msg_at("pop", sprintf("Columns changed: %s = Sequence, %s = Population%s. %d of %d sequences matched.", cols$key, cols$pop, if (is.na(cols$grp)) "" else paste0(", ", cols$grp, " = Group"), as.integer(m), as.integer(rv$sets$n)), "info")
    }, ignoreInit = TRUE)
    shiny::observeEvent(input$clear_csv, { rv$sets <- NULL; if (identical(nz(input$pop_src, "names"), "file")) shiny::updateRadioButtons(session, "pop_src", selected = "names"); regen(); set_msg_at("pop", "The loaded table was removed; populations now come from the rule.", "info") })
    shiny::observeEvent(names(prep()$aln), regen(), ignoreNULL = FALSE)
    shiny::observeEvent(input$apply_rule, regen())
    load_pop_csv <- function(f) {
      d0 <- rv$pop
      if (is.null(d0)) { set_msg_at("pop", "Load the alignment first (step 1); then load the table of sequence sets.", "warning"); return() }
      d <- pg_read_sets(f$datapath, f$name)
      if (is.character(d)) { set_msg_at("pop", d, "warning"); return() }
      if (ncol(d) < 2) { set_msg_at("pop", "The table needs at least two columns: Sequence and Population (and optionally Group).", "warning"); return() }
      cols <- pg_pick_set_cols(names(d)); r <- apply_sets(d0, d, cols); mt <- r$mt; ok <- r$ok
      if (!any(ok)) {
        set_msg_at("pop", sprintf("None of the %d sequences matched the table. First keys in column '%s': %s. First sequence names: %s.", nrow(d0), r$kcol,
                        paste(utils::head(d[[r$kcol]], 3), collapse = ", "), paste(utils::head(d0$Sequence, 3), collapse = ", ")), "warning"); return()
      }
      rv$sets <- list(table = d, name = f$name, cols = cols, matched = sum(ok), n = nrow(d0))
      d0 <- r$d0; rv$pop <- d0; rv$pop_how <- "file"; bump()
      tb <- table(factor(mt$how[ok], levels = c("exact", "prefix", "accession")))
      cn <- local({ u <- unique(d0$Population); g <- split(u, tolower(u)); g <- g[vapply(g, length, 1L) > 1]
                    if (length(g)) paste0(" Populations that differ only in capitals (rename one in the Populations table to merge them): ", paste(vapply(g, paste, "", collapse = " / "), collapse = "; "), ".") else "" })
      set_msg_at("pop", paste0(sprintf("%d of %d sequences matched the table (exact: %d, name starts with the key: %d, accession number: %d). Columns used: %s = Sequence, %s = Population%s - change them below if they are not the right ones.%s",
                      sum(ok), nrow(d0), tb[["exact"]], tb[["prefix"]], tb[["accession"]], cols$key, cols$pop, if (is.na(cols$grp)) "" else paste0(", ", cols$grp, " = Group"),
                      if (any(!ok)) paste0(" Not found in the CSV (keep their previous population): ", paste(utils::head(d0$Sequence[!ok], 5), collapse = ", "), if (sum(!ok) > 5) ", ..." else "") else ""), cn),
              if (any(!ok) || nzchar(cn)) "warning" else "success")
    }
    shiny::observeEvent(input$pop_csv, load_pop_csv(input$pop_csv))
    # Qt dialog: paths chosen on this computer (the file boxes of the page cannot be filled from outside)
    pending_csv <- shiny::reactiveVal(NULL)   # the table may arrive before the alignment is processed: wait for it
    shiny::observeEvent(input$qt_popfile, {
      f <- list(datapath = as.character(input$qt_popfile), name = basename(as.character(input$qt_popfile)))
      if (is.null(rv$pop)) pending_csv(f) else load_pop_csv(f)
    })
    shiny::observeEvent(rv$pop, { f <- pending_csv(); if (!is.null(f) && !is.null(rv$pop)) { pending_csv(NULL); load_pop_csv(f) } }, ignoreNULL = TRUE)

    output$seq_tbl <- DT::renderDataTable({
      rv$pop_ver; d <- shiny::isolate(rv$pop); shiny::req(d)
      dtt(d, rownames = FALSE, selection = "multiple", editable = list(target = "cell", disable = list(columns = 0)),
                    options = list(pageLength = 10, lengthMenu = LENGTHS, dom = "lftp", scrollX = TRUE))
    })
    seq_proxy <- DT::dataTableProxy("seq_tbl", session = session)
    shiny::observeEvent(input$seq_tbl_cell_edit, { d <- rv$pop; shiny::req(d); rv$pop <- DT::editData(d, input$seq_tbl_cell_edit, rownames = FALSE); to_manual(); bump_groups() })
    shiny::observeEvent(input$sel_all, DT::selectRows(seq_proxy, input$seq_tbl_rows_all))
    shiny::observeEvent(input$sel_none, DT::selectRows(seq_proxy, NULL))
    shiny::observe({                                         # keep the drop-downs in step with the table
      d <- rv$pop; shiny::req(d)
      shiny::updateSelectizeInput(session, "set_pop", choices = sort(unique(d$Population)), selected = shiny::isolate(input$set_pop))
      shiny::updateSelectizeInput(session, "set_grp", choices = sort(unique(d$Group)), selected = shiny::isolate(input$set_grp))
    })
    canon <- function(nm, existing) { hit <- existing[tolower(existing) == tolower(nm)]; if (length(hit)) hit[1] else nm }
    shiny::observeEvent(input$assign_pop, {
      s <- input$seq_tbl_rows_selected; nm <- trimws(nz(input$set_pop, ""))
      if (!length(s) || !nzchar(nm)) { set_msg_at("pop", "Select sequences in the table first, then pick (or type) the population.", "warning"); return() }
      d <- rv$pop; nm2 <- canon(nm, unique(d$Population))
      if (!identical(nm2, nm)) shiny::showNotification(sprintf("'%s' is the same as the existing population '%s' (capitals ignored): '%s' is used.", nm, nm2, nm2), type = "message")
      d$Population[s] <- nm2; rv$pop <- d; to_manual(); bump()
    })
    shiny::observeEvent(input$assign_grp, {
      s <- input$seq_tbl_rows_selected; nm <- trimws(nz(input$set_grp, ""))
      if (!length(s) || !nzchar(nm)) { set_msg_at("pop", "Select sequences in the table first, then pick (or type) the group.", "warning"); return() }
      d <- rv$pop; nm2 <- canon(nm, unique(d$Group))
      if (!identical(nm2, nm)) shiny::showNotification(sprintf("'%s' is the same as the existing group '%s' (capitals ignored): '%s' is used.", nm, nm2, nm2), type = "message")
      d$Group[s] <- nm2; rv$pop <- d; to_manual(); bump()
    })
    grp_df <- shiny::reactive({
      d <- rv$pop; shiny::req(d)
      u <- unique(d$Population)
      data.frame(Population = u, n = as.integer(table(factor(d$Population, levels = u))),
                 Group = vapply(u, function(p) d$Group[d$Population == p][1], ""), stringsAsFactors = FALSE, row.names = NULL)
    })
    output$grp_tbl <- DT::renderDataTable({
      grp_ver(); g <- shiny::isolate(grp_df())
      dtt(g, rownames = FALSE, editable = list(target = "cell", disable = list(columns = 1)),
                    options = list(pageLength = 10, lengthMenu = LENGTHS, dom = "lftp", scrollX = TRUE))
    })
    shiny::observeEvent(input$grp_tbl_cell_edit, {
      info <- input$grp_tbl_cell_edit; g <- grp_df(); d <- rv$pop
      p <- g$Population[info$row]; val <- trimws(as.character(info$value))
      if (!nzchar(val)) return()
      if (info$col == 0) d$Population[d$Population == p] <- val else d$Group[d$Population == p] <- val
      rv$pop <- d; to_manual()
      if (info$col == 0) bump()
    })
    pop_for <- function(M) {
      d <- rv$pop; i <- match(rownames(M), d$Sequence)
      pop <- trimws(d$Population[i]); grp <- trimws(d$Group[i])
      pop[is.na(pop) | !nzchar(pop)] <- "unassigned"; grp[is.na(grp) | !nzchar(grp)] <- "all"
      list(pop = pop, grp = grp)
    }
    output$gf_pops_ui <- shiny::renderUI({
      d <- rv$pop; shiny::req(d); u <- unique(d$Population)
      shiny::checkboxGroupInput(session$ns("gf_pops"), "Populations included (differentiation)", u, selected = u)
    })

    # ---- haplotypes, diversity, neutrality, differentiation --------------------------------------------------------
    shiny::observeEvent(input$run_div, {
      cl <- tryCatch(clean(), error = function(e) NULL)
      if (is.null(cl)) { set_msg_at("run", "Provide an alignment with at least 3 sequences first.", "warning"); return() }
      M <- cl$M
      if (ncol(M) == 0) { set_msg_at("run", "No site is free of gaps / N / IUPAC codes in every sequence. Trim the alignment further first.", "warning"); return() }
      pop <- pop_for(M)$pop
      pa <- prep()$aln; dmode <- nz(input$div_gap, "all"); dmiss <- nz(input$miss, 0)
      # the matrix of a set of sequences (rows ii) under the chosen rule for gaps / missing data
      Xu <- function(ii) {
        X <- switch(dmode, all = M[ii, , drop = FALSE], subset = pg_prepare(pa[ii], "level", dmiss)$M, pairwise = pg_prepare(pa[ii], "pairwise")$M)
        if (ncol(X) == 0) X <- M[ii, , drop = FALSE]
        X
      }
      Lden <- if (isTRUE(input$div_inv)) { Mw <- pg_prepare(pa, "pairwise")$M
        sum(apply(Mw, 2, function(x) length(unique(x[x %in% c("A", "C", "G", "T")])) > 1)) } else NULL
      prog <- .pg_progress(session, "Population genetics"); on.exit(prog$close(), add = TRUE)
      units <- c(list(All = seq_len(nrow(M))), split(seq_len(nrow(M)), pop))
      nsim <- nz(input$nsim, 1000); rows <- list(); ewrows <- list(); tjl <- list(); fuli <- list(); ewl <- list()
      for (u in seq_along(units)) {
        ii <- units[[u]]; Xi <- Xu(ii); st <- pg_stats(Xi, Lden)
        prog$set(70 * (u - 1) / length(units), paste("statistics:", names(units)[u]))
        pfs <- if (nsim > 0 && st$n >= 3 && is.finite(st$Fs))
          pg_fs_pvalue(st$n, st$k, st$Fs, nsim, progress = function(f) prog$set(70 * ((u - 1) + f) / length(units), sprintf("Fu's Fs simulations: %s", names(units)[u]))) else NA_real_
        hn <- as.integer(table(pg_haplotypes(Xi)$hap)); hn <- hn[hn > 0]
        tt <- if (st$n >= 2) pg_thetas(st$n, st$S, hn, st$k) else c(theta_Hom = NA, theta_S = NA, theta_k = NA, theta_Pi = NA)
        ewl[[u]] <- data.frame(Population = names(units)[u], t(pg_ewens_tests(st$n, st$Nh, tt[["theta_Pi"]], tt[["theta_Hom"]])), check.names = FALSE, stringsAsFactors = FALSE)
        ew <- pg_ew_test(hn, as.integer(nz(input$ew_nsim, 1000)))
        ewrows[[u]] <- data.frame(Population = names(units)[u], n = ew$n, h = ew$k, F_obs = ew$F, F_mean_simulated = ew$mean, p_F_sim_le_F_obs = ew$p_lower,
                                  p_F_sim_ge_F_obs = ew$p_upper, simulated = ew$nsim, stringsAsFactors = FALSE)
        sdv <- pg_arl_sd(st$n, st$k)
        tjs <- if (nsim > 0 && st$n >= 4 && st$S >= 1) pg_tajima_sim(st$n, st$S, st$k, nsim) else .pg_tj_na
        tjl[[names(units)[u]]] <- tjs
        fl <- tryCatch(pg_fuli(Xi, nsim, progress = function(f) prog$set(70 * ((u - 1) + f) / length(units), sprintf("Fu and Li's D*, F*, R2: %s", names(units)[u]))), error = function(e) .pg_fuli_na)
        fuli[[u]] <- data.frame(Population = names(units)[u], t(fl), check.names = FALSE, stringsAsFactors = FALSE)
        rows[[u]] <- data.frame(Population = names(units)[u], n = st$n, Sites = st$L, S = st$S, Eta = st$eta, h = st$Nh, Hd = st$Hd, SD_Hd = st$SD_Hd,
                                k = st$k, Pi = st$pi, theta_Hom = tt[["theta_Hom"]], theta_S = tt[["theta_S"]], theta_k = tt[["theta_k"]], theta_Pi = tt[["theta_Pi"]],
                                Tajima_D = st$D, p_D = st$p_D, Fs = st$Fs, p_Fs = pfs,
                                SD_k = sdv[["SD_k"]], SD_Pi = sdv[["SD_k"]] / (if (is.null(Lden)) st$L else Lden), SD_theta_Pi = sdv[["SD_theta_Pi"]], p_D_sim = tjs[["p"]], stringsAsFactors = FALSE)
      }
      hp <- pg_haplotypes(M)
      cnt <- table(hap = hp$hap, pop = factor(pop))
      counts <- matrix(as.integer(cnt), nrow = nrow(cnt), dimnames = list(rownames(cnt), colnames(cnt)))
      members <- tapply(rownames(M), hp$hap, function(x) paste(x, collapse = ", "))
      haptab <- cbind(data.frame(Haplotype = hp$id, Count = hp$n, Frequency = round(hp$n / nrow(M), 4), stringsAsFactors = FALSE),
                      as.data.frame.matrix(counts), Sequences = as.character(members[hp$id]))
      # genetic differentiation (DnaSP)
      gf <- NULL; keep <- nz(input$gf_pops, unique(pop)); gap <- nz(input$gf_gap, "excluded")
      prog$set(72, "genetic differentiation")
      Mg <- pg_prepare(pa, gap)$M; pg_pop <- pop_for(Mg)$pop; inc <- pg_pop %in% keep
      if (length(unique(pg_pop[inc])) >= 2 && sum(inc) >= 3 && ncol(Mg) > 0) {
        gf <- tryCatch(pg_geneflow(Mg[inc, , drop = FALSE], pg_pop[inc], gap, nz(input$gf_genome, "diploid"),
                                   as.integer(nz(input$gf_nperm, 0)), if (is.na(nz(input$gf_seed, NA))) NULL else input$gf_seed,
                                   progress = function(f) prog$set(72 + 28 * f, "permutations")),
                       error = function(e) { set_msg_at("run", paste("Genetic differentiation:", conditionMessage(e)), "warning"); NULL })
      } else set_msg_at("run", "Genetic differentiation needs at least 2 populations (with 3 sequences in total) among the included populations.", "info")
      ex <- NULL; exsteps <- as.integer(nz(input$ex_steps, 0)); exdem <- as.integer(nz(input$ex_demem, 10000)); exseed <- if (is.na(nz(input$gf_seed, NA))) NULL else input$gf_seed
      if (exsteps > 0 && ncol(counts) >= 2) {
        prog$set(90, "exact test of differentiation")
        ex <- pg_exact_test(counts, exsteps, exdem, exseed)
        if (isTRUE(input$ex_pair) && ncol(counts) >= 3)
          ex$pairwise <- pg_exact_pairwise(counts, max(2000L, exsteps %/% 5L), max(500L, exdem %/% 5L), exseed, progress = function(f) prog$set(90 + 9 * f, "exact test, pairs"))
      }
      rv$div <- list(M = M, cols = cl$cols, ew = do.call(rbind, ewrows), exact = ex, total = cl$total, stats = do.call(rbind, rows), hap = hp, counts = counts, haptab = haptab, gf = gf,
                     mm = lapply(units, function(ii) pg_mismatch(M[ii, , drop = FALSE])), pop = pop, aln = pa,
                     det = pg_arl_details(M, pop), tj = tjl, nsim_fs = nsim, fuli = do.call(rbind, fuli), ewens = do.call(rbind, ewl),
                     divmode = list(mode = dmode, miss = dmiss, Lden = Lden))
      prog$set(100, "done"); if (!is.null(gf)) set_msg_at("run", NULL)
    })
    fmt <- function(d, cols, digits) { for (cc in intersect(cols, names(d))) d[[cc]] <- round(d[[cc]], digits); d }
    output$sum_tbl <- DT::renderDataTable({
      d <- rv$div; shiny::req(d)
      s <- fmt(fmt(d$stats, c("Hd", "SD_Hd", "k", "Tajima_D", "Fs", "theta_Hom", "theta_S", "theta_k", "theta_Pi"), 3), c("Pi", "SD_Pi"), 5)
      dtt(fmt(fmt(s, c("p_D", "p_Fs", "p_D_sim"), 4), c("SD_k", "SD_theta_Pi"), 3), rownames = FALSE, options = list(pageLength = 15, dom = "tp", scrollX = TRUE))
    })
    output$ew_tbl <- DT::renderDataTable({
      d <- rv$div; shiny::req(d$ew)
      dtt(fmt(d$ew, c("F_obs", "F_mean_simulated", "p_F_sim_le_F_obs", "p_F_sim_ge_F_obs"), 4), rownames = FALSE, options = list(dom = "t", scrollX = TRUE))
    })
    output$ewens_tbl <- DT::renderDataTable({
      d <- rv$div; shiny::req(d$ewens); e <- d$ewens
      names(e) <- c("Population", "n genes", "k haplotypes", "Theta(Pi)", "Strobeck S: P(K <= k)", "P(K = k)", "Theta(Hom)", "Chakraborty: P(K >= k)")
      dtt(fmt(e, setdiff(names(e), c("Population", "n genes", "k haplotypes")), 4), rownames = FALSE, options = list(dom = "t", scrollX = TRUE))
    })
    output$fuli_tbl <- DT::renderDataTable({
      d <- rv$div; shiny::req(d$fuli)
      f <- d$fuli; names(f) <- sub("^D_star", "D*", sub("^F_star", "F*", gsub("Y_star", "Y*", gsub("D_star", "D*", gsub("F_star", "F*", names(f))))))
      f$n <- as.integer(f$n)
      dtt(fmt(f, setdiff(names(f), c("Population", "n", "S_biallelic", "singletons", "nsim")), 4), rownames = FALSE, options = list(dom = "t", scrollX = TRUE))
    })
    # ---- more tests: SFS, LD, Mantel, Hardy-Weinberg, EM phasing, Fst outliers ----------------------------------------------------------
    rx <- shiny::reactiveValues(ld = NULL, sfs = NULL, mn = NULL, hw = NULL, em = NULL, fo = NULL)
    output$mt_unit_ui <- shiny::renderUI({ d <- rv$div; shiny::selectInput(session$ns("mt_unit"), "Sequences", c("All", if (!is.null(d)) unique(d$pop)), width = "200px") })
    shiny::observeEvent(rv$div, { rx$ld <- NULL; rx$sfs <- NULL; rx$fo <- NULL }, ignoreNULL = FALSE)
    shiny::observeEvent(input$run_ld, {
      d <- rv$div; if (is.null(d)) { set_msg_at("more", "Run 'Haplotypes, diversity, differentiation' first.", "warning"); return() }
      u <- nz(input$mt_unit, "All"); ii <- if (identical(u, "All")) seq_len(nrow(d$M)) else which(d$pop == u)
      X <- d$M[ii, , drop = FALSE]
      if (length(ii) < 4) { set_msg_at("more", "At least 4 sequences are needed.", "warning"); return() }
      rx$sfs <- pg_sfs(X)
      rx$ld <- tryCatch(pg_ld(X, d$cols, nz(input$ld_max, 300), isTRUE(input$ld_fisher)), error = function(e) { set_msg_at("more", paste("LD:", conditionMessage(e)), "danger"); NULL })
      if (is.null(rx$ld)) set_msg_at("more", "Fewer than 2 biallelic sites: no linkage disequilibrium to compute.", "info") else set_msg_at("more", NULL)
    })
    output$sfs_tbl <- DT::renderDataTable({ s <- shiny::req(rx$sfs); dtt(fmt(s, c("Expected_neutral", "Observed_proportion", "Expected_proportion"), 4), rownames = FALSE, options = list(dom = "tp", pageLength = 8, scrollX = TRUE)) })
    output$sfs_plot <- shiny::renderPlot({
      s <- shiny::req(rx$sfs); o <- graphics::par(mar = c(4.2, 4.2, 1, 1))
      b <- graphics::barplot(s$Observed_sites, names.arg = s$Copies_of_minor_variant, col = "#5b8def", xlab = "Copies of the minor variant (folded spectrum)", ylab = "Number of sites")
      graphics::lines(b, s$Expected_neutral, type = "b", pch = 19, col = "#d62728", lwd = 2); graphics::legend("topright", c("observed", "neutral expectation"), fill = c("#5b8def", NA), border = NA, lty = c(NA, 1), pch = c(NA, 19), col = c("#5b8def", "#d62728"), bty = "n")
      graphics::par(o)
    })
    output$ld_sum <- DT::renderDataTable({ l <- shiny::req(rx$ld); dtt(fmt(l$summary, c("ZnS", "Za", "ZZ"), 4), rownames = FALSE, options = list(dom = "t", scrollX = TRUE)) })
    output$ld_tbl <- DT::renderDataTable({ l <- shiny::req(rx$ld); dtt(fmt(l$pairs, c("D", "Dprime", "r2", "P_Fisher", "P_Bonferroni"), 4), rownames = FALSE, options = list(pageLength = 10, dom = "tp", scrollX = TRUE)) })
    ld_draw <- function() {
      l <- shiny::req(rx$ld); p <- l$pairs; o <- graphics::par(mar = c(4.2, 4.2, 1, 1))
      graphics::plot(p$Distance, p$r2, pch = 19, cex = 0.7, col = grDevices::adjustcolor("#2a6fdb", 0.5), xlab = "Distance between the sites (alignment columns)", ylab = expression(r^2), ylim = c(0, 1))
      if (nrow(p) > 10) graphics::lines(stats::lowess(p$Distance, p$r2), col = "#d62728", lwd = 2); graphics::par(o)
    }
    output$ld_plot <- shiny::renderPlot({ ld_draw() })
    output$ld_png <- shiny::downloadHandler("ld_decay.png", function(file) { grDevices::png(file, 1800, 1200, res = 220); on.exit(grDevices::dev.off()); ld_draw() })
    output$ld_csv <- shiny::downloadHandler("ld_pairs.csv", function(file) { l <- shiny::req(rx$ld); utils::write.csv(l$pairs, file, row.names = FALSE) })
    output$sfs_csv <- shiny::downloadHandler("sfs.csv", function(file) { s <- shiny::req(rx$sfs); utils::write.csv(s, file, row.names = FALSE) })

    shiny::observeEvent(input$run_mn, {
      a <- rv$amova; cd <- rv_coords()
      if (is.null(a)) { set_msg_at("more", "Press 'AMOVA and pairwise FST' first (tab AMOVA).", "warning"); return() }
      if (is.null(cd)) { set_msg_at("more", "Load the coordinates of the populations first (tab Map).", "warning"); return() }
      ph <- a$pw$phi; ph[is.na(ph)] <- t(ph)[is.na(ph)]; ph[ph < 0] <- 0; nm <- rownames(ph); m <- match(tolower(nm), tolower(cd$Population)); ok <- !is.na(m)
      if (sum(ok) < 4) { set_msg_at("more", "At least 4 populations with coordinates and pairwise FST are needed.", "warning"); return() }
      G <- ph[ok, ok]; diag(G) <- 0; if (identical(input$mn_gen, "rousset")) G <- G / pmax(1 - G, 1e-9)
      Dg <- pg_haversine(cd$lon[m[ok]], cd$lat[m[ok]]); if (identical(input$mn_geo, "ln")) { Dg <- log(pmax(Dg, 1)) }
      dimnames(G) <- dimnames(Dg) <- list(nm[ok], nm[ok])
      r <- pg_mantel(G, Dg, nz(input$mn_perm, 9999)); if (is.null(r)) { set_msg_at("more", "The Mantel test could not be computed.", "warning"); return() }
      rx$mn <- list(r = r, G = G, Dg = Dg, gen = nz(input$mn_gen, "fst"), geo = nz(input$mn_geo, "km")); set_msg_at("more", NULL)
    })
    output$mn_text <- shiny::renderText({ m <- rx$mn; if (is.null(m)) return("Not run yet."); r <- m$r
      sprintf("Populations: %d (%d pairs)\nMantel r = %.4f\np (one-sided, r > 0; isolation by distance) = %.4f\np (two-sided) = %.4f\nPermutations: %d", r$npop, r$npairs, r$r, r$p_right, r$p_two, r$nperm) })
    mn_draw <- function() {
      m <- shiny::req(rx$mn); lt <- lower.tri(m$G); o <- graphics::par(mar = c(4.5, 4.5, 1, 1))
      graphics::plot(m$Dg[lt], m$G[lt], pch = 19, col = "#2a6fdb", xlab = if (m$geo == "ln") "ln(distance in km)" else "Geographic distance (km)", ylab = if (m$gen == "rousset") "FST / (1 - FST)" else "Pairwise FST")
      if (sum(lt) > 2) graphics::abline(stats::lm(m$G[lt] ~ m$Dg[lt]), col = "#d62728", lwd = 2)
      graphics::legend("topleft", sprintf("r = %.3f, p = %.4f", m$r$r, m$r$p_right), bty = "n"); graphics::par(o)
    }
    output$mn_plot <- shiny::renderPlot({ mn_draw() })
    output$mn_png <- shiny::downloadHandler("mantel.png", function(file) { grDevices::png(file, 1600, 1200, res = 220); on.exit(grDevices::dev.off()); mn_draw() })

    shiny::observeEvent(input$run_hw, {
      pa <- tryCatch(prep()$aln, error = function(e) NULL); if (is.null(pa)) { set_msg_at("more", "Provide an alignment first.", "warning"); return() }
      M <- pg_prepare(pa, "pairwise")$M; pop <- pop_for(M)$pop
      r <- pg_hw_sites(M, pop, minn = nz(input$hw_min, 5))
      if (is.null(r)) { set_msg_at("more", "No site has heterozygous IUPAC codes (R, Y, S, W, K, M) with enough individuals: the sequences look haploid, so there are no genotypes to test.", "info"); rx$hw <- NULL; return() }
      rx$hw <- r; set_msg_at("more", sprintf("%d site x population tests.", nrow(r)), "success")
    })
    output$hw_tbl <- DT::renderDataTable({ r <- shiny::req(rx$hw); dtt(fmt(r, c("Freq_first", "Ho", "He", "Fis"), 4), rownames = FALSE, options = list(pageLength = 10, dom = "lftp", scrollX = TRUE)) })
    output$hw_csv <- shiny::downloadHandler("hardy_weinberg_sites.csv", function(file) { r <- shiny::req(rx$hw); utils::write.csv(r, file, row.names = FALSE) })

    shiny::observeEvent(input$run_em, {
      pa <- tryCatch(prep()$aln, error = function(e) NULL); if (is.null(pa)) { set_msg_at("more", "Provide an alignment first.", "warning"); return() }
      M <- pg_prepare(pa, "pairwise")$M
      if (!any(M %in% names(.pg_iupac2))) { set_msg_at("more", "No heterozygous IUPAC code (R, Y, S, W, K, M) in the alignment: nothing to phase (the sequences look haploid).", "info"); rx$em <- NULL; return() }
      r <- tryCatch(pg_phase_em(M, rownames(M), nz(input$em_max, 12)), error = function(e) { set_msg_at("more", paste("Phasing:", conditionMessage(e)), "danger"); NULL })
      if (is.null(r)) return(); rx$em <- r
      set_msg_at("more", sprintf("%d individuals phased, %d haplotypes, %d EM iterations%s%s.", r$n_ind, nrow(r$freq), r$iterations,
                                 if (length(r$skipped)) sprintf("; %d individuals with too many heterozygous sites were left out", length(r$skipped)) else "",
                                 if (r$dropped_columns > 0) sprintf("; %d columns with gaps / N / other codes were not used", r$dropped_columns) else ""), "success")
    })
    output$em_text <- shiny::renderText({ r <- rx$em; if (is.null(r)) return("Not run yet."); sprintf("Log-likelihood = %.3f after %d iterations; %d individuals; %d haplotypes.", r$loglik, r$iterations, r$n_ind, nrow(r$freq)) })
    output$em_freq <- DT::renderDataTable({ r <- shiny::req(rx$em); dtt(fmt(r$freq, c("Frequency", "Copies"), 4), rownames = FALSE, options = list(pageLength = 8, dom = "tp", scrollX = TRUE)) })
    output$em_ind <- DT::renderDataTable({ r <- shiny::req(rx$em); dtt(fmt(r$ind, "Posterior", 4), rownames = FALSE, options = list(pageLength = 8, dom = "tp", scrollX = TRUE)) })
    output$em_fa <- shiny::downloadHandler("phased_haplotypes.fasta", function(file) { r <- shiny::req(rx$em); writeLines(as.vector(rbind(paste0(">", rownames(r$phased)), apply(r$phased, 1, paste, collapse = ""))), file) })
    output$em_csv <- shiny::downloadHandler("phased_haplotype_frequencies.csv", function(file) { r <- shiny::req(rx$em); utils::write.csv(r$freq, file, row.names = FALSE) })
    output$em_ind_csv <- shiny::downloadHandler("phased_individuals.csv", function(file) { r <- shiny::req(rx$em); utils::write.csv(r$ind, file, row.names = FALSE) })

    shiny::observeEvent(input$run_fo, {
      d <- rv$div; if (is.null(d)) { set_msg_at("more", "Run 'Haplotypes, diversity, differentiation' first.", "warning"); return() }
      prog <- .pg_progress(session, "Fst outliers"); on.exit(prog$close(), add = TRUE)
      sd <- nz(input$fo_seed, NA)
      r <- tryCatch(pg_fst_outlier(d$M, d$pop, d$cols, as.integer(nz(input$fo_sim, 10000)), nz(input$fo_mac, 2), progress = function(f, m = "") prog$set(100 * f, m), seed = if (is.na(sd)) NULL else sd),
                    error = function(e) { set_msg_at("more", paste("Fst outliers:", conditionMessage(e)), "danger"); NULL })
      if (is.null(r)) { if (is.null(mg$more)) set_msg_at("more", "Not enough biallelic sites or populations (at least 2 populations with 2 sequences).", "info"); rx$fo <- NULL; return() }
      rx$fo <- r; set_msg_at("more", NULL)
    })
    output$fo_text <- shiny::renderText({ r <- rx$fo; if (is.null(r)) return("Not run yet."); t <- r$table
      sprintf("%d populations (%s); %d biallelic sites; mean Fst = %.4f; simulated migration rate 4Nm = %.3f; %d simulated loci.\nSites with P_high < 0.05: %d; with P_low < 0.05: %d (expected by chance: about %d each).", r$npop, paste(r$pops, collapse = ", "), nrow(t), r$mean_fst, r$migration, r$nsim, sum(t$P_high < 0.05), sum(t$P_low < 0.05), round(0.05 * nrow(t))) })
    fo_draw <- function() {
      r <- shiny::req(rx$fo); t <- r$table; o <- graphics::par(mar = c(4.5, 4.5, 1, 1)); e <- r$envelope
      graphics::plot(r$sim$He, r$sim$Fst, pch = 16, cex = 0.35, col = grDevices::adjustcolor("grey60", 0.25), xlab = "Heterozygosity (He)", ylab = "Fst", xlim = c(0, 0.5), ylim = c(min(0, min(t$Fst, na.rm = TRUE)), max(t$Fst, r$sim$Fst, na.rm = TRUE)))
      k <- is.finite(e$q95); graphics::lines(e$He[k], e$q95[k], col = "#e08a00", lwd = 2); k <- is.finite(e$q99); graphics::lines(e$He[k], e$q99[k], col = "#d62728", lwd = 2); k <- is.finite(e$median); graphics::lines(e$He[k], e$median[k], col = "#2a6fdb", lwd = 2)
      graphics::points(t$He, t$Fst, pch = 19, cex = 1.1, col = ifelse(t$P_high < 0.05, "#d62728", ifelse(t$P_low < 0.05, "#7b3fbf", "black")))
      graphics::legend("topright", c("observed sites", "median / 95% / 99% of the neutral simulation"), pch = c(19, NA), lty = c(NA, 1), col = c("black", "#e08a00"), bty = "n"); graphics::par(o)
    }
    output$fo_plot <- shiny::renderPlot({ fo_draw() })
    output$fo_png <- shiny::downloadHandler("fst_outliers.png", function(file) { grDevices::png(file, 1800, 1300, res = 220); on.exit(grDevices::dev.off()); fo_draw() })
    output$fo_csv <- shiny::downloadHandler("fst_outliers.csv", function(file) { r <- shiny::req(rx$fo); utils::write.csv(r$table, file, row.names = FALSE) })
    output$fo_tbl <- DT::renderDataTable({ r <- shiny::req(rx$fo); dtt(fmt(r$table, c("Fst", "He", "P_high", "P_low"), 4), rownames = FALSE, options = list(pageLength = 10, dom = "lftp", scrollX = TRUE)) })

    # ---- sliding window -------------------------------------------------------------------------------------------------------
    shiny::observeEvent(input$run_sl, {
      d <- rv$div
      if (is.null(d)) { set_msg_at("slide", "Run 'Haplotypes, diversity, differentiation' (section 4) first: the window uses the same sites and populations.", "warning"); return() }
      w <- nz(input$sl_w, 100); st <- nz(input$sl_s, 25)
      if (w < 10 || st < 1) { set_msg_at("slide", "The window must have at least 10 columns and the step at least 1.", "warning"); return() }
      tab <- tryCatch(pg_sliding(d$M, d$cols, d$total, d$pop, w, st, nz(input$sl_min, 10), isTRUE(input$sl_pop)), error = function(e) { set_msg_at("slide", paste("Sliding window:", conditionMessage(e)), "danger"); NULL })
      if (is.null(tab) || !any(is.finite(tab$Pi))) { set_msg_at("slide", "No window has enough sites. Lower the minimum number of sites or enlarge the window.", "warning"); rv$slide <- NULL; return() }
      rv$slide <- list(tab = tab, w = w, st = st, total = d$total)
      set_msg_at("slide", sprintf("%d windows of %d columns (step %d) in %d columns of the alignment.", length(unique(tab$Start)), w, st, d$total), "success")
    })
    sl_draw <- function() {
      r <- shiny::req(rv$slide); tab <- r$tab; v <- nz(input$sl_stat, "Pi")
      lab <- c(Pi = "Nucleotide diversity (Pi per site)", ThetaW = "Watterson's theta per site", TajimaD = "Tajima's D", S = "Segregating sites (S)", Hd = "Haplotype diversity (Hd)", h = "Number of haplotypes")[[v]]
      us <- unique(tab$Population); col <- c("#222222", grDevices::hcl.colors(max(1, length(us) - 1), "Dark 3"))[seq_along(us)]
      y <- tab[[v]]; ok <- is.finite(y)
      if (!any(ok)) { graphics::plot.new(); graphics::text(0.5, 0.5, "This statistic is not defined in any window."); return(invisible()) }
      graphics::par(mar = c(4.2, 4.6, 1.2, if (length(us) > 1) 9.5 else 1), xpd = FALSE)
      graphics::plot(NA, xlim = c(1, r$total), ylim = range(y[ok]), xlab = "Alignment column (centre of the window)", ylab = lab, las = 1)
      if (v == "TajimaD") graphics::abline(h = 0, col = "grey60", lty = 2)
      for (i in seq_along(us)) { z <- tab[tab$Population == us[i], ]; graphics::lines(z$Mid, z[[v]], col = col[i], lwd = if (i == 1) 2.2 else 1.5, type = "o", pch = 16, cex = 0.5) }
      if (length(us) > 1) graphics::legend(graphics::par("usr")[2], graphics::par("usr")[4], us, col = col, lwd = 2, bty = "n", cex = 0.85, xpd = NA)
    }
    output$sl_plot <- shiny::renderPlot({ sl_draw() }, height = function() nz(input$sl_h, 340), res = 96)
    output$sl_png <- shiny::downloadHandler("sliding_window.png", function(file) { grDevices::png(file, width = 9, height = nz(input$sl_h, 340) / 96 * 1.0, units = "in", res = 200); on.exit(grDevices::dev.off()); sl_draw() })
    output$sl_csv <- shiny::downloadHandler("sliding_window.csv", function(file) { r <- shiny::req(rv$slide); utils::write.csv(r$tab, file, row.names = FALSE) })
    output$sl_tbl <- DT::renderDataTable({
      r <- shiny::req(rv$slide)
      dtt(fmt(fmt(r$tab, c("Pi", "ThetaW"), 5), c("TajimaD", "Hd"), 3), rownames = FALSE, options = list(pageLength = 10, dom = "tp", scrollX = TRUE))
    })
    output$exact_text <- shiny::renderText({
      d <- rv$div; if (is.null(d)) return("Run 'Haplotypes, diversity, differentiation' first.")
      e <- d$exact; if (is.null(e) || is.na(e$p)) return("Not computed (set the number of Markov chain steps above 0, with at least 2 populations).")
      sprintf("All populations: p = %.5f (s.e. %.5f); %d steps of the Markov chain after the dememorization.", e$p, e$se, e$steps)
    })
    output$exact_pw <- DT::renderDataTable({
      d <- rv$div; shiny::req(d$exact$pairwise)
      m <- formatC(d$exact$pairwise, digits = 4, format = "f"); m[row(m) == col(m)] <- "-"
      dtt(cbind(Population = rownames(m), as.data.frame(m, stringsAsFactors = FALSE)), rownames = FALSE, options = list(dom = "t", scrollX = TRUE))
    })
    output$hap_tbl <- DT::renderDataTable({
      d <- rv$div; shiny::req(d)
      dtt(d$haptab, rownames = FALSE, options = list(pageLength = 10, lengthMenu = LENGTHS, dom = "lftp", scrollX = TRUE))
    })
    output$hap_text <- shiny::renderText({
      d <- rv$div; if (is.null(d)) return("Run 'Haplotypes, diversity, differentiation' first.")
      paste(pg_hap_text(d$hap$hap, rownames(d$M)), collapse = "\n")
    })

    # ---- genetic differentiation text and tables (DnaSP layout) ----------------------------------------------------------
    dnasp_text <- function(d) {
      g <- d$gf; if (is.null(g)) return("Genetic differentiation was not computed (see the message above).")
      pops <- g$pops; fm <- function(x, k = 5) formatC(x, digits = k, format = "f")
      out <- c("Gene Flow and Genetic Differentiation", sprintf("Number of populations included: %d", nrow(pops) - 1),
               sprintf("Sites with alignment gaps are: %s", c(excluded = "Excluded", fifth = "Considered (fifth state)", pairwise = "Excluded only in pairwise comparisons")[[g$gap]]),
               sprintf("Total sites used: %d", g$sites), sprintf("Genome: %s", c(diploid = "Diploid-Autosome", haploid = "Haploid", X = "X-linked", Y = "Y-linked")[[g$genome]]), "")
      for (i in seq_len(nrow(pops))) {
        p <- pops[i, ]
        out <- c(out, if (p$Population == "Total") "Total Data Estimates" else sprintf("Population %d: %s", i, p$Population),
                 sprintf("  Number of sequences: %d", p$n), sprintf("  Number of segregating sites, S: %d", p$S), sprintf("  Number of haplotypes, h: %d", p$h),
                 sprintf("  Haplotype diversity, Hd: %s", fm(p$Hd)), sprintf("  Average number of differences, K: %s", fm(p$K)),
                 sprintf("  Nucleotide diversity, Pi: %s", fm(p$Pi)), sprintf("  Nucleotide diversity with JC, PiJC: %s", fm(p$PiJC)), "")
      }
      te <- g$tests
      out <- c(out, "Genetic differentiation estimates",
               sprintf("  Chi-square (table) Chi2: %s  P-value (chi-square distribution, df = %d): %s%s", fm(te$Value[1], 3), g$df, fm(te$P_chi2_distribution[1], 4),
                       if (!is.na(te$P_permutation[1])) sprintf("   permutation P: %s", fm(te$P_permutation[1], 4)) else ""),
               sprintf("  Hs: %s  Hst: %s%s", fm(te$Value[2]), fm(te$Value[3]), if (!is.na(te$P_permutation[3])) sprintf("   permutation P (Hst): %s", fm(te$P_permutation[3], 4)) else ""),
               sprintf("  Ks: %s  Kst: %s%s", fm(te$Value[4]), fm(te$Value[5]), if (!is.na(te$P_permutation[5])) sprintf("   permutation P (Kst): %s", fm(te$P_permutation[5], 4)) else ""),
               sprintf("  Snn: %s%s", fm(te$Value[6]), if (!is.na(te$P_permutation[6])) sprintf("   permutation P: %s", fm(te$P_permutation[6], 4)) else ""), "",
               "Gene flow estimates",
               sprintf("  Gst (Hudson et al. 1992 form): %s   Nm: %s", fm(te$Value[7]), fm(te$Nm[7], 2)),
               sprintf("  Fst (Hudson, Slatkin & Maddison 1992): %s   Nm: %s", fm(te$Value[8]), fm(te$Nm[8], 2)),
               sprintf("  Nst (Lynch & Crease 1990, JC): %s   Nm: %s", fm(te$Value[9]), fm(te$Nm[9], 2)))
      out
    }
    output$gf_text <- shiny::renderText({ d <- rv$div; if (is.null(d)) "Run 'Haplotypes, diversity, differentiation' first." else paste(dnasp_text(d), collapse = "\n") })
    output$gf_tests <- DT::renderDataTable({ d <- rv$div; shiny::req(d$gf); dtt(fmt(fmt(d$gf$tests, c("Value", "Nm"), 5), c("P_chi2_distribution", "P_permutation"), 4), rownames = FALSE, options = list(dom = "t")) })
    output$gf_pair <- DT::renderDataTable({
      d <- rv$div; shiny::req(d$gf)
      dtt(fmt(d$gf$pairwise, setdiff(names(d$gf$pairwise), c("Population1", "Population2")), 5), rownames = FALSE, options = list(pageLength = 15, dom = "tp", scrollX = TRUE))
    })

    # ---- mismatch ----------------------------------------------------------------------------------------------
    output$mm_unit_ui <- shiny::renderUI({ d <- rv$div; shiny::req(d); shiny::selectInput(session$ns("mm_unit"), "Population", names(d$mm)[!vapply(d$mm, is.null, TRUE)], width = "260px") })
    output$mm_plot <- shiny::renderPlot({
      d <- rv$div; shiny::req(d, input$mm_unit); m <- d$mm[[input$mm_unit]]; shiny::req(m)
      graphics::barplot(m$freq, names.arg = m$diff, col = "#8db4e2", border = "#5f86b8", xlab = "Pairwise differences", ylab = "Frequency",
                        main = sprintf("%s: mean %.2f differences, raggedness r = %.4f", input$mm_unit, m$mean, m$raggedness))
    })

    # ---- network ----------------------------------------------------------------------------------------------------
    shiny::observeEvent(input$run_net, {
      d <- rv$div
      if (is.null(d)) { set_msg_at("net", "Run 'Haplotypes, diversity, differentiation' first.", "warning"); return() }
      hp <- d$hap; ids <- hp$id
      if (length(ids) < 2) { set_msg_at("net", "A network needs at least 2 haplotypes.", "warning"); return() }
      prog <- .pg_progress(session, "Haplotype network"); on.exit(prog$close(), add = TRUE); prog$set(15, "building")
      method <- nz(input$net_method, "msn"); note <- NULL
      Mh <- do.call(rbind, strsplit(hp$seq, "")); rownames(Mh) <- ids
      nodes_seq <- hp$seq; node_lab <- ids; is_med <- rep(FALSE, length(ids)); edges <- NULL
      if (method == "mjn") {
        mj <- tryCatch(pg_net_mj(hp$seq, nz(input$net_eps, 0)), error = function(e) { note <<- paste0("The median-joining network could not be built (", conditionMessage(e), "): the minimum spanning network is shown instead."); NULL })
        if (!is.null(mj)) {
          mvn <- 0L; nl <- vapply(seq_along(mj$med), function(i) if (mj$med[i]) { mvn <<- mvn + 1L; paste0("mv", mvn) } else ids[i], "")
          nodes_seq <- mj$seq; node_lab <- nl; is_med <- mj$med; edges <- mj$edges
        }
      } else if (method == "haplo") {
        if (!requireNamespace("pegas", quietly = TRUE)) { note <- "The pegas package is not installed (install.packages(\"pegas\")): the minimum spanning network is shown instead." }
        else {
          pgres <- tryCatch({
            h <- pegas::haplotype(.pg_dnabin(Mh)); hn <- pegas::haploNet(h)
            em <- as.matrix(hn); idx <- attr(h, "index")
            if (ncol(em) < 2) stop("unexpected structure")
            map <- vapply(idx, function(v) ids[as.integer(v[1])], "")
            if (any(is.na(map)) || any(em[, 1:2] < 1) || any(em[, 1:2] > length(map))) stop("link indices out of range")
            ord <- match(ids, map)
            list(seq = hp$seq, lab = ids, med = rep(FALSE, length(ids)),
                 edges = data.frame(a = ord[as.integer(em[, 1])], b = ord[as.integer(em[, 2])], stringsAsFactors = FALSE))
          }, error = function(e) { note <<- paste0("The pegas network could not be read (", conditionMessage(e), "): the minimum spanning network is shown instead."); NULL })
          if (!is.null(pgres)) { nodes_seq <- pgres$seq; node_lab <- pgres$lab; is_med <- pgres$med; edges <- pgres$edges }
        }
      }
      if (is.null(edges)) edges <- pg_net_msn(hp$seq)
      edges <- edges[edges$a != edges$b, , drop = FALSE]
      S <- do.call(rbind, strsplit(nodes_seq, ""))
      steps <- vapply(seq_len(nrow(edges)), function(i) sum(S[edges$a[i], ] != S[edges$b[i], ]), 1L)
      posl <- lapply(seq_len(nrow(edges)), function(i) as.integer(d$cols[which(S[edges$a[i], ] != S[edges$b[i], ])]))
      counts <- d$counts; pops <- colnames(counts); cols <- .pg_cols(length(pops))
      nodes <- lapply(seq_along(node_lab), function(i) {
        k <- match(node_lab[i], rownames(counts))
        list(label = node_lab[i], n = if (is_med[i] || is.na(k)) 0L else sum(counts[k, ]), counts = if (is_med[i] || is.na(k)) rep(0L, length(pops)) else as.integer(counts[k, ]), median = is_med[i])
      })
      prog$set(80, "drawing")
      session$sendCustomMessage("netv_load", list(id = session$ns("net"), nodes = nodes,
        edges = lapply(seq_len(nrow(edges)), function(i) list(a = edges$a[i] - 1L, b = edges$b[i] - 1L, steps = steps[i], pos = posl[[i]])),
        pops = as.list(pops), colors = as.list(cols)))
      prog$set(100, "done"); rv$net_built <- TRUE; set_msg_at("net", note, "info")
    })

    # ---- AMOVA and pairwise FST (Arlequin) -------------------------------------------------------------------------------
    shiny::observeEvent(input$run_amova, {
      cl <- tryCatch(clean(), error = function(e) NULL)
      if (is.null(cl)) { set_msg_at("run", "Provide an alignment with at least 3 sequences first.", "warning"); return() }
      ms <- nz(input$miss, 0)
      M <- if (ms > 0) pg_prepare(prep()$aln, "level", ms)$M else cl$M
      if (ncol(M) == 0) { set_msg_at("run", "No site is free of gaps / N / IUPAC codes in every sequence (or within the allowed level of missing data). Trim the alignment further first.", "warning"); return() }
      pp <- pop_for(M); pop <- factor(pp$pop)
      if (nlevels(pop) < 2) { set_msg_at("run", "Define at least 2 populations (section 2).", "warning"); return() }
      if (length(pop) <= nlevels(pop)) { set_msg_at("run", sprintf("The AMOVA needs more sequences than populations (%d sequences, %d populations). If every sequence is its own population, the sets were not applied: check the population table.", length(pop), nlevels(pop)), "danger"); return() }
      prog <- .pg_progress(session, "AMOVA"); on.exit(prog$close(), add = TRUE); prog$set(3, "distances")
      g <- nz(input$gamma, 0)
      D <- tryCatch(as.matrix(ape::dist.dna(.pg_dnabin(M), model = nz(input$dmodel, "N"), gamma = if (g > 0 && nz(input$dmodel, "N") != "N") g else FALSE, pairwise.deletion = ms > 0, as.matrix = TRUE)), error = function(e) NULL)
      if (is.null(D) || any(!is.finite(D))) { set_msg_at("run", "Some distances cannot be computed with this method (saturated sequences). Use 'Pairwise difference' or Jukes-Cantor.", "danger"); return() }
      nperm <- as.integer(nz(input$nperm, 1000))
      am <- tryCatch(pg_amova(D, pop, pp$grp, nperm, progress = function(f) prog$set(5 + 55 * f, "AMOVA permutations")), error = function(e) { set_msg_at("run", conditionMessage(e), "danger"); NULL })
      if (is.null(am)) return()
      pw <- pg_pairwise_phist(D, pop, as.integer(nz(input$nperm_pair, 200)), progress = function(f) prog$set(60 + 40 * f, "pairwise FST"))
      rv$amova <- list(am = am, pw = pw, nperm_pair = as.integer(nz(input$nperm_pair, 200)), model = nz(input$dmodel, "N"), gamma = g, n = nrow(M), sites = ncol(M), npop = nlevels(pop), alpha = nz(input$alpha, 0.05), slatkin = isTRUE(input$slatkin),
                       reyn = isTRUE(input$reyn), pipop = if (isTRUE(input$pipop)) pg_pop_pi(D, pop) else NULL, miss = ms)
      prog$set(100, "done"); set_msg_at("run", NULL)
    })
    amova_df <- function(am) {
      o <- am$obs
      if (am$levels == 2) { sg <- unname(o$sigma[c("a", "w")]); df <- unname(o$df[c("a", "w")]); ssd <- unname(o$ssd[c("a", "w")]); src <- c("Among populations", "Within populations") }
      else { sg <- unname(o$sigma[c("a", "b", "c")]); df <- unname(o$df[c("ag", "ap", "wp")]); ssd <- unname(o$ssd[c("ag", "ap", "wp")])
             src <- c("Among groups", "Among populations within groups", "Within populations") }
      tot <- sum(sg)
      data.frame(Source_of_variation = c(src, "Total"), d.f. = c(df, sum(df)), Sum_of_squares = round(c(ssd, sum(ssd)), 4),
                 Variance_components = round(c(sg, tot), 5), Percentage_of_variation = round(100 * c(sg, tot) / tot, 2), stringsAsFactors = FALSE)
    }
    phi_df <- function(am) {
      o <- am$obs
      if (am$levels == 2) data.frame(Fixation_index = "FST", Value = round(o$phiST, 5), P = round(am$p[["phiST"]], 4), stringsAsFactors = FALSE)
      else data.frame(Fixation_index = c("FSC", "FST", "FCT"), Value = round(c(o$phiSC, o$phiST, o$phiCT), 5),
                      P = round(c(am$p[["phiSC"]], am$p[["phiST"]], am$p[["phiCT"]]), 4), stringsAsFactors = FALSE)
    }
    pw_df <- function(a) {
      pw <- a$pw; ph <- pw$phi; if (isTRUE(a$slatkin)) ph <- ph / (1 - ph)
      m <- matrix("", nrow(ph), ncol(ph), dimnames = dimnames(ph)); lo <- lower.tri(m); up <- upper.tri(m)
      m[lo] <- formatC(ph[lo], digits = 5, format = "f"); m[up] <- paste0(formatC(pw$p[up], digits = 4, format = "f"), ifelse(!is.na(pw$p[up]) & pw$p[up] < a$alpha, " +", ""))
      m[row(m) == col(m)] <- "-"
      cbind(Population = rownames(m), as.data.frame(m, stringsAsFactors = FALSE))
    }
    output$amova_head <- shiny::renderUI({
      a <- rv$amova; if (is.null(a)) return(shiny::helpText("Press 'AMOVA and pairwise FST'."))
      shiny::helpText(sprintf("%d sequences, %d sites, %d populations; %d-level AMOVA; distance: %s%s; %d permutations.", a$n, a$sites, a$npop, a$am$levels,
                              c(N = "pairwise difference", K80 = "Kimura 2-parameter", TN93 = "Tamura-Nei", JC69 = "Jukes-Cantor")[[a$model]], if (a$gamma > 0 && a$model != "N") sprintf(", gamma a = %s", a$gamma) else "", a$am$nperm))
    })
    output$amova_tbl <- DT::renderDataTable({ a <- rv$amova; shiny::req(a); dtt(amova_df(a$am), rownames = FALSE, options = list(dom = "t")) })
    output$phi_tbl <- DT::renderDataTable({ a <- rv$amova; shiny::req(a); dtt(phi_df(a$am), rownames = FALSE, options = list(dom = "t")) })
    output$amova_extra <- shiny::renderUI({
      a <- rv$amova; shiny::req(a)
      mt <- function(m, f) { m <- as.data.frame(apply(m, 2, f), stringsAsFactors = FALSE); cbind(Population = rownames(a$pw$phi), m) }
      tagl <- list()
      if (isTRUE(a$reyn)) { r <- pg_reynolds(a$pw$phi); r[upper.tri(r, TRUE)] <- NA
        tagl <- c(tagl, list(shiny::h5(style = "margin-top:14px", "Reynolds' distance, -ln(1 - FST)"), shiny::HTML(.pg_html_table(mt(r, function(x) ifelse(is.na(x), "", formatC(x, digits = 5, format = "f"))))))) }
      if (!is.null(a$pipop)) { pp <- a$pipop; m <- pp$net; m[upper.tri(m)] <- pp$xy[upper.tri(m)]; diag(m) <- diag(pp$xy)
        tagl <- c(tagl, list(shiny::h5(style = "margin-top:14px", "Pairwise differences: within populations (diagonal), between populations PiXY (above), corrected PiXY - (PiX + PiY) / 2 (below)"),
                             shiny::HTML(.pg_html_table(mt(m, function(x) formatC(x, digits = 5, format = "f")))))) }
      if (isTRUE(a$miss > 0)) tagl <- c(tagl, list(shiny::helpText(sprintf("Allowed level of missing data %g: %d sites used for the distances; missing entries are skipped pair by pair.", a$miss, a$sites))))
      do.call(shiny::tagList, tagl)
    })
    output$pw_tbl <- DT::renderDataTable({ a <- rv$amova; shiny::req(a); dtt(pw_df(a), rownames = FALSE, options = list(dom = "t", scrollX = TRUE)) })


    # ---- Structure: K clusters, K by population and K by individual -------------------------------------------------------------
    cur_pop <- function(ids) {
      d <- rv$pop; p <- if (is.null(d)) rep(NA_character_, length(ids)) else trimws(d$Population[match(ids, d$Sequence)])
      p[is.na(p) | !nzchar(p)] <- "unassigned"; p
    }
    shiny::observeEvent(input$run_struct, {
      cl <- tryCatch(clean(), error = function(e) NULL)
      if (is.null(cl)) { set_msg_at("run", "Provide an alignment with at least 3 sequences first.", "warning"); return() }
      M <- cl$M
      if (ncol(M) == 0) { set_msg_at("run", "No site is free of gaps / N / IUPAC codes in every sequence. Trim the alignment further first.", "warning"); return() }
      if (nrow(M) < 4) { set_msg_at("run", "Structure needs at least 4 sequences.", "warning"); return() }
      k1 <- max(1L, as.integer(nz(input$st_kmin, 1))); k2 <- min(as.integer(nz(input$st_kmax, 5)), nrow(M) - 1L)
      if (k2 < k1) { set_msg_at("run", sprintf("Choose K from 1 to at most %d (the number of sequences minus 1).", nrow(M) - 1L), "warning"); return() }
      prog <- .pg_progress(session, "Structure"); on.exit(prog$close(), add = TRUE)
      res <- tryCatch(pg_structure(M, k1:k2, max(1L, as.integer(nz(input$st_rep, 5))), nz(input$st_model, "admix"), as.integer(nz(input$st_iter, 300)),
                                   if (is.na(nz(input$st_seed, NA))) NULL else input$st_seed, progress = function(f, txt) prog$set(100 * f, txt),
                                   cv_frac = max(0, min(30, nz(input$st_cv, 5))) / 100,
                                   method = nz(input$st_method, "em"), burnin = max(100L, as.integer(nz(input$st_burn, 2000))), nrep = max(100L, as.integer(nz(input$st_nrep, 3000))),
                                   alpha = max(0.01, min(10, nz(input$st_alpha, 1))), infer_alpha = !identical(input$st_alpha_inf, FALSE)),
                      error = function(e) { set_msg_at("run", conditionMessage(e), "danger"); NULL })
      if (is.null(res)) return()
      rv$struct <- list(res = res, ext = FALSE, sites = res$S, n = nrow(M))
      prog$set(100, "done"); set_msg_at("run", NULL)
      shiny::updateTabsetPanel(session, "tabs", selected = "Structure")
    })
    shiny::observeEvent(input$st_q_csv, {
      p <- tryCatch(prep(), error = function(e) NULL)
      if (is.null(p)) { set_msg_at("struct", "Provide the alignment first (section 1): its sequence names are used to match the rows of the Q file.", "warning"); return() }
      f <- input$st_q_csv$datapath
      d <- NULL
      for (sep in c(",", "\t", ";", " ")) {
        d <- tryCatch(utils::read.table(f, sep = sep, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, fileEncoding = "UTF-8-BOM", comment.char = "", quote = "\"", strip.white = TRUE), error = function(e) NULL)
        if (!is.null(d) && ncol(d) >= 3) break
      }
      if (is.null(d) || ncol(d) < 3) { set_msg_at("struct", "The Q file needs a first column with the sequence names and at least two columns of numbers (one per cluster).", "warning"); return() }
      num <- vapply(d[-1], is.numeric, TRUE)
      if (sum(num) < 2) { set_msg_at("struct", "The Q file needs at least two numeric columns (one per cluster).", "warning"); return() }
      mt <- pg_match_keys(names(p$aln), d[[1]]); ok <- !is.na(mt$index)
      if (!any(ok)) { set_msg_at("struct", sprintf("None of the %d sequences matched the Q file. First keys of the file: %s. First sequence names: %s.", length(p$aln),
                                      paste(utils::head(d[[1]], 3), collapse = ", "), paste(utils::head(names(p$aln), 3), collapse = ", ")), "warning"); return() }
      Q <- as.matrix(d[-1][mt$index[ok], num, drop = FALSE]); Q[!is.finite(Q)] <- 0
      Q <- Q / pmax(rowSums(Q), 1e-12)
      dimnames(Q) <- list(names(p$aln)[ok], paste0("Q", seq_len(ncol(Q))))
      rv$struct <- list(res = list(Q = stats::setNames(list(Q), ncol(Q)), ll = NULL, runs = NULL, model = "file", n = nrow(Q), S = NA_integer_), ext = TRUE, sites = NA_integer_, n = nrow(Q))
      set_msg_at("struct", sprintf("Q file read: %d of %d sequences matched, K = %d.", sum(ok), length(p$aln), ncol(Q)), if (all(ok)) "success" else "warning")
      shiny::updateTabsetPanel(session, "tabs", selected = "Structure")
    })
    output$st_head <- shiny::renderUI({
      s <- rv$struct
      if (is.null(s)) return(shiny::helpText("Press 'Structure (K clusters)' (left, section 4), or load a Q table from STRUCTURE / ADMIXTURE below."))
      if (s$ext) return(shiny::helpText(sprintf("Q table from a file: %d sequences, K = %d.", s$n, ncol(s$res$Q[[1]]))))
      r <- s$res; ks <- as.integer(names(r$Q))
      shiny::helpText(sprintf("%d sequences, %d variable sites; %s model, %s; K = %d to %d; %d run(s) per K (the run with the best likelihood is drawn).", s$n, s$sites,
                              if (r$model == "admix") "admixture" else "no-admixture",
                              if (identical(r$method, "mcmc")) sprintf("Bayesian MCMC (%d burn-in + %d sweeps)", r$burnin, r$nrep) else "maximum likelihood (EM)", min(ks), max(ks), r$reps))
    })
    output$st_k_ui <- shiny::renderUI({
      s <- shiny::req(rv$struct); ks <- as.integer(names(s$res$Q))
      ce <- if (is.null(s$res$ll)) NA else s$res$ll$cross_entropy; dk <- if (is.null(s$res$ll)) NA else s$res$ll$deltaK
      def <- if (any(is.finite(ce))) ks[which.min(ce)] else if (any(is.finite(dk))) ks[which.max(dk)] else max(ks)
      old <- shiny::isolate(input$st_k)
      shiny::selectInput(session$ns("st_k"), "K to show", ks, selected = if (!is.null(old) && as.integer(old) %in% ks) as.integer(old) else def)
    })
    st_cur <- shiny::reactive({
      s <- shiny::req(rv$struct); K <- shiny::req(input$st_k); Q <- s$res$Q[[as.character(K)]]; shiny::req(!is.null(Q))
      list(Q = Q, pop = cur_pop(rownames(Q)), K = as.integer(K))
    })
    st_draw <- function() {
      x <- st_cur(); od <- trimws(strsplit(nz(input$st_order, ""), ",")[[1]]); od <- od[nzchar(od)]
      pg_structure_plot(x$Q, x$pop, pg_str_palette(x$K, nz(input$st_palette, "figure"), nz(input$st_colors, "")), nz(input$st_mode, "both"), od,
                        nz(input$st_title, ""), nz(input$st_cex, 0.85), isTRUE(input$st_sort), isTRUE(input$st_names))
    }
    output$st_plot <- shiny::renderPlot({ st_draw() }, height = function() as.integer(nz(input$st_h, 460)), res = 96)
    st_file <- function(dev, ext) shiny::downloadHandler(function() paste0("structure_K", nz(input$st_k, ""), ext), function(file) {
      h <- max(2.5, nz(input$st_h, 460) / 100)
      if (dev == "png") grDevices::png(file, width = 11, height = h, units = "in", res = 200)
      else if (dev == "pdf") grDevices::pdf(file, width = 11, height = h)
      else grDevices::svg(file, width = 11, height = h)
      on.exit(grDevices::dev.off()); st_draw()
    })
    output$st_png <- st_file("png", ".png"); output$st_pdf <- st_file("pdf", ".pdf"); output$st_svg <- st_file("svg", ".svg")
    output$st_ll_plot <- shiny::renderPlot(height = function() as.integer(nz(input$st_llh, 270)), expr = {
      s <- shiny::req(rv$struct); ll <- s$res$ll; shiny::req(!is.null(ll), nrow(ll) > 0)
      op <- graphics::par(mfrow = c(1, 4), mar = c(4, 4.4, 2.2, 1)); on.exit(graphics::par(op))
      one <- function(y, ylab, main, best = NULL, sd = NULL) {
        if (!any(is.finite(y))) { graphics::plot.new(); graphics::text(0.5, 0.5, "not available", cex = 0.9); graphics::title(main); return() }
        graphics::plot(ll$K, y, type = "b", pch = 19, col = "#1f4e8c", xlab = "Number of clusters K", ylab = ylab, main = main, xaxt = "n"); graphics::axis(1, ll$K)
        if (!is.null(sd)) { w <- which(is.finite(sd) & sd > 0.004 * max(diff(range(y, na.rm = TRUE)), 1e-12)); if (length(w)) graphics::arrows(ll$K[w], y[w] - sd[w], ll$K[w], y[w] + sd[w], angle = 90, code = 3, length = 0.04, col = "#1f4e8c") }
        if (!is.null(best) && is.finite(best)) graphics::points(ll$K[best], y[best], pch = 1, cex = 2.2, col = "#c0392b")
      }
      one(ll$cross_entropy, "Cross-entropy (lower = better)", "Cross-entropy", if (any(is.finite(ll$cross_entropy))) which.min(ll$cross_entropy) else NULL, ll$sd_cross_entropy)
      one(ll$mean_lnL, "ln L (mean of the runs)", "Likelihood", NULL, ll$sd_lnL)
      one(ll$BIC, "BIC (lower = better)", "BIC", if (any(is.finite(ll$BIC))) which.min(ll$BIC) else NULL)
      if (any(is.finite(ll$deltaK))) one(ll$deltaK, "delta K", "Evanno's delta K", which.max(ll$deltaK))
      else { graphics::plot.new(); graphics::text(0.5, 0.5, "delta K needs at least 3\nconsecutive K and 2 runs per K", cex = 0.85); graphics::title("Evanno's delta K") }
    }, res = 96)
    rv_coords <- shiny::reactiveVal(NULL)
    read_coords <- function(path, name, area) {
      d <- NULL
      if (grepl("\\.xlsx?$", name, ignore.case = TRUE)) { sh <- tryCatch(pg_read_xlsx(path), error = function(e) NULL); if (!is.null(sh) && length(sh)) d <- sh[[1]] }
      else for (sep in c(",", "\t", ";")) {
        d <- tryCatch(utils::read.table(path, sep = sep, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, fileEncoding = "UTF-8-BOM", comment.char = "", quote = "\"", strip.white = TRUE), error = function(e) NULL)
        if (!is.null(d) && ncol(d) >= 3) break
      }
      if (is.null(d) || ncol(d) < 3) { set_msg_at(area, "The coordinates file needs the columns Population, Latitude and Longitude.", "warning"); return() }
      nm <- tolower(trimws(names(d)))
      pc <- which(nm %in% c("population", "pop", "site", "locality", "name"))[1]; if (is.na(pc)) pc <- 1
      la <- which(nm %in% c("latitude", "lat", "y"))[1]; lo <- which(nm %in% c("longitude", "lon", "long", "lng", "x"))[1]
      if (is.na(la) || is.na(lo)) { set_msg_at(area, "Could not find the columns Latitude and Longitude in the coordinates file.", "warning"); return() }
      out <- data.frame(Population = trimws(as.character(d[[pc]])), lon = suppressWarnings(as.numeric(gsub(",", ".", d[[lo]]))), lat = suppressWarnings(as.numeric(gsub(",", ".", d[[la]]))), stringsAsFactors = FALSE)
      out <- out[is.finite(out$lon) & is.finite(out$lat), , drop = FALSE]
      if (!nrow(out)) { set_msg_at(area, "No row of the coordinates file has numeric Latitude and Longitude.", "warning"); return() }
      rv_coords(out); set_msg_at(area, sprintf("Coordinates read for %d population(s).", nrow(out)), "success")
    }
    shiny::observeEvent(input$st_coord, read_coords(input$st_coord$datapath, input$st_coord$name, "struct"))
    shiny::observeEvent(input$map_coord, read_coords(input$map_coord$datapath, input$map_coord$name, "map"))
    # ---- map of the haplotypes ----------------------------------------------------------------------------------------------------
    map_msg <- function() {
      if (is.null(rv$div)) return(list(id = session$ns("hmap"), empty = "Run 'Haplotypes, diversity, differentiation' first."))
      cd <- rv_coords(); if (is.null(cd)) return(list(id = session$ns("hmap"), empty = "Load the coordinates of the populations (CSV or Excel) to draw the map."))
      cnt <- rv$div$counts; hn <- rownames(cnt); pn <- colnames(cnt)
      cols <- .pg_hapcols(nrow(cnt)); uc <- trimws(unlist(strsplit(nz(input$map_cols, ""), "[,;\\s]+"))); uc <- uc[nzchar(uc)]
      if (length(uc)) { k <- seq_len(min(length(uc), length(cols))); cols[k] <- uc[k] }
      m <- match(tolower(pn), tolower(cd$Population)); ok <- !is.na(m)
      if (!any(ok)) return(list(id = session$ns("hmap"), empty = "None of the populations has coordinates: the names in the file must be the same as in section 2."))
      if (any(!ok)) set_msg_at("map", sprintf("No coordinates for: %s (not drawn).", paste(pn[!ok], collapse = ", ")), "warning") else set_msg_at("map", NULL)
      pops <- lapply(which(ok), function(j) list(name = pn[j], lon = cd$lon[m[j]], lat = cd$lat[m[j]], counts = as.list(as.numeric(cnt[, j])), n = sum(cnt[, j])))
      pops <- Filter(function(p) p$n > 0, pops)
      list(id = session$ns("hmap"), pops = unname(pops), haps = as.list(hn), colors = as.list(cols),
           opt = list(base = nz(input$map_base, "sat"), pie = nz(input$map_pie, 5), height = as.integer(nz(input$map_h, 560)), names = !identical(input$map_lab, FALSE),
                      hapLegend = !identical(input$map_hleg, FALSE), sizeLegend = !identical(input$map_sleg, FALSE), scale = !identical(input$map_scale, FALSE),
                      title = nz(input$map_title, ""), land = nz(input$map_land, "#e8e4d8"), sea = nz(input$map_sea, "#dbe9f4")))
    }
    shiny::observe({
      msg <- map_msg(); session$sendCustomMessage("hmap_load", jsonlite::fromJSON(jsonlite::toJSON(msg, auto_unbox = TRUE, null = "null"), simplifyVector = FALSE))
    })
    # the files of the map (PNG + SVG made in the browser) when everything is saved
    shiny::observeEvent(input$hmap_files, {
      x <- input$hmap_files; o <- rv$last_out; if (is.null(o) || is.null(x$svg)) return()
      f1 <- file.path(o$dir, paste0(o$name, "_haplotype_map.svg")); f2 <- file.path(o$dir, paste0(o$name, "_haplotype_map.png"))
      tryCatch(writeLines(x$svg, f1, useBytes = TRUE), error = function(e) NULL)
      if (!is.null(x$png)) tryCatch(writeBin(jsonlite::base64_dec(sub("^data:image/png;base64,", "", x$png)), f2), error = function(e) NULL)
    })
    output$map_tpl <- shiny::downloadHandler("coordinates_template.csv", function(file) {
      pp <- if (!is.null(rv$div)) colnames(rv$div$counts) else if (!is.null(rv$pop)) unique(rv$pop$Population) else character()
      cd <- rv_coords(); la <- lo <- rep(NA_real_, length(pp)); if (!is.null(cd)) { m <- match(tolower(pp), tolower(cd$Population)); la <- cd$lat[m]; lo <- cd$lon[m] }
      utils::write.csv(data.frame(Population = pp, Latitude = la, Longitude = lo), file, row.names = FALSE, na = "")
    })
    st_map_draw <- function() {
      x <- st_cur(); cd <- shiny::req(rv_coords()); pops <- unique(x$pop)
      pg_structure_map(pg_str_bypop(x$Q, x$pop, pops), as.integer(table(factor(x$pop, levels = pops))), cd, pg_str_palette(x$K, nz(input$st_palette, "figure"), nz(input$st_colors, "")),
                       nz(input$st_pie, 3), isTRUE(input$st_pie_n), !identical(input$st_pie_lab, FALSE), nz(input$st_cex, 0.85), sprintf("K=%d", x$K))
    }
    output$st_map <- shiny::renderPlot({
      if (is.null(rv_coords())) { graphics::plot.new(); graphics::text(0.5, 0.5, "Load the coordinates of the populations to draw the map", cex = 1.1, col = "grey40"); return(invisible()) }
      tryCatch(st_map_draw(), error = function(e) { graphics::plot.new(); graphics::text(0.5, 0.5, conditionMessage(e), cex = 1, col = "#c0392b") })
    }, height = function() as.integer(nz(input$st_mh, 460)), res = 96)
    output$st_map_png <- shiny::downloadHandler(function() paste0("structure_map_K", nz(input$st_k, ""), ".png"), function(file) { grDevices::png(file, width = 9, height = 6, units = "in", res = 200); on.exit(grDevices::dev.off()); st_map_draw() })
    # ---- Arlequin settings (.ars) ---------------------------------------------------------------------------------
    rv_ars <- shiny::reactiveVal(NULL)
    shiny::observeEvent(input$ars, {
      f <- input$ars$datapath
      txt <- tryCatch(rawToChar(readBin(f, "raw", file.size(f))), error = function(e) NULL)
      kv <- if (is.null(txt)) NULL else tryCatch(pg_ars_parse(iconv(txt, "latin1", "UTF-8", sub = "")), error = function(e) NULL)
      if (is.null(kv) || !length(kv) || is.null(kv[["TaskNumber"]]) && is.null(kv[["NumPermutationsAMOVA"]])) {
        rv_ars(NULL); set_msg_at("set", "That file does not look like an Arlequin settings file (.ars).", "warning"); return()
      }
      r <- pg_ars_to_settings(kv); v <- r$values
      for (id in c("nperm", "nperm_pair", "alpha", "gamma", "nsim", "ew_nsim", "ex_steps", "ex_demem", "miss"))
        if (!is.null(v[[id]])) shiny::updateNumericInput(session, id, value = v[[id]])
      if (!is.null(v$dmodel)) shiny::updateSelectInput(session, "dmodel", selected = v$dmodel)
      if (!is.null(v$slatkin)) shiny::updateCheckboxInput(session, "slatkin", value = isTRUE(v$slatkin))
      if (!is.null(v$reyn)) shiny::updateCheckboxInput(session, "reyn", value = isTRUE(v$reyn))
      if (!is.null(v$pipop)) shiny::updateCheckboxInput(session, "pipop", value = isTRUE(v$pipop))
      rv_ars(list(text = iconv(txt, "latin1", "UTF-8", sub = ""), name = input$ars$name, applied = r$applied, notes = r$notes))
    })
    output$ars_msg <- shiny::renderUI({
      a <- rv_ars(); if (is.null(a)) return(NULL)
      shiny::div(class = "alert alert-success", style = "padding:8px 10px;font-size:13px",
        shiny::strong(sprintf("%s: %d setting(s) applied to the form.", a$name, nrow(a$applied))),
        if (nrow(a$applied)) shiny::tags$ul(style = "margin:4px 0 0 0;padding-left:18px",
          lapply(seq_len(nrow(a$applied)), function(i) shiny::tags$li(sprintf("%s = %s", a$applied$Setting[i], a$applied$Value[i])))),
        if (length(a$notes)) shiny::tagList(shiny::strong(style = "display:block;margin-top:6px", "Not applied / to know:"),
          shiny::tags$ul(style = "margin:4px 0 0 0;padding-left:18px", lapply(a$notes, shiny::tags$li))))
    })
    output$ars_save <- shiny::renderUI({ if (is.null(rv_ars())) NULL else shiny::div(style = "margin-bottom:8px", shiny::downloadButton(session$ns("ars_dl"), "Save the settings (.ars) for Arlequin")) })
    output$ars_dl <- shiny::downloadHandler(
      filename = function() { a <- rv_ars(); if (is.null(a)) "arlequin_settings.ars" else sub("(\\.ars)?$", "_SangeRly.ars", a$name, ignore.case = TRUE) },
      content = function(file) {
        a <- shiny::req(rv_ars())
        out <- pg_ars_write(a$text, list(nperm = nz(input$nperm, 1000), nperm_pair = nz(input$nperm_pair, 1000), alpha = nz(input$alpha, 0.05), dmodel = nz(input$dmodel, "N"),
                                         gamma = nz(input$gamma, 0), nsim = nz(input$nsim, 0), ew_nsim = nz(input$ew_nsim, 0), ex_steps = nz(input$ex_steps, 0),
                                         ex_demem = nz(input$ex_demem, 10000), slatkin = isTRUE(input$slatkin),
                                         reyn = isTRUE(input$reyn), pipop = isTRUE(input$pipop), miss = nz(input$miss, 0)))
        writeBin(charToRaw(iconv(out, "UTF-8", "latin1", sub = "")), file)
      })
    output$st_ll_tbl <- DT::renderDataTable({
      s <- shiny::req(rv$struct); ll <- s$res$ll; shiny::req(!is.null(ll))
      ll[] <- lapply(ll, function(x) if (is.numeric(x) && !is.integer(x)) signif(x, 6) else x)
      dtt(ll, rownames = FALSE, options = list(dom = "t", scrollX = TRUE))
    })
    st_qpop <- function() { x <- st_cur(); od <- trimws(strsplit(nz(input$st_order, ""), ",")[[1]]); od <- od[nzchar(od)]
                            pops <- unique(x$pop); if (length(od)) pops <- c(od[od %in% pops], setdiff(pops, od))
                            m <- pg_str_bypop(x$Q, x$pop, pops); cbind(Population = rownames(m), n = as.integer(table(factor(x$pop, levels = pops))), as.data.frame(round(m, 4), row.names = NULL)) }
    st_qind <- function() { x <- st_cur(); cbind(Sequence = rownames(x$Q), Population = x$pop, as.data.frame(round(x$Q, 4), row.names = NULL)) }
    output$st_qpop_tbl <- DT::renderDataTable(dtt(st_qpop(), rownames = FALSE, options = list(dom = "t", scrollX = TRUE, pageLength = 50)))
    output$st_q_tbl <- DT::renderDataTable(dtt(st_qind(), rownames = FALSE, options = list(pageLength = 10, lengthMenu = LENGTHS, dom = "lftp", scrollX = TRUE)))
    output$st_dl_q <- shiny::downloadHandler(function() paste0("structure_Q_K", nz(input$st_k, ""), ".csv"), function(file) utils::write.csv(st_qind(), file, row.names = FALSE))
    output$st_dl_qpop <- shiny::downloadHandler(function() paste0("structure_Qpop_K", nz(input$st_k, ""), ".csv"), function(file) utils::write.csv(st_qpop(), file, row.names = FALSE))
    output$st_dl_ll <- shiny::downloadHandler("structure_likelihood.csv", function(file) {
      s <- shiny::req(rv$struct); shiny::req(!is.null(s$res$ll)); utils::write.csv(merge(s$res$runs, s$res$ll[, c("K", "BIC", "deltaK")], by = "K"), file, row.names = FALSE)
    })

    # ---- reports -------------------------------------------------------------------------------------------------------
    report_html <- function() {
      d <- shiny::req(rv$div); a <- rv$amova
      sec <- list(
        list(h = "Run information", pre = c(sprintf("Sequences: %d; sites used (no gaps / N / IUPAC): %d of %d", nrow(d$M), ncol(d$M), d$total),
                                            sprintf("Populations: %s", paste(sort(unique(d$pop)), collapse = ", ")), sprintf("Created by SangeRly on %s", format(Sys.time(), "%Y-%m-%d %H:%M")))),
        list(h = "Samples", tbl = grp_df()),
        list(h = "Molecular diversity indices and neutrality tests", tbl = d$stats, digits = 5,
             p = paste0(sprintf("Sites used: %s%s. ", c(all = "complete in all the sequences", subset = "complete in each population separately", pairwise = "all sites, pairwise deletion for k and pi")[[nz(d$divmode$mode, "all")]],
                          if (is.null(d$divmode$Lden)) "" else sprintf("; per-site values divided by %d variable sites", d$divmode$Lden)),
                   "S segregating sites; h haplotypes; Hd haplotype diversity; k mean pairwise differences (SD_k: Tajima 1993); Pi per site; theta (Hom, S, k, Pi; SD_theta_Pi: Tajima 1983). Tajima's D p: beta distribution; p_D_sim: coalescent simulations as Arlequin, P(D sim < D obs). Fu's Fs p: coalescent simulations.")),
        list(h = "Haplotype distribution", pre = pg_hap_text(d$hap$hap, rownames(d$M))))
      if (!is.null(d$ew)) sec <- c(sec, list(list(h = "Ewens-Watterson neutrality test", tbl = d$ew, digits = 4,
                                                  p = "Infinite alleles model; F = homozygosity; null distribution simulated given the number of haplotypes.")))
      if (!is.null(d$fuli)) sec <- c(sec, list(list(h = "Fu and Li's D* and F*, Ramos-Onsins and Rozas' R2", tbl = d$fuli, digits = 4,
                                                    p = "Biallelic sites only. p_..._le / _ge: fraction of neutral coalescent samples (same n and S) with a value <= / >= the observed one.")))
      if (!is.null(d$exact) && !is.na(d$exact$p)) {
        sec <- c(sec, list(list(h = "Exact test of population differentiation", pre = sprintf("All populations: p = %.5f (s.e. %.5f); %d Markov chain steps.", d$exact$p, d$exact$se, d$exact$steps))))
        if (!is.null(d$exact$pairwise)) { m <- formatC(d$exact$pairwise, digits = 4, format = "f"); m[row(m) == col(m)] <- "-"
                                          sec <- c(sec, list(list(h = "Exact test, p-values for every pair", tbl = cbind(Population = rownames(m), as.data.frame(m, stringsAsFactors = FALSE))))) }
      }
      if (!is.null(d$gf)) sec <- c(sec, list(list(h = "Genetic differentiation and gene flow", pre = dnasp_text(d)),
                                             list(h = "Tests", tbl = d$gf$tests, digits = 5), list(h = "Pairwise comparisons", tbl = d$gf$pairwise, digits = 5)))
      if (!is.null(a)) sec <- c(sec, list(list(h = "AMOVA", tbl = amova_df(a$am), p = sprintf("%d-level AMOVA, %d permutations.", a$am$levels, a$am$nperm)),
                                          list(h = "Fixation indices", tbl = phi_df(a$am)), list(h = "Pairwise FST (below the diagonal) and p-values (above)", tbl = pw_df(a), p = sprintf("+ : p < %s", a$alpha))))
      st <- rv$struct
      if (!is.null(st) && !st$ext) {
        ks <- as.integer(names(st$res$Q)); kk <- nz(input$st_k, max(ks)); kk <- if (as.integer(kk) %in% ks) as.integer(kk) else max(ks)
        sec <- c(sec, list(list(h = "Structure: likelihood by K", tbl = st$res$ll, digits = 4,
                                p = sprintf("%s model; %d variable sites treated as independent; best of %d run(s) per K. delta K: Evanno et al. 2005.", if (st$res$model == "admix") "admixture" else "no-admixture", st$sites, st$res$reps)),
                          list(h = sprintf("Structure: mean ancestry of every population, K = %d", kk),
                               tbl = { Qk <- st$res$Q[[as.character(kk)]]; m <- pg_str_bypop(Qk, cur_pop(rownames(Qk))); cbind(Population = rownames(m), as.data.frame(round(m, 4), row.names = NULL)) })))
      }
      pg_html("SangeRly - population genetics", sec, "Report in the style of the Arlequin results")
    }
    output$rep_tbl <- shiny::downloadHandler("population_genetics_tables.html", function(file) writeLines(report_html(), file, useBytes = TRUE))
    # the result in the layout of the Arlequin result file
    arl_rep <- function() {
      d <- shiny::req(rv$div); a <- rv$amova
      dm <- c(N = "Pairwise difference", K80 = "Kimura 2-parameter", TN93 = "Tamura-Nei", JC69 = "Jukes-Cantor")
      settings <- c(sprintf(" Molecular distance (AMOVA, pairwise FST): %s; gamma a = %s", dm[[nz(input$dmodel, "N")]], format(nz(input$gamma, 0))),
                    sprintf(" No. of permutations: AMOVA = %d; pairwise FST = %d", as.integer(nz(input$nperm, 1000)), as.integer(nz(input$nperm_pair, 200))),
                    sprintf(" Tajima's D and Fu's Fs: no. of simulated samples = %d", as.integer(nz(d$nsim_fs, nz(input$nsim, 1000)))),
                    sprintf(" Ewens-Watterson: no. of simulated samples = %d", as.integer(nz(input$ew_nsim, 1000))),
                    sprintf(" Exact test of sample differentiation: Markov chain steps = %d; dememorization steps = %d", as.integer(nz(input$ex_steps, 0)), as.integer(nz(input$ex_demem, 10000))),
                    sprintf(" Allowed level of missing data (AMOVA, pairwise FST, pi) = %s", format(nz(input$miss, 0))),
                    sprintf(" Significance level = %s", format(nz(input$alpha, 0.05))),
                    " Sites with a gap, N or IUPAC code in any sequence are left out (complete deletion)")
      groups <- NULL
      if (!is.null(a)) { pp <- pop_for(d$M); g <- tapply(pp$pop, pp$grp, function(x) sort(unique(x))); if (length(g) > 1) groups <- g }
      pg_arl_report(d, a, list(title = "SangeRly results", settings = settings, nsim_fs = nz(d$nsim_fs, 1000), alpha = nz(input$alpha, 0.05), groups = groups,
                               nperm_pair = if (is.null(a$nperm_pair)) nz(input$nperm_pair, 200) else a$nperm_pair, dnasp = if (!is.null(d$gf)) dnasp_text(d)))
    }
    output$rep_html <- shiny::downloadHandler("sangerly_results.html", function(file) writeLines(pg_arl_single("SangeRly results", arl_rep()), file, useBytes = TRUE))
    rv_saved <- shiny::reactiveVal(NULL)
    shiny::observeEvent(input$rep_browse, {
      sg_dirpick(session, function(p) shiny::updateTextInput(session, "rep_dir", value = gsub("\\\\", "/", p)),
                 "Folder for the result", nz(input$rep_dir, ""))
    })

    # ---- which steps are done, and which results are outdated ---------------------------------------------------------------
    aln_stamp <- function() { p <- prep(); if (is.null(p)) "0" else paste(length(p$aln), nchar(p$aln[1]), sum(nchar(names(p$aln)))) }
    all_stamp <- function() paste(grp_ver(), aln_stamp())
    shiny::observeEvent(rv$div, rv$st_div <- all_stamp(), ignoreNULL = TRUE)
    shiny::observeEvent(rv$div, rv$slide <- NULL, ignoreNULL = FALSE)
    shiny::observeEvent(rv$amova, rv$st_amova <- all_stamp(), ignoreNULL = TRUE)
    shiny::observeEvent(rv$struct, rv$st_struct <- aln_stamp(), ignoreNULL = TRUE)
    output$steps <- shiny::renderUI({
      p <- prep(); sa <- all_stamp(); sl <- aln_stamp(); d <- rv$pop; np <- if (is.null(d)) 0L else length(unique(d$Population))
      state <- function(x, stamp, now) if (is.null(x)) "todo" else if (!identical(stamp, now)) "stale" else "ok"
      s3 <- state(rv$div, rv$st_div, sa); s4 <- state(rv$amova, rv$st_amova, sa); s5 <- state(rv$struct, rv$st_struct, sl)
      lab <- c(ok = "done", todo = "not run yet", stale = "OUTDATED - run again")
      chip <- function(n, title, st, detail) shiny::div(class = paste0("pg-chip pg-", st), shiny::strong(sprintf("%d. %s", n, title)), shiny::br(), detail)
      how <- c(names = "from the names", file = "from the file", manual = "by hand")[[if (rv$pop_how %in% c("names", "file", "manual")) rv$pop_how else "names"]]
      c1 <- if (is.null(p)) chip(1, "Alignment", "todo", "load an alignment") else chip(1, "Alignment", "ok", sprintf("%d sequences, %d columns", length(p$aln), nchar(p$aln[1])))
      c2 <- if (is.null(d)) chip(2, "Populations", "todo", "waiting for the alignment") else if (np < 2) chip(2, "Populations", "warn", "only one population: AMOVA and FST need two or more") else chip(2, "Populations", "ok", sprintf("%d populations, %s", np, how))
      stale <- c(if (s3 == "stale") "Haplotypes", if (s4 == "stale") "AMOVA", if (s5 == "stale") "Structure")
      shiny::tagList(shiny::div(class = "pg-steps", c1, c2, chip(3, "Haplotypes, diversity", s3, lab[[s3]]), chip(4, "AMOVA, pairwise FST", s4, lab[[s4]]), chip(5, "Structure", s5, lab[[s5]])),
        if (length(stale)) shiny::div(class = "alert alert-warning", style = "padding:6px 10px;font-size:13px",
          sprintf("The results of %s were computed before you changed the alignment or the populations, so the tables and figures below may not match them. Press the buttons in section 4 again.", paste(stale, collapse = ", "))))
    })

    # ---- save everything to one folder ------------------------------------------------------------------------------------------
    save_all <- function() {
      dir <- trimws(nz(input$rep_dir, "")); nm <- gsub("[^A-Za-z0-9._-]", "_", nz(input$rep_name, "sangerly")); if (!nzchar(nm)) nm <- "sangerly"
      if (!nzchar(dir)) { rv_saved(list(type = "warning", txt = "Type or choose a folder.")); return() }
      if (is.null(rv$div) && is.null(rv$amova) && is.null(rv$struct)) { rv_saved(list(type = "warning", txt = "Nothing to save yet: run one of the analyses in section 4 first.")); return() }
      ok <- tryCatch({ dir.create(dir, recursive = TRUE, showWarnings = FALSE); dir.exists(dir) }, error = function(e) FALSE)
      if (!ok) { rv_saved(list(type = "danger", txt = paste("The folder could not be created:", dir))); return() }
      done <- character(); fail <- character()
      put <- function(file, fun) { r <- tryCatch({ fun(file.path(dir, file)); TRUE }, error = function(e) { message("save_all: ", file, ": ", conditionMessage(e)); FALSE }); if (r) done <<- c(done, file) else fail <<- c(fail, file) }
      csv <- function(x) function(f) utils::write.csv(x, f, row.names = FALSE)
      if (!is.null(rv$pop)) put(paste0(nm, "_population_table.csv"), csv(rv$pop))
      d <- rv$div
      if (!is.null(d)) {
        put(paste0(nm, "_haplotypes.csv"), csv(d$haptab)); put(paste0(nm, "_summary.csv"), csv(d$stats))
        if (!is.null(d$fuli)) put(paste0(nm, "_fu_li_R2.csv"), csv(d$fuli))
        if (!is.null(rv$slide)) put(paste0(nm, "_sliding_window.csv"), csv(rv$slide$tab))
        put(paste0(nm, "_report_tables.html"), function(f) writeLines(report_html(), f, useBytes = TRUE))
        res <- tryCatch(pg_arl_save(dir, nm, arl_rep()), error = function(e) NULL)
        if (is.null(res)) fail <- c(fail, paste0(nm, ".res (Arlequin-style folder)")) else done <- c(done, paste0(nm, ".res/ (Arlequin-style: ", basename(res), ")"))
      }
      a <- rv$amova
      if (!is.null(a)) {
        put(paste0(nm, "_amova.csv"), function(f) { ph <- phi_df(a$am)
          out <- rbind(cbind(Table = "AMOVA", amova_df(a$am), P = NA),
                       data.frame(Table = "Fixation indices", Source_of_variation = ph$Fixation_index, d.f. = NA, Sum_of_squares = NA, Variance_components = ph$Value, Percentage_of_variation = NA, P = ph$P, stringsAsFactors = FALSE))
          utils::write.csv(out, f, row.names = FALSE) })
        put(paste0(nm, "_pairwise_fst.csv"), csv(pw_df(a)))
      }
      s5 <- rv$struct
      if (!is.null(s5)) {
        put(paste0(nm, "_structure_Q_individuals.csv"), csv(st_qind())); put(paste0(nm, "_structure_Q_populations.csv"), csv(st_qpop()))
        if (!is.null(s5$res$ll)) put(paste0(nm, "_structure_likelihood.csv"), csv(merge(s5$res$runs, s5$res$ll[, c("K", "BIC", "deltaK")], by = "K")))
      }
      rv$last_out <- list(dir = dir, name = nm, done = done, fail = fail)
      if (isTRUE(rv$net_built)) session$sendCustomMessage("netv_request", list(id = session$ns("net"), input = session$ns("net_files")))
      if (!is.null(rv_coords())) session$sendCustomMessage("hmap_request", list(id = session$ns("hmap"), input = session$ns("hmap_files")))
      rv_saved(list(type = if (length(fail)) "warning" else "success",
                    txt = paste0(sprintf("Saved %d item(s) in %s: ", length(done), dir), paste(done, collapse = "; "),
                                 if (length(fail)) paste0(". Could not write: ", paste(fail, collapse = "; ")) else "", if (isTRUE(rv$net_built)) ". The network (SVG + PNG) follows in a moment." else ". (No network was built, so there is no network file.)")))
      tryCatch(utils::browseURL(dir), error = function(e) NULL, warning = function(w) NULL)
    }
    shiny::observeEvent(input$rep_all, save_all())
    shiny::observeEvent(input$net_save, {
      dir <- trimws(nz(input$rep_dir, "")); nm <- gsub("[^A-Za-z0-9._-]", "_", nz(input$rep_name, "sangerly"))
      if (!isTRUE(rv$net_built)) { set_msg_at("net", "Build the network first.", "warning"); return() }
      if (!nzchar(dir)) { set_msg_at("net", "Choose the output folder in the Result tab first.", "warning"); return() }
      dir.create(dir, recursive = TRUE, showWarnings = FALSE)
      rv$last_out <- list(dir = dir, name = nm, done = character(), fail = character(), net_only = TRUE)
      session$sendCustomMessage("netv_request", list(id = session$ns("net"), input = session$ns("net_files")))
    })
    shiny::observeEvent(input$net_files, {
      x <- input$net_files; o <- rv$last_out; if (is.null(o) || is.null(x$svg)) return()
      f1 <- file.path(o$dir, paste0(o$name, "_haplotype_network.svg")); f2 <- file.path(o$dir, paste0(o$name, "_haplotype_network.png"))
      okS <- tryCatch({ writeLines(x$svg, f1, useBytes = TRUE); TRUE }, error = function(e) FALSE)
      okP <- tryCatch({ writeBin(jsonlite::base64_dec(sub("^data:image/png;base64,", "", x$png)), f2); TRUE }, error = function(e) FALSE)
      if (isTRUE(o$net_only)) { set_msg_at("net", if (okS && okP) sprintf("Saved %s and %s in %s.", basename(f1), basename(f2), o$dir) else "The network could not be written to the folder.", if (okS && okP) "success" else "warning"); tryCatch(utils::browseURL(o$dir), error = function(e) NULL, warning = function(w) NULL) }
      else if (okS) rv_saved(list(type = "success", txt = paste0("Saved in ", o$dir, ": ", paste(c(o$done, basename(f1), if (okP) basename(f2)), collapse = "; "), if (length(o$fail)) paste0(". Could not write: ", paste(o$fail, collapse = "; ")) else "", ".")))
    })
    shiny::observeEvent(input$rep_folder, {
      if (is.null(rv$div)) { rv_saved(list(type = "warning", txt = "Run 'Haplotypes, diversity, differentiation' first.")); return() }
      dir <- trimws(nz(input$rep_dir, "")); if (!nzchar(dir)) { rv_saved(list(type = "warning", txt = "Type or choose a folder.")); return() }
      res <- tryCatch(pg_arl_save(dir, nz(input$rep_name, "sangerly"), arl_rep()), error = function(e) e)
      if (inherits(res, "error")) rv_saved(list(type = "danger", txt = paste("The result could not be saved:", conditionMessage(res))))
      else { rv_saved(list(type = "success", txt = sprintf("Saved. Open this file: %s", res))); tryCatch(utils::browseURL(res), error = function(e) NULL, warning = function(w) NULL) }
    })
    output$rep_saved <- shiny::renderUI({ x <- rv_saved(); if (is.null(x)) NULL else shiny::div(class = paste0("alert alert-", x$type), style = "padding:6px 10px;margin:4px 0;font-size:13px;word-break:break-all", x$txt) })
    output$rep_txt <- shiny::downloadHandler("population_genetics_report.txt", function(file) {
      d <- shiny::req(rv$div); writeLines(c(pg_hap_text(d$hap$hap, rownames(d$M)), "", dnasp_text(d)), file, useBytes = TRUE)
    })

    # ---- files for DnaSP, Arlequin, PopART ---------------------------------------------------------------------------------
    txt_dl <- function(name, fun) shiny::downloadHandler(name, function(file) { shiny::req(rv$div); writeLines(fun(), file, useBytes = TRUE) })
    csv_dl <- function(name, fun) shiny::downloadHandler(name, function(file) { shiny::req(rv$div); utils::write.csv(fun(), file, row.names = FALSE) })
    seqs_of <- function() { M <- rv$div$M; list(seq = apply(M, 1, paste, collapse = ""), nm = rownames(M), M = M) }
    cur_counts <- function() {
      d <- rv$div; pp <- pop_for(d$M); cnt <- table(hap = d$hap$hap, pop = factor(pp$pop))
      matrix(as.integer(cnt), nrow = nrow(cnt), dimnames = list(rownames(cnt), colnames(cnt)))
    }
    group_vec <- function(M) { pp <- pop_for(M); g <- tapply(pp$grp, pp$pop, function(x) x[1]); if (length(unique(g)) > 1) g else NULL }
    # haplotypes with the options of DnaSP's "Haplotype/DNA Sequences Data File"
    hap_file_data <- function() {
      d <- shiny::req(rv$div); pa <- d$aln
      cl <- pg_prepare(pa, if (identical(input$hf_gaps, "yes")) "pairwise" else "excluded"); M <- cl$M
      if (identical(input$hf_invar, "removed")) {
        keep <- vapply(seq_len(ncol(M)), function(j) length(unique(M[, j])) > 1, TRUE)
        if (!any(keep)) stop("Every site is invariable.")
        M <- M[, keep, drop = FALSE]
      }
      hp <- pg_haplotypes(M); pp <- pop_for(M); cnt <- table(hap = hp$hap, pop = factor(pp$pop))
      list(hp = hp, M = M, counts = matrix(as.integer(cnt), nrow = nrow(cnt), dimnames = list(rownames(cnt), colnames(cnt))), g = group_vec(M))
    }
    output$dl_hapfile <- shiny::downloadHandler(
      filename = function() switch(nz(input$hf_format, "nexus"), nexus = "haplotypes.nex", arlequin = "haplotypes.arp", fasta = "haplotypes.fasta"),
      content = function(file) {
        h <- tryCatch(hap_file_data(), error = function(e) { shiny::showNotification(conditionMessage(e), type = "error"); NULL }); shiny::req(h)
        txt <- switch(nz(input$hf_format, "nexus"),
          nexus = pg_nexus_popart(h$hp$seq, h$hp$id, h$counts),
          arlequin = pg_arlequin(h$hp$seq, h$hp$id, h$counts, h$g, "SangeRly"),
          fasta = pg_fasta(h$hp$seq, h$hp$id))
        writeLines(txt, file, useBytes = TRUE)
      })
    output$dl_hap <- csv_dl("haplotypes.csv", function() rv$div$haptab)
    output$dl_sum <- csv_dl("population_genetics_summary.csv", function() rv$div$stats)
    output$ex_fasta <- txt_dl("sequences_used.fasta", function() { s <- seqs_of(); pg_fasta(s$seq, s$nm) })
    output$ex_nexus <- txt_dl("sequences.nex", function() { s <- seqs_of(); pg_nexus(s$seq, pg_safe_names(s$nm), pop_for(s$M)$pop) })
    output$ex_phylip <- txt_dl("sequences.phy", function() { s <- seqs_of(); pg_phylip(s$seq, s$nm)$text })
    output$ex_phymap <- csv_dl("phylip_names.csv", function() { s <- seqs_of(); pg_phylip(s$seq, s$nm)$map })
    output$ex_mega <- txt_dl("sequences.meg", function() { s <- seqs_of(); pg_mega(s$seq, pg_safe_names(s$nm)) })
    output$ex_popart <- txt_dl("haplotypes_popart.nex", function() { d <- rv$div; pg_nexus_popart(d$hap$seq, d$hap$id, cur_counts()) })
    output$ex_pop <- shiny::downloadHandler("population_table.csv", function(file) { shiny::req(rv$pop); utils::write.csv(rv$pop, file, row.names = FALSE) })
    output$dl_amova <- shiny::downloadHandler("amova.csv", function(file) {
      a <- shiny::req(rv$amova); ph <- phi_df(a$am)
      out <- rbind(cbind(Table = "AMOVA", amova_df(a$am), P = NA),
                   data.frame(Table = "Fixation indices", Source_of_variation = ph$Fixation_index, d.f. = NA, Sum_of_squares = NA, Variance_components = ph$Value, Percentage_of_variation = NA, P = ph$P, stringsAsFactors = FALSE))
      utils::write.csv(out, file, row.names = FALSE)
    })
    output$dl_pw <- shiny::downloadHandler("pairwise_fst.csv", function(file) { a <- shiny::req(rv$amova); utils::write.csv(pw_df(a), file, row.names = FALSE) })

    # ---- save / open the sets and settings of this step -----------------------------------------------------------------------
    sids <- c("src", "pop_method", "pop_sep", "pop_nfield", "pop_nchars", "pop_regex", "nsim", "ew_nsim", "ex_steps", "ex_demem", "ex_pair", "div_gap", "div_inv", "gf_gap", "gf_genome", "gf_nperm", "gf_seed", "dmodel", "gamma",
              "nperm", "nperm_pair", "alpha", "slatkin", "reyn", "pipop", "miss", "net_method", "net_eps", "hf_format", "hf_gaps", "hf_invar",
              "st_method", "st_model", "st_kmin", "st_kmax", "st_rep", "st_iter", "st_burn", "st_nrep", "st_alpha", "st_alpha_inf", "st_seed", "st_mode", "st_palette", "st_title", "st_order", "st_colors", "st_cex", "st_sort", "st_names", "st_h", "st_llh", "st_mh", "pop_src", "st_cv", "st_pie", "st_pie_n", "st_pie_lab")
    get_state <- function() list(version = 1L, pop = rv$pop, settings = stats::setNames(lapply(sids, function(i) input[[i]]), sids))
    set_state <- function(s) {
      if (!is.list(s) || is.null(s$pop)) { set_msg_at("run", "This file is not a SangeRly population-genetics file.", "danger"); return(invisible()) }
      rv$pop <- s$pop; rv$pop_how <- "manual"; bump()
      for (i in intersect(names(s$settings), sids)) {
        v <- s$settings[[i]]; if (is.null(v)) next
        if (i == "pop_src") v <- "manual"                    # the saved table is used as it is
        tryCatch({
          if (i %in% c("src", "pop_src", "gf_gap", "net_method", "hf_format", "hf_gaps", "hf_invar")) shiny::updateRadioButtons(session, i, selected = v)
          else if (i %in% c("pop_method", "gf_genome", "dmodel", "st_method", "st_model", "st_mode", "st_palette")) shiny::updateSelectInput(session, i, selected = v)
          else if (i %in% c("pop_sep", "pop_regex", "st_title", "st_order", "st_colors")) shiny::updateTextInput(session, i, value = v)
          else if (i %in% c("slatkin", "reyn", "pipop", "ex_pair", "st_sort", "st_names", "st_pie_n", "st_pie_lab", "st_alpha_inf")) shiny::updateCheckboxInput(session, i, value = isTRUE(v))
          else if (i %in% c("st_cex", "st_h", "st_pie", "st_llh", "st_mh")) shiny::updateSliderInput(session, i, value = v)
          else shiny::updateNumericInput(session, i, value = v)
        }, error = function(e) NULL)
      }
      set_msg_at("run", "Sets and settings restored.", "success")
    }
    output$proj_save <- shiny::downloadHandler("sangerly_popgen.rds", function(file) saveRDS(get_state(), file))
    shiny::observeEvent(input$proj_load, { s <- tryCatch(readRDS(input$proj_load$datapath), error = function(e) NULL); set_state(s) })

    # for other tabs (Barcode gap): the population table (Sequence, Population, Group) and the ancestry cluster of every sequence
    struct_groups <- shiny::reactive({
      x <- tryCatch(st_cur(), error = function(e) NULL)
      if (is.null(x) || is.null(rownames(x$Q)) || ncol(x$Q) < 2) return(NULL)
      stats::setNames(sprintf("Cluster_%d", max.col(x$Q, ties.method = "first")), rownames(x$Q))
    })
    list(div = shiny::reactive(rv$div), amova = shiny::reactive(rv$amova), struct = shiny::reactive(rv$struct),
         pop = shiny::reactive(rv$pop), struct_groups = struct_groups, state = get_state, restore = set_state)
  })
}
