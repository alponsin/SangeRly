# ==============================================================================
# SangeRly - BLAST module (R Shiny)
# ------------------------------------------------------------------------------
# Runs NCBI BLAST (blastn) through the BLAST URL API with the same options as
# the website: database, program (megablast / discontiguous / blastn),
# organism include & exclude (via Entrez query), e-value, word size,
# low-complexity filter, identity cutoff and number of hits.
#
# Required packages: shiny, httr, xml2, jsonlite, DT
#   install.packages(c("shiny", "httr", "xml2", "jsonlite", "DT"))
#
# Use inside the SangeRly app:
#   source("R/blast_module.R")
#   ui:     blastUI("blast")
#   server: blast_out <- blastServer("blast", consensus = reactive(consensus_seq()),
#                                    email = "you@example.org")
#
# Use without Shiny (plain script):
#   res <- blast_run_sync(">q\nACGT...", db = "core_nt",
#                         include = "Mammalia", exclude = "Homo sapiens")
#
# NCBI rules: no more than 1 request per 10 seconds, poll an RID no more than
# once per minute, and include a contact email (parameter `email`). Results
# are kept on the NCBI server for about 24 hours.
# ==============================================================================

BLAST_URL <- "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
TOOL_NAME <- "Sangerly"

# Nucleotide databases commonly used on the BLAST website ---------------------
blast_databases <- c(
  "Core nucleotide (core_nt)"                         = "core_nt",
  "Nucleotide collection (nt)"                        = "nt",
  "RefSeq Select RNA (refseq_select_rna)"             = "refseq_select_rna",
  "RefSeq RNA (refseq_rna)"                           = "refseq_rna",
  "RefSeq representative genomes"                     = "refseq_representative_genomes",
  "RefSeq genomes"                                    = "refseq_genomes",
  "Whole-genome shotgun contigs (wgs)"                = "wgs",
  "Transcriptome shotgun assembly (tsa_nt)"           = "tsa_nt",
  "16S ribosomal RNA (Bacteria & Archaea)"            = "16S_ribosomal_RNA",
  "18S ribosomal RNA (SSU eukaryote)"                 = "SSU_eukaryote_rRNA",
  "28S ribosomal RNA (LSU eukaryote)"                 = "LSU_eukaryote_rRNA",
  "ITS RefSeq Fungi"                                  = "ITS_RefSeq_Fungi",
  "ITS eukaryote sequences"                           = "ITS_eukaryote_sequences",
  "Environmental samples (env_nt)"                    = "env_nt"
)

blast_programs <- c(
  "Highly similar sequences (megablast)"              = "megablast",
  "More dissimilar sequences (discontiguous megablast)" = "dc-megablast",
  "Somewhat similar sequences (blastn)"               = "blastn"
)

# ------------------------------------------------------------------------------
# Input helpers
# ------------------------------------------------------------------------------

# Take the first record from FASTA text / a raw sequence, clean and validate it.
blast_clean_query <- function(txt) {
  lines <- trimws(unlist(strsplit(as.character(txt), "\r?\n")))
  lines <- lines[nzchar(lines)]
  if (!length(lines)) stop("The query sequence is empty.")
  hdr_idx <- grep("^>", lines)
  header <- "query"
  if (length(hdr_idx)) {
    header <- sub("^>\\s*", "", lines[hdr_idx[1]])
    end <- if (length(hdr_idx) > 1) hdr_idx[2] - 1 else length(lines)
    lines <- lines[(hdr_idx[1] + 1):end]
  }
  seq <- toupper(gsub("[^A-Za-z-]", "", paste(lines, collapse = "")))
  seq <- gsub("-", "", seq)
  if (nchar(seq) < 25) stop("The sequence is too short (minimum 25 nt).")
  bad <- gsub("[ACGTURYSWKMBDHVN]", "", seq)
  if (nzchar(bad)) {
    stop("Invalid characters in the sequence: ",
         paste(unique(strsplit(bad, "")[[1]]), collapse = ", "))
  }
  list(header = header, seq = gsub("U", "T", seq))
}

# Check / clean ONE sequence; `label` is used in the error message.
.blast_check_seq <- function(seq, label = "query") {
  seq <- toupper(gsub("[^A-Za-z-]", "", seq))
  seq <- gsub("-", "", seq)
  if (nchar(seq) < 25) stop(label, ": the sequence is too short (minimum 25 nt).", call. = FALSE)
  bad <- gsub("[ACGTURYSWKMBDHVN]", "", seq)
  if (nzchar(bad)) {
    stop(label, ": invalid characters in the sequence: ",
         paste(unique(strsplit(bad, "")[[1]]), collapse = ", "), call. = FALSE)
  }
  gsub("U", "T", seq)
}

# Split FASTA text (or one raw sequence) into ALL its records: data.frame(id, header, seq).
blast_parse_records <- function(txt) {
  none <- data.frame(id = character(), header = character(), seq = character(),
                     stringsAsFactors = FALSE)
  lines <- trimws(unlist(strsplit(as.character(txt), "\r?\n")))
  lines <- lines[nzchar(lines)]
  if (!length(lines)) return(none)
  hdr_idx <- grep("^>", lines)
  if (!length(hdr_idx)) {
    return(data.frame(id = "query", header = "query",
                      seq = toupper(gsub("[^A-Za-z-]", "", paste(lines, collapse = ""))),
                      stringsAsFactors = FALSE))
  }
  ends <- c(hdr_idx[-1] - 1, length(lines))
  hdr <- sub("^>\\s*", "", lines[hdr_idx])
  seqs <- mapply(function(a, b) if (b > a) paste(lines[(a + 1):b], collapse = "") else "",
                 hdr_idx, ends, USE.NAMES = FALSE)
  ids <- vapply(seq_along(hdr), function(i) {
    w <- strsplit(hdr[i], "\\s+")[[1]][1]
    w <- if (is.na(w)) "" else gsub("[^A-Za-z0-9._-]", "_", w)
    if (!nzchar(w)) paste0("seq", i) else w
  }, character(1))
  data.frame(id = make.unique(ids, sep = "_"), header = hdr,
             seq = toupper(gsub("[^A-Za-z-]", "", seqs)), stringsAsFactors = FALSE)
}

# Split organism input (one per line / comma / semicolon separated)
.split_terms <- function(x) {
  if (is.null(x) || !nzchar(trimws(x))) return(character())
  t <- trimws(unlist(strsplit(x, "[\n,;]+")))
  unique(t[nzchar(t)])
}

# Turn one taxon term into Entrez syntax.
#  - "9606" atau "txid9606"  -> txid9606[Organism]
#  - "Panthera leo"          -> "Panthera leo"[Organism]
# An [Organism] search automatically includes all child taxa (explode),
# like the "Organism" box on the BLAST website.
.organism_term <- function(t) {
  t <- trimws(t)
  if (grepl("^(txid)?[0-9]+$", t, ignore.case = TRUE)) {
    sprintf("txid%s[Organism]", sub("^txid", "", t, ignore.case = TRUE))
  } else {
    sprintf('"%s"[Organism]', t)
  }
}

# Build the Entrez query from the include / exclude lists.
blast_build_entrez <- function(include = character(),
                               exclude = character(),
                               drop_env = FALSE,
                               extra = "") {
  inc <- vapply(include, .organism_term, character(1), USE.NAMES = FALSE)
  exc <- vapply(exclude, .organism_term, character(1), USE.NAMES = FALSE)
  if (drop_env) {
    exc <- c(exc, "environmental samples[filter]", "metagenomes[orgn]")
  }
  inc_part <- if (length(inc)) paste0("(", paste(inc, collapse = " OR "), ")") else "all[filter]"
  q <- inc_part
  if (length(exc)) q <- paste0(q, " NOT (", paste(exc, collapse = " OR "), ")")
  extra <- trimws(extra %||% "")
  if (nzchar(extra)) q <- paste0("(", q, ") AND (", extra, ")")
  if (identical(q, "all[filter]")) "" else q
}

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || (is.character(a) && !nzchar(a[1]))) b else a

# Validate taxon names against NCBI Taxonomy (like the website autocomplete).
blast_resolve_taxa <- function(terms, email = "") {
  if (!length(terms)) return(data.frame())
  out <- lapply(terms, function(t) {
    Sys.sleep(0.4)  # respect the limit of 3 requests/second without an API key
    if (grepl("^(txid)?[0-9]+$", t, ignore.case = TRUE)) {
      id <- sub("^txid", "", t, ignore.case = TRUE)
    } else {
      r <- httr::GET("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
                     query = list(db = "taxonomy", term = paste0('"', t, '"[All Names]'),
                                  retmode = "json", tool = TOOL_NAME, email = email),
                     httr::timeout(30))
      j <- jsonlite::fromJSON(httr::content(r, as = "text", encoding = "UTF-8"))
      ids <- j$esearchresult$idlist
      if (!length(ids)) {
        return(data.frame(Input = t, TaxID = NA, Scientific_name = "NOT FOUND",
                          Rank = NA, stringsAsFactors = FALSE))
      }
      id <- ids[1]
    }
    Sys.sleep(0.4)
    r2 <- httr::GET("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi",
                    query = list(db = "taxonomy", id = id, retmode = "json",
                                 tool = TOOL_NAME, email = email),
                    httr::timeout(30))
    j2 <- jsonlite::fromJSON(httr::content(r2, as = "text", encoding = "UTF-8"))
    rec <- j2$result[[id]]
    data.frame(Input = t, TaxID = id,
               Scientific_name = rec$scientificname %||% NA,
               Rank = rec$rank %||% NA, stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}

# ------------------------------------------------------------------------------
# Communication with the NCBI BLAST URL API
# ------------------------------------------------------------------------------

blast_submit <- function(seq, db = "core_nt", program = "megablast",
                         entrez = "", evalue = 0.05, hitlist = 50,
                         word_size = NA, low_complexity = TRUE,
                         perc_ident = 0, email = "") {
  body <- list(
    CMD = "Put", PROGRAM = "blastn", DATABASE = db, QUERY = seq,
    EXPECT = evalue, HITLIST_SIZE = hitlist,
    FILTER = if (isTRUE(low_complexity)) "L" else "F",
    tool = TOOL_NAME
  )
  if (nzchar(email)) body$email <- email
  if (nzchar(entrez)) body$ENTREZ_QUERY <- entrez
  if (!is.na(perc_ident) && perc_ident > 0) body$PERC_IDENT <- perc_ident
  if (program == "megablast") {
    body$MEGABLAST <- "on"
  } else if (program == "dc-megablast") {
    body$MEGABLAST <- "on"; body$TEMPLATE_TYPE <- "0"; body$TEMPLATE_LENGTH <- "18"
  }
  if (!is.na(word_size)) body$WORD_SIZE <- word_size

  r <- httr::POST(BLAST_URL, body = body, encode = "form", httr::timeout(90),
                  httr::user_agent(paste0(TOOL_NAME, "/0.1")))
  httr::stop_for_status(r)
  txt <- httr::content(r, as = "text", encoding = "UTF-8")
  rid  <- regmatches(txt, regexpr("(?<=RID = )\\S+", txt, perl = TRUE))
  rtoe <- regmatches(txt, regexpr("(?<=RTOE = )[0-9]+", txt, perl = TRUE))
  if (!length(rid)) {
    stop("NCBI did not return an RID. Check the Entrez query / database. ",
         "Response excerpt: ", substr(gsub("\\s+", " ", gsub("<[^>]+>", " ", txt)), 1, 300))
  }
  list(rid = rid, rtoe = if (length(rtoe)) as.numeric(rtoe) else 15)
}

blast_status <- function(rid) {
  r <- httr::GET(BLAST_URL,
                 query = list(CMD = "Get", FORMAT_OBJECT = "SearchInfo", RID = rid),
                 httr::timeout(60), httr::user_agent(paste0(TOOL_NAME, "/0.1")))
  httr::stop_for_status(r)
  txt <- httr::content(r, as = "text", encoding = "UTF-8")
  st <- regmatches(txt, regexpr("(?<=Status=)\\w+", txt, perl = TRUE))
  list(status = if (length(st)) st else "UNKNOWN",
       hits = grepl("ThereAreHits=yes", txt, fixed = TRUE))
}

blast_fetch_xml <- function(rid) {
  r <- httr::GET(BLAST_URL,
                 query = list(CMD = "Get", RID = rid, FORMAT_TYPE = "XML",
                              ALIGNMENTS = 500, DESCRIPTIONS = 500),
                 httr::timeout(180), httr::user_agent(paste0(TOOL_NAME, "/0.1")))
  httr::stop_for_status(r)
  httr::content(r, as = "text", encoding = "UTF-8")
}

# Species name taken from the first two words of the hit definition (an estimate!).
.species_from_def <- function(def) {
  vapply(def, function(d) {
    d <- sub("^(PREDICTED|UNVERIFIED|MAG|TPA_asm|TPA)[:_ ]*\\s*", "", d, ignore.case = TRUE)
    w <- strsplit(trimws(d), "\\s+")[[1]]
    if (length(w) >= 2) paste(w[1], w[2]) else d
  }, character(1), USE.NAMES = FALSE)
}

BLAST_MAX_FULL <- 2e6   # FASTA complete sequence: hits longer than this (bp) are skipped

.empty_blast_xml <- paste0("<BlastOutput><BlastOutput_iterations><Iteration>",
  "<Iteration_query-len>1</Iteration_query-len><Iteration_hits/></Iteration>",
  "</BlastOutput_iterations></BlastOutput>")

blast_parse_xml <- function(xml_txt, query_names = NULL) {
  empty <- data.frame(Query = character(), Rank = integer(), Species = character(),
                      Accession = character(),
                      Description = character(), Identity = numeric(),
                      QueryCover = numeric(), Evalue = numeric(),
                      BitScore = numeric(), TotalScore = numeric(), AlignLen = integer(),
                      Length = numeric(), HitFrom = numeric(), HitTo = numeric(),
                      Strand = character(), Hseq = character(),
                      stringsAsFactors = FALSE)
  doc <- xml2::read_xml(xml_txt)
  iters <- xml2::xml_find_all(doc, "//Iteration")
  if (!length(iters)) return(empty)

  num <- function(node, tag) as.numeric(xml2::xml_text(xml2::xml_find_all(node, paste0(".//", tag))))
  txt1 <- function(node, tag) xml2::xml_text(xml2::xml_find_first(node, paste0("./", tag)))

  per_query <- lapply(seq_along(iters), function(k) {
    it <- iters[[k]]
    qlen <- as.numeric(xml2::xml_text(xml2::xml_find_first(it, "./Iteration_query-len")))
    hits <- xml2::xml_find_all(it, "./Iteration_hits/Hit")
    if (!length(hits)) return(NULL)
    qname <- if (!is.null(query_names) && length(query_names) >= k) query_names[k] else
      if (length(iters) > 1) paste0("query", k) else "query"

    rows <- lapply(seq_along(hits), function(i) {
      h <- hits[[i]]
      qf <- num(h, "Hsp_query-from"); qt <- num(h, "Hsp_query-to")
      ident <- num(h, "Hsp_identity"); alen <- num(h, "Hsp_align-len")
      bits <- num(h, "Hsp_bit-score"); ev <- num(h, "Hsp_evalue")
      b <- which.max(bits)
      hf <- num(h, "Hsp_hit-from"); ht <- num(h, "Hsp_hit-to")
      qfr <- num(h, "Hsp_query-frame"); hfr <- num(h, "Hsp_hit-frame")
      hseq <- xml2::xml_text(xml2::xml_find_all(h, ".//Hsp_hseq"))
      strand <- paste0(if (length(qfr) >= b && !is.na(qfr[b]) && qfr[b] < 0) "Minus" else "Plus", "/",
                       if (length(hfr) >= b && !is.na(hfr[b]) && hfr[b] < 0) "Minus" else "Plus")
      covered <- unique(unlist(Map(function(a, z) seq.int(min(a, z), max(a, z)), qf, qt)))
      def <- txt1(h, "Hit_def")
      data.frame(
        Query = qname,
        Rank = i,
        Species = .species_from_def(def),
        Accession = txt1(h, "Hit_accession"),
        Description = def,
        Identity = round(100 * ident[b] / alen[b], 2),
        QueryCover = round(100 * length(covered) / qlen, 1),
        Evalue = ev[b],
        BitScore = round(bits[b], 1),
        TotalScore = round(sum(bits), 1),
        AlignLen = as.integer(alen[b]),
        Length = suppressWarnings(as.numeric(txt1(h, "Hit_len"))),
        HitFrom = hf[b], HitTo = ht[b], Strand = strand,
        Hseq = if (length(hseq) >= b) hseq[b] else "",
        stringsAsFactors = FALSE
      )
    })
    do.call(rbind, rows)
  })
  per_query <- per_query[!vapply(per_query, is.null, logical(1))]
  if (!length(per_query)) return(empty)
  do.call(rbind, per_query)
}

# All HSPs (alignment segments) of all hits: used by Graphic Summary and Alignments.
blast_parse_hsps <- function(xml_txt, query_names = NULL) {
  doc <- xml2::read_xml(xml_txt)
  iters <- xml2::xml_find_all(doc, "//Iteration")
  out <- list()
  for (k in seq_along(iters)) {
    it <- iters[[k]]
    qlen <- as.numeric(xml2::xml_text(xml2::xml_find_first(it, "./Iteration_query-len")))
    qname <- if (!is.null(query_names) && length(query_names) >= k) query_names[k] else
      if (length(iters) > 1) paste0("query", k) else "query"
    hits <- xml2::xml_find_all(it, "./Iteration_hits/Hit")
    for (i in seq_along(hits)) {
      h <- hits[[i]]
      hsps <- xml2::xml_find_all(h, ".//Hsp")
      if (!length(hsps)) next
      g <- function(tag) xml2::xml_text(xml2::xml_find_first(hsps, paste0("./", tag)))
      out[[length(out) + 1]] <- data.frame(
        Query = qname, QLen = qlen, Rank = i,
        Accession = xml2::xml_text(xml2::xml_find_first(h, "./Hit_accession")),
        Description = xml2::xml_text(xml2::xml_find_first(h, "./Hit_def")),
        QFrom = as.numeric(g("Hsp_query-from")), QTo = as.numeric(g("Hsp_query-to")),
        HFrom = as.numeric(g("Hsp_hit-from")), HTo = as.numeric(g("Hsp_hit-to")),
        Bits = as.numeric(g("Hsp_bit-score")), Score = as.numeric(g("Hsp_score")),
        Evalue = as.numeric(g("Hsp_evalue")),
        Identity = as.numeric(g("Hsp_identity")), Gaps = as.numeric(g("Hsp_gaps")),
        AlignLen = as.numeric(g("Hsp_align-len")),
        QFrame = as.numeric(g("Hsp_query-frame")), HFrame = as.numeric(g("Hsp_hit-frame")),
        Qseq = g("Hsp_qseq"), Hseq = g("Hsp_hseq"), Midline = g("Hsp_midline"),
        stringsAsFactors = FALSE)
    }
  }
  if (!length(out)) return(NULL)
  do.call(rbind, out)
}

# Text alignment in the NCBI style (Query / midline / Sbjct, 60 columns, with coordinates).
blast_format_alignment <- function(q, h, mid, qfrom, hfrom, qminus = FALSE, hminus = FALSE, width = 60) {
  n <- nchar(q)
  if (is.na(n) || n == 0) return("")
  qs <- strsplit(q, "")[[1]]; hs_ <- strsplit(h, "")[[1]]
  ms <- if (is.na(mid) || !nzchar(mid)) rep(" ", n) else strsplit(mid, "")[[1]]
  if (length(ms) < n) ms <- c(ms, rep(" ", n - length(ms)))
  qstep <- if (isTRUE(qminus)) -1 else 1
  hstep <- if (isTRUE(hminus)) -1 else 1
  fmt <- function(x) formatC(x, format = "d", big.mark = "")
  qpos <- qfrom; hpos <- hfrom
  w <- nchar(fmt(max(abs(c(qfrom, hfrom)), na.rm = TRUE) + n))
  out <- character()
  for (st in seq(1, n, by = width)) {
    en <- min(n, st + width - 1)
    qseg <- qs[st:en]; hseg <- hs_[st:en]; mseg <- ms[st:en]
    nq <- sum(qseg != "-"); nh <- sum(hseg != "-")
    out <- c(out,
      sprintf("Query  %*s  %s  %s", w, fmt(qpos), paste(qseg, collapse = ""), fmt(qpos + qstep * (nq - 1))),
      sprintf("       %*s  %s", w, "", paste(mseg, collapse = "")),
      sprintf("Sbjct  %*s  %s  %s", w, fmt(hpos), paste(hseg, collapse = ""), fmt(hpos + hstep * (nh - 1))),
      "")
    qpos <- qpos + qstep * nq; hpos <- hpos + hstep * nh
  }
  paste(out, collapse = "\n")
}

# Lineage (superkingdom ... genus) of species names through NCBI Taxonomy (needs internet).
blast_taxonomy_lineage <- function(species, email = "") {
  species <- unique(species[!is.na(species) & nzchar(species)])
  ranks <- c("superkingdom", "kingdom", "phylum", "class", "order", "family", "genus")
  ids <- vapply(species, function(sp) {
    Sys.sleep(0.4)   # NCBI: at most 3 requests/second without an API key
    r <- tryCatch(httr::GET("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
                            query = list(db = "taxonomy", term = paste0('"', sp, '"[All Names]'),
                                         retmode = "json", tool = TOOL_NAME, email = email),
                            httr::timeout(30)), error = function(e) NULL)
    if (is.null(r)) return(NA_character_)
    j <- tryCatch(jsonlite::fromJSON(httr::content(r, as = "text", encoding = "UTF-8")),
                  error = function(e) NULL)
    idl <- j$esearchresult$idlist
    if (length(idl)) idl[1] else NA_character_
  }, character(1), USE.NAMES = FALSE)
  res <- data.frame(Species = species, TaxID = ids, stringsAsFactors = FALSE)
  for (rk in ranks) res[[rk]] <- NA_character_
  ok <- !is.na(ids)
  if (!any(ok)) return(res)
  r <- httr::POST("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
                  body = list(db = "taxonomy", id = paste(unique(ids[ok]), collapse = ","),
                              retmode = "xml", tool = TOOL_NAME, email = email),
                  encode = "form", httr::timeout(120))
  httr::stop_for_status(r)
  doc <- xml2::read_xml(httr::content(r, as = "text", encoding = "UTF-8"))
  tx_all <- xml2::xml_find_all(doc, "/TaxaSet/Taxon")
  for (ii in seq_along(tx_all)) {
    tx <- tx_all[[ii]]
    id <- xml2::xml_text(xml2::xml_find_first(tx, "./TaxId"))
    lin <- xml2::xml_find_all(tx, "./LineageEx/Taxon")
    nm <- xml2::xml_text(xml2::xml_find_first(lin, "./ScientificName"))
    rkn <- xml2::xml_text(xml2::xml_find_first(lin, "./Rank"))
    rows <- which(res$TaxID == id)
    for (rk in ranks) {
      v <- nm[rkn == rk]
      if (length(v) && length(rows)) res[rows, rk] <- v[1]
    }
  }
  res
}

# ---- run on the NCBI website instead of inside the app ----------------------------------------
# Parameters of the web form (pre-filled; the user presses BLAST on the website).
blast_web_params <- function(fasta, db, program, entrez = "", evalue = 0.05, hitlist = 50,
                             word_size = NA, low_complexity = TRUE) {
  p <- list(PAGE_TYPE = "BlastSearch", PROGRAM = "blastn",
            BLAST_PROGRAMS = switch(program, "dc-megablast" = "discoMegablast", blastn = "blastn", "megaBlast"),
            DATABASE = db, QUERY = fasta, LINK_LOC = "blasthome",
            MAX_NUM_SEQ = hitlist, EXPECT = evalue)
  if (isTRUE(low_complexity)) p$FILTER <- "L"
  if (!is.na(word_size)) p$WORD_SIZE <- word_size
  if (nzchar(entrez)) p$EQ_TEXT <- entrez
  p
}

# Open the web form in the browser. On the desktop / local R the browser is opened by R (through a
# small local page that sends the form, so long sequences are fine); on a web server the browser of
# the visitor opens it (JavaScript).
blast_open_site <- function(params, session = NULL) {
  params <- lapply(params, as.character)
  if (nzchar(Sys.getenv("SANGERLY_DESKTOP")) || interactive() || is.null(session)) {
    esc <- function(x) gsub(">", "&gt;", gsub("<", "&lt;", gsub("&", "&amp;", x, fixed = TRUE), fixed = TRUE), fixed = TRUE)
    fields <- paste0('<textarea name="', names(params), '" style="display:none">',
                     vapply(params, esc, character(1)), "</textarea>", collapse = "\n")
    html <- paste0('<!DOCTYPE html><html><head><meta charset="utf-8"><title>Opening NCBI BLAST</title></head><body>',
                   "<p>Opening NCBI BLAST... if nothing happens, press the button.</p>",
                   '<form id="f" method="get" action="', BLAST_URL, '">', fields,
                   '<button type="submit">Open NCBI BLAST</button></form>',
                   '<script>document.getElementById("f").submit();</script></body></html>')
    f <- tempfile("sangerly_blast_", fileext = ".html")
    writeLines(html, f, useBytes = TRUE)
    # Dipanggil lewat get() agar pemeriksa RStudio tidak menandainya saat publish.
    # Cabang ini hanya jalan di desktop / R lokal, bukan di server.
    open_fn <- get("browseURL", envir = asNamespace("utils"))
    open_fn(f)
  } else {
    session$sendCustomMessage("sangerly_open_blast", list(action = BLAST_URL, params = params))
  }
  invisible(TRUE)
}

# FASTA "aligned sequences": the part of the subject aligned to the query (gaps removed), like
# Download > FASTA (aligned sequences) on the BLAST website.
blast_fasta_aligned <- function(df) {
  df <- df[!is.na(df$Hseq) & nzchar(gsub("-", "", df$Hseq)), , drop = FALSE]
  if (!nrow(df)) return("")
  seqs <- gsub("-", "", df$Hseq)
  wrap <- function(s) paste(substring(s, seq(1, nchar(s), by = 70), pmin(nchar(s), seq(70, nchar(s) + 69, by = 70))),
                            collapse = "\n")
  hdr <- sprintf(">%s %s [%d-%d]", df$Accession, df$Description, as.integer(df$HitFrom), as.integer(df$HitTo))
  paste0(hdr, "\n", vapply(seqs, wrap, "", USE.NAMES = FALSE), collapse = "\n")
}

# FASTA "complete sequence": full sequence of every accession via NCBI efetch
blast_fetch_fasta_complete <- function(accs, email = "", chunk = 50) {
  accs <- unique(accs[!is.na(accs) & nzchar(accs)])
  if (!length(accs)) return("")
  out <- character()
  for (i in seq(1, length(accs), by = chunk)) {
    ids <- accs[i:min(length(accs), i + chunk - 1)]
    r <- httr::POST("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
                    body = list(db = "nuccore", id = paste(ids, collapse = ","), rettype = "fasta",
                                retmode = "text", tool = TOOL_NAME, email = email),
                    encode = "form", httr::timeout(300))
    httr::stop_for_status(r)
    out <- c(out, httr::content(r, as = "text", encoding = "UTF-8"))
    Sys.sleep(0.5)   # NCBI limit without an API key: 3 requests/second
  }
  paste(out, collapse = "")
}

# Synchronous version (for scripts/tests, no Shiny). Blocks until finished.
blast_run_sync <- function(query, db = "core_nt", program = "megablast",
                           include = character(), exclude = character(),
                           drop_env = FALSE, extra_entrez = "",
                           evalue = 0.05, hitlist = 50, word_size = NA,
                           low_complexity = TRUE, perc_ident = 0,
                           email = "", timeout_min = 20, verbose = TRUE) {
  q <- blast_clean_query(query)
  entrez <- blast_build_entrez(include, exclude, drop_env, extra_entrez)
  if (verbose) message("Entrez query: ", if (nzchar(entrez)) entrez else "(none)")
  sub <- blast_submit(q$seq, db, program, entrez, evalue, hitlist, word_size,
                      low_complexity, perc_ident, email)
  if (verbose) message("RID: ", sub$rid, " (estimated ", sub$rtoe, " s)")
  Sys.sleep(max(sub$rtoe, 10))
  t0 <- Sys.time()
  repeat {
    st <- blast_status(sub$rid)
    if (st$status == "READY") break
    if (st$status %in% c("FAILED", "UNKNOWN")) stop("BLAST failed / RID expired: ", st$status)
    if (as.numeric(difftime(Sys.time(), t0, units = "mins")) > timeout_min)
      stop("Timed out after ", timeout_min, " minutes. RID: ", sub$rid)
    if (verbose) message("Waiting...")
    Sys.sleep(60)
  }
  if (!st$hits) return(list(rid = sub$rid, entrez = entrez, hits = blast_parse_xml(
    "<BlastOutput><BlastOutput_iterations><Iteration><Iteration_query-len>1</Iteration_query-len><Iteration_hits/></Iteration></BlastOutput_iterations></BlastOutput>")))
  list(rid = sub$rid, entrez = entrez, hits = blast_parse_xml(blast_fetch_xml(sub$rid)))
}

# ------------------------------------------------------------------------------
# Interpretation of results (simple rules, NOT automatic species assignment)
# ------------------------------------------------------------------------------
blast_interpret <- function(df, delta = 1) {
  if (is.null(df) || !nrow(df)) {
    return("No hits passed the filters. Relax the e-value, the identity cutoff or the organism restriction.")
  }
  d <- df[order(-df$Identity, -df$QueryCover, -df$BitScore), , drop = FALSE]
  best <- d[1, ]
  per_sp <- d[!duplicated(d$Species), , drop = FALSE]
  close <- per_sp[per_sp$Species != best$Species &
                    per_sp$Identity >= best$Identity - delta, , drop = FALSE]
  msg <- sprintf("Best match: %s (%s) - identity %.2f%%, query cover %.1f%%.",
                 best$Species, best$Accession, best$Identity, best$QueryCover)
  if (nrow(close)) {
    msg <- c(msg, sprintf(
      "WARNING: %d other species are within %.1f%% of the top identity (%s). Species-level identification may not be certain.",
      nrow(close), delta, paste(utils::head(close$Species, 5), collapse = ", ")))
  } else {
    msg <- c(msg, "No other species is very close in the hits shown.")
  }
  if (best$QueryCover < 90) {
    msg <- c(msg, "Query cover < 90%: the hit may cover only part of the sequence (check trimming / possible chimera).")
  }
  c(msg, paste0("Note: species names are taken from the first two words of the hit description and can be wrong ",
                "(e.g. 'uncultured', misannotation). Identity thresholds for species boundaries differ between taxa ",
                "and genes; also check geographic distribution, reference taxonomy and database record quality."))
}

# ------------------------------------------------------------------------------
# Shiny module
# ------------------------------------------------------------------------------
blastUI <- function(id) {
  ns <- shiny::NS(id)
  shiny::fluidRow(
    shiny::column(
      5,
      shiny::tags$script(shiny::HTML(paste0(
        "(function(){function reg(){Shiny.addCustomMessageHandler('sangerly_open_blast',function(m){",
        "var f=document.createElement('form');f.method='GET';f.action=m.action;f.target='_blank';",
        "f.style.display='none';Object.keys(m.params).forEach(function(k){",
        "var t=document.createElement('textarea');t.name=k;t.value=m.params[k];f.appendChild(t);});",
        "document.body.appendChild(f);f.submit();document.body.removeChild(f);});}",
        "if(window.Shiny&&Shiny.addCustomMessageHandler){reg();}else{$(document).on('shiny:connected',reg);}})();"))),
      shiny::h4("Query sequence"),
      shiny::textAreaInput(ns("query"), NULL, rows = 6, width = "100%",
                           placeholder = "Paste FASTA / sequence, or use the button below to take the consensus"),
      shiny::actionButton(ns("use_consensus"), "Use consensus from the assembly"),
      shiny::fileInput(ns("fasta_file"), NULL, accept = c(".fa", ".fasta", ".fas", ".txt"),
                       buttonLabel = "Upload FASTA", placeholder = "(optional)"),
      shiny::uiOutput(ns("seq_select")),
      shiny::h4("Database & program"),
      shiny::selectInput(ns("db"), "Database", choices = blast_databases, selected = "core_nt"),
      shiny::selectInput(ns("program"), "Program", choices = blast_programs),
      shiny::h4("Restrict organisms"),
      shiny::textAreaInput(ns("include"), "Include (one per line: scientific name or taxid)",
                           rows = 3, width = "100%", placeholder = "Mammalia\n9606"),
      shiny::textAreaInput(ns("exclude"), "Exclude (one per line)",
                           rows = 3, width = "100%", placeholder = "Homo sapiens"),
      shiny::checkboxInput(ns("drop_env"), "Exclude environmental samples / metagenome (uncultured)", FALSE),
      shiny::actionButton(ns("check_taxa"), "Validate taxon names"),
      shiny::tableOutput(ns("taxa_table")),
      shiny::h4("Advanced parameters"),
      shiny::fluidRow(
        shiny::column(6, shiny::numericInput(ns("hitlist"), "Max target sequences", 50, min = 5, max = 500, step = 5)),
        shiny::column(6, shiny::numericInput(ns("evalue"), "E-value", 0.05, min = 0, step = 0.01))
      ),
      shiny::fluidRow(
        shiny::column(6, shiny::numericInput(ns("word_size"), "Word size (empty = default)", NA, min = 4, step = 1)),
        shiny::column(6, shiny::numericInput(ns("perc_ident"), "Min % identity (server side)", 0, min = 0, max = 100))
      ),
      shiny::checkboxInput(ns("low_complexity"), "Filter low-complexity regions", TRUE),
      shiny::textInput(ns("extra"), "Additional Entrez query (optional)", "", width = "100%"),
      shiny::strong("Entrez query that will be sent:"),
      shiny::verbatimTextOutput(ns("entrez_preview"), placeholder = TRUE),
      shiny::radioButtons(ns("mode"), "Where to run BLAST",
        choices = c("In this app (results appear here)" = "app",
                    "On the NCBI website (opens your browser)" = "web"), selected = "app"),
      shiny::actionButton(ns("run"), "RUN BLAST", class = "btn-primary")
    ),
    shiny::column(
      7,
      shiny::h4("Status"),
      shiny::verbatimTextOutput(ns("status"), placeholder = TRUE),
      shiny::div(style = "margin-bottom:10px;",
        shiny::actionButton(ns("stop"), "STOP", icon = shiny::icon("stop"), class = "btn-danger"),
        " ",
        shiny::actionButton(ns("refresh"), "Refresh / reload", icon = shiny::icon("arrows-rotate"))),
      shiny::fluidRow(
        shiny::column(8, shiny::textInput(ns("rid_load"), NULL, width = "100%",
                                          placeholder = "Load results from an existing RID (valid ~24 h)")),
        shiny::column(4, shiny::actionButton(ns("load_rid"), "Load RID"))),
      shiny::uiOutput(ns("query_filter")),
      shiny::fluidRow(
        shiny::column(6, shiny::sliderInput(ns("min_ident"), "Show identity >=", 0, 100, 0, step = 0.5)),
        shiny::column(6, shiny::sliderInput(ns("min_qcov"), "Show query cover >=", 0, 100, 0, step = 1))
      ),
      shiny::tabsetPanel(id = ns("rtabs"),
        shiny::tabPanel("Descriptions", shiny::br(), DT::DTOutput(ns("table"))),
        shiny::tabPanel("Graphic Summary", shiny::br(),
          shiny::plotOutput(ns("gsummary"), height = "560px"),
          shiny::tags$small("Each bar is one alignment on the query; rows are hits in order of score. ",
                            "Colours follow the alignment score, as on the NCBI website.")),
        shiny::tabPanel("Alignments", shiny::br(),
          shiny::numericInput(ns("aln_n"), "Hits to show (used when no row is selected in Descriptions)",
                              10, min = 1, max = 100, step = 1, width = "320px"),
          shiny::uiOutput(ns("alignments"))),
        shiny::tabPanel("Taxonomy", shiny::br(),
          shiny::actionButton(ns("tax_fetch"), "Add lineage from NCBI Taxonomy (needs internet)"),
          shiny::tags$small(shiny::tags$br(),
            "Species are the first two words of each hit title (an estimate). Counts follow the filters above."),
          shiny::br(), shiny::br(), DT::DTOutput(ns("tax_table")))),
      shiny::br(),
      shiny::tags$div(
        shiny::strong("Download"), shiny::tags$small(" (selected rows; if none is selected, all rows shown): "),
        shiny::br(),
        shiny::downloadButton(ns("dl_fasta_full"), "FASTA (complete sequence)"),
        shiny::downloadButton(ns("dl_fasta_aln"), "FASTA (aligned sequences)"),
        shiny::downloadButton(ns("dl"), "Hit table (CSV)"),
        shiny::downloadButton(ns("dl_xml"), "XML")),
      shiny::div(style = "display:flex;gap:8px;align-items:center;flex-wrap:wrap;margin:10px 0",
        shiny::strong("Send to the MSA tab:"),
        shiny::span("the top"),
        shiny::numericInput(ns("n_msa"), NULL, 10, min = 1, max = 500, step = 1, width = "90px"),
        shiny::span("hits"),
        shiny::checkboxInput(ns("msa_query"), "with my query sequence(s)", TRUE),
        shiny::actionButton(ns("to_msa"), "Send to MSA", icon = shiny::icon("share-nodes"), class = "btn-primary"),
        shiny::tags$small("(if rows are selected in the table, those rows are sent; the part of each hit that aligns to the query is used)")),
      shiny::uiOutput(ns("ncbi_link")),
      shiny::h4("Preliminary interpretation"),
      shiny::uiOutput(ns("interp"))
    )
  )
}

blastServer <- function(id, consensus = shiny::reactive(NULL), email = "") {
  shiny::moduleServer(id, function(input, output, session) {

    em <- function() if (is.function(email)) email() else email

    rv <- shiny::reactiveValues(stage = "idle", rid = NULL, rtoe = NULL,
                                t_submit = NULL, next_poll = NULL, msg = "Ready.",
                                results = NULL, entrez = "", xml = NULL,
                                qnames = NULL, hitlist = 50, last_poll = NULL, tax = NULL, qrec = NULL)

    # --- sequence input ------------------------------------------------------
    shiny::observeEvent(input$use_consensus, {
      s <- consensus()
      if (is.null(s) || !nzchar(paste(s, collapse = ""))) {
        shiny::showNotification("There is no consensus from the assembly yet.", type = "warning")
      } else {
        # several contigs with names -> multi-FASTA, so that each one can be ticked separately
        txt <- if (length(s) > 1 && !is.null(names(s)) && !any(grepl("^>", s)))
          paste0(">", names(s), "\n", s, collapse = "\n") else paste(s, collapse = "\n")
        shiny::updateTextAreaInput(session, "query", value = txt)
      }
    })
    shiny::observeEvent(input$fasta_file, {
      shiny::updateTextAreaInput(session, "query",
        value = paste(readLines(input$fasta_file$datapath, warn = FALSE), collapse = "\n"))
    })

    # --- choose which sequences to BLAST (tick boxes) ----------------------------
    records <- shiny::debounce(shiny::reactive({
      if (!nzchar(trimws(input$query %||% ""))) return(NULL)
      tryCatch(blast_parse_records(input$query), error = function(e) NULL)
    }), 500)

    output$seq_select <- shiny::renderUI({
      rec <- records()
      shiny::req(rec, nrow(rec) > 0)
      keys <- rec$id
      prev <- shiny::isolate(input$sel)
      chosen_now <- if (length(prev) && any(prev %in% keys)) intersect(prev, keys) else keys
      labs <- sprintf("%s (%d nt)", rec$id, nchar(gsub("[^A-Za-z]", "", rec$seq)))
      shiny::tagList(
        shiny::checkboxGroupInput(session$ns("sel"),
          sprintf("Sequences to BLAST (%d found)", nrow(rec)),
          choices = stats::setNames(keys, labs), selected = chosen_now),
        shiny::actionLink(session$ns("sel_all"), "Select all"), " | ",
        shiny::actionLink(session$ns("sel_none"), "Clear"),
        shiny::br(), shiny::br())
    })
    shiny::observeEvent(input$sel_all, {
      rec <- records(); shiny::req(rec)
      shiny::updateCheckboxGroupInput(session, "sel", selected = rec$id)
    })
    shiny::observeEvent(input$sel_none, {
      shiny::updateCheckboxGroupInput(session, "sel", selected = character(0))
    })

    # --- Entrez query -------------------------------------------------------
    entrez <- shiny::reactive({
      blast_build_entrez(.split_terms(input$include), .split_terms(input$exclude),
                         isTRUE(input$drop_env), input$extra)
    })
    output$entrez_preview <- shiny::renderText({
      e <- entrez(); if (nzchar(e)) e else "(no organism restriction)"
    })

    shiny::observeEvent(input$check_taxa, {
      terms <- c(.split_terms(input$include), .split_terms(input$exclude))
      if (!length(terms)) return()
      res <- tryCatch(blast_resolve_taxa(terms, em()),
                      error = function(e) data.frame(Input = "Error", Scientific_name = conditionMessage(e)))
      output$taxa_table <- shiny::renderTable(res)
    })

    # --- submit job ----------------------------------------------------------
    shiny::observeEvent(input$run, {
      if (!identical(input$mode, "web") && rv$stage %in% c("submitting", "waiting")) {
        shiny::showNotification("A BLAST job is still running.", type = "warning"); return()
      }
      rec <- tryCatch(blast_parse_records(input$query), error = function(e) NULL)
      if (is.null(rec) || !nrow(rec)) {
        shiny::showNotification("The query sequence is empty.", type = "error"); return()
      }
      sel <- input$sel
      drawn <- !is.null(records()) && nrow(records()) > 0
      if (is.null(sel)) sel <- if (drawn) character(0) else rec$id   # nothing ticked vs. list not drawn yet
      if (!length(sel)) {
        shiny::showNotification("Tick at least one sequence to BLAST.", type = "warning"); return()
      }
      rec <- rec[rec$id %in% sel, , drop = FALSE]
      q <- tryCatch({
        rec$seq <- vapply(seq_len(nrow(rec)),
                          function(i) .blast_check_seq(rec$seq[i], rec$id[i]), character(1))
        rec
      }, error = function(e) {
        shiny::showNotification(conditionMessage(e), type = "error"); NULL })
      if (is.null(q)) return()
      fasta <- paste0(">", q$id, "\n", q$seq, collapse = "\n")
      hl <- input$hitlist
      if (is.null(hl) || is.na(hl) || hl < 1) hl <- 50

      if (identical(input$mode, "web")) {            # run on the NCBI website
        ws0 <- if (is.na(input$word_size)) NA else input$word_size
        params <- blast_web_params(fasta, input$db, input$program, entrez(), input$evalue, hl, ws0,
                                   isTRUE(input$low_complexity))
        ok <- tryCatch({ blast_open_site(params, session); TRUE }, error = function(e) {
          shiny::showNotification(paste("Could not open the browser:", conditionMessage(e)), type = "error"); FALSE })
        if (ok) rv$msg <- sprintf(paste0("Opened the NCBI BLAST website in your browser (%d sequence(s)). ",
                                         "Check the form there and press BLAST."), nrow(q))
        return()
      }

      rv$stage <- "submitting"; rv$msg <- "Sending to NCBI..."; rv$results <- NULL
      rv$entrez <- entrez(); rv$qnames <- q$id; rv$hitlist <- hl; rv$tax <- NULL
      rv$params <- list(db = input$db, program = input$program, evalue = input$evalue, hitlist = hl,
                        date = format(Sys.Date(), "%d %B %Y"), n_query = nrow(q))   # kept for the Methods paragraph
      rv$qrec <- data.frame(id = q$id, seq = q$seq, stringsAsFactors = FALSE)
      ws <- if (is.na(input$word_size)) NA else input$word_size
      sub <- tryCatch(
        blast_submit(fasta, input$db, input$program, rv$entrez, input$evalue,
                     hl, ws, isTRUE(input$low_complexity),
                     input$perc_ident %||% 0, em()),
        error = function(e) { rv$stage <- "error"; rv$msg <- paste("ERROR:", conditionMessage(e)); NULL })
      if (is.null(sub)) return()
      rv$rid <- sub$rid; rv$rtoe <- sub$rtoe
      rv$t_submit <- Sys.time()
      rv$next_poll <- Sys.time() + max(sub$rtoe, 10)
      rv$stage <- "waiting"
      rv$msg <- sprintf("Job submitted (%d sequence(s)). RID: %s", nrow(q), sub$rid)
    })

    # --- polling (at most once a minute after the initial wait) -------------------------
    shiny::observe({
      shiny::req(rv$stage == "waiting")
      shiny::invalidateLater(2000, session)
      shiny::isolate({
        elapsed <- as.numeric(difftime(Sys.time(), rv$t_submit, units = "secs"))
        if (elapsed > 20 * 60) {
          rv$stage <- "error"
          rv$msg <- sprintf("Exceeded 20 minutes. Check manually at NCBI with RID %s.", rv$rid); return()
        }
        if (Sys.time() < rv$next_poll) {
          rv$msg <- sprintf("Waiting for NCBI... %d s elapsed (RID %s)", round(elapsed), rv$rid)
          return()
        }
        rv$last_poll <- Sys.time()
        st <- tryCatch(blast_status(rv$rid), error = function(e) {
          list(status = "RETRY", hits = FALSE, err = conditionMessage(e)) })
        if (st$status == "READY") {
          if (!st$hits) {
            rv$results <- blast_parse_xml(paste0("<BlastOutput><BlastOutput_iterations><Iteration>",
              "<Iteration_query-len>1</Iteration_query-len><Iteration_hits/></Iteration>",
              "</BlastOutput_iterations></BlastOutput>"))
            rv$stage <- "done"; rv$msg <- "Finished: no hits."
          } else {
            xml_txt <- tryCatch(blast_fetch_xml(rv$rid),
                                error = function(e) { rv$msg <- paste("ERROR downloading:", conditionMessage(e)); NULL })
            res <- if (is.null(xml_txt)) NULL else
              tryCatch(blast_parse_xml(xml_txt, rv$qnames),
                       error = function(e) { rv$msg <- paste("ERROR reading results:", conditionMessage(e)); NULL })
            if (!is.null(res)) {
              rv$xml <- xml_txt
              rv$results <- res; rv$stage <- "done"
              rv$msg <- sprintf("Finished: %d hit(s) (RID %s).", nrow(res), rv$rid)
            } else rv$stage <- "error"
          }
        } else if (st$status %in% c("FAILED", "UNKNOWN")) {
          rv$stage <- "error"; rv$msg <- paste("BLAST failed / unknown RID:", st$status)
        } else {
          rv$next_poll <- Sys.time() + 60
          rv$msg <- sprintf("Status %s - checking again in 60 s (RID %s)", st$status, rv$rid)
        }
      })
    })

    output$status <- shiny::renderText(rv$msg)

    # --- STOP: stop waiting for the running job (the job itself may still finish at NCBI) ----
    shiny::observeEvent(input$stop, {
      if (identical(rv$stage, "waiting")) {
        rv$stage <- "stopped"
        rv$msg <- sprintf(paste0("Stopped. The job (RID %s) may still run at NCBI; use Refresh / reload ",
                                 "or Load RID to get its results later (valid ~24 h)."), rv$rid %||% "-")
      } else {
        shiny::showNotification("No BLAST job is waiting for results.", type = "message")
      }
    })

    # --- load / reload the results of an RID --------------------------------------------------
    reload_rid <- function(rid) {
      rid <- trimws(rid %||% "")
      if (!nzchar(rid)) { shiny::showNotification("Enter an RID first.", type = "warning"); return(invisible()) }
      if (identical(rv$stage, "waiting") && !identical(rid, rv$rid)) {
        shiny::showNotification("A job is still running. Press STOP first.", type = "warning"); return(invisible())
      }
      if (!identical(rid, rv$rid)) { rv$qnames <- NULL; rv$xml <- NULL; rv$qrec <- NULL }
      rv$last_poll <- Sys.time()
      st <- tryCatch(blast_status(rid), error = function(e)
        list(status = "ERROR", hits = FALSE, err = conditionMessage(e)))
      if (identical(st$status, "WAITING")) {
        rv$rid <- rid; rv$t_submit <- Sys.time(); rv$next_poll <- Sys.time() + 60; rv$stage <- "waiting"
        rv$msg <- sprintf("Job %s is still running at NCBI; waiting for it...", rid)
        return(invisible())
      }
      if (!identical(st$status, "READY")) {
        rv$msg <- sprintf("Cannot load RID %s: status %s%s", rid, st$status,
                          if (!is.null(st$err)) paste0(" (", st$err, ")") else "")
        return(invisible())
      }
      xml_txt <- if (st$hits) tryCatch(blast_fetch_xml(rid), error = function(e) {
        rv$msg <- paste("ERROR downloading:", conditionMessage(e)); NULL }) else .empty_blast_xml
      if (is.null(xml_txt)) return(invisible())
      res <- tryCatch(blast_parse_xml(xml_txt, rv$qnames), error = function(e) {
        rv$msg <- paste("ERROR reading results:", conditionMessage(e)); NULL })
      if (is.null(res)) return(invisible())
      rv$rid <- rid; rv$xml <- xml_txt; rv$results <- res; rv$stage <- "done"
      if (nrow(res)) rv$hitlist <- max(table(res$Query))
      rv$msg <- if (nrow(res)) sprintf("Loaded: %d hit(s) (RID %s).", nrow(res), rid) else
        sprintf("Loaded: no hits (RID %s).", rid)
      invisible()
    }

    shiny::observeEvent(input$refresh, {
      if (is.null(rv$rid)) {
        shiny::showNotification("There is no BLAST job yet. Run one, or use Load RID.", type = "message"); return()
      }
      if (identical(rv$stage, "waiting")) {
        since <- if (is.null(rv$last_poll)) Inf else as.numeric(difftime(Sys.time(), rv$last_poll, units = "secs"))
        if (since < 60) {
          shiny::showNotification(sprintf(
            "NCBI asks for at most one status check per minute. Try again in %d s.", ceiling(60 - since)),
            type = "message")
        } else rv$next_poll <- Sys.time()      # check now
      } else reload_rid(rv$rid)
    })
    shiny::observeEvent(input$load_rid, reload_rid(input$rid_load))

    # --- results --------------------------------------------------------------
    filtered <- shiny::reactive({
      r <- rv$results
      if (is.null(r)) return(NULL)
      r <- r[r$Identity >= input$min_ident & r$QueryCover >= input$min_qcov, , drop = FALSE]
      sq <- input$show_query
      if (!is.null(sq) && !identical(sq, "(all)") && sq %in% r$Query) r <- r[r$Query == sq, , drop = FALSE]
      r
    })

    output$query_filter <- shiny::renderUI({
      r <- rv$results; shiny::req(r)
      qs <- unique(r$Query)
      shiny::req(length(qs) > 1)
      shiny::selectInput(session$ns("show_query"), "Show results for", choices = c("(all)", qs),
                         selected = "(all)")
    })

    output$table <- DT::renderDT({
      r <- filtered(); shiny::req(r)
      r$Hseq <- NULL
      if (length(unique(r$Query)) < 2) r$Query <- NULL
      r$Accession <- sprintf('<a href="https://www.ncbi.nlm.nih.gov/nucleotide/%s" target="_blank">%s</a>',
                             r$Accession, r$Accession)
      # one page = Max target sequences (e.g. 50 -> 50 rows per page)
      page_len <- as.integer(rv$hitlist %||% 50)
      menu <- sort(unique(c(10, 25, 50, 100, page_len)))
      DT::datatable(r, escape = FALSE, rownames = FALSE, selection = "multiple",
                    options = list(pageLength = page_len,
                                   lengthMenu = list(c(menu, -1), c(menu, "All")),
                                   scrollX = TRUE))
    })

    # rows selected in the table; if none, all rows shown
    chosen <- shiny::reactive({
      r <- filtered(); shiny::req(r)
      s <- input$table_rows_selected
      if (length(s)) r[s, , drop = FALSE] else r
    })

    output$interp <- shiny::renderUI({
      r <- filtered(); shiny::req(r)
      qs <- unique(r$Query)
      if (length(qs) < 2) return(shiny::tags$ul(lapply(blast_interpret(r), shiny::tags$li)))
      shiny::tagList(lapply(qs, function(q) shiny::tagList(
        shiny::tags$strong(q),
        shiny::tags$ul(lapply(blast_interpret(r[r$Query == q, , drop = FALSE]), shiny::tags$li)))))
    })

    output$ncbi_link <- shiny::renderUI({
      shiny::req(rv$rid, rv$stage == "done")
      shiny::tags$p(shiny::tags$a("Open the results on the NCBI website (valid ~24 h)", target = "_blank",
        href = paste0(BLAST_URL, "?CMD=Get&RID=", rv$rid)))
    })

    output$dl <- shiny::downloadHandler(
      filename = function() paste0("blast_", rv$rid %||% "results", "_hits.csv"),
      content = function(file) {
        d <- chosen(); d$Hseq <- NULL
        utils::write.csv(d, file, row.names = FALSE)
      }
    )

    output$dl_fasta_aln <- shiny::downloadHandler(
      filename = function() paste0("blast_", rv$rid %||% "results", "_aligned.fasta"),
      content = function(file) writeLines(blast_fasta_aligned(chosen()), file)
    )

    output$dl_fasta_full <- shiny::downloadHandler(
      filename = function() paste0("blast_", rv$rid %||% "results", "_complete.fasta"),
      content = function(file) {
        d <- chosen()
        big <- !is.na(d$Length) & d$Length > BLAST_MAX_FULL
        if (any(big))
          shiny::showNotification(sprintf("%d hit(s) skipped because they are longer than %s bp; use FASTA (aligned sequences) for those.",
                                          sum(big), format(BLAST_MAX_FULL, big.mark = ",")), type = "warning", duration = 15)
        txt <- tryCatch(blast_fetch_fasta_complete(d$Accession[!big], em()),
          error = function(e) {
            shiny::showNotification(paste("Could not fetch sequences from NCBI:", conditionMessage(e)), type = "error", duration = 15)
            paste("# ERROR:", conditionMessage(e))
          })
        writeLines(txt, file)
      }
    )

    output$dl_xml <- shiny::downloadHandler(
      filename = function() paste0("blast_", rv$rid %||% "results", ".xml"),
      content = function(file) writeLines(rv$xml %||% "", file, useBytes = TRUE)
    )

    shiny::observeEvent(input$mode, {
      shiny::updateActionButton(session, "run",
        label = if (identical(input$mode, "web")) "OPEN NCBI BLAST" else "RUN BLAST")
    })

    # --- Graphic Summary / Alignments / Taxonomy (like the NCBI result tabs) -----------------------
    hsps <- shiny::reactive({
      shiny::req(rv$xml)
      tryCatch(blast_parse_hsps(rv$xml, rv$qnames), error = function(e) NULL)
    })

    output$gsummary <- shiny::renderPlot({
      r <- filtered(); shiny::req(r, nrow(r) > 0)
      hs <- hsps(); shiny::req(hs, nrow(hs) > 0)
      hs <- hs[paste(hs$Query, hs$Rank) %in% paste(r$Query, r$Rank), , drop = FALSE]
      shiny::req(nrow(hs) > 0)
      hs$Row <- stats::ave(hs$Rank, hs$Query, FUN = function(x) match(x, sort(unique(x))))
      hs$Band <- cut(hs$Bits, c(-Inf, 40, 50, 80, 200, Inf), right = FALSE,
                     labels = c("<40", "40-50", "50-80", "80-200", ">=200"))
      qbar <- unique(hs[, c("Query", "QLen")])
      p <- ggplot2::ggplot(hs) +
        ggplot2::geom_segment(ggplot2::aes(x = pmin(QFrom, QTo), xend = pmax(QFrom, QTo),
                                           y = Row, yend = Row, colour = Band), linewidth = 2.2) +
        ggplot2::geom_segment(data = qbar, ggplot2::aes(x = 1, xend = QLen, y = 0, yend = 0),
                              colour = "#2c7be5", linewidth = 4, inherit.aes = FALSE) +
        ggplot2::scale_colour_manual(name = "Alignment score (bits)", drop = FALSE,
          values = c("<40" = "black", "40-50" = "blue", "50-80" = "darkgreen",
                     "80-200" = "magenta", ">=200" = "red")) +
        ggplot2::scale_y_reverse() +
        ggplot2::labs(x = "Position on the query (nt)", y = "Hit rank") +
        ggplot2::theme_minimal(base_size = 13) +
        ggplot2::theme(legend.position = "top")
      if (length(unique(hs$Query)) > 1) p <- p + ggplot2::facet_wrap(ggplot2::vars(Query), ncol = 1, scales = "free_x")
      p
    })

    output$alignments <- shiny::renderUI({
      r <- chosen(); shiny::req(r, nrow(r) > 0)
      hs <- hsps(); shiny::req(hs, nrow(hs) > 0)
      if (!length(input$table_rows_selected)) {
        n <- input$aln_n
        if (is.null(n) || is.na(n) || n < 1) n <- 10
        r <- utils::head(r, n)
      }
      multi <- length(unique(rv$results$Query)) > 1
      blocks <- lapply(seq_len(nrow(r)), function(i) {
        h <- hs[hs$Query == r$Query[i] & hs$Rank == r$Rank[i], , drop = FALSE]
        h <- h[order(-h$Bits), , drop = FALSE]
        shiny::tagList(
          shiny::tags$h5(
            if (multi) shiny::tags$span(paste0("[", r$Query[i], "] ")),
            shiny::tags$a(r$Accession[i], target = "_blank",
                          href = paste0("https://www.ncbi.nlm.nih.gov/nucleotide/", r$Accession[i])),
            " ", r$Description[i]),
          lapply(seq_len(nrow(h)), function(j) {
            x <- h[j, ]
            gp <- if (is.na(x$Gaps)) 0 else x$Gaps
            strand <- paste0(if (!is.na(x$QFrame) && x$QFrame < 0) "Minus" else "Plus", "/",
                             if (!is.na(x$HFrame) && x$HFrame < 0) "Minus" else "Plus")
            stat <- sprintf("Score %.0f bits (%.0f), Expect %s, Identities %.0f/%.0f (%.0f%%), Gaps %.0f/%.0f (%.0f%%), Strand %s",
                            x$Bits, x$Score, format(x$Evalue, digits = 3, scientific = TRUE),
                            x$Identity, x$AlignLen, 100 * x$Identity / x$AlignLen,
                            gp, x$AlignLen, 100 * gp / x$AlignLen, strand)
            shiny::tagList(
              shiny::tags$div(style = "font-size:90%;", stat),
              shiny::tags$pre(style = "font-size:12px; overflow-x:auto;",
                blast_format_alignment(x$Qseq, x$Hseq, x$Midline, x$QFrom, x$HFrom,
                                       qminus = !is.na(x$QFrame) && x$QFrame < 0,
                                       hminus = !is.na(x$HFrame) && x$HFrame < 0)))
          }))
      })
      shiny::tagList(blocks)
    })

    tax_base <- shiny::reactive({
      r <- filtered(); shiny::req(r, nrow(r) > 0)
      d <- do.call(rbind, lapply(unique(r$Species), function(sp) {
        x <- r[r$Species == sp, , drop = FALSE]
        data.frame(Species = sp, Hits = nrow(x), BestIdentity = max(x$Identity, na.rm = TRUE),
                   BestAccession = x$Accession[which.max(x$Identity)], stringsAsFactors = FALSE)
      }))
      d[order(-d$Hits, -d$BestIdentity), , drop = FALSE]
    })

    shiny::observeEvent(input$tax_fetch, {
      b <- tax_base()
      sp <- utils::head(b$Species, 60)
      shiny::withProgress(message = "Fetching lineage from NCBI Taxonomy...", value = 0.3, {
        rv$tax <- tryCatch(blast_taxonomy_lineage(sp, em()), error = function(e) {
          shiny::showNotification(paste("Taxonomy lookup failed:", conditionMessage(e)), type = "error", duration = 15)
          NULL })
      })
    })

    output$tax_table <- DT::renderDT({
      b <- tax_base()
      if (!is.null(rv$tax)) {
        b <- merge(b, rv$tax, by = "Species", all.x = TRUE, sort = FALSE)
        b <- b[order(-b$Hits, -b$BestIdentity), , drop = FALSE]
        b$TaxID <- ifelse(is.na(b$TaxID), NA_character_, sprintf(
          '<a href="https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=%s" target="_blank">%s</a>',
          b$TaxID, b$TaxID))
      }
      b$BestAccession <- sprintf('<a href="https://www.ncbi.nlm.nih.gov/nucleotide/%s" target="_blank">%s</a>',
                                 b$BestAccession, b$BestAccession)
      DT::datatable(b, escape = FALSE, rownames = FALSE, options = list(pageLength = 25, scrollX = TRUE))
    })

    # --- send the top hits (and the query) to the MSA tab ---------------------------------------------
    outbox <- shiny::reactiveVal(NULL)
    shiny::observeEvent(input$to_msa, {
      r <- filtered()
      if (is.null(r) || !nrow(r)) { shiny::showNotification("There are no BLAST hits yet.", type = "warning"); return() }
      sel <- input$table_rows_selected
      d <- if (length(sel)) r[sel, , drop = FALSE] else {
        n <- input$n_msa
        if (is.null(n) || is.na(n) || n < 1) n <- 10
        do.call(rbind, lapply(split(r, factor(r$Query, levels = unique(r$Query))), utils::head, n))
      }
      d <- d[!is.na(d$Hseq) & nzchar(gsub("-", "", d$Hseq)), , drop = FALSE]
      d <- d[!duplicated(d$Accession), , drop = FALSE]
      if (!nrow(d)) { shiny::showNotification("The selected hits have no aligned sequence to send.", type = "warning"); return() }
      sp <- gsub("^_+|_+$", "", gsub("[^A-Za-z0-9]+", "_", d$Species))
      recs <- data.frame(id = ifelse(nzchar(sp), paste0(d$Accession, "_", sp), d$Accession),
                         seq = toupper(gsub("-", "", d$Hseq)), stringsAsFactors = FALSE)
      nq <- 0L
      if (isTRUE(input$msa_query) && !is.null(rv$qrec) && nrow(rv$qrec)) {
        recs <- rbind(data.frame(id = rv$qrec$id, seq = rv$qrec$seq, stringsAsFactors = FALSE), recs)
        nq <- nrow(rv$qrec)
      }
      outbox(list(recs = recs, nq = nq, nh = nrow(d), stamp = Sys.time()))
    })

    # ---- project file: what is saved / restored (used by the Project tab) -----------------------------
    get_state <- function() shiny::isolate(list(results = rv$results, qrec = rv$qrec, rid = rv$rid, xml = rv$xml,
                                                entrez = rv$entrez, params = rv$params, query = input$query))
    set_state <- function(s) {
      rv$results <- s$results; rv$qrec <- s$qrec; rv$rid <- s$rid; rv$xml <- s$xml
      rv$entrez <- if (is.null(s$entrez)) "" else s$entrez
      rv$params <- s$params
      rv$stage <- if (!is.null(s$results)) "done" else "idle"
      rv$msg <- if (!is.null(s$results)) "Loaded from a project file." else "Ready."
      if (!is.null(s$query)) shiny::updateTextAreaInput(session, "query", value = s$query)
      invisible()
    }
    params <- function() {
      p <- shiny::isolate(rv$params)
      if (is.null(p)) p <- shiny::isolate(list(db = input$db, program = input$program, evalue = input$evalue,
                                               hitlist = rv$hitlist, date = NA_character_, n_query = NA_integer_))
      p$rid <- shiny::isolate(rv$rid)
      p
    }

    list(get_state = get_state, set_state = set_state, params = params,
         results = filtered, rid = shiny::reactive(rv$rid),
         entrez = shiny::reactive(rv$entrez), to_msa = shiny::reactive(outbox()))
  })
}
