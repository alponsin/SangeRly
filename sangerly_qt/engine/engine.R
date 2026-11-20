# SangeRly engine: a long-running R worker that the Qt front-end talks to.
# Protocol: one JSON object per line on stdin  {"id":1,"action":"assemble","args":{...}}
#           one JSON object per line on stdout {"id":1,"ok":true,"result":{...}}  or  {"id":1,"ok":false,"error":"..."}
#           progress lines                     {"event":"progress","id":1,"i":2,"n":5,"name":"S2"}
# Everything R prints with message()/warning() goes to stderr (shown in the Log panel).
# Start:  Rscript engine.R <project folder that contains R/>
args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) args[1] else Sys.getenv("SANGERLY_ROOT", getwd())
root <- normalizePath(root, winslash = "/", mustWork = TRUE)
suppressPackageStartupMessages({ library(jsonlite) })
old_wd <- setwd(root)
for (f in c("assembly_core.R", "seqtools_core.R", "stats_core.R", "pileup_core.R", "batch_core.R"))
  if (file.exists(file.path("R", f))) source(file.path("R", f), local = FALSE)
setwd(old_wd)

out <- function(x) { writeLines(enc2utf8(as.character(toJSON(x, auto_unbox = TRUE, null = "null", na = "null", digits = NA))), stdout(), useBytes = TRUE); flush(stdout()) }   # raw UTF-8 bytes (Windows would re-encode with cat)
.num <- function(x, d) if (is.null(x) || length(x) == 0 || is.na(suppressWarnings(as.numeric(x)))) d else as.numeric(x)

make_params <- function(p) {
  P <- default_params()
  if (is.null(p)) return(P)
  for (k in names(p)) if (k %in% names(P)) P[[k]] <- if (is.character(P[[k]])) as.character(p[[k]]) else as.numeric(p[[k]])
  P
}

df_rows <- function(d) if (is.null(d) || !nrow(d)) list() else lapply(seq_len(nrow(d)), function(i) as.list(d[i, , drop = FALSE]))

act_ping <- function(a) list(r = R.version.string, root = root,
  has_biostrings = requireNamespace("Biostrings", quietly = TRUE), pid = Sys.getpid())

# file list -> sample / direction using the same rules as the Shiny app
act_scan <- function(a) {
  paths <- unlist(a$paths); mode <- if (is.null(a$naming)) "name" else a$naming
  names <- basename(paths)
  sm <- switch(mode,
    one_per_file = tools::file_path_sans_ext(names),
    custom = rep(if (is.null(a$custom) || !nzchar(a$custom)) "sample" else a$custom, length(names)),
    group_label(names, mode))
  list(files = lapply(seq_along(paths), function(i)
    list(name = names[i], path = paths[i], sample = sm[i], direction = ab1_direction(names[i]))))
}

act_assemble <- function(a) {
  fl <- a$files
  files <- data.frame(name = vapply(fl, function(x) x$name, ""), path = vapply(fl, function(x) x$path, ""),
                      sample = vapply(fl, function(x) x$sample, ""), stringsAsFactors = FALSE)
  P <- make_params(a$params); mode <- if (identical(a$mode, "reference")) "ref" else "denovo"
  refs <- NULL
  if (mode == "ref") {
    if (is.null(a$reference) || !file.exists(a$reference)) stop("Reference FASTA not found.")
    refs <- read_fasta_all(a$reference)
  }
  id <- a$.id
  prog <- function(i, n, name) out(list(event = "progress", id = id, i = i, n = n, name = name))
  res <- build_batch(files, mode, P, refs, progress = prog)
  outdir <- a$out_dir; what <- unlist(a$save)
  if (!is.null(outdir) && nzchar(outdir)) dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  saved <- character(); samples <- list(); fasta_all <- character()
  for (nm in names(res)) {
    aln <- res[[nm]]
    if (!isTRUE(aln$ok)) { samples[[length(samples) + 1]] <- list(sample = nm, ok = FALSE, msg = aln$msg); next }
    cf <- aln_final(aln, aln$cols)
    cons <- paste(cf$Final[cf$Final != "-"], collapse = "")
    vr <- aln_variants(aln, cf)
    nrev <- sum(cf$Review & !cf$Manual)
    fasta_all <- c(fasta_all, paste0(">", nm), cons)
    if (!is.null(outdir) && nzchar(outdir)) {
      sfx <- function(s) file.path(outdir, paste0(gsub("[^A-Za-z0-9._-]+", "_", nm), s))
      if ("consensus" %in% what) { writeLines(c(paste0(">", nm), cons), f <- sfx("_consensus.fasta")); saved <- c(saved, f) }
      if ("alignment" %in% what) { writeLines(aln_alignment_fasta(aln, cf), f <- sfx("_alignment.fasta")); saved <- c(saved, f) }
      if ("decisions" %in% what) { utils::write.csv(cf[, intersect(c("Col", "RefPos", "RefBase", "Cov", "Type", "Auto", "Final", "Supported", "Recommendation", "Reason", "Review"), names(cf))],
                                                    f <- sfx("_decisions.csv"), row.names = FALSE); saved <- c(saved, f) }
      if ("variants" %in% what) { utils::write.csv(vr, f <- sfx("_variants.csv"), row.names = FALSE); saved <- c(saved, f) }
    }
    samples[[length(samples) + 1]] <- list(sample = nm, ok = TRUE, mode = aln$mode, notes = unlist(aln$notes),
      length = nchar(cons), columns = nrow(cf), reads = unlist(aln$labels), review = nrev,
      consensus = cons, variants = df_rows(vr),
      ambiguous = sum(!strsplit(cons, "")[[1]] %in% c("A", "C", "G", "T")))
  }
  if (!is.null(outdir) && nzchar(outdir) && "consensus_all" %in% what && length(fasta_all)) {
    writeLines(fasta_all, f <- file.path(outdir, "all_consensus.fasta")); saved <- c(saved, f)
  }
  list(samples = samples, saved = as.list(saved))
}

handlers <- list(ping = act_ping, scan = act_scan, assemble = act_assemble)
con <- file("stdin", "r")
out(list(event = "ready", pid = Sys.getpid(), r = R.version.string))
repeat {
  ln <- readLines(con, n = 1, warn = FALSE, encoding = "UTF-8")
  if (!length(ln)) break
  if (!nzchar(trimws(ln))) next
  req <- tryCatch(fromJSON(ln, simplifyVector = FALSE), error = function(e) NULL)
  if (is.null(req)) { out(list(ok = FALSE, error = "Bad JSON request")); next }
  if (identical(req$action, "quit")) break
  h <- handlers[[req$action]]
  if (is.null(h)) { out(list(id = req$id, ok = FALSE, error = paste("Unknown action:", req$action))); next }
  a <- req$args; if (is.null(a)) a <- list(); a$.id <- req$id
  r <- tryCatch(list(ok = TRUE, result = h(a)), error = function(e) list(ok = FALSE, error = conditionMessage(e)))
  r$id <- req$id; out(r)
}
