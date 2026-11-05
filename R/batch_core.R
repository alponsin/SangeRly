# ==============================================================================
# Sangerly - batch helpers (no Shiny): many .ab1 files -> samples -> one alignment per sample
# Needs (sourced first): assembly_core.R, pileup_core.R
# ==============================================================================

# Direction token in a file name: _F / -R / .fwd / 27F / 1492R / forward / reverse ...
.dir_rx <- "(^|[^A-Za-z0-9])[0-9]*(forward|reverse|fwd|rev|F|R)(?![A-Za-z])"

ab1_direction <- function(name) {
  s <- tools::file_path_sans_ext(basename(name))
  m <- gregexpr(.dir_rx, s, perl = TRUE, ignore.case = TRUE)[[1]]
  if (m[1] < 0) return(NA_character_)
  tok <- tolower(regmatches(s, list(m))[[1]])
  tok <- tok[length(tok)]
  tok <- gsub("[^a-z]", "", tok)
  if (tok %in% c("f", "fwd", "forward")) "F" else "R"
}

# Sample name = file name without the direction token (and, optionally, without long run numbers,
# e.g. the 1st BASE / Macrogen run IDs that differ between F and R of the same sample)
ab1_sample <- function(name, ignore_run = FALSE) {
  s <- tools::file_path_sans_ext(basename(name))
  m <- gregexpr(.dir_rx, s, perl = TRUE, ignore.case = TRUE)[[1]]
  if (m[1] >= 0) {
    k <- length(m)
    a <- m[k]; l <- attr(m, "match.length")[k]
    tok <- substr(s, a, a + l - 1L)
    rest <- substr(s, a + l, nchar(s))
    # primer code right after the tag (1st BASE: ..._F-1119 / ..._R-2092) differs between F and R: drop it too
    rest <- sub("^[-_.][0-9]{3,5}(?![0-9A-Za-z])", "", rest, perl = TRUE)
    s <- paste0(substr(s, 1, a - 1L), rest)
  }
  if (ignore_run) s <- gsub("[0-9]{6,}", "", s)
  s <- gsub("[_ .-]{2,}", "_", s)
  s <- gsub("^[_ .-]+|[_ .-]+$", "", s)
  if (!nzchar(s)) s <- tools::file_path_sans_ext(basename(name))
  s
}

# files: data.frame(name, path, ...). mode: "name" | "name_norun" | "all"
group_label <- function(names, mode = "name") {
  if (identical(mode, "all")) return(rep("sample", length(names)))
  vapply(names, ab1_sample, "", ignore_run = identical(mode, "name_norun"), USE.NAMES = FALSE)
}

# all records of a FASTA file as list(list(name, seq), ...)
read_fasta_all <- function(path) {
  x <- Biostrings::readDNAStringSet(path)
  if (!length(x)) stop("The FASTA file is empty.")
  lapply(seq_along(x), function(i) list(name = names(x)[i], seq = toupper(as.character(x[[i]]))))
}

# Best reference record for a sample (first usable read, both strands, local alignment score)
pick_reference <- function(raws, refs, P) {
  if (length(refs) == 1) return(1L)
  mat <- .pw("nucleotideSubstitutionMatrix")(match = 2, mismatch = -3, baseOnly = FALSE)
  rd <- NULL
  for (r in raws) { rd <- tryCatch(prep_read(r, P), error = function(e) NULL); if (!is.null(rd)) break }
  if (is.null(rd)) return(1L)
  q <- rd$seq; if (nchar(q) > 600) q <- substr(q, 1, 600)
  qr <- revcomp_chr(q)
  sc <- vapply(refs, function(r) {
    s <- substr(r$seq, 1, 30000)
    max(as.numeric(.pw("score")(.map_run(q, s, mat))), as.numeric(.pw("score")(.map_run(qr, s, mat))))
  }, 1)
  which.max(sc)
}

# One sample -> alignment object (or list(ok = FALSE, msg)).
# files: rows (name, path) of this sample. mode "denovo" needs exactly one F and one R.
build_sample <- function(files, mode, P, refs = NULL, reader = read_ab1) {
  load_notes <- character(); raws <- list(); labs <- character()
  for (i in seq_len(nrow(files))) {
    r <- tryCatch(reader(files$path[i], files$name[i]),
                  error = function(e) { load_notes <<- c(load_notes, conditionMessage(e)); NULL })
    if (!is.null(r)) { raws[[length(raws) + 1]] <- r; labs <- c(labs, tools::file_path_sans_ext(files$name[i])) }
  }
  if (!length(raws)) return(list(ok = FALSE, msg = paste("No AB1 file could be read.", paste(load_notes, collapse = " "))))
  if (identical(mode, "denovo")) {
    if (length(raws) != 2)
      return(list(ok = FALSE, msg = sprintf(paste0("De novo mode needs exactly 2 reads (Forward + Reverse) per sample; this sample has %d. ",
        "Fix the Sample column, or use reference-guided mode."), length(raws))))
    d <- vapply(files$name[seq_len(nrow(files))][seq_along(raws)], ab1_direction, "")
    if (!anyNA(d) && d[1] == "R" && d[2] == "F") { raws <- raws[2:1]; labs <- labs[2:1] }
    aln <- build_aln_denovo(raws[[1]], raws[[2]], P, labels = c("Forward", "Reverse"))
    if (isTRUE(aln$ok)) {
      aln$labels <- make.unique(labs)
      if (!anyNA(d) && d[1] == d[2]) aln$notes <- c(sprintf("Both files look like %s reads by name; the first was used as Forward.", d[1]), aln$notes)
    }
  } else {
    if (is.null(refs) || !length(refs)) return(list(ok = FALSE, msg = "No reference loaded."))
    k <- pick_reference(raws, refs, P)
    aln <- build_aln_ref(raws, refs[[k]], P, labs)
    if (isTRUE(aln$ok) && length(refs) > 1)
      aln$notes <- c(sprintf("Reference chosen automatically: %s.", substr(refs[[k]]$name, 1, 60)), aln$notes)
  }
  if (isTRUE(aln$ok)) aln$notes <- c(load_notes, aln$notes)
  aln
}

# Several samples; one failing sample does not stop the others. progress(i, n, name) is optional.
build_batch <- function(files, mode, P, refs = NULL, progress = NULL, reader = read_ab1) {
  sm <- unique(files$sample)
  out <- list()
  for (i in seq_along(sm)) {
    if (is.function(progress)) progress(i, length(sm), sm[i])
    f <- files[files$sample == sm[i], , drop = FALSE]
    out[[sm[i]]] <- tryCatch(build_sample(f, mode, P, refs, reader),
                             error = function(e) list(ok = FALSE, msg = conditionMessage(e)))
  }
  out
}
