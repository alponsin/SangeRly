# ==============================================================================
# SangeRly - Sanger assembly core (no Shiny)
# ------------------------------------------------------------------------------
# Flow: read AB1 -> quality trimming -> orient R (automatic) -> overlap alignment
# of F vs R. Per-column decisions live in pileup_core.R (rule-based, transparent:
# Phred and secondary-peak ratio).
#
# Packages: sangerseqR (Bioconductor), Biostrings, pwalign (Bioconductor >= 3.19;
# older Biostrings is used if pwalign is missing), BiocGenerics.
#   if (!requireNamespace("BiocManager")) install.packages("BiocManager")
#   BiocManager::install(c("sangerseqR", "Biostrings", "pwalign", "BiocGenerics"))
# ==============================================================================

.nz <- function(a, b) if (is.null(a) || length(a) == 0) b else a

COMP_CHAR <- c(A = "T", C = "G", G = "C", T = "A", N = "N", R = "Y", Y = "R",
               S = "S", W = "W", K = "M", M = "K", B = "V", V = "B", D = "H", H = "D")

revcomp_chr <- function(s) {
  ch <- strsplit(s, "")[[1]]
  out <- unname(COMP_CHAR[ch])
  out[is.na(out)] <- "N"
  paste(rev(out), collapse = "")
}

# True if a package is installed, WITHOUT loading it (requireNamespace() loads the whole package, which is
# slow for DECIPHER / msa / ape and was done at every app start just to fill a dropdown).
.has_pkg <- function(p) nzchar(system.file(package = p))

# Alignment functions come from pwalign (new) or Biostrings (old)
.pw <- function(fn) {
  pkg <- if (requireNamespace("pwalign", quietly = TRUE)) "pwalign" else "Biostrings"
  getExportedValue(pkg, fn)
}

default_params <- function() {
  list(trim_method = "M2",   # "M2" = sliding window (Trimmomatic-like), "M1" = modified Mott
       q_trim = 20,          # M2: minimum mean Phred in the window
       win = 10,             # M2: window size
       mott_cutoff = 1e-4,   # M1: error cutoff (sangeranalyseR default)
       q_good = 30,          # Phred >= this = good
       q_poor = 20,          # Phred <  this = poor
       ratio = 0.33,         # secondary/primary peak ratio >= this = genuine double peak
       dq_pref = 15,         # minimum Phred-support difference to favour one base
       min_overlap = 30,     # minimum F/R overlap (bp)
       min_ident = 0.85,     # minimum identity in the overlap
       type4 = "N")          # type 4 conflict: "N", "iupac", or "F" (follow best-supported base)
}

# ------------------------------------------------------------------------------
# Reading AB1
# ------------------------------------------------------------------------------
# Minimal ABIF reader in base R (does not depend on sangerseqR's handling of data types).
# Returns a list keyed "TAG.number" with numeric vectors (types 1, 3, 4, 5, 7, 8, 10, 11, 12, 13) or a
# character string (types 2, 18, 19, 20, 21).
.abif_read <- function(path) {
  con <- file(path, "rb"); on.exit(close(con))
  raw <- readBin(con, "raw", n = file.size(path))
  if (length(raw) < 40 || rawToChar(raw[1:4]) != "ABIF") stop("not an ABIF file")
  i32 <- function(o) readBin(raw[o + 0:3], "integer", n = 1, size = 4, endian = "big")
  u16 <- function(o) readBin(raw[o + 0:1], "integer", n = 1, size = 2, signed = FALSE, endian = "big")
  nent <- i32(19); doff <- i32(27)
  out <- list()
  for (e in seq_len(nent)) {
    o <- doff + (e - 1L) * 28L + 1L
    tag <- rawToChar(raw[o:(o + 3L)]); num <- i32(o + 4L)
    typ <- u16(o + 8L); cnt <- i32(o + 12L); size <- i32(o + 16L)
    start <- if (size <= 4) o + 20L else i32(o + 20L) + 1L
    if (size <= 0 || cnt <= 0) next
    bytes <- raw[start:(start + size - 1L)]
    val <- switch(as.character(typ),
      "1"  = as.integer(bytes),
      "2"  = as.integer(bytes),                                    # bytes kept as integers; text made below
      "3"  = as.integer(bytes),
      "4"  = readBin(bytes, "integer", n = cnt, size = 2, signed = TRUE,  endian = "big"),
      "5"  = readBin(bytes, "integer", n = cnt, size = 4, endian = "big"),
      "18" = , "19" = , "10" = , "11" = , "12" = , "13" = , "7" = , "8" = NULL,
      NULL)
    if (!is.null(val)) out[[paste0(tag, ".", num)]] <- list(type = typ, val = val)
  }
  out
}

read_ab1 <- function(path, name = basename(path)) {
  d <- tryCatch(.abif_read(path), error = function(e) stop(name, ": ", conditionMessage(e)))
  get <- function(...) { for (k in c(...)) if (!is.null(d[[k]])) return(d[[k]]$val); NULL }
  seq  <- toupper(rawToChar(as.raw(get("PBAS.2", "PBAS.1"))))
  qual <- as.integer(get("PCON.2", "PCON.1"))
  ploc <- as.integer(get("PLOC.2", "PLOC.1"))
  fwo  <- toupper(rawToChar(as.raw(.nz(get("FWO_.1"), utf8ToInt("GATC")))))
  if (nchar(fwo) != 4) fwo <- "GATC"
  ch <- lapply(9:12, function(k) as.numeric(.nz(get(paste0("DATA.", k)), numeric(0))))
  if (any(vapply(ch, length, 1L) == 0)) stop(name, ": trace data (DATA.9-12) not found.")
  traces <- do.call(cbind, ch)
  colnames(traces) <- strsplit(fwo, "")[[1]]
  n <- nchar(seq)
  if (n == 0 || length(qual) != n || length(ploc) != n)
    stop(name, sprintf(": base calls / quality / peak positions are inconsistent (non-standard AB1 file?) [bases=%d, quality=%d, peaks=%d]",
                       n, length(qual), length(ploc)))
  list(name = name, seq = seq, qual = qual, ploc = ploc, traces = traces)
}

# ------------------------------------------------------------------------------
# Trimming
# ------------------------------------------------------------------------------
# M2: first to last window whose mean Phred >= q
trim_range <- function(qual, q = 20, win = 10) {
  n <- length(qual)
  if (n < win) return(c(NA_integer_, NA_integer_))
  cs <- c(0, cumsum(qual))
  m <- (cs[(win + 1):(n + 1)] - cs[1:(n - win + 1)]) / win
  ok <- which(m >= q)
  if (!length(ok)) return(c(NA_integer_, NA_integer_))
  c(min(ok), max(ok) + win - 1L)
}

# M1: modified Mott (as in Phred/Phrap and Biopython).
# Base score = cutoff - 10^(-Q/10); keep the segment with the maximum summed score.
trim_mott <- function(qual, cutoff = 1e-4) {
  sc <- cutoff - 10^(-qual / 10)
  cur <- 0; start <- 1L; best <- 0; bs <- NA_integer_; be <- NA_integer_
  for (i in seq_along(sc)) {
    cur <- cur + sc[i]
    if (cur <= 0) { cur <- 0; start <- i + 1L }
    else if (cur > best) { best <- cur; bs <- start; be <- i }
  }
  c(bs, be)
}

# P: parameter list (default_params()). Raw data (original orientation) is kept too.
prep_read <- function(raw, P = default_params()) {
  method <- .nz(P$trim_method, "M2")
  rng <- if (identical(method, "M1")) trim_mott(raw$qual, .nz(P$mott_cutoff, 1e-4))
         else trim_range(raw$qual, P$q_trim, P$win)
  if (anyNA(rng)) stop(raw$name, ": no region passes trimming (",
                       if (identical(method, "M1")) "Mott" else paste0("mean Q", P$q_trim),
                       "). Relax the trimming parameters or resequence.")
  idx <- rng[1]:rng[2]
  list(name = raw$name,
       seq = paste(strsplit(raw$seq, "")[[1]][idx], collapse = ""),
       qual = raw$qual[idx], ploc = raw$ploc[idx], idx = idx,
       traces = raw$traces, rc = FALSE,
       raw = list(seq = strsplit(raw$seq, "")[[1]], qual = raw$qual, ploc = raw$ploc),
       qc = list(qual_full = raw$qual, from = rng[1], to = rng[2], n_full = length(raw$qual)))
}

# Reverse-complement a read (sequence, quality, peak positions, trace channel names)
rc_read <- function(rd) {
  rd$seq <- revcomp_chr(rd$seq)
  rd$qual <- rev(rd$qual); rd$ploc <- rev(rd$ploc); rd$idx <- rev(rd$idx)
  colnames(rd$traces) <- unname(COMP_CHAR[colnames(rd$traces)])
  rd$rc <- !rd$rc
  rd
}

# ------------------------------------------------------------------------------
# Overlap alignment F vs R
# ------------------------------------------------------------------------------
align_overlap <- function(s1, s2) {
  mat <- .pw("nucleotideSubstitutionMatrix")(match = 2, mismatch = -3, baseOnly = FALSE)
  aln <- .pw("pairwiseAlignment")(pattern = s1, subject = s2, type = "overlap",
                                  substitutionMatrix = mat, gapOpening = 8, gapExtension = 4)
  pw <- .pw("pattern")(aln); sj <- .pw("subject")(aln)
  a1 <- as.character(.pw("alignedPattern")(aln))
  a2 <- as.character(.pw("alignedSubject")(aln))
  p1 <- c(BiocGenerics::start(pw), BiocGenerics::end(pw))
  p2 <- c(BiocGenerics::start(sj), BiocGenerics::end(sj))
  c1 <- strsplit(a1, "")[[1]]; c2 <- strsplit(a2, "")[[1]]
  both <- c1 != "-" & c2 != "-"
  list(score = as.numeric(.pw("score")(aln)), a1 = a1, a2 = a2, p1 = p1, p2 = p2,
       overlap_len = length(c1),
       identity = if (any(both)) mean(c1[both] == c2[both]) else 0)
}

# Build the two full alignment rows (including the non-overlapping ends)
build_rows <- function(s1, s2, al) {
  lead1 <- if (al$p1[1] > 1) substr(s1, 1, al$p1[1] - 1) else ""
  lead2 <- if (al$p2[1] > 1) substr(s2, 1, al$p2[1] - 1) else ""
  tail1 <- substr(s1, al$p1[2] + 1, nchar(s1))
  tail2 <- substr(s2, al$p2[2] + 1, nchar(s2))
  g <- function(n) strrep("-", n)
  row1 <- paste0(lead1, g(nchar(lead2)), al$a1, tail1, g(nchar(tail2)))
  row2 <- paste0(g(nchar(lead1)), lead2, al$a2, g(nchar(tail1)), tail2)
  ov_start <- nchar(lead1) + nchar(lead2) + 1
  list(row1 = row1, row2 = row2, ov = c(ov_start, ov_start + nchar(al$a1) - 1))
}

# IUPAC code for a set of bases
.iupac <- function(x) {
  s <- sort(unique(x[x %in% c("A", "C", "G", "T")]))
  if (length(s) == 1) return(s)
  tbl <- c(AG = "R", CT = "Y", CG = "S", AT = "W", GT = "K", AC = "M",
           ACG = "V", ACT = "H", AGT = "D", CGT = "B")
  out <- unname(tbl[paste(s, collapse = "")])
  if (length(out) == 0 || is.na(out)) "N" else out
}

# Orient R, align the pair, return the two aligned rows plus per-column read indices.
pair_alignment <- function(f_raw, r_raw, P = default_params()) {
  f  <- prep_read(f_raw, P)
  r0 <- prep_read(r_raw, P)
  r1 <- rc_read(r0)

  a_rc <- align_overlap(f$seq, r1$seq)
  a_as <- align_overlap(f$seq, r0$seq)
  use_rc <- a_rc$score >= a_as$score
  r  <- if (use_rc) r1 else r0
  al <- if (use_rc) a_rc else a_as

  out <- list(ok = FALSE, f = f, r = r, r_was_rc = use_rc, P = P,
              overlap_len = al$overlap_len, overlap_identity = al$identity)
  if (al$overlap_len < P$min_overlap || al$identity < P$min_ident) {
    out$msg <- sprintf(paste0("F and R do not overlap convincingly (overlap %d bp, identity %.1f%%; ",
                              "required: >= %d bp and >= %.0f%%). Check that the file pair is correct, ",
                              "or relax the trimming parameters."),
                       al$overlap_len, 100 * al$identity, P$min_overlap, 100 * P$min_ident)
    return(out)
  }

  rows <- build_rows(f$seq, r$seq, al)
  c1 <- strsplit(rows$row1, "")[[1]]; c2 <- strsplit(rows$row2, "")[[1]]
  out$ok <- TRUE
  out$cols <- data.frame(Fwd = c1, Rev = c2,
                         Fwd_idx = ifelse(c1 == "-", NA_integer_, cumsum(c1 != "-")),
                         Rev_idx = ifelse(c2 == "-", NA_integer_, cumsum(c2 != "-")),
                         stringsAsFactors = FALSE)
  out
}

# Apply manual decisions (ov: named character vector; names = column numbers)
final_cols <- function(cols, ov = character()) {
  cols$Final <- cols$Auto
  cols$Manual <- FALSE
  if (length(ov)) {
    k <- as.integer(names(ov))
    cols$Final[k] <- unname(ov)
    cols$Manual[k] <- TRUE
  }
  cols
}

# ------------------------------------------------------------------------------
# Reference sequences
# ------------------------------------------------------------------------------
read_fasta_first <- function(path) {
  x <- Biostrings::readDNAStringSet(path)
  if (!length(x)) stop("The FASTA file is empty.")
  list(name = names(x)[1], seq = toupper(as.character(x[[1]])))
}

fetch_ncbi_fasta <- function(accession, email = "") {
  r <- httr::GET("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
                 query = list(db = "nuccore", id = accession, rettype = "fasta", retmode = "text",
                              tool = "Sangerly", email = email), httr::timeout(60))
  httr::stop_for_status(r)
  txt <- httr::content(r, as = "text", encoding = "UTF-8")
  lines <- strsplit(txt, "\r?\n")[[1]]
  if (!length(lines) || !grepl("^>", lines[1])) stop("Accession not found: ", accession)
  list(name = sub("^>", "", lines[1]), seq = toupper(paste(lines[-1], collapse = "")))
}
