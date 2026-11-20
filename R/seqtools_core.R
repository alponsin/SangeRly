# SangeRly - small sequence tools used by the "Sequence analysis" tab (no Shiny):
# nucleotide composition and restriction-site search.

# ---- common restriction enzymes: name, recognition site (IUPAC), cut after the n-th base of the site
.sa_enzymes <- data.frame(
  Enzyme = c("AflII", "AgeI", "AluI", "ApaI", "AvaI", "BamHI", "BclI", "BglII", "BstEII", "ClaI",
             "DdeI", "DpnII", "DraI", "EagI", "EcoRI", "EcoRV", "FspI", "HaeIII", "HhaI", "HincII",
             "HindIII", "HinfI", "HpaI", "HpaII", "KpnI", "MboI", "MfeI", "MluI", "MscI", "MspI",
             "NaeI", "NcoI", "NdeI", "NheI", "NotI", "NruI", "NsiI", "PstI", "PvuI", "PvuII",
             "RsaI", "SacI", "SacII", "Sau3AI", "SalI", "ScaI", "SmaI", "SnaBI", "SphI", "SpeI",
             "SspI", "StuI", "TaqI", "XbaI", "XhoI", "XmaI"),
  Site   = c("CTTAAG", "ACCGGT", "AGCT", "GGGCCC", "CYCGRG", "GGATCC", "TGATCA", "AGATCT", "GGTNACC", "ATCGAT",
             "CTNAG", "GATC", "TTTAAA", "CGGCCG", "GAATTC", "GATATC", "TGCGCA", "GGCC", "GCGC", "GTYRAC",
             "AAGCTT", "GANTC", "GTTAAC", "CCGG", "GGTACC", "GATC", "CAATTG", "ACGCGT", "TGGCCA", "CCGG",
             "GCCGGC", "CCATGG", "CATATG", "GCTAGC", "GCGGCCGC", "TCGCGA", "ATGCAT", "CTGCAG", "CGATCG", "CAGCTG",
             "GTAC", "GAGCTC", "CCGCGG", "GATC", "GTCGAC", "AGTACT", "CCCGGG", "TACGTA", "GCATGC", "ACTAGT",
             "AATATT", "AGGCCT", "TCGA", "TCTAGA", "CTCGAG", "CCCGGG"),
  Cut    = c(1, 1, 2, 5, 1, 1, 1, 1, 1, 2,
             1, 0, 3, 1, 1, 3, 3, 2, 3, 3,
             1, 1, 3, 1, 5, 0, 1, 1, 3, 1,
             3, 1, 2, 1, 2, 3, 5, 5, 4, 3,
             2, 5, 4, 0, 1, 3, 3, 3, 5, 1,
             3, 3, 1, 1, 1, 1),
  stringsAsFactors = FALSE)

.sa_composition <- function(id, seq) {
  b <- strsplit(seq, "")[[1]]
  n <- length(b)
  a <- sum(b == "A"); cc <- sum(b == "C"); g <- sum(b == "G"); t <- sum(b == "T")
  acgt <- a + cc + g + t
  data.frame(Sequence = id, Length = n, A = a, C = cc, G = g, T = t, Other = n - acgt,
             GC_pct = round(100 * (g + cc) / acgt, 2), AT_pct = round(100 * (a + t) / acgt, 2),
             GC_skew = round((g - cc) / (g + cc), 3), stringsAsFactors = FALSE)
}

.sa_iupac_re <- function(site) {
  map <- c(A = "A", C = "C", G = "G", T = "T", R = "[AG]", Y = "[CT]", S = "[GC]", W = "[AT]",
           K = "[GT]", M = "[AC]", B = "[CGT]", D = "[AGT]", H = "[ACT]", V = "[ACG]", N = "[ACGT]")
  paste0(unname(map[strsplit(site, "")[[1]]]), collapse = "")
}

# All sites of the given enzymes (all the listed sites are palindromic: one strand is enough).
.sa_find_sites <- function(seq, enz) {
  rows <- lapply(seq_len(nrow(enz)), function(i) {
    re <- paste0("(?=", .sa_iupac_re(enz$Site[i]), ")")
    m <- gregexpr(re, seq, perl = TRUE)[[1]]
    if (m[1] == -1) return(NULL)
    data.frame(Enzyme = enz$Enzyme[i], Site = enz$Site[i], Start = as.integer(m),
               CutAfter = as.integer(m) + enz$Cut[i] - 1L, stringsAsFactors = FALSE)
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (!length(rows)) return(NULL)
  do.call(rbind, rows)
}
