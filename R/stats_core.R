# Sangerly - statistics and export helpers (no Shiny):
#   * alignment site classes like MEGA (conserved, variable, parsimony-informative, singleton)
#   * genetic distance tables (matrix / pairs) and a heat map (lower triangle, numbers in the cells)
#   * Excel export (package writexl)

.sg_xlsx <- function(sheets, file) {
  if (!requireNamespace("writexl", quietly = TRUE))
    stop("the 'writexl' package is not installed: run install.packages(\"writexl\") and build again")
  writexl::write_xlsx(sheets, file)
}

# ---- alignment sites -------------------------------------------------------------------------------------
# M: character matrix (sequences x columns). One row per alignment column.
# complete_deletion = TRUE : a column that has a gap or an ambiguous base in any sequence is left out (as MEGA's
#                            "complete deletion"); FALSE: gaps / ambiguous bases are ignored and the column is
#                            classified from the A, C, G, T that are present.
.site_stats <- function(M, complete_deletion = TRUE) {
  M <- toupper(M)
  n <- nrow(M); L <- ncol(M)
  bases <- c("A", "C", "G", "T")
  cnt <- vapply(bases, function(b) as.numeric(colSums(M == b)), numeric(L))
  if (is.null(dim(cnt))) cnt <- matrix(cnt, nrow = L, dimnames = list(NULL, bases))
  gaps <- as.numeric(colSums(M == "-"))
  nbase <- rowSums(cnt)
  other <- n - nbase - gaps                               # N and IUPAC codes
  types <- rowSums(cnt > 0)
  informative <- rowSums(cnt >= 2) >= 2
  top <- apply(cnt, 1, max)
  topbase <- ifelse(nbase > 0, bases[max.col(cnt, ties.method = "first")], NA_character_)
  excluded <- if (complete_deletion) (gaps > 0 | other > 0) else (nbase == 0)
  excl_label <- if (complete_deletion) "Excluded (gap or missing)" else "No base"
  cls <- ifelse(gaps == n, "All gaps",
         ifelse(excluded, excl_label,
         ifelse(types <= 1, "Conserved",
         ifelse(informative, "Variable, parsimony-informative", "Variable, singleton"))))
  states <- apply(cnt > 0, 1, function(v) paste(bases[v], collapse = "/"))
  data.frame(Site = seq_len(L), Class = cls, Most_common_base = topbase,
             A = cnt[, "A"], C = cnt[, "C"], G = cnt[, "G"], T = cnt[, "T"],
             Gaps = gaps, Other = other, Sequences_with_base = nbase,
             Most_common_pct = ifelse(nbase > 0, round(100 * top / nbase, 1), NA_real_),
             States = states, stringsAsFactors = FALSE, row.names = NULL)
}

.site_summary <- function(df) {
  used <- !(df$Class %in% c("All gaps", "Excluded (gap or missing)", "No base"))
  nU <- sum(used)
  pct <- function(k) if (nU > 0) round(100 * k / nU, 2) else NA_real_
  cons <- sum(df$Class == "Conserved")
  pinf <- sum(df$Class == "Variable, parsimony-informative")
  sing <- sum(df$Class == "Variable, singleton")
  data.frame(
    Measure = c("Sites in the alignment", "Sites used", "Sites left out (gaps / missing / no base)",
                "Conserved sites (C)", "Variable sites (V)", "Parsimony-informative sites (Pi)", "Singleton sites (S)"),
    Count = c(nrow(df), nU, nrow(df) - nU, cons, pinf + sing, pinf, sing),
    Percent_of_sites_used = c(NA, NA, NA, pct(cons), pct(pinf + sing), pct(pinf), pct(sing)),
    stringsAsFactors = FALSE)
}

.site_composition <- function(df) {
  used <- !(df$Class %in% c("All gaps", "Excluded (gap or missing)", "No base"))
  tot <- colSums(df[used, c("A", "C", "G", "T"), drop = FALSE])
  s <- sum(tot)
  data.frame(Base = names(tot), Count = as.numeric(tot),
             Percent = if (s > 0) round(100 * as.numeric(tot) / s, 2) else NA_real_, stringsAsFactors = FALSE)
}

.site_colours <- c("Conserved" = "#2c7be5", "Variable, parsimony-informative" = "#e5484d",
                   "Variable, singleton" = "#f5a524", "Excluded (gap or missing)" = "#c9ced6",
                   "No base" = "#c9ced6", "All gaps" = "#7b8794")

# ---- genetic distances ---------------------------------------------------------------------------------------
.dist_matrix_df <- function(D) {
  D <- as.matrix(D)
  data.frame(Taxon = rownames(D), D, check.names = FALSE, row.names = NULL, stringsAsFactors = FALSE)
}

# one row per pair of sequences ("column" format)
.dist_pairs <- function(D) {
  D <- as.matrix(D); nm <- rownames(D)
  ix <- which(lower.tri(D), arr.ind = TRUE)
  data.frame(Sequence_1 = nm[ix[, 2]], Sequence_2 = nm[ix[, 1]], Distance = D[ix], stringsAsFactors = FALSE)
}

# Size of the heat map (inches): grows with the number of sequences, with the number of decimals
# (a cell must be wide enough for its number) and with the length of the names.
.dist_dims <- function(D, digits = 2, show_num = TRUE) {
  D <- as.matrix(D); n <- nrow(D)
  lab <- rownames(D); if (is.null(lab)) lab <- as.character(seq_len(n))
  if (is.null(digits) || is.na(digits)) digits <- 2
  num_pt <- 8; ax_pt <- 9
  cell <- if (isTRUE(show_num)) max(0.30, (digits + 2) * num_pt * 0.62 / 72 + 0.12) else 0.22
  lab_in <- max(nchar(lab), 1) * ax_pt * 0.55 / 72
  list(n = n, cell = cell, num_pt = num_pt, ax_pt = ax_pt, lab_in = lab_in,
       w = n * cell + lab_in + 0.9, h = n * cell + 0.75 * lab_in + 1.1)
}

# Lower-triangle heat map: blue = small distance, red = large, numbers in the cells, grey upper triangle,
# names on the right (rows) and on top (columns, tilted 45 degrees).
#   scale_by = "matrix": one colour scale for all cells; the legend shows the real minimum, middle and maximum
#                        distance of the data (the diagonal, always 0, is not counted).
#   scale_by = "row"   : every row is scaled from its own smallest to its own largest distance (diagonal not counted).
.dist_heatmap <- function(D, digits = 2, show_num = TRUE, scale_by = "matrix", title = NULL) {
  D <- as.matrix(D); n <- nrow(D)
  lab <- rownames(D); if (is.null(lab)) lab <- as.character(seq_len(n))
  if (is.null(digits) || is.na(digits)) digits <- 2
  dm <- .dist_dims(D, digits, show_num)
  off <- row(D) != col(D)
  z <- matrix(0, n, n)
  if (identical(scale_by, "row")) {
    for (i in seq_len(n)) {
      v <- D[i, off[i, ]]
      r <- if (length(v) && any(is.finite(v))) range(v, na.rm = TRUE) else c(0, 0)
      z[i, ] <- if (diff(r) > 0) (D[i, ] - r[1]) / diff(r) else 0
    }
    brk <- c(0, 1); leg <- c("row min", "row max")
  } else {
    r <- suppressWarnings(range(D[off], na.rm = TRUE)); if (!all(is.finite(r))) r <- c(0, 1)
    z[] <- if (diff(r) > 0) (D - r[1]) / diff(r) else 0
    brk <- c(0, 0.5, 1); leg <- formatC(r[1] + brk * diff(r), format = "f", digits = max(digits, 2))
  }
  z <- pmin(1, pmax(0, z)); z[!off] <- 0
  lo <- which(lower.tri(D, diag = TRUE), arr.ind = TRUE)
  up <- which(upper.tri(D), arr.ind = TRUE)
  dl <- data.frame(x = lo[, 2], y = lo[, 1], d = D[lo], z = z[lo])
  dl$dark <- dl$z < 0.25 | dl$z > 0.75
  du <- data.frame(x = up[, 2], y = up[, 1])
  p <- ggplot2::ggplot() +
    ggplot2::geom_tile(data = du, ggplot2::aes(x = x, y = y), fill = "grey75", colour = "white", linewidth = 0.2) +
    ggplot2::geom_tile(data = dl, ggplot2::aes(x = x, y = y, fill = z), colour = "white", linewidth = 0.2) +
    ggplot2::scale_fill_gradient2(low = "#0000FF", mid = "#FFFFFF", high = "#FF0000", midpoint = 0.5,
      limits = c(0, 1), breaks = brk, labels = leg, name = NULL) +
    ggplot2::scale_x_continuous(position = "top", breaks = seq_len(n), labels = lab, expand = c(0, 0)) +
    ggplot2::scale_y_reverse(breaks = seq_len(n), labels = lab, position = "right", expand = c(0, 0)) +
    ggplot2::coord_fixed(clip = "off") +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(), axis.title = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(angle = 45, hjust = 0, vjust = 0, size = dm$ax_pt),
                   axis.text.y = ggplot2::element_text(size = dm$ax_pt),
                   legend.position = "top", legend.justification = "left",
                   legend.key.width = grid::unit(2.4, "cm"), legend.key.height = grid::unit(0.3, "cm"),
                   legend.text = ggplot2::element_text(size = 9),
                   plot.margin = ggplot2::margin(4, 4, 4, 4))
  if (isTRUE(show_num)) {
    p <- p + ggplot2::geom_text(data = dl, ggplot2::aes(x = x, y = y, colour = dark,
                                  label = formatC(d, format = "f", digits = digits)),
                                size = dm$num_pt / 2.845, show.legend = FALSE) +
      ggplot2::scale_colour_manual(values = c("TRUE" = "white", "FALSE" = "black"))
  }
  if (!is.null(title)) p <- p + ggplot2::labs(title = title)
  p
}

# write a plot to a file: type = "svg", "jpg", "png" or "pdf"; D, digits and show_num give the size
.dist_save <- function(p, file, type, D, digits = 2, show_num = TRUE) {
  dm <- .dist_dims(D, digits, show_num); w <- dm$w; h <- dm$h
  res <- max(72, min(300, floor(30000 / max(w, h))))        # keep very large matrices below the pixel limit
  switch(type,
    svg = grDevices::svg(file, width = w, height = h),
    jpg = grDevices::jpeg(file, width = round(w * res), height = round(h * res), res = res, quality = 95),
    png = grDevices::png(file, width = round(w * res), height = round(h * res), res = res),
    pdf = grDevices::pdf(file, width = w, height = h),
    stop("unknown figure type"))
  on.exit(grDevices::dev.off(), add = TRUE)
  print(p)
  invisible(file)
}
