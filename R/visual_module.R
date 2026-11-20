# SangeRly - "Visualisation" tab
#   * a tree editor inside the program (no website needed):
#       - open a tree file in any common format (Newick, IQ-TREE .treefile / .contree, RAxML, NEXUS ...)
#       - layout: rectangular phylogram / cladogram, circular (fan), unrooted, radial
#       - root: on an outgroup, midpoint, or RIGHT-click a tip / node in the picture and choose "Reroot"
#       - names: replace by a column of your Excel table, add a second column, or edit any name by hand
#       - groups from an Excel / CSV table: symbols (shape + colour), colour blocks behind the names
#         (5 width modes, optional gradient), coloured names, each with its own legend
#       - symbols: before or after the names, aligned in a column, offset - / +
#       - branch support as numbers or as circles (size = support) with a legend
#       - ONE picture size, the same on the screen and in the downloaded files; the tree and the legends are
#         moved by holding the LEFT mouse button and dragging
#       - export: SVG, PNG, PDF and Newick
#   * the links to iTOL and ChiPlot and the Newick download stay in the tab "External viewers".
#
# Packages: ape, DT (already used by the program); readxl to open .xlsx files (CSV / TSV work without it):
#           install.packages("readxl")      optional: phangorn (midpoint rooting)

# ---------------------------------------------------------------------------------------------------
# Core (no Shiny)
# ---------------------------------------------------------------------------------------------------

.vz_or <- function(a, b) {
  if (is.null(a) || length(a) == 0L) return(b)
  if (length(a) == 1L && is.na(a)) return(b)
  a
}

.vz_pal <- list(
  distinct = c("#E6194B", "#3CB44B", "#4363D8", "#F58231", "#911EB4", "#42D4F4", "#F032E6", "#BFEF45", "#FABED4",
               "#469990", "#DCBEFF", "#9A6324", "#FFE119", "#800000", "#AAFFC3", "#808000", "#000075", "#A9A9A9",
               "#FFD8B1", "#000000"),
  okabe    = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999"),
  set2     = c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3"),
  bright   = c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD", "#8C564B", "#E377C2", "#7F7F7F", "#BCBD22", "#17BECF"),
  grey     = c("#000000", "#444444", "#777777", "#999999", "#BBBBBB", "#DDDDDD"))
.vz_pal_choices <- c("Distinct (20 colours)" = "distinct", "Colour-blind safe (Okabe-Ito)" = "okabe",
                     "Soft (Set2)" = "set2", "Bright" = "bright", "Greys" = "grey")

# symbol shapes (R plotting symbols); 99 = hidden (the category gets no symbol)
.vz_shapes <- c("Circle" = 21L, "Square" = 22L, "Diamond" = 23L, "Triangle up" = 24L, "Triangle down" = 25L,
                "Circle (open)" = 1L, "Square (open)" = 0L, "Triangle (open)" = 2L, "Diamond (open)" = 5L,
                "Plus" = 3L, "Cross" = 4L, "Asterisk" = 8L, "Hidden" = 99L)
.vz_default_shapes <- c(21L, 22L, 24L, 23L, 25L, 3L, 4L, 8L, 1L, 0L, 2L, 5L)

.vz_layouts <- list(
  phylogram = list(type = "phylogram", len = TRUE,  rect = TRUE),
  cladogram = list(type = "cladogram", len = FALSE, rect = TRUE),
  fan       = list(type = "fan",       len = TRUE,  rect = FALSE),
  fan_clado = list(type = "fan",       len = FALSE, rect = FALSE),
  unrooted  = list(type = "unrooted",  len = TRUE,  rect = FALSE),
  radial    = list(type = "radial",    len = TRUE,  rect = FALSE))
.vz_layout_choices <- c("Rectangular phylogram" = "phylogram", "Rectangular cladogram" = "cladogram",
                        "Circular phylogram (fan)" = "fan", "Circular cladogram" = "fan_clado",
                        "Unrooted" = "unrooted", "Radial" = "radial")

# how far a colour block reaches (the same choices as in iTOL, plus "behind each name only")
.vz_blk_modes <- c("Behind each name only" = "tip_label",
                   "From the branch tip to a uniform end" = "tip_end",
                   "Uniform start and end" = "uniform",
                   "From the LCA to the longest tip of the group" = "lca_max",
                   "From the LCA to a uniform end" = "lca_end")
# annotation layers: the types that can be added (as often as wanted)
.vz_layer_types <- c("Symbols" = "sym", "Colour blocks" = "blk", "Name colour" = "lc")
.vz_fonts <- c("Sans (Arial / Helvetica)" = "sans", "Serif (Times)" = "serif", "Monospace (Courier)" = "mono")
.vz_leg_what <- c("Symbol / colour and name" = "both", "Symbol / colour only" = "mark", "Name only (in the colour of the group)" = "name")
.vz_leg_positions <- c("Right of the tree (outside)" = "right", "Left of the tree (outside)" = "left",
                       "Inside, top right" = "tr", "Inside, top left" = "tl",
                       "Inside, bottom right" = "br", "Inside, bottom left" = "bl")
# lines from the tips to the aligned names (value = the R line type; 0 = no line)
.vz_guide_choices <- c("None" = "0", "Dotted" = "3", "Dashed" = "2", "Solid" = "1")

# names in the table are matched with the tree names after the same cleaning the Phylogeny tab applies
.vz_clean <- function(x) gsub("[^A-Za-z0-9._-]", "_", trimws(as.character(x)))

# the column of the table that holds the most of the tree names
.vz_guess_id <- function(df, tips) {
  sc <- vapply(names(df), function(n) sum(.vz_clean(df[[n]]) %in% tips), numeric(1))
  names(df)[which.max(sc)]
}

.vz_key <- function(x) gsub("[^A-Za-z0-9]", "", x)

# Read a tree file in any of the usual formats. The format is recognised from the CONTENT, not from the name:
#   Newick : .nwk .newick .tre .tree .txt, IQ-TREE .treefile and .contree, RAxML bipartitions / bestTree ...
#   NEXUS  : .nex .nexus .trees, MrBayes .con.tre, BEAST / PAUP* trees (with or without a translate table)
# A file may hold many trees; the result is a list of trees (at least 3 tips each). Comments in square brackets
# (for example [&prob=0.9]) are dropped if the first attempt to read the file fails because of them.
.vz_read_trees <- function(path) {
  txt <- tryCatch(readLines(path, warn = FALSE), error = function(e) character())
  if (!length(txt) || !any(nzchar(trimws(txt)))) stop("the file is empty")
  all_txt <- paste(txt, collapse = "\n")
  nexus <- grepl("^\\s*#NEXUS", all_txt, ignore.case = TRUE) || grepl("begin\\s+trees", all_txt, ignore.case = TRUE)
  read_one <- function(f) if (nexus) ape::read.nexus(f) else ape::read.tree(f)
  trees <- tryCatch(read_one(path), error = function(e) NULL)
  if (is.null(trees)) {
    tmp <- tempfile(fileext = if (nexus) ".nex" else ".tre")
    on.exit(unlink(tmp), add = TRUE)
    writeLines(gsub("\\[[^]]*\\]", "", all_txt), tmp)
    trees <- tryCatch(read_one(tmp), error = function(e) NULL)
  }
  if (is.null(trees)) stop(if (nexus) "this NEXUS file has no readable tree block" else "this is not a Newick or NEXUS tree file")
  lst <- if (inherits(trees, "multiPhylo")) lapply(seq_along(trees), function(i) trees[[i]]) else list(trees)
  lst <- Filter(function(t) inherits(t, "phylo") && ape::Ntip(t) >= 3, lst)
  if (!length(lst)) stop("no tree with at least 3 tips was found in the file")
  lst
}

# height of the tree in the units of the plot (used to size the gaps between tips, symbols and names)
.vz_height <- function(tr, use_len) {
  h <- tryCatch({
    if (use_len) max(ape::node.depth.edgelength(tr)[seq_len(ape::Ntip(tr))])
    else max(ape::node.depth(tr, method = 2)) - 1
  }, error = function(e) NA_real_)
  if (!is.finite(h) || h <= 0) 1 else h
}

# branch support of every internal node on a 0-100 scale (NA where the label is not a number);
# NULL if the tree has no node labels. Proportions (posterior probabilities, 0-1) are turned into percentages.
.vz_support <- function(tr) {
  nl <- tr$node.label
  if (is.null(nl) || length(nl) != tr$Nnode) return(NULL)
  num <- suppressWarnings(as.numeric(sub(".*/", "", as.character(nl))))
  fin <- is.finite(num)
  if (any(fin) && max(num[fin]) <= 1) num <- num * 100
  num
}

# Limits of the drawing area for the rectangular layouts, so that every name fits in the plot box (nothing is cut off).
# ape guesses the room for the names, and the guess can be too small (italic / bold / other fonts, long names). The real width
# of every name is measured (inches) and the right limit is solved so that the farthest name just fits.
# Call it after par(plt = ..., family = ...) is set. Returns c(left, right) or NULL (then ape decides).
.vz_rect_limits <- function(tr2, ly, use_len, a, lab, off, fnt, aligned, wide = FALSE) {
  n <- ape::Ntip(tr2)
  old <- grDevices::dev.cur()
  grDevices::pdf(NULL)
  lp <- tryCatch({
    do.call(ape::plot.phylo, list(tr2, type = ly$type, use.edge.length = use_len, show.tip.label = FALSE))
    env <- tryCatch(ape::.PlotPhyloEnv, error = function(e) get(".PlotPhyloEnv", envir = asNamespace("ape")))
    get("last_plot.phylo", envir = env)
  }, error = function(e) NULL)
  grDevices::dev.off()
  if (old > 1L) try(grDevices::dev.set(old), silent = TRUE)
  if (is.null(lp)) return(NULL)
  xx <- lp$xx[seq_len(n)]
  pin <- graphics::par("pin")
  if (length(pin) != 2L || any(!is.finite(pin)) || any(pin <= 0)) return(NULL)
  w_in <- graphics::strwidth(lab, units = "inches", cex = a$cex, font = fnt)
  if (wide) w_in <- rep(max(w_in), n)
  fr <- w_in / pin[1]                                            # share of the box width that every name needs
  if (!all(is.finite(fr)) || any(fr >= 0.9)) return(NULL)
  x0 <- min(0, lp$xx)
  xs <- if (aligned) rep(max(xx), n) else xx
  len <- max((xs + 1.5 * off - x0) / (1 - fr))                   # the name ends (plus a little room) at the right limit
  if (!is.finite(len) || len <= 0) return(NULL)
  c(x0, x0 + len * 1.003)
}

# Limits of the drawing area for circular, radial and unrooted trees, so that the tree AND every name fit in the picture.
# ape only guesses the room for the names from the NUMBER of letters, so long names (for example "Species | country | host")
# were cut off at the edge of the picture. Here the real width of every name is measured (in inches, so the result is
# the same on the screen and in SVG / PNG / PDF) and the limits are solved so that the farthest name just fits.
# Call it after par(plt = ...) is set (it uses the size of the plot box). Returns list(x, y) or NULL (then ape decides).
.vz_polar_limits <- function(tr2, ly, use_len, a, lab, off, fnt, is_fan, pad_u = 0, pad_in = 0, wide = FALSE) {
  n <- ape::Ntip(tr2)
  # 1) where the nodes are: this does not depend on the device, so draw once without names on a throw-away device
  old <- grDevices::dev.cur()
  grDevices::pdf(NULL)
  lp <- tryCatch({
    pa <- list(tr2, type = ly$type, use.edge.length = use_len, show.tip.label = FALSE)
    if (is_fan) { pa$rotate.tree <- a$rotate; pa$open.angle <- a$open }
    do.call(ape::plot.phylo, pa)
    env <- tryCatch(ape::.PlotPhyloEnv, error = function(e) get(".PlotPhyloEnv", envir = asNamespace("ape")))
    get("last_plot.phylo", envir = env)
  }, error = function(e) NULL)
  grDevices::dev.off()
  if (old > 1L) try(grDevices::dev.set(old), silent = TRUE)
  if (is.null(lp)) return(NULL)
  allx <- lp$xx; ally <- lp$yy
  xx <- allx[seq_len(n)]; yy <- ally[seq_len(n)]
  pin <- graphics::par("pin")                                   # size of the plot box, inches
  if (length(pin) != 2L || any(!is.finite(pin)) || any(pin <= 0)) return(NULL)
  w_in <- graphics::strwidth(lab, units = "inches", cex = a$cex, font = fnt)
  h_in <- graphics::strheight("Ag", units = "inches", cex = a$cex, font = fnt)
  if (wide) w_in <- rep(max(w_in), n)                          # colour blocks reach the same end for every name
  unrooted <- identical(ly$type, "unrooted")
  ctr <- if (unrooted) c(mean(range(xx)), mean(range(yy))) else c(allx[n + 1L], ally[n + 1L])
  th <- atan2(yy - ctr[2], xx - ctr[1]); ca <- cos(th); sa <- sin(th)
  # 2) room taken by the names when one inch is s units of the tree drawing
  extent <- function(s) {
    wu <- (w_in + pad_in) * s + pad_u
    hu <- h_in * s
    if (unrooted) {                                             # names are horizontal: to the right of the tips on the right side, else to the left
      right <- ca >= 0
      x0 <- xx + ifelse(right, off, -off); x1 <- x0 + ifelse(right, wu, -wu)
      list(x = c(pmin(x0, x1), pmax(x0, x1)), y = c(yy - hu / 2, yy + hu / 2))
    } else {                                                    # names run outwards along the radius
      x0 <- xx + ca * off; y0 <- yy + sa * off
      x1 <- x0 + ca * wu;  y1 <- y0 + sa * wu
      px <- abs(sa) * hu / 2; py <- abs(ca) * hu / 2
      list(x = c(pmin(x0, x1) - px, pmax(x0, x1) + px), y = c(pmin(y0, y1) - py, pmax(y0, y1) + py))
    }
  }
  s0 <- max(diff(range(allx)) / pin[1], diff(range(ally)) / pin[2], 1e-9)
  s <- s0
  for (it in seq_len(80)) {
    e <- extent(s)
    rx <- diff(range(c(e$x, allx))); ry <- diff(range(c(e$y, ally)))
    s2 <- max(rx / pin[1], ry / pin[2]) * 1.004
    if (!is.finite(s2) || s2 > 1e3 * s0) return(NULL)           # the names are wider than the picture: leave it to ape
    done <- abs(s2 - s) < 1e-5 * s
    s <- s2
    if (done) break
  }
  e <- extent(s)
  cx <- mean(range(c(e$x, allx))); cy <- mean(range(c(e$y, ally)))
  list(x = cx + c(-1, 1) * s * pin[1] / 2, y = cy + c(-1, 1) * s * pin[2] / 2)
}

# Draw the tree with base graphics (ape) and add blocks, symbols, support, scale bar and legends.
#   tr  : the tree (original tip names; this fixes the order of the tips)
#   lab : the names to show, one per tip, in the order of tr$tip.label
#   a   : list of settings (see plot_args() in the server)
#         a$plt = c(x1, x2, y1, y2): the box of the tree inside the picture, as fractions of the picture
#         a$W, a$H = size of the picture in pixels
# Positions that leave this function (legend boxes, nodes) are fractions of the picture: x to the right,
# y UPWARDS (0 = bottom), the same system that the browser script uses.
.vz_draw <- function(tr, lab, a) {
  ly <- .vz_layouts[[a$layout]]
  if (is.null(ly)) ly <- .vz_layouts$phylogram
  use_len <- ly$len && !is.null(tr$edge.length)
  n <- ape::Ntip(tr)
  H <- .vz_height(tr, use_len)
  rect <- isTRUE(ly$rect)
  is_fan <- identical(ly$type, "fan")
  lays <- if (is.null(a$layers)) list() else a$layers
  syl <- Filter(function(z) identical(z$type, "sym"), lays)  # symbol layers, in the order of the list
  bkl <- Filter(function(z) identical(z$type, "blk"), lays)  # colour block layers, in the order of the list
  has_sym <- length(syl) > 0L
  # Every symbol layer owns a "cell" next to the tip (where = tip) or next to the end of the name (where = end); the cells of
  # several layers follow one another. gap = width of the cell, off = shift of the symbol (negative = left / inwards).
  cell <- function(z) z$gap * H + max(0, z$off / 100 * H)
  inl <- Filter(function(z) identical(z$where, "tip"), syl)
  enl <- Filter(function(z) identical(z$where, "end"), syl)
  inl_al <- any(vapply(inl, function(z) isTRUE(z$align), logical(1)))
  st_in <- sum(vapply(inl, cell, numeric(1)))                # width of the cells before the names
  st_end <- sum(vapply(enl, cell, numeric(1)))               # width of the cells after the names
  gap1 <- if (length(enl)) enl[[1]]$gap * H else if (has_sym) syl[[1]]$gap * H else 0.03 * H
  # distance between a tip and the start of its name
  off <- if (length(inl)) st_in
         else if (has_sym || length(bkl)) max(0.01 * H, 0.5 * gap1)
         else 0.01 * H
  fnt <- 1L + (if (isTRUE(a$bold)) 1L else 0L) + (if (isTRUE(a$italic)) 2L else 0L)   # 1 plain, 2 bold, 3 italic, 4 bold italic
  fam <- .vz_or(a$family, "sans")
  tr2 <- tr
  tr2$tip.label <- lab
  tipcol <- if (is.null(a$tipcol)) "black" else a$tipcol
  aligned <- rect && (isTRUE(a$align) || inl_al)

  # the tree lives in a box inside the picture (the legends and the symbols after the names are outside the box)
  graphics::par(plt = a$plt, xpd = NA, family = fam)

  args <- list(tr2, type = ly$type, use.edge.length = use_len, show.tip.label = TRUE,
               tip.color = tipcol, font = fnt, cex = a$cex, label.offset = off,
               edge.width = a$lwd, edge.color = a$edge_col)
  if (aligned) args$align.tip.label <- as.integer(.vz_or(a$guide, 0L))   # the number is the line type; 0 = no line
  if (is_fan) { args$rotate.tree <- a$rotate; args$open.angle <- a$open }
  if (rect) {
    wide_r <- any(vapply(bkl, function(z) !identical(z$mode, "tip_label"), logical(1)))
    lim <- tryCatch(.vz_rect_limits(tr2, ly, use_len, a, lab, off, fnt, aligned, wide_r), error = function(e) NULL)
    if (!is.null(lim)) { args$x.lim <- lim; graphics::par(xaxs = "i") }
  }
  if (!rect) {
    # circular / radial / unrooted: fit the tree and the REAL width of every name in the picture (nothing is cut off)
    pad_u <- 0.5 * off + (if (length(enl)) st_end - 0.5 * enl[[length(enl)]]$gap * H else 0)
    pad_in <- if (length(enl)) 0.09 * max(vapply(enl, function(z) z$size, numeric(1))) else 0
    wide <- any(vapply(bkl, function(z) !identical(z$mode, "tip_label"), logical(1)))
    lim <- tryCatch(.vz_polar_limits(tr2, ly, use_len, a, lab, off, fnt, is_fan, pad_u, pad_in, wide), error = function(e) NULL)
    if (!is.null(lim)) { args$x.lim <- lim$x; args$y.lim <- lim$y; graphics::par(xaxs = "i", yaxs = "i") }
    graphics::par(xpd = NA)
  }
  do.call(ape::plot.phylo, args)

  env <- tryCatch(ape::.PlotPhyloEnv, error = function(e) get(".PlotPhyloEnv", envir = asNamespace("ape")))
  lp <- get("last_plot.phylo", envir = env)
  xx <- lp$xx; yy <- lp$yy
  tips <- seq_len(n)
  ndc_x <- function(x) graphics::grconvertX(x, "user", "ndc")
  ndc_y <- function(y) graphics::grconvertY(y, "user", "ndc")
  nx <- ndc_x(xx); ny <- ndc_y(yy)                            # every node, as fractions of the picture (for the right-click)
  cx <- xx[n + 1L]; cy <- yy[n + 1L]                          # the root (the centre of the circular layouts)
  lw <- graphics::strwidth(lab, units = "user", cex = a$cex, font = fnt)
  dx <- xx[tips] - cx; dy <- yy[tips] - cy
  th <- atan2(dy, dx)                                          # angle of every tip (circular layouts)
  dth <- (2 * pi - a$open * pi / 180) / n
  # position of the nodes along the main axis: x for rectangular trees, distance from the centre for circular ones
  p_all <- if (rect) xx else sqrt((xx - cx)^2 + (yy - cy)^2)
  ptip <- p_all[tips]
  lab_start <- (if (aligned) rep(max(ptip), n) else ptip) + off
  lab_end <- lab_start + lw
  common_end <- max(lab_end) + 0.5 * off
  tip_end <- if (rect) ndc_x(lab_end) else NULL               # where each name ends (for the right-click on a name)
  ty_sorted <- sort(ny[tips])
  half_row <- if (length(ty_sorted) > 1L) 0.5 * stats::median(diff(ty_sorted)) else 0.02

  # colour blocks behind the names (semi-transparent, drawn over the tree, so lines and names stay readable)
  for (bz in bkl) if (rect || is_fan) {
    mode <- .vz_or(bz$mode, "tip_label")
    ord <- tryCatch({ e <- ape::reorder.phylo(tr, "cladewise")$edge[, 2]; e[e <= n] }, error = function(e) tips)
    rl <- rle(as.character(bz$grp[ord]))                    # runs of tips of the same group, in the order of the picture
    rid <- rep(seq_along(rl$lengths), rl$lengths)
    run_of <- integer(n); run_of[ord] <- rid
    lca_p <- numeric(length(rl$lengths)); run_max <- numeric(length(rl$lengths))
    for (k in seq_along(lca_p)) {
      tk <- ord[rid == k]
      run_max[k] <- max(ptip[tk])
      lca_p[k] <- if (length(tk) < 2L) min(ptip[tk]) else {
        m <- tryCatch(ape::getMRCA(tr, tk), error = function(e) NA)
        if (is.null(m) || is.na(m)) min(ptip[tk]) else p_all[m]
      }
    }
    p0 <- switch(mode, uniform = rep(max(ptip), n), lca_max = , lca_end = lca_p[run_of], ptip)
    p1 <- switch(mode, tip_label = lab_end + 0.5 * off, lca_max = run_max[run_of] + 0.5 * off, rep(common_end, n))
    p1 <- pmax(p1, p0 + 0.02 * H)
    nsl <- if (isTRUE(bz$grad)) (if (rect) 30L else 14L) else 1L
    for (s in seq_len(nsl)) {
      f0 <- (s - 1) / nsl; f1 <- s / nsl; u <- (f0 + f1) / 2
      al <- if (nsl == 1L) bz$alpha else bz$alpha * (if (isTRUE(bz$rev)) 0.1 + 0.9 * u else 1 - 0.9 * u)
      cs <- grDevices::adjustcolor(bz$col, alpha.f = al)
      ra <- p0 + f0 * (p1 - p0); rb <- p0 + f1 * (p1 - p0)
      if (rect) {
        graphics::rect(ra, yy[tips] - 0.46, rb, yy[tips] + 0.46, col = cs, border = NA)
      } else {
        for (i in tips) {
          tt <- seq(th[i] - 0.46 * dth, th[i] + 0.46 * dth, length.out = 5)
          graphics::polygon(cx + c(ra[i] * cos(tt), rev(rb[i] * cos(tt))), cy + c(ra[i] * sin(tt), rev(rb[i] * sin(tt))),
                            col = cs[i], border = NA)
        }
      }
    }
  }

  # symbols (drawn outside the box where needed)
  graphics::par(xpd = NA)
  if (has_sym) {
    if (rect || is_fan) { base_in <- ptip; base_en <- lab_end } else { base_in <- rep(0, n); base_en <- off + lw }
    if (!rect && !is_fan) {                                    # unrooted / radial: symbols lie on the line from the centre through the tip
      ux <- xx[tips] - mean(xx); uy <- yy[tips] - mean(yy)
      dd <- sqrt(ux^2 + uy^2); dd[dd == 0] <- 1
    }
    st_i <- 0                                                  # width of the cells before the names that are already used
    edge <- base_en                                            # where the next cell after the names starts, for every tip
    for (z in syl) {
      gz <- z$gap * H; oz <- z$off / 100 * H
      if (identical(z$where, "tip")) {
        start <- if (isTRUE(z$align) || (rect && inl_al)) rep(max(base_in), n) else base_in
        p <- start + st_i + 0.5 * gz + oz
        st_i <- st_i + gz + max(0, oz)
      } else {
        start <- if (isTRUE(z$align)) rep(max(edge), n) else edge
        p <- start + 0.5 * gz + oz
        edge <- start + gz + max(0, oz)
      }
      if (rect) { sx <- p; sy <- yy[tips] }
      else if (is_fan) { sx <- cx + p * cos(th); sy <- cy + p * sin(th) }
      else { sx <- xx[tips] + p * ux / dd; sy <- yy[tips] + p * uy / dd }
      keep <- which(z$pch != 99L)
      if (length(keep)) graphics::points(sx[keep], sy[keep], pch = z$pch[keep], col = z$col[keep], bg = z$col[keep], cex = z$size)
    }
  }

  # branch support: numbers or circles (size grows with the support)
  sup <- .vz_support(tr)
  if (!is.null(sup) && !identical(a$sup_style, "none")) {
    nodes <- n + seq_along(sup)
    ok <- is.finite(sup) & sup >= a$min_sup
    if (identical(a$sup_style, "text")) {
      nl <- as.character(tr$node.label); nl[!ok] <- ""
      if (any(nzchar(nl))) ape::nodelabels(nl, frame = "none", adj = c(1.15, -0.35), cex = a$cex * 0.75, col = "#c0392b")
    } else if (any(ok)) {
      graphics::points(xx[nodes[ok]], yy[nodes[ok]], pch = 21, col = "#222222", bg = a$sup_col,
                       cex = a$sup_cex * (0.35 + 0.65 * pmin(100, sup[ok]) / 100))
    }
  }
  if (isTRUE(a$scale_bar) && identical(ly$type, "phylogram") && use_len) ape::add.scale.bar(cex = a$cex * 0.8)
  if (nzchar(a$title))
    graphics::text(graphics::grconvertX(0.5, "ndc", "user"), graphics::grconvertY(1 - 16 / a$H, "ndc", "user"),
                   a$title, font = 2, cex = 1.2, xpd = NA)

  # the node selected with the right mouse button
  if (!is.null(a$sel) && a$sel >= 1 && a$sel <= length(xx))
    graphics::points(xx[a$sel], yy[a$sel], pch = 21, cex = 2.2, col = "#e5484d", bg = NA, lwd = 2)

  # legends: stacked from the chosen corner / margin; a legend that was dragged keeps its own place.
  # All places are fractions of the picture (x to the right, y upwards).
  leg_boxes <- list()
  if (length(a$legends)) {
    W <- a$W; Hp <- a$H
    bx0 <- a$plt[1]; bx1 <- a$plt[2]; by0 <- a$plt[3]; by1 <- a$plt[4]
    an <- switch(a$leg_pos,
      left  = list(x = 10 / W,                          xj = 0, y = by1,         yj = 1, dir = -1),
      tr    = list(x = bx1 - 0.01,                      xj = 1, y = by1 - 0.01,  yj = 1, dir = -1),
      tl    = list(x = bx0 + 0.01,                      xj = 0, y = by1 - 0.01,  yj = 1, dir = -1),
      br    = list(x = bx1 - 0.01,                      xj = 1, y = by0 + 0.01,  yj = 0, dir = 1),
      bl    = list(x = bx0 + 0.01,                      xj = 0, y = by0 + 0.01,  yj = 0, dir = 1),
      list(x = bx1 + (.vz_or(a$sym_out, 0) + 14) / W,   xj = 0, y = by1,         yj = 1, dir = -1))
    y <- an$y
    for (L in a$legends) {
      xy <- a$leg_xy[[L$id]]
      if (!is.null(xy)) { xn <- xy[1]; yn <- xy[2]; xj <- 0; yj <- 1 }
      else { xn <- an$x; yn <- y; xj <- an$xj; yj <- an$yj }
      what <- .vz_or(L$what, "both")
      graphics::par(family = .vz_or(L$family, "sans"))
      la <- list(x = graphics::grconvertX(xn, "ndc", "user"), y = graphics::grconvertY(yn, "ndc", "user"),
                 legend = if (identical(what, "mark")) rep("", length(L$labels)) else L$labels,
                 title = if (is.null(L$title) || !nzchar(L$title)) NULL else L$title,
                 bty = "n", cex = .vz_or(L$cex, a$leg_cex), xjust = xj, yjust = yj, title.adj = 0, y.intersp = 1.1,
                 text.font = .vz_or(L$font, 1L))
      if (!identical(what, "name")) { la$pch <- L$pch; la$col <- L$col; la$pt.bg <- L$bg; la$pt.cex <- L$ptcex }
      else la$text.col <- L$col                                  # names only: the names carry the colours
      if ("title.font" %in% names(formals(graphics::legend))) la$title.font <- .vz_or(L$font, 1L)
      r <- do.call(graphics::legend, la)
      rc <- r$rect
      lx <- ndc_x(c(rc$left, rc$left + rc$w)); lyv <- ndc_y(c(rc$top - rc$h, rc$top))
      leg_boxes[[length(leg_boxes) + 1L]] <- list(id = L$id, x0 = min(lx), x1 = max(lx), y0 = min(lyv), y1 = max(lyv))
      if (is.null(xy)) y <- y + an$dir * (abs(diff(lyv)) + 0.012)
    }
  }
  invisible(list(n = n, nx = nx, ny = ny, rect = rect, tip_end = tip_end, half_row = half_row,
                 leg_boxes = leg_boxes,
                 tree_box = list(x0 = a$plt[1], x1 = a$plt[2], y0 = a$plt[3], y1 = a$plt[4])))
}

# ---------------------------------------------------------------------------------------------------
# UI helpers
# ---------------------------------------------------------------------------------------------------

# Look of the annotation layers (a list of small cards inside the "2. Groups" card)
.vz_lay_css <- function() {
  shiny::tags$style(shiny::HTML("
    .vz-lay-title{font-weight:600;margin:12px 0 4px}
    .vz-lay{border:1px solid #dfe5ec;border-radius:8px;margin:6px 0;background:#fff}
    .vz-lay>summary{list-style:none}
    .vz-lay>summary::-webkit-details-marker{display:none}
    .vz-lay-sum{display:flex;align-items:center;gap:6px;padding:5px 6px 5px 8px;cursor:pointer;font-size:13px;color:#1f3c5f;background:#f7f9fc;border-radius:8px}
    .vz-lay[open]>.vz-lay-sum{border-bottom:1px solid #e3e8ee;border-radius:8px 8px 0 0}
    .vz-lay-sum::before{content:'\\25B8';font-size:12px;transition:transform .15s ease;flex:none}
    .vz-lay[open]>.vz-lay-sum::before{transform:rotate(90deg)}
    .vz-lay-type{font-weight:600;flex:none}
    .vz-lay-col{font-weight:500;color:#52606d;flex:1;min-width:0;overflow:hidden;text-overflow:ellipsis;white-space:nowrap}
    .vz-lay-btns{display:flex;gap:1px;flex:none}
    .vz-lay-btn{border:1px solid transparent;background:transparent;border-radius:5px;width:24px;height:24px;line-height:1;padding:0;cursor:pointer;color:#52606d;font-size:11px}
    .vz-lay-btn:hover{background:#e6edf5;border-color:#c9d2db}
    .vz-lay-btn:disabled{opacity:.3;cursor:default}
    .vz-lay-btn.off{opacity:.35;text-decoration:line-through}
    .vz-lay-body{padding:8px 10px 4px}
    .vz-lay-sub{font-weight:600;font-size:11.5px;color:#52606d;margin:10px 0 4px;text-transform:uppercase;letter-spacing:.04em}
    .vz-lay-add{display:flex;gap:6px;align-items:flex-end;margin:8px 0 4px}
    .vz-lay-add .form-group{margin-bottom:0;flex:1}
  "))
}

# Layout of the Visualisation page: the settings (left) scroll by themselves and the picture (right) stays in view, so the page
# is never longer than the window. The heights are set in the browser from the free space under the top of the page.
.vz_fit <- function() {
  shiny::tagList(
    shiny::tags$style(shiny::HTML("
      @media (min-width: 768px){
        .vz-left,.vz-right{overflow-y:auto;overflow-x:hidden;overscroll-behavior:contain}
        .vz-left{padding-right:10px}
        .vz-left::-webkit-scrollbar,.vz-right::-webkit-scrollbar{width:9px}
        .vz-left::-webkit-scrollbar-thumb,.vz-right::-webkit-scrollbar-thumb{background:#c5ccd6;border-radius:6px}
      }")),
    shiny::tags$script(shiny::HTML("
(function(){
  if (window.__sgVzFit) return;
  window.__sgVzFit = true;
  var fitOn = true;
  try { fitOn = localStorage.getItem('sg_vz_fit') !== '0'; } catch(e) {}
  function fit(){
    var L = document.querySelector('.vz-left'), R = document.querySelector('.vz-right');
    if (!L || !R || !L.offsetParent) return;
    var pb = R.querySelector('.vz-plotbox');
    if (pb && !document.getElementById('vz_fitwrap')) {
      var lb = document.createElement('label');
      lb.id = 'vz_fitwrap';
      lb.style.cssText = 'display:flex;align-items:center;gap:6px;margin:0 0 6px;font-weight:500;font-size:12.5px;cursor:pointer';
      lb.innerHTML = '<input type=\"checkbox\" id=\"vz_fitbox\"> Fit the picture to the window (untick for the full size)';
      pb.parentElement.insertBefore(lb, pb);
      var cb = lb.querySelector('input');
      cb.checked = fitOn;
      cb.addEventListener('change', function(){
        fitOn = cb.checked;
        try { localStorage.setItem('sg_vz_fit', fitOn ? '1' : '0'); } catch(e) {}
        fit();
      });
    }
    if (window.innerWidth < 768) { L.style.height = ''; R.style.height = ''; if (pb) pb.style.maxHeight = '84vh'; return; }
    var emb = document.documentElement.classList.contains('sg-embed');
    var top = L.getBoundingClientRect().top + (window.pageYOffset || 0);
    var avail = Math.max(320, Math.floor(window.innerHeight - top - (emb ? 36 : 16)));
    L.style.height = avail + 'px';
    R.style.height = avail + 'px';
    if (pb && pb.offsetParent) {
      var last = pb.parentElement.lastElementChild;
      var below = Math.max(0, last.getBoundingClientRect().bottom - pb.getBoundingClientRect().bottom);
      var above = pb.getBoundingClientRect().top - R.getBoundingClientRect().top + R.scrollTop;
      var mh = Math.max(220, avail - above - below - 16);
      pb.style.maxHeight = mh + 'px';
      var img = pb.querySelector('img');
      if (img) {
        var w0 = parseInt(img.getAttribute('width')), h0 = parseInt(img.getAttribute('height'));
        if (fitOn && w0 > 0 && h0 > 0) {
          var k = Math.min(1, (pb.clientWidth - 4) / w0, (mh - 4) / h0);
          img.style.width = Math.floor(w0 * k) + 'px';
          img.style.height = Math.floor(h0 * k) + 'px';
        } else { img.style.width = ''; img.style.height = ''; }
      }
    }
  }
  window.addEventListener('resize', fit);
  setInterval(fit, 500);
  document.addEventListener('shown.bs.tab', function(){ setTimeout(fit, 50); }, true);
  if (window.jQuery) jQuery(document).on('shiny:value shiny:visualchange', function(){ setTimeout(fit, 30); });
})();"))
  )
}


# a colour box (browser colour picker); it reports its value to Shiny when the user picks a colour
.vz_color <- function(id_full, value, title = NULL) {
  shiny::tags$input(type = "color", id = id_full, value = value, title = title,
                    style = "width:34px;height:26px;padding:1px;border:1px solid #b8c2cc;border-radius:5px;cursor:pointer;background:#fff;flex:none",
                    onchange = sprintf("Shiny.setInputValue('%s', this.value)", id_full))
}

.vz_shape_select <- function(id_full, selected) {
  shiny::tags$select(id = id_full,
    style = "height:26px;border:1px solid #b8c2cc;border-radius:5px;background:#fff;font-size:12px;flex:none",
    onchange = sprintf("Shiny.setInputValue('%s', this.value)", id_full),
    lapply(seq_along(.vz_shapes), function(j)
      shiny::tags$option(value = unname(.vz_shapes[j]), names(.vz_shapes)[j],
                         selected = if (as.integer(selected) == unname(.vz_shapes[j])) NA else NULL)))
}

.vz_section <- function(title, ..., open = FALSE) {
  shiny::tags$details(open = if (open) NA else NULL,
    shiny::tags$summary(style = "cursor:pointer;font-weight:600;margin:8px 0", title), ...)
}

# Mouse handling on the picture (runs in the browser, so dragging is instant):
#   * hold the LEFT button on a legend or on the tree and drag: a dashed frame follows the mouse; when the button is
#     released the server redraws the picture with the legend / tree in the new place;
#   * the RIGHT button opens a small menu (reroot ...) at the mouse; a plain left click does nothing.
# The server sends the boxes of the legends and of the tree after every drawing ('<output>_layout') and the content
# of the menu ('<output>_menu'); the browser answers with the inputs '<ns>drag', '<ns>ctx' and '<ns>menu_act'.
.vz_js <- function(out_id, ns_prefix) {
  js <- "
(function(){
  var OUT = '__OUT__', NS = '__NS__';
  window.__sgVz = window.__sgVz || {};
  if (window.__sgVz[OUT]) return;
  window.__sgVz[OUT] = true;
  var lay = {legends: [], tree: null}, drag = null, ghost = null, menu = null;
  function send(name, val){ Shiny.setInputValue(NS + name, val, {priority: 'event'}); }
  function imgOf(e){
    var c = document.getElementById(OUT);
    if (!c || !e.target || e.target.tagName !== 'IMG' || !c.contains(e.target)) return null;
    return e.target;
  }
  function frac(img, e){
    var r = img.getBoundingClientRect();
    return {r: r, fx: (e.clientX - r.left) / r.width, fy: 1 - (e.clientY - r.top) / r.height};
  }
  function inBox(b, fx, fy){ return !!b && fx >= b.x0 && fx <= b.x1 && fy >= b.y0 && fy <= b.y1; }
  function hit(fx, fy){
    var L = lay.legends || [];
    for (var i = L.length - 1; i >= 0; i--) {
      if (inBox(L[i], fx, fy)) return {kind: 'legend', box: L[i], id: L[i].id};
    }
    if (inBox(lay.tree, fx, fy)) return {kind: 'tree', box: lay.tree};
    return null;
  }
  function closeMenu(){ if (menu) { menu.remove(); menu = null; } }
  function ghostAt(r, box, dx, dy){
    if (!ghost) {
      ghost = document.createElement('div');
      ghost.style.cssText = 'position:fixed;border:2px dashed #e5484d;background:rgba(229,72,77,0.08);pointer-events:none;z-index:99999;border-radius:3px';
      document.body.appendChild(ghost);
    }
    ghost.style.left = (r.left + box.x0 * r.width + dx) + 'px';
    ghost.style.top = (r.top + (1 - box.y1) * r.height + dy) + 'px';
    ghost.style.width = ((box.x1 - box.x0) * r.width) + 'px';
    ghost.style.height = ((box.y1 - box.y0) * r.height) + 'px';
  }
  function showMenu(m){
    closeMenu();
    menu = document.createElement('div');
    menu.style.cssText = 'position:fixed;z-index:100000;background:#fff;border:1px solid #c9d2db;border-radius:8px;box-shadow:0 6px 22px rgba(0,0,0,0.18);padding:6px;min-width:210px;max-width:340px;font-size:13px';
    var t = document.createElement('div');
    t.textContent = m.title || '';
    t.style.cssText = 'padding:4px 8px 6px;color:#52606d;border-bottom:1px solid #eef1f4;margin-bottom:4px;word-break:break-word';
    menu.appendChild(t);
    (m.items || []).forEach(function(it){
      var b = document.createElement('div');
      b.textContent = it.label;
      b.style.cssText = 'padding:6px 10px;cursor:pointer;border-radius:5px';
      b.onmouseenter = function(){ b.style.background = '#eef3f8'; };
      b.onmouseleave = function(){ b.style.background = 'transparent'; };
      b.onclick = function(ev){ ev.stopPropagation(); closeMenu(); send('menu_act', {act: it.act, id: it.id || null}); };
      menu.appendChild(b);
    });
    document.body.appendChild(menu);
    var w = menu.offsetWidth, h = menu.offsetHeight;
    menu.style.left = Math.max(4, Math.min(m.x, window.innerWidth - w - 6)) + 'px';
    menu.style.top = Math.max(4, Math.min(m.y, window.innerHeight - h - 6)) + 'px';
  }
  document.addEventListener('mousedown', function(e){
    if (menu && !menu.contains(e.target)) closeMenu();
    if (e.button !== 0) return;
    var img = imgOf(e);
    if (!img) return;
    var p = frac(img, e), h = hit(p.fx, p.fy);
    if (!h) return;
    e.preventDefault();
    drag = {h: h, r: p.r, sx: e.clientX, sy: e.clientY, moved: false};
  }, true);
  document.addEventListener('mousemove', function(e){
    if (drag) {
      var dx = e.clientX - drag.sx, dy = e.clientY - drag.sy;
      if (!drag.moved && Math.abs(dx) + Math.abs(dy) < 4) return;
      drag.moved = true;
      document.body.style.cursor = 'grabbing';
      ghostAt(drag.r, drag.h.box, dx, dy);
      return;
    }
    var img = imgOf(e);
    if (!img) return;
    var p = frac(img, e);
    img.style.cursor = hit(p.fx, p.fy) ? 'grab' : 'default';
  }, true);
  document.addEventListener('mouseup', function(e){
    if (!drag) return;
    var d = drag;
    drag = null;
    document.body.style.cursor = '';
    if (ghost) { ghost.remove(); ghost = null; }
    if (!d.moved || e.button !== 0) return;
    var dx = e.clientX - d.sx, dy = e.clientY - d.sy;
    if (d.h.kind === 'legend') {
      var b = d.h.box, w = b.x1 - b.x0, ht = b.y1 - b.y0;
      var nx = Math.min(Math.max(b.x0 + dx / d.r.width, 0), 1 - w);
      var ny = Math.min(Math.max(b.y1 - dy / d.r.height, ht), 1);
      send('drag', {kind: 'legend', id: d.h.id, x: nx, y: ny});
    } else {
      send('drag', {kind: 'tree', dx: dx, dy: dy});
    }
  }, true);
  document.addEventListener('dragstart', function(e){ if (imgOf(e)) e.preventDefault(); }, true);
  document.addEventListener('contextmenu', function(e){
    var img = imgOf(e);
    if (!img) return;
    e.preventDefault();
    closeMenu();
    var p = frac(img, e), h = hit(p.fx, p.fy);
    if (h && h.kind === 'legend') {
      showMenu({x: e.clientX, y: e.clientY, title: 'Legend',
                items: [{label: 'Put this legend back in its starting place', act: 'leg_reset_one', id: h.id}]});
      return;
    }
    send('ctx', {fx: p.fx, fy: p.fy, cx: e.clientX, cy: e.clientY});
  }, true);
  document.addEventListener('keydown', function(e){ if (e.key === 'Escape') closeMenu(); });
  (function wait(){
    if (!(window.Shiny && Shiny.addCustomMessageHandler)) return setTimeout(wait, 30);
    Shiny.addCustomMessageHandler(OUT + '_layout', function(m){ lay = m || lay; });
    Shiny.addCustomMessageHandler(OUT + '_menu', showMenu);
    // a colour box (<input type=color>) set from R, for example when a project is opened
    Shiny.addCustomMessageHandler('sg_setcolor', function(m){
      var el = document.getElementById(m.id); if (!el) return;
      el.value = m.value; Shiny.setInputValue(m.id, m.value);
    });
  })();
})();"
  js <- gsub("__OUT__", out_id, js, fixed = TRUE)
  js <- gsub("__NS__", ns_prefix, js, fixed = TRUE)
  shiny::tags$script(shiny::HTML(js))
}

visualUI <- function(id) {
  ns <- shiny::NS(id)
  dlb <- function(x, lab) shiny::downloadButton(ns(x), lab)
  site <- function(title, url, text) {
    .card(title,
      shiny::p(text),
      shiny::tags$a(href = url, target = "_blank", rel = "noopener noreferrer", class = "btn btn-primary",
                    shiny::icon("arrow-up-right-from-square"), paste(" Open", url)))
  }
  shiny::div(class = "container-fluid",
    .vz_js(ns("tree_plot"), ns("")),
    .vz_fit(),
    .vz_lay_css(), .sg_ui_head(),
    shiny::fluidRow(class = "sg-siderow",
      shiny::column(4, shiny::div(class = "vz-left",
        .card("Project",
          shiny::helpText("Save the tree, the table of groups, the colours and shapes, the edited names, the place of the tree and ",
                          "legends and all the settings in one file. Open it later to go on editing, without starting again."),
          shiny::div(style = "display:flex;flex-direction:column;gap:2px",
            dlb("dl_proj", "Save project (.json)"),
            shiny::fileInput(ns("proj_file"), NULL, buttonLabel = "Open project...", placeholder = "", accept = c(".json")))),
        .card("Tree",
          shiny::radioButtons(ns("src"), NULL,
            c("NJ tree (Phylogeny tab)" = "nj", "IQ-TREE tree (Phylogeny tab)" = "iq", "A tree file" = "file")),
          shiny::conditionalPanel(sprintf("input['%s'] == 'file'", ns("src")),
            shiny::fileInput(ns("nwk"), NULL, buttonLabel = "Open tree file...", placeholder = ""),
            shiny::uiOutput(ns("file_info")),
            shiny::helpText("Any common tree file: Newick (.nwk .newick .tre .tree), IQ-TREE (.treefile .contree), ",
                            "RAxML (bipartitions, bestTree), NEXUS (.nex .nexus .trees, MrBayes .con.tre, BEAST). ",
                            "The format is recognised from the content. A file with several trees lets you choose one.")),
          shiny::selectInput(ns("layout"), "Layout", .vz_layout_choices),
          shiny::conditionalPanel(
            sprintf("input['%s'] == 'fan' || input['%s'] == 'fan_clado'", ns("layout"), ns("layout")),
            shiny::sliderInput(ns("rotate"), "Rotate (degrees)", 0, 360, 0, 5),
            shiny::sliderInput(ns("open"), "Open angle (degrees)", 0, 180, 0, 5)),
          shiny::radioButtons(ns("root"), "Root",
            c("As computed" = "none", "On the outgroup" = "outgroup", "Midpoint (needs phangorn)" = "midpoint")),
          shiny::selectizeInput(ns("outgroup"), "Outgroup", choices = NULL, multiple = TRUE),
          shiny::uiOutput(ns("root_note")),
          shiny::helpText("Or right-click a tip or a node in the picture and choose 'Reroot'."),
          shiny::radioButtons(ns("order"), "Order of the branches",
            c("As computed" = "none", "Ladderize" = "up", "Ladderize (reverse)" = "down")),
          shiny::checkboxInput(ns("last_only"),
            "Keep only the last support value (IQ-TREE 'SH-aLRT/UFBoot' becomes 'UFBoot')", TRUE)),
        .card("Groups from an Excel / CSV table",
          shiny::fileInput(ns("meta"), NULL, buttonLabel = "Open table...", placeholder = "",
                           accept = c(".xlsx", ".xls", ".csv", ".tsv", ".txt")),
          shiny::uiOutput(ns("meta_ui")),
          shiny::selectInput(ns("pal"), "Default colours", .vz_pal_choices),
          shiny::div(class = "vz-lay-title", "Annotation layers"),
          shiny::helpText("Every layer uses one column of the table. Add as many layers as you like, also of the same type ",
                          "(for example two layers of symbols: one for the Location and one for the Group). ",
                          "Use the arrows to change the order, the eye to hide a layer and the cross to remove it."),
          shiny::uiOutput(ns("layers_ui")),
          shiny::div(class = "vz-lay-add",
            shiny::selectInput(ns("add_type"), "Add a layer", .vz_layer_types, width = "100%"),
            shiny::actionButton(ns("add_layer"), "Add", icon = shiny::icon("plus"), style = "margin-bottom:15px")),
          shiny::actionButton(ns("reset_cols"), "Reset group colours and shapes", icon = shiny::icon("rotate-left"))),
        .card("Appearance",
          .vz_section("Names and branches",
            shiny::checkboxInput(ns("underscore"), "Replace _ with a space in the names", TRUE),
            shiny::selectInput(ns("fam"), "Font of the names", .vz_fonts),
            shiny::div(style = "display:flex;gap:18px;flex-wrap:wrap",
              shiny::checkboxInput(ns("italic"), "Italic names", TRUE),
              shiny::checkboxInput(ns("bold"), "Bold names", FALSE)),
            shiny::sliderInput(ns("cex"), "Name size", 0.4, 3, 0.9, 0.1),
            shiny::sliderInput(ns("lwd"), "Branch width", 0.5, 4, 1.5, 0.5),
            shiny::checkboxInput(ns("align"), "Align the names (rectangular layouts)", FALSE),
            shiny::selectInput(ns("guide"), "Lines from the tips to aligned names", .vz_guide_choices, selected = "0"),
            shiny::checkboxInput(ns("scale_bar"), "Scale bar (rectangular phylogram)", TRUE),
            shiny::textInput(ns("title"), "Title (optional)", "", width = "100%"), open = TRUE),
          .vz_section("Picture and tree size",
            shiny::uiOutput(ns("size_txt")),
            shiny::checkboxInput(ns("auto_size"), "The picture follows the tree (untick to set the picture size by hand)", TRUE),
            shiny::conditionalPanel(sprintf("!input['%s']", ns("auto_size")),
              shiny::numericInput(ns("man_w"), "Picture width (pixels)", 900, min = 300, max = 8000, step = 50),
              shiny::numericInput(ns("man_h"), "Picture height (pixels)", 700, min = 300, max = 12000, step = 50)),
            shiny::sliderInput(ns("zoom_x"), "Tree width (%)", 20, 300, 100, 5),
            shiny::sliderInput(ns("zoom_y"), "Tree height (%)", 20, 300, 100, 5),
            shiny::helpText("With 'the picture follows the tree', 100% is the normal size of the tree and the picture grows or ",
                            "shrinks with it. With a picture size set by hand, 100% means the tree fills the picture, and a ",
                            "smaller value leaves room around it. The downloaded files have exactly this size."),
            shiny::actionButton(ns("tree_reset"), "Put the tree back in its place", icon = shiny::icon("rotate-left"))),
          .vz_section("Branch support",
            shiny::selectInput(ns("sup_style"), "Show the support as",
              c("Numbers" = "text", "Circles (size = support)" = "circle", "Nothing" = "none")),
            shiny::sliderInput(ns("min_sup"), "Hide support below", 0, 100, 0, 5),
            shiny::conditionalPanel(sprintf("input['%s'] == 'circle'", ns("sup_style")),
              shiny::sliderInput(ns("sup_cex"), "Circle size", 0.4, 3, 1.2, 0.1),
              shiny::div(style = "display:flex;align-items:center;gap:8px;margin-bottom:8px",
                .vz_color(ns("sup_col"), "#555555", "Circle colour"), shiny::span("Circle colour")),
              shiny::checkboxInput(ns("sup_leg"), "Legend for the circles", TRUE))),
          .vz_section("Legends",
            shiny::checkboxInput(ns("show_leg"), "Show the legends", TRUE),
            shiny::sliderInput(ns("leg_cex"), "Legend text size", 0.5, 1.6, 0.9, 0.1),
            shiny::selectInput(ns("leg_pos"), "Starting position", .vz_leg_positions),
            shiny::helpText("Then hold the left mouse button on a legend and drag it where you want. ",
                            "Right-click a legend to send it back to its starting place."),
            shiny::actionButton(ns("leg_reset"), "Reset the legend positions", icon = shiny::icon("rotate-left")))))),
      shiny::column(8, shiny::div(class = "vz-right",
        shiny::uiOutput(ns("msg")),
        shiny::tabsetPanel(id = ns("tabs"),
          shiny::tabPanel("Tree", shiny::br(),
            shiny::div(class = "vz-plotbox", style = "overflow:auto;max-height:84vh;border:1px solid #e3e8ee;border-radius:6px;background:#fff",
              shiny::plotOutput(ns("tree_plot"), width = "auto", height = "auto")),
            shiny::helpText("Left mouse button: hold it on a legend or on the tree and drag to move it. ",
                            "Right mouse button: click a tip or a node for a menu (reroot)."),
            shiny::uiOutput(ns("sel_ui")),
            shiny::div(style = "display:flex;gap:6px;flex-wrap:wrap;align-items:flex-end;margin:10px 0",
              dlb("dl_svg", "SVG"), dlb("dl_png", "PNG"), dlb("dl_pdf", "PDF"),
              dlb("dl_nwk", "Newick"), dlb("dl_nwk_lab", "Newick (names as shown)"),
              shiny::div(style = "width:150px",
                shiny::selectInput(ns("dl_dpi"), "PNG resolution", c("150 dpi" = "150", "300 dpi" = "300", "600 dpi" = "600"),
                                   selected = "300")))),
          shiny::tabPanel("Edit names", shiny::br(),
            shiny::helpText("Double-click a name in the second column to change it; press Ctrl+Enter to confirm. ",
                            "Changes appear in the tree at once. An empty name goes back to the automatic name."),
            DT::dataTableOutput(ns("names_tbl")),
            shiny::actionButton(ns("reset_edits"), "Undo all name edits", icon = shiny::icon("rotate-left"))),
          shiny::tabPanel("External viewers", shiny::br(),
            shiny::fluidRow(
              shiny::column(6, site("iTOL - Interactive Tree Of Life", "https://itol.embl.de/",
                "Display, annotate and export phylogenetic trees (Newick files).")),
              shiny::column(6, site("ChiPlot", "https://chiplot.online/",
                "Online tool for drawing and styling trees and other charts."))),
            .card("Newick files for the web viewers",
              shiny::uiOutput(ns("status")),
              shiny::div(style = "display:flex;gap:8px;flex-wrap:wrap",
                shiny::downloadButton(ns("dl_nj"), "NJ tree (Newick)"),
                shiny::downloadButton(ns("dl_iq"), "IQ-TREE tree (Newick)")),
              shiny::helpText("These are the trees exactly as computed (not rerooted here); they use the option 'Keep only the last ",
                              "support value' of the Tree card. Sequence names were cleaned in the Phylogeny tab (spaces and symbols ",
                              "became '_'). Please read the terms of use of each website before uploading data."))))))))
}

# ---------------------------------------------------------------------------------------------------
# Server
# ---------------------------------------------------------------------------------------------------

visualServer <- function(id, nj = shiny::reactive(NULL), iq = shiny::reactive(NULL)) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns
    rv <- shiny::reactiveValues(file_trees = NULL, meta = NULL, sheets = NULL, edits = character(),
                                epoch = 0L, tbl_epoch = 0L, sel = NULL, msg = NULL,
                                leg_xy = list(), tx = 0, ty = 0, layer_ov = list(),
                                lay = list(list(id = "sy", type = "sym"), list(id = "bk", type = "blk"), list(id = "lc", type = "lc")),
                                lay_n = 0L, vis = list())
    pend <- new.env()                   # values waiting for a control that does not exist yet (project opened): cols, outgroup, skip_*
    coords <- new.env()                 # nodes of the last drawing (for the right-click)
    memo <- new.env()                   # size of the picture when it follows the tree; open = layers shown open; made = layers with outputs
    memo$open <- "sy"; memo$made <- character()
    get_tree <- function(r) if (is.null(r)) NULL else r$tree

    # ---- the tree ---------------------------------------------------------------------------------------
    shiny::observeEvent(input$nwk, {
      res <- tryCatch(.vz_read_trees(input$nwk$datapath), error = function(e) conditionMessage(e))
      if (is.character(res)) {
        rv$file_trees <- NULL
        rv$msg <- list(txt = paste("The tree file could not be read:", res), type = "warning")
      } else {
        rv$file_trees <- res; rv$msg <- NULL
        shiny::updateRadioButtons(session, "src", selected = "file")
      }
    })
    output$file_info <- shiny::renderUI({
      tl <- rv$file_trees
      if (is.null(tl)) return(NULL)
      shiny::tagList(
        shiny::helpText(sprintf("%d tree%s read.", length(tl), if (length(tl) == 1L) "" else "s")),
        if (length(tl) > 1L)
          shiny::selectInput(ns("tree_pick"), "Tree to show",
            stats::setNames(seq_along(tl), sprintf("Tree %d of %d (%d tips)", seq_along(tl), length(tl),
                                                   vapply(tl, ape::Ntip, numeric(1)))),
            selected = shiny::isolate(input$tree_pick)))
    })
    src_tree <- shiny::reactive({
      switch(.vz_or(input$src, "nj"), nj = get_tree(nj()), iq = get_tree(iq()),
             file = {
               tl <- rv$file_trees
               if (is.null(tl)) NULL else tl[[min(length(tl), max(1L, as.integer(.vz_or(input$tree_pick, 1L))))]]
             })
    })
    base_tree <- shiny::reactive({
      tr <- src_tree(); shiny::req(tr)
      if (isTRUE(input$last_only) && !is.null(tr$node.label))
        tr$node.label <- sub(".*/", "", as.character(tr$node.label))
      tr
    })
    shiny::observe({
      tr <- tryCatch(base_tree(), error = function(e) NULL)
      sel <- shiny::isolate(input$outgroup)
      if (!is.null(tr) && !is.null(pend$outgroup) && identical(input$src, "file")) {     # outgroup of an opened project
        sel <- intersect(pend$outgroup, tr$tip.label); pend$outgroup <- NULL
      }
      shiny::updateSelectizeInput(session, "outgroup", choices = if (is.null(tr)) character() else tr$tip.label,
                                  selected = sel, server = FALSE)
    })

    rooted_res <- shiny::reactive({
      tr <- base_tree(); mode <- .vz_or(input$root, "none"); note <- NULL
      if (identical(mode, "outgroup")) {
        og <- intersect(input$outgroup, tr$tip.label)
        if (!length(og)) note <- "Choose the outgroup (or right-click a tip or node in the picture)."
        else if (length(og) >= ape::Ntip(tr)) note <- "The outgroup cannot be all the tips."
        else {
          r <- tryCatch(ape::root(tr, outgroup = og, resolve.root = TRUE, edgelabel = TRUE),
                        error = function(e) tryCatch(ape::root(tr, outgroup = og, resolve.root = TRUE), error = function(e2) NULL))
          if (is.null(r)) note <- "This outgroup could not be used (the chosen tips do not form one group in the tree); the tree is left as computed."
          else return(list(tree = r, note = NULL))
        }
      } else if (identical(mode, "midpoint")) {
        if (!requireNamespace("phangorn", quietly = TRUE)) note <- "Midpoint rooting needs the 'phangorn' package: install.packages(\"phangorn\")"
        else {
          r <- tryCatch(phangorn::midpoint(tr), error = function(e) NULL)
          if (is.null(r)) note <- "Midpoint rooting failed for this tree."
          else return(list(tree = r, note = NULL))
        }
      }
      list(tree = tr, note = note)
    })
    drawn <- shiny::reactive({
      tr <- rooted_res()$tree
      o <- .vz_or(input$order, "none")
      if (identical(o, "up")) tr <- ape::ladderize(tr, right = TRUE)
      else if (identical(o, "down")) tr <- ape::ladderize(tr, right = FALSE)
      tr
    })
    output$root_note <- shiny::renderUI({
      nt <- tryCatch(rooted_res()$note, error = function(e) NULL)
      if (!is.null(nt)) shiny::div(class = "alert alert-warning", style = "padding:6px 10px;font-size:13px", nt)
    })
    shiny::observeEvent(drawn(), { rv$sel <- NULL }, ignoreInit = TRUE)

    # ---- table with the groups (Excel / CSV) ----------------------------------------------------------------
    shiny::observeEvent(list(input$meta, input$sheet), {
      f <- input$meta
      if (is.null(f)) return()
      ext <- tolower(tools::file_ext(f$name))
      sheets <- NULL
      df <- tryCatch({
        if (ext %in% c("xlsx", "xls", "xlsm")) {
          if (!requireNamespace("readxl", quietly = TRUE))
            stop("the 'readxl' package is needed to open Excel files: install.packages(\"readxl\") (or save the table as CSV)")
          sheets <- readxl::excel_sheets(f$datapath)
          sh <- input$sheet
          if (is.null(sh) || !(sh %in% sheets)) sh <- sheets[1]
          as.data.frame(readxl::read_excel(f$datapath, sheet = sh, col_types = "text"), stringsAsFactors = FALSE)
        } else if (ext %in% c("tsv", "txt")) {
          utils::read.delim(f$datapath, stringsAsFactors = FALSE, check.names = FALSE, colClasses = "character")
        } else {
          utils::read.csv(f$datapath, stringsAsFactors = FALSE, check.names = FALSE, colClasses = "character")
        }
      }, error = function(e) {
        rv$msg <- list(txt = paste("The table could not be read:", conditionMessage(e)), type = "danger")
        NULL
      })
      if (is.null(df)) return()
      nm <- names(df)
      bad <- is.na(nm) | !nzchar(trimws(nm))
      nm[bad] <- paste0("Column", seq_along(nm))[bad]
      names(df) <- make.unique(nm)
      df[] <- lapply(df, function(x) trimws(as.character(x)))
      rv$meta <- df
      rv$sheets <- sheets
      rv$msg <- NULL
    })
    output$meta_ui <- shiny::renderUI({
      m <- rv$meta
      if (is.null(m)) return(shiny::helpText("Optional. One row per sequence; the first row holds the column names (for example ",
                                             "ID, Species, Location, Group). The sequence names must be the same as in the tree."))
      cols <- names(m)
      tr <- tryCatch(shiny::isolate(base_tree()), error = function(e) NULL)
      guess <- if (is.null(tr)) cols[1] else .vz_guess_id(m, tr$tip.label)
      ch <- c("(none)", cols)
      pk <- function(id, choices, default = NULL) {         # a choice saved in a project wins while the project is being opened
        v <- pend$cols[[id]]
        if (!is.null(v) && length(v) == 1L && v %in% choices) v else default
      }
      shiny::onFlushed(function() pend$cols <- NULL, once = TRUE)
      shiny::tagList(
        if (length(rv$sheets) > 1) shiny::selectInput(ns("sheet"), "Sheet", rv$sheets, selected = shiny::isolate(input$sheet)),
        shiny::selectInput(ns("id_col"), "Column with the sequence names", cols, selected = pk("id_col", cols, guess)),
        shiny::uiOutput(ns("match_info")),
        shiny::selectInput(ns("lab_col"), "Names shown in the tree: from column", c("(name in the tree)", cols),
                           selected = pk("lab_col", c("(name in the tree)", cols))),
        shiny::selectInput(ns("lab_app"), "Add after the name", ch, selected = pk("lab_app", ch)))
    })
    match_idx <- shiny::reactive({
      m <- rv$meta; shiny::req(m)
      idc <- input$id_col; shiny::req(idc, idc %in% names(m))
      match(drawn()$tip.label, .vz_clean(m[[idc]]))
    })
    aligned <- shiny::reactive({            # one row of the table per tip, in the order of the tips (NA row = not found)
      m <- rv$meta
      if (is.null(m)) return(NULL)
      idx <- tryCatch(match_idx(), error = function(e) NULL)
      if (is.null(idx)) return(NULL)
      m[idx, , drop = FALSE]
    })
    output$match_info <- shiny::renderUI({
      idx <- tryCatch(match_idx(), error = function(e) NULL)
      if (is.null(idx)) return(NULL)
      k <- sum(!is.na(idx)); n <- length(idx)
      if (k == n) return(shiny::div(class = "alert alert-success", style = "padding:6px 10px;font-size:13px",
                                    sprintf("All %d sequences of the tree were found in the table.", n)))
      miss <- drawn()$tip.label[is.na(idx)]
      shiny::div(class = "alert alert-warning", style = "padding:6px 10px;font-size:13px",
                 sprintf("%d of %d sequences were found in the table. Not found: %s%s. Choose another name column or check the spelling.",
                         k, n, paste(utils::head(miss, 5), collapse = ", "), if (length(miss) > 5) ", ..." else ""))
    })

    # ---- groups: levels, colours, shapes -----------------------------------------------------------------------
    col_vals <- function(col) {
      al <- aligned()
      if (is.null(al) || is.null(col) || !(col %in% names(al))) return(NULL)
      v <- trimws(as.character(al[[col]]))
      v[is.na(v) | !nzchar(v)] <- "(no data)"
      v
    }
    col_levels <- function(col, vals) {      # order of the table, then anything else; "(no data)" last
      lv <- unique(trimws(as.character(rv$meta[[col]])))
      lv <- lv[!is.na(lv) & nzchar(lv)]
      lv <- lv[lv %in% vals]
      lv <- c(lv, setdiff(unique(vals), c(lv, "(no data)")))
      if ("(no data)" %in% vals) lv <- c(lv, "(no data)")
      lv
    }
    layer_defaults <- function(prefix, col) {
      if (is.null(col) || identical(col, "(none)")) return(NULL)
      v <- col_vals(col)
      if (is.null(v)) return(NULL)
      lv <- col_levels(col, v)
      pal <- .vz_pal[[.vz_or(input$pal, "distinct")]]
      nd <- lv == "(no data)"
      cols <- pal[(seq_along(lv) - 1L) %% length(pal) + 1L]; cols[nd] <- "#BDBDBD"
      pch <- .vz_default_shapes[(seq_along(lv) - 1L) %% length(.vz_default_shapes) + 1L]; pch[nd] <- 99L
      ov <- rv$layer_ov[[prefix]]            # colours and shapes of an opened project (same column, palette and reset counter)
      if (!is.null(ov) && identical(ov$col, col) && identical(ov$pal, .vz_or(input$pal, "distinct")) && identical(ov$epoch, rv$epoch)) {
        h <- match(lv, ov$lv)
        cols[!is.na(h)] <- ov$cols[h[!is.na(h)]]
        pch[!is.na(h)] <- ov$pch[h[!is.na(h)]]
      }
      list(lv = lv, vals = v, cols = cols, pch = as.integer(pch), col = col,
           key = paste(.vz_key(col), .vz_or(input$pal, "distinct"), rv$epoch, sep = "_"))
    }
    layer_map <- function(prefix, col) {       # defaults, overridden by what the user picked
      d <- layer_defaults(prefix, col)
      if (is.null(d)) return(NULL)
      for (i in seq_along(d$lv)) {
        d$cols[i] <- .vz_or(input[[sprintf("%s_c_%s_%d", prefix, d$key, i)]], d$cols[i])
        d$pch[i] <- as.integer(.vz_or(input[[sprintf("%s_s_%s_%d", prefix, d$key, i)]], d$pch[i]))
      }
      d
    }
    layer_rows <- function(prefix, col, with_shape) {
      d <- layer_defaults(prefix, col)
      if (is.null(d)) return(shiny::helpText("Choose a column above."))
      shiny::tagList(lapply(seq_along(d$lv), function(i) {
        cid <- sprintf("%s_c_%s_%d", prefix, d$key, i)
        sid <- sprintf("%s_s_%s_%d", prefix, d$key, i)
        cur_c <- shiny::isolate(.vz_or(input[[cid]], d$cols[i]))
        cur_s <- shiny::isolate(.vz_or(input[[sid]], d$pch[i]))
        shiny::div(style = "display:flex;align-items:center;gap:8px;margin:3px 0",
          .vz_color(ns(cid), cur_c, d$lv[i]),
          if (with_shape) .vz_shape_select(ns(sid), cur_s),
          shiny::span(style = "font-size:13px", d$lv[i]))
      }))
    }
    shiny::observeEvent(input$reset_cols, { rv$epoch <- rv$epoch + 1L })

    # ---- annotation layers: any number, any type, the same type as often as wanted ---------------------------------
    # A layer is list(id, type). Its controls have the ids "<id>_col" (the column), "<id>_where" ... and, per group,
    # "<id>_c_<key>_<i>" (colour) and "<id>_s_<key>_<i>" (shape). The first three layers are "sy", "bk" and "lc".
    lid <- function(L, nm) paste0(L$id, "_", nm)
    lay_name <- function(type) { i <- match(type, .vz_layer_types); if (is.na(i)) type else names(.vz_layer_types)[i] }
    lay_in <- function(L, nm, default) {         # value of a control while the list is redrawn: what it shows, else the project, else the default
      v <- shiny::isolate(input[[lid(L, nm)]])
      if (is.null(v)) v <- pend$lay[[L$id]][[nm]]
      .vz_or(v, default)
    }
    lay_now <- function(L, nm, default) .vz_or(input[[lid(L, nm)]], .vz_or(pend$lay[[L$id]][[nm]], default))
    lay_ui <- function(L, k, n, ch, vis, open) {
      btn <- function(act, label, title, off = FALSE, dis = FALSE) shiny::tags$button(type = "button",
        class = paste("vz-lay-btn", if (off) "off" else ""), title = title, disabled = if (dis) NA else NULL,
        onclick = sprintf("event.preventDefault();event.stopPropagation();%sShiny.setInputValue('%s',{id:'%s',act:'%s',t:Date.now()},{priority:'event'})",
                          if (identical(act, "vis")) "this.classList.toggle('off');" else "", ns("lay_act"), L$id, act),
        shiny::HTML(label))
      head <- shiny::tags$summary(class = "vz-lay-sum",
        shiny::span(class = "vz-lay-type", lay_name(L$type)),
        shiny::span(class = "vz-lay-col", shiny::textOutput(ns(lid(L, "ttl")), inline = TRUE)),
        shiny::span(class = "vz-lay-btns",
          btn("vis", "&#128065;&#65038;", "Show / hide this layer", off = !vis),
          btn("up", "&#9650;", "Move this layer up", dis = k == 1L),
          btn("down", "&#9660;", "Move this layer down", dis = k == n),
          btn("del", "&#10005;", "Remove this layer")))
      sel <- lay_in(L, "col", "(none)"); if (!(sel %in% ch)) sel <- "(none)"
      sets <- switch(L$type,
        sym = shiny::tagList(
          shiny::radioButtons(ns(lid(L, "where")), "Position",
            c("Before the name (at the tip)" = "tip", "After the name (column at the right)" = "end"),
            selected = lay_in(L, "where", "tip")),
          shiny::checkboxInput(ns(lid(L, "align")), "Align the symbols in one column", isTRUE(lay_in(L, "align", FALSE))),
          shiny::numericInput(ns(lid(L, "off")), "Offset (- to the left, + to the right; % of the tree height)",
                              as.numeric(lay_in(L, "off", 0)), min = -100, max = 300, step = 1),
          shiny::sliderInput(ns(lid(L, "size")), "Symbol size", 0.4, 3, as.numeric(lay_in(L, "size", 1.2)), 0.1),
          shiny::sliderInput(ns(lid(L, "gap")), "Space around the symbols (% of the tree height)", 0, 10,
                             as.numeric(lay_in(L, "gap", 3)), 0.5)),
        blk = shiny::tagList(
          shiny::selectInput(ns(lid(L, "mode")), "Width of the blocks", .vz_blk_modes, selected = lay_in(L, "mode", "tip_label")),
          shiny::checkboxInput(ns(lid(L, "grad")), "Colour gradient (fades out)", isTRUE(lay_in(L, "grad", FALSE))),
          shiny::conditionalPanel(sprintf("input['%s']", ns(lid(L, "grad"))),
            shiny::checkboxInput(ns(lid(L, "rev")), "Reverse the gradient", isTRUE(lay_in(L, "rev", FALSE)))),
          shiny::sliderInput(ns(lid(L, "alpha")), "Colour strength", 0.1, 1, as.numeric(lay_in(L, "alpha", 0.35)), 0.05),
          shiny::helpText("Blocks are drawn in the rectangular and circular layouts. The 'LCA' widths start at the common ",
                          "ancestor of the neighbouring tips that belong to the same group.")),
        NULL)
      shiny::tags$details(class = "vz-lay", open = if (open) NA else NULL,
        ontoggle = sprintf("Shiny.setInputValue('%s',{id:'%s',open:this.open,t:Date.now()},{priority:'event'})", ns("lay_open"), L$id),
        head,
        shiny::div(class = "vz-lay-body",
          shiny::selectInput(ns(lid(L, "col")), "Column", ch, selected = sel),
          shiny::uiOutput(ns(lid(L, "levels"))),
          if (!is.null(sets)) shiny::tagList(shiny::div(class = "vz-lay-sub", "Settings of this layer"), sets),
          shiny::div(class = "vz-lay-sub", "Legend of this layer"),
          shiny::checkboxInput(ns(lid(L, "lg_show")), "Show this legend", isTRUE(lay_in(L, "lg_show", TRUE))),
          shiny::selectInput(ns(lid(L, "lg_what")), "Show in the legend", .vz_leg_what, selected = lay_in(L, "lg_what", "both")),
          shiny::checkboxInput(ns(lid(L, "lg_tt")), "Show the title", isTRUE(lay_in(L, "lg_tt", TRUE))),
          shiny::textInput(ns(lid(L, "lg_title")), "Title (empty = the name of the column)", as.character(lay_in(L, "lg_title", "")), width = "100%"),
          shiny::selectInput(ns(lid(L, "lg_fam")), "Font", .vz_fonts, selected = lay_in(L, "lg_fam", "sans")),
          shiny::div(style = "display:flex;gap:18px;flex-wrap:wrap",
            shiny::checkboxInput(ns(lid(L, "lg_bold")), "Bold", isTRUE(lay_in(L, "lg_bold", FALSE))),
            shiny::checkboxInput(ns(lid(L, "lg_ital")), "Italic", isTRUE(lay_in(L, "lg_ital", FALSE)))),
          shiny::sliderInput(ns(lid(L, "lg_cex")), "Text size (1 = the general legend size)", 0.5, 2.5,
                             as.numeric(lay_in(L, "lg_cex", 1)), 0.1)))
    }
    output$layers_ui <- shiny::renderUI({
      m <- rv$meta; lay <- rv$lay
      if (is.null(m)) return(shiny::helpText("Open a table above first: every layer takes the groups from one column of the table."))
      if (!length(lay)) return(shiny::helpText("No layer yet. Choose a type below and press Add."))
      ch <- c("(none)", names(m)); vis <- shiny::isolate(rv$vis); n <- length(lay)
      shiny::tagList(lapply(seq_len(n), function(k) {
        L <- lay[[k]]
        lay_ui(L, k, n, ch, !isFALSE(vis[[L$id]]), L$id %in% memo$open)
      }))
    })
    shiny::observe({                      # the group list and the title of every layer (made once per layer)
      for (L in rv$lay) if (!(L$id %in% memo$made)) {
        memo$made <- c(memo$made, L$id)
        local({
          LL <- L
          output[[lid(LL, "levels")]] <- shiny::renderUI(layer_rows(LL$id, lay_now(LL, "col", NULL), identical(LL$type, "sym")))
          output[[lid(LL, "ttl")]] <- shiny::renderText({
            cl <- lay_now(LL, "col", NULL)
            if (is.null(cl) || identical(cl, "(none)")) "" else paste0("- ", cl)
          })
        })
      }
    })
    shiny::observeEvent(input$lay_open, {
      o <- input$lay_open
      memo$open <- if (isTRUE(o$open)) union(memo$open, o$id) else setdiff(memo$open, o$id)
    })
    shiny::observeEvent(input$add_layer, {
      ty <- .vz_or(input$add_type, "sym")
      rv$lay_n <- rv$lay_n + 1L
      id <- paste0("ly", rv$lay_n)
      memo$open <- c(memo$open, id)
      rv$lay <- c(rv$lay, list(list(id = id, type = ty)))
    })
    shiny::observeEvent(input$lay_act, {
      a <- input$lay_act; id <- a$id; lay <- rv$lay
      k <- which(vapply(lay, function(L) identical(L$id, id), logical(1)))
      if (length(k) != 1L) return()
      if (identical(a$act, "vis")) {
        v <- rv$vis; v[[id]] <- isFALSE(v[[id]]); rv$vis <- v         # hidden -> shown, anything else -> hidden
      } else if (identical(a$act, "up") && k > 1L) {
        lay[c(k - 1L, k)] <- lay[c(k, k - 1L)]; rv$lay <- lay
      } else if (identical(a$act, "down") && k < length(lay)) {
        lay[c(k, k + 1L)] <- lay[c(k + 1L, k)]; rv$lay <- lay
      } else if (identical(a$act, "del")) {
        lay[[k]] <- NULL; rv$lay <- lay
        xy <- rv$leg_xy; xy[[id]] <- NULL; rv$leg_xy <- xy
      }
    })
    layers_all <- shiny::reactive({       # every layer that has data, with its settings (hidden ones too)
      out <- list(); vis <- rv$vis
      for (L in rv$lay) {
        d <- layer_map(L$id, lay_now(L, "col", NULL))
        if (is.null(d)) next
        d$id <- L$id; d$type <- L$type; d$show <- !isFALSE(vis[[L$id]])
        d$set <- switch(L$type,
          sym = list(where = as.character(lay_now(L, "where", "tip")), align = isTRUE(lay_now(L, "align", FALSE)),
                     off = as.numeric(lay_now(L, "off", 0)), size = as.numeric(lay_now(L, "size", 1.2)),
                     gap = as.numeric(lay_now(L, "gap", 3))),
          blk = list(mode = as.character(lay_now(L, "mode", "tip_label")), grad = isTRUE(lay_now(L, "grad", FALSE)),
                     rev = isTRUE(lay_now(L, "rev", FALSE)), alpha = as.numeric(lay_now(L, "alpha", 0.35))),
          list())
        d$leg <- list(show = !isFALSE(lay_now(L, "lg_show", TRUE)), tt = !isFALSE(lay_now(L, "lg_tt", TRUE)),
                      title = as.character(lay_now(L, "lg_title", "")), what = as.character(lay_now(L, "lg_what", "both")),
                      bold = isTRUE(lay_now(L, "lg_bold", FALSE)), ital = isTRUE(lay_now(L, "lg_ital", FALSE)),
                      fam = as.character(lay_now(L, "lg_fam", "sans")), cex = as.numeric(lay_now(L, "lg_cex", 1)))
        out[[length(out) + 1L]] <- d
      }
      out
    })
    layers <- shiny::reactive(Filter(function(d) isTRUE(d$show), layers_all()))   # the layers that are drawn

    # ---- names -----------------------------------------------------------------------------------------------------
    auto_labels <- shiny::reactive({
      tr <- drawn(); lab <- tr$tip.label; al <- aligned()
      lc <- input$lab_col; la <- input$lab_app
      if (!is.null(al) && !is.null(lc) && lc %in% names(al)) {
        v <- al[[lc]]; ok <- !is.na(v) & nzchar(v); lab[ok] <- v[ok]
      }
      if (!is.null(al) && !is.null(la) && la %in% names(al)) {
        v <- al[[la]]; ok <- !is.na(v) & nzchar(v); lab[ok] <- paste0(lab[ok], " | ", v[ok])
      }
      if (isTRUE(input$underscore)) lab <- gsub("_", " ", lab)
      lab
    })
    labels <- shiny::reactive({
      tr <- drawn(); lab <- auto_labels(); ed <- rv$edits
      if (length(ed)) {
        h <- match(tr$tip.label, names(ed)); ok <- !is.na(h)
        lab[ok] <- unname(ed[h[ok]])
      }
      lab
    })
    output$names_tbl <- DT::renderDataTable({
      rv$tbl_epoch
      tr <- drawn(); lab <- auto_labels(); ed <- shiny::isolate(rv$edits)
      if (length(ed)) {
        h <- match(tr$tip.label, names(ed)); ok <- !is.na(h)
        lab[ok] <- unname(ed[h[ok]])
      }
      DT::datatable(data.frame(ID = tr$tip.label, `Name shown in the tree` = lab, check.names = FALSE, stringsAsFactors = FALSE),
                    rownames = FALSE, selection = "none", editable = list(target = "cell", disable = list(columns = 0)),
                    options = list(pageLength = 15, dom = "tp"))
    })
    shiny::observeEvent(input$names_tbl_cell_edit, {
      info <- input$names_tbl_cell_edit
      if (info$col == 0) return()
      tr <- shiny::isolate(drawn())
      id <- tr$tip.label[info$row]
      v <- trimws(as.character(info$value))
      ed <- rv$edits
      if (!nzchar(v)) ed <- ed[names(ed) != id] else ed[id] <- v
      rv$edits <- ed
    })
    shiny::observeEvent(input$reset_edits, { rv$edits <- character(); rv$tbl_epoch <- rv$tbl_epoch + 1L })

    # ---- size of the picture and place of the tree in it ---------------------------------------------------------
    # ONE size for the screen and for the downloaded files. The tree sits in a box (left, top, tw, th, in pixels);
    # the legends and the symbols after the names are outside the box.
    legend_dims <- function(L) {
      ch <- 0; rows <- 0; k <- 0
      if (!isTRUE(input$show_leg)) return(list(chars = 0, rows = 0, k = 0))
      for (z in L) {
        lg <- z$leg
        if (!is.null(lg) && !isTRUE(lg$show)) next
        lvv <- if (identical(z$type, "sym")) z$lv[z$pch != 99L] else z$lv
        if (!length(lvv)) next
        what <- if (is.null(lg)) "both" else lg$what
        m <- if (is.null(lg)) 1 else lg$cex
        f <- (if (!is.null(lg) && isTRUE(lg$bold)) 1.12 else 1) * (if (!is.null(lg) && identical(lg$fam, "mono")) 1.25 else 1)
        ttl <- if (!is.null(lg) && nzchar(lg$title)) lg$title else z$col
        ch <- max(ch, (if (identical(what, "mark")) 0 else max(nchar(lvv))) * m * f, nchar(ttl) * m * f)
        rows <- rows + (length(lvv) + 1) * m; k <- k + 1
      }
      if (identical(input$sup_style, "circle") && isTRUE(input$sup_leg) && !is.null(.vz_support(drawn()))) {
        ch <- max(ch, 14); rows <- rows + 4; k <- k + 1
      }
      list(chars = ch, rows = rows, k = k)
    }
    size_info <- shiny::reactive({
      tr <- drawn(); lab <- labels(); L <- layers()
      ly <- .vz_layouts[[.vz_or(input$layout, "phylogram")]]
      rect <- isTRUE(ly$rect)
      n <- ape::Ntip(tr)
      cex <- .vz_or(input$cex, 0.9)
      zx <- .vz_or(input$zoom_x, 100) / 100; zy <- .vz_or(input$zoom_y, 100) / 100
      line_px <- 19 * cex                                         # height of one row of names
      wf <- (if (isTRUE(input$bold)) 1.1 else 1) * (if (identical(input$fam, "mono")) 1.2 else 1)   # bold and monospace names are wider
      lab_px <- max(nchar(lab), 1L) * 7.6 * cex * wf + 12         # width of the longest name
      body <- 520 * zx                                            # width of the branches of a rectangular tree
      sym_px <- 0; inline_px <- 0                                 # room taken by the symbol layers: after the names / before the names
      for (z in L) if (identical(z$type, "sym") && any(z$pch != 99L)) {
        s <- z$set
        if (identical(s$where, "end")) sym_px <- sym_px + 14 * s$size + 10 + max(0, s$off) / 100 * body
        else inline_px <- inline_px + 14 * s$size + (s$gap + max(0, s$off)) / 100 * body
      }
      if (rect) {
        tw <- body + lab_px + inline_px
        th <- n * line_px * zy + 20
        out_r <- sym_px                                           # symbols after the names stick out at the right
      } else {
        fan <- identical(ly$type, "fan")
        op <- if (fan) .vz_or(input$open, 0) else 0
        R0 <- if (fan) max(180, n * line_px / (2 * pi * max(0.2, 1 - op / 360))) else max(220, 0.45 * n * line_px)
        tw <- 2 * (R0 * zx + lab_px + inline_px + sym_px) + 20
        th <- 2 * (R0 * zy + lab_px + inline_px + sym_px) + 20
        out_r <- 0
      }
      ld <- legend_dims(L)
      leg_cex <- .vz_or(input$leg_cex, 0.9)
      pos <- .vz_or(input$leg_pos, "right")
      leg_px <- if (ld$chars > 0) ld$chars * 7.4 * leg_cex + 45 else 0
      leg_r <- if (identical(pos, "right")) leg_px else 0
      leg_l <- if (identical(pos, "left")) leg_px else 0
      leg_h <- if (pos %in% c("right", "left")) ld$rows * 17 * leg_cex + ld$k * 14 else 0
      has_title <- nzchar(trimws(.vz_or(input$title, "")))
      padL <- 15 + leg_l; padT <- 12 + (if (has_title) 30 else 0); padR <- 15; padB <- 20
      left <- max(5, padL + rv$tx); top <- max(5, padT + rv$ty)   # where the tree box starts (dragging moves it)
      W_auto <- min(8000, max(350, left + tw + out_r + leg_r + padR))
      H_auto <- min(9000, max(300, top + th + padB, top + leg_h + padB))
      memo$W <- W_auto; memo$H <- H_auto
      auto <- !identical(input$auto_size, FALSE)
      if (auto) {
        W <- W_auto; H <- H_auto
      } else {
        W <- max(300, .vz_or(input$man_w, 900)); H <- max(300, .vz_or(input$man_h, 700))
        avail_w <- max(100, W - left - out_r - leg_r - padR); avail_h <- max(100, H - top - padB)
        tw <- min(avail_w, avail_w * zx); th <- min(avail_h, avail_h * zy)
        left <- min(left, max(5, W - tw - 5)); top <- min(top, max(5, H - th - 5))
      }
      plt <- pmin(1, pmax(0, c(left / W, (left + tw) / W, 1 - (top + th) / H, 1 - top / H)))
      if (!all(is.finite(plt)) || plt[1] >= plt[2] || plt[3] >= plt[4]) plt <- c(0.05, 0.95, 0.05, 0.95)
      list(W = W, H = H, plt = plt, padL = padL, padT = padT, sym_out = out_r, auto = auto)
    })
    output$size_txt <- shiny::renderUI({
      si <- tryCatch(size_info(), error = function(e) NULL)
      if (is.null(si)) return(NULL)
      shiny::helpText(sprintf("Picture: %d x %d pixels (%s).", round(si$W), round(si$H),
                              if (si$auto) "follows the tree" else "set by hand"))
    })
    # when the user takes over the picture size, start from the size it has now
    shiny::observeEvent(input$auto_size, {
      if (isTRUE(pend$skip_auto)) { pend$skip_auto <- NULL; return() }     # the size comes from an opened project
      if (identical(input$auto_size, FALSE) && !is.null(memo$W)) {
        shiny::updateNumericInput(session, "man_w", value = round(memo$W))
        shiny::updateNumericInput(session, "man_h", value = round(memo$H))
      }
    }, ignoreInit = TRUE)
    shiny::observeEvent(input$tree_reset, { rv$tx <- 0; rv$ty <- 0 })

    # ---- drawing -----------------------------------------------------------------------------------------------------
    plot_args <- shiny::reactive({
      L <- layers(); si <- size_info()
      sup_style <- .vz_or(input$sup_style, "text")
      sup_col <- .vz_or(input$sup_col, "#555555"); sup_cex <- .vz_or(input$sup_cex, 1.2)
      legs <- list(); dl <- list()
      for (z in L) {
        pos <- match(z$vals, z$lv)
        if (identical(z$type, "sym") && any(z$pch != 99L)) {
          dl[[length(dl) + 1L]] <- list(type = "sym", pch = z$pch[pos], col = z$cols[pos], size = z$set$size, gap = z$set$gap / 100,
                                        where = z$set$where, align = z$set$align, off = z$set$off)
        } else if (identical(z$type, "blk")) {
          dl[[length(dl) + 1L]] <- list(type = "blk", col = z$cols[pos], grp = z$vals, mode = z$set$mode, grad = z$set$grad,
                                        rev = z$set$rev, alpha = z$set$alpha)
        }
      }
      if (isTRUE(input$show_leg)) {
        for (z in L) {
          lg <- z$leg
          if (!isTRUE(lg$show)) next
          ttl <- if (!isTRUE(lg$tt)) NULL else if (nzchar(trimws(lg$title))) trimws(lg$title) else z$col
          sty <- list(what = lg$what, font = 1L + (if (isTRUE(lg$bold)) 1L else 0L) + (if (isTRUE(lg$ital)) 2L else 0L),
                      family = lg$fam, cex = .vz_or(input$leg_cex, 0.9) * lg$cex)
          if (identical(z$type, "sym")) {
            keep <- z$pch != 99L
            if (any(keep)) legs[[length(legs) + 1L]] <- c(list(id = z$id, title = ttl, labels = z$lv[keep], pch = z$pch[keep],
                                                               col = z$cols[keep], bg = z$cols[keep], ptcex = 1.3), sty)
          } else if (identical(z$type, "blk")) {
            legs[[length(legs) + 1L]] <- c(list(id = z$id, title = ttl, labels = z$lv, pch = rep(22L, length(z$lv)),
                                                col = z$cols, bg = grDevices::adjustcolor(z$cols, alpha.f = z$set$alpha), ptcex = 1.8), sty)
          } else {
            legs[[length(legs) + 1L]] <- c(list(id = z$id, title = ttl, labels = z$lv, pch = rep(15L, length(z$lv)),
                                                col = z$cols, bg = z$cols, ptcex = 1.5), sty)
          }
        }
        if (identical(sup_style, "circle") && isTRUE(input$sup_leg) && !is.null(.vz_support(drawn()))) {
          v <- c(100, 80, 60)
          legs[[length(legs) + 1L]] <- list(id = "sup", title = "Branch support", labels = as.character(v), pch = rep(21L, 3),
                                            col = rep("#222222", 3), bg = rep(sup_col, 3), ptcex = sup_cex * (0.35 + 0.65 * v / 100))
        }
      }
      lcl <- Filter(function(z) identical(z$type, "lc"), L)       # only one colour per name: the last "Name colour" layer wins
      lc <- if (length(lcl)) lcl[[length(lcl)]] else NULL
      list(layout = .vz_or(input$layout, "phylogram"),
           cex = .vz_or(input$cex, 0.9), lwd = .vz_or(input$lwd, 1.5), edge_col = "#222222",
           italic = isTRUE(input$italic), bold = isTRUE(input$bold), family = .vz_or(input$fam, "sans"),
           tipcol = if (is.null(lc)) NULL else lc$cols[match(lc$vals, lc$lv)],
           layers = dl, sym_out = si$sym_out,
           align = isTRUE(input$align), guide = .vz_or(input$guide, "0"),
           sup_style = sup_style, min_sup = .vz_or(input$min_sup, 0), sup_col = sup_col, sup_cex = sup_cex,
           scale_bar = isTRUE(input$scale_bar), title = trimws(.vz_or(input$title, "")),
           legends = legs, leg_cex = .vz_or(input$leg_cex, 0.9), leg_pos = .vz_or(input$leg_pos, "right"),
           leg_xy = rv$leg_xy, rotate = .vz_or(input$rotate, 0), open = .vz_or(input$open, 0),
           plt = si$plt, W = si$W, H = si$H)
    })
    output$tree_plot <- shiny::renderPlot({
      tr <- drawn(); a <- plot_args(); lab <- labels()
      a$sel <- rv$sel
      res <- .vz_draw(tr, lab, a)
      for (nm in c("nx", "ny", "rect", "tip_end", "half_row")) assign(nm, res[[nm]], envir = coords)
      coords$W <- a$W; coords$H <- a$H
      # tell the browser where the legends and the tree are, so it can drag them without asking the server
      session$sendCustomMessage(paste0(ns("tree_plot"), "_layout"),
                                list(legends = unname(res$leg_boxes), tree = res$tree_box))
    }, width = function() tryCatch(as.integer(size_info()$W), error = function(e) 400L),
       height = function() tryCatch(as.integer(size_info()$H), error = function(e) 60L), res = 96)
    output$msg <- shiny::renderUI({
      if (!is.null(rv$msg)) return(shiny::div(class = paste0("alert alert-", rv$msg$type), rv$msg$txt))
      tr <- tryCatch(src_tree(), error = function(e) NULL)
      if (is.null(tr)) shiny::helpText("No tree yet: build one in the Phylogeny tab (NJ or IQ-TREE), or choose 'A tree file' on the left.")
    })

    # ---- dragging (left button): a legend or the tree ------------------------------------------------------------------
    shiny::observeEvent(input$leg_reset, { rv$leg_xy <- list() })
    shiny::observeEvent(input$leg_pos, {
      if (isTRUE(pend$skip_leg)) { pend$skip_leg <- NULL; return() }       # the legend places come from an opened project
      rv$leg_xy <- list()
    }, ignoreInit = TRUE)
    shiny::observeEvent(input$drag, {
      d <- input$drag
      if (identical(d$kind, "legend") && !is.null(d$id)) {
        xy <- rv$leg_xy; xy[[d$id]] <- c(d$x, d$y); rv$leg_xy <- xy
      } else if (identical(d$kind, "tree")) {
        si <- size_info()
        rv$tx <- max(rv$tx + d$dx, 5 - si$padL)                   # the tree cannot leave the picture at the left / top
        rv$ty <- max(rv$ty + d$dy, 5 - si$padT)
      }
    })

    # ---- right button: select a tip or a node, reroot -----------------------------------------------------------------
    describe_node <- function(tr, k, lab) {
      n <- ape::Ntip(tr)
      if (k <= n) return(paste0("Tip: ", lab[k]))
      if (k == n + 1L) return("The root of the tree")
      tips <- ape::extract.clade(tr, k)$tip.label
      nm <- lab[match(tips, tr$tip.label)]
      sprintf("Node with %d tips (%s%s)", length(tips), paste(utils::head(nm, 3), collapse = ", "),
              if (length(tips) > 3) ", ..." else "")
    }
    do_reroot <- function() {
      tr <- drawn(); k <- rv$sel
      shiny::req(tr, k)
      n <- ape::Ntip(tr)
      if (k == n + 1L) return()
      og <- if (k <= n) tr$tip.label[k] else ape::extract.clade(tr, k)$tip.label
      shiny::updateRadioButtons(session, "root", selected = "outgroup")
      shiny::updateSelectizeInput(session, "outgroup", selected = og)
      rv$sel <- NULL
    }
    shiny::observeEvent(input$ctx, {
      cl <- input$ctx
      tr <- tryCatch(drawn(), error = function(e) NULL)
      if (is.null(tr) || is.null(coords$nx) || is.null(cl$fx)) return()
      n <- ape::Ntip(tr); tips <- seq_len(n)
      dxp <- (coords$nx - cl$fx) * coords$W                       # distance from the click to every node, in pixels
      dyp <- (coords$ny - cl$fy) * coords$H
      k <- NULL
      if (isTRUE(coords$rect) && !is.null(coords$tip_end)) {     # a click on a name selects that tip
        cand <- which(abs(dyp[tips]) <= coords$half_row * coords$H &
                        (cl$fx - coords$nx[tips]) * coords$W >= -8 & cl$fx <= coords$tip_end[tips])
        if (length(cand)) k <- cand[which.min(abs(dyp[cand]))]
      }
      if (is.null(k)) k <- which.min(dxp^2 + dyp^2)
      rv$sel <- k
      items <- list()
      if (k != n + 1L)
        items[[1L]] <- list(label = if (k <= n) "Reroot on this tip (it becomes the outgroup)"
                                    else "Reroot on this node (its tips become the outgroup)", act = "reroot")
      items[[length(items) + 1L]] <- list(label = "Clear the selection", act = "clear")
      session$sendCustomMessage(paste0(ns("tree_plot"), "_menu"),
                                list(x = cl$cx, y = cl$cy, title = describe_node(tr, k, labels()), items = items))
    })
    shiny::observeEvent(input$menu_act, {
      m <- input$menu_act
      if (identical(m$act, "reroot")) do_reroot()
      else if (identical(m$act, "clear")) rv$sel <- NULL
      else if (identical(m$act, "leg_reset_one") && !is.null(m$id)) { xy <- rv$leg_xy; xy[[m$id]] <- NULL; rv$leg_xy <- xy }
    })
    output$sel_ui <- shiny::renderUI({     # the selection is also shown here, in case the menu was closed
      tr <- tryCatch(drawn(), error = function(e) NULL)
      k <- rv$sel
      if (is.null(tr) || is.null(k)) return(NULL)
      n <- ape::Ntip(tr)
      shiny::div(class = "alert alert-info", style = "margin-top:8px;padding:8px 12px",
        shiny::strong(describe_node(tr, k, shiny::isolate(labels()))), shiny::tags$br(),
        if (k != n + 1L) shiny::actionButton(ns("reroot_sel"), "Reroot on the selection", icon = shiny::icon("code-branch")),
        shiny::actionButton(ns("clear_sel"), "Clear the selection"))
    })
    shiny::observeEvent(input$clear_sel, { rv$sel <- NULL })
    shiny::observeEvent(input$reroot_sel, { do_reroot() })

    # ---- downloads (the picture has the same size and layout as on the screen) --------------------------------------------
    save_fig <- function(type) function(file) {
      tr <- drawn(); a <- plot_args(); lab <- labels(); si <- size_info()
      shiny::req(tr)
      a$sel <- NULL
      w <- si$W / 96; h <- si$H / 96
      dpi <- min(as.integer(.vz_or(input$dl_dpi, "300")), floor(20000 / max(w, h)))
      switch(type,
        svg = grDevices::svg(file, width = w, height = h),
        png = grDevices::png(file, width = round(w * dpi), height = round(h * dpi), res = dpi),
        pdf = grDevices::pdf(file, width = w, height = h))
      on.exit(grDevices::dev.off(), add = TRUE)
      .vz_draw(tr, lab, a)
    }
    output$dl_svg <- shiny::downloadHandler(function() "tree_figure.svg", save_fig("svg"))
    output$dl_png <- shiny::downloadHandler(function() "tree_figure.png", save_fig("png"))
    output$dl_pdf <- shiny::downloadHandler(function() "tree_figure.pdf", save_fig("pdf"))
    output$dl_nwk <- shiny::downloadHandler(function() "tree_edited.nwk",
      function(file) { tr <- drawn(); shiny::req(tr); ape::write.tree(tr, file) })
    output$dl_nwk_lab <- shiny::downloadHandler(function() "tree_edited_names.nwk",
      function(file) { tr <- drawn(); shiny::req(tr); tr$tip.label <- labels(); ape::write.tree(tr, file) })

    # ---- external viewers (Newick as computed) ---------------------------------------------------------------------------
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

    # ---- project file: save everything needed to go on editing, and open it again ---------------------------------------
    proj_ids <- list(
      radio  = c("src", "root", "order"),
      select = c("layout", "pal", "guide", "fam", "sup_style", "leg_pos", "dl_dpi"),
      check  = c("last_only", "underscore", "italic", "bold", "align", "scale_bar", "auto_size",
                 "sup_leg", "show_leg"),
      slider = c("rotate", "open", "cex", "lwd", "zoom_x", "zoom_y", "min_sup", "sup_cex", "leg_cex"),
      number = c("man_w", "man_h"),
      text   = "title")
    output$dl_proj <- shiny::downloadHandler(
      function() paste0("sangerly_visualisation_", format(Sys.time(), "%Y%m%d"), ".json"),
      function(file) {
        tr <- tryCatch(src_tree(), error = function(e) NULL)
        vals <- list()
        for (k in names(proj_ids)) for (i in proj_ids[[k]]) { v <- input[[i]]; if (!is.null(v)) vals[[i]] <- v }
        m <- rv$meta
        tab <- if (is.null(m)) NULL else {
          m[is.na(m)] <- ""
          list(names = names(m), columns = unname(lapply(m, as.character)))
        }
        cols <- list(); for (i in c("id_col", "lab_col", "lab_app")) if (!is.null(input[[i]])) cols[[i]] <- input[[i]]
        ly <- tryCatch(layers_all(), error = function(e) list())
        annot <- lapply(ly, function(d) list(id = d$id, type = d$type, col = d$col, show = d$show, set = d$set, legend = d$leg,
                                             colours = list(col = d$col, pal = .vz_or(input$pal, "distinct"), lv = d$lv, cols = d$cols, pch = d$pch)))
        ed <- rv$edits
        out <- list(app = "SangeRly", type = "visualisation", format = 2L, saved = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                    tree = if (is.null(tr)) NULL else list(newick = paste(ape::write.tree(tr), collapse = ""), from = .vz_or(input$src, "nj")),
                    table = tab, inputs = vals, columns = cols, outgroup = as.list(as.character(input$outgroup)),
                    sup_col = input$sup_col, annot = unname(annot),
                    edits = if (length(ed)) as.list(ed) else list(),
                    tree_xy = c(rv$tx, rv$ty), legends = rv$leg_xy)
        writeLines(jsonlite::toJSON(out, auto_unbox = TRUE, null = "null", na = "null", pretty = TRUE, digits = NA), file, useBytes = TRUE)
      })
    shiny::observeEvent(input$proj_file, {
      f <- input$proj_file
      p <- tryCatch(jsonlite::fromJSON(f$datapath, simplifyVector = FALSE), error = function(e) NULL)
      if (is.null(p) || !identical(p$app, "SangeRly") || !identical(p$type, "visualisation")) {
        shiny::showNotification("This is not a SangeRly Visualisation project file.", type = "error"); return()
      }
      tr <- NULL
      if (!is.null(p$tree$newick)) {
        tr <- tryCatch(ape::read.tree(text = p$tree$newick), error = function(e) NULL)
        if (inherits(tr, "multiPhylo")) tr <- tr[[1]]
      }
      if (is.null(tr)) { shiny::showNotification("The project holds no tree that can be read.", type = "error"); return() }
      chr <- function(x) as.character(unlist(x))
      rv$file_trees <- list(tr)                                     # the tree travels inside the project
      if (!is.null(p$table)) {
        cl <- lapply(p$table$columns, chr); names(cl) <- chr(p$table$names)
        rv$meta <- as.data.frame(cl, stringsAsFactors = FALSE, check.names = FALSE, optional = TRUE); names(rv$meta) <- chr(p$table$names)
        rv$sheets <- NULL
      } else { rv$meta <- NULL; rv$sheets <- NULL }
      pend$cols <- lapply(p$columns, function(x) chr(x)[1])
      pend$outgroup <- if (length(p$outgroup)) chr(p$outgroup) else NULL
      rv$epoch <- rv$epoch + 1L                                     # new keys: no colour of the old picture is kept
      # the annotation layers (a project of an older version holds one layer of each kind: symbols, blocks, name colour)
      ivp <- p$inputs
      num1 <- function(v, d) { v <- suppressWarnings(as.numeric(chr(v)[1])); if (is.finite(v)) v else d }
      specs <- list()
      if (length(p$annot)) {
        for (a in p$annot) {
          ty <- chr(a$type)[1]
          if (is.na(ty) || !(ty %in% .vz_layer_types)) next
          specs[[length(specs) + 1L]] <- list(old = chr(a$id)[1], type = ty, col = chr(a$col)[1], show = !identical(a$show, FALSE),
                                              set = a$set, ov = a$colours, leg = a$legend)
        }
      } else {
        ci <- function(k) { v <- p$columns[[k]]; if (is.null(v)) "(none)" else chr(v)[1] }
        specs <- list(
          list(old = "sy", type = "sym", col = ci("sym_col"), show = TRUE, ov = p$layers[["sy"]],
               set = list(where = ivp$sym_where, align = ivp$sym_align, off = ivp$sym_off, size = ivp$sym_size, gap = ivp$gap)),
          list(old = "bk", type = "blk", col = ci("blk_col"), show = TRUE, ov = p$layers[["bk"]],
               set = list(mode = ivp$blk_mode, grad = ivp$blk_grad, rev = ivp$blk_rev, alpha = ivp$blk_alpha)),
          list(old = "lc", type = "lc", col = ci("lc_col"), show = TRUE, ov = p$layers[["lc"]], set = list()))
      }
      idmap <- list(); newlay <- list(); newvis <- list(); pl <- list(); ovl <- list()
      for (sp in specs) {
        rv$lay_n <- rv$lay_n + 1L
        nid <- paste0("ly", rv$lay_n)                              # new ids: nothing of the old controls is kept
        if (length(sp$old) == 1L && !is.na(sp$old)) idmap[[sp$old]] <- nid
        newlay[[length(newlay) + 1L]] <- list(id = nid, type = sp$type)
        newvis[[nid]] <- isTRUE(sp$show)
        st <- sp$set
        s1 <- list(col = sp$col)
        if (identical(sp$type, "sym"))
          s1 <- c(s1, list(where = if (identical(chr(st$where)[1], "end")) "end" else "tip", align = isTRUE(st$align),
                           off = num1(st$off, 0), size = num1(st$size, 1.2), gap = num1(st$gap, 3)))
        if (identical(sp$type, "blk"))
          s1 <- c(s1, list(mode = if (length(st$mode)) chr(st$mode)[1] else "tip_label", grad = isTRUE(st$grad), rev = isTRUE(st$rev),
                           alpha = num1(st$alpha, 0.35)))
        lgs <- sp$leg
        if (!is.null(lgs)) s1 <- c(s1, list(lg_show = !identical(lgs$show, FALSE), lg_tt = !identical(lgs$tt, FALSE),
                                            lg_title = if (length(lgs$title)) chr(lgs$title)[1] else "",
                                            lg_what = if (length(lgs$what)) chr(lgs$what)[1] else "both",
                                            lg_bold = isTRUE(lgs$bold), lg_ital = isTRUE(lgs$ital),
                                            lg_fam = if (length(lgs$fam)) chr(lgs$fam)[1] else "sans", lg_cex = num1(lgs$cex, 1)))
        pl[[nid]] <- s1
        ov <- sp$ov
        if (!is.null(ov)) ovl[[nid]] <- list(col = chr(ov$col)[1], pal = chr(ov$pal)[1], lv = chr(ov$lv), cols = chr(ov$cols),
                                             pch = as.integer(chr(ov$pch)), epoch = rv$epoch)
      }
      rv$layer_ov <- ovl; rv$vis <- newvis; pend$lay <- pl
      memo$open <- if (length(newlay)) newlay[[1]]$id else character()
      rv$lay <- newlay
      ed <- unlist(p$edits); rv$edits <- if (length(ed)) ed else character(); rv$tbl_epoch <- rv$tbl_epoch + 1L
      xy <- as.numeric(unlist(p$tree_xy)); if (length(xy) == 2L && all(is.finite(xy))) { rv$tx <- xy[1]; rv$ty <- xy[2] }
      lg <- lapply(p$legends, function(v) as.numeric(unlist(v)))
      newxy <- list()
      for (nm in names(lg)) newxy[[if (!is.null(idmap[[nm]])) idmap[[nm]] else nm]] <- lg[[nm]]   # legends follow their layers
      rv$leg_xy <- newxy
      rv$sel <- NULL; rv$msg <- NULL
      iv <- p$inputs; iv$src <- "file"
      if (!is.null(iv$leg_pos) && !identical(iv$leg_pos, input$leg_pos)) pend$skip_leg <- TRUE
      if (identical(iv$auto_size, FALSE) && !identical(input$auto_size, FALSE)) pend$skip_auto <- TRUE
      num <- function(v) suppressWarnings(as.numeric(v))
      for (i in proj_ids$radio)  if (!is.null(iv[[i]])) shiny::updateRadioButtons(session, i, selected = chr(iv[[i]])[1])
      for (i in proj_ids$select) if (!is.null(iv[[i]])) shiny::updateSelectInput(session, i, selected = chr(iv[[i]])[1])
      for (i in proj_ids$check)  if (!is.null(iv[[i]])) shiny::updateCheckboxInput(session, i, value = isTRUE(iv[[i]]))
      for (i in proj_ids$slider) if (!is.null(iv[[i]]) && is.finite(num(iv[[i]]))) shiny::updateSliderInput(session, i, value = num(iv[[i]]))
      for (i in proj_ids$number) if (!is.null(iv[[i]]) && is.finite(num(iv[[i]]))) shiny::updateNumericInput(session, i, value = num(iv[[i]]))
      for (i in proj_ids$text)   if (!is.null(iv[[i]])) shiny::updateTextInput(session, i, value = chr(iv[[i]])[1])
      if (!is.null(p$sup_col)) session$sendCustomMessage("sg_setcolor", list(id = ns("sup_col"), value = chr(p$sup_col)[1]))
      shiny::showNotification(sprintf("Project opened (%d tips%s). The tree is stored in the project file.", ape::Ntip(tr),
                                      if (!is.null(p$saved)) paste0(", saved ", chr(p$saved)[1]) else ""), type = "message")
    })
  })
}
