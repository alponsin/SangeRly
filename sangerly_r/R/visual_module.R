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
  has_sym <- !is.null(a$sym) && any(a$sym$pch != 99L)
  sym_end <- has_sym && identical(a$sym_where, "end")       # symbols after the names (a column at the right)
  gapH <- a$gap * H
  offH <- a$sym_off / 100 * H                               # symbol offset: negative = left / inwards, positive = right / outwards
  # distance between a tip and the start of its name
  off <- if (has_sym && !sym_end) gapH + max(0, offH)
         else if (has_sym || !is.null(a$blk)) max(0.01 * H, 0.5 * gapH)
         else 0.01 * H
  fnt <- if (isTRUE(a$italic)) 3 else 1
  tr2 <- tr
  tr2$tip.label <- lab
  tipcol <- if (is.null(a$tipcol)) "black" else a$tipcol
  aligned <- rect && (isTRUE(a$align) || (has_sym && !sym_end && isTRUE(a$sym_align)))

  # the tree lives in a box inside the picture (the legends and the symbols after the names are outside the box)
  graphics::par(plt = a$plt, xpd = FALSE)

  args <- list(tr2, type = ly$type, use.edge.length = use_len, show.tip.label = TRUE,
               tip.color = tipcol, font = fnt, cex = a$cex, label.offset = off,
               edge.width = a$lwd, edge.color = a$edge_col)
  if (aligned) args$align.tip.label <- as.integer(.vz_or(a$guide, 0L))   # the number is the line type; 0 = no line
  if (is_fan) { args$rotate.tree <- a$rotate; args$open.angle <- a$open }
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
  if (!is.null(a$blk) && (rect || is_fan)) {
    mode <- .vz_or(a$blk_mode, "tip_label")
    ord <- tryCatch({ e <- ape::reorder.phylo(tr, "cladewise")$edge[, 2]; e[e <= n] }, error = function(e) tips)
    rl <- rle(as.character(a$blk$grp[ord]))                    # runs of tips of the same group, in the order of the picture
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
    nsl <- if (isTRUE(a$blk_grad)) (if (rect) 30L else 14L) else 1L
    for (s in seq_len(nsl)) {
      f0 <- (s - 1) / nsl; f1 <- s / nsl; u <- (f0 + f1) / 2
      al <- if (nsl == 1L) a$blk_alpha else a$blk_alpha * (if (isTRUE(a$blk_rev)) 0.1 + 0.9 * u else 1 - 0.9 * u)
      cs <- grDevices::adjustcolor(a$blk$col, alpha.f = al)
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
    keep <- which(a$sym$pch != 99L)
    if (rect || is_fan) {
      ref <- (if (sym_end) lab_end else ptip) + 0.5 * gapH
      if (isTRUE(a$sym_align)) ref <- rep(max(ref), n)         # all symbols in one column
      pp <- ref + offH
      if (rect) { sx <- pp; sy <- yy[tips] } else { sx <- cx + pp * cos(th); sy <- cy + pp * sin(th) }
    } else {
      ux <- xx[tips] - mean(xx); uy <- yy[tips] - mean(yy)
      d <- sqrt(ux^2 + uy^2); d[d == 0] <- 1
      pp <- (if (sym_end) off + lw else 0) + 0.5 * gapH + offH
      sx <- xx[tips] + pp * ux / d; sy <- yy[tips] + pp * uy / d
    }
    graphics::points(sx[keep], sy[keep], pch = a$sym$pch[keep], col = a$sym$col[keep], bg = a$sym$col[keep],
                     cex = a$sym_size)
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
      r <- graphics::legend(graphics::grconvertX(xn, "ndc", "user"), graphics::grconvertY(yn, "ndc", "user"),
                            legend = L$labels, title = L$title, pch = L$pch, col = L$col, pt.bg = L$bg,
                            pt.cex = L$ptcex, bty = "n", cex = a$leg_cex, xjust = xj, yjust = yj, title.adj = 0,
                            y.intersp = 1.1)
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
    shiny::fluidRow(
      shiny::column(4,
        .card("1. Tree",
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
        .card("2. Groups from an Excel / CSV table",
          shiny::fileInput(ns("meta"), NULL, buttonLabel = "Open table...", placeholder = "",
                           accept = c(".xlsx", ".xls", ".csv", ".tsv", ".txt")),
          shiny::uiOutput(ns("meta_ui")),
          shiny::selectInput(ns("pal"), "Default colours", .vz_pal_choices),
          .vz_section("Symbols: shape and colour per group", shiny::uiOutput(ns("sym_levels")), open = TRUE),
          .vz_section("Colour blocks: colour per group", shiny::uiOutput(ns("blk_levels"))),
          .vz_section("Name colours: colour per group", shiny::uiOutput(ns("lc_levels"))),
          shiny::actionButton(ns("reset_cols"), "Reset group colours and shapes", icon = shiny::icon("rotate-left"))),
        .card("3. Appearance",
          .vz_section("Names and branches",
            shiny::checkboxInput(ns("underscore"), "Replace _ with a space in the names", TRUE),
            shiny::checkboxInput(ns("italic"), "Italic names", TRUE),
            shiny::sliderInput(ns("cex"), "Name size", 0.4, 2, 0.9, 0.1),
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
          .vz_section("Symbols",
            shiny::radioButtons(ns("sym_where"), "Position",
              c("Before the name (at the tip)" = "tip", "After the name (column at the right)" = "end")),
            shiny::checkboxInput(ns("sym_align"), "Align the symbols in one column", FALSE),
            shiny::numericInput(ns("sym_off"), "Offset (- to the left, + to the right; % of the tree height)",
                                0, min = -100, max = 300, step = 1),
            shiny::sliderInput(ns("sym_size"), "Symbol size", 0.4, 3, 1.2, 0.1),
            shiny::sliderInput(ns("gap"), "Space around the symbols (% of the tree height)", 0, 10, 3, 0.5)),
          .vz_section("Colour blocks",
            shiny::selectInput(ns("blk_mode"), "Width of the blocks", .vz_blk_modes),
            shiny::checkboxInput(ns("blk_grad"), "Colour gradient (fades out)", FALSE),
            shiny::conditionalPanel(sprintf("input['%s']", ns("blk_grad")),
              shiny::checkboxInput(ns("blk_rev"), "Reverse the gradient", FALSE)),
            shiny::sliderInput(ns("blk_alpha"), "Colour strength", 0.1, 1, 0.35, 0.05),
            shiny::helpText("Blocks are drawn in the rectangular and circular layouts. The 'LCA' widths start at the common ",
                            "ancestor of the neighbouring tips that belong to the same group.")),
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
            shiny::actionButton(ns("leg_reset"), "Reset the legend positions", icon = shiny::icon("rotate-left"))))),
      shiny::column(8,
        shiny::uiOutput(ns("msg")),
        shiny::tabsetPanel(id = ns("tabs"),
          shiny::tabPanel("Tree", shiny::br(),
            shiny::div(style = "overflow:auto;max-height:84vh;border:1px solid #e3e8ee;border-radius:6px;background:#fff",
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
                              "became '_'). Please read the terms of use of each website before uploading data.")))))))
}

# ---------------------------------------------------------------------------------------------------
# Server
# ---------------------------------------------------------------------------------------------------

visualServer <- function(id, nj = shiny::reactive(NULL), iq = shiny::reactive(NULL)) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns
    rv <- shiny::reactiveValues(file_trees = NULL, meta = NULL, sheets = NULL, edits = character(),
                                epoch = 0L, tbl_epoch = 0L, sel = NULL, msg = NULL,
                                leg_xy = list(), tx = 0, ty = 0)
    coords <- new.env()                 # nodes of the last drawing (for the right-click)
    memo <- new.env()                   # size of the picture when it follows the tree
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
      shiny::updateSelectizeInput(session, "outgroup", choices = if (is.null(tr)) character() else tr$tip.label,
                                  selected = shiny::isolate(input$outgroup), server = FALSE)
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
      shiny::tagList(
        if (length(rv$sheets) > 1) shiny::selectInput(ns("sheet"), "Sheet", rv$sheets, selected = shiny::isolate(input$sheet)),
        shiny::selectInput(ns("id_col"), "Column with the sequence names", cols, selected = guess),
        shiny::uiOutput(ns("match_info")),
        shiny::selectInput(ns("lab_col"), "Names shown in the tree: from column", c("(name in the tree)", cols)),
        shiny::selectInput(ns("lab_app"), "Add after the name", ch),
        shiny::selectInput(ns("sym_col"), "Symbols by", ch),
        shiny::selectInput(ns("blk_col"), "Colour blocks by", ch),
        shiny::selectInput(ns("lc_col"), "Name colour by", ch))
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
    output$sym_levels <- shiny::renderUI(layer_rows("sy", input$sym_col, TRUE))
    output$blk_levels <- shiny::renderUI(layer_rows("bk", input$blk_col, FALSE))
    output$lc_levels  <- shiny::renderUI(layer_rows("lc", input$lc_col, FALSE))
    shiny::observeEvent(input$reset_cols, { rv$epoch <- rv$epoch + 1L })
    layers <- shiny::reactive(list(sy = layer_map("sy", input$sym_col),
                                   bk = layer_map("bk", input$blk_col),
                                   lc = layer_map("lc", input$lc_col)))

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
      if (!is.null(L$sy) && any(L$sy$pch != 99L)) {
        keep <- L$sy$pch != 99L
        ch <- max(ch, nchar(L$sy$lv[keep]), nchar(L$sy$col)); rows <- rows + sum(keep) + 1; k <- k + 1
      }
      for (x in list(L$bk, L$lc)) if (!is.null(x)) {
        ch <- max(ch, nchar(x$lv), nchar(x$col)); rows <- rows + length(x$lv) + 1; k <- k + 1
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
      lab_px <- max(nchar(lab), 1L) * 7.6 * cex + 12              # width of the longest name
      has_sym <- !is.null(L$sy) && any(L$sy$pch != 99L)
      sym_end <- has_sym && identical(input$sym_where, "end")
      sym_size <- .vz_or(input$sym_size, 1.2)
      off_pct <- .vz_or(input$sym_off, 0); gap <- .vz_or(input$gap, 3)
      body <- 520 * zx                                            # width of the branches of a rectangular tree
      sym_px <- if (sym_end) 14 * sym_size + 10 + max(0, off_pct) / 100 * body else 0
      inline_px <- if (has_sym && !sym_end) 14 * sym_size + (gap + max(0, off_pct)) / 100 * body else 0
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
      if (identical(input$auto_size, FALSE) && !is.null(memo$W)) {
        shiny::updateNumericInput(session, "man_w", value = round(memo$W))
        shiny::updateNumericInput(session, "man_h", value = round(memo$H))
      }
    }, ignoreInit = TRUE)
    shiny::observeEvent(input$tree_reset, { rv$tx <- 0; rv$ty <- 0 })

    # ---- drawing -----------------------------------------------------------------------------------------------------
    plot_args <- shiny::reactive({
      L <- layers(); sy <- L$sy; bk <- L$bk; lc <- L$lc; si <- size_info()
      alpha <- .vz_or(input$blk_alpha, 0.35)
      sup_style <- .vz_or(input$sup_style, "text")
      sup_col <- .vz_or(input$sup_col, "#555555"); sup_cex <- .vz_or(input$sup_cex, 1.2)
      legs <- list()
      if (isTRUE(input$show_leg)) {
        if (!is.null(sy)) {
          keep <- sy$pch != 99L
          if (any(keep)) legs[[length(legs) + 1L]] <- list(id = "sy", title = sy$col, labels = sy$lv[keep], pch = sy$pch[keep],
                                                           col = sy$cols[keep], bg = sy$cols[keep], ptcex = 1.3)
        }
        if (!is.null(bk)) legs[[length(legs) + 1L]] <- list(id = "bk", title = bk$col, labels = bk$lv, pch = rep(22L, length(bk$lv)),
                                                            col = bk$cols, bg = grDevices::adjustcolor(bk$cols, alpha.f = alpha), ptcex = 1.8)
        if (!is.null(lc)) legs[[length(legs) + 1L]] <- list(id = "lc", title = lc$col, labels = lc$lv, pch = rep(15L, length(lc$lv)),
                                                            col = lc$cols, bg = lc$cols, ptcex = 1.5)
        if (identical(sup_style, "circle") && isTRUE(input$sup_leg) && !is.null(.vz_support(drawn()))) {
          v <- c(100, 80, 60)
          legs[[length(legs) + 1L]] <- list(id = "sup", title = "Branch support", labels = as.character(v), pch = rep(21L, 3),
                                            col = rep("#222222", 3), bg = rep(sup_col, 3), ptcex = sup_cex * (0.35 + 0.65 * v / 100))
        }
      }
      list(layout = .vz_or(input$layout, "phylogram"),
           cex = .vz_or(input$cex, 0.9), lwd = .vz_or(input$lwd, 1.5), edge_col = "#222222",
           italic = isTRUE(input$italic),
           tipcol = if (is.null(lc)) NULL else lc$cols[match(lc$vals, lc$lv)],
           sym = if (is.null(sy)) NULL else list(pch = sy$pch[match(sy$vals, sy$lv)], col = sy$cols[match(sy$vals, sy$lv)]),
           blk = if (is.null(bk)) NULL else list(col = bk$cols[match(bk$vals, bk$lv)], grp = bk$vals),
           blk_alpha = alpha, blk_mode = .vz_or(input$blk_mode, "tip_label"),
           blk_grad = isTRUE(input$blk_grad), blk_rev = isTRUE(input$blk_rev),
           sym_size = .vz_or(input$sym_size, 1.2), gap = .vz_or(input$gap, 3) / 100,
           sym_where = .vz_or(input$sym_where, "tip"), sym_align = isTRUE(input$sym_align),
           sym_off = .vz_or(input$sym_off, 0), sym_out = si$sym_out,
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
    shiny::observeEvent(input$leg_pos, { rv$leg_xy <- list() }, ignoreInit = TRUE)
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
  })
}
