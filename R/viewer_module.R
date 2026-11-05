# ==============================================================================
# Sangerly - interactive alignment viewer (R side)
# ------------------------------------------------------------------------------
# The viewer itself is JavaScript (www/sangerly_viewer.js): it draws the alignment and
# the chromatograms on a canvas, so panning (drag), zooming and selecting are instant.
# This file only (a) provides the HTML skeleton and CSS and (b) turns the alignment into
# the two messages the viewer understands:
#   sv_static : everything that does not change while editing (bases, Phred, traces)
#   sv_dyn    : everything that changes when you delete or decide (consensus, flags, masks)
# Needs BASE_COLS / TRACE_COLS from assembly_module.R (looked up when called).
# ==============================================================================

# The viewer script is embedded INLINE in the page (read from www/sangerly_viewer.js), so it no
# longer depends on Shiny serving the www/ folder (a common cause of "viewer script did not
# load": app started before www/ existed, app.R run from another working directory, file named
# sangerly_viewer.js.txt because Windows hides extensions, or a cached 404 in the RStudio viewer).
.sv_js_paths <- function()
  c("www/sangerly_viewer.js", "sangerly_viewer.js", "R/sangerly_viewer.js", "../www/sangerly_viewer.js")

.sv_js_tag <- function() {
  p <- .sv_js_paths(); p <- p[file.exists(p)]
  if (length(p)) {
    js <- paste(readLines(p[1], warn = FALSE, encoding = "UTF-8"), collapse = "\n")
    js <- gsub("</script", "<\\/script", js, fixed = TRUE)
    return(tags$script(HTML(js)))
  }
  message("Sangerly: sangerly_viewer.js not found (looked in: ", paste(.sv_js_paths(), collapse = ", "),
          "; working directory: ", getwd(), "). Falling back to src=.")
  tags$script(src = "sangerly_viewer.js")
}

svViewerHead <- function() {
  tagList(
    tags$style(HTML("
      .sv-root{background:#fff;border-radius:8px;box-shadow:0 1px 3px rgba(0,0,0,.12);margin-bottom:12px;overflow:hidden}
      .sv-toolbar{display:flex;flex-wrap:wrap;gap:6px 18px;align-items:center;padding:6px 12px;border-bottom:1px solid #e3e8ee;font-size:12.5px;color:#1f3c5f}
      .sv-toolbar label{display:flex;align-items:center;gap:6px;margin:0;font-weight:500}
      .sv-toolbar input[type=range]{width:110px}
      .sv-info{color:#6b7785;flex:1 1 320px;min-width:0;white-space:nowrap;overflow:hidden;text-overflow:ellipsis}
      .sv-scroll{display:none;position:relative;overflow-x:auto;overflow-y:hidden;outline:none;cursor:grab;
                 -webkit-user-select:none;user-select:none}
      .sv-scroll:focus{box-shadow:inset 0 0 0 2px rgba(44,123,229,.35)}
      .sv-canvas{display:block;position:sticky;left:0;top:0}
      .sv-spacer{position:absolute;left:0;top:0;height:1px;pointer-events:none}
      .sv-empty{padding:28px;text-align:center;color:#6b7785}
    ")),
    .sv_js_tag(),
    # If the script file is missing or blocked, say so in the viewer box instead of staying empty.
    tags$script(HTML("
      setTimeout(function () {
        if (window.SangerlyViewer) return;
        var els = document.querySelectorAll('.sv-empty');
        for (var i = 0; i < els.length; i++) {
          els[i].style.color = '#c0392b';
          els[i].innerHTML = '<b>The viewer script did not load.</b> Make sure <code>www/sangerly_viewer.js</code> exists ' +
            'next to app.R (the file name must end in .js, not .js.txt), restart the app, and open it in a ' +
            'normal browser (Open in Browser) with Ctrl+F5. Details: press F12 and open the Console tab.';
        }
      }, 2500);
    "))
  )
}

svViewerUI <- function(id) {
  div(id = id, class = "sv-root",
      div(class = "sv-toolbar",
          tags$label("Zoom", tags$input(type = "range", class = "sv-zoom", min = 4, max = 48, step = 1, value = 22)),
          tags$label(tags$input(type = "checkbox", class = "sv-trace", checked = NA), "Chromatograms"),
          tags$label("Trace height", tags$input(type = "range", class = "sv-th", min = 30, max = 180, step = 2, value = 72)),
          span(class = "sv-info", "")),
      div(class = "sv-scroll", tabindex = 0, tags$canvas(class = "sv-canvas"), div(class = "sv-spacer")),
      div(class = "sv-empty", "Assemble first to see the alignment here."))
}

# One read: trace samples placed on the alignment-column axis (as in the old ggplot viewer),
# thinned to ~1 sample per pixel; channels scaled to 0-255.
sv_trace <- function(aln, i) {
  rd <- aln$reads[[i]]
  kk <- which(!is.na(aln$idx[i, ]))
  if (length(kk) < 2) return(NULL)
  s <- rd$ploc[aln$idx[i, kk]]
  o <- order(s); s <- s[o]; ko <- kk[o]
  if (length(unique(s)) < 2) return(NULL)
  sp <- stats::median(abs(diff(s))); if (!is.finite(sp) || sp <= 0) sp <- 10
  sc <- seq(floor(min(s) - sp / 2), ceiling(max(s) + sp / 2))
  sc <- sc[sc >= 1 & sc <= nrow(rd$traces)]
  if (!length(sc)) return(NULL)
  xs <- suppressWarnings(stats::approx(s, ko, xout = sc, rule = 1, ties = mean)$y)
  slope <- (ko[length(ko)] - ko[1]) / (s[length(s)] - s[1])
  lo <- is.na(xs) & sc < s[1]; hi <- is.na(xs) & sc > s[length(s)]
  xs[lo] <- ko[1] + (sc[lo] - s[1]) * slope
  xs[hi] <- ko[length(ko)] + (sc[hi] - s[length(s)]) * slope
  keep <- if (sp >= 8) seq(1, length(sc), by = 2) else seq_along(sc)
  sc <- sc[keep]; xs <- xs[keep]
  mx <- max(rd$traces[sc, , drop = FALSE]); if (!is.finite(mx) || mx <= 0) mx <- 1
  out <- list(x = I(round(xs, 2)))
  for (ch in c("A", "C", "G", "T")) {
    v <- if (ch %in% colnames(rd$traces)) rd$traces[sc, ch] else rep(0, length(sc))
    out[[ch]] <- I(as.integer(round(255 * v / mx)))
  }
  out
}

sv_static <- function(aln, id) {
  nc <- ncol(aln$base); nr <- nrow(aln$base)
  reads <- lapply(seq_len(nr), function(i) {
    b <- aln$base[i, ]; b[is.na(b)] <- " "
    q <- aln$q[i, ]; q[is.na(q)] <- -1
    list(label = aln$labels[i], base = paste(b, collapse = ""), q = I(as.integer(round(q))), tr = sv_trace(aln, i))
  })
  has_ref <- !is.null(aln$ref)
  list(id = id, nc = nc, nr = nr, has_ref = has_ref,
       ref = if (has_ref) paste(aln$refbase, collapse = "") else "",
       refpos = if (has_ref && !anyNA(aln$refpos)) I(as.integer(aln$refpos)) else NULL,
       colors = list(base = as.list(BASE_COLS), trace = as.list(TRACE_COLS)),
       reads = reads)
}

sv_dyn <- function(cf, mask, delcol, id) {
  flag <- function(x) paste(as.integer(x), collapse = "")
  list(id = id,
       final = paste(cf$Final, collapse = ""),
       review = flag(cf$Review & !cf$Manual),
       manual = flag(cf$Manual),
       delcol = flag(delcol),
       mask = I(unname(vapply(seq_len(nrow(mask)), function(i) flag(mask[i, ]), ""))))
}
