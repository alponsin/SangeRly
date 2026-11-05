# Sangerly desktop launcher (runs inside the portable R; started by Sangerly.vbs)
# Opens a splash window at once, then the app in its own window (Edge/Chrome "app mode"),
# and stops R when that window is closed.
root <- normalizePath(getwd(), winslash = "/")
logdir <- file.path(Sys.getenv("LOCALAPPDATA", tempdir()), "Sangerly")
dir.create(logdir, showWarnings = FALSE, recursive = TRUE)
errfile <- file.path(logdir, "last_error.txt")
if (file.exists(errfile)) unlink(errfile)

# Log of what happens while the app runs (read it if the app closes by itself)
logfile <- file.path(logdir, "session.log")
if (file.exists(logfile) && file.size(logfile) > 2e5) unlink(logfile)
sglog <- function(...) try(cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0("[", Sys.getpid(), "]"), ..., "\n",
                               file = logfile, append = TRUE), silent = TRUE)

.libPaths(file.path(root, "R", "library"))
# Keep the compiled Bootstrap theme between launches (only the first start is slow)
options(sass.cache = sass::sass_file_cache(file.path(logdir, "sass_cache")))
Sys.setenv(SANGERLY_DESKTOP = "1")
options(shiny.maxRequestSize = 100 * 1024^2, warn = 1)
options(shiny.error = function() sglog("ERROR in the app:", geterrmessage()))
sglog("---- launcher started;", R.version.string)

# The app calls this (if it exists) at the start of every session. When the last window has been
# closed for 5 seconds, R stops. A page refresh (F5) does not stop it.
.sg <- new.env(); .sg$n <- 0L; .sg$seen <- FALSE
assign("sangerly_register_session", function(session) {
  if (isTRUE(session$userData$sg_registered)) return(invisible())   # already counted
  session$userData$sg_registered <- TRUE
  .sg$n <- .sg$n + 1L; .sg$seen <- TRUE
  sglog("window connected")
  session$onSessionEnded(function() {
    .sg$n <- .sg$n - 1L
    sglog("window disconnected (R stops in 10 minutes unless a window comes back)")
    # 10 minutes of grace: a refresh, a sleeping tab, a sleeping PC or an error does not kill R.
    # Opening Sangerly again within that time reuses this R (instant start).
    later::later(function() if (.sg$n <= 0L) { sglog("no window came back: stopping R"); shiny::stopApp() }, 600)
  })
}, envir = globalenv())

open_window <- function(url) {
  pf <- c(Sys.getenv("ProgramFiles(x86)"), Sys.getenv("ProgramFiles"), Sys.getenv("LOCALAPPDATA"))
  pf <- pf[nzchar(pf)]
  cands <- c(file.path(pf, "Microsoft/Edge/Application/msedge.exe"),
             file.path(pf, "Google/Chrome/Application/chrome.exe"))
  hit <- cands[file.exists(cands)]
  if (length(hit)) system2(hit[1], c(paste0("--app=", url), "--window-size=1400,900"), wait = FALSE)
  else get("browseURL", envir = asNamespace("utils"))(url)
}

# ---- one copy at a time -----------------------------------------------------------------------------------
# A second start (for example a double click while the first copy is still loading) must not start another R.
inst_file <- file.path(logdir, "instance.txt")
port_alive <- function(p) tryCatch({
  con <- suppressWarnings(socketConnection("127.0.0.1", p, open = "r+b", blocking = TRUE, timeout = 1))
  close(con); TRUE }, error = function(e) FALSE)
Sys.sleep(stats::runif(1, 0, 0.8))      # two clicks at the same moment: let one of them go first
if (file.exists(inst_file)) {
  info <- suppressWarnings(readLines(inst_file, warn = FALSE))
  p0 <- suppressWarnings(as.integer(info[1])); t0 <- suppressWarnings(as.numeric(info[2]))
  starting <- !is.na(t0) && (as.numeric(Sys.time()) - t0) < 60
  if (!is.na(p0) && port_alive(p0)) {
    sglog("another copy is running on port", p0, ": opening a window on it")
    open_window(paste0("http://127.0.0.1:", p0, "/"))
    quit(save = "no", status = 0)
  } else if (starting) {
    sglog("another copy is still starting: this one exits")
    quit(save = "no", status = 0)
  }
}

# ---- splash: shown immediately, jumps to the app as soon as the app answers ---------------
port <- httpuv::randomPort()
writeLines(c(as.character(port), as.character(as.numeric(Sys.time()))), inst_file)
app_url <- paste0("http://127.0.0.1:", port, "/")
file_url <- function(p) paste0("file:///", utils::URLencode(normalizePath(p, winslash = "/"), reserved = FALSE))
icon_url <- if (file.exists(file.path(root, "Sangerly.ico"))) file_url(file.path(root, "Sangerly.ico")) else ""
splash <- '<!DOCTYPE html><html><head><meta charset="utf-8"><title>Sangerly</title>
<link rel="icon" href="@ICON@">
<style>
html,body{height:100%;margin:0}
body{background:#2c7be5;color:#fff;font-family:"Segoe UI",system-ui,Arial,sans-serif;
     display:flex;flex-direction:column;align-items:center;justify-content:center}
img{width:96px;height:96px;margin-bottom:18px}
h1{margin:0 0 6px;font-weight:600;font-size:34px}
p{margin:0 0 26px;opacity:.9}
.s{width:34px;height:34px;border:4px solid rgba(255,255,255,.35);border-top-color:#fff;
   border-radius:50%;animation:r 1s linear infinite}
@keyframes r{to{transform:rotate(360deg)}}
</style></head><body>
<img src="@ICON@" alt="" onerror="this.style.display=\'none\'">
<h1>Sangerly</h1><p>Loading, please wait...</p><div class="s"></div>
<script>
var target="@URL@", t0=Date.now();
function go(){location.replace(target);}
function poll(){
  fetch(target,{mode:"no-cors",cache:"no-store"}).then(go).catch(function(){
    if(Date.now()-t0>90000) go(); else setTimeout(poll,500);});
}
poll();
</script></body></html>'
splash <- gsub("@ICON@", icon_url, splash, fixed = TRUE)
splash <- gsub("@URL@", app_url, splash, fixed = TRUE)
splash_file <- file.path(tempdir(), "sangerly_splash.html")
writeLines(splash, splash_file, useBytes = TRUE)
open_window(file_url(splash_file))

# If the splash window is closed before the app ever started, do not leave R running.
later::later(function() if (!.sg$seen) shiny::stopApp(), 180)

res <- tryCatch({
  library(shiny)
  sglog("launcher: shiny loaded")
  app_dir <- file.path(root, "app")
  # Count every session here, so the launcher does not depend on a line inside app.R
  app <- tryCatch({
    a <- shiny::shinyAppDir(app_dir)
    orig <- a$serverFuncSource
    if (!is.function(orig)) stop("serverFuncSource not found")
    a$serverFuncSource <- function() {
      f <- orig()
      function(input, output, session) {
        sangerly_register_session(session)
        if (length(formals(f)) >= 3) f(input, output, session) else f(input, output)
      }
    }
    a
  }, error = function(e) {
    sglog("could not wrap the app (", conditionMessage(e), "); running it as it is")
    NULL })
  if (is.null(app)) { .sg$seen <- TRUE; app <- app_dir }   # no wrapper: never stop R on a timer
  setwd(app_dir)
  sglog("runApp starting")
  runApp(app, port = port, host = "127.0.0.1", launch.browser = FALSE)
  NULL
}, error = function(e) conditionMessage(e))

sglog("runApp ended", if (!is.null(res)) paste("with error:", res))
try({ if (file.exists(inst_file) && identical(readLines(inst_file, warn = FALSE)[1], as.character(port))) unlink(inst_file) }, silent = TRUE)
if (!is.null(res)) {
  writeLines(c(format(Sys.time()), res, "", paste("R_HOME:", R.home()),
               paste("libPaths:", paste(.libPaths(), collapse = "; ")), paste("working dir:", getwd())), errfile)
  quit(save = "no", status = 1)
}
