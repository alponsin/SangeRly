# SangeRly desktop launcher (runs inside the portable R; started by Sangerly.vbs)
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
Sys.setenv(SANGERLY_DESKTOP = "1")
options(shiny.maxRequestSize = 100 * 1024^2, warn = 1)
options(shiny.error = function() sglog("ERROR in the app:", geterrmessage()))
sglog("---- launcher started;", R.version.string)

# ---- folder dialog that always comes to the front ---------------------------------------------------------------------
# "Browse folder..." uses utils::choose.dir(). In the packaged program R runs without a console window (Sangerly.vbs starts it
# hidden) and the app lives in an Edge/Chrome window, so that dialog can open BEHIND the app window or not open at all, and the
# button seems to do nothing. "Browse files..." (utils::choose.files) is a different Windows dialog and is not affected.
# Here choose.dir is replaced (only in this launcher) by the standard Windows folder dialog opened from PowerShell with a
# top-most owner window, so it is always visible. If PowerShell cannot be used, the original function is called.
sg_ps_choose_dir <- function(default = "", caption = "Select folder") {
  ps <- file.path(Sys.getenv("SystemRoot", "C:\\Windows"), "System32", "WindowsPowerShell", "v1.0", "powershell.exe")
  if (!file.exists(ps)) return("__FAIL__")
  out <- tempfile("sg_dir_", fileext = ".txt"); scr <- tempfile("sg_dir_", fileext = ".ps1")
  writeLines(c(
    "Add-Type -AssemblyName System.Windows.Forms; Add-Type -AssemblyName System.Drawing",
    "[System.Windows.Forms.Application]::EnableVisualStyles()",
    "$f = New-Object System.Windows.Forms.Form",
    "$f.TopMost = $true; $f.ShowInTaskbar = $false; $f.StartPosition = 'Manual'",
    "$f.Location = New-Object System.Drawing.Point(-3000, -3000); $f.Size = New-Object System.Drawing.Size(1, 1)",
    "$f.Show(); $f.Activate()",
    "$d = New-Object System.Windows.Forms.FolderBrowserDialog",
    "$d.Description = $env:SG_CAPTION; $d.ShowNewFolderButton = $true",
    "if ($env:SG_DEFAULT -and (Test-Path -LiteralPath $env:SG_DEFAULT)) { $d.SelectedPath = $env:SG_DEFAULT }",
    "$r = $d.ShowDialog($f)",
    "$u = New-Object System.Text.UTF8Encoding($false)",
    "if ($r -eq [System.Windows.Forms.DialogResult]::OK) { [System.IO.File]::WriteAllText($env:SG_OUT, $d.SelectedPath, $u) } else { [System.IO.File]::WriteAllText($env:SG_OUT, '__CANCEL__', $u) }",
    "$f.Close()"), scr, useBytes = TRUE)
  on.exit(unlink(c(out, scr)), add = TRUE)
  st <- tryCatch(suppressWarnings(system2(ps, c("-NoProfile", "-NonInteractive", "-STA", "-ExecutionPolicy", "Bypass", "-WindowStyle", "Hidden", "-File", shQuote(scr)),
                                          stdout = FALSE, stderr = FALSE, wait = TRUE,
                                          env = c(paste0("SG_CAPTION=", if (is.null(caption)) "Select folder" else caption),
                                                  paste0("SG_DEFAULT=", if (is.null(default)) "" else default), paste0("SG_OUT=", out)))),
                 error = function(e) { sglog("folder dialog: PowerShell failed:", conditionMessage(e)); 1L })
  if (!file.exists(out)) { sglog("folder dialog: no answer from PowerShell (status", st, ")"); return("__FAIL__") }
  r <- readLines(out, warn = FALSE, encoding = "UTF-8")
  r <- if (length(r)) r[1] else ""
  if (identical(r, "__CANCEL__") || !nzchar(r)) return(NA_character_)
  gsub("/", "\\\\", r)
}
if (.Platform$OS.type == "windows") try({
  .sg_orig_choose_dir <- utils::choose.dir
  .sg_new_choose_dir <- function(default = "", caption = NULL) {
    r <- sg_ps_choose_dir(default, if (is.null(caption)) "Select folder" else caption)
    if (identical(r, "__FAIL__")) .sg_orig_choose_dir(default, caption) else r
  }
  for (e in list(asNamespace("utils"), as.environment("package:utils"))) {
    unlockBinding("choose.dir", e); assign("choose.dir", .sg_new_choose_dir, envir = e); lockBinding("choose.dir", e)
  }
  sglog("folder dialog: utils::choose.dir replaced by the top-most PowerShell dialog")
}, silent = TRUE)

# The app calls this (if it exists) at the start of every session. When the last window has been
# closed for 5 seconds, R stops. A page refresh (F5) does not stop it.
.sg <- new.env(); .sg$n <- 0L; .sg$seen <- FALSE; .sg$t0 <- Sys.time()
assign("sangerly_register_session", function(session) {
  if (isTRUE(session$userData$sg_registered)) return(invisible())   # already counted
  session$userData$sg_registered <- TRUE
  .sg$n <- .sg$n + 1L; .sg$seen <- TRUE
  sglog(sprintf("window connected (%.0f s after the launcher started)", as.numeric(difftime(Sys.time(), .sg$t0, units = "secs"))))
  session$onSessionEnded(function() {
    .sg$n <- .sg$n - 1L
    sglog("window disconnected (R stops in 2 minutes unless a window comes back)")
    # 2 minutes of grace: a page refresh or a short sleep does not kill R; a closed window frees the memory soon.
    later::later(function() if (.sg$n <= 0L) { sglog("no window came back: stopping R"); shiny::stopApp() }, 120)
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

# ---- splash window: shown immediately, jumps to the app as soon as the app answers -------------------------
file_url <- function(p) paste0("file:///", utils::URLencode(normalizePath(p, winslash = "/"), reserved = FALSE))
icon_url <- if (file.exists(file.path(root, "Sangerly.ico"))) file_url(file.path(root, "Sangerly.ico")) else ""
splash_tpl <- '<!DOCTYPE html><html><head><meta charset="utf-8"><title>SangeRly</title>
<link rel="icon" href="@ICON@">
<style>
html,body{height:100%;margin:0}
body{background:radial-gradient(1200px 700px at 50% 38%,#ffffff 0%,#f4f8ff 55%,#e9f1fd 100%);color:#1f3c5f;
     font-family:"Segoe UI",system-ui,Arial,sans-serif;display:flex;flex-direction:column;align-items:center;justify-content:center;overflow:hidden}
.logo{width:min(34vh,260px);height:min(34vh,260px);
      filter:drop-shadow(0 16px 26px rgba(11,61,145,.30));animation:pop .7s ease-out both,float 3.2s ease-in-out .7s infinite}
.logo img{width:100%;height:100%;display:block}
@keyframes pop{from{opacity:0;transform:scale(.82)}to{opacity:1;transform:scale(1)}}
@keyframes float{0%,100%{transform:translateY(0)}50%{transform:translateY(-8px)}}
h1{margin:34px 0 4px;font-weight:700;font-size:clamp(34px,6vh,52px);letter-spacing:.5px;
   background:linear-gradient(90deg,#0b3d91,#14a3b8);-webkit-background-clip:text;background-clip:text;color:transparent}
.tag{margin:0;font-size:clamp(14px,2.2vh,18px);color:#4a6382;text-align:center;padding:0 20px}
.bar{width:min(380px,70vw);height:6px;background:#dbe6f5;border-radius:6px;margin:34px 0 12px;overflow:hidden;position:relative}
.bar i{position:absolute;top:0;left:-40%;width:40%;height:100%;border-radius:6px;background:linear-gradient(90deg,#0b3d91,#14a3b8,#22b07d);animation:slide 1.3s ease-in-out infinite}
@keyframes slide{0%{left:-40%}100%{left:100%}}
#msg{font-size:14px;color:#6b7f99;min-height:20px;transition:opacity .25s}
.foot{position:fixed;bottom:16px;left:0;right:0;text-align:center;font-size:12px;color:#8a9bb2}
</style></head><body>
<div class="logo"><img src="@ICON@" alt="" onerror="this.parentNode.style.display=\'none\'"></div>
<h1>SangeRly</h1>
<p class="tag">Sanger assembly &middot; BLAST &middot; alignment &middot; phylogeny &middot; population genetics</p>
<div class="bar"><i></i></div>
<div id="msg">Starting, please wait...</div>
<div class="foot">Works offline &middot; the first start can take a little longer</div>
<script>
var target="@URL@", t0=Date.now();
var msgs=["Starting, please wait...","Loading the analysis packages...","Preparing the tools...","Almost ready..."];
setInterval(function(){var m=document.getElementById("msg"),i=Math.min(msgs.length-1,Math.floor((Date.now()-t0)/3500));if(m.textContent!==msgs[i]){m.style.opacity=0;setTimeout(function(){m.textContent=msgs[i];m.style.opacity=1;},250);}},500);
function go(){location.replace(target);}
function poll(){
  fetch(target,{mode:"no-cors",cache:"no-store"}).then(go).catch(function(){
    if(Date.now()-t0>240000) go(); else setTimeout(poll,500);});
}
poll();
</script></body></html>'
open_splash <- function(p) {
  html <- gsub("@ICON@", icon_url, splash_tpl, fixed = TRUE)
  html <- gsub("@URL@", paste0("http://127.0.0.1:", p, "/"), html, fixed = TRUE)
  f <- tempfile("sangerly_splash_", fileext = ".html")
  writeLines(html, f, useBytes = TRUE)
  open_window(file_url(f))
}

# ---- double click protection ------------------------------------------------------------------------------
# Only ONE thing is checked: a copy that was started a few seconds ago (the same double click). Its splash window
# is already on the screen, so this second start simply ends. Every other start is a fresh, independent copy
# (starting takes only a few seconds): no window is ever sent to an old copy, and no old file can block a start.
inst_file <- file.path(logdir, "instance.txt")       # line 1: port, line 2: start time, line 3: process id
Sys.sleep(stats::runif(1, 0, 0.3))
if (file.exists(inst_file)) {
  info <- suppressWarnings(readLines(inst_file, warn = FALSE))
  age <- as.numeric(Sys.time()) - suppressWarnings(as.numeric(info[2]))
  if (!is.na(age) && age >= 0 && age < 8) {
    sglog("a copy was started", round(age, 1), "s ago (double click): this start ends")
    quit(save = "no", status = 0)
  }
}

port <- httpuv::randomPort()
writeLines(c(as.character(port), as.character(as.numeric(Sys.time())), as.character(Sys.getpid())), inst_file)
open_splash(port)

# If the splash window is closed before the app ever started, do not leave R running.
later::later(function() if (!.sg$seen) shiny::stopApp(), 900)   # 15 min: a slow first start must not be cut short

# Keep the compiled Bootstrap theme between launches (only the first start is slow)
try(options(sass.cache = sass::sass_file_cache(file.path(logdir, "sass_cache"))), silent = TRUE)

res <- tryCatch({
  library(shiny)
  sglog(sprintf("shiny loaded (%.0f s after the launcher started)", as.numeric(difftime(Sys.time(), .sg$t0, units = "secs"))))
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
