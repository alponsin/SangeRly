# Starts the Shiny app of SangeRly for the Qt window:  Rscript shiny_host.R <folder with app.R> <port>
# The app is the same one as in the Shiny program (every tool); the Qt window shows it with ?embed=1.
args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(args[1], winslash = "/", mustWork = TRUE)
port <- as.integer(args[2])
Sys.setenv(SANGERLY_DESKTOP = "1")
options(shiny.maxRequestSize = Inf, warn = 1)
logdir <- file.path(Sys.getenv("LOCALAPPDATA", tempdir()), "Sangerly"); dir.create(logdir, showWarnings = FALSE, recursive = TRUE)
logfile <- file.path(logdir, "session.log")
sglog <- function(...) try(cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0("[", Sys.getpid(), "]"), ..., "\n", file = logfile, append = TRUE), silent = TRUE)
options(shiny.error = function() sglog("ERROR in the app:", geterrmessage()))
# keep the compiled Bootstrap theme between launches (only the first start is slow)
try(options(sass.cache = sass::sass_file_cache(file.path(logdir, "sass_cache"))), silent = TRUE)

# "Browse folder..." (utils::choose.dir) opens behind the window when R runs hidden: use a top-most Windows dialog instead
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


setwd(root)
sglog("shiny_host: starting the app on port", port, "from", root)
shiny::runApp(root, port = port, host = "127.0.0.1", launch.browser = FALSE)
