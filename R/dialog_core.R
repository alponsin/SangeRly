# SangeRly - native "choose a folder / choose files" windows that appear IN FRONT of the app.
#
# utils::choose.dir() / utils::choose.files() open their window without a visible owner. In the desktop version R runs
# hidden (Rscript started by Sangerly.vbs), so Windows keeps that dialog hidden or behind the app window. Here the
# dialog is opened by PowerShell with a visible, invisible-looking top-most owner window, so it comes to the front.
# At the end of this file utils::choose.dir / utils::choose.files are REPLACED (for this R session only) by these
# functions, so every "Browse..." button of the program uses them without changing the other files.
# What happens is written to  %LOCALAPPDATA%\Sangerly\dialog.log  (read it if a window still does not appear).

.dlg <- new.env(parent = emptyenv())

.dlg_log <- function(...) {
  d <- file.path(Sys.getenv("LOCALAPPDATA", tempdir()), "Sangerly")
  try({ dir.create(d, showWarnings = FALSE, recursive = TRUE)
        cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), ..., "\n", file = file.path(d, "dialog.log"), append = TRUE) }, silent = TRUE)
}

.ps_run <- function(body, env) {
  if (.Platform$OS.type != "windows") return(NULL)
  ps <- tempfile("sg_dialog_", fileext = ".ps1"); errf <- tempfile("sg_dialog_err_", fileext = ".txt")
  writeLines(c(
    "$ErrorActionPreference = 'Stop'",
    "try {",
    "  Add-Type -AssemblyName System.Windows.Forms",
    "  Add-Type -AssemblyName System.Drawing",
    "  [System.Windows.Forms.Application]::EnableVisualStyles()",
    "  [Console]::OutputEncoding = New-Object System.Text.UTF8Encoding $false",
    "  $owner = New-Object System.Windows.Forms.Form",
    "  $owner.TopMost = $true; $owner.ShowInTaskbar = $false; $owner.FormBorderStyle = 'None'",
    "  $owner.StartPosition = 'CenterScreen'; $owner.Size = New-Object System.Drawing.Size(1,1); $owner.Opacity = 0",
    "  $owner.Show(); $owner.Activate()",
    paste0("  ", body),
    "} catch { [System.IO.File]::WriteAllText($env:SG_ERR, $_.Exception.ToString()) }",
    "finally { if ($owner) { $owner.Close() } }"), ps, useBytes = TRUE)
  env <- c(env, SG_ERR = errf)
  old <- Sys.getenv(names(env), unset = NA)
  do.call(Sys.setenv, as.list(env))
  on.exit({
    for (i in seq_along(old)) if (is.na(old[i])) Sys.unsetenv(names(old)[i]) else do.call(Sys.setenv, stats::setNames(list(old[[i]]), names(old)[i]))
    unlink(c(ps, errf))
  }, add = TRUE)
  t0 <- Sys.time()
  out <- tryCatch(suppressWarnings(system2("powershell.exe",
           c("-NoProfile", "-STA", "-ExecutionPolicy", "Bypass", "-WindowStyle", "Hidden", "-File", shQuote(ps, type = "cmd")),
           stdout = TRUE, stderr = FALSE, invisible = FALSE)), error = function(e) { .dlg_log("powershell could not start:", conditionMessage(e)); NULL })
  er <- if (file.exists(errf)) paste(readLines(errf, warn = FALSE), collapse = " ") else ""
  .dlg_log(sprintf("powershell finished in %.1f s; %d line(s) of output%s", as.numeric(difftime(Sys.time(), t0, units = "secs")),
                   length(out), if (nzchar(er)) paste0("; ERROR: ", substr(er, 1, 300)) else ""))
  if (is.null(out) || nzchar(er)) return(NULL)
  st <- attr(out, "status"); if (!is.null(st) && st != 0) return(NULL)
  enc2native(out)
}

# a folder; NA when cancelled
.pick_dir <- function(caption = "Choose a folder", default = "", ...) {
  start <- if (!is.null(default) && nzchar(default) && dir.exists(default)) normalizePath(default, winslash = "\\") else ""
  out <- .ps_run(c(
    "$d = New-Object System.Windows.Forms.FolderBrowserDialog",
    "$d.Description = $env:SG_CAPTION; $d.ShowNewFolderButton = $true",
    "if ($env:SG_START -and (Test-Path -LiteralPath $env:SG_START)) { $d.SelectedPath = $env:SG_START }",
    "if ($d.ShowDialog($owner) -eq [System.Windows.Forms.DialogResult]::OK) { Write-Output $d.SelectedPath }"),
    c(SG_CAPTION = if (is.null(caption)) "Choose a folder" else caption, SG_START = start))
  if (is.null(out)) {                                    # PowerShell not usable: the old dialog
    f <- .dlg$orig_dir; if (is.null(f)) f <- utils::choose.dir
    return(tryCatch(f(caption = caption, default = default), error = function(e) NA_character_))
  }
  out <- out[nzchar(out)]
  if (length(out)) out[1] else NA_character_
}

# one or several files; character(0) when cancelled. filters: like choose.files (a matrix: description, pattern)
.pick_files <- function(caption = "Choose files", filters = NULL, multi = TRUE, ...) {
  flt <- if (is.null(filters)) "All files (*.*)|*.*" else
    paste(apply(matrix(filters, ncol = 2), 1, function(r) paste(r, collapse = "|")), collapse = "|")
  out <- .ps_run(c(
    "$d = New-Object System.Windows.Forms.OpenFileDialog",
    "$d.Title = $env:SG_CAPTION; $d.Filter = $env:SG_FILTER; $d.Multiselect = ($env:SG_MULTI -eq '1')",
    "if ($d.ShowDialog($owner) -eq [System.Windows.Forms.DialogResult]::OK) { $d.FileNames | ForEach-Object { Write-Output $_ } }"),
    c(SG_CAPTION = if (is.null(caption)) "Choose files" else caption, SG_FILTER = flt, SG_MULTI = if (isTRUE(multi)) "1" else "0"))
  if (is.null(out)) {
    f <- .dlg$orig_files; if (is.null(f)) f <- utils::choose.files
    return(tryCatch(if (is.null(filters)) f(caption = caption, multi = multi) else f(caption = caption, multi = multi, filters = matrix(filters, ncol = 2)),
                    error = function(e) character()))
  }
  out[nzchar(out)]
}

# replace the two utils functions in this R session (the modules call utils::choose.dir / utils::choose.files)
.install_dialog_shim <- function() {
  if (.Platform$OS.type != "windows" || isTRUE(.dlg$installed)) return(invisible(isTRUE(.dlg$installed)))
  .dlg$orig_dir <- utils::choose.dir; .dlg$orig_files <- utils::choose.files
  ok <- tryCatch({
    utils::assignInNamespace("choose.dir", function(default = "", caption = NULL) .pick_dir(caption = caption, default = default), "utils")
    utils::assignInNamespace("choose.files", function(default = "", caption = "Select files", multi = TRUE, filters = NULL, index = 1)
      .pick_files(caption = caption, filters = filters, multi = multi), "utils")
    TRUE }, error = function(e) { .dlg_log("could not replace the utils dialogs:", conditionMessage(e)); FALSE })
  .dlg$installed <- ok
  .dlg_log("dialog shim installed:", ok)
  invisible(ok)
}
.install_dialog_shim()
