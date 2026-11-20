# Replaces utils::choose.dir() / utils::choose.files() by the dialogs of R/dialog_core.R in the R files of the app,
# so that the "Browse..." windows open IN FRONT of the app (desktop version).
# Run it from the project folder (the one with app.R):
#   source("fix_dialogs.R")                              # dry run: only lists the lines
#   dry_run <- FALSE; source("fix_dialogs.R")            # applies it (originals are kept in _backup_dialogs/)
if (!exists("dry_run")) dry_run <- TRUE
files <- list.files("R", pattern = "\\.R$", full.names = TRUE)
files <- files[!grepl("dialog_core\\.R$", files)]
n <- 0L
for (f in files) {
  x <- readLines(f, warn = FALSE)
  y <- gsub("utils::choose.dir(", ".pick_dir(", x, fixed = TRUE)
  y <- gsub("utils::choose.files(", ".pick_files(", y, fixed = TRUE)
  chg <- which(y != x)
  if (!length(chg)) next
  n <- n + length(chg)
  cat(sprintf("%s: %d line(s)\n", f, length(chg)))
  if (!dry_run) {
    bk <- file.path("_backup_dialogs", f); dir.create(dirname(bk), recursive = TRUE, showWarnings = FALSE)
    file.copy(f, bk, overwrite = FALSE); writeLines(y, f, useBytes = TRUE)
  }
}
cat(sprintf("\n%s %d line(s).%s\n", if (dry_run) "Would change" else "Changed", n,
            if (dry_run) " Nothing written; set dry_run <- FALSE to apply." else ""))
if (!file.exists("R/dialog_core.R")) cat("WARNING: R/dialog_core.R is missing; copy it to the R folder.\n")
if (!any(grepl("dialog_core", readLines("app.R", warn = FALSE)))) cat("Remember to add  source(\"R/dialog_core.R\")  to app.R (before the modules).\n")
