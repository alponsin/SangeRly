# ==============================================================================
# Build the offline Windows package of the app (portable R + packages + app + launcher).
# Run ONCE on a Windows PC that has internet and in which the app already works:
#   setwd("D:/alponsin/sangerly")          # the folder that contains app.R
#   source("desktop/build_portable.R")
# It also bundles trimAl and IQ-TREE (see "External programs" below).
# Result: dist/Sangerly/   (double-click Sangerly.vbs; copy the folder anywhere)
#         dist/installer.iss  (optional: compile it with Inno Setup to get Sangerly_Setup.exe)
# ==============================================================================
stopifnot(.Platform$OS.type == "windows")
app_dir <- normalizePath(getwd(), winslash = "/")
if (!file.exists(file.path(app_dir, "app.R"))) stop("Run this from the folder that contains app.R.")
tpl  <- file.path(app_dir, "desktop")
out  <- file.path(app_dir, "dist")
root <- file.path(out, "Sangerly")
# ---- External programs to bundle (edit here if you move them) -----------------------------------------
# Defaults: inside the project folder. If the file is not there, the script searches the folder.
find_exe <- function(default, pattern) {
  if (file.exists(default)) return(normalizePath(default, winslash = "/"))
  f <- list.files(app_dir, pattern, recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
  f <- f[!dir.exists(f) & !grepl("/(dist|desktop)/", f)]
  if (length(f)) normalizePath(sort(f, decreasing = TRUE)[1], winslash = "/") else ""
}
trimal_exe <- find_exe(file.path(app_dir, "trimAL", "trimal.exe"), "^trimal\\.exe$")
iqtree_exe <- find_exe(file.path(app_dir, "iqtree-3.0.1-Windows", "bin", "iqtree3.exe"), "^iqtree[0-9]*\\.exe$")
cat("trimAl :", if (nzchar(trimal_exe)) trimal_exe else "NOT FOUND (the app will have no trimAl)", "\n")
cat("IQ-TREE:", if (nzchar(iqtree_exe)) iqtree_exe else "NOT FOUND (the app will have no IQ-TREE)", "\n")

for (f in c("launch.R", "Sangerly.vbs", "Sangerly_debug.bat", "installer.iss"))
  if (!file.exists(file.path(tpl, f))) stop("Missing desktop/", f)

# ---- which packages does the app use? ----------------------------------------------------
files <- c(file.path(app_dir, "app.R"), list.files(file.path(app_dir, "R"), "\\.R$", full.names = TRUE))
txt <- unlist(lapply(files, readLines, warn = FALSE))
used <- unique(c(unlist(regmatches(txt, gregexpr("[A-Za-z][A-Za-z0-9.]*(?=::)", txt, perl = TRUE))),
                 unlist(regmatches(txt, gregexpr("(?<=library\\()[A-Za-z0-9.]+", txt, perl = TRUE)))))
need <- unique(c("shiny", "bslib", "ggplot2", "DT", "httr", "xml2", "jsonlite", "htmltools",
                 "sangerseqR", "Biostrings", "pwalign", "BiocGenerics", "fontawesome", "later", "httpuv", used))
ip <- installed.packages()
base_pk <- rownames(ip)[ip[, "Priority"] %in% c("base", "recommended")]
need <- setdiff(need, base_pk)
missing <- need[!need %in% rownames(ip)]
# names found by the scan that are not packages (e.g. local objects before '::') are ignored;
# but the core ones must exist
core <- c("shiny", "bslib", "ggplot2", "DT", "httr", "xml2", "jsonlite", "sangerseqR", "Biostrings", "pwalign", "BiocGenerics")
opt <- c("msa", "DECIPHER", "ape", "phangorn", "processx", "writexl")
if (any(opt %in% missing))
  warning("Not installed, so NOT bundled: ", paste(intersect(opt, missing), collapse = ", "),
          "\n  Install them (BiocManager::install(c(\"msa\", \"DECIPHER\")); install.packages(c(\"ape\", \"phangorn\", \"processx\"))) and build again.",
          call. = FALSE)
if (any(core %in% missing)) stop("Install these packages first (BiocManager::install): ", paste(intersect(core, missing), collapse = ", "))
need <- setdiff(need, missing)
deps <- unlist(tools::package_dependencies(need, db = ip, which = c("Depends", "Imports", "LinkingTo"), recursive = TRUE))
all_pk <- setdiff(unique(c(need, deps)), base_pk)
all_pk <- all_pk[all_pk %in% rownames(ip)]
cat(length(all_pk), "packages will be bundled.\n")

# ---- 1. copy R -------------------------------------------------------------------------------
unlink(out, recursive = TRUE); dir.create(root, recursive = TRUE)
rhome <- normalizePath(R.home(), winslash = "/")
cat("1/5 Copying R", as.character(getRversion()), "from", rhome, "...\n")
dir.create(file.path(root, "R"))
file.copy(list.files(rhome, full.names = TRUE), file.path(root, "R"), recursive = TRUE)
unlink(file.path(root, "R", c("doc", "tests")), recursive = TRUE)

# ---- 2. copy packages ------------------------------------------------------------------------
cat("2/5 Copying packages ...\n")
lib <- file.path(root, "R", "library")
for (p in all_pk) {
  if (!dir.exists(file.path(lib, p))) file.copy(find.package(p), lib, recursive = TRUE)
}

# ---- 3. copy the app -------------------------------------------------------------------------------
cat("3/5 Copying the app ...\n")
dir.create(file.path(root, "app"))
items <- setdiff(list.files(app_dir), c("dist", "desktop", "rsconnect", "renv", "installer"))
items <- items[!grepl("\\.Rproj$|^\\.", items)]
items <- items[!grepl("^(iqtree|trimal)", items, ignore.case = TRUE)]   # bundled separately below
items <- items[!grepl("^make_logos", items)]                            # one-time helper, not needed in the package
file.copy(file.path(app_dir, items), file.path(root, "app"), recursive = TRUE)
unlink(file.path(root, "app", "www", "img", "original"), recursive = TRUE)   # full-size logo originals are not needed
if (!file.exists(file.path(root, "app", "www", "img", "logo_assembly.png")))
  warning("www/img/logo_*.png were not copied: the logos will be missing in the package.", call. = FALSE)

# ---- 4. launcher ---------------------------------------------------------------------------------------
cat("4/5 Writing the launcher ...\n")
file.copy(file.path(tpl, c("launch.R", "Sangerly.vbs", "Sangerly_debug.bat")), root)
file.copy(file.path(tpl, "installer.iss"), out)
if (file.exists(file.path(tpl, "Sangerly.ico"))) file.copy(file.path(tpl, "Sangerly.ico"), root)

# ---- External programs: copied into the app folder, where the modules look for them -------------------------
#   trimAl  -> app/trimAL/trimal.exe          IQ-TREE -> app/IQTREE/bin/iqtree3.exe (+ its DLLs)
app_root <- file.path(root, "app")
if (nzchar(trimal_exe)) {
  dir.create(file.path(app_root, "trimAL"), showWarnings = FALSE)
  src <- list.files(dirname(trimal_exe), pattern = "^trimal\\.exe$|\\.dll$|^LICENSE|^COPYING", full.names = TRUE, ignore.case = TRUE)
  file.copy(src, file.path(app_root, "trimAL"), overwrite = TRUE)
}
if (nzchar(iqtree_exe)) {
  dir.create(file.path(app_root, "IQTREE", "bin"), recursive = TRUE, showWarnings = FALSE)
  file.copy(list.files(dirname(iqtree_exe), full.names = TRUE), file.path(app_root, "IQTREE", "bin"),
            recursive = TRUE, overwrite = TRUE)
  lic <- list.files(dirname(dirname(iqtree_exe)), pattern = "^(LICENSE|COPYING|README)", full.names = TRUE, ignore.case = TRUE)
  if (length(lic)) file.copy(lic, file.path(app_root, "IQTREE"), overwrite = TRUE)
}

# ---- 5. smoke test with the portable R only ---------------------------------------------------------------
cat("5/5 Testing the portable R (no internet needed) ...\n")
rs <- file.path(root, "R", "bin", "Rscript.exe")
pk_test <- c("shiny", "bslib", "ggplot2", "DT", "httr", "xml2", "jsonlite", "sangerseqR", "Biostrings", "pwalign",
             "BiocGenerics", "later", "httpuv", intersect(c("msa", "DECIPHER", "ape", "phangorn", "processx", "writexl"), all_pk))
tf <- tempfile(fileext = ".R")
writeLines(c(sprintf('.libPaths("%s")', file.path(root, "R", "library")),
             sprintf('for (p in %s) {', paste(deparse(pk_test), collapse = "")),
             '  r <- try(library(p, character.only = TRUE), silent = TRUE)',
             '  if (inherits(r, "try-error")) cat("FAILED:", p, conditionMessage(attr(r, "condition")), "\\n")',
             '}',
             'cat("SMOKE_OK\\n")'), tf)
res <- suppressWarnings(system2(rs, c("--vanilla", shQuote(tf, type = "cmd")), stdout = TRUE, stderr = TRUE,
                                env = c("R_LIBS=", "R_LIBS_USER=")))
if (any(grepl("SMOKE_OK", res)) && !any(grepl("^FAILED:", res))) cat("OK: the portable R loads every package.\n") else {
  cat("PROBLEM - the portable R could not load all packages:\n"); cat(utils::tail(res, 15), sep = "\n")
}
# do the bundled programs start? (a missing DLL would show here)
try_exe <- function(label, exe, arg) {
  if (!file.exists(exe)) { cat(sprintf("%-8s not bundled\n", label)); return(invisible()) }
  r <- tryCatch(suppressWarnings(system2(exe, arg, stdout = TRUE, stderr = TRUE)), error = function(e) conditionMessage(e))
  cat(sprintf("%-8s %s\n", label, paste(utils::head(r[nzchar(r)], 2), collapse = " | ")))
}
try_exe("trimAl", file.path(app_root, "trimAL", "trimal.exe"), "--version")
try_exe("IQ-TREE", file.path(app_root, "IQTREE", "bin", basename(iqtree_exe)), "--version")
sz <- sum(file.info(list.files(root, recursive = TRUE, full.names = TRUE))$size, na.rm = TRUE)
cat(sprintf("\nDone: %s  (%.0f MB)\nTest: switch the internet OFF and double-click dist/Sangerly/Sangerly.vbs\n", root, sz / 1024^2))
