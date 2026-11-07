# SangeRly

**SangeRly** (Sanger + R) is an R/Shiny application for Sanger sequencing data. It takes you from raw chromatograms (`.ab1`) to
assembled consensus sequences, species identification, multiple alignment, phylogenetic trees and publication-ready
figures, in one window. It runs as a normal Shiny app, and as an **offline Windows desktop program**
(portable R + installer, no R installation needed).

> Status: research software under active development. Please check results (especially base calls at conflicting
> positions) against the chromatograms before using them in a publication.

<!-- Add screenshots here, for example: ![Home](docs/home.png) -->

## Workflow

| Step | What it does |
|---|---|
| **1. Assembly** | Reads `.ab1` files, trims low-quality ends (modified Mott or sliding window), assembles forward + reverse reads (or many reads / tiled amplicons against a reference), shows the alignment with chromatograms, lets you edit or delete bases, calls SNPs, translates, and exports consensus, reports and tables. Includes composition and restriction-site analysis. |
| **2. BLAST** | Sends the consensus to NCBI BLAST (or opens the NCBI website with the form filled in), with Descriptions, Graphic Summary, Alignments and Taxonomy tabs. The top hits can be sent straight to the MSA tab. |
| **3. MSA and trimming** | Aligns sequences (Clustal Omega, MUSCLE, ClustalW or DECIPHER), lets you edit the alignment, trims it with trimAl, and reports conserved / variable / parsimony-informative / singleton sites (Excel export). |
| **4. Phylogeny** | Neighbor-Joining trees with several distance models and Maximum Likelihood trees with IQ-TREE (ModelFinder, ultrafast bootstrap, SH-aLRT). Distance matrices export to Excel, with heat maps as SVG / JPG / PNG / PDF. |
| **5. Visualization** | Tree and data figures in the style of iTOL / ChiPlot. |

Everything works **offline** except BLAST, fetching a reference by NCBI accession, taxon-name validation and the
NCBI Taxonomy lookup, which need an internet connection.

## Install (Windows desktop version)

1. Download `Sangerly_Setup.exe` (or the portable zip) from the **Releases** page.
2. Run it. It installs for the current user (no administrator rights) and adds a Start-menu entry.
3. Windows SmartScreen may warn about an unsigned program: choose **More info -> Run anyway**.
4. Start SangeRly and wait for the splash window (the first start takes 20-60 s).

Requirements: 64-bit Windows 10/11 and Microsoft Edge or Chrome (Edge is built into Windows).

## Run from source

Install R (4.4 or newer) and the packages:

```r
install.packages(c("shiny", "bslib", "DT", "ggplot2", "httr", "xml2", "jsonlite", "htmltools",
                   "ape", "phangorn", "processx", "writexl", "BiocManager"))
BiocManager::install(c("sangerseqR", "Biostrings", "pwalign", "BiocGenerics", "msa", "DECIPHER"))
```

Then, from the project folder:

```r
shiny::runApp()
```

### External programs (optional)

Two programs are used when present. They are **not** part of this repository; download them from their official
sites and put them inside the project folder:

| Program | Used for | Where to put it |
|---|---|---|
| [trimAl](https://github.com/inab/trimal) | alignment trimming | `trimAL/trimal.exe` |
| [IQ-TREE](https://iqtree.org) | Maximum Likelihood trees, ModelFinder | `IQTREE/bin/iqtree3.exe` (any folder named `iqtree*` also works) |

Without them, the app still works: the built-in trimmers and Neighbor-Joining remain available. In the IQ-TREE panel
you can also type the full path of `iqtree3.exe`.

## Build the Windows desktop version

1. Put trimAl and IQ-TREE in the project folder as above (they are copied into the build).
2. Run from R:
   ```r
   setwd("path/to/sangerly"); source("desktop/build_portable.R")
   ```
   This copies R, the packages, the app and the two programs into `dist/Sangerly/` and tests that everything loads.
3. Test `dist/Sangerly/Sangerly.vbs` with the internet switched off.
4. Optional installer: open `dist/installer.iss` in [Inno Setup](https://jrsoftware.org/isinfo.php) (free) and press
   Compile. The result is `dist/Output/Sangerly_Setup.exe`.

Details are in [`desktop/README_desktop.md`](desktop/README_desktop.md).

## Troubleshooting (desktop version)

Log files are in `%LOCALAPPDATA%\Sangerly\`:

* `session.log`: when windows connect and disconnect, and any error inside the app. Read this first.
* `last_error.txt`: why the app could not start.
* `Sangerly_debug.bat` (in the install folder) runs the app with a console window that stays open.

If a window does not open, wait for the splash (the first start on a new PC can take a minute) instead of clicking
again. If an old R process is stuck, end `Rscript.exe` in Task Manager and start SangeRly again.

## Credits and licenses

SangeRly is released under the **GNU General Public License v3.0** (see [`LICENSE`](LICENSE)).
Third-party software (R packages, trimAl, IQ-TREE) and their licenses are listed in
[`THIRD_PARTY_NOTICES.md`](THIRD_PARTY_NOTICES.md).

SangeRly reads chromatograms with **sangerseqR**, aligns reads with **Biostrings / pwalign**, and uses
**msa**, **DECIPHER**, **ape**, **phangorn**, **trimAl** and **IQ-TREE** for alignment and phylogeny. Its trimming
options and HTML report are inspired by **sangeranalyseR**, which SangeRly does not use as a dependency.
If you use SangeRly in a publication, please also cite the tools you used through it.

## How to cite

> Alponsin and Putra, A.E. SangeRly: an R/Shiny application for Sanger sequence assembly, identification and phylogeny
> version 1.0, 2026. https://github.com/alponsin/sangerly

(Replace with the final reference when a paper or a Zenodo DOI is available.)
