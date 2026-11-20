# SangeRly

**SangeRly** is a free desktop program for Sanger sequencing and DNA-barcode work, from the raw `.ab1` chromatograms to a phylogenetic tree and population-genetic statistics. It is written in R (Shiny) and comes in a native Windows window (Qt).

Version **1.0.0** · [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24126657.svg)](https://doi.org/10.5281/zenodo.24126657) · License: GPL-3.0-or-later

<!-- Add a screenshot here: ![SangeRly](docs/screenshot.png) -->

## Features

| Tool | What it does |
|---|---|
| **Assembly** | Assemble forward and reverse `.ab1` reads (de novo, or against a reference), inspect chromatograms, delete low-Phred bases, call SNPs, edit and export the consensus. |
| **BLAST** | Identify the species by searching the consensus against the NCBI databases. |
| **MSA and trimming** | Align many sequences (DECIPHER, Clustal Omega, ClustalW, MUSCLE, reference-guided, or already aligned) and trim poorly aligned columns. |
| **Phylogeny** | Neighbor-joining and maximum-likelihood (IQ-TREE) trees, ModelFinder, an interactive viewer and tree comparison. |
| **Visualisation** | Prepare Newick files for iTOL and ChiPlot and draw figures. |
| **Population structure** | Haplotypes, diversity, Tajima's D, Fu's Fs, networks, maps, AMOVA / PhiST, Structure (K clusters), sliding window, and export for DnaSP, Arlequin and PopART. |
| **Genetic distances** | Distance matrices, barcode gap, PCoA, statistical tests and species identification. |
| **Primer mining and design** | Find target sequences in BOLD and NCBI, design primers, check their quality and BLAST them. |
| **More** | Open and view files, project and report, GenBank export (BankIt / table2asn), tutorial, help, citation. |

Everything works offline except BLAST, fetching a reference by NCBI accession, validating taxon names, and the primer-mining searches.

## Download (Windows)

Get `SangeRly_Setup.exe` from the [**Releases**](https://github.com/alponsin/sangerly/releases) page and run it. It installs for the current user (no administrator rights needed) and adds a Start-menu entry and an optional desktop shortcut. R is included, so you do not need to install anything else.

> **Unsigned installer.** Windows SmartScreen may warn because the installer is not code-signed yet: click **More info**, then **Run anyway**. The SHA-256 checksum of each release is given on the Releases page.

**Platform.** The installer and the portable R are built for 64-bit Windows 10 and 11 only. macOS, Linux and 32-bit Windows are not supported: running from source may work on macOS and Linux but has not been tested.

## Run from source

You need [R](https://www.r-project.org/) 4.3 or newer.

### The Qt window (recommended)

```
cd sangerly_qt
pip install -r requirements.txt
set SANGERLY_ROOT=D:\path\to\SangeRly          :: Windows, the folder with app.R (the repository folder)
export SANGERLY_ROOT=/path/to/SangeRly         # macOS / Linux
python run.py
```

On the first start R loads its packages, which takes longer. If the log says "Rscript was not found", open **Settings** in the window and choose `Rscript` and the folder with `app.R`.

### The Shiny app only (in the browser)

```r
setwd("path/to/SangeRly")        # the repository folder, with app.R
shiny::runApp()
```

### R packages

```r
install.packages(c("shiny", "bslib", "ggplot2", "DT", "httr", "xml2", "jsonlite",
                   "readxl", "writexl", "processx", "BiocManager"))
BiocManager::install(c("sangerseqR", "Biostrings", "pwalign", "BiocGenerics"))
BiocManager::install(c("DECIPHER", "msa"))                  # MSA (at least one)
install.packages(c("ape", "phangorn", "pegas", "maps", "mapdata"))   # Phylogeny, Population structure (optional)
```

### External programs (optional)

| Program | Used for |
|---|---|
| [IQ-TREE](http://www.iqtree.org/) | maximum-likelihood trees, ModelFinder |
| [trimAl](https://github.com/inab/trimal) | trimming alignments |
| [minimap2](https://github.com/lh3/minimap2) | reference-guided alignment (without it, a built-in R aligner is used) |

Put them in the `app` folder next to the program (for example `app\IQTREE\bin\iqtree3.exe`), or type the path in the tool. The Windows installer already contains them when they were present at build time.

## Build the Windows installer

Run once on a Windows PC with internet where the program already works.

1. **Portable R.** In R, from the folder with `app.R`: `source("desktop/build_portable.R")`. It ends with `OK: the portable R loads every package`. The result is `dist/Sangerly`.
2. **The exe.** From `sangerly_qt`, run `packaging\build_all.bat`, or:
   ```
   python packaging/build_exe.py --portable "D:/path/SangeRly/dist/Sangerly" --app "D:/path/SangeRly" --tools "D:/path/SangeRly"
   ```
   The result is `sangerly_qt/dist_qt/SangeRly/` (test `SangeRly.exe` there). IQ-TREE, trimAl and minimap2 are copied from the folders next to `app.R` when present.
3. **The installer.** Open `sangerly_qt/packaging/installer_qt.iss` in [Inno Setup](https://jrsoftware.org/isinfo.php) 6.3+ and compile. The result is `sangerly_qt/Output/SangeRly_Setup.exe`. To sign it, see [`sangerly_qt/packaging/SIGNING.md`](sangerly_qt/packaging/SIGNING.md).

`SangeRly.exe` alone does not run: it needs the folders next to it, so share the installer (or a zip of the whole `SangeRly` folder).

## Repository layout

```
app.R, R/, www/    the R code (Shiny app)
Data_example/      example data for every tool
desktop/           builds the portable R
sangerly_qt/       the Qt window: run.py, sangerly/, engine/, packaging/ (and CHANGES.md)
LICENSE            GNU GPL v3
CITATION.cff       how to cite
THIRD_PARTY_NOTICES.md, licenses/   licenses of the programs shipped in the installer
CONTRIBUTING.md, CODE_OF_CONDUCT.md, .github/   issue and pull-request templates
```

Example data for every tool is in `Data_example/`.

## macOS, Linux and 32-bit Windows

Not supported. The Qt code looks for R in the usual macOS and Linux locations, so running from source may work, but this has not been tested. The portable build and the installer are for 64-bit Windows only. 32-bit Windows cannot be supported: R 4.2 and later and PySide6 (Qt) have no 32-bit Windows builds.

## Citation

If you use SangeRly in your work, please cite it (GitHub also shows a "Cite this repository" button, made from `CITATION.cff`):

> Alponsin, Putra, A. E. (2026). *Sangerly: R Shiny application for Sanger sequence assembly, BLAST identification, alignment and phylogeny* (Version 1.0.0) [Software]. https://doi.org/10.5281/zenodo.24126657

The **About and citation** page in the program lists the references of the methods it uses (and of IQ-TREE, trimAl and minimap2, which you should cite when you use their results).

## License

SangeRly is free software, licensed under the **GNU General Public License v3.0 or later** (see [`LICENSE`](LICENSE)).

The Windows installer also contains other programs (R, IQ-TREE, trimAl, minimap2, Qt and others) that keep their own licenses. They are listed in [`THIRD_PARTY_NOTICES.md`](THIRD_PARTY_NOTICES.md), and the license texts of IQ-TREE, trimAl and minimap2 are in [`licenses/`](licenses/).

## Contributing

Bug reports, ideas and pull requests are welcome: see [`CONTRIBUTING.md`](CONTRIBUTING.md) and the [Code of Conduct](CODE_OF_CONDUCT.md). Please report bugs and questions in [Issues](https://github.com/alponsin/sangerly/issues).
