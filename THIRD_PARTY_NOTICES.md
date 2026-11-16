# SangeRly: licence and third-party software

SangeRly (c) 2026 Alponsin and Andani Eka Putra, licensed under the GNU General Public License,
version 3 or (at your option) any later version. See the file `LICENSE`.

SangeRly runs the following programs as separate executables. They are NOT part of the SangeRly
source code and keep their own licences. They are included only in the Windows desktop installer.
The full licence text of each one is in the folder [`licenses/`](licenses/).

| Program  | Licence      | Licence text | Source |
|----------|--------------|--------------|--------|
| trimAl   | GPL-3.0      | [licenses/trimAl-GPL-3.0.txt](licenses/trimAl-GPL-3.0.txt) | https://github.com/inab/trimal |
| IQ-TREE  | GPL-2.0      | [licenses/IQ-TREE-GPL-2.0.txt](licenses/IQ-TREE-GPL-2.0.txt) | https://github.com/iqtree/iqtree3 |
| minimap2 | MIT          | [licenses/minimap2-MIT.txt](licenses/minimap2-MIT.txt) | https://github.com/lh3/minimap2 |

Please cite them when you use their results (see the About page in the app).
SangeRly also uses R packages (shiny, ape, Biostrings, msa, DECIPHER and others), each under its own licence.

## Also included in the Windows installer

| Component | Licence | Source |
|-----------|---------|--------|
| R (portable copy) | GPL-2 or GPL-3 | https://www.r-project.org/ and https://cran.r-project.org/ |
| R packages (copied from the build computer) | each package's own licence, see the `DESCRIPTION` file of the package in `R/library/` | https://cran.r-project.org/ , https://bioconductor.org/ |
| Qt and PySide6 (the window, with QtWebEngine / Chromium) | LGPL-3.0 and other open-source licences of Qt and Chromium | https://www.qt.io/ , https://doc.qt.io/qtforpython/ |
| PyInstaller (builds `SangeRly.exe`) | GPL-2.0 with a bootloader exception that allows any program | https://pyinstaller.org/ |

The source code of the unmodified versions of these programs is available from the addresses above.
If you cannot get it from there, open an issue at https://github.com/alponsin/sangerly/issues.
