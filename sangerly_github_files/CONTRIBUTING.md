# Contributing to SangeRly

Thank you for your interest. Bug reports, questions, ideas and pull requests are welcome.

## Report a bug or ask a question
Open an [issue](https://github.com/alponsin/sangerly/issues) and tell us:
- the SangeRly version (shown in the status bar of the window, or the About page),
- your operating system,
- what you did, what you expected and what happened,
- if possible, a small example file (an `.ab1` or FASTA) and the text of the Log panel.

## Suggest a change
Open an issue first for larger changes, so we can agree on the idea before you spend time on it.

## Pull requests
1. Fork the repository and create a branch (`git checkout -b my-change`).
2. Keep a change small and focused. Do not mix formatting changes with new features.
3. Run the program and try the tool you changed with the files in `sangerly_r/Data_example/`.
4. Describe in `CHANGES.md` what changed and why.
5. Open the pull request and explain what you tested.

## Do not commit
Build results and programs: `dist/`, `dist_qt/`, `build_qt/`, `Output/`, portable R, and the `.exe` files of IQ-TREE, trimAl and minimap2 (the `.gitignore` already skips most of them). Installers go to [Releases](https://github.com/alponsin/sangerly/releases).

## Licence
By contributing you agree that your contribution is licensed under the GNU General Public License v3 or later, like the rest of SangeRly (see `LICENSE`).

Please follow the [Code of Conduct](CODE_OF_CONDUCT.md).
