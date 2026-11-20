# Code signing policy

Free code signing for the Windows releases of SangeRly is provided by [SignPath.io](https://about.signpath.io/), certificate by the [SignPath Foundation](https://signpath.org/).

- **What is signed:** `SangeRly_Setup.exe` of each release (built from the tagged source in this repository by the GitHub Actions workflow [`release.yml`](.github/workflows/release.yml); nothing is built or signed on a personal computer).
- **Roles:** Committers and reviewers: the repository owner and collaborators listed under *Contributors*. Approvers of a signing request: the repository owner (**Alponsin**). Every release is approved by hand in SignPath before it is signed.
- **Rules:** only code from this repository is signed; the bundled third-party programs (R, IQ-TREE, trimAl, minimap2, Qt and others) are unmodified binaries of their own projects, listed in [`THIRD_PARTY_NOTICES.md`](THIRD_PARTY_NOTICES.md). The program does not contain malware, adware or code that changes system settings beyond its own install folder.
- **Privacy:** see [`PRIVACY.md`](PRIVACY.md).
- **Verify:** right-click the installer, *Properties*, *Digital Signatures*; the SHA-256 checksum is on the release page.
