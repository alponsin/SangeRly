# Privacy

SangeRly runs on your own computer. It has no account, no analytics, no telemetry and no server of its own. Your reads, alignments, trees and projects stay in the folders you choose.

The program contacts other services **only when you press the button that needs them**, and sends only what that function needs:

| Function | Service | What is sent |
|---|---|---|
| BLAST identification | NCBI BLAST (`blast.ncbi.nlm.nih.gov`) | the query sequence and the BLAST settings you chose |
| Find sequences, GenBank records and annotation, primer-BLAST | NCBI (`eutils.ncbi.nlm.nih.gov`, `www.ncbi.nlm.nih.gov`) | your search terms or accession numbers; for primer-BLAST, the primers and the template |
| Barcode (BOLD) search | BOLD (`portal.boldsystems.org`) | your search terms or the query sequence |
| Find primers | PR2 Primers (`app.pr2-primers.org`), qPrimerDB (`qprimerdb.biodb.org`), PrimerBank (`pga.mgh.harvard.edu`) | your search terms |
| Links such as iTOL or ChiPlot | your browser opens the website | nothing is sent automatically; you upload the file yourself |

These services have their own privacy policies and terms of use. **Do not send unpublished or confidential sequences to a public service** unless you accept that.

Inside the program, the R process and the Qt window talk to each other only through `127.0.0.1` (your own computer); nothing listens on the network.

The program does not check for updates and does not send crash reports. Log files are written to your computer only.
