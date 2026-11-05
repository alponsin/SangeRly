# ==============================================================================
# Sangerly - Help page (installation, use, problems). Texts follow the EN | ID switch.
# ==============================================================================

helpUI <- function(id) {
  t <- function(en, idn) shiny::tags$span(class = "i18n", `data-en` = en, `data-id` = idn, en)
  sec <- function(title_en, title_id, ..., open = FALSE)
    shiny::tags$details(class = "sg-help", open = if (open) NA else NULL,
      shiny::tags$summary(t(title_en, title_id)), shiny::div(class = "sg-help-body", ...))
  p <- function(en, idn) shiny::p(t(en, idn))
  li <- function(en, idn) shiny::tags$li(t(en, idn))
  code <- function(x) shiny::tags$code(x)

  shiny::div(class = "container", style = "max-width:900px",
    shiny::tags$style(shiny::HTML("
      details.sg-help{background:#fff;border:1px solid #d9e2ec;border-radius:10px;margin:10px 0;padding:0}
      details.sg-help>summary{cursor:pointer;padding:12px 16px;font-weight:700;color:#1f3c5f;list-style:none}
      details.sg-help>summary::before{content:'\\25B8';display:inline-block;width:18px;color:#14a3b8}
      details.sg-help[open]>summary::before{content:'\\25BE'}
      .sg-help-body{padding:2px 18px 14px;color:#334e68}
      .sg-help-body ol,.sg-help-body ul{padding-left:20px} .sg-help-body li{margin:4px 0}
      .sg-help-body code{background:#f0f4f8;padding:1px 5px;border-radius:4px}")),
    shiny::div(class = "sg-hero", shiny::h1(shiny::icon("circle-question"), " Help"),
               shiny::p(t("Installation, how to use Sangerly, and what to do when something does not work.",
                          "Instalasi, cara memakai Sangerly, dan apa yang dilakukan bila ada yang tidak berjalan."))),

    sec("Installing (Windows)", "Instalasi (Windows)", open = TRUE,
      shiny::tags$ol(
        li("Run the installer (Sangerly-Setup.exe). Windows 10 or 11, 64-bit. Administrator rights are not needed if you install in your own folder, for example D:\\Biotool\\sangerly.",
           "Jalankan installer (Sangerly-Setup.exe). Windows 10 atau 11, 64-bit. Hak administrator tidak diperlukan bila dipasang di folder sendiri, misalnya D:\\Biotool\\sangerly."),
        li("If Windows shows 'Windows protected your PC' (SmartScreen), click 'More info' and then 'Run anyway'. This appears because the installer is not digitally signed.",
           "Jika Windows menampilkan 'Windows protected your PC' (SmartScreen), klik 'More info' lalu 'Run anyway'. Ini muncul karena installer belum ditandatangani secara digital."),
        li("Install on a local drive (SSD if possible), not on a USB stick or a OneDrive folder.",
           "Pasang di drive lokal (sebaiknya SSD), bukan di flashdisk atau folder OneDrive."),
        li("Add the Sangerly folder to the antivirus exclusions (see 'Antivirus' below). This makes the start much faster.",
           "Tambahkan folder Sangerly ke pengecualian antivirus (lihat 'Antivirus' di bawah). Ini membuat program jauh lebih cepat dibuka."),
        li("Start Sangerly from the Start menu or the desktop icon. The first start after a restart can take up to a minute.",
           "Buka Sangerly dari menu Start atau ikon desktop. Buka pertama setelah komputer dihidupkan bisa sampai satu menit."))),

    sec("Quick start", "Mulai cepat",
      shiny::tags$ol(
        li("Assembly: load the Forward and Reverse .ab1 files of a sample, check the chromatograms and edit the consensus.",
           "Assembly: muat berkas .ab1 Forward dan Reverse, periksa kromatogram, dan edit konsensus."),
        li("BLAST: search the consensus in NCBI (internet needed) and send the hits you want to the alignment with 'Send to MSA'.",
           "BLAST: cari konsensus di NCBI (butuh internet) dan kirim hasil yang diinginkan ke alignment dengan 'Send to MSA'."),
        li("MSA and trimming: align, check the alignment in the editor, then trim the poorly aligned columns.",
           "MSA dan trimming: lakukan alignment, periksa di editor, lalu pangkas kolom yang sejajarnya buruk."),
        li("Phylogeny: build a neighbour-joining tree, and a maximum-likelihood tree with IQ-TREE.",
           "Filogeni: bangun pohon neighbour-joining dan pohon maximum-likelihood dengan IQ-TREE."),
        li("Visualisation: export Newick for iTOL, ChiPlot or tvBOT, or draw the figure here.",
           "Visualisasi: ekspor Newick untuk iTOL, ChiPlot atau tvBOT, atau gambar figurnya di sini."),
        li("More > Project and report: save the whole project (.rds) and download the report with the Methods paragraph. Use 'Load example data' to try the program.",
           "More > Project and report: simpan seluruh proyek (.rds) dan unduh laporan beserta paragraf Methods. Pakai 'Load example data' untuk mencoba program."))),

    sec("Antivirus and Windows Defender", "Antivirus dan Windows Defender",
      p("Sangerly contains thousands of small files (R, packages, IQ-TREE, trimAl). Antivirus programs scan each of them, which makes the start slow. Some antivirus programs also report a false alarm on programs such as IQ-TREE or trimAl, because they are unsigned command-line tools.",
        "Sangerly berisi ribuan berkas kecil (R, paket, IQ-TREE, trimAl). Antivirus memindai satu per satu sehingga program lambat dibuka. Sebagian antivirus juga salah mengira program seperti IQ-TREE atau trimAl berbahaya, karena keduanya program baris perintah tanpa tanda tangan digital."),
      shiny::tags$ol(
        li("Open Windows Security > Virus & threat protection > Manage settings > Exclusions > Add or remove exclusions.",
           "Buka Windows Security > Virus & threat protection > Manage settings > Exclusions > Add or remove exclusions."),
        li("Click 'Add an exclusion' > Folder and choose the Sangerly folder (for example D:\\Biotool\\sangerly). Add the folder %LOCALAPPDATA%\\Sangerly as well.",
           "Klik 'Add an exclusion' > Folder dan pilih folder Sangerly (misalnya D:\\Biotool\\sangerly). Tambahkan juga folder %LOCALAPPDATA%\\Sangerly."),
        li("Or, in PowerShell opened as Administrator: Add-MpPreference -ExclusionPath \"D:\\Biotool\\sangerly\".",
           "Atau, di PowerShell yang dibuka sebagai Administrator: Add-MpPreference -ExclusionPath \"D:\\Biotool\\sangerly\"."),
        li("If a file was removed or quarantined: Windows Security > Protection history > choose the item > Actions > Restore (or 'Allow on device'), then reinstall if needed.",
           "Jika ada berkas yang dihapus atau dikarantina: Windows Security > Protection history > pilih item > Actions > Restore (atau 'Allow on device'), lalu pasang ulang bila perlu."),
        li("Other antivirus (Kaspersky, Avast, ESET, and so on): add the same folders in that program.",
           "Antivirus lain (Kaspersky, Avast, ESET, dan lain-lain): tambahkan folder yang sama di program tersebut."),
        li("If the exclusion button is grey or you get 'not enough permissions', the computer is managed (campus or office). Ask the IT department to allow the folder.",
           "Jika tombol pengecualian berwarna abu-abu atau muncul 'not enough permissions', komputer dikelola (kampus atau kantor). Minta bagian IT mengizinkan folder tersebut.")),
      p("An exclusion means the folder is not scanned: keep only Sangerly files in it.",
        "Pengecualian berarti folder tidak dipindai: simpan hanya berkas Sangerly di dalamnya.")),

    sec("Sangerly is slow to open or the window stays empty", "Sangerly lama dibuka atau jendelanya kosong",
      shiny::tags$ul(
        li("Click the icon only once and wait up to a minute. Every extra click opens another empty window on the same program.",
           "Klik ikon satu kali saja dan tunggu sampai satu menit. Setiap klik tambahan membuka jendela kosong baru pada program yang sama."),
        li("Add the antivirus exclusion (above). It is the most common cause of a slow start.",
           "Tambahkan pengecualian antivirus (di atas). Ini penyebab paling umum dari buka yang lambat."),
        li("The second start is faster: Sangerly keeps running for 10 minutes after you close the window.",
           "Buka kedua lebih cepat: Sangerly tetap berjalan 10 menit setelah jendela ditutup."),
        li("Empty window: press Ctrl+Shift+R in it. If it stays empty, close all windows, end the processes Rscript.exe and Rterm.exe in Task Manager, delete %LOCALAPPDATA%\\Sangerly\\instance.txt and start again.",
           "Jendela kosong: tekan Ctrl+Shift+R di jendela itu. Bila tetap kosong, tutup semua jendela, akhiri proses Rscript.exe dan Rterm.exe di Task Manager, hapus %LOCALAPPDATA%\\Sangerly\\instance.txt lalu buka lagi."),
        li("Messages and times of every start are written to %LOCALAPPDATA%\\Sangerly\\session.log; a fatal error is written to last_error.txt in the same folder. Send these files when you ask for help.",
           "Pesan dan waktu setiap pembukaan ditulis di %LOCALAPPDATA%\\Sangerly\\session.log; galat fatal ditulis di last_error.txt pada folder yang sama. Kirim berkas ini bila meminta bantuan."))),

    sec("Other problems", "Masalah lain",
      shiny::tags$ul(
        li("'IQ-TREE or trimAl not found': they are inside the Sangerly folder (iqtree and trimAL). If an antivirus removed them, restore them (see above) or type the path of the program in the card.",
           "'IQ-TREE atau trimAl tidak ditemukan': keduanya ada di folder Sangerly (iqtree dan trimAL). Jika dihapus antivirus, pulihkan (lihat di atas) atau ketik lokasi programnya di kartu yang bersangkutan."),
        li("BLAST does not answer: it needs the internet and NCBI can be busy; wait and try again, or lower the number of hits.",
           "BLAST tidak menjawab: butuh internet dan NCBI kadang sibuk; tunggu lalu coba lagi, atau kurangi jumlah hit."),
        li("A project file (.rds) does not open: open only files made by Sangerly, from the same or a newer version. Raw chromatograms are not stored in the project.",
           "Berkas proyek (.rds) tidak terbuka: buka hanya berkas buatan Sangerly, dari versi yang sama atau lebih baru. Kromatogram mentah tidak disimpan di proyek."),
        li("The program closes by itself: read %LOCALAPPDATA%\\Sangerly\\session.log.",
           "Program menutup sendiri: baca %LOCALAPPDATA%\\Sangerly\\session.log."))),

    sec("Getting help and reporting a problem", "Bantuan dan melaporkan masalah",
      shiny::p(t("Report a problem or a suggestion at ", "Laporkan masalah atau saran di "),
               shiny::tags$a(href = paste0(SANGERLY_URL, "/issues"), target = "_blank", paste0(SANGERLY_URL, "/issues")),
               t(". Include what you did, the message you saw, the version (About page) and session.log.",
                 ". Sertakan apa yang Anda lakukan, pesan yang muncul, versi (halaman About) dan session.log."))))
}
