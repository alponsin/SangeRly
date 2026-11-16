# ==============================================================================
# SangeRly - Tutorial page (More > Tutorial)
#   * "What do you want to do?" cards -> jump to the right lessons
#   * one lesson per step (and per extra tool): numbered steps, tips, "Open this step" button
#   * progress bar ("x of y lessons read", remembered in the browser), search box, expand / collapse all
#   * a "has results" badge on the lessons whose step already has results (the sg_done message of app.R)
# All texts follow the EN | ID switch. No server part is needed: the buttons use the same "go_tab"
# input as the cards of the Home page.
# ==============================================================================

tutorialUI <- function(id) {
  tg  <- shiny::tags
  t   <- function(en, idn) tg$span(class = "i18n", `data-en` = en, `data-id` = idn, en)
  img <- function(p) if (exists("img_src", mode = "function")) img_src(p) else p

  li   <- function(en, idn) tg$li(t(en, idn))
  steps_ol <- function(...) tg$ol(class = "tt-steps", ...)
  tip  <- function(en, idn, kind = "tip")
    shiny::div(class = paste("tt-note", kind),
               shiny::icon(if (kind == "warn") "triangle-exclamation" else "lightbulb"),
               shiny::div(t(en, idn)))
  need <- function(en, idn)
    shiny::div(class = "tt-need", shiny::icon("box-open"), tg$b(t("You need: ", "Yang diperlukan: ")), t(en, idn))
  goal <- function(en, idn)
    shiny::div(class = "tt-goal", shiny::icon("flag-checkered"), tg$b(t("Result: ", "Hasil: ")), t(en, idn))

  # one lesson (a collapsible block)
  lesson <- function(key, num, title_en, title_id, mins, ..., go = NULL, go_en = NULL, go_id = NULL,
                     logo = NULL, icon = NULL, open = FALSE, count = TRUE) {
    tg$details(class = "tt-lesson", id = paste0("tt_", key), `data-lesson` = key,
               `data-go` = go, `data-count` = if (count) "1" else NULL, open = if (open) NA else NULL,
      tg$summary(
        shiny::span(class = "tt-num", if (!is.null(num)) num else shiny::icon("star")),
        if (!is.null(logo))
          tg$img(class = "tt-logo", src = img(logo), alt = "", onerror = "this.style.display='none'")
        else if (!is.null(icon)) shiny::span(class = "tt-ico", shiny::icon(icon)),
        shiny::span(class = "tt-title", t(title_en, title_id)),
        shiny::span(class = "tt-meta",
                    shiny::span(class = "tt-badge tt-b-res", shiny::icon("circle-check"), t("has results", "ada hasil")),
                    shiny::span(class = "tt-badge tt-b-read", shiny::icon("check"), t("read", "sudah dibaca")),
                    shiny::span(class = "tt-time", shiny::icon("clock"), mins))),
      shiny::div(class = "tt-body", ...,
        shiny::div(class = "tt-actions",
          if (!is.null(go))
            tg$button(type = "button", class = "btn btn-primary btn-sm tt-go", `data-go` = go,
                      t(go_en, go_id), " ", shiny::icon("arrow-right")),
          tg$button(type = "button", class = "btn btn-outline-secondary btn-sm tt-markread", `data-lesson` = key,
                    shiny::span(class = "tt-off", t("Mark as read and continue", "Tandai sudah dibaca dan lanjut")),
                    shiny::span(class = "tt-on",  t("Read (click to undo)", "Sudah dibaca (klik untuk batal)")))))
    )
  }

  # a "what do you want to do" card: title + a path of lessons
  goal_card <- function(icon, en, idn, path) {
    shiny::div(class = "tt-gcard",
      shiny::div(class = "tt-gico", shiny::icon(icon)),
      tg$h5(t(en, idn)),
      shiny::div(class = "tt-path",
        lapply(seq_along(path), function(i) {
          p <- path[[i]]
          shiny::tagList(
            if (i > 1) shiny::span(class = "tt-arrow", shiny::icon("chevron-right")),
            tg$button(type = "button", class = "tt-chip tt-jump", `data-lesson` = p[[1]], t(p[[2]], p[[3]])))
        })))
  }

  shiny::div(class = "container", style = "max-width:980px",
    tg$style(shiny::HTML("
      /* hero + progress */
      .tt-wrap .sg-hero{padding-bottom:4px}
      .tt-prog{background:#fff;border:1px solid #d9e2ec;border-radius:12px;padding:12px 16px;margin:14px 0 6px}
      .tt-prog-top{display:flex;justify-content:space-between;align-items:center;gap:10px;flex-wrap:wrap;color:#1f3c5f;font-weight:600}
      .tt-bar{height:10px;background:#e4ebf2;border-radius:6px;overflow:hidden;margin-top:8px}
      .tt-bar>div{height:100%;width:0;background:linear-gradient(90deg,#0b3d91,#14a3b8,#22b07d);transition:width .35s ease}
      .tt-prog .btn-link{font-size:12px;padding:0;color:#6b7785}
      /* toolbar */
      .tt-tools{display:flex;gap:8px;align-items:center;flex-wrap:wrap;margin:14px 0 4px}
      .tt-tools input{flex:1 1 240px;border:1px solid #b8c2cc;border-radius:8px;padding:7px 12px;min-width:0}
      .tt-tools input:focus{outline:none;border-color:#14a3b8;box-shadow:0 0 0 3px rgba(20,163,184,.18)}
      .tt-h{color:#1f3c5f;font-weight:700;margin:26px 0 10px;font-size:19px}
      /* goal cards */
      .tt-goals{display:grid;grid-template-columns:repeat(auto-fit,minmax(290px,1fr));gap:14px}
      .tt-gcard{background:#fff;border:1px solid #d9e2ec;border-radius:12px;padding:14px 16px;box-shadow:0 1px 3px rgba(0,0,0,.07)}
      .tt-gcard h5{margin:6px 0 10px;color:#1f3c5f;font-weight:700;font-size:15.5px}
      .tt-gico{width:36px;height:36px;border-radius:10px;background:linear-gradient(135deg,#0b3d91,#14a3b8);color:#fff;
               display:flex;align-items:center;justify-content:center;font-size:17px}
      .tt-path{display:flex;flex-wrap:wrap;align-items:center;gap:4px}
      .tt-arrow{color:#9fb3c8;font-size:11px}
      .tt-chip{border:1px solid #14a3b8;background:#f0fbfd;color:#0b5f70;border-radius:999px;padding:3px 11px;font-size:12.5px;font-weight:600;cursor:pointer}
      .tt-chip:hover{background:#14a3b8;color:#fff}
      /* example-data box */
      .tt-try{display:flex;gap:14px;align-items:center;flex-wrap:wrap;background:#eef7ff;border:1px solid #b9d9f7;border-radius:12px;padding:12px 16px;margin-top:16px;color:#1f3c5f}
      .tt-try .ico{font-size:24px;color:#2c7be5}
      .tt-try .txt{flex:1 1 280px}
      /* lessons */
      details.tt-lesson{background:#fff;border:1px solid #d9e2ec;border-radius:12px;margin:10px 0;box-shadow:0 1px 3px rgba(0,0,0,.05);
                        transition:border-color .15s ease, box-shadow .15s ease}
      details.tt-lesson[open]{border-color:#14a3b8;box-shadow:0 4px 14px rgba(20,163,184,.15)}
      details.tt-lesson>summary{cursor:pointer;list-style:none;display:flex;align-items:center;gap:12px;padding:11px 16px;color:#1f3c5f}
      details.tt-lesson>summary::-webkit-details-marker{display:none}
      details.tt-lesson>summary::after{content:'\\25BE';margin-left:auto;color:#14a3b8;transition:transform .15s ease;order:9}
      details.tt-lesson[open]>summary::after{transform:rotate(180deg)}
      .tt-num{flex:0 0 30px;width:30px;height:30px;border-radius:50%;background:linear-gradient(135deg,#0b3d91,#14a3b8);color:#fff;
              display:flex;align-items:center;justify-content:center;font-weight:700;font-size:14px}
      .tt-logo{height:38px;width:60px;object-fit:contain;flex:0 0 60px}
      .tt-ico{flex:0 0 34px;text-align:center;font-size:22px;color:#14a3b8}
      .tt-title{font-weight:700;font-size:16.5px;flex:1 1 auto}
      .tt-meta{display:flex;gap:8px;align-items:center;flex-wrap:wrap;justify-content:flex-end}
      .tt-time{color:#6b7785;font-size:12.5px;white-space:nowrap}
      .tt-time .fa,.tt-time .fas{margin-right:4px}
      .tt-badge{display:none;font-size:11.5px;font-weight:700;border-radius:999px;padding:2px 9px;white-space:nowrap}
      .tt-badge .fa,.tt-badge .fas{margin-right:4px}
      .tt-b-res{background:#e6f6ef;color:#15805a}
      .tt-b-read{background:#e8eefb;color:#0b3d91}
      details.tt-lesson.tt-hasres .tt-b-res{display:inline-block}
      details.tt-lesson.tt-read .tt-b-read{display:inline-block}
      details.tt-lesson.tt-read .tt-num{background:#22b07d}
      .tt-body{padding:2px 20px 16px 58px;color:#334e68}
      @media (max-width:640px){.tt-body{padding-left:18px}.tt-logo{display:none}}
      .tt-need,.tt-goal{background:#f6f8fb;border-radius:8px;padding:7px 12px;margin:8px 0;font-size:14px}
      .tt-need .fa,.tt-need .fas,.tt-goal .fa,.tt-goal .fas{color:#14a3b8;margin-right:6px}
      ol.tt-steps{list-style:none;counter-reset:tt;padding-left:0;margin:10px 0}
      ol.tt-steps>li{counter-increment:tt;position:relative;padding:5px 0 5px 36px;margin:2px 0;min-height:30px}
      ol.tt-steps>li::before{content:counter(tt);position:absolute;left:0;top:4px;width:24px;height:24px;border-radius:50%;
                             background:#e1f3f7;color:#0b5f70;font-weight:700;font-size:12.5px;display:flex;align-items:center;justify-content:center}
      .tt-note{display:flex;gap:10px;align-items:flex-start;border-radius:8px;padding:8px 12px;margin:8px 0;font-size:13.5px}
      .tt-note .fa,.tt-note .fas{margin-top:3px}
      .tt-note.tip{background:#fff8e1;border-left:4px solid #f5b301}
      .tt-note.tip .fa,.tt-note.tip .fas{color:#c98f00}
      .tt-note.warn{background:#fdecea;border-left:4px solid #e0554a}
      .tt-note.warn .fa,.tt-note.warn .fas{color:#c0392b}
      .tt-actions{display:flex;gap:8px;flex-wrap:wrap;margin-top:14px;padding-top:12px;border-top:1px dashed #d9e2ec}
      .tt-markread .tt-on{display:none}
      .tt-markread.on .tt-on{display:inline}
      .tt-markread.on .tt-off{display:none}
      .tt-markread.on{background:#22b07d;border-color:#22b07d;color:#fff}
      .tt-empty{display:none;text-align:center;color:#6b7785;padding:18px}
      .tt-faq summary{font-weight:600}
      .tt-hide{display:none !important}
    ")),

    tg$script(shiny::HTML("
      (function(){
        var KEY = 'sg_tut_read';
        function getRead(){ try { return JSON.parse(localStorage.getItem(KEY) || '[]'); } catch(e) { return []; } }
        function setRead(a){ try { localStorage.setItem(KEY, JSON.stringify(a)); } catch(e) {} }
        function lessons(){ return Array.prototype.slice.call(document.querySelectorAll('details.tt-lesson')); }
        function refresh(){
          var r = getRead(), total = 0, done = 0;
          lessons().forEach(function(l){
            var k = l.getAttribute('data-lesson'), on = r.indexOf(k) >= 0;
            l.classList.toggle('tt-read', on);
            var b = l.querySelector('.tt-markread'); if (b) b.classList.toggle('on', on);
            if (l.hasAttribute('data-count')) { total++; if (on) done++; }
          });
          var n = document.getElementById('tt_n'), t = document.getElementById('tt_t'), f = document.getElementById('tt_fill');
          if (n) n.textContent = done; if (t) t.textContent = total;
          if (f) f.style.width = (total ? (100 * done / total) : 0) + '%';
          var ok = document.getElementById('tt_alldone'); if (ok) ok.style.display = (total && done === total) ? 'inline' : 'none';
        }
        function openLesson(key, scroll){
          var el = document.getElementById('tt_' + key); if (!el) return;
          el.open = true;
          if (scroll !== false) setTimeout(function(){ el.scrollIntoView({behavior: 'smooth', block: 'start'}); }, 40);
        }
        $(document).on('click', '.tt-go', function(e){
          e.preventDefault();
          if (window.Shiny && Shiny.setInputValue) Shiny.setInputValue('go_tab', this.getAttribute('data-go'), {priority: 'event'});
          window.scrollTo({top: 0});
        });
        $(document).on('click', '.tt-jump', function(){ openLesson(this.getAttribute('data-lesson')); });
        $(document).on('click', '.tt-markread', function(){
          var k = this.getAttribute('data-lesson'), r = getRead(), i = r.indexOf(k), wasOn = i >= 0;
          if (wasOn) r.splice(i, 1); else r.push(k);
          setRead(r); refresh();
          if (!wasOn) {                                     // go on to the next lesson
            var ls = lessons(), pos = ls.findIndex(function(l){ return l.getAttribute('data-lesson') === k; });
            ls[pos].open = false;
            var nx = null;
            for (var j = pos + 1; j < ls.length; j++) { if (!ls[j].classList.contains('tt-hide')) { nx = ls[j]; break; } }
            if (nx) { nx.open = true; setTimeout(function(){ nx.scrollIntoView({behavior: 'smooth', block: 'start'}); }, 60); }
          }
        });
        $(document).on('click', '#tt_expand',   function(){ lessons().forEach(function(l){ if (!l.classList.contains('tt-hide')) l.open = true; }); });
        $(document).on('click', '#tt_collapse', function(){ lessons().forEach(function(l){ l.open = false; }); });
        $(document).on('click', '#tt_reset',    function(){ setRead([]); refresh(); });
        // search: looks in both languages
        var cache = null;
        function build(){ cache = {}; lessons().forEach(function(l){
          var s = l.textContent;
          l.querySelectorAll('.i18n').forEach(function(e){ s += ' ' + (e.getAttribute('data-en') || '') + ' ' + (e.getAttribute('data-id') || ''); });
          cache[l.getAttribute('data-lesson')] = s.toLowerCase(); }); }
        $(document).on('input', '#tt_search', function(){
          if (!cache) build();
          var q = this.value.trim().toLowerCase(), shown = 0;
          lessons().forEach(function(l){
            var hit = !q || cache[l.getAttribute('data-lesson')].indexOf(q) >= 0;
            l.classList.toggle('tt-hide', !hit); if (hit) shown++;
            if (q && hit) l.open = true;
          });
          var em = document.getElementById('tt_empty'); if (em) em.style.display = shown ? 'none' : 'block';
        });
        $(function(){ refresh(); });
      })();
    ")),

    shiny::div(class = "tt-wrap",
      shiny::div(class = "sg-hero",
        shiny::h1(shiny::icon("graduation-cap"), " ", t("Tutorial", "Tutorial")),
        shiny::p(t("New to SangeRly? Pick what you want to do, then follow the short lessons. Each one tells you what you need, what to click and what you get.",
                   "Baru memakai SangeRly? Pilih apa yang ingin Anda kerjakan, lalu ikuti pelajaran singkatnya. Setiap pelajaran menjelaskan apa yang diperlukan, apa yang diklik, dan hasilnya."))),

      shiny::div(class = "tt-prog",
        shiny::div(class = "tt-prog-top",
          shiny::span(shiny::icon("list-check"), " ", tg$b(id = "tt_n", "0"), " / ", tg$span(id = "tt_t", "0"), " ",
                      t("lessons read", "pelajaran dibaca"),
                      shiny::span(id = "tt_alldone", style = "display:none;color:#15805a", " - ", t("well done!", "kerja bagus!"))),
          tg$button(type = "button", id = "tt_reset", class = "btn btn-link", t("Reset progress", "Atur ulang progres"))),
        shiny::div(class = "tt-bar", shiny::div(id = "tt_fill"))),

      # ---- what do you want to do? ------------------------------------------------------------------------
      tg$h4(class = "tt-h", t("What do you want to do?", "Apa yang ingin Anda kerjakan?")),
      shiny::div(class = "tt-goals",
        goal_card("dna", "Identify a species from my .ab1 files", "Mengidentifikasi spesies dari berkas .ab1 saya",
                  list(list("assembly", "1. Assembly", "1. Perakitan"), list("blast", "2. BLAST", "2. BLAST"))),
        goal_card("diagram-project", "Build a phylogenetic tree", "Membuat pohon filogenetik",
                  list(list("assembly", "Assembly", "Perakitan"), list("msa", "MSA", "MSA"),
                       list("tree", "Phylogeny", "Filogeni"), list("vis", "Visualisation", "Visualisasi"))),
        goal_card("circle-nodes", "Haplotypes, diversity and networks", "Haplotipe, keragaman, dan jaringan",
                  list(list("msa", "MSA", "MSA"), list("popgen", "1. Data and run", "1. Data dan jalankan"), list("popgen_res", "2. Results", "2. Hasil"),
                       list("popgen_net", "3. Network and map", "3. Jaringan dan peta"), list("popgen_more", "4. More tests", "4. More tests"), list("popgen_save", "5. Save and compare", "5. Simpan dan bandingkan"))),
        goal_card("ruler-horizontal", "Barcode gap and species limits", "Barcode gap dan batas spesies",
                  list(list("msa", "MSA", "MSA"), list("barcode", "Genetic distances", "Jarak genetik"))),
        goal_card("scissors", "Design primers", "Merancang primer",
                  list(list("primer", "Primer design", "Desain primer"))),
        goal_card("file-export", "Submit sequences to GenBank", "Mengirim sekuens ke GenBank",
                  list(list("msa", "MSA", "MSA"), list("genbank", "GenBank export", "Ekspor GenBank"))),
        goal_card("eye", "Just look at a file", "Hanya ingin melihat sebuah berkas",
                  list(list("open", "Open and view files", "Buka dan lihat berkas")))),

      shiny::div(class = "tt-try",
        shiny::span(class = "ico", shiny::icon("flask")),
        shiny::div(class = "txt", tg$b(t("No data yet? Practise with an example.", "Belum punya data? Berlatih dengan data contoh.")), tg$br(),
                   t("Use 'Load the synthetic demo' in Project and report, or 'Load example data' in Population structure. Nothing you do there changes your own files.",
                     "Pakai 'Load the synthetic demo' di Proyek dan laporan, atau 'Load example data' di Struktur populasi. Apa pun yang Anda lakukan di sana tidak mengubah berkas Anda.")),
        tg$button(type = "button", class = "btn btn-primary btn-sm tt-go", `data-go` = "project", t("Project and report", "Proyek dan laporan")),
        tg$button(type = "button", class = "btn btn-outline-primary btn-sm tt-go", `data-go` = "popgen", t("Population structure", "Struktur populasi"))),

      # ---- lessons ----------------------------------------------------------------------------------------
      tg$h4(class = "tt-h", t("Lessons", "Pelajaran")),
      shiny::div(class = "tt-tools",
        tg$input(type = "search", id = "tt_search", placeholder = "Search the tutorial... (e.g. BLAST, trimming, FST, primer)", autocomplete = "off"),
        tg$button(type = "button", id = "tt_expand", class = "btn btn-outline-secondary btn-sm", t("Expand all", "Buka semua")),
        tg$button(type = "button", id = "tt_collapse", class = "btn btn-outline-secondary btn-sm", t("Collapse all", "Tutup semua"))),

      lesson("start", NULL, "Start here: how SangeRly works", "Mulai di sini: cara kerja SangeRly", "2 min",
             count = FALSE, open = TRUE, icon = "compass",
        steps_ol(
          li("Choose a step on the Home page, or click its logo in the top menu. Every step is a tab and you can open any of them at any time.",
             "Pilih satu tahap di halaman Home, atau klik logonya di menu atas. Setiap tahap adalah satu tab dan bisa dibuka kapan saja."),
          li("The tabs pass data to each other. The consensus of Assembly goes to BLAST and MSA; the alignment of MSA goes to Phylogeny, Population structure and Genetic distances. Buttons such as 'Send to MSA' move sequences for you.",
             "Antar tab saling meneruskan data. Konsensus dari Perakitan dipakai BLAST dan MSA; alignment dari MSA dipakai Filogeni, Struktur populasi, dan Jarak genetik. Tombol seperti 'Send to MSA' memindahkan sekuens untuk Anda."),
          li("A green tick on a card or on a logo means that this step already has results.",
             "Tanda centang hijau pada kartu atau logo berarti tahap itu sudah punya hasil."),
          li("Use the EN | ID switch at the top right to change the language at any time.",
             "Pakai tombol EN | ID di kanan atas untuk mengganti bahasa kapan saja."),
          li("Save your work before you close the program: More > Project and report > 'Save project (.rds)'. Next time use 'Open a project (.rds)'.",
             "Simpan pekerjaan sebelum menutup program: More > Proyek dan laporan > 'Save project (.rds)'. Lain kali pakai 'Open a project (.rds)'.")),
        tip("Hover over a logo in the top menu to see the name of the step.",
            "Arahkan kursor ke logo di menu atas untuk melihat nama tahapnya.")),

      lesson("assembly", "1", "Assembly: from .ab1 reads to a consensus", "Perakitan: dari read .ab1 menjadi konsensus", "10 min",
             go = "assembly", go_en = "Open Assembly", go_id = "Buka Perakitan", logo = "img/logo_assembly.png",
        need("The .ab1 files of your samples (Forward and Reverse). A reference FASTA is optional.",
             "Berkas .ab1 sampel Anda (Forward dan Reverse). FASTA referensi bersifat opsional."),
        steps_ol(
          li("Add your files with 'Add .ab1 files', 'Browse files...' or 'Add folder'. You can add files many times, also from different folders.",
             "Tambahkan berkas dengan 'Add .ab1 files', 'Browse files...' atau 'Add folder'. Boleh berkali-kali, juga dari folder yang berbeda."),
          li("Check 'How files are grouped into samples'. Forward and Reverse files of one sample must end up in the same sample (for example S1_F and S1_R).",
             "Periksa 'How files are grouped into samples'. Berkas Forward dan Reverse dari satu sampel harus berada di sampel yang sama (misalnya S1_F dan S1_R)."),
          li("Optional: load 'Reference region(s)' (FASTA) when you have many reads or tiled regions, and choose how the reads are placed on it.",
             "Opsional: muat 'Reference region(s)' (FASTA) bila Anda punya banyak read atau wilayah yang bertumpuk, lalu pilih cara read ditempatkan."),
          li("Choose the 'Trimming method' (it removes low-Phred ends) and press 'ASSEMBLE / RECALCULATE'.",
             "Pilih 'Trimming method' (membuang ujung dengan Phred rendah) lalu tekan 'ASSEMBLE / RECALCULATE'."),
          li("Go through the samples. The buttons '< Needs review' and 'Needs review >' jump to the ones where the reads disagree. Look at the chromatogram: if a base is wrong, delete it with 'Delete this read's base here' or 'Delete whole position (all reads + consensus)'. 'Restore everything' undoes all deletions.",
             "Telusuri sampelnya. Tombol '< Needs review' dan 'Needs review >' melompat ke sampel yang read-nya tidak cocok. Lihat kromatogramnya: bila ada basa yang salah, hapus dengan 'Delete this read's base here' atau 'Delete whole position (all reads + consensus)'. 'Restore everything' membatalkan semua penghapusan."),
          li("Download 'Consensus FASTA' (one sample) or 'All consensus sequences (FASTA)'. 'Save session (.rds)' keeps everything to continue later.",
             "Unduh 'Consensus FASTA' (satu sampel) atau 'All consensus sequences (FASTA)'. 'Save session (.rds)' menyimpan semuanya untuk dilanjutkan nanti.")),
        tg$h6(class = "mt-3", t("Reading the Alignment tab", "Membaca tab Alignment")),
        shiny::div(class = "mb-2",
          lapply(list(
            list("background:#f7e08a", "needs review", "perlu diperiksa"),
            list("background:#bfe6cf", "decided manually", "diputuskan manual"),
            list("background:#9aa5b1", "position deleted", "posisi dihapus"),
            list("border:2px solid #c0392b", "read differs from consensus", "read berbeda dari konsensus"),
            list("background:#eee", "x = base deleted", "x = basa dihapus"),
            list("background:#eee", "base transparency = Phred", "transparansi basa = Phred")),
            function(ch) tg$span(class = "legend-chip", style = ch[[1]], t(ch[[2]], ch[[3]])))),
        tg$ul(class = "tt-steps",
          li("'Jump to base number' goes to a column (de novo mode) or to a reference base number (reference mode). Extra columns for insertions or read ends that run past the reference have no number.",
             "'Jump to base number' menuju sebuah kolom (mode de novo) atau nomor basa referensi (mode referensi). Kolom tambahan untuk insersi atau ujung read yang melewati referensi tidak bernomor."),
          li("Drag the alignment to move along it. Click a base to select it; click the Consensus row to select the whole position. Shift+click or Shift+arrows extend the selection. Ctrl+wheel zooms.",
             "Seret alignment untuk bergerak. Klik sebuah basa untuk memilihnya; klik baris Consensus untuk memilih seluruh posisi. Shift+klik atau Shift+panah memperluas pilihan. Ctrl+roda mouse untuk zoom."),
          tg$li(t("Delete / Backspace removes the selection: on the Consensus row the whole position (all reads and the consensus), on a read row only that read's bases. Ctrl+Z undoes.",
                  "Delete / Backspace menghapus pilihan: pada baris Consensus seluruh posisi (semua read dan konsensus), pada baris read hanya basa read itu. Ctrl+Z membatalkan.")),
          tg$li(tg$b(t("Edit the consensus from the keyboard: ", "Edit konsensus dari keyboard: ")),
                t("select a position (or a range) and type A, C, G, T or N (also IUPAC letters, or - ) to set the consensus base; = (or . / Enter) uses the reference base; 0 (or Esc) returns to the rule decision. The selection moves on to the next position, so you can keep typing. Typing always edits the consensus, whichever row is selected.",
                  "pilih sebuah posisi (atau rentang) lalu ketik A, C, G, T atau N (juga huruf IUPAC, atau - ) untuk mengatur basa konsensus; = (atau . / Enter) memakai basa referensi; 0 (atau Esc) kembali ke keputusan aturan. Pilihan berpindah ke posisi berikutnya, jadi Anda bisa terus mengetik. Mengetik selalu mengedit konsensus, baris mana pun yang dipilih.")),
          li("Click the viewer first so that it receives the keys.",
             "Klik viewer terlebih dahulu agar menerima tombol keyboard.")),
        goal("A clean consensus sequence for each sample, ready for BLAST and MSA.",
             "Sekuens konsensus yang bersih untuk setiap sampel, siap untuk BLAST dan MSA."),
        tip("Always look at the chromatogram where Forward and Reverse disagree. A clear single peak is more reliable than the base call.",
            "Selalu lihat kromatogram di posisi yang Forward dan Reverse-nya berbeda. Satu puncak yang jelas lebih dapat dipercaya daripada hasil base call.")),

      lesson("blast", "2", "BLAST: identify the species", "BLAST: identifikasi spesies", "5 min",
             go = "blast", go_en = "Open BLAST", go_id = "Buka BLAST", logo = "img/logo_blast.png",
        need("A consensus from Assembly (or a sequence you paste). Internet is needed to search NCBI.",
             "Konsensus dari Perakitan (atau sekuens yang Anda tempel). Butuh internet untuk mencari di NCBI."),
        steps_ol(
          li("Type your e-mail at the top of the tab (recommended by NCBI; it is not stored).",
             "Ketik e-mail Anda di bagian atas tab (disarankan NCBI; tidak disimpan)."),
          li("Press 'Use consensus from the assembly', or paste a sequence.",
             "Tekan 'Use consensus from the assembly', atau tempel sebuah sekuens."),
          li("Choose the 'Database' and the 'Program'. To search only inside a group of organisms, add taxon names and press 'Validate taxon names'.",
             "Pilih 'Database' dan 'Program'. Untuk mencari hanya di satu kelompok organisme, tambahkan nama takson lalu tekan 'Validate taxon names'."),
          li("Choose 'Where to run BLAST': in this app (results appear here) or on the NCBI website. Press 'RUN BLAST' and wait; it can take a few minutes. 'STOP' cancels it and 'Load RID' reopens an earlier search.",
             "Pilih 'Where to run BLAST': di aplikasi ini (hasil muncul di sini) atau di situs NCBI. Tekan 'RUN BLAST' dan tunggu; bisa beberapa menit. 'STOP' membatalkan dan 'Load RID' membuka kembali pencarian sebelumnya."),
          li("Read the hit table: the best hits have a high identity and a high coverage. 'Add lineage from NCBI Taxonomy' adds the taxonomy of each hit.",
             "Baca tabel hit: hit terbaik memiliki identity dan coverage yang tinggi. 'Add lineage from NCBI Taxonomy' menambahkan taksonomi tiap hit."),
          li("Download the 'Hit table (CSV)' or the FASTA files, or press 'Send to MSA' to take the hits to the alignment.",
             "Unduh 'Hit table (CSV)' atau berkas FASTA, atau tekan 'Send to MSA' untuk membawa hit ke alignment.")),
        goal("The most likely species of each sample, and reference sequences for your tree.",
             "Spesies yang paling mungkin untuk tiap sampel, dan sekuens referensi untuk pohon Anda."),
        tip("If BLAST does not answer, NCBI is probably busy. Wait a little and try again, or ask for fewer hits.",
            "Bila BLAST tidak menjawab, kemungkinan NCBI sedang sibuk. Tunggu sebentar lalu coba lagi, atau minta jumlah hit yang lebih sedikit.")),

      lesson("msa", "3", "MSA and trimming: align your sequences", "MSA dan trimming: menyejajarkan sekuens", "8 min",
             go = "msa", go_en = "Open MSA", go_id = "Buka MSA", logo = "img/logo_msa.png",
        need("Two or more sequences: consensus sequences, BLAST hits, FASTA files or pasted text.",
             "Dua sekuens atau lebih: konsensus, hit BLAST, berkas FASTA, atau teks yang ditempel."),
        steps_ol(
          li("Step 1, in the window choose Tools > MSA and trimming > Multiple sequence alignment > Align with ...: DECIPHER, Clustal Omega, ClustalW, MUSCLE, 'Align to a reference' (like ViralMSA), or 'I already have an alignment' to go straight to trimming. A window opens, like in UGENE.",
             "Langkah 1, di jendela pilih Tools > MSA and trimming > Multiple sequence alignment > Align with ...: DECIPHER, Clustal Omega, ClustalW, MUSCLE, 'Align to a reference' (seperti ViralMSA), atau 'I already have an alignment' untuk langsung ke trimming. Sebuah jendela terbuka, seperti di UGENE."),
          li("Step 2, in that window give the input: the sequences already in the list (for example hits sent from BLAST), FASTA files or a folder, the consensus of the Assembly, or pasted FASTA. 'Output file' (optional) also saves the alignment as FASTA. Press 'Align'; the alignment opens in the Editor tab.",
             "Langkah 2, di jendela itu isi inputnya: sekuens yang sudah ada di daftar (misalnya hit yang dikirim dari BLAST), berkas FASTA atau folder, konsensus dari Perakitan, atau FASTA yang ditempel. 'Output file' (opsional) juga menyimpan alignment sebagai FASTA. Tekan 'Align'; alignment terbuka di tab Editor."),
          li("Check the alignment in the editor. Ragged ends or columns full of gaps are normal for Sanger reads.",
             "Periksa alignment di editor. Ujung yang tidak rata atau kolom penuh gap adalah hal biasa pada read Sanger."),
          li("Step 3, trimming: choose Tools > MSA and trimming > Trim the alignment. In the window choose the 'Trimming engine' and the 'Rule', then press 'Trim' to remove the poorly aligned columns. The result shows in the 'Trimmed result' tab.",
             "Langkah 3, trimming: pilih Tools > MSA and trimming > Trim the alignment. Di jendelanya pilih 'Trimming engine' dan 'Rule', lalu tekan 'Trim' untuk membuang kolom yang sejajarnya buruk. Hasilnya tampil di tab 'Trimmed result'."),
          li("Download 'Alignment (FASTA)' or 'Trimmed alignment (FASTA)'. The other tabs use the trimmed alignment automatically.",
             "Unduh 'Alignment (FASTA)' atau 'Trimmed alignment (FASTA)'. Tab lain memakai alignment yang sudah dipangkas secara otomatis.")),
        goal("A trimmed alignment that feeds Phylogeny, Population structure and Genetic distances.",
             "Alignment terpangkas yang dipakai oleh Filogeni, Struktur populasi, dan Jarak genetik."),
        tip("The tab can also list the mutations of every sequence against a reference ('Compare every sequence to') and save them as Excel, HTML or CSV.",
            "Tab ini juga dapat mendaftar mutasi setiap sekuens terhadap sebuah referensi ('Compare every sequence to') dan menyimpannya sebagai Excel, HTML, atau CSV.")),

      lesson("tree", "4", "Phylogeny: build a tree", "Filogeni: membangun pohon", "6 min",
             go = "tree", go_en = "Open Phylogeny", go_id = "Buka Filogeni", logo = "img/logo_phylo.png",
        need("An alignment from the MSA tab (the trimmed one is best).",
             "Alignment dari tab MSA (yang sudah dipangkas paling baik)."),
        steps_ol(
          li("In the window choose Tools > Phylogeny > Build a tree. A window opens, like 'Build Phylogenetic Tree' in UGENE: pick the 'Tree building method' on top, then the tabs 'Model', 'Bootstrapping and Consensus tree' and 'Display Options'.",
             "Di jendela pilih Tools > Phylogeny > Build a tree. Sebuah jendela terbuka, seperti 'Build Phylogenetic Tree' di UGENE: pilih 'Tree building method' di atas, lalu tab 'Model', 'Bootstrapping and Consensus tree' dan 'Display Options'."),
          li("Quick look first: Neighbor-Joining with a 'Distance matrix model'. It takes seconds.",
             "Lihat cepat dulu: Neighbor-Joining dengan 'Distance matrix model'. Hanya butuh beberapa detik."),
          li("For a maximum-likelihood tree choose Maximum Likelihood (IQ-TREE), the 'Substitution model' (ModelFinder can pick the best one) and the 'Branch support', then press 'Build'.",
             "Untuk pohon maximum-likelihood pilih Maximum Likelihood (IQ-TREE), 'Substitution model' (ModelFinder dapat memilih yang terbaik) dan 'Branch support', lalu tekan 'Build'."),
          li("Use 'Root the drawn tree' to root on your outgroup, and 'Layout' to change the shape.",
             "Pakai 'Root the drawn tree' untuk mengakarkan pohon pada outgroup Anda, dan 'Layout' untuk mengubah bentuknya."),
          li("Save the tree files with 'Save to this folder' (use 'Browse folder...' to choose where).",
             "Simpan berkas pohon dengan 'Save to this folder' (pakai 'Browse folder...' untuk memilih lokasi).")),
        goal("An NJ tree and/or an IQ-TREE tree with branch support.",
             "Pohon NJ dan/atau pohon IQ-TREE dengan dukungan cabang."),
        tip("Use NJ for a first look and ML (IQ-TREE) for the tree that goes in a paper. IQ-TREE is a separate program that comes with the installer.",
            "Pakai NJ untuk gambaran awal dan ML (IQ-TREE) untuk pohon yang masuk ke naskah. IQ-TREE adalah program terpisah yang ikut dalam installer.")),

      lesson("vis", "5", "Visualisation: make the tree figure", "Visualisasi: membuat figur pohon", "8 min",
             go = "vis", go_en = "Open Visualisation", go_id = "Buka Visualisasi", logo = "img/logo_visual.png",
        need("A tree from the Phylogeny tab (or a Newick file). A table with information about your sequences is optional.",
             "Pohon dari tab Filogeni (atau berkas Newick). Tabel berisi informasi sekuens Anda bersifat opsional."),
        steps_ol(
          li("Choose the 'Tree to show' (NJ or IQ-TREE), or open a tree file with 'Open tree file...'.",
             "Pilih 'Tree to show' (NJ atau IQ-TREE), atau buka berkas pohon dengan 'Open tree file...'."),
          li("To colour the tree by species, locality and so on, press 'Open table...' and load a table. Then choose 'Column with the sequence names' and 'Names shown in the tree: from column'.",
             "Untuk mewarnai pohon berdasarkan spesies, lokasi, dan sebagainya, tekan 'Open table...' dan muat sebuah tabel. Lalu pilih 'Column with the sequence names' dan 'Names shown in the tree: from column'."),
          li("Add colour with 'Symbols by', 'Colour blocks by' and 'Name colour by'. Choose the 'Layout', the 'Root' and the 'Show the support as' style.",
             "Tambahkan warna dengan 'Symbols by', 'Colour blocks by', dan 'Name colour by'. Pilih 'Layout', 'Root', dan gaya 'Show the support as'."),
          li("Drag the tree and the legends with the mouse to place them. 'Put the tree back in its place' and 'Reset the legend positions' undo this.",
             "Seret pohon dan legenda dengan mouse untuk menempatkannya. 'Put the tree back in its place' dan 'Reset the legend positions' membatalkannya."),
          li("Choose the 'PNG resolution' and download the figure. 'NJ tree (Newick)' and 'IQ-TREE tree (Newick)' give files for iTOL, ChiPlot or tvBOT.",
             "Pilih 'PNG resolution' lalu unduh figurnya. 'NJ tree (Newick)' dan 'IQ-TREE tree (Newick)' memberi berkas untuk iTOL, ChiPlot, atau tvBOT.")),
        goal("A publication figure of your tree, or Newick files for online tree viewers.",
             "Figur pohon siap terbit, atau berkas Newick untuk penampil pohon daring."),
        tip("Names in the table must be the same as the names in the tree. If nothing is coloured, check the 'Column with the sequence names'.",
            "Nama pada tabel harus sama dengan nama pada pohon. Jika tidak ada yang berwarna, periksa 'Column with the sequence names'.", "warn")),

      lesson("popgen", "6", "Population structure 1: data, populations, settings, run", "Struktur populasi 1: data, populasi, pengaturan, menjalankan", "15 min",
             go = "popgen", go_en = "Open Population structure", go_id = "Buka Struktur populasi", logo = "img/logo_popup.png",
        need("An alignment of ONE gene or region (from the MSA tab or an aligned FASTA file; all sequences the same length, at least 3) and the population of each sequence.",
             "Alignment SATU gen atau wilayah (dari tab MSA atau berkas FASTA ter-align; semua sekuens sama panjang, minimal 3) dan populasi dari setiap sekuens."),
        steps_ol(
          li("Section 1 - Alignment. Choose 'MSA tab: kept columns (after trimming / editing)', 'MSA tab: whole alignment' or 'An aligned FASTA file' (then 'Open aligned FASTA...'). The line under the choice tells how many sequences and columns there are and how many sites have no gap, N or IUPAC code. If more than half of the columns are lost, trim the ragged ends in the MSA tab first.",
             "Bagian 1 - Alignment. Pilih 'MSA tab: kept columns (after trimming / editing)', 'MSA tab: whole alignment' atau 'An aligned FASTA file' (lalu 'Open aligned FASTA...'). Baris di bawah pilihan menunjukkan jumlah sekuens dan kolom serta berapa situs yang bebas gap, N, atau kode IUPAC. Bila lebih dari separuh kolom hilang, pangkas dulu ujung yang tidak rata di tab MSA."),
          li("No data at hand? Press 'Load example data': 60 ITS sequences of Durio zibethinus with their locations. Everything below can be tried on it. Ready-made example files are in the 'contoh' folder (contoh_alignment_Durio_ITS.fasta, contoh_populasi_Durio.csv / .xlsx, contoh_populasi_hanya_nomor_akses.csv, contoh_koordinat_Durio.csv).",
             "Belum punya data? Tekan 'Load example data': 60 sekuens ITS Durio zibethinus beserta lokasinya. Semua langkah di bawah bisa dicoba dengan data ini. Berkas contoh ada di folder 'contoh' (contoh_alignment_Durio_ITS.fasta, contoh_populasi_Durio.csv / .xlsx, contoh_populasi_hanya_nomor_akses.csv, contoh_koordinat_Durio.csv)."),
          li("Section 2 - Populations. Choose ONE source; the label under the choice always shows which one is active. (a) 'From the sequence names': a rule on the name, for example the first field before '_' (Siberut_01, Siberut_02 -> population 'Siberut'); you can set the separator, the field number, the number of characters or a regular expression.",
             "Bagian 2 - Populasi. Pilih SATU sumber; label di bawah pilihan selalu menunjukkan sumber yang aktif. (a) 'From the sequence names': aturan pada nama, misalnya field pertama sebelum '_' (Siberut_01, Siberut_02 -> populasi 'Siberut'); pemisah, nomor field, jumlah karakter, atau regular expression dapat diatur."),
          li("(b) 'From a CSV / Excel file': one row per sequence, with the columns Sequence (the full name as in the alignment, or only the accession such as MT230650.1), Population and, optionally, Group (region or island, used by the two-level AMOVA). The first row holds the column names. Press 'Table with my sequence names (CSV)' to download a file already filled with your names, add the Population column in Excel, then open it with 'Open CSV / Excel...'. If the program picks the wrong column, change it in 'Sequence name in column', 'Population in column' and 'Group in column'.",
             "(b) 'From a CSV / Excel file': satu baris per sekuens, dengan kolom Sequence (nama lengkap seperti di alignment, atau hanya nomor akses seperti MT230650.1), Population dan, opsional, Group (wilayah atau pulau, dipakai AMOVA dua tingkat). Baris pertama berisi nama kolom. Tekan 'Table with my sequence names (CSV)' untuk mengunduh berkas yang sudah berisi nama sekuens Anda, tambahkan kolom Population di Excel, lalu buka dengan 'Open CSV / Excel...'. Bila program salah memilih kolom, ubah di 'Sequence name in column', 'Population in column', dan 'Group in column'."),
          li("(c) 'By hand': select sequences in the list, type the population (or group) name and press 'Update population' / 'Update group'. You can also double-click a cell in the table of populations to rename a population or to give it a group. The table under the list shows how many sequences each population has. A population of 1 or 2 sequences is kept but gives few statistics.",
             "(c) 'By hand': pilih sekuens di daftar, ketik nama populasi (atau grup), lalu tekan 'Update population' / 'Update group'. Anda juga dapat mengklik dua kali sel pada tabel populasi untuk mengganti nama populasi atau memberinya grup. Tabel di bawah daftar menunjukkan jumlah sekuens tiap populasi. Populasi dengan 1 atau 2 sekuens tetap dipakai tetapi memberi sedikit statistik."),
          li("Section 3 - Settings. The defaults suit a first look. The important ones: 'Sites used for diversity' (complete in all sequences = default; complete in each population = DnaSP default; all sites with pairwise deletion = DnaSP 'pairwise'), 'Divide pi by the variable sites' (only to copy Arlequin when its file came from DnaSP with the invariable sites removed), 'Simulated samples for the p-value of Fu's Fs' and 'Ewens-Watterson test: simulated samples' (1000 to explore, 10 000 or more for the final p-values), 'Genome' (diploid, haploid, X or Y-linked for the DnaSP gene-flow estimates), the permutations of the differentiation test, the 'Exact test of differentiation' steps (use 1 000 000 for a final value), the AMOVA distance method, gamma, permutations and significance level.",
             "Bagian 3 - Pengaturan. Nilai awal cukup untuk gambaran pertama. Yang penting: 'Sites used for diversity' (complete in all sequences = default; complete in each population = default DnaSP; all sites with pairwise deletion = 'pairwise' DnaSP), 'Divide pi by the variable sites' (hanya untuk meniru Arlequin bila berkasnya berasal dari DnaSP dengan situs invarian dibuang), 'Simulated samples for the p-value of Fu's Fs' dan 'Ewens-Watterson test: simulated samples' (1000 untuk eksplorasi, 10 000 atau lebih untuk nilai p akhir), 'Genome' (diploid, haploid, X atau Y untuk estimasi aliran gen DnaSP), jumlah permutasi uji diferensiasi, langkah 'Exact test of differentiation' (pakai 1 000 000 untuk nilai akhir), metode jarak AMOVA, gamma, permutasi, dan taraf nyata."),
          li("Have an Arlequin settings file? Press 'Open Arlequin settings (.ars)...' and the form is filled with the same values, so the analysis uses the same settings as in Arlequin.",
             "Punya berkas pengaturan Arlequin? Tekan 'Open Arlequin settings (.ars)...' dan formulir terisi dengan nilai yang sama, sehingga analisis memakai pengaturan seperti di Arlequin."),
          li("Section 4 - Run. 'Haplotypes, diversity, differentiation' first (it feeds the Result, Haplotype, Network, Map, Sliding window and More tests tabs). Then 'AMOVA and pairwise FST' (needs at least 2 populations; a Group column gives the hierarchical version) and, if you want, 'Structure (K clusters)'. The coloured boxes at the top show which steps are done, not run yet or outdated.",
             "Bagian 4 - Jalankan. Mulai dengan 'Haplotypes, diversity, differentiation' (hasilnya dipakai tab Result, Haplotype, Network, Map, Sliding window, dan More tests). Lalu 'AMOVA and pairwise FST' (minimal 2 populasi; kolom Group memberi versi hierarkis) dan, bila perlu, 'Structure (K clusters)'. Kotak berwarna di bagian atas menunjukkan tahap yang selesai, belum dijalankan, atau kedaluwarsa.")),
        goal("Haplotype and diversity tables, neutrality tests, FST / AMOVA, ready for the tabs of the next lessons.",
             "Tabel haplotipe dan keragaman, uji netralitas, FST / AMOVA, siap untuk tab pada pelajaran berikutnya."),
        tip("Use one gene per analysis. Change the populations or the alignment after a run and the results stay on the screen but are marked 'outdated': run them again.",
            "Pakai satu gen per analisis. Bila populasi atau alignment diubah setelah dijalankan, hasil tetap tampil tetapi ditandai 'outdated': jalankan lagi."),
        tip("Populations with fewer than 4 sequences are left out of Tajima's D, Fu & Li and R2, and a population with 1 sequence has no diversity to measure. FST with populations of 1-3 sequences is not reliable.",
            "Populasi dengan kurang dari 4 sekuens tidak dimasukkan dalam Tajima's D, Fu & Li, dan R2, dan populasi dengan 1 sekuens tidak punya keragaman yang dapat diukur. FST dengan populasi 1-3 sekuens tidak dapat diandalkan.", "warn")),

      lesson("popgen_res", "6b", "Population structure 2: read the results", "Struktur populasi 2: membaca hasil", "15 min",
             go = "popgen", go_en = "Open Population structure", go_id = "Buka Struktur populasi", icon = "table",
        need("The run of 'Haplotypes, diversity, differentiation' (and, for the last tabs, 'AMOVA and pairwise FST' and 'Structure').",
             "Hasil 'Haplotypes, diversity, differentiation' (dan, untuk tab terakhir, 'AMOVA and pairwise FST' serta 'Structure')."),
        steps_ol(
          li("Result tab, diversity table (one row for all sequences and one per population): n, sites, S (segregating sites), Eta, h (haplotypes), Hd (haplotype diversity, with SD), k (mean pairwise differences), Pi (per site), the four thetas (Hom, S, k, Pi), Tajima's D, Fu's Fs and the p-values. Hover over a column title to read what it means. Hd and Pi high = much diversity; both low = small, young or bottlenecked population.",
             "Tab Result, tabel keragaman (satu baris untuk semua sekuens dan satu per populasi): n, situs, S (situs segregasi), Eta, h (haplotipe), Hd (keragaman haplotipe, dengan SD), k (rerata perbedaan berpasangan), Pi (per situs), empat theta (Hom, S, k, Pi), Tajima's D, Fs Fu, dan nilai p. Arahkan kursor ke judul kolom untuk membaca artinya. Hd dan Pi tinggi = keragaman besar; keduanya rendah = populasi kecil, muda, atau mengalami bottleneck."),
          li("Neutrality tests: Tajima's D, Fu's Fs, Fu & Li's D*, F* and Y*, Ramos-Onsins & Rozas' R2, Ewens-Watterson, Strobeck's S and Chakraborty's test. Negative D, Fs, D*, F* and a small R2 (with a small p) mean an excess of rare variants: population expansion, a selective sweep or purifying selection. Positive values mean a deficit of rare variants: population structure, a recent bottleneck or balancing selection. One locus cannot separate demography from selection. Strobeck's S gives the probability of k or fewer haplotypes; Chakraborty's test gives the probability of k or more haplotypes given Theta(Hom), and a small value suggests a sample that mixes populations.",
             "Uji netralitas: Tajima's D, Fs Fu, D*, F*, dan Y* Fu & Li, R2 Ramos-Onsins & Rozas, Ewens-Watterson, S Strobeck, dan uji Chakraborty. D, Fs, D*, F* negatif dan R2 kecil (dengan p kecil) berarti kelebihan varian langka: ekspansi populasi, selective sweep, atau seleksi pemurni. Nilai positif berarti kekurangan varian langka: struktur populasi, bottleneck baru-baru ini, atau seleksi penyeimbang. Satu lokus tidak dapat memisahkan demografi dari seleksi. S Strobeck memberi peluang k haplotipe atau kurang; uji Chakraborty memberi peluang k haplotipe atau lebih pada Theta(Hom), dan nilai kecil menunjukkan sampel gabungan beberapa populasi."),
          li("Differentiation (DnaSP style): Hs, Hst, Ks, Kst, Fst, Nst, Snn, the chi-square of the haplotype table, gene flow (Nm) and the exact test of population differentiation (p-values by permutation or by the Markov chain). FST / PhiST: 0 = populations alike, 1 = completely different; p < 0.05 means the difference is unlikely by chance.",
             "Diferensiasi (gaya DnaSP): Hs, Hst, Ks, Kst, Fst, Nst, Snn, khi-kuadrat tabel haplotipe, aliran gen (Nm), dan exact test diferensiasi populasi (nilai p dengan permutasi atau rantai Markov). FST / PhiST: 0 = populasi serupa, 1 = benar-benar berbeda; p < 0,05 berarti perbedaan sulit terjadi secara kebetulan."),
          li("Haplotype tab: the table of haplotypes with the count per population and the names of the sequences in each. Mismatch distribution (Result tab): the histogram of pairwise differences with the expectation of a constant and of an expanding population.",
             "Tab Haplotype: tabel haplotipe dengan jumlah per populasi dan nama sekuens di tiap haplotipe. Mismatch distribution (tab Result): histogram perbedaan berpasangan beserta harapan populasi konstan dan populasi yang berekspansi."),
          li("AMOVA tab: variance among groups / among populations within groups / within populations, the fixation indices (PhiCT, PhiSC, PhiST) with p-values by permutation, and the pairwise FST (below the diagonal) with p-values (above); a + marks significance. Options: Slatkin's linearized FST, Reynolds' distance and the pi between populations.",
             "Tab AMOVA: ragam antar grup / antar populasi dalam grup / dalam populasi, indeks fiksasi (PhiCT, PhiSC, PhiST) dengan nilai p hasil permutasi, dan FST berpasangan (di bawah diagonal) dengan nilai p (di atas); tanda + berarti signifikan. Opsi: FST linier Slatkin, jarak Reynolds, dan pi antarpopulasi."),
          li("Sliding window tab: set the window (columns), the step and the minimum number of sites, press 'Run the sliding window' and choose what to draw (Pi, Watterson's theta, Tajima's D, S, Hd, h), also for each population. It shows where along the gene the variation is.",
             "Tab Sliding window: atur jendela (kolom), langkah, dan jumlah situs minimum, tekan 'Run the sliding window', lalu pilih yang digambar (Pi, theta Watterson, Tajima's D, S, Hd, h), juga untuk tiap populasi. Tab ini menunjukkan di mana variasi berada sepanjang gen."),
          li("Structure tab: choose the method (maximum-likelihood EM like ADMIXTURE / LEA, or the Gibbs sampler of STRUCTURE), the model (admixture or no admixture), K from-to, runs per K and press 'Structure (K clusters)'. K is chosen with the cross-entropy, BIC and Evanno's delta K. You can also load a Q table of STRUCTURE / ADMIXTURE (first column the sequence name, then one column per cluster) and the same figure is drawn.",
             "Tab Structure: pilih metode (maximum-likelihood EM seperti ADMIXTURE / LEA, atau Gibbs sampler STRUCTURE), model (admixture atau no admixture), K dari-sampai, jumlah run per K, lalu tekan 'Structure (K clusters)'. K dipilih dengan cross-entropy, BIC, dan delta K Evanno. Anda juga dapat memuat tabel Q dari STRUCTURE / ADMIXTURE (kolom pertama nama sekuens, lalu satu kolom per klaster) dan figur yang sama digambar.")),
        goal("You can say how diverse each population is, whether the data fit a neutral constant-size population and how much the populations differ.",
             "Anda dapat menyatakan seberapa beragam tiap populasi, apakah data sesuai populasi netral berukuran konstan, dan seberapa berbeda antarpopulasi."),
        tip("p-values from simulations or permutations change a little between runs. Use more simulations (10 000 or more) for the number you publish.",
            "Nilai p dari simulasi atau permutasi sedikit berubah antar run. Pakai lebih banyak simulasi (10 000 atau lebih) untuk angka yang dipublikasikan.")),

      lesson("popgen_net", "6c", "Population structure 3: haplotype network and map", "Struktur populasi 3: jaringan haplotipe dan peta", "15 min",
             go = "popgen", go_en = "Open Population structure", go_id = "Buka Struktur populasi", icon = "circle-nodes",
        need("The run of 'Haplotypes, diversity, differentiation'. For the map: a CSV or Excel file with the coordinates of the populations.",
             "Hasil 'Haplotypes, diversity, differentiation'. Untuk peta: berkas CSV atau Excel berisi koordinat populasi."),
        steps_ol(
          li("Network tab: choose the type under 'Network' - 'Median-joining (as NETWORK / DnaSP)' (it adds median vectors and leaves no dangling ones), 'Minimum spanning network' or 'Haplotype network (pegas::haploNet)' - and press 'Build network'. Circle size = number of sequences, colours = populations; the figure is drawn in your browser.",
             "Tab Network: pilih jenis di 'Network' - 'Median-joining (as NETWORK / DnaSP)' (menambahkan median vector tanpa yang menggantung), 'Minimum spanning network' atau 'Haplotype network (pegas::haploNet)' - lalu tekan 'Build network'. Ukuran lingkaran = jumlah sekuens, warna = populasi; figur digambar di browser Anda."),
          li("Mutations on the links: choose 'along the link, one number per position' so that every mutated position (e.g. 234, 256, 278) has its own place on the link, with a small tick; the other choices ('beside the link', 'on the link', 'near the first / second node') put all the numbers in one text. 'Max. positions shown' limits long lists (0 = all).",
             "Mutasi pada garis penghubung: pilih 'along the link, one number per position' agar setiap posisi bermutasi (mis. 234, 256, 278) punya tempat sendiri pada garis, dengan tanda kecil; pilihan lain ('beside the link', 'on the link', 'near the first / second node') menaruh semua angka dalam satu teks. 'Max. positions shown' membatasi daftar yang panjang (0 = semua)."),
          li("Edit the figure by mouse: drag a node, drag a label, drag a number of a link (each one moves alone), 'Undo' and 'Reset view' bring things back. Tick 'Size legend (circles)' for the nested circles that show how many sequences a circle stands for, and the legend of the populations.",
             "Edit figur dengan mouse: seret node, seret label, seret angka pada garis (masing-masing bergerak sendiri); 'Undo' dan 'Reset view' mengembalikan. Centang 'Size legend (circles)' untuk lingkaran bersarang yang menunjukkan jumlah sekuens satu lingkaran, dan legenda populasi."),
          li("Save the network as SVG (editable in Inkscape / Illustrator, labels stay objects) or PNG; 'Save SVG + PNG in the output folder' (and 'Save all results') writes both.",
             "Simpan jaringan sebagai SVG (dapat diedit di Inkscape / Illustrator, label tetap berupa objek) atau PNG; 'Save SVG + PNG in the output folder' (dan 'Save all results') menulis keduanya."),
          li("Map tab (as PopART): press 'Coordinates (CSV / Excel)...' and open a file with one row per population and the columns Population, Latitude and Longitude in decimal degrees (south and west negative). The names must be the same as the populations of section 2. The button 'Coordinates template (CSV)' gives a file with your populations to fill in; see contoh_koordinat_Durio.csv. The same coordinates are used by the Structure map and the Mantel test.",
             "Tab Map (seperti PopART): tekan 'Coordinates (CSV / Excel)...' dan buka berkas dengan satu baris per populasi dan kolom Population, Latitude, dan Longitude dalam derajat desimal (selatan dan barat negatif). Nama harus sama dengan populasi pada bagian 2. Tombol 'Coordinates template (CSV)' memberi berkas berisi populasi Anda untuk diisi; lihat contoh_koordinat_Durio.csv. Koordinat yang sama dipakai peta Structure dan uji Mantel."),
          li("Choose the 'Base map': satellite (Esri World Imagery), topographic, ocean, streets (OpenStreetMap), light grey, dark, terrain (OpenTopoMap), 'Coastlines (offline)' or a plain colour. Satellite and the other pictures are loaded from public servers, so they need internet; without internet the coastlines are shown and the pies are still drawn.",
             "Pilih 'Base map': satelit (Esri World Imagery), topografi, laut, jalan (OpenStreetMap), abu terang, gelap, terrain (OpenTopoMap), 'Coastlines (offline)' atau warna polos. Satelit dan peta gambar lain dimuat dari server publik sehingga butuh internet; tanpa internet garis pantai tampil dan pie tetap digambar."),
          li("On the map: drag to move, wheel or double-click to zoom, 'Fit' returns to all populations, drag a pie to a better place (e.g. when two overlap) and 'Reset pies' undoes it. Choose the size of the largest pie, the height, the names, the haplotype legend, the size legend, the scale bar, a title and your own haplotype colours. The buttons PNG and SVG above the map save what you see (the SVG has the base map as one picture; pies, names and legends are editable objects). Keep the credit of the base map written in the corner when you publish.",
             "Pada peta: seret untuk menggeser, roda mouse atau klik dua kali untuk zoom, 'Fit' kembali ke semua populasi, seret pie ke tempat yang lebih baik (mis. bila dua pie bertumpuk) dan 'Reset pies' membatalkannya. Atur ukuran pie terbesar, tinggi, nama, legenda haplotipe, legenda ukuran, skala, judul, dan warna haplotipe sendiri. Tombol PNG dan SVG di atas peta menyimpan yang Anda lihat (SVG memuat peta dasar sebagai satu gambar; pie, nama, dan legenda berupa objek yang dapat diedit). Pertahankan kredit peta dasar yang tertulis di sudut saat dipublikasikan.")),
        goal("A haplotype network and a haplotype map ready for a paper (SVG / PNG).",
             "Jaringan haplotipe dan peta haplotipe siap untuk naskah (SVG / PNG)."),
        tip("The network is built from the sites used in the run. If it looks too crowded, build it from a population subset or lower 'Max. positions shown'.",
            "Jaringan dibuat dari situs yang dipakai pada run. Bila terlalu ramai, buat dari subset populasi atau turunkan 'Max. positions shown'.")),

      lesson("popgen_more", "6d", "Population structure 4: More tests (LD, Mantel, HWE, phasing, outliers)", "Struktur populasi 4: More tests (LD, Mantel, HWE, fasa, outlier)", "15 min",
             go = "popgen", go_en = "Open Population structure", go_id = "Buka Struktur populasi", icon = "flask-vial",
        need("The run of 'Haplotypes, diversity, differentiation'. Mantel: also the AMOVA run and the coordinates (Map tab). Hardy-Weinberg and phasing: diploid sequences with IUPAC codes (R, Y, S, W, K, M) at the heterozygous sites.",
             "Hasil 'Haplotypes, diversity, differentiation'. Mantel: juga hasil AMOVA dan koordinat (tab Map). Hardy-Weinberg dan fasa: sekuens diploid dengan kode IUPAC (R, Y, S, W, K, M) pada situs heterozigot."),
        steps_ol(
          li("SFS and linkage (LD): choose the sequences (all or one population), the maximum number of sites and press 'Run SFS and LD'. The site frequency spectrum (folded) is compared with the neutral expectation (an excess of singletons = expansion / purifying selection). LD gives D, D', r2 and Fisher's exact test for every pair of biallelic sites, ZnS (mean r2), Za (adjacent sites), ZZ and Rm (minimum number of recombination events, four-gamete test), plus the plot of r2 against distance.",
             "SFS dan linkage (LD): pilih sekuens (semua atau satu populasi), jumlah situs maksimum, lalu tekan 'Run SFS and LD'. Spektrum frekuensi situs (terlipat) dibandingkan dengan harapan netral (kelebihan singleton = ekspansi / seleksi pemurni). LD memberi D, D', r2, dan uji eksak Fisher untuk setiap pasang situs biallelic, ZnS (rerata r2), Za (situs berdekatan), ZZ, dan Rm (jumlah minimum peristiwa rekombinasi, uji empat gamet), serta grafik r2 terhadap jarak."),
          li("Mantel test (isolation by distance): needs at least 4 populations with coordinates and pairwise FST (run 'AMOVA and pairwise FST' first). Choose the genetic distance (FST or FST / (1 - FST)), the geographic distance (km or ln km) and the permutations, press 'Run the Mantel test'. A positive r with a small one-sided p means that populations far apart are more different.",
             "Uji Mantel (isolation by distance): butuh minimal 4 populasi dengan koordinat dan FST berpasangan (jalankan 'AMOVA and pairwise FST' dulu). Pilih jarak genetik (FST atau FST / (1 - FST)), jarak geografis (km atau ln km), dan jumlah permutasi, lalu tekan 'Run the Mantel test'. r positif dengan p satu sisi kecil berarti populasi yang berjauhan lebih berbeda."),
          li("Hardy-Weinberg: for diploid individuals, every site with two bases and their heterozygous code (e.g. A, G and R) is tested with the exact test of Wigginton et al. (2005); the table gives the genotype counts, Ho, He, Fis and the p-value (also corrected by Bonferroni). Sequences without any IUPAC heterozygous code are haploid (mtDNA, cpDNA): there is nothing to test and the program says so.",
             "Hardy-Weinberg: untuk individu diploid, setiap situs dengan dua basa dan kode heterozigotnya (mis. A, G, dan R) diuji dengan exact test Wigginton dkk. (2005); tabel memberi jumlah genotipe, Ho, He, Fis, dan nilai p (juga terkoreksi Bonferroni). Sekuens tanpa kode heterozigot IUPAC bersifat haploid (mtDNA, cpDNA): tidak ada yang diuji dan program menyatakannya."),
          li("Phasing (EM): for diploid Sanger sequences with IUPAC codes, set the maximum number of heterozygous sites per individual (12 = 2048 possible pairs) and press 'Infer the haplotypes'. The EM algorithm of Excoffier and Slatkin (1995) returns the haplotype frequencies, the best pair of haplotypes of each individual with its posterior probability, and a FASTA with two phased sequences per individual (ID_a, ID_b). Open that FASTA as an alignment (sequence names ID_a, ID_b) to build networks and diversity with haplotypes instead of genotypes.",
             "Fasa (EM): untuk sekuens Sanger diploid dengan kode IUPAC, atur jumlah maksimum situs heterozigot per individu (12 = 2048 pasangan mungkin) lalu tekan 'Infer the haplotypes'. Algoritma EM Excoffier dan Slatkin (1995) menghasilkan frekuensi haplotipe, pasangan haplotipe terbaik tiap individu beserta probabilitas posterior, dan FASTA berisi dua sekuens terfasa per individu (ID_a, ID_b). Buka FASTA itu sebagai alignment (nama sekuens ID_a, ID_b) untuk membuat jaringan dan keragaman dengan haplotipe, bukan genotipe."),
          li("Fst outliers (FDIST idea): for biallelic sites (SNPs) the Fst and the heterozygosity of every site are compared with a neutral island model simulated by the coalescent with the same sample sizes; the migration rate is tuned so that the mean simulated Fst equals the observed one. Set the number of simulated loci and the minimum minor-allele count and press 'Detect outlier sites'. P_high small = more differentiated than neutral (divergent selection); P_low small = less differentiated (balancing selection). About 5% of the sites are expected to be 'significant' by chance.",
             "Outlier Fst (ide FDIST): untuk situs biallelic (SNP), Fst dan heterozigositas tiap situs dibandingkan dengan model pulau netral yang disimulasikan secara koalesen dengan ukuran sampel yang sama; laju migrasi diatur agar rerata Fst simulasi sama dengan yang teramati. Atur jumlah lokus simulasi dan jumlah minimum alel minor lalu tekan 'Detect outlier sites'. P_high kecil = lebih berbeda daripada netral (seleksi divergen); P_low kecil = kurang berbeda (seleksi penyeimbang). Sekitar 5% situs diharapkan 'signifikan' secara kebetulan.")),
        goal("LD and recombination, isolation by distance, genotype tests, phased haplotypes and candidate outlier sites.",
             "LD dan rekombinasi, isolation by distance, uji genotipe, haplotipe terfasa, dan situs outlier kandidat."),
        tip("Sites of one gene are linked, so the SNP-based tests (LD, Fst outliers) are exploratory for a single gene. They are meant for many independent SNPs.",
            "Situs dalam satu gen saling terkait, sehingga uji berbasis SNP (LD, outlier Fst) bersifat eksploratif untuk satu gen. Uji ini dirancang untuk banyak SNP independen.", "warn")),

      lesson("popgen_save", "6e", "Population structure 5: save, and compare with Arlequin and DnaSP", "Struktur populasi 5: menyimpan, dan membandingkan dengan Arlequin dan DnaSP", "10 min",
             go = "popgen", go_en = "Open Population structure", go_id = "Buka Struktur populasi", icon = "scale-balanced",
        steps_ol(
          li("Result tab: choose a folder and press 'Save all results': the tables (CSV), the HTML report (Arlequin-like), the text report (DnaSP-like), the Arlequin-style folder, the Structure files, the sliding window, the network and the map (SVG + PNG) are written there and the folder opens. 'Save settings and sets (.rds)' / 'Open saved...' keep the populations and the settings of this step in a file.",
             "Tab Result: pilih folder lalu tekan 'Save all results': tabel (CSV), laporan HTML (mirip Arlequin), laporan teks (mirip DnaSP), folder gaya Arlequin, berkas Structure, sliding window, jaringan, dan peta (SVG + PNG) ditulis di sana dan folder terbuka. 'Save settings and sets (.rds)' / 'Open saved...' menyimpan populasi dan pengaturan tahap ini dalam satu berkas."),
          li("The same data give the same values as Arlequin 3.5 and DnaSP 6 for: segregating sites, haplotypes, haplotype diversity (and SD), mean number of pairwise differences (and SD), Theta(Hom), Theta(S), Theta(k), Theta(Pi), Tajima's D, Fu's Fs, the Ewens-Watterson F, the AMOVA variance components and the (pairwise) FST. This was checked with a published DnaSP + Arlequin analysis of 60 Durio zibethinus ITS sequences. Fu & Li's D* and F* follow DnaSP 6: D* with the formula of Fu and Li (1993), F* with the general variances of Achaz (2009), on the biallelic sites; Y* is Achaz (2008).",
             "Data yang sama memberi nilai yang sama dengan Arlequin 3.5 dan DnaSP 6 untuk: situs segregasi, haplotipe, keragaman haplotipe (dan SD), rerata perbedaan berpasangan (dan SD), Theta(Hom), Theta(S), Theta(k), Theta(Pi), Tajima's D, Fs Fu, F Ewens-Watterson, komponen ragam AMOVA, dan FST (berpasangan). Ini telah diperiksa dengan analisis DnaSP + Arlequin yang dipublikasikan atas 60 sekuens ITS Durio zibethinus. D* dan F* Fu & Li mengikuti DnaSP 6: D* dengan rumus Fu dan Li (1993), F* dengan ragam umum Achaz (2009), pada situs biallelic; Y* menurut Achaz (2008)."),
          li("To get exactly the same numbers: use the same populations (the same number of sequences in every sample); tick 'Divide pi by the variable sites' when Arlequin read a DnaSP file with the invariable sites removed; choose the gap rule you used in DnaSP ('complete in each population' = DnaSP default, or pairwise deletion) and set the same 'Allowed level of missing data per site' as in the Arlequin settings.",
             "Agar angkanya persis sama: pakai populasi yang sama (jumlah sekuens sama di setiap sampel); centang 'Divide pi by the variable sites' bila Arlequin membaca berkas DnaSP dengan situs invarian dibuang; pilih aturan gap yang Anda pakai di DnaSP ('complete in each population' = default DnaSP, atau pairwise deletion) dan atur 'Allowed level of missing data per site' yang sama dengan pengaturan Arlequin."),
          li("What is never exactly identical: p-values from simulations or permutations (Tajima's D, Fu's Fs, Fu & Li, R2, Ewens-Watterson, exact test, permutation tests) change a little between runs and programs (a few hundredths with 1000 simulations). Use more simulations for a final value, and 1 000 000 steps for the exact test.",
             "Yang tidak pernah persis sama: nilai p dari simulasi atau permutasi (Tajima's D, Fs Fu, Fu & Li, R2, Ewens-Watterson, exact test, uji permutasi) sedikit berubah antar run dan antar program (beberapa per seratus dengan 1000 simulasi). Pakai lebih banyak simulasi untuk nilai akhir, dan 1 000 000 langkah untuk exact test."),
          li("Not included: Arlequin's hierarchical FDIST2 with several groups, and Hardy-Weinberg / LD tests of multilocus genotypes.",
             "Tidak termasuk: FDIST2 hierarkis Arlequin dengan beberapa grup, serta uji Hardy-Weinberg / LD untuk genotipe multilokus.")),
        goal("Result files for a paper, and a clear idea of which numbers must match DnaSP / Arlequin and which vary a little.",
             "Berkas hasil untuk naskah, dan gambaran jelas angka mana yang harus cocok dengan DnaSP / Arlequin dan mana yang sedikit bervariasi.")),

      lesson("barcode", "7", "Genetic distances: barcode gap and species limits", "Jarak genetik: barcode gap dan batas spesies", "8 min",
             go = "barcode", go_en = "Open Genetic distances", go_id = "Buka Jarak genetik", logo = "img/logo_barcode.png",
        need("An alignment (from MSA, or 'Open aligned FASTA...') and a way to tell the groups apart, such as the species name in the sequence names.",
             "Alignment (dari MSA, atau 'Open aligned FASTA...') dan cara membedakan kelompoknya, misalnya nama spesies pada nama sekuens."),
        steps_ol(
          li("Choose the 'Alignment from'.",
             "Pilih 'Alignment from'."),
          li("Define the groups in 'Groups (species) from': from the sequence names (with the rule 'Species = ') or from a group table. Press 'Apply groups'.",
             "Tentukan kelompok di 'Groups (species) from': dari nama sekuens (dengan aturan 'Species = ') atau dari tabel kelompok. Tekan 'Apply groups'."),
          li("Optional: select the outgroup and press 'Set as outgroup'.",
             "Opsional: pilih outgroup lalu tekan 'Set as outgroup'."),
          li("Choose the 'Model' and press 'Calculate'.",
             "Pilih 'Model' lalu tekan 'Calculate'."),
          li("Open the result tabs: distance distributions (barcode gap), box plots, PCoA, the distance matrix and a tree.",
             "Buka tab hasil: sebaran jarak (barcode gap), box plot, PCoA, matriks jarak, dan pohon."),
          li("Download 'Excel (all tables)', or the figure of the open tab with 'Download the figure of the open tab'.",
             "Unduh 'Excel (all tables)', atau figur dari tab yang terbuka dengan 'Download the figure of the open tab'.")),
        goal("Within- and between-group distances, a barcode gap check and species identification statistics.",
             "Jarak dalam dan antar kelompok, pemeriksaan barcode gap, dan statistik identifikasi spesies.")),

      lesson("primer", "8", "Primer mining and design", "Penambangan dan desain primer", "12 min",
             go = "primer", go_en = "Open Primer design", go_id = "Buka Desain primer", logo = "img/logo_primer.png",
        need("A target: a gene or marker name, or sequences you already have. Internet is needed to search BOLD and NCBI.",
             "Target: nama gen atau marker, atau sekuens yang sudah Anda punya. Butuh internet untuk mencari di BOLD dan NCBI."),
        steps_ol(
          li("Find sequences of your target with 'Search / load' (BOLD or NCBI), or use sequences from the MSA tab or the assembly consensus.",
             "Cari sekuens target dengan 'Search / load' (BOLD atau NCBI), atau pakai sekuens dari tab MSA atau konsensus perakitan."),
          li("Select the ones you want and press 'Send selected to Design'. For degenerate primers, press 'Align (selected or all)' first to see the conserved regions.",
             "Pilih yang Anda inginkan lalu tekan 'Send selected to Design'. Untuk primer degenerate, tekan 'Align (selected or all)' dulu untuk melihat daerah terkonservasi."),
          li("To work on one gene or exons, press 'Get GenBank annotation', choose the gene and press 'Use selected gene in Design' or 'Design exon-exon junction primers for it'.",
             "Untuk bekerja pada satu gen atau ekson, tekan 'Get GenBank annotation', pilih gennya, lalu tekan 'Use selected gene in Design' atau 'Design exon-exon junction primers for it'."),
          li("Set the options and press 'Design primers'. Download the table as CSV.",
             "Atur opsinya lalu tekan 'Design primers'. Unduh tabelnya sebagai CSV."),
          li("Select good primers, press 'Send selected to Quality' and then 'Analyse' to check dimers, hairpins and melting temperature.",
             "Pilih primer yang bagus, tekan 'Send selected to Quality' lalu 'Analyse' untuk memeriksa dimer, hairpin, dan suhu leleh."),
          li("Finish with 'BLAST primers' to see whether they bind only where you want. 'Test connection' checks that NCBI answers.",
             "Akhiri dengan 'BLAST primers' untuk melihat apakah primer hanya menempel di tempat yang Anda inginkan. 'Test connection' memeriksa apakah NCBI menjawab.")),
        goal("Primer pairs with a quality report and a specificity check.",
             "Pasangan primer dengan laporan kualitas dan pemeriksaan spesifisitas.")),

      lesson("open", NULL, "More tools: open and view files", "Alat tambahan: buka dan lihat berkas", "3 min",
             go = "open", go_en = "Open 'Open and view files'", go_id = "Buka 'Buka dan lihat berkas'",
             icon = "eye", count = FALSE,
        steps_ol(
          li("Choose files with 'Choose one or more files', 'Browse files...' or 'Browse folder...'. AB1, FASTA, FASTQ, GenBank, Newick trees and tables are supported.",
             "Pilih berkas dengan 'Choose one or more files', 'Browse files...' atau 'Browse folder...'. Format AB1, FASTA, FASTQ, GenBank, pohon Newick, dan tabel didukung."),
          li("Move between files with '< Previous' and 'Next >'. For sequence files you can 'Edit sequence...', 'Crop...', 'Add sequences...' and 'Reverse complement selected'; 'Undo' and 'Save' are at the top.",
             "Berpindah antarberkas dengan '< Previous' dan 'Next >'. Pada berkas sekuens Anda dapat 'Edit sequence...', 'Crop...', 'Add sequences...', dan 'Reverse complement selected'; 'Undo' dan 'Save' ada di atas."),
          li("'Send to the MSA tab' sends the sequences to the alignment.",
             "'Send to the MSA tab' mengirim sekuens ke alignment.")),
        tip("Nothing here needs an assembly, so it is a quick way to check a chromatogram or a tree file.",
            "Di sini tidak perlu perakitan, jadi ini cara cepat untuk memeriksa kromatogram atau berkas pohon.")),

      lesson("project", NULL, "More tools: project and report", "Alat tambahan: proyek dan laporan", "3 min",
             go = "project", go_en = "Open Project and report", go_id = "Buka Proyek dan laporan",
             icon = "folder-open", count = FALSE,
        steps_ol(
          li("Choose 'What did you sequence?' (the marker, for example COI). It fills in the defaults for the other tabs.",
             "Pilih 'What did you sequence?' (marker, misalnya COI). Ini mengisi nilai awal untuk tab lain."),
          li("'Save project (.rds)' saves your whole work; 'Open a project (.rds)' brings it back. Raw chromatograms are not stored in the project.",
             "'Save project (.rds)' menyimpan seluruh pekerjaan; 'Open a project (.rds)' memulihkannya. Kromatogram mentah tidak disimpan dalam proyek."),
          li("'Download the report (HTML)' makes a report and 'Methods + references (.txt)' gives a Methods paragraph with the references of the methods you used.",
             "'Download the report (HTML)' membuat laporan dan 'Methods + references (.txt)' memberi paragraf Methods beserta pustaka metode yang Anda pakai."),
          li("'Load the synthetic demo' fills the program with example data so you can practise.",
             "'Load the synthetic demo' mengisi program dengan data contoh agar Anda bisa berlatih."))),

      lesson("genbank", NULL, "More tools: GenBank export", "Alat tambahan: ekspor GenBank", "5 min",
             go = "genbank", go_en = "Open GenBank export", go_id = "Buka Ekspor GenBank",
             icon = "file-export", count = FALSE,
        need("Sequences in the list of the MSA tab.",
             "Sekuens di daftar pada tab MSA."),
        steps_ol(
          li("Press 'Load / refresh from the MSA list'.",
             "Tekan 'Load / refresh from the MSA list'."),
          li("Fill the table with the information of each sequence (organism, voucher, locality and so on). 'Fill the empty cells with these defaults' saves typing.",
             "Isi tabel dengan informasi tiap sekuens (organisme, voucher, lokasi, dan sebagainya). 'Fill the empty cells with these defaults' menghemat pengetikan."),
          li("Choose the 'Feature' (and the genetic code for a CDS), then download 'FASTA with source modifiers' and 'Feature table (.tbl)'.",
             "Pilih 'Feature' (dan kode genetik untuk CDS), lalu unduh 'FASTA with source modifiers' dan 'Feature table (.tbl)'."),
          li("Upload the files to NCBI BankIt, or use them with table2asn.",
             "Unggah berkas ke NCBI BankIt, atau pakai dengan table2asn.")),
        goal("Files ready for submission to GenBank.", "Berkas siap dikirim ke GenBank.")),

      lesson("faq", NULL, "Quick answers and common problems", "Jawaban cepat dan masalah umum", "3 min",
             go = "help", go_en = "Open Help", go_id = "Buka Bantuan", icon = "circle-question", count = FALSE,
        tg$ul(
          li("A button does nothing or a tab is empty: the step before it has no data yet. Check the green ticks on Home, or the 'You need' line of the lesson.",
             "Tombol tidak bereaksi atau tab kosong: tahap sebelumnya belum punya data. Periksa tanda centang hijau di Home, atau baris 'Yang diperlukan' pada pelajaran."),
          li("Which steps need internet? BLAST on NCBI and the BOLD / NCBI searches in Primer design. The rest works offline.",
             "Tahap mana yang butuh internet? BLAST di NCBI serta pencarian BOLD / NCBI pada Desain primer. Selebihnya bisa offline."),
          li("I closed the window and lost my work: save often with 'Save project (.rds)' (More > Project and report) or 'Save session (.rds)' in Assembly.",
             "Saya menutup jendela dan pekerjaan hilang: simpan sering dengan 'Save project (.rds)' (More > Proyek dan laporan) atau 'Save session (.rds)' di Perakitan."),
          li("The program is slow to start, or IQ-TREE / trimAl is not found: see the Help page (antivirus and other problems).",
             "Program lambat dibuka, atau IQ-TREE / trimAl tidak ditemukan: lihat halaman Bantuan (antivirus dan masalah lain)."))),

      shiny::div(id = "tt_empty", class = "tt-empty", shiny::icon("magnifying-glass"), " ",
                 t("No lesson matches your search.", "Tidak ada pelajaran yang cocok dengan pencarian Anda.")),
      shiny::div(style = "height:28px")
    )
  )
}
