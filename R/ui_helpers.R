# Small building blocks that keep the pages short. The look (CSS) and the behaviour (JavaScript) are embedded below, so nothing
# has to be copied to the www folder: the first helper that is placed on a page adds them once (shiny::singleton).
#
#   sg_dlgroup(..., label = "Save")   download buttons that are shown as ONE menu "Save v" (click to open, click outside to close).
#                                      Without JavaScript the buttons are simply shown one after the other.
#   sg_fold(title, ..., open = FALSE)  a card that can be folded (click on its title). Everything inside keeps working while it is folded.
#   sg_opts(title, ..., open = FALSE)  a light folding block for a group of options.
#   sg_hint(text)                      a small (i) icon; the long explanation appears when the mouse rests on it.
#   sg_findbtn(ns)                     small "Find  Ctrl+F" button of the strip above an editor.
#   sg_findbar(ns)                     the find bar: hidden until Ctrl+F; text box, mismatches, Find, previous / next hit, "3 / 12", close.
#                                      Needs sg_find_server() in the module server. Input ids: find_pat, find_mm, find_go, find_prev,
#                                      find_next, find_close; outputs: find_count, find_tbl.
#   sg_find_server(...)                the search itself (both strands) for sg_findbar().

.sg_ui_css <- r"-----(/* SangeRly - small interface helpers: Save menu, folding cards, find bar (Ctrl+F). Appearance only. */

/* "Save v" menu (the script builds it from .sg-dlgroup) */
.sg-dlgroup{display:inline-flex;flex-wrap:wrap;gap:6px;margin:0 6px 4px 0}
.sg-dlgroup.sg-save{display:inline-block;position:relative}
.sg-save-btn{font-weight:600}
.sg-save-menu{min-width:230px;font-size:13px;padding:5px;border-radius:10px;border:1px solid #d8e0ea;box-shadow:0 8px 24px rgba(15,23,42,.14);z-index:2000}
.sg-save-menu .dropdown-item{border-radius:6px;padding:6px 10px;cursor:pointer;white-space:normal;display:flex;align-items:center;gap:8px}
.sg-save-menu .dropdown-item i{width:16px;color:#64748b;flex:none}
.sg-save-menu .dropdown-item-text{padding:4px 10px}
.sg-save-menu .form-group{margin-bottom:4px}
.sg-save-menu .form-group label{font-size:12px;margin-bottom:2px;color:#52606d}

/* bar with a few small tools (Save / Open) */
.sg-bar{display:flex;align-items:center;gap:10px;flex-wrap:wrap;background:#f6f8fa;border:1px solid #e6ebf1;border-radius:8px;padding:5px 10px;margin-bottom:8px}
.sg-bar .sg-dlgroup{margin:0}
.sg-open .form-group{margin:0}
.sg-open .form-control{display:none}
.sg-open .input-group{display:block}
.sg-open .btn{padding:.25rem .6rem;font-size:.82rem}
.sg-open .progress{margin:4px 0 0;height:6px}

/* a card that folds */
details.sg-fold>summary.card-header{cursor:pointer;list-style:none;display:flex;align-items:center;gap:8px;user-select:none}
details.sg-fold>summary.card-header::-webkit-details-marker{display:none}
details.sg-fold>summary.card-header::after{content:'\25BE';margin-left:auto;color:#64748b;transition:transform .15s ease}
details.sg-fold:not([open])>summary.card-header::after{transform:rotate(-90deg)}
details.sg-fold:not([open]){margin-bottom:.6rem!important}
details.sg-fold:not([open])>summary.card-header{border-bottom:0!important;border-radius:var(--sg-r,8px)!important}

/* a light folding block of options */
details.sg-opts{border:1px solid #e1e7ee;border-radius:8px;padding:2px 10px;margin:6px 0 8px;background:#fbfcfe}
details.sg-opts>summary{cursor:pointer;font-weight:600;font-size:13px;color:#1f3c5f;padding:5px 0;list-style:none;user-select:none}
details.sg-opts>summary::-webkit-details-marker{display:none}
details.sg-opts>summary::before{content:'\25B8';display:inline-block;margin-right:6px;transition:transform .15s ease}
details.sg-opts[open]>summary::before{transform:rotate(90deg)}
details.sg-opts[open]{padding-bottom:6px}

/* small (i) icon with a long explanation */
.sg-hint{color:#64748b;cursor:help;margin-left:6px;font-size:.92em}
.sg-hint:hover,.sg-hint:focus{color:#2563eb}

/* thin strip above an editor: Reload, (i) and the Find button */
.sg-edbar{display:flex;align-items:center;gap:8px;flex-wrap:wrap;margin-bottom:6px}
.sg-edbar .btn{padding:.15rem .6rem;font-size:.82rem}
.sg-edbar .sg-hint{margin-left:0}
.sg-kbd{font-size:10.5px;color:#64748b;border:1px solid #cbd5e1;border-radius:4px;padding:0 4px;margin-left:6px;background:#f8fafc}

/* find bar: hidden until Ctrl+F (or the Find button), Esc closes it */
.sg-findwrap{display:none;margin-bottom:6px}
.sg-findwrap.open{display:block}
.sg-find{display:flex;align-items:center;gap:8px;flex-wrap:wrap;border:1px solid #cfd8e3;border-radius:8px;padding:4px 10px;background:#fff;box-shadow:0 2px 8px rgba(15,23,42,.06)}
.sg-find .form-group{margin:0}
.sg-find .form-control,.sg-find .form-select,.sg-find select{height:30px;padding:2px 8px;font-size:13px}
.sg-find .selectize-input{min-height:30px;padding:3px 8px;font-size:13px}
.sg-find .sg-grow{flex:1}
.sg-find .btn{padding:.2rem .7rem;font-size:.84rem}
.sg-find .btn-outline-secondary{padding:.2rem .5rem}
.sg-find-count{font-size:12.5px;color:#52606d;min-width:56px;text-align:center;font-variant-numeric:tabular-nums}
.sg-find-x{border:0;background:transparent;font-size:22px;line-height:1;color:#64748b;padding:0 6px;cursor:pointer}
.sg-find-x:hover{color:#dc2626}
.sg-find-res:empty{display:none}
.sg-findwrap.sg-findcard>.card{margin-bottom:0}

/* a side column (options) becomes a rail of tabs on the right (like the Options Panel of UGENE): click a tab to open it, click again to close */
.sg-siderow{flex-wrap:nowrap!important;align-items:flex-start}
.sg-siderow>.sg-maincol{order:1;flex:1 1 0!important;width:auto!important;max-width:none!important;min-width:0}
.sg-siderow>.sg-sidecol{order:2;display:none;flex:0 0 370px!important;width:370px!important;max-width:370px!important;max-height:calc(100vh - 175px);overflow-y:auto;border-left:1px solid #d5dce6;background:#f8fafc;padding:6px!important}
.sg-siderow.sg-side-open>.sg-sidecol{display:block}
.sg-sidecol .sg-pane:not(.sg-on){display:none!important}
.sg-maincol .dataTables_wrapper{overflow-x:auto}
.sg-rail{order:3;flex:0 0 34px;display:flex;flex-direction:column;gap:4px;padding:6px 3px;background:#eef1f5;border-left:1px solid #d5dce6;max-height:calc(100vh - 175px);overflow:auto;align-self:stretch}
.sg-railbtn{flex:0 0 auto;width:26px;height:26px;padding:0;border:1px solid transparent;background:transparent;border-radius:5px;cursor:pointer;color:#1f3c5f;font-size:14px;display:flex;align-items:center;justify-content:center}
.sg-sidehead{position:sticky;top:-6px;z-index:2;margin:-6px -6px 6px;padding:6px 10px;background:#e9eef5;border-bottom:1px solid #d5dce6;font-size:11px;font-weight:700;letter-spacing:.06em;text-transform:uppercase;color:#475569}
.sg-railbtn:hover{background:#e0e7f0}
.sg-railbtn.active{background:#fff;border-color:#c6d1de;color:#2563eb}

/* hold the left mouse button on a table and drag to move it */
.dataTables_scrollBody,.sg-pan{cursor:grab}
body.sg-panning,body.sg-panning *{cursor:grabbing!important;user-select:none!important}

/* desktop density (inside the Qt window): smaller type, tighter boxes, flat panels */
html.sg-embed{font-size:13px;--sg-r:5px}
html.sg-embed body{font-size:13px}
html.sg-embed body>.container-fluid{padding:8px 10px 14px!important}
html.sg-embed .card,html.sg-embed .well,html.sg-embed .ibox{box-shadow:none!important;border-radius:5px!important}
html.sg-embed .card{margin-bottom:8px}
html.sg-embed .card-header{padding:5px 10px 5px 14px;font-size:12px;border-radius:5px 5px 0 0!important}
html.sg-embed .card-header::before{top:5px;bottom:5px;width:3px}
html.sg-embed .card-body{padding:8px 10px}
html.sg-embed .form-control,html.sg-embed .form-select{font-size:12.5px;padding:2px 7px;min-height:26px;height:auto;border-radius:4px}
html.sg-embed .selectize-input{min-height:26px;padding:3px 7px;font-size:12.5px;border-radius:4px}
html.sg-embed .btn{font-size:12.5px;padding:3px 10px;border-radius:4px}
html.sg-embed .btn-sm{padding:2px 8px;font-size:12px}
html.sg-embed label,html.sg-embed .control-label{font-size:12px;margin-bottom:2px}
html.sg-embed .form-group,html.sg-embed .shiny-input-container:not(.shiny-input-container-inline){margin-bottom:6px}
html.sg-embed .form-check,html.sg-embed .radio,html.sg-embed .checkbox{margin:1px 0;min-height:0;font-size:12.5px}
html.sg-embed h1{font-size:19px;margin:2px 0 4px}
html.sg-embed h2{font-size:17px}html.sg-embed h3{font-size:15px}html.sg-embed h4{font-size:14px}html.sg-embed h5{font-size:13px}
html.sg-embed .help-block,html.sg-embed .form-text,html.sg-embed p{font-size:12px}
html.sg-embed .nav-tabs .nav-link{padding:4px 12px;font-size:12.5px}
html.sg-embed table.dataTable{font-size:12px}
html.sg-embed table.dataTable td,html.sg-embed table.dataTable th{padding:3px 8px}
html.sg-embed .alert{padding:6px 10px;font-size:12.5px}
)-----"

.sg_ui_js <- r"-----(/* SangeRly - small interface helpers (no tool is changed).
   1) .sg-dlgroup: the download buttons inside become ONE menu (Save v). The buttons themselves are moved, not copied,
      so the downloads work exactly as before.
   2) .sg-findwrap: a find bar that is hidden until Ctrl+F (or its Find button) is pressed; Esc or the cross closes it.
      In the text box: Enter = find (a second Enter on the same text = next hit), Shift+Enter = previous hit. */
(function () {
  'use strict';
  if (window.__sgUi) return;
  window.__sgUi = true;

  /* ---------------------------------------------------------------- 1) Save menu */
  function build(g) {
    if (g.getAttribute('data-sg-done')) return;
    var items = Array.prototype.slice.call(g.children);
    if (!items.length) return;
    g.setAttribute('data-sg-done', '1');
    var label = g.getAttribute('data-label') || 'Save';
    var ic = g.getAttribute('data-icon') || 'floppy-disk';
    g.classList.add('dropdown', 'sg-save');
    g.removeAttribute('style');

    var btn = document.createElement('button');
    btn.type = 'button';
    btn.className = 'btn btn-sm btn-outline-secondary dropdown-toggle sg-save-btn';
    btn.setAttribute('data-bs-toggle', 'dropdown');
    btn.setAttribute('data-bs-auto-close', 'outside');
    btn.setAttribute('aria-expanded', 'false');
    btn.innerHTML = '<i class="fa fa-' + ic + '" role="presentation"></i> <span class="sg-save-lbl"></span>';
    btn.querySelector('.sg-save-lbl').textContent = label;

    var ul = document.createElement('ul');
    ul.className = 'dropdown-menu sg-save-menu';
    items.forEach(function (el) {
      var li = document.createElement('li');
      var isLink = el.matches && el.matches('a, button');
      if (isLink) {
        el.className = 'dropdown-item ' + (el.classList.contains('shiny-download-link') ? 'shiny-download-link' : '') +
                       (el.classList.contains('action-button') ? ' action-button' : '');
        li.appendChild(el);
      } else {
        li.className = 'dropdown-item-text';
        li.appendChild(el);
      }
      ul.appendChild(li);
    });
    g.innerHTML = '';
    g.appendChild(btn);
    g.appendChild(ul);

    try {   // the menu may sit inside a box that scrolls: let it float above everything
      if (window.bootstrap && window.bootstrap.Dropdown) {
        new window.bootstrap.Dropdown(btn, { popperConfig: function (d) { d.strategy = 'fixed'; return d; } });
      }
    } catch (e) { /* the plain menu still works */ }

    ul.addEventListener('click', function (e) {   // after a download or an action the menu closes
      var a = e.target.closest ? e.target.closest('a.dropdown-item, button.dropdown-item') : null;
      if (!a) return;
      setTimeout(function () {
        try { window.bootstrap.Dropdown.getOrCreateInstance(btn).hide(); } catch (err) { /* ignore */ }
      }, 150);
    });
  }

  var pending = false;
  function scan() {
    pending = false;
    var gs = document.querySelectorAll('.sg-dlgroup:not([data-sg-done])');
    for (var i = 0; i < gs.length; i++) build(gs[i]);
  }
  function schedule() {
    if (pending) return;
    pending = true;
    window.setTimeout(scan, 150);
  }
  function start() {
    scan();
    new MutationObserver(schedule).observe(document.documentElement, { childList: true, subtree: true });
  }
  if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', start); else start();

  /* ---------------------------------------------------------------- 2) find bar (Ctrl+F) */
  function wraps() { return Array.prototype.slice.call(document.querySelectorAll('.sg-findwrap')); }
  function shown(w) { var p = w.parentElement; return !!p && p.getClientRects().length > 0; }   // is the tab that holds it on screen?
  function focusFind(w) {
    var i = w.querySelector('input[type=text]');
    if (i) { i.focus(); i.select(); }
  }
  function openFind(w) { w.classList.add('open'); focusFind(w); }
  function closeFind(w) {
    w.classList.remove('open');
    var id = w.getAttribute('data-close');       // the server removes the colours of the hits
    if (id && window.Shiny && window.Shiny.setInputValue) window.Shiny.setInputValue(id, Date.now(), { priority: 'event' });
  }

  document.addEventListener('keydown', function (e) {
    var k = (e.key || '').toLowerCase();
    if ((e.ctrlKey || e.metaKey) && !e.shiftKey && !e.altKey && k === 'f') {
      var ws = wraps().filter(shown);
      if (!ws.length) return;                    // no find bar on this page: the browser keeps its own Ctrl+F
      var act = document.activeElement;
      var w = ws.filter(function (x) { return x.contains(act); })[0] || ws[0];
      e.preventDefault(); e.stopPropagation();
      openFind(w);
      return;
    }
    var t = e.target && e.target.closest ? e.target.closest('.sg-findwrap') : null;
    if (!t) return;
    if (k === 'escape' && t.classList.contains('open')) {
      e.preventDefault(); e.stopPropagation();
      closeFind(t);
      return;
    }
    if (k === 'enter' && e.target.matches && e.target.matches('input[type=text]')) {
      e.preventDefault(); e.stopPropagation();
      var inp = e.target, sel = t.querySelector('select');
      var key = inp.value + '|' + (sel ? sel.value : '');
      var go = t.querySelector('button.btn-primary');
      var nxt = t.querySelector('.sg-find-next'), prv = t.querySelector('.sg-find-prev');
      if (window.Shiny && window.Shiny.setInputValue && inp.id) window.Shiny.setInputValue(inp.id, inp.value, { priority: 'event' });
      if (t.getAttribute('data-last') === key && (nxt || prv)) {        // same text again: go to the next / previous hit
        var b = e.shiftKey ? prv : nxt;
        if (b) b.click();
      } else {
        t.setAttribute('data-last', key);
        setTimeout(function () { if (go) go.click(); }, 60);
      }
    }
  }, true);

  document.addEventListener('click', function (e) {
    var el = e.target && e.target.closest ? e.target : null;
    if (!el) return;
    var b = el.closest('.sg-find-open');
    if (b) {
      var w = document.getElementById(b.getAttribute('data-find'));
      if (w) { if (w.classList.contains('open')) closeFind(w); else openFind(w); }
      return;
    }
    var x = el.closest('.sg-find-x');
    if (x) {
      var w2 = x.closest('.sg-findwrap');
      if (w2) closeFind(w2);
    }
    var go = el.closest('.sg-findwrap button.btn-primary');           // a new search: next Enter on the same text goes to the next hit
    if (go) { var w3 = go.closest('.sg-findwrap'); var i3 = w3.querySelector('input[type=text]'), s3 = w3.querySelector('select');
              if (i3) w3.setAttribute('data-last', i3.value + '|' + (s3 ? s3.value : '')); }
  });

  /* ---------------------------------------------------------------- 3) side column -> rail of tabs (UGENE style) */
  function tx(t) { return (t || '').replace(/\s+/g, ' ').trim(); }
  function isTag(el, re) { return re.test(el.tagName); }
  function title(el) {
    var h = el.querySelector(':scope > .card-header, :scope > summary, :scope > .card > .card-header');
    return tx(h ? h.textContent : '').replace(/^\d+\.\s*/, '');
  }
  function host(col) {
    var c = col;
    for (var d = 0; d < 3; d++) {
      var k = [].filter.call(c.children, function (x) { return !isTag(x, /^(SCRIPT|STYLE)$/); });
      if (k.length === 1 && !k[0].classList.contains('card')) c = k[0]; else break;
    }
    return c;
  }
  var ICONS = [[/align|sequence record|sequences$/i, 'align-left'], [/taxa|group|population|set/i, 'users'], [/tree|neighbo|iq-?tree|phylo/i, 'sitemap'],
               [/display|appearance|layout|look/i, 'sliders'], [/setting|preference|option|advanced/i, 'gear'], [/run|analysis|compute/i, 'play'],
               [/project|file|input|open/i, 'folder-open'], [/trim/i, 'scissors'], [/session|history/i, 'clock'], [/database|blast/i, 'database'],
               [/primer|design|template|what/i, 'dna'], [/result|table|export|save/i, 'table-cells'], [/colou?r|palette/i, 'palette'], [/search|find|mining|mine/i, 'magnifying-glass'],
               [/map/i, 'earth-americas'], [/quality|check/i, 'circle-check']];
  var SPARE = ['list', 'table', 'flag', 'star', 'bookmark', 'cube', 'layer-group'];
  function iconFor(t, used) {
    var i, c;
    for (i = 0; i < ICONS.length; i++) if (ICONS[i][0].test(t) && used.indexOf(ICONS[i][1]) < 0) return ICONS[i][1];
    for (i = 0; i < SPARE.length; i++) if (used.indexOf(SPARE[i]) < 0) return SPARE[i];
    return 'bars';
  }
  function showPane(row, rail, i) {
    var hd = row.querySelector('.sg-sidehead');
    if (hd) hd.textContent = i >= 0 && rail.children[i] ? (rail.children[i].getAttribute('data-t') || '') : '';
    var all = row.querySelectorAll('.sg-sidecol .sg-pane');
    [].forEach.call(all, function (p) { p.classList.toggle('sg-on', i >= 0 && p.getAttribute('data-sgp') === String(i)); });
    [].forEach.call(rail.children, function (b) { b.classList.toggle('active', b.getAttribute('data-i') === String(i)); });
    row.classList.toggle("sg-side-open", i >= 0);
    try { window.dispatchEvent(new Event('resize')); } catch (e) { /* ignore */ }
    setTimeout(function () { try { window.dispatchEvent(new Event('resize')); } catch (e) { /* ignore */ } }, 120);
  }
  function sideRows() {
    var rows = document.querySelectorAll('.row:not([data-sg-side])');
    for (var n = 0; n < rows.length; n++) {
      var r = rows[n], k = r.children;
      var ok = k.length >= 2 && r.classList.contains('sg-siderow');        // only the rows that ask for it (Assembly and Open keep their own menu)
      if (!ok) continue;
      var hst = host(k[0]);
      var items = [].filter.call(hst.children, function (x) { return !isTag(x, /^(SCRIPT|STYLE)$/); });
      var panes = [], byT = {}, cur = null;
      items.forEach(function (el) {
        if (el.classList.contains('card')) {
          var t = title(el) || 'Options';
          if (byT[t]) byT[t].els.push(el); else { byT[t] = { t: t, els: [el] }; panes.push(byT[t]); }
          cur = null;
        } else {
          if (!cur) { cur = byT['Options'] || { t: 'Options', els: [] }; if (!byT['Options']) { byT['Options'] = cur; panes.push(cur); } }
          cur.els.push(el);
        }
      });
      if (!panes.length) continue;                       // not drawn yet: try again later
      r.setAttribute('data-sg-side', '1');
      r.classList.add('sg-siderow');
      k[0].classList.add('sg-sidecol'); k[1].classList.add('sg-maincol');
      var rail = document.createElement('div'); rail.className = 'sg-rail';
      var used = [];
      panes.forEach(function (p, i) {
        p.els.forEach(function (el) { el.classList.add('sg-pane'); el.setAttribute('data-sgp', i); });
        var b = document.createElement('button'); b.type = 'button'; b.className = 'sg-railbtn';
        b.setAttribute('data-i', i); b.setAttribute('data-t', p.t); b.title = p.t; b.setAttribute('aria-label', p.t);
        var ic = iconFor(p.t, used); used.push(ic); b.innerHTML = '<i class="fa fa-' + ic + '" role="presentation"></i>';
        rail.appendChild(b);
      });
      var hd = document.createElement('div'); hd.className = 'sg-sidehead'; k[0].insertBefore(hd, k[0].firstChild);
      r.appendChild(rail);
      rail.addEventListener('click', function (ev) {
        var b = ev.target.closest('.sg-railbtn'); if (!b) return;
        var rl = ev.currentTarget, row = rl.parentNode, was = b.classList.contains('active');
        showPane(row, rl, was ? -1 : parseInt(b.getAttribute('data-i'), 10));
      });
      showPane(r, rail, 0);
    }
  }
  function unnumber() {
    var h = document.querySelectorAll('.card-header:not([data-sg-nn]),details.sg-fold>summary:not([data-sg-nn])');
    for (var i = 0; i < h.length; i++) { h[i].setAttribute('data-sg-nn', '1');
      var t = h[i].firstChild; if (t && t.nodeType === 3) t.nodeValue = t.nodeValue.replace(/^\s*\d\.\s+/, ''); }
  }
  var sideTimer = null;      // at most 5 scans per second, however busy the page is (a design run updates the page all the time)
  new MutationObserver(function () { if (sideTimer) return; sideTimer = setTimeout(function () { sideTimer = null; unnumber(); sideRows(); }, 200); })
    .observe(document.documentElement, { childList: true, subtree: true });
  if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', sideRows); else sideRows();

  /* ---------------------------------------------------------------- 4) hold the left mouse button on a table and drag to move it */
  var pan = null;
  document.addEventListener('mousedown', function (e) {
    if (e.button !== 0 || !e.target.closest) return;
    if (e.target.closest('input,select,textarea,button,a,label,canvas,svg,img,.msaed-root,.sg-sl,[draggable=true],.selectize-control,.leaflet-container')) return;
    var box = e.target.closest('.dataTables_scrollBody,.dataTables_wrapper,.table-responsive,.sg-pan');
    if (!box) return;
    var sc = box.classList.contains('dataTables_wrapper') ? (box.querySelector('.dataTables_scrollBody') || box) : box;
    if (sc.scrollWidth <= sc.clientWidth + 1 && sc.scrollHeight <= sc.clientHeight + 1) return;
    pan = { s: sc, x: e.clientX, y: e.clientY, l: sc.scrollLeft, t: sc.scrollTop, moved: false };
  }, true);
  window.addEventListener('mousemove', function (e) {
    if (!pan) return;
    var dx = e.clientX - pan.x, dy = e.clientY - pan.y;
    if (!pan.moved && Math.abs(dx) + Math.abs(dy) > 5) { pan.moved = true; document.body.classList.add('sg-panning'); }
    if (pan.moved) { pan.s.scrollLeft = pan.l - dx; pan.s.scrollTop = pan.t - dy; e.preventDefault(); }
  });
  window.addEventListener('mouseup', function () {
    if (pan && pan.moved) {
      document.body.classList.remove('sg-panning');
      var stop = function (ev) { ev.stopPropagation(); ev.preventDefault(); };      // the drag is not a click on a row
      window.addEventListener('click', stop, { capture: true, once: true });
      setTimeout(function () { window.removeEventListener('click', stop, true); }, 50);
    }
    pan = null;
  });
})();
)-----"

.sg_ui_head <- function() {
  shiny::tagList(
    shiny::singleton(shiny::tags$style(shiny::HTML(.sg_ui_css))),
    shiny::singleton(shiny::tags$script(shiny::HTML(.sg_ui_js))))
}

sg_dlgroup <- function(..., label = "Save", icon = "floppy-disk", class = NULL) {
  shiny::tagList(.sg_ui_head(),
    shiny::tags$div(class = paste(c("sg-dlgroup", class), collapse = " "), `data-label` = label, `data-icon` = icon, ...))
}

sg_fold <- function(title, ..., open = FALSE, class = NULL) {
  shiny::tagList(.sg_ui_head(),
    shiny::tags$details(class = paste(c("sg-fold card mb-3", class), collapse = " "), open = if (isTRUE(open)) NA else NULL,
      shiny::tags$summary(class = "card-header fw-semibold", title),
      shiny::div(class = "card-body", ...)))
}

sg_opts <- function(title, ..., open = FALSE) {
  shiny::tagList(.sg_ui_head(),
    shiny::tags$details(class = "sg-opts", open = if (isTRUE(open)) NA else NULL,
      shiny::tags$summary(title), ...))
}

sg_hint <- function(text) {
  shiny::tagList(.sg_ui_head(),
    shiny::tags$span(class = "sg-hint", title = text, tabindex = "0", shiny::icon("circle-info")))
}

sg_findbtn <- function(ns, id = "findwrap") {
  shiny::tagList(.sg_ui_head(),
    shiny::tags$button(type = "button", class = "btn btn-sm btn-outline-secondary sg-find-open", `data-find` = ns(id),
      title = "Find a sequence (Ctrl+F)", shiny::icon("magnifying-glass"), " Find", shiny::span(class = "sg-kbd", "Ctrl+F")))
}

sg_findbar <- function(ns, hint = NULL, id = "findwrap") {
  if (is.null(hint)) hint <- paste("Finds the sequence and its reverse complement in every sequence (strand + = as typed,",
                                   "strand - = reverse complement). The hits are coloured in the editor; click a row of the result,",
                                   "or use the arrows, to go to a hit. Enter = find (again = next hit), Shift+Enter = previous hit, Esc = close.")
  shiny::tagList(.sg_ui_head(),
    shiny::div(id = ns(id), class = "sg-findwrap", `data-close` = ns("find_close"),
      shiny::div(class = "sg-find",
        shiny::icon("magnifying-glass"),
        shiny::textInput(ns("find_pat"), NULL, "", placeholder = "Find a sequence (IUPAC letters), both strands", width = "320px"),
        shiny::selectInput(ns("find_mm"), NULL, c("exact" = 0, "1 mismatch" = 1, "2 mismatches" = 2, "3 mismatches" = 3), width = "130px"),
        shiny::actionButton(ns("find_go"), "Find", class = "btn-primary btn-sm"),
        shiny::actionButton(ns("find_prev"), shiny::icon("chevron-up"), class = "btn-sm btn-outline-secondary sg-find-prev", title = "Previous hit (Shift+Enter)"),
        shiny::actionButton(ns("find_next"), shiny::icon("chevron-down"), class = "btn-sm btn-outline-secondary sg-find-next", title = "Next hit (Enter)"),
        shiny::span(class = "sg-find-count", shiny::textOutput(ns("find_count"), inline = TRUE)),
        sg_hint(hint),
        shiny::span(class = "sg-grow"),
        shiny::tags$button(type = "button", class = "sg-find-x", title = "Close (Esc)", shiny::HTML("&times;"))),
      shiny::div(class = "sg-find-res", DT::DTOutput(ns("find_tbl")))))
}

# The search of sg_findbar(). get_aln() returns a named character vector (the rows of the alignment, gaps as "-") or NULL
# (then "Align the sequences first" / "Open a sequence file first" is shown with msg_none). editor_id = full id of the editor
# (ns("editor")). Returns the reactiveVal that holds the hits, so that the module can clear it (for example when another file opens).
sg_find_server <- function(input, output, session, get_aln, editor_id, msg_none = "Open or align sequences first.") {
  hits <- shiny::reactiveVal(NULL)
  idx <- shiny::reactiveVal(0L)
  goto <- function(h, i) {
    session$sendCustomMessage("msaed_goto", list(id = editor_id, name = h$Sequence[i], row = h$Row[i],
                                                 c0 = h$Start[i] - 1L, c1 = h$End[i] - 1L))
  }
  shiny::observeEvent(input$find_go, {
    aln <- get_aln()
    if (is.null(aln) || !length(aln)) { shiny::showNotification(msg_none, type = "warning"); return() }
    pat <- input$find_pat
    if (is.null(pat) || nchar(gsub("[^A-Za-z]", "", pat)) < 4) {
      shiny::showNotification("Type at least 4 bases.", type = "warning"); hits(NULL); idx(0L); return()
    }
    mm <- suppressWarnings(as.integer(input$find_mm)); if (is.na(mm)) mm <- 0L
    res <- lapply(seq_along(aln), function(i) {
      h <- .find_oligo(aln[[i]], pat, mm)
      if (is.null(h)) return(NULL)
      cbind(Sequence = names(aln)[i], Row = i - 1L, h, stringsAsFactors = FALSE)
    })
    res <- res[!vapply(res, is.null, logical(1))]
    h <- if (length(res)) utils::head(do.call(rbind, res), 1000) else NULL
    hits(h); idx(if (is.null(h)) 0L else 1L)
    if (is.null(h)) {
      session$sendCustomMessage("msaed_unmark", list(id = editor_id))
      shiny::showNotification("Not found (neither the sequence nor its reverse complement).", type = "message")
    } else {
      marks <- lapply(seq_len(nrow(h)), function(j)
        list(name = h$Sequence[j], c0 = as.integer(h$Start[j]) - 1L, c1 = as.integer(h$End[j]) - 1L, s = as.character(h$Strand[j])))
      session$sendCustomMessage("msaed_marks", list(id = editor_id, marks = marks))
      goto(h, 1L)
      shiny::showNotification(sprintf("%d hit(s): coloured in the editor (yellow = as typed, green = reverse complement).", nrow(h)), type = "message")
    }
  })
  output$find_count <- shiny::renderText({
    h <- hits(); i <- idx()
    if (is.null(h) || i < 1L) "" else sprintf("%d / %d", i, nrow(h))
  })
  output$find_tbl <- DT::renderDT({
    h <- hits(); shiny::req(h)
    DT::datatable(h[, c("Sequence", "Strand", "Start", "End", "Mismatches", "Matched")], rownames = FALSE,
                  selection = "single", options = list(pageLength = 6, dom = "tip", scrollX = TRUE))
  })
  shiny::observeEvent(input$find_tbl_rows_selected, {
    h <- hits(); i <- input$find_tbl_rows_selected
    if (is.null(h) || length(i) != 1L) return()
    idx(i); goto(h, i)
  })
  step <- function(d) {
    h <- hits(); if (is.null(h) || !nrow(h)) return()
    i <- ((idx() - 1L + d) %% nrow(h)) + 1L
    idx(i); goto(h, i)
  }
  shiny::observeEvent(input$find_next, step(1L))
  shiny::observeEvent(input$find_prev, step(-1L))
  shiny::observeEvent(input$find_close, {       # the find bar was closed: take the colours off
    hits(NULL); idx(0L)
    session$sendCustomMessage("msaed_unmark", list(id = editor_id))
  })
  hits
}
