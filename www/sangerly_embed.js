/* SangeRly inside the Qt window: switches the compact look on (?embed=1) and gives the window a small control object.
   Outside the Qt window (normal browser) this file does nothing. */
(function () {
  var q = new URLSearchParams(location.search);
  if (q.get('embed') !== '1') return;
  document.documentElement.classList.add('sg-embed');
  var st = { saved: null, notes: [], msa: null, nav: null };
  var VERSION = 16;                                   // the Qt window checks this number (sangerly/mainwindow.py: NEED_JS)

  // click the tab titled `titles[i]` inside the tool `tab`; tries for ~6 s because the tool is drawn after the tab changes
  // MSA: "method:<code>" chooses the alignment method (Tools > Multiple sequence alignment > Align with ...), "trim" goes to the trimming step
  function msaMethod(code) {
    if (!(window.Shiny && Shiny.setInputValue)) return false;
    Shiny.setInputValue('msa-set_method', code, { priority: 'event' });
    return true;
  }
  // Genetic distances: "analysis:<code>" chooses the analysis (Tools > Genetic distances > Compute ...)
  function distAnalysis(code) {
    if (!(window.Shiny && Shiny.setInputValue)) return false;
    Shiny.setInputValue('barcode-set_analysis', code, { priority: 'event' });
    // the preferences table is the next step: bring the Compute button into view
    setTimeout(function () { var b = document.getElementById('barcode-run'); if (b && b.scrollIntoView) b.scrollIntoView({ block: 'center' }); }, 900);
    return true;
  }
  // the barcode result tabs (Distribution, Box plot ...) exist only for the barcode analysis: switch to it when one of them is asked for
  function barcodeNeeds(want) {
    var cur = document.querySelector('input[name="barcode-analysis"]:checked');
    if (want === 'Result' || (cur && cur.value === 'barcode')) return false;
    distAnalysis('barcode'); return true;
  }
  function msaTrim() {
    var c = document.getElementById('msa-trim_card');
    if (c && c.offsetParent) { c.scrollIntoView({ behavior: 'smooth', block: 'start' }); return true; }
    if (window.Shiny && Shiny.setInputValue) Shiny.setInputValue('msa-trim_hint', Math.random(), { priority: 'event' });   // "align first"
    return false;
  }
  function openSub(tab, titles, i, tries) {
    var want = titles[i];
    if (/^view:/.test(want)) return;
    if (/^method:/.test(want)) { msaMethod(want.slice(7)); return; }
    if (/^analysis:/.test(want)) { setTimeout(function () { distAnalysis(want.slice(9)); }, 400); return; }
    if (want === 'trim') { setTimeout(msaTrim, 400); return; }
    if (tab === 'barcode' && i === 0 && !openSub._bc && barcodeNeeds(want)) {
      openSub._bc = true; setTimeout(function () { openSub._bc = false; openSub(tab, titles, i, 0); }, 1200); return;
    }
    var pane = document.querySelector('.tab-pane[data-value="' + tab + '"]');
    var hit = null;
    if (pane) [].forEach.call(pane.querySelectorAll('ul.nav a.nav-link, ul.nav a[data-value]'), function (a) {
      if (!hit && (a.getAttribute('data-value') === want || a.textContent.trim() === want)) hit = a;
    });
    if (hit) {
      hit.click();
      if (i + 1 < titles.length) setTimeout(function () { openSub(tab, titles, i + 1, 0); }, 250);
    } else if (tries < 40) {
      setTimeout(function () { openSub(tab, titles, i, tries + 1); }, 150);
    }
  }

  // ---- sliders become number boxes: the slider keeps working underneath, the person types the value ----
  function decimals(step) { var t = String(step); var k = t.indexOf('.'); return k < 0 ? 0 : t.length - k - 1; }
  function numify(el) {
    var $ = window.jQuery; if (!$) return;
    var inst = $(el).data('ionRangeSlider'); if (!inst) return;
    var o = inst.options || {};
    if (el.getAttribute('data-sg-num')) return;
    el.setAttribute('data-sg-num', '1');
    var dt = el.getAttribute('data-data-type');
    if ((dt && dt !== 'number') || (o.values && o.values.length)) return;   // dates and named steps keep their slider
    var dbl = o.type === 'double', step = o.step || 1, dec = decimals(step);
    var box = document.createElement('span'); box.className = 'sg-numbox';
    function mk(key) {
      var n = document.createElement('input'); n.type = 'number'; n.className = 'form-control sg-num';
      n.min = o.min; n.max = o.max; n.step = step; n.setAttribute('data-key', key);
      n.value = (+inst.result[key]).toFixed(dec);
      n.addEventListener('change', function () {
        var v = parseFloat(n.value); if (isNaN(v)) { n.value = (+inst.result[key]).toFixed(dec); return; }
        v = Math.min(o.max, Math.max(o.min, v));
        var cur = { from: +inst.result.from, to: +inst.result.to };
        if (dbl) { if (key === 'from') cur.from = Math.min(v, cur.to); else cur.to = Math.max(v, cur.from); } else cur.from = v;
        $(el).data('ionRangeSlider').update(dbl ? { from: cur.from, to: cur.to } : { from: cur.from });
        $(el).trigger('change');
      });
      return n;
    }
    box.appendChild(mk('from'));
    if (dbl) { var d = document.createElement('span'); d.className = 'sg-numdash'; d.textContent = '–'; box.appendChild(d); box.appendChild(mk('to')); }
    var unit = (o.postfix || '').trim(); if (unit) { var u = document.createElement('span'); u.className = 'sg-numunit'; u.textContent = unit; box.appendChild(u); }
    el.classList.add('sg-numified');
    el.parentNode.insertBefore(box, el.nextSibling);
  }
  function scanSliders() {
    var $ = window.jQuery; if (!$) return;
    [].forEach.call(document.querySelectorAll('input.js-range-slider'), function (el) {
      if (!el.getAttribute('data-sg-num')) numify(el);
      var box = el.nextSibling; if (!(box && box.classList && box.classList.contains('sg-numbox'))) return;
      var inst = $(el).data('ionRangeSlider'); if (!inst) return;
      var o = inst.options || {}, dec = decimals(o.step || 1);       // follows changes made from R (updateSliderInput)
      [].forEach.call(box.querySelectorAll('input'), function (n) {
        if (document.activeElement !== n) { var v = (+inst.result[n.getAttribute('data-key')]).toFixed(dec); if (n.value !== v) n.value = v; }
      });
    });
  }
  setInterval(scanSliders, 600);

  // ---- Assembly tab: the "1. Input" panel is replaced by a card that opens the Qt assembly window; trimming / session become "Advanced" ----
  function qt(cmd) { console.info('SGQT:' + cmd); }            // the Qt window listens to the console (see sangerly/webhost.py)
  function decorateAssembly() {
    var pane = document.querySelector('.tab-pane[data-value="assembly"]'); if (!pane) return;
    var secs = pane.querySelectorAll('.sidebar-sec'); if (!secs.length) return;
    var input = secs[0];
    if (!input.getAttribute('data-sg-card')) {
      input.setAttribute('data-sg-card', '1'); input.classList.add('sg-input-orig');
      var card = document.createElement('div'); card.className = 'sidebar-sec sg-input-card';
      card.innerHTML = '<div class="sg-ic-title">Input</div>' +
        '<div class="sg-ic-text">Choose the mode (de novo or reference-guided), the .ab1 files and the sample names in the assembly window. You can also drop .ab1 files or folders anywhere on this window.</div>' +
        '<button type="button" class="btn btn-primary sg-ic-btn" data-cmd="assembly">Add .ab1 files…</button> ';
      [].forEach.call(card.querySelectorAll('button'), function (b) { b.addEventListener('click', function () { qt(b.getAttribute('data-cmd')); }); });
      input.parentNode.insertBefore(card, input);
    }
    // "How files are grouped into samples" + its explanation (the assembly window groups the files): hide the nearest row only
    var gm = document.getElementById('asm-group_mode');
    if (gm && !gm.getAttribute('data-sg-hide')) { gm.setAttribute('data-sg-hide', '1'); var row = gm.closest('.row'); if (row) row.classList.add('sg-hide'); }
    // the card keeps one button only (the assembly window does the rest)
    [].forEach.call(document.querySelectorAll('.sg-input-card .sg-ic-btn'), function (b, k) { if (k > 0) b.style.display = 'none'; });
    for (var i = 1; i < secs.length; i++) {
      if (secs[i].getAttribute('data-sg-adv')) continue;
      secs[i].setAttribute('data-sg-adv', '1'); secs[i].classList.add('sg-adv');
      if (i === 2) { var sm = secs[i].querySelector('summary'); if (sm && /^[23]\./.test(sm.textContent.trim())) sm.textContent = 'Session'; }   // no numbers in the titles
      if (i === 1) secs[i].classList.add('sg-hide');      // "2. Trimming and thresholds": these settings are in the assembly window
      var d = secs[i].querySelector('details'); if (d && i === 1) d.removeAttribute('open');   // trimming: closed until needed
    }
  }
  setInterval(decorateAssembly, 600);

  // ---- the Qt dialogs set the controls of a tool and press its button: {set: {id: value}, click: id, wait: ms} -----------------------------
  function setVal(id, v) {
    var el = document.getElementById(id); if (!el) return false;
    var ev = function (t) { el.dispatchEvent(new Event(t, { bubbles: true })); };
    if (el.classList.contains('shiny-input-radiogroup')) {
      var hit = false;
      [].forEach.call(el.querySelectorAll('input[type=radio]'), function (r) {
        var on = r.value === String(v); r.checked = on; if (on) { hit = true; r.dispatchEvent(new Event('change', { bubbles: true })); }
      });
      return hit;
    }
    if (el.tagName === 'SELECT') { if (el.selectize) el.selectize.setValue(String(v), false); else { el.value = String(v); ev('change'); } return true; }
    if ((el.type || '').toLowerCase() === 'checkbox') { el.checked = !!v; ev('change'); return true; }
    if (el.classList.contains('js-range-slider') && window.jQuery && jQuery(el).data('ionRangeSlider')) {
      jQuery(el).data('ionRangeSlider').update({ from: +v }); jQuery(el).trigger('change'); return true;
    }
    el.value = (v === null || v === undefined) ? '' : String(v); ev('input'); ev('change'); return true;
  }
  function popReady() {
    var r = document.querySelector('#popgen-seq_tbl tbody tr td');
    return !!(r && !/dataTables_empty/.test(r.className));
  }
  function applyInputs(o) {
    o = o || {}; var set = o.set || {};
    if (o.fasta) shinySet('tree-qt_fasta', o.fasta);
    Object.keys(set).forEach(function (k) { setVal(k, set[k]); });
    if (o.click) setTimeout(function () { var b = document.getElementById(o.click); if (b) b.click(); }, o.wait || 450);
    (o.send || []).forEach(function (x) { setTimeout(function () { shinySet(x.name, x.value); }, x.wait || 0); });
    (o.clicks || []).forEach(function (x) {
      var press = function () { setTimeout(function () { var b = document.getElementById(x.id); if (b) b.click(); }, x.wait || 450); };
      if (!x.ready) { press(); return; }
      var n = 0;                                              // wait (max 30 s) until the sequence table of the populations card has rows
      (function poll() { if (popReady() || n++ > 100) press(); else setTimeout(poll, 300); })();
    });
    return true;
  }
  function shinySet(name, obj) {
    if (!(window.Shiny && Shiny.setInputValue)) return false;
    Shiny.setInputValue(name, obj, { priority: 'event' }); return true;
  }

  // ---- the tabs inside a tool are chosen in the Qt window (document tabs), so the strip of the page is hidden ----------------------------
  var TOOLS = ['assembly', 'blast', 'msa', 'tree', 'vis', 'popgen', 'barcode', 'primer'];
  function mainNav(pane) { return pane.querySelector('ul.nav.nav-tabs'); }
  function decorateTabs() {
    TOOLS.forEach(function (t) {
      var pane = document.querySelector('.tab-pane[data-value="' + t + '"]'); if (!pane) return;
      var nav = mainNav(pane); if (nav) nav.classList.add('sg-maintabs');
    });
  }
  function activeSub() {
    var a = document.querySelector('.navbar-nav a.nav-link.active[data-value]'); if (!a) return null;
    var pane = document.querySelector('.tab-pane[data-value="' + a.getAttribute('data-value') + '"]'); if (!pane) return null;
    var nav = mainNav(pane); if (!nav) return null;
    var l = nav.querySelector('a.nav-link.active, a.active'); if (!l) return null;
    return l.getAttribute('data-value') || l.textContent.trim();
  }

  // ---- MSA: the method, the sequences and the trimming are chosen in the Qt dialogs; the page keeps the results only -------------------------
  function decorateMsa() {
    var pane = document.querySelector('.tab-pane[data-value="msa"]'); if (!pane) return;
    var anchor = document.getElementById('msa-trim_card') || pane.querySelector('.msa-methods') || document.getElementById('msa-run');
    if (!anchor) return;
    var col = anchor.closest('[class*="col-sm-"]'); if (!col || col.getAttribute('data-sg-hide')) return;
    col.setAttribute('data-sg-hide', '1'); col.classList.add('sg-hide2');
    var next = col.nextElementSibling; if (next) next.classList.add('sg-full');
  }
  // ---- Phylogeny: "Build tree" dialog replaces the cards 1-3; the display card stays; IQ-TREE path box only when IQ-TREE is not found ----
  function decorateTree() {
    var pane = document.querySelector('.tab-pane[data-value="tree"]'); if (!pane) return;
    var cards = pane.querySelectorAll('.card');
    [].forEach.call(cards, function (c) {
      var h = c.querySelector('.card-header'); var t = h ? h.textContent.trim() : '';
      if (/^1\.\s*Alignment/.test(t) || /^2\.\s*Neighbo/.test(t)) c.classList.add('sg-hide2');
      if (/^3\.\s*IQ-TREE/.test(t)) {
        var p = document.getElementById('tree-iq_path_ui'), f = document.getElementById('tree-iq_found');
        if (p && p.parentNode) [].forEach.call(p.parentNode.children, function (ch) { if (ch !== p && ch !== f) ch.classList.add('sg-hide2'); });
        var need = !!(p && p.textContent.trim());
        c.classList.toggle('sg-hide2', !need);
        if (h && need && h.textContent !== 'IQ-TREE program') h.textContent = 'IQ-TREE program';
      }
      if (/^4\.\s*Display/.test(t) && h) h.textContent = 'Display options';
    });
  }
  // ---- numbered cards "2. ...", "3. ..." fold with a triangle (card 1 stays open); the state is remembered -------------------------------------
  // Folded = the body is clipped to zero height (see the CSS), not display:none, so Shiny keeps the outputs inside it alive.
  function foldKey(card, title) {
    var tool = card.closest(TOOLS.map(function (t) { return '.tab-pane[data-value="' + t + '"]'; }).join(','));
    return 'sgFold:' + (tool ? tool.getAttribute('data-value') : '') + ':' + title.replace(/\s+/g, ' ');
  }
  function foldGet(k) { try { return localStorage.getItem(k) === 'open'; } catch (e) { return false; } }
  function foldSet(k, open) { try { localStorage.setItem(k, open ? 'open' : 'closed'); } catch (e) { /* no storage: fine */ } }
  function foldCards() {
    [].forEach.call(document.querySelectorAll('.tab-pane .card > .card-header'), function (h) {
      var card = h.parentNode, t = h.textContent.trim(), m = /^(\d+)\.\s/.exec(t);
      if (card.getAttribute('data-sg-fold')) {
        if (!h.querySelector('.sg-tri') && !m) { card.removeAttribute('data-sg-fold'); card.classList.remove('sg-fold', 'sg-collapsed'); }
        return;
      }
      if (!m || m[1] === '1' || card.classList.contains('sg-nofold') || card.closest('.tab-pane[data-value="tree"]')) return;
      var key = foldKey(card, t);
      card.setAttribute('data-sg-fold', '1'); card.classList.add('sg-fold');
      card.classList.toggle('sg-collapsed', !foldGet(key));
      var tri = document.createElement('span'); tri.className = 'sg-tri'; tri.setAttribute('aria-hidden', 'true'); h.appendChild(tri);
      h.setAttribute('role', 'button'); h.setAttribute('tabindex', '0');
      h.setAttribute('aria-expanded', card.classList.contains('sg-collapsed') ? 'false' : 'true');
      var toggle = function () {
        var c = card.classList.toggle('sg-collapsed');
        h.setAttribute('aria-expanded', c ? 'false' : 'true'); foldSet(key, !c);
        setTimeout(function () { window.dispatchEvent(new Event('resize')); }, 30);   // tables and plots inside measure themselves again
      };
      h.addEventListener('click', toggle);
      h.addEventListener('keydown', function (e) { if (e.key === 'Enter' || e.key === ' ') { e.preventDefault(); toggle(); } });
    });
  }
  setInterval(function () { decorateTabs(); decorateMsa(); decorateTree(); foldCards(); }, 1000);
  window.sgQt = {
    // open a tab by its id: home, open, assembly, blast, msa, tree, vis, popgen, barcode, primer, project, genbank, tutorial, help, about
    // sub: optional title of a tab inside the tool (or a list of titles for nested tabs), e.g. go('tree', 'NJ tree')
    go: function (tab, sub) {
      if (window.Shiny && Shiny.setInputValue) Shiny.setInputValue('go_tab', tab, { priority: 'event' });
      if (sub) openSub(tab, [].concat(sub), 0, 0);
    },
    // result of the "New assembly" dialog: {mode, reference, files:[{name,path,sample}], params:{...}, out_dir, save:[...]}
    assemble: function (obj) {
      if (!(window.Shiny && Shiny.setInputValue)) return false;
      Shiny.setInputValue('go_tab', 'assembly', { priority: 'event' });
      Shiny.setInputValue('asm-qt_cmd', obj, { priority: 'event' });
      return true;
    },
    // read by the window every second: steps that have results (green ticks) and the open tab
    status: function () {
      var ids = [].map.call(document.querySelectorAll('.sg-card.done[data-go]'), function (c) { return c.getAttribute('data-go'); });
      var a = document.querySelector('.navbar-nav a.nav-link.active[data-value]');
      var sel = document.getElementById('asm-sample_sel');
      var samples = sel ? [].map.call(sel.options, function (o) { return o.value; }) : [];
      return JSON.stringify({ ready: !!(window.Shiny && Shiny.shinyapp && Shiny.shinyapp.isConnected && Shiny.shinyapp.isConnected()),
                              done: ids, tab: a ? a.getAttribute('data-value') : null, saved: st.saved,
                              samples: samples, cur: sel ? sel.value : null, sub: activeSub(), msa: st.msa, nav: st.nav, ver: VERSION });
    },
    // choose a sample of the open assembly (the selector of the Assembly tab)
    msaMethod: msaMethod, msaTrim: msaTrim, distAnalysis: distAnalysis,
    openPaths: function (a) { return shinySet('open-qt_open', { paths: a }); },
    // File menu of the Qt window: the buttons of the Open tool (Save, Save as, Undo, Send to MSA, Rename ...) and Import
    openCmd: function (id) { var b = document.getElementById('open-' + id); if (!b) return false; b.click(); return true; },
    importSeqs: function (o) { return shinySet('open-qt_import', o); },
    alnPath: function (p) { if (!(window.Shiny && Shiny.setInputValue)) return false; Shiny.setInputValue('barcode-qt_aln_path', p, { priority: 'event' }); return true; },
    // dialogs of the Qt window
    msaRun: function (o) { return shinySet('msa-qt_align', o); },
    msaTrimRun: function (o) { return shinySet('msa-qt_trim', o); },
    apply: applyInputs,
    classic: function (on) { document.documentElement.classList.toggle('sg-classic', !!on); return true; },
    showTabs: function (on) { document.documentElement.classList.toggle('sg-showtabs', !!on); return true; },
    pickSample: function (name) {
      var sel = document.getElementById('asm-sample_sel'); if (!sel) return false;
      if (sel.selectize) sel.selectize.setValue(name); else { sel.value = name; window.jQuery && jQuery(sel).trigger('change'); }
      return true;
    }
  };
  (function wait() {
    if (!(window.Shiny && Shiny.addCustomMessageHandler)) return setTimeout(wait, 30);
    Shiny.addCustomMessageHandler('sg_saved', function (m) { st.saved = m; });
    Shiny.addCustomMessageHandler('sg_msa_state', function (m) { st.msa = m; });
    Shiny.addCustomMessageHandler('sg_nav', function (m) { st.nav = m; });
    Shiny.addCustomMessageHandler('sg_click', function (m) {     // press a button of the page (e.g. 'Design primers' after 'Send to Primer design')
      setTimeout(function () { var b = document.getElementById(m.id); if (b) b.click(); }, m.wait || 300);
    });
  })();
  // links that open a new window go to the system browser (the Qt side handles it); nothing to do here
})();
