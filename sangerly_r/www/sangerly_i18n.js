/* SangeRly - Indonesian for the WHOLE program (the EN | ID switch of the menu bar).
   The texts of the program are English; when "ID" is chosen, every text found in www/i18n_id.json is shown in
   Indonesian: labels, buttons, tab titles, notes, notifications, table headers, tooltips and placeholders, also the ones
   created later. A long paragraph is translated sentence by sentence when the whole paragraph is not in the file.
   Texts that are not in the file stay in English.

   To find what is still missing: choose ID, visit the tabs, then press  Ctrl + Alt + Shift + T
   -> the browser saves "sangerly_untranslated.txt" with every English text that was shown but not translated.
   (Send that file to get the translations; they are added to www/i18n_id.json, no code change needed.) */
(function () {
  "use strict";
  var DICT = null, PATS = [], lang = "en", obs = null;
  var track = [];                                      // translated text nodes: {n: node, o: original, t: translated}
  var missing = Object.create(null);                   // English texts shown in ID mode without a translation
  var SKIP = "script,style,textarea,pre,code,kbd,noscript,.i18n,.no-i18n,.msaed-canvas,canvas,.op-seq,.sg-seq";
  var ATTRS = ["placeholder", "title", "aria-label", "data-original-title"];

  /* a validation message inside a table output used to be drawn over the text below it */
  (function css() {
    var s = document.createElement("style");
    s.textContent = ".datatables.html-widget{min-height:1.9em}" +
      ".shiny-output-error-validation{position:static!important;display:block;padding:4px 0;color:#6c757d}" +
      ".datatables .shiny-output-error{position:static!important}";
    document.head.appendChild(s);
  })();

  function lookup(t) {
    if (DICT && Object.prototype.hasOwnProperty.call(DICT, t)) return DICT[t];
    for (var i = 0; i < PATS.length; i++) { if (PATS[i][0].test(t)) return t.replace(PATS[i][0], PATS[i][1]); }
    return undefined;
  }
  function lookupSmart(t) {
    var m = lookup(t); if (m !== undefined) return m;
    var parts = t.split(/(?<=[.!?])\s+(?=[A-Z(<'"])/);          // sentence by sentence
    if (parts.length < 2) return undefined;
    var out = [], hit = 0;
    parts.forEach(function (p) { var r = lookup(p.trim()); if (r !== undefined) { hit++; out.push(r); } else out.push(p); });
    return hit ? out.join(" ") : undefined;
  }
  function skipEl(el) { return !el || (el.closest && el.closest(SKIP)); }
  function noteMissing(t, node) {
    if (t.length < 3 || t.length > 600 || !/[A-Za-z]{3}/.test(t)) return;
    if (/^[\d\s.,:;%()+\-/*=<>[\]{}|#@_]+$/.test(t)) return;
    var p = node && node.parentNode; if (p && p.closest && p.closest("td,tbody")) return;     // table data, not interface text
    missing[t] = 1;
  }

  function doText(n) {
    var p = n.parentNode; if (!p || skipEl(p)) return;
    var v = n.nodeValue, t = v.trim(); if (!t) return;
    var m = lookupSmart(t);
    if (m === undefined) { noteMissing(t, n); return; }
    if (m === t) return;
    var lead = v.slice(0, v.indexOf(t)), trail = v.slice(v.indexOf(t) + t.length);
    var nv = lead + m + trail;
    track.push({ n: n, o: v, t: nv });
    n.nodeValue = nv;
  }
  function doAttrs(el) {
    if (!el.getAttribute || skipEl(el)) return;
    for (var i = 0; i < ATTRS.length; i++) {
      var a = ATTRS[i], v = el.getAttribute(a); if (!v) continue;
      var key = "data-sg-o-" + a, orig = el.hasAttribute(key) ? el.getAttribute(key) : v;
      var m = lookupSmart(orig.trim());
      if (m === undefined) { noteMissing(orig.trim(), null); continue; }
      if (!el.hasAttribute(key)) el.setAttribute(key, v);
      el.setAttribute(a, m);
    }
  }
  function walk(root) {
    if (!root) return;
    if (root.nodeType === 3) { doText(root); return; }
    if (root.nodeType !== 1) return;
    doAttrs(root);
    var tw = document.createTreeWalker(root, NodeFilter.SHOW_TEXT | NodeFilter.SHOW_ELEMENT, null);
    var n; while ((n = tw.nextNode())) { if (n.nodeType === 3) doText(n); else doAttrs(n); }
  }
  function restore() {
    track = track.filter(function (r) { return r.n.isConnected; });
    track.forEach(function (r) { if (r.n.nodeValue === r.t) r.n.nodeValue = r.o; });
    track = [];
    document.querySelectorAll("[data-sg-o-placeholder],[data-sg-o-title],[data-sg-o-aria-label],[data-sg-o-data-original-title]").forEach(function (el) {
      ATTRS.forEach(function (a) { var k = "data-sg-o-" + a; if (el.hasAttribute(k)) { el.setAttribute(a, el.getAttribute(k)); el.removeAttribute(k); } });
    });
  }
  function pause() { if (obs) obs.disconnect(); }
  function resume() { if (obs) obs.observe(document.body, { childList: true, subtree: true, characterData: true, attributes: false }); }
  function apply() {
    pause();
    try { if (lang === "id") { document.documentElement.setAttribute("lang", "id"); walk(document.body); }
          else { restore(); document.documentElement.setAttribute("lang", "en"); } } finally { resume(); }
  }
  function onMutations(ms) {
    if (lang !== "id") return;
    pause();
    try {
      ms.forEach(function (m) {
        if (m.type === "characterData") { doText(m.target); }
        else { for (var i = 0; i < m.addedNodes.length; i++) walk(m.addedNodes[i]); }
      });
    } finally { resume(); }
  }
  function readLang() { try { return localStorage.getItem("sg_lang") || "en"; } catch (e) { return "en"; } }

  function exportMissing() {
    var list = Object.keys(missing).sort();
    var blob = new Blob([list.join("\n") + "\n"], { type: "text/plain;charset=utf-8" });
    var a = document.createElement("a");
    a.href = URL.createObjectURL(blob); a.download = "sangerly_untranslated.txt";
    document.body.appendChild(a); a.click(); setTimeout(function () { URL.revokeObjectURL(a.href); a.remove(); }, 500);
    return list.length;
  }
  window.SG_I18N = { missing: function () { return Object.keys(missing).sort(); }, export: exportMissing };

  function start() {
    obs = new MutationObserver(onMutations);
    lang = readLang();
    resume();
    if (lang === "id") apply();
    document.addEventListener("click", function (e) {            // the EN | ID buttons are handled by the app; this runs right after them
      var b = e.target.closest && e.target.closest("#sg_lang button"); if (!b) return;
      setTimeout(function () { lang = readLang(); apply(); }, 30);
    }, true);
    document.addEventListener("keydown", function (e) {          // Ctrl + Alt + Shift + T: save the list of untranslated texts
      if (e.ctrlKey && e.altKey && e.shiftKey && (e.key === "T" || e.key === "t")) { e.preventDefault(); exportMissing(); }
    });
  }
  function load() {
    fetch("i18n_id.json", { cache: "no-store" }).then(function (r) { return r.json(); }).then(function (j) {
      DICT = j.dict || {};
      PATS = (j.patterns || []).map(function (p) { return [new RegExp(p[0]), p[1]]; });
      start();
    }).catch(function () { /* no dictionary: everything stays English */ });
  }
  if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", load); else load();
})();
