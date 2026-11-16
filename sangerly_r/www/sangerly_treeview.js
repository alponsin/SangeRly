/* SangeRly - interactive tree viewer for the Phylogeny tab (plain SVG, no library).
   Layouts: rectangular, circular, radial (unrooted). Support as coloured dots or numbers, tip colours by the first word of the name,
   search, collapse a clade (click its node), flip a clade (double-click), zoom (wheel), pan (drag), export SVG / PNG. */
(function () {
  'use strict';
  var NS = 'http://www.w3.org/2000/svg';
  var PAL = ['#2563eb', '#16a34a', '#dc2626', '#d97706', '#7c3aed', '#0891b2', '#db2777', '#65a30d', '#ea580c', '#475569', '#0d9488', '#9333ea'];

  // ---- Newick parser (quoted names, node labels, branch lengths) -------------------------------------------------
  function parseNewick(txt) {
    var s = String(txt || '').trim(), i = 0;
    function name() {
      var out = '';
      if (s[i] === "'") { i++; while (i < s.length && !(s[i] === "'" && s[i + 1] !== "'")) { if (s[i] === "'" && s[i + 1] === "'") { out += "'"; i += 2; } else out += s[i++]; } i++; return out; }
      while (i < s.length && ':,();'.indexOf(s[i]) < 0) out += s[i++];
      return out.trim().replace(/_/g, ' ');
    }
    function node() {
      var n = { ch: [], name: '', len: 0, label: '', collapsed: false };
      if (s[i] === '(') {
        i++;
        do { if (s[i] === ',') i++; var c = node(); c.parent = n; n.ch.push(c); } while (s[i] === ',');
        i++;                                               // ')'
        n.label = name();
      } else n.name = name();
      if (s[i] === ':') { i++; var t = ''; while (i < s.length && ',();'.indexOf(s[i]) < 0) t += s[i++]; n.len = parseFloat(t) || 0; }
      return n;
    }
    var r = node();
    return r;
  }
  function leaves(n, out) { out = out || []; if (!n.ch.length || n.collapsed) out.push(n); else n.ch.forEach(function (c) { leaves(c, out); }); return out; }
  function allLeaves(n, out) { out = out || []; if (!n.ch.length) out.push(n); else n.ch.forEach(function (c) { allLeaves(c, out); }); return out; }
  function supportOf(n) {                                  // "SH-aLRT/UFBoot": the last number counts
    if (!n.label) return NaN;
    var p = n.label.split('/'); return parseFloat(p[p.length - 1]);
  }
  function genus(nm) { return (nm.split(/[ _|]/)[0] || nm); }

  // ---- layouts: every node gets x, y (cartesian) -------------------------------------------------------------------------
  function depthOf(n, d, useLen, W) { n.d = d; W.max = Math.max(W.max, d); n.ch.forEach(function (c) { depthOf(c, d + (useLen ? Math.max(c.len, 0) : 1), useLen, W); }); }
  function layoutRect(root, useLen, align) {
    var W = { max: 0 }; depthOf(root, 0, useLen, W);
    var k = 0;
    (function y(n) {
      if (!n.ch.length || n.collapsed) { n.y = k++; return; }
      n.ch.forEach(y); n.y = (n.ch[0].y + n.ch[n.ch.length - 1].y) / 2;
    })(root);
    if (!useLen) { // cladogram: tips on the right edge
      (function lvl(n) { if (!n.ch.length || n.collapsed) { n.x = W.max; return; } n.ch.forEach(lvl); })(root);
      (function h(n) { if (!n.ch.length || n.collapsed) return; n.ch.forEach(h); n.x = Math.min.apply(null, n.ch.map(function (c) { return c.x; })) - 1; })(root);
      return { nl: k, w: W.max || 1, kind: 'rect', align: true };
    }
    (function xs(n) { n.x = n.d; n.ch.forEach(xs); })(root);
    return { nl: k, w: W.max || 1, kind: 'rect', align: align };
  }
  function layoutCirc(root, useLen) {
    var r = layoutRect(root, useLen, false), nl = r.nl;
    (function c(n) { var a = (n.y / nl) * 2 * Math.PI; n.a = a; n.r = n.x; n.cx = n.r * Math.cos(a); n.cy = n.r * Math.sin(a); n.ch.forEach(c); })(root);
    return { nl: nl, w: r.w, kind: 'circ' };
  }
  function layoutRadial(root, useLen) {          // equal-angle algorithm
    function cnt(n) { n.nl = (!n.ch.length || n.collapsed) ? 1 : n.ch.reduce(function (a, c) { return a + cnt(c); }, 0); return n.nl; }
    var tot = cnt(root), W = { max: 0 };
    root.cx = 0; root.cy = 0;
    (function go(n, a0, a1) {
      var a = a0;
      if (!n.ch.length || n.collapsed) return;
      n.ch.forEach(function (c) {
        var span = (a1 - a0) * c.nl / n.nl, mid = a + span / 2, L = useLen ? Math.max(c.len, 0) : 1;
        c.cx = n.cx + L * Math.cos(mid); c.cy = n.cy + L * Math.sin(mid); c.ang = mid;
        W.max = Math.max(W.max, Math.hypot(c.cx, c.cy));
        go(c, a, a + span); a += span;
      });
    })(root, 0, 2 * Math.PI);
    return { nl: tot, w: W.max || 1, kind: 'radial' };
  }

  // ---- viewer ---------------------------------------------------------------------------------------------------------------
  function Viewer(el) {
    var me = this; me.el = el; me.trees = {}; me.cur = null; me.root = null;
    me.o = { layout: 'rect', len: true, labels: true, size: 12, sup: 'dot', minsup: 0, hi: 95, mid: 70, color: false, align: false, q: '' };
    el.innerHTML =
      '<div class="sgt-bar">' +
      '<label>Tree <select class="sgt-tree"></select></label>' +
      '<label>Layout <select class="sgt-layout"><option value="rect">Rectangular</option><option value="circ">Circular</option><option value="radial">Radial (unrooted)</option></select></label>' +
      '<label><input type="checkbox" class="sgt-len" checked> Branch lengths</label>' +
      '<label><input type="checkbox" class="sgt-align"> Align tips</label>' +
      '<label><input type="checkbox" class="sgt-labels" checked> Names</label>' +
      '<label>Size <input type="range" class="sgt-size" min="6" max="26" value="12"></label>' +
      '<label>Support <select class="sgt-sup"><option value="dot">Dots</option><option value="num">Numbers</option><option value="both">Both</option><option value="none">None</option></select></label>' +
      '<label><input type="checkbox" class="sgt-color"> Colour by first word</label>' +
      '<input type="search" class="sgt-q" placeholder="Find a name...">' +
      '<span class="sgt-btns"><button type="button" class="sgt-fit">Fit</button><button type="button" class="sgt-exp">Expand all</button>' +
      '<button type="button" class="sgt-ladd">Ladderize</button><button type="button" class="sgt-svg">SVG</button><button type="button" class="sgt-png">PNG</button></span>' +
      '</div><div class="sgt-legend"></div><div class="sgt-stage"><svg class="sgt-svg-el" xmlns="' + NS + '"></svg><div class="sgt-tip"></div>' +
      '<div class="sgt-help">Click an inner node: collapse / expand. Double-click: flip the two sides. Wheel: zoom. Drag: move.</div></div>';
    me.svg = el.querySelector('.sgt-svg-el'); me.tip = el.querySelector('.sgt-tip');
    me.vb = { x: 0, y: 0, w: 1000, h: 600 };
    var q = function (c) { return el.querySelector(c); };
    q('.sgt-tree').onchange = function () { me.setTree(this.value); };
    q('.sgt-layout').onchange = function () { me.o.layout = this.value; me.draw(true); };
    q('.sgt-len').onchange = function () { me.o.len = this.checked; me.draw(true); };
    q('.sgt-align').onchange = function () { me.o.align = this.checked; me.draw(); };
    q('.sgt-labels').onchange = function () { me.o.labels = this.checked; me.draw(); };
    q('.sgt-size').oninput = function () { me.o.size = +this.value; me.draw(); };
    q('.sgt-sup').onchange = function () { me.o.sup = this.value; me.draw(); };
    q('.sgt-color').onchange = function () { me.o.color = this.checked; me.draw(); };
    q('.sgt-q').oninput = function () { me.o.q = this.value.toLowerCase(); me.draw(); };
    q('.sgt-fit').onclick = function () { me.draw(true); };
    q('.sgt-exp').onclick = function () { if (me.root) { (function e(n) { n.collapsed = false; n.ch.forEach(e); })(me.root); me.draw(true); } };
    q('.sgt-ladd').onclick = function () { if (me.root) { me.ladderize(); me.draw(); } };
    q('.sgt-svg').onclick = function () { me.exportSVG(); };
    q('.sgt-png').onclick = function () { me.exportPNG(); };
    // zoom and pan
    me.svg.addEventListener('wheel', function (e) {
      e.preventDefault();
      var r = me.svg.getBoundingClientRect(), f = e.deltaY < 0 ? 0.87 : 1.15;
      var mx = me.vb.x + (e.clientX - r.left) / r.width * me.vb.w, my = me.vb.y + (e.clientY - r.top) / r.height * me.vb.h;
      me.vb.x = mx - (mx - me.vb.x) * f; me.vb.y = my - (my - me.vb.y) * f; me.vb.w *= f; me.vb.h *= f; me.applyVB();
    }, { passive: false });
    var drag = null;
    me.svg.addEventListener('mousedown', function (e) { if (e.target.closest && e.target.closest('.sgt-node')) return; drag = { x: e.clientX, y: e.clientY, vx: me.vb.x, vy: me.vb.y }; me.svg.style.cursor = 'grabbing'; });
    window.addEventListener('mousemove', function (e) {
      if (!drag) return; var r = me.svg.getBoundingClientRect();
      me.vb.x = drag.vx - (e.clientX - drag.x) / r.width * me.vb.w; me.vb.y = drag.vy - (e.clientY - drag.y) / r.height * me.vb.h; me.applyVB();
    });
    window.addEventListener('mouseup', function () { drag = null; me.svg.style.cursor = ''; });
  }
  Viewer.prototype.applyVB = function () { this.svg.setAttribute('viewBox', [this.vb.x, this.vb.y, this.vb.w, this.vb.h].join(' ')); };
  Viewer.prototype.setData = function (trees) {
    var me = this, sel = me.el.querySelector('.sgt-tree'), names = Object.keys(trees || {}), keep = me.cur;
    me.trees = trees || {}; me.parsed = {};
    sel.innerHTML = names.map(function (n) { return '<option>' + n + '</option>'; }).join('');
    if (!names.length) { me.root = null; me.svg.innerHTML = ''; me.el.querySelector('.sgt-legend').textContent = 'Build a tree (NJ or IQ-TREE) first; it appears here.'; return; }
    me.setTree(names.indexOf(keep) >= 0 ? keep : names[0]);
  };
  Viewer.prototype.setTree = function (nm) {
    this.cur = nm; this.el.querySelector('.sgt-tree').value = nm;
    try { this.root = parseNewick(this.trees[nm]); } catch (e) { this.root = null; }
    this.draw(true);
  };
  Viewer.prototype.ladderize = function () {
    (function sz(n) { n.size = n.ch.length ? n.ch.reduce(function (a, c) { return a + sz(c); }, 0) : 1; return n.size; })(this.root);
    (function ord(n) { n.ch.sort(function (a, b) { return a.size - b.size; }); n.ch.forEach(ord); })(this.root);
  };
  function el(tag, attrs, parent) { var e = document.createElementNS(NS, tag); for (var k in attrs) e.setAttribute(k, attrs[k]); if (parent) parent.appendChild(e); return e; }
  Viewer.prototype.draw = function (fit) {
    var me = this, root = me.root, o = me.o; if (!root) return;
    var svg = me.svg; while (svg.firstChild) svg.removeChild(svg.firstChild);
    var useLen = o.len && (function a(n) { return n.len > 0 || n.ch.some(a); })(root);
    var L = o.layout === 'circ' ? layoutCirc(root, useLen) : o.layout === 'radial' ? layoutRadial(root, useLen) : layoutRect(root, useLen, o.align);
    var tips = allLeaves(root), groups = {};
    tips.forEach(function (t) { groups[genus(t.name)] = 0; });
    var gnames = Object.keys(groups).sort(); gnames.forEach(function (g, i) { groups[g] = PAL[i % PAL.length]; });
    var leg = me.el.querySelector('.sgt-legend');
    leg.innerHTML = o.color && gnames.length > 1 && gnames.length <= 40 ? gnames.map(function (g) { return '<span><i style="background:' + groups[g] + '"></i>' + g + '</span>'; }).join('') : (tips.length + ' tips');
    var step = 18 * (o.size / 12), sx, W, H;
    var g = el('g', {}, svg), fs = o.size, tipFont = 'font-family="Segoe UI,Arial,sans-serif"';
    function hit(t) { return o.q && t.name.toLowerCase().indexOf(o.q) >= 0; }
    function col(t) { return o.color ? groups[genus(t.name)] : '#0f172a'; }
    function P(n) { // pixel position
      if (L.kind === 'rect') return [n.x * sx, n.y * step];
      return [n.cx * sx, n.cy * sx];
    }
    if (L.kind === 'rect') { W = 760; sx = W / (L.w || 1); H = (L.nl + 1) * step; }
    else { sx = 300 / (L.w || 1); W = 700; H = 700; }
    var lines = el('g', { stroke: '#475569', 'stroke-width': 1.5, fill: 'none', 'stroke-linecap': 'round', 'stroke-linejoin': 'round' }, g);
    var dots = el('g', {}, g), labs = el('g', {}, g);
    var lblX = L.kind === 'rect' && L.align ? W + 8 : null;
    function branch(n) {
      n.ch.forEach(function (c) {
        var p = P(n), q = P(c), cp;
        if (L.kind === 'rect') {
          el('path', { d: 'M' + p[0] + ',' + p[1] + 'V' + q[1] + 'H' + q[0] }, lines);
        } else if (L.kind === 'circ') {
          var r1 = n.r * sx, a = c.a, sweep = (c.a > n.a) ? 1 : 0, ax = n.cx * sx, ay = n.cy * sx;
          var mx = r1 * Math.cos(a), my = r1 * Math.sin(a);
          el('path', { d: 'M' + ax + ',' + ay + 'A' + r1 + ',' + r1 + ' 0 0 ' + sweep + ' ' + mx + ',' + my + 'L' + q[0] + ',' + q[1] }, lines);
        } else el('path', { d: 'M' + p[0] + ',' + p[1] + 'L' + q[0] + ',' + q[1] }, lines);
        if (c.ch.length && !c.collapsed) branch(c);
      });
    }
    if (root.ch.length && !root.collapsed) branch(root);
    // nodes, supports, tip labels
    (function nodes(n) {
      var p = P(n), isTip = !n.ch.length;
      if (!isTip) {
        var s = supportOf(n), show = !isNaN(s) && s >= o.minsup && o.sup !== 'none';
        var grp = el('g', { 'class': 'sgt-node', style: 'cursor:pointer' }, dots);
        if (n.collapsed) {
          var lv = allLeaves(n), w = 26, hh = Math.min(26, step * 0.9);
          el('path', { d: 'M' + p[0] + ',' + p[1] + 'l' + w + ',' + (-hh / 2) + 'v' + hh + 'z', fill: '#94a3b8', stroke: '#475569' }, grp);
          var tx = el('text', { x: p[0] + w + 6, y: p[1] + fs / 3, 'font-size': fs, fill: '#334155' }, grp); tx.textContent = lv.length + ' tips: ' + lv[0].name + ' ...';
        } else if (show && (o.sup === 'dot' || o.sup === 'both')) {
          el('circle', { cx: p[0], cy: p[1], r: Math.max(3, fs * 0.33), fill: s >= o.hi ? '#16a34a' : s >= o.mid ? '#d97706' : '#dc2626', stroke: '#fff', 'stroke-width': 1 }, grp);
        } else el('circle', { cx: p[0], cy: p[1], r: 3, fill: '#fff', stroke: '#64748b', 'stroke-width': 1.2 }, grp);
        if (show && (o.sup === 'num' || o.sup === 'both') && !n.collapsed) {
          var t2 = el('text', { x: p[0] - 4, y: p[1] - 5, 'text-anchor': 'end', 'font-size': Math.max(7, fs * 0.72), fill: '#b91c1c' }, labs); t2.textContent = n.label;
        }
        var ti = document.createElementNS(NS, 'title'); ti.textContent = (n.label ? 'support ' + n.label + '\n' : '') + allLeaves(n).length + ' tips' + (n.len ? '\nbranch ' + n.len : ''); grp.appendChild(ti);
        grp.addEventListener('click', function (e) { e.stopPropagation(); n.collapsed = !n.collapsed; me.draw(); });
        grp.addEventListener('dblclick', function (e) { e.stopPropagation(); n.collapsed = false; n.ch.reverse(); me.draw(); });
      }
      if (isTip && o.labels) {
        var x = lblX !== null ? lblX : p[0], deg = 0, anc = 'start', dx = 7;
        if (lblX !== null) el('path', { d: 'M' + p[0] + ',' + p[1] + 'H' + lblX, stroke: '#cbd5e1', 'stroke-dasharray': '2 3', fill: 'none' }, lines);
        var T = el('text', { x: x + dx, y: p[1] + fs / 3, 'font-size': fs, fill: col(n), 'font-weight': hit(n) ? 700 : 400 }, labs);
        if (L.kind !== 'rect') {
          var a = L.kind === 'circ' ? n.a : (n.ang || 0), flip = Math.cos(a) < 0;
          deg = a * 180 / Math.PI + (flip ? 180 : 0); anc = flip ? 'end' : 'start';
          T.setAttribute('x', (flip ? -dx : dx)); T.setAttribute('y', fs / 3); T.setAttribute('text-anchor', anc);
          T.setAttribute('transform', 'translate(' + p[0] + ',' + p[1] + ') rotate(' + deg + ')');
        }
        T.textContent = n.name;
        if (hit(n)) { var b = el('circle', { cx: p[0], cy: p[1], r: fs * 0.7, fill: 'none', stroke: '#f59e0b', 'stroke-width': 2.5 }, dots); }
      } else if (isTip) el('circle', { cx: p[0], cy: p[1], r: 2.5, fill: col(n) }, dots);
      if (!n.collapsed) n.ch.forEach(nodes);
    })(root);
    // scale bar (rectangular and circular with real branch lengths)
    if (useLen && L.kind !== 'radial' || (useLen && L.kind === 'radial')) {
      var raw = L.w / 5, mag = Math.pow(10, Math.floor(Math.log10(raw))), unit = [1, 2, 5, 10].map(function (m) { return m * mag; }).filter(function (v) { return v <= raw * 1.5; }).pop() || mag;
      var bx = L.kind === 'rect' ? 0 : -W / 2 + 20, by = L.kind === 'rect' ? H + 14 : H / 2 - 20;
      el('path', { d: 'M' + bx + ',' + by + 'h' + unit * sx, stroke: '#0f172a', 'stroke-width': 2 }, g);
      var st = el('text', { x: bx, y: by + 16, 'font-size': 11, fill: '#334155' }, g); st.textContent = unit.toPrecision(2) + ' substitutions/site';
    }
    if (fit) {
      var bb = g.getBBox(), pad = 30;
      me.vb = { x: bb.x - pad, y: bb.y - pad, w: bb.width + 2 * pad, h: bb.height + 2 * pad };
      var stage = me.el.querySelector('.sgt-stage'); stage.style.height = Math.max(420, Math.min(900, me.vb.h * (stage.clientWidth || 800) / me.vb.w)) + 'px';
      // keep the aspect of the stage so that the tree fills it
      var r = me.svg.getBoundingClientRect(), ar = (r.width || 800) / (r.height || 600), cur = me.vb.w / me.vb.h;
      if (cur < ar) { var nw = me.vb.h * ar; me.vb.x -= (nw - me.vb.w) / 2; me.vb.w = nw; } else { var nh = me.vb.w / ar; me.vb.y -= (nh - me.vb.h) / 2; me.vb.h = nh; }
    }
    me.applyVB();
  };
  Viewer.prototype.svgText = function () {
    var c = this.svg.cloneNode(true); c.setAttribute('xmlns', NS);
    var v = this.vb; c.setAttribute('viewBox', [v.x, v.y, v.w, v.h].join(' ')); c.setAttribute('width', 1600); c.setAttribute('height', Math.round(1600 * v.h / v.w));
    c.insertAdjacentHTML('afterbegin', '<rect x="' + v.x + '" y="' + v.y + '" width="' + v.w + '" height="' + v.h + '" fill="#fff"/>');
    return '<?xml version="1.0" encoding="UTF-8"?>\n' + new XMLSerializer().serializeToString(c);
  };
  function save(blob, name) { var a = document.createElement('a'); a.href = URL.createObjectURL(blob); a.download = name; document.body.appendChild(a); a.click(); setTimeout(function () { URL.revokeObjectURL(a.href); a.remove(); }, 1500); }
  Viewer.prototype.exportSVG = function () { save(new Blob([this.svgText()], { type: 'image/svg+xml' }), 'tree_' + (this.cur || 'tree') + '.svg'); };
  Viewer.prototype.exportPNG = function () {
    var v = this.vb, W = 2400, H = Math.round(2400 * v.h / v.w), img = new Image(), me = this;
    img.onload = function () { var c = document.createElement('canvas'); c.width = W; c.height = H; var x = c.getContext('2d'); x.fillStyle = '#fff'; x.fillRect(0, 0, W, H); x.drawImage(img, 0, 0, W, H); c.toBlob(function (b) { save(b, 'tree_' + (me.cur || 'tree') + '.png'); }); };
    img.src = 'data:image/svg+xml;charset=utf-8,' + encodeURIComponent(this.svgText());
  };

  // ---- hook to Shiny ----------------------------------------------------------------------------------------------------------------
  var viewers = {}, pending = {};
  function ensure(id) {
    if (viewers[id]) return viewers[id];
    var e = document.getElementById(id); if (!e) return null;
    viewers[id] = new Viewer(e); return viewers[id];
  }
  (function wait() {
    if (!(window.Shiny && Shiny.addCustomMessageHandler)) return setTimeout(wait, 40);
    Shiny.addCustomMessageHandler('sg_tree_data', function (m) {
      var v = ensure(m.id); if (v) v.setData(m.trees || {}); else pending[m.id] = m.trees;
    });
    document.addEventListener('shown.bs.tab', function () { Object.keys(viewers).forEach(function (k) { if (viewers[k].root && viewers[k].el.offsetParent) viewers[k].draw(true); }); });
  })();
  document.addEventListener('DOMContentLoaded', function () { Object.keys(pending).forEach(function (id) { var v = ensure(id); if (v) v.setData(pending[id]); }); });
})();
