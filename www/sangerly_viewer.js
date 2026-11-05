/* ============================================================================
   Sangerly - alignment + chromatogram viewer (canvas, runs in the browser)
   ----------------------------------------------------------------------------
   Mouse : drag = pan | click a read base = select that base | click the consensus
           row = select the whole position | Shift+click = extend to a range |
           Ctrl+wheel = zoom | click/drag the overview strip = jump
   Keys  : Delete / Backspace = delete selection | Ctrl+Z = undo | A C G T N - = set consensus base |
           Left/Right = move | Up/Down = change row   (click the viewer first)
   Talks to Shiny through custom messages  sv_static / sv_dyn / sv_goto / sv_clear
   and the inputs <id>_sel, <id>_del, <id>_undo, <id>_view.
   ========================================================================== */
(function () {
  "use strict";

  var LW = 124;        // width of the fixed label column
  var MINI_H = 34;     // overview strip height
  var MINI_GAP = 6;
  var RULER_H = 18;
  var ROW_H = 22;      // height of one row of base tiles
  var GAP = 8;         // space between reads
  var CHANNELS = ["A", "C", "G", "T"];
  var REG = {};

  function arr(x) { return Array.isArray(x) ? x : (x === null || x === undefined ? [] : [x]); }
  function clamp(v, a, b) { return Math.max(a, Math.min(b, v)); }
  function lowerBound(a, v) {
    var lo = 0, hi = a.length;
    while (lo < hi) { var m = (lo + hi) >> 1; if (a[m] < v) lo = m + 1; else hi = m; }
    return lo;
  }
  function send(id, name, value) {
    if (window.Shiny && window.Shiny.setInputValue)
      window.Shiny.setInputValue(id + "_" + name, value, { priority: "event" });
  }

  function Viewer(id, root) {
    this.id = id; this.root = root;
    this.scroll = root.querySelector(".sv-scroll");
    this.canvas = root.querySelector(".sv-canvas");
    this.spacer = root.querySelector(".sv-spacer");
    this.empty = root.querySelector(".sv-empty");
    this.info = root.querySelector(".sv-info");
    this.zoomEl = root.querySelector(".sv-zoom");
    this.traceEl = root.querySelector(".sv-trace");
    this.thEl = root.querySelector(".sv-th");
    this.ctx = this.canvas.getContext("2d");
    this.cw = parseInt(this.zoomEl.value, 10) || 22;
    this.showTrace = this.traceEl.checked;
    this.bandH = parseInt(this.thEl.value, 10) || 72;
    this.d = null; this.dyn = null;
    this.selCol = 0; this.selTo = 0; this.selRead = 0;
    this.W = 0; this.H = 0; this.dpr = 1;
    this.mini = null; this.miniDirty = true;
    this.drag = null; this.raf = 0; this.viewTimer = 0; this.lastView = "";
    this.pendingGoto = 0; this.savedSl = 0; this.wasVisible = false;
    this.wire();
  }

  Viewer.prototype.wire = function () {
    var self = this, sc = this.scroll;
    sc.addEventListener("scroll", function () {
      if (self.W >= 50) self.savedSl = sc.scrollLeft;
      self.schedule(); self.reportView();
    });
    this.zoomEl.addEventListener("input", function () { self.setZoom(+this.value); });
    this.traceEl.addEventListener("change", function () { self.showTrace = this.checked; self.relayout(); });
    this.thEl.addEventListener("input", function () { self.bandH = +this.value; self.relayout(); });
    if (window.ResizeObserver) new ResizeObserver(function () { self.relayout(); }).observe(sc);
    window.addEventListener("resize", function () { self.relayout(); });
    sc.addEventListener("mousedown", function (e) { self.onDown(e); });
    sc.addEventListener("keydown", function (e) { self.onKey(e); });
    sc.addEventListener("wheel", function (e) {
      if (!(e.ctrlKey || e.metaKey)) return;
      e.preventDefault();
      var r = self.canvas.getBoundingClientRect();
      self.setZoom(self.cw + (e.deltaY < 0 ? 2 : -2), e.clientX - r.left);
    }, { passive: false });
  };

  /* ---- data -------------------------------------------------------------- */
  Viewer.prototype.setStatic = function (m) {
    var d = { nc: +m.nc, nr: +m.nr, hasRef: !!m.has_ref, colors: m.colors || {}, reads: [] };
    d.ref = m.ref || "";
    d.refpos = m.refpos ? arr(m.refpos) : null;
    arr(m.reads).forEach(function (r) {
      var o = { label: String(r.label), base: String(r.base), q: arr(r.q), tr: null, first: 0, last: -1 };
      for (var k = 0; k < o.base.length; k++) if (o.base.charAt(k) !== " ") { o.first = k + 1; break; }
      for (var j = o.base.length - 1; j >= 0; j--) if (o.base.charAt(j) !== " ") { o.last = j + 1; break; }
      if (r.tr && r.tr.x) {
        var x = arr(r.tr.x), t = { x: null };
        var rev = x.length > 1 && x[0] > x[x.length - 1];     // reverse-complemented reads run backwards
        var order = function (a) { a = arr(a); return rev ? a.slice().reverse() : a; };
        t.x = Float32Array.from(order(x));
        CHANNELS.forEach(function (c) { t[c] = Uint8Array.from(order(r.tr[c])); });
        o.tr = t;
      }
      d.reads.push(o);
    });
    this.d = d;
    this.selCol = this.selTo = this.selRead = 0;
    this.empty.style.display = "none";
    this.scroll.style.display = "block";
    this.miniDirty = true;
    this.updateInfo();
    this.relayout();
  };

  Viewer.prototype.setDyn = function (m) {
    this.dyn = { final: String(m.final), review: String(m.review), manual: String(m.manual),
                 delcol: String(m.delcol), mask: arr(m.mask).map(String) };
    this.miniDirty = true;
    this.schedule();
  };

  Viewer.prototype.clear = function () {
    this.d = null; this.dyn = null;
    this.scroll.style.display = "none";
    this.empty.style.display = "block";
  };

  Viewer.prototype.valid = function () {
    var d = this.d, y = this.dyn;
    return d && y && y.final.length === d.nc && y.mask.length === d.nr;
  };

  /* ---- layout ------------------------------------------------------------- */
  Viewer.prototype.relayout = function () {
    var W = this.scroll.clientWidth;
    this.W = W;
    var d = this.d;
    if (!d || W < 50) { this.wasVisible = false; return; }
    var y = MINI_H + MINI_GAP;
    this.yRuler = y; y += RULER_H;
    if (d.hasRef) { this.yRef = y; y += ROW_H; }
    this.yCons = y; y += ROW_H + GAP;
    this.yRead = []; this.hRead = ROW_H + (this.showTrace ? this.bandH : 0);
    for (var i = 0; i < d.nr; i++) { this.yRead.push(y); y += this.hRead + GAP; }
    this.H = y + 2;
    this.dpr = window.devicePixelRatio || 1;
    var pw = Math.round(W * this.dpr), ph = Math.round(this.H * this.dpr);
    if (this.canvas.width !== pw) this.canvas.width = pw;
    if (this.canvas.height !== ph) this.canvas.height = ph;
    this.canvas.style.width = W + "px";
    this.canvas.style.height = this.H + "px";
    this.spacer.style.width = (LW + d.nc * this.cw) + "px";
    this.miniDirty = true;
    if (!this.wasVisible) {                       // just became visible (tab shown / first data)
      this.wasVisible = true;
      if (this.pendingGoto) { var c = this.pendingGoto; this.pendingGoto = 0; this.scrollToCol(c, true); }
      else if (this.savedSl) this.scroll.scrollLeft = this.savedSl;
    }
    this.schedule();
  };

  Viewer.prototype.dataWidth = function () { return Math.max(10, this.W - LW); };

  Viewer.prototype.setZoom = function (v, anchorX) {
    v = clamp(Math.round(v), 4, 48);
    if (v === this.cw) return;
    var Wd = this.dataWidth(), sl = this.scroll.scrollLeft;
    var ax = anchorX === undefined ? LW + Wd / 2 : anchorX;
    var colAt = (sl + (ax - LW)) / this.cw;
    this.cw = v; this.zoomEl.value = v;
    this.relayout();
    this.scroll.scrollLeft = colAt * v - (ax - LW);
  };

  Viewer.prototype.scrollToCol = function (col, force) {
    if (this.W < 50) { this.pendingGoto = col; return; }
    var cw = this.cw, Wd = this.dataWidth(), sl = this.scroll.scrollLeft;
    var left = (col - 1) * cw, right = col * cw;
    if (force || left < sl || right > sl + Wd) this.scroll.scrollLeft = left + cw / 2 - Wd / 2;
  };

  Viewer.prototype.ensureVisible = function (col) {
    var cw = this.cw, Wd = this.dataWidth(), sl = this.scroll.scrollLeft;
    var left = (col - 1) * cw, right = col * cw;
    if (left < sl) this.scroll.scrollLeft = Math.max(0, left - 20);
    else if (right > sl + Wd) this.scroll.scrollLeft = right - Wd + 20;
  };

  Viewer.prototype.reportView = function () {
    var self = this;
    if (this.viewTimer || !this.d || this.W < 50) return;
    this.viewTimer = setTimeout(function () {
      self.viewTimer = 0;
      var cw = self.cw, sl = self.scroll.scrollLeft, Wd = self.dataWidth();
      var first = Math.max(1, Math.floor(sl / cw) + 1), last = Math.min(self.d.nc, Math.ceil((sl + Wd) / cw));
      var key = first + "-" + last;
      if (key !== self.lastView) { self.lastView = key; send(self.id, "view", { first: first, last: last }); }
    }, 200);
  };

  /* ---- drawing ------------------------------------------------------------- */
  Viewer.prototype.schedule = function () {
    var self = this;
    if (!this.raf) this.raf = requestAnimationFrame(function () {
      try { self.draw(); } catch (e) { self.fail(e); }
    });
  };

  // Show the reason in the viewer box instead of staying silently empty.
  Viewer.prototype.fail = function (e) {
    if (window.console) console.error("Sangerly viewer:", e);
    this.scroll.style.display = "none";
    this.empty.style.display = "block";
    this.empty.style.color = "#c0392b";
    this.empty.textContent = "Viewer error: " + (e && e.message ? e.message : e) +
      " (press F12 and open the Console tab for details)";
  };

  function fillTile(ctx, x, y, w, h, color, alpha) {
    ctx.globalAlpha = alpha; ctx.fillStyle = color;
    ctx.fillRect(x + 0.5, y + 1, w - 1, h - 2);
    ctx.globalAlpha = 1;
  }

  Viewer.prototype.buildMini = function (Wd) {
    var d = this.d, dy = this.dyn, nc = d.nc, nr = d.nr, dpr = this.dpr;
    var m = document.createElement("canvas");
    m.width = Math.max(1, Math.ceil(Wd * dpr)); m.height = Math.ceil(MINI_H * dpr);
    var c = m.getContext("2d"); c.scale(dpr, dpr);
    c.fillStyle = "#eef2f6"; c.fillRect(0, 0, Wd, MINI_H);
    var s = Wd / nc, lh = Math.max(1, Math.min(3, 18 / nr));
    c.fillStyle = "#9db4c8";
    for (var i = 0; i < nr; i++) {
      var r = d.reads[i];
      if (r.last >= r.first) c.fillRect((r.first - 1) * s, 3 + i * (lh + 0.5), Math.max(1, (r.last - r.first + 1) * s), lh);
    }
    for (var k = 0; k < nc; k++) {
      var col = dy.delcol.charAt(k) === "1" ? "#2c3e50" : dy.manual.charAt(k) === "1" ? "#27ae60" :
                dy.review.charAt(k) === "1" ? "#e67e22" : null;
      if (col) { c.fillStyle = col; c.fillRect(k * s, MINI_H - 12, Math.max(1, s), 9); }
    }
    this.mini = m; this.miniDirty = false;
  };

  Viewer.prototype.draw = function () {
    this.raf = 0;
    if (!this.valid() || this.W < 50) return;
    var d = this.d, dy = this.dyn, ctx = this.ctx, cw = this.cw, nc = d.nc, nr = d.nr;
    var W = this.W, Wd = this.dataWidth(), sl = this.scroll.scrollLeft;
    var first = Math.max(1, Math.floor(sl / cw) + 1), last = Math.min(nc, Math.ceil((sl + Wd) / cw));
    var col = d.colors || {}, bc = col.base || {}, tc = col.trace || {};
    var showLetters = cw >= 9, fs = Math.min(14, Math.max(8, Math.round(cw * 0.62)));
    var X = function (k) { return LW + (k - 1) * cw - sl; };
    var baseColor = function (ch) { return bc[ch] || "#caa6e8"; };
    var selA = Math.min(this.selCol, this.selTo), selB = Math.max(this.selCol, this.selTo);

    ctx.setTransform(this.dpr, 0, 0, this.dpr, 0, 0);
    ctx.clearRect(0, 0, W, this.H);

    /* overview strip */
    if (this.miniDirty || !this.mini) this.buildMini(Wd);
    ctx.drawImage(this.mini, 0, 0, this.mini.width, this.mini.height, LW, 0, Wd, MINI_H);
    var ms = Wd / nc, vx = LW + (sl / cw) * ms, vw = Math.min(Wd, (Wd / cw) * ms);
    ctx.fillStyle = "rgba(41,128,185,.16)"; ctx.fillRect(vx, 0, Math.max(3, vw), MINI_H);
    ctx.strokeStyle = "#2980b9"; ctx.lineWidth = 1; ctx.strokeRect(vx + .5, .5, Math.max(3, vw) - 1, MINI_H - 1);

    /* data region (clipped so nothing runs under the label column) */
    ctx.save();
    ctx.beginPath(); ctx.rect(LW, MINI_H + MINI_GAP, Wd, this.H); ctx.clip();
    var topCol = this.yRuler + RULER_H - 2;

    for (var k = first; k <= last; k++) {                       // column highlights
      var x = X(k), fill = null;
      if (dy.delcol.charAt(k - 1) === "1") fill = "rgba(44,62,80,.30)";
      else if (dy.manual.charAt(k - 1) === "1") fill = "rgba(39,174,96,.18)";
      else if (dy.review.charAt(k - 1) === "1") fill = "rgba(241,196,15,.28)";
      if (fill) { ctx.fillStyle = fill; ctx.fillRect(x, topCol, cw, this.H - topCol); }
      if (this.selCol && k >= selA && k <= selB) {
        ctx.fillStyle = "rgba(41,128,185,.10)"; ctx.fillRect(x, topCol, cw, this.H - topCol);
      }
    }

    ctx.textAlign = "center"; ctx.textBaseline = "middle";
    ctx.font = "600 " + fs + "px system-ui, Segoe UI, Arial, sans-serif";

    var step = [1, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000].filter(function (s) { return s * cw >= 40; })[0] || 5000;
    ctx.fillStyle = "#555"; ctx.strokeStyle = "#b8c2cc"; ctx.lineWidth = 1;
    ctx.font = "11px system-ui, Segoe UI, Arial, sans-serif";
    for (k = first; k <= last; k++) {                           // ruler
      var lab = d.refpos ? d.refpos[k - 1] : k;
      if (lab > 0 && lab % step === 0) {
        ctx.fillText(String(lab), X(k) + cw / 2, this.yRuler + 7);
        ctx.beginPath(); ctx.moveTo(X(k) + cw / 2 + .5, this.yRuler + 13); ctx.lineTo(X(k) + cw / 2 + .5, this.yRuler + RULER_H - 1); ctx.stroke();
      }
    }
    ctx.font = "600 " + fs + "px system-ui, Segoe UI, Arial, sans-serif";

    if (d.hasRef) {                                             // reference row
      for (k = first; k <= last; k++) {
        var rc = d.ref.charAt(k - 1);
        fillTile(ctx, X(k), this.yRef, cw, ROW_H, baseColor(rc), 1);
        if (showLetters) { ctx.fillStyle = "#1c1c1c"; ctx.fillText(rc, X(k) + cw / 2, this.yRef + ROW_H / 2 + 1); }
      }
    }
    for (k = first; k <= last; k++) {                           // consensus row
      var cc = dy.final.charAt(k - 1), deleted = dy.delcol.charAt(k - 1) === "1";
      fillTile(ctx, X(k), this.yCons, cw, ROW_H, deleted ? "#9aa5b1" : baseColor(cc), 1);
      if (showLetters && !deleted) { ctx.fillStyle = "#1c1c1c"; ctx.fillText(cc, X(k) + cw / 2, this.yCons + ROW_H / 2 + 1); }
      if (cw >= 6) { ctx.strokeStyle = "#2c3e50"; ctx.lineWidth = .6; ctx.strokeRect(X(k) + .5, this.yCons + 1, cw - 1, ROW_H - 2); }
    }

    for (var i = 0; i < nr; i++) {                              // reads
      var rd = d.reads[i], yT = this.yRead[i], mask = dy.mask[i];
      for (k = Math.max(first, rd.first); k <= Math.min(last, rd.last); k++) {
        var ch = rd.base.charAt(k - 1);
        if (ch === " ") continue;
        var xk = X(k);
        if (mask.charAt(k - 1) === "1") {
          fillTile(ctx, xk, yT, cw, ROW_H, "#efefef", 1);
          if (showLetters) { ctx.fillStyle = "#9a9a9a"; ctx.fillText("x", xk + cw / 2, yT + ROW_H / 2 + 1); }
          continue;
        }
        var q = rd.q[k - 1]; if (q === undefined || q < 0) q = 20;
        fillTile(ctx, xk, yT, cw, ROW_H, baseColor(ch), 0.35 + 0.65 * Math.min(q, 40) / 40);
        if (showLetters) { ctx.fillStyle = "#1c1c1c"; ctx.fillText(ch, xk + cw / 2, yT + ROW_H / 2 + 1); }
        if (dy.delcol.charAt(k - 1) !== "1" && ch !== dy.final.charAt(k - 1)) {
          ctx.strokeStyle = "#c0392b"; ctx.lineWidth = 1.6; ctx.strokeRect(xk + 1, yT + 1.5, cw - 2, ROW_H - 3);
        }
      }
      if (this.showTrace && rd.tr) this.drawTrace(ctx, rd.tr, yT + ROW_H + 2, this.bandH - 4, sl, first, last, tc);
      if (this.showTrace) {
        ctx.strokeStyle = "#e3e8ee"; ctx.lineWidth = 1;
        ctx.beginPath(); ctx.moveTo(LW, yT + this.hRead + GAP / 2 + .5); ctx.lineTo(W, yT + this.hRead + GAP / 2 + .5); ctx.stroke();
      }
    }

    if (this.selCol) {                                          // selection outline
      var xa = X(selA), wsel = (selB - selA + 1) * cw;
      ctx.setLineDash([4, 3]); ctx.strokeStyle = "#2980b9"; ctx.lineWidth = 1.2;
      ctx.strokeRect(xa + .5, topCol + .5, wsel - 1, this.H - topCol - 2);
      ctx.setLineDash([]);
      if (this.selRead >= 1 && this.selRead <= nr) {
        ctx.strokeStyle = "#1f5fa8"; ctx.lineWidth = 2.4;
        ctx.strokeRect(xa + 1, this.yRead[this.selRead - 1] + 1, wsel - 2, ROW_H - 2);
      } else {
        ctx.strokeStyle = "#1f5fa8"; ctx.lineWidth = 2.4;
        ctx.strokeRect(xa + 1, this.yCons + 1, wsel - 2, ROW_H - 2);
      }
    }
    ctx.restore();

    /* label column */
    ctx.fillStyle = "#ffffff"; ctx.fillRect(0, MINI_H + MINI_GAP, LW, this.H);
    ctx.strokeStyle = "#d5dbe2"; ctx.lineWidth = 1;
    ctx.beginPath(); ctx.moveTo(LW - .5, MINI_H + MINI_GAP); ctx.lineTo(LW - .5, this.H); ctx.stroke();
    ctx.textAlign = "left"; ctx.textBaseline = "middle";
    ctx.font = "11px system-ui, Segoe UI, Arial, sans-serif"; ctx.fillStyle = "#6b7785";
    ctx.fillText("Overview", 8, MINI_H / 2);
    ctx.fillText(d.refpos ? "Reference pos." : "Column", 8, this.yRuler + 8);
    var fit = function (t) {
      var maxw = LW - 14;
      if (ctx.measureText(t).width <= maxw) return t;
      while (t.length > 1 && ctx.measureText(t + "\u2026").width > maxw) t = t.slice(0, -1);
      return t + "\u2026";
    };
    ctx.font = "12px system-ui, Segoe UI, Arial, sans-serif"; ctx.fillStyle = "#1f3c5f";
    if (d.hasRef) ctx.fillText("Reference", 8, this.yRef + ROW_H / 2);
    ctx.font = "600 12px system-ui, Segoe UI, Arial, sans-serif";
    ctx.fillText("Consensus", 8, this.yCons + ROW_H / 2);
    for (i = 0; i < nr; i++) {
      ctx.font = (this.selRead === i + 1 ? "700 " : "") + "12px system-ui, Segoe UI, Arial, sans-serif";
      ctx.fillStyle = this.selRead === i + 1 ? "#1f5fa8" : "#1f3c5f";
      ctx.fillText(fit(d.reads[i].label), 8, this.yRead[i] + ROW_H / 2);
    }
  };

  Viewer.prototype.drawTrace = function (ctx, tr, yb, bh, sl, first, last, tc) {
    var cw = this.cw, xs = tr.x;
    var i0 = Math.max(0, lowerBound(xs, first - 2) - 1), i1 = Math.min(xs.length - 1, lowerBound(xs, last + 2) + 1);
    if (i1 <= i0) return;
    ctx.lineWidth = 1.1; ctx.lineJoin = "round";
    for (var c = 0; c < CHANNELS.length; c++) {
      var name = CHANNELS[c], v = tr[name];
      ctx.strokeStyle = tc[name] || "#444";
      ctx.beginPath();
      for (var i = i0; i <= i1; i++) {
        var px = LW + (xs[i] - 0.5) * cw - sl, py = yb + bh - (v[i] / 255) * bh;
        if (i === i0) ctx.moveTo(px, py); else ctx.lineTo(px, py);
      }
      ctx.stroke();
    }
  };

  /* ---- mouse ---------------------------------------------------------------- */
  Viewer.prototype.hit = function (x, y) {
    var d = this.d;
    if (y < MINI_H + 2) return { zone: "mini" };
    if (x < LW) return null;
    var col = Math.floor((x - LW + this.scroll.scrollLeft) / this.cw) + 1;
    if (col < 1 || col > d.nc) return null;
    if (y >= this.yRuler && y < this.yCons + ROW_H + GAP / 2) return { zone: "cons", col: col, read: 0 };
    for (var i = 0; i < d.nr; i++)
      if (y >= this.yRead[i] && y < this.yRead[i] + this.hRead + GAP / 2) return { zone: "read", col: col, read: i + 1 };
    return null;
  };

  Viewer.prototype.miniTo = function (x) {
    var d = this.d, Wd = this.dataWidth();
    var col = clamp((x - LW) / (Wd / d.nc), 0, d.nc);
    this.scroll.scrollLeft = col * this.cw - Wd / 2;
  };

  Viewer.prototype.onDown = function (e) {
    if (e.button !== 0 || e.target !== this.canvas || !this.valid()) return;
    var self = this, r = this.canvas.getBoundingClientRect();
    var x = e.clientX - r.left, y = e.clientY - r.top;
    this.scroll.focus({ preventScroll: true });
    var h = this.hit(x, y);
    if (h && h.zone === "mini") { this.drag = { mode: "mini" }; this.miniTo(x); }
    else this.drag = { mode: "pending", x0: e.clientX, sl0: this.scroll.scrollLeft, hit: h, shift: e.shiftKey, moved: false };
    var move = function (ev) {
      var dr = self.drag; if (!dr) return;
      if (dr.mode === "mini") { self.miniTo(ev.clientX - self.canvas.getBoundingClientRect().left); return; }
      var dx = ev.clientX - dr.x0;
      if (!dr.moved && Math.abs(dx) > 4) { dr.moved = true; self.scroll.style.cursor = "grabbing"; }
      if (dr.moved) self.scroll.scrollLeft = dr.sl0 - dx;
    };
    var up = function () {
      window.removeEventListener("mousemove", move); window.removeEventListener("mouseup", up);
      var dr = self.drag; self.drag = null; self.scroll.style.cursor = "";
      if (dr && dr.mode === "pending" && !dr.moved && dr.hit && dr.hit.col) self.select(dr.hit, dr.shift);
    };
    window.addEventListener("mousemove", move); window.addEventListener("mouseup", up);
    e.preventDefault();
  };

  /* ---- selection and keys ------------------------------------------------------- */
  Viewer.prototype.select = function (h, extend) {
    if (extend && this.selCol) this.selTo = h.col;
    else { this.selCol = h.col; this.selTo = h.col; this.selRead = h.read; }
    this.afterSelect();
  };

  Viewer.prototype.afterSelect = function () {
    send(this.id, "sel", { col: this.selCol, from: Math.min(this.selCol, this.selTo),
                           to: Math.max(this.selCol, this.selTo), read: this.selRead });
    this.updateInfo(); this.schedule();
  };

  Viewer.prototype.updateInfo = function () {
    if (!this.info) return;
    var d = this.d;
    if (!d || !this.selCol) { this.info.textContent = "Click a base to select it. Drag to move along the sequence."; return; }
    var a = Math.min(this.selCol, this.selTo), b = Math.max(this.selCol, this.selTo);
    var pos = function (k) { return d.refpos ? "ref " + d.refpos[k - 1] : "col " + k; };
    var where = a === b ? "position " + a + " (" + pos(a) + ")" : "positions " + a + "-" + b + " (" + pos(a) + " to " + pos(b) + ")";
    var who = this.selRead >= 1 ? "read " + d.reads[this.selRead - 1].label : "consensus";
    var act = this.selRead >= 1 ? "Delete = remove " + (a === b ? "this base" : "these bases") + " from this read only"
                                : "Delete = remove " + (a === b ? "this whole position" : "these whole positions") + " from all reads and the consensus";
    this.info.textContent = "Selected: " + who + ", " + where + ". " + act + ". Type A/C/G/T/N to set the consensus base, = for the reference base, 0 to reset. Ctrl+Z undoes.";
    this.info.title = this.info.textContent;     // full text on hover; the line itself never wraps (no layout jump)
  };

  Viewer.prototype.onKey = function (e) {
    var d = this.d; if (!d) return;
    var k = e.key;
    if (k === "Delete" || k === "Backspace") {
      e.preventDefault();
      if (this.selCol) send(this.id, "del", { from: Math.min(this.selCol, this.selTo), to: Math.max(this.selCol, this.selTo), read: this.selRead });
      return;
    }
    if ((e.ctrlKey || e.metaKey) && (k === "z" || k === "Z")) { e.preventDefault(); send(this.id, "undo", Date.now()); return; }
    // Type a base to set the consensus base at the selected position(s), then move on to the next one.
    //   A C G T N (and IUPAC R Y S W K M B D H V) or "-" = that base | "=" "." Enter = reference base | 0 Esc = back to the rule
    if (this.selCol && !e.ctrlKey && !e.metaKey && !e.altKey) {
      var up = k.length === 1 ? k.toUpperCase() : k, tok = null;
      if (k.length === 1 && "ACGTNRYSWKMBDHV-".indexOf(up) >= 0) tok = up;
      else if (k === "=" || k === "." || k === "Enter") tok = "REF";
      else if (k === "0" || k === "Escape") tok = "RULE";
      if (tok) {
        e.preventDefault();
        var sa = Math.min(this.selCol, this.selTo), sb = Math.max(this.selCol, this.selTo);
        send(this.id, "set", { from: sa, to: sb, base: tok, nonce: Date.now() });
        var nx2 = Math.min(d.nc, sb + 1);
        this.selCol = this.selTo = nx2;
        this.ensureVisible(nx2); this.afterSelect();
        return;
      }
    }
    if (k === "ArrowLeft" || k === "ArrowRight") {
      e.preventDefault();
      var cur = this.selCol || Math.max(1, Math.floor(this.scroll.scrollLeft / this.cw) + 1);
      var nx = clamp(cur + (k === "ArrowRight" ? 1 : -1), 1, d.nc);
      if (e.shiftKey && this.selCol) this.selTo = nx; else { this.selCol = this.selTo = nx; }
      this.ensureVisible(nx); this.afterSelect(); return;
    }
    if (k === "ArrowUp" || k === "ArrowDown") {
      e.preventDefault();
      if (!this.selCol) return;
      this.selRead = clamp(this.selRead + (k === "ArrowDown" ? 1 : -1), 0, d.nr);
      this.afterSelect();
    }
  };

  Viewer.prototype.goto = function (col) {
    this.selCol = this.selTo = col; this.selRead = 0;
    this.updateInfo();
    this.scrollToCol(col, false);
    this.schedule();
  };

  /* ---- Shiny wiring ----------------------------------------------------------------- */
  function get(id) {
    if (REG[id]) return REG[id];
    var root = document.getElementById(id);
    if (!root) return null;
    REG[id] = new Viewer(id, root);
    return REG[id];
  }
  function register() {
    var S = window.Shiny;
    S.addCustomMessageHandler("sv_static", function (m) {
      var v = get(m.id); if (!v) return;
      try { v.setStatic(m); } catch (e) { v.fail(e); }
    });
    S.addCustomMessageHandler("sv_dyn", function (m) {
      var v = get(m.id); if (!v) return;
      try { v.setDyn(m); } catch (e) { v.fail(e); }
    });
    S.addCustomMessageHandler("sv_goto", function (m) { var v = get(m.id); if (v && v.d) v.goto(+m.col); });
    S.addCustomMessageHandler("sv_clear", function (m) { var v = get(m.id); if (v) v.clear(); });
  }
  (function wait() {
    if (window.Shiny && window.Shiny.addCustomMessageHandler) register();
    else setTimeout(wait, 30);
  })();

  window.SangerlyViewer = { get: get, _registry: REG };   // exposed for debugging / tests
})();
