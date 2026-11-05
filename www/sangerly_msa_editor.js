/* ============================================================================
   Sangerly - alignment editor (AliView-like), runs in the browser
   ----------------------------------------------------------------------------
   Mouse : click a base = select | drag or Shift+click = rectangle |
           click the ruler / consensus row = select whole columns |
           click a name = select whole sequences (drag = several) | double-click a name = rename |
           drag a SELECTED name up/down = move the sequence(s) | right-click = menu (rename, move, delete, ...) |
           Ctrl+wheel = zoom
   Keys  : letters (A C G T N and IUPAC) or "-" = set the selected bases |
           Space = insert a gap before the selection (shifts the sequence right) |
           Delete = columns: remove them from the result, sequences: delete them, bases: as chosen in "Delete key"
                    (default: remove their column) | Shift+Delete = replace the bases by gaps |
           Backspace = bases: delete and shift left |  Ctrl+Z / Ctrl+Y = undo / redo | arrows = move
           Alt+Up / Alt+Down = move the selected sequences up / down
           "=" = copy the reference sequence's bases onto the selected bases (set the reference first)
   The keys work without clicking the sequence area first (after a click anywhere inside the editor).
   Removed columns are NOT lost: they are only left out of the result and can be restored.
   Talks to Shiny through the custom message  msaed_load  and the input  <id>_state  (Apply edits).
   ========================================================================== */
(function () {
  "use strict";

  var GAP = 45;                       // "-"
  var LW = 170, RULER_H = 20, CONS_H = 22, BAR_H = 34, ROW_H = 20;
  var MAX_UNDO = 40;

  function arr(x) { return Array.isArray(x) ? x : (x === null || x === undefined ? [] : [x]); }
  function clamp(v, a, b) { return Math.max(a, Math.min(b, v)); }
  function u8ToStr(u8) {
    var s = "", CH = 8192;
    for (var i = 0; i < u8.length; i += CH) s += String.fromCharCode.apply(null, u8.subarray(i, i + CH));
    return s;
  }

  /* ======================================================================================
     Model: no DOM, so it can be tested on its own.
     rows[r] is a Uint8Array of character codes; keep[c] = 1 (column in the result) or 0.
     ====================================================================================== */
  function Model(names, seqs, keep) {
    var i;
    this.names = names.slice();
    this.rows = seqs.map(function (s) {
      var a = new Uint8Array(s.length);
      for (var j = 0; j < s.length; j++) a[j] = s.charCodeAt(j);
      return a;
    });
    this.nc = this.rows.length ? this.rows[0].length : 0;
    this.keep = new Uint8Array(this.nc);
    for (i = 0; i < this.nc; i++) this.keep[i] = (typeof keep === "string" && keep.charAt(i) === "0") ? 0 : 1;
    this.undoStack = []; this.redoStack = [];
    this.rev = 0; this.savedRev = 0;
    this._cnt = new Int32Array(128);
  }

  Model.prototype._push = function (op) {
    op.apply();
    this.undoStack.push(op);
    if (this.undoStack.length > MAX_UNDO) this.undoStack.shift();
    this.redoStack = [];
    this.rev++;
  };
  Model.prototype.undo = function () {
    var op = this.undoStack.pop(); if (!op) return false;
    op.revert(); this.redoStack.push(op); this.rev++; return true;
  };
  Model.prototype.redo = function () {
    var op = this.redoStack.pop(); if (!op) return false;
    op.apply(); this.undoStack.push(op); this.rev++; return true;
  };

  // Set the bases of rows r0..r1 at the given (original) columns to the character ch.
  Model.prototype.setCells = function (r0, r1, cols, ch) {
    var code = ch.charCodeAt(0), changes = [], r, k, c, self = this;
    for (r = r0; r <= r1; r++)
      for (k = 0; k < cols.length; k++) {
        c = cols[k];
        if (this.rows[r][c] !== code) changes.push([r, c, this.rows[r][c]]);
      }
    if (!changes.length) return 0;
    this._push({
      apply: function () { for (var i = 0; i < changes.length; i++) self.rows[changes[i][0]][changes[i][1]] = code; },
      revert: function () { for (var i = 0; i < changes.length; i++) self.rows[changes[i][0]][changes[i][1]] = changes[i][2]; }
    });
    return changes.length;
  };

  // Copy the bases of the reference row `ref` onto rows r0..r1 at the given columns. onlyN: only change N (or n)
  // into a real base (A, C, G, T) of the reference.
  Model.prototype.setCellsFromRef = function (r0, r1, cols, ref, onlyN) {
    var self = this, ch = [], refRow = this.rows[ref], r, k, c, nv, ov;
    if (!refRow) return 0;
    for (r = r0; r <= r1; r++) {
      if (r === ref) continue;
      for (k = 0; k < cols.length; k++) {
        c = cols[k]; nv = refRow[c]; ov = this.rows[r][c];
        if (nv === ov) continue;
        if (onlyN && (!(ov === 78 || ov === 110) || !(nv === 65 || nv === 67 || nv === 71 || nv === 84))) continue;
        ch.push([r, c, ov, nv]);
      }
    }
    if (!ch.length) return 0;
    this._push({
      apply: function () { for (var i = 0; i < ch.length; i++) self.rows[ch[i][0]][ch[i][1]] = ch[i][3]; },
      revert: function () { for (var i = 0; i < ch.length; i++) self.rows[ch[i][0]][ch[i][1]] = ch[i][2]; }
    });
    return ch.length;
  };

  // Leave columns out of (val = 0) or put them back into (val = 1) the result.
  Model.prototype.setKeep = function (cols, val) {
    var self = this, changed = [], i;
    for (i = 0; i < cols.length; i++) if (this.keep[cols[i]] !== val) changed.push(cols[i]);
    if (!changed.length) return 0;
    this._push({
      apply: function () { for (var i = 0; i < changed.length; i++) self.keep[changed[i]] = val; },
      revert: function () { for (var i = 0; i < changed.length; i++) self.keep[changed[i]] = 1 - val; }
    });
    return changed.length;
  };

  // Delete the sequences r0..r1.
  Model.prototype.delRows = function (r0, r1) {
    var self = this, n = r1 - r0 + 1, nm, rw;
    this._push({
      apply: function () { nm = self.names.splice(r0, n); rw = self.rows.splice(r0, n); },
      revert: function () {
        Array.prototype.splice.apply(self.names, [r0, 0].concat(nm));
        Array.prototype.splice.apply(self.rows, [r0, 0].concat(rw));
      }
    });
    return n;
  };

  // Move the block of sequences r0..r1 so that it starts at index `to` (an index of the list WITHOUT the block).
  Model.prototype.moveBlock = function (r0, r1, to) {
    var self = this, n = r1 - r0 + 1, oldNames, oldRows;
    this._push({
      apply: function () {
        oldNames = self.names; oldRows = self.rows;
        var idx = [], i; for (i = 0; i < oldNames.length; i++) idx.push(i);
        var blk = idx.splice(r0, n);
        Array.prototype.splice.apply(idx, [to, 0].concat(blk));
        self.names = idx.map(function (k) { return oldNames[k]; });
        self.rows = idx.map(function (k) { return oldRows[k]; });
      },
      revert: function () { self.names = oldNames; self.rows = oldRows; }
    });
    return n;
  };

  Model.prototype.rename = function (r, name) {
    var self = this, old = this.names[r];
    if (old === name) return 0;
    this._push({ apply: function () { self.names[r] = name; }, revert: function () { self.names[r] = old; } });
    return 1;
  };

  // Insert n gaps at column p in the sequences r0..r1 (their bases move right); every other sequence gets
  // n gaps at the end so that all sequences keep the same length. The new end columns are kept.
  Model.prototype.insertGap = function (r0, r1, p, n) {
    var self = this, oldRows, oldKeep, oldNc;
    this._push({
      apply: function () {
        oldRows = self.rows; oldKeep = self.keep; oldNc = self.nc;
        var nr = [], i, a, b;
        for (i = 0; i < oldRows.length; i++) {
          a = oldRows[i]; b = new Uint8Array(oldNc + n);
          if (i >= r0 && i <= r1) { b.set(a.subarray(0, p), 0); b.fill(GAP, p, p + n); b.set(a.subarray(p), p + n); }
          else { b.set(a, 0); b.fill(GAP, oldNc, oldNc + n); }
          nr.push(b);
        }
        var nk = new Uint8Array(oldNc + n); nk.set(oldKeep, 0); nk.fill(1, oldNc);
        self.rows = nr; self.keep = nk; self.nc = oldNc + n;
      },
      revert: function () { self.rows = oldRows; self.keep = oldKeep; self.nc = oldNc; }
    });
    return n;
  };

  // Delete the bases at columns p0..p1 in the sequences r0..r1 and shift the rest left (gaps are added at the end).
  Model.prototype.delShift = function (r0, r1, p0, p1) {
    var self = this, n = p1 - p0 + 1, saved = {};
    this._push({
      apply: function () {
        for (var r = r0; r <= r1; r++) {
          var a = self.rows[r], b = new Uint8Array(self.nc);
          b.set(a.subarray(0, p0), 0); b.set(a.subarray(p1 + 1), p0); b.fill(GAP, self.nc - n, self.nc);
          saved[r] = a; self.rows[r] = b;
        }
      },
      revert: function () { for (var r = r0; r <= r1; r++) self.rows[r] = saved[r]; }
    });
    return n;
  };

  Model.prototype.allGapCols = function () {
    var out = [], c, r, all;
    for (c = 0; c < this.nc; c++) {
      all = true;
      for (r = 0; r < this.rows.length; r++) if (this.rows[r][c] !== GAP) { all = false; break; }
      if (all) out.push(c);
    }
    return out;
  };

  // Most common base of a column (gaps ignored); GAP if the column has no base.
  Model.prototype.consensusCode = function (c) {
    var cnt = this._cnt, r, code, best = 0, bestN = 0, k;
    cnt.fill(0);
    for (r = 0; r < this.rows.length; r++) {
      code = this.rows[r][c];
      if (code !== GAP && code < 128) { k = ++cnt[code]; if (k > bestN) { bestN = k; best = code; } }
    }
    return bestN ? best : GAP;
  };

  Model.prototype.keptCount = function () { var n = 0; for (var i = 0; i < this.nc; i++) n += this.keep[i]; return n; };

  Model.prototype.rowString = function (r, onlyKept) {
    var a = this.rows[r], out;
    if (!onlyKept) return u8ToStr(a);
    out = new Uint8Array(this.keptCount());
    for (var i = 0, j = 0; i < this.nc; i++) if (this.keep[i]) out[j++] = a[i];
    return u8ToStr(out);
  };

  // mode: "all" (every column) | "kept" (columns of the result) | "ungapped" (kept columns, gaps removed)
  Model.prototype.toFasta = function (mode) {
    var out = [], r, s, i;
    for (r = 0; r < this.rows.length; r++) {
      s = this.rowString(r, mode !== "all");
      if (mode === "ungapped") s = s.replace(/-/g, "");
      out.push(">" + this.names[r]);
      for (i = 0; i < s.length; i += 60) out.push(s.substring(i, i + 60));
      if (!s.length) out.push("");
    }
    return out.join("\n") + "\n";
  };

  Model.prototype.state = function () {
    var names = this.names.slice(), seqs = [], keep = "", r, i;
    for (r = 0; r < this.rows.length; r++) seqs.push(u8ToStr(this.rows[r]));
    for (i = 0; i < this.nc; i++) keep += this.keep[i] ? "1" : "0";
    return { names: names, seqs: seqs, keep: keep };
  };

  if (typeof module !== "undefined" && module.exports) module.exports = { Model: Model };
  if (typeof window === "undefined") return;

  /* ======================================================================================
     Editor: canvas drawing + mouse + keyboard
     ====================================================================================== */
  var REG = {};
  var ACTIVE = null;                 // the editor the user touched last (receives the keys even if the focus drifted)
  var FILL = {};
  FILL[65] = "#c6f0d6"; FILL[67] = "#cfe4fb"; FILL[71] = "#e3e3e3"; FILL[84] = "#ffd0d0"; FILL[85] = "#ffd0d0";

  function Editor(id, root) {
    var self = this;
    this.id = id; this.root = root;
    this.scroll = root.querySelector(".msaed-scroll");
    this.canvas = root.querySelector(".msaed-canvas");
    this.spacer = root.querySelector(".msaed-spacer");
    this.empty = root.querySelector(".msaed-empty");
    this.infoEl = root.querySelector(".msaed-info");
    this.ctx = this.canvas.getContext("2d");
    this.cw = 12; this.showRm = true; this.dots = false;
    this.model = null; this.cols = []; this.sel = null; this.drag = null;
    this.msg = ""; this.W = 0; this.H = 0; this.raf = 0;
    this.ro = root.getAttribute("data-readonly") === "1";            // read-only view
    this.refName = null; this.delMode = "col";
    var dm = root.querySelector(".msaed-delmode");
    if (dm) { this.delMode = dm.value || "col"; dm.addEventListener("change", function () { self.delMode = dm.value; self.scroll.focus(); }); }
    root.addEventListener("mousedown", function () { ACTIVE = self; }, true);
    root.addEventListener("focusin", function () { ACTIVE = self; });
    this.barMode = "cons"; this.setHead();
    this.ov = root.querySelector(".msaed-ov"); this.ovBox = root.querySelector(".msaed-ovbox");
    this.ovCache = null; this.statCache = null; this.ovDrag = false;
    var bs = root.querySelector(".msaed-bars");
    if (bs) {
      this.barMode = bs.value || "cons"; this.setHead();
      bs.addEventListener("change", function () { self.barMode = bs.value; self.setHead(); self.layout(); });
    }
    if (this.ov) {
      this.ov.addEventListener("mousedown", function (e) { self.ovDrag = true; self.ovGo(e); e.preventDefault(); });
      window.addEventListener("mousemove", function (e) { if (self.ovDrag) self.ovGo(e); });
      window.addEventListener("mouseup", function () { self.ovDrag = false; });
    }

    var z = root.querySelector(".msaed-zoom");
    if (z) { this.cw = +z.value || 12; z.addEventListener("input", function () { self.cw = +z.value; self.layout(); }); }
    var sr = root.querySelector(".msaed-showrm");
    if (sr) { this.showRm = sr.checked; sr.addEventListener("change", function () { self.showRm = sr.checked; self.sel = null; self.rebuildCols(); self.layout(); }); }
    var dt = root.querySelector(".msaed-dots");
    if (dt) { this.dots = dt.checked; dt.addEventListener("change", function () { self.dots = dt.checked; self.schedule(); }); }
    var btns = root.querySelectorAll("[data-act]");
    for (var i = 0; i < btns.length; i++) (function (b) {
      b.addEventListener("click", function () { self.act(b.getAttribute("data-act")); self.scroll.focus(); });
    })(btns[i]);

    this.scroll.addEventListener("scroll", function () { self.schedule(); });
    this.scroll.addEventListener("keydown", function (e) { e._msaed = 1; self.onKey(e); });
    this.scroll.addEventListener("wheel", function (e) {
      if (!e.ctrlKey) return;
      e.preventDefault();
      self.cw = clamp(self.cw + (e.deltaY < 0 ? 1 : -1), 4, 30);
      if (z) z.value = Math.min(self.cw, +z.max);
      self.layout();
    }, { passive: false });
    this.canvas.addEventListener("mousedown", function (e) { self.onDown(e); });
    this.canvas.addEventListener("dblclick", function (e) { self.onDbl(e); });
    this.canvas.addEventListener("contextmenu", function (e) { self.onMenu(e); });
    window.addEventListener("mousemove", function (e) { if (self.drag) self.onMove(e); });
    window.addEventListener("mouseup", function () { self.onUp(); });
    if (window.ResizeObserver) new ResizeObserver(function () { if (self.model) self.layout(); }).observe(this.scroll);
    else window.addEventListener("resize", function () { if (self.model) self.layout(); });
  }

  Editor.prototype.refIndex = function () {
    return (this.model && this.refName !== null) ? this.model.names.indexOf(this.refName) : -1;
  };

  Editor.prototype.setHead = function () {
    this.bar = this.barMode === "none" ? 0 : BAR_H;
    this.head = this.bar + RULER_H + CONS_H;
  };

  Editor.prototype.clear = function () {
    this.model = null; this.sel = null; this.statCache = null; this.ovCache = null;
    this.scroll.style.display = "none"; if (this.ovBox) this.ovBox.style.display = "none";
    this.empty.style.display = "block"; this.empty.style.color = "";
    if (this.infoEl) this.infoEl.textContent = "";
  };

  Editor.prototype.fail = function (e) {
    if (this.ovBox) this.ovBox.style.display = "none";
    this.empty.style.display = "block"; this.empty.style.color = "#c0392b";
    this.empty.textContent = "The editor failed: " + (e && e.message ? e.message : e);
    this.scroll.style.display = "none";
  };

  Editor.prototype.setData = function (m) {
    var names = arr(m.names).map(String), seqs = arr(m.seqs).map(String);
    var keep = typeof m.keep === "string" ? m.keep : (Array.isArray(m.keep) && m.keep.length ? String(m.keep[0]) : null);
    this.model = new Model(names, seqs, keep);
    this.sel = null; this.msg = "";
    this.statCache = null; this.ovCache = null;
    this.empty.style.display = "none"; this.scroll.style.display = "block";
    if (this.ovBox) this.ovBox.style.display = "block";
    this.rebuildCols(); this.scroll.scrollLeft = 0; this.scroll.scrollTop = 0;
    this.layout();
  };

  Editor.prototype.rebuildCols = function () {
    var m = this.model, out = [];
    for (var c = 0; c < m.nc; c++) if (this.showRm || m.keep[c]) out.push(c);
    this.cols = out;
  };

  Editor.prototype.layout = function () {
    if (!this.model) return;
    var m = this.model, HEAD = this.head;
    this.spacer.style.width = (LW + this.cols.length * this.cw + 20) + "px";
    this.spacer.style.height = (HEAD + m.rows.length * ROW_H + 20) + "px";
    this.sizeCanvas();
    this.sizeOv();
    this.updateInfo();
    this.schedule();
  };
  Editor.prototype.sizeOv = function () {
    if (!this.ov || !this.model) return;
    var w = this.scroll.clientWidth, nr = this.model.rows.length, dpr = window.devicePixelRatio || 1;
    this.ovPH = 44; this.ovMapH = clamp(nr * 3, 10, 90); this.ovMapY = this.ovPH + 6;
    var h = this.ovMapY + this.ovMapH + 4;
    if (w < 2) return;
    if (w === this.ovW && h === this.ovH && dpr === this.ovDpr) return;
    this.ovW = w; this.ovH = h; this.ovDpr = dpr; this.ovCache = null;
    this.ov.width = Math.round(w * dpr); this.ov.height = Math.round(h * dpr);
    this.ov.style.width = w + "px"; this.ov.style.height = h + "px";
  };
  Editor.prototype.sizeCanvas = function () {
    var w = this.scroll.clientWidth, h = this.scroll.clientHeight, dpr = window.devicePixelRatio || 1;
    if (w < 2 || h < 2) return;
    if (w === this.W && h === this.H && dpr === this.dpr) return;
    this.W = w; this.H = h; this.dpr = dpr;
    this.canvas.width = Math.round(w * dpr); this.canvas.height = Math.round(h * dpr);
    this.canvas.style.width = w + "px"; this.canvas.style.height = h + "px";
  };
  Editor.prototype.schedule = function () {
    var self = this;
    if (this.raf) return;
    this.raf = requestAnimationFrame(function () { self.raf = 0; try { self.draw(); } catch (e) { self.fail(e); } });
  };

  /* ---- selection helpers -------------------------------------------------------------- */
  Editor.prototype.norm = function () {
    var s = this.sel; if (!s || !this.model) return null;
    var nr = this.model.rows.length, nc = this.cols.length;
    var r0 = Math.min(s.ar, s.er), r1 = Math.max(s.ar, s.er), c0 = Math.min(s.ac, s.ec), c1 = Math.max(s.ac, s.ec);
    if (s.mode === "cols") { r0 = 0; r1 = nr - 1; }
    if (s.mode === "rows") { c0 = 0; c1 = nc - 1; }
    r0 = clamp(r0, 0, nr - 1); r1 = clamp(r1, 0, nr - 1); c0 = clamp(c0, 0, nc - 1); c1 = clamp(c1, 0, nc - 1);
    if (!nr || !nc) return null;
    return { mode: s.mode, r0: r0, r1: r1, c0: c0, c1: c1 };
  };
  Editor.prototype.selCols = function (n) {
    var out = []; for (var c = n.c0; c <= n.c1; c++) out.push(this.cols[c]); return out;
  };
  Editor.prototype.hit = function (e) {
    var HEAD = this.head;
    var rect = this.canvas.getBoundingClientRect(), x = e.clientX - rect.left, y = e.clientY - rect.top;
    var sx = this.scroll.scrollLeft, sy = this.scroll.scrollTop, nr = this.model.rows.length, nc = this.cols.length;
    var ci = clamp(Math.floor((x - LW + sx) / this.cw), 0, Math.max(0, nc - 1));
    var r = clamp(Math.floor((y - HEAD + sy) / ROW_H), 0, Math.max(0, nr - 1));
    if (y < HEAD && x < LW) return { type: "all" };
    if (y < HEAD) return { type: "col", ci: ci };
    if (x < LW) return { type: "row", r: r };
    return { type: "cell", r: r, ci: ci };
  };

  Editor.prototype.onDown = function (e) {
    if (!this.model || e.button !== 0) return;
    var nr = this.model.rows.length, nc = this.cols.length;
    if (!nr || !nc) return;
    var h = this.hit(e), s = this.sel;
    this.scroll.focus();
    if (h.type === "all") { this.sel = { mode: "cells", ar: 0, er: nr - 1, ac: 0, ec: nc - 1 }; this.drag = null; }
    else if (h.type === "cell") {
      if (e.shiftKey && s && s.mode === "cells") { s.er = h.r; s.ec = h.ci; }
      else this.sel = { mode: "cells", ar: h.r, er: h.r, ac: h.ci, ec: h.ci };
      this.drag = "cells";
    } else if (h.type === "col") {
      if (e.shiftKey && s && s.mode === "cols") s.ec = h.ci;
      else this.sel = { mode: "cols", ar: 0, er: nr - 1, ac: h.ci, ec: h.ci };
      this.drag = "cols";
    } else {
      var n0 = this.norm();
      if (!this.ro && !e.shiftKey && s && s.mode === "rows" && n0 && h.r >= n0.r0 && h.r <= n0.r1) {
        // a name that is already selected: dragging it moves the selected sequences
        this.drag = "rowmove";
        this.mv = { r0: n0.r0, r1: n0.r1, y0: e.clientY, moved: false, at: null, hitr: h.r };
      } else {
        if (e.shiftKey && s && s.mode === "rows") s.er = h.r;
        else this.sel = { mode: "rows", ar: h.r, er: h.r, ac: 0, ec: nc - 1 };
        this.drag = "rows";
      }
    }
    e.preventDefault();
    this.updateInfo(); this.schedule();
  };
  Editor.prototype.onUp = function () {
    var d = this.drag; this.drag = null;
    if (d === "rowmove" && this.mv) {
      var mv = this.mv; this.mv = null;
      if (!mv.moved) { this.sel = { mode: "rows", ar: mv.hitr, er: mv.hitr, ac: 0, ec: Math.max(0, this.cols.length - 1) }; this.updateInfo(); }
      else if (mv.at !== null) this.moveRowsTo(mv.r0, mv.r1, mv.at);
      this.schedule();
    }
  };
  // Put the sequences r0..r1 before the sequence that is now at index `at` (at = number of rows: to the end).
  Editor.prototype.moveRowsTo = function (r0, r1, at) {
    var m = this.model, n = r1 - r0 + 1;
    if (at >= r0 && at <= r1 + 1) return;                       // nothing changes
    var to = at > r1 ? at - n : at;
    m.moveBlock(r0, r1, to);
    this.sel = { mode: "rows", ar: to, er: to + n - 1, ac: 0, ec: Math.max(0, this.cols.length - 1) };
    this.say(n + " sequence(s) moved (Ctrl+Z undoes).");
    this.afterEdit(false);
  };
  Editor.prototype.onMove = function (e) {
    if (this.drag === "rowmove" && this.mv) {
      var mv = this.mv, rect = this.canvas.getBoundingClientRect(), y = e.clientY - rect.top;
      if (!mv.moved && Math.abs(e.clientY - mv.y0) > 4) mv.moved = true;
      if (mv.moved) {
        var nr = this.model.rows.length, sc = this.scroll;
        if (y > this.H - 18) sc.scrollTop += 12; else if (y < this.head + 18) sc.scrollTop -= 12;
        mv.at = clamp(Math.round((y - this.head + sc.scrollTop) / ROW_H), 0, nr);
        this.schedule();
      }
      return;
    }
    var h = this.hit(e), s = this.sel; if (!s) return;
    if (this.drag === "cells") { s.er = h.r !== undefined ? h.r : s.er; s.ec = h.ci !== undefined ? h.ci : s.ec; }
    else if (this.drag === "cols") { if (h.ci !== undefined) s.ec = h.ci; }
    else if (this.drag === "rows") { if (h.r !== undefined) s.er = h.r; }
    this.updateInfo(); this.schedule();
  };
  Editor.prototype.onDbl = function (e) {
    if (!this.model || this.ro) return;
    var h = this.hit(e);
    if (h.type !== "row") return;
    this.renameRow(h.r);
  };
  Editor.prototype.renameRow = function (r) {
    var m = this.model;
    var nm = window.prompt("New name for this sequence:", m.names[r]);
    if (nm === null) return;
    nm = nm.replace(/[\r\n>]/g, " ").trim();
    if (!nm || nm === m.names[r]) return;
    if (m.names.indexOf(nm) >= 0) { this.say("The name '" + nm + "' is already used by another sequence."); return; }
    if (this.refName === m.names[r]) this.refName = nm;
    m.rename(r, nm);
    this.updateInfo(); this.schedule();
  };

  /* ---- keyboard -------------------------------------------------------------------------- */
  Editor.prototype.move = function (dr, dc, extend) {
    var nr = this.model.rows.length, nc = this.cols.length;
    if (!nr || !nc) return;
    if (!this.sel) this.sel = { mode: "cells", ar: 0, er: 0, ac: 0, ec: 0 };
    var s = this.sel;
    if (s.mode === "cols") { s.mode = "cells"; s.ar = 0; s.er = 0; }
    if (s.mode === "rows") { s.mode = "cells"; s.ac = 0; s.ec = 0; }
    s.er = clamp(s.er + dr, 0, nr - 1); s.ec = clamp(s.ec + dc, 0, nc - 1);
    if (!extend) { s.ar = s.er; s.ac = s.ec; }
    this.ensureVisible(s.er, s.ec);
    this.updateInfo(); this.schedule();
  };
  Editor.prototype.ensureVisible = function (r, ci) {
    var sc = this.scroll, x = ci * this.cw, y = r * ROW_H, HEAD = this.head;
    if (x < sc.scrollLeft) sc.scrollLeft = x;
    else if (x + this.cw > sc.scrollLeft + this.W - LW) sc.scrollLeft = x + this.cw - (this.W - LW);
    if (y < sc.scrollTop) sc.scrollTop = y;
    else if (y + ROW_H > sc.scrollTop + this.H - HEAD) sc.scrollTop = y + ROW_H - (this.H - HEAD);
  };
  Editor.prototype.onKey = function (e) {
    if (!this.model) return;
    var k = e.key, handled = true, ctrl = e.ctrlKey || e.metaKey;
    if (this.ro && !(k.indexOf("Arrow") === 0 || k === "Escape" || (ctrl && (k === "a" || k === "A")))) return;
    if (ctrl && (k === "z" || k === "Z")) this.act(e.shiftKey ? "redo" : "undo");
    else if (ctrl && (k === "y" || k === "Y")) this.act("redo");
    else if (ctrl && (k === "a" || k === "A")) {
      var nr = this.model.rows.length, nc = this.cols.length;
      if (nr && nc) { this.sel = { mode: "cells", ar: 0, er: nr - 1, ac: 0, ec: nc - 1 }; this.updateInfo(); this.schedule(); }
    }
    else if (e.altKey && !ctrl && (k === "ArrowUp" || k === "ArrowDown")) this.act(k === "ArrowUp" ? "moveup" : "movedown");
    else if (ctrl || e.altKey) handled = false;
    else if (k === "ArrowLeft") this.move(0, -1, e.shiftKey);
    else if (k === "ArrowRight") this.move(0, 1, e.shiftKey);
    else if (k === "ArrowUp") this.move(-1, 0, e.shiftKey);
    else if (k === "ArrowDown") this.move(1, 0, e.shiftKey);
    else if (k === "Delete") this.act(e.shiftKey ? "gap" : "delete");
    else if (k === "=") this.act("copyref");
    else if (k === "Backspace") this.act("backspace");
    else if (k === " ") this.act("insgap");
    else if (k === "Escape") { this.sel = null; this.updateInfo(); this.schedule(); }
    else if (k.length === 1 && /[A-Za-z\-.]/.test(k)) this.typeBase(k);
    else handled = false;
    if (handled) e.preventDefault();
  };
  Editor.prototype.typeBase = function (k) {
    var n = this.norm(); if (!n) return;
    if (n.mode !== "cells") { this.say("Select bases (not whole columns or sequences) to type a letter."); return; }
    var ch = (k === "." ? "-" : k.toUpperCase());
    var single = n.r0 === n.r1 && n.c0 === n.c1;
    this.model.setCells(n.r0, n.r1, this.selCols(n), ch);
    if (single) this.move(0, 1, false);       // keep typing: the selection moves on
    this.updateInfo(); this.schedule();
  };

  /* ---- actions (buttons and keys) --------------------------------------------------------- */
  Editor.prototype.say = function (t) {
    var self = this; this.msg = t; this.updateInfo();
    clearTimeout(this.msgT); this.msgT = setTimeout(function () { self.msg = ""; self.updateInfo(); }, 5000);
  };
  Editor.prototype.act = function (a) {
    var m = this.model; if (!m) return;
    var n = this.norm(), cols;
    var editing = { copyref: 1, nref: 1, undo: 1, redo: 1, delcols: 1, restore: 1, restoreall: 1, delrows: 1, gap: 1, insgap: 1, shiftdel: 1,
                    "delete": 1, backspace: 1, allgap: 1, commit: 1,
                    rename: 1, moveup: 1, movedown: 1, movetop: 1, movebottom: 1 };
    if (this.ro && editing[a]) { this.say("This view is read-only. Edit in the Editor tab."); return; }
    var needSel = { rename: 1, moveup: 1, movedown: 1, movetop: 1, movebottom: 1, delcols: 1, restore: 1, delrows: 1, gap: 1, insgap: 1, shiftdel: 1, "delete": 1, backspace: 1 };
    if (needSel[a] && !n) { this.say("Select something first (click a base, a column or a sequence name)."); return; }
    if (a === "delete") a = n.mode === "cols" ? "delcols" : (n.mode === "rows" ? "delrows"
      : (this.delMode === "gap" ? "gap" : (this.delMode === "shift" ? "shiftdel" : "delcols")));
    if (a === "backspace") a = n.mode === "cols" ? "delcols" : (n.mode === "rows" ? "delrows" : "shiftdel");

    if (a === "rename") { this.renameRow(n.r0); return; }
    if (a === "moveup" || a === "movedown" || a === "movetop" || a === "movebottom") {
      var nrw = m.rows.length;
      if (a === "moveup") { if (n.r0 > 0) this.moveRowsTo(n.r0, n.r1, n.r0 - 1); }
      else if (a === "movedown") { if (n.r1 < nrw - 1) this.moveRowsTo(n.r0, n.r1, n.r1 + 2); }
      else if (a === "movetop") this.moveRowsTo(n.r0, n.r1, 0);
      else this.moveRowsTo(n.r0, n.r1, nrw);
      return;
    }
    if (a === "undo") { if (!m.undo()) this.say("Nothing to undo."); this.afterEdit(true); return; }
    if (a === "redo") { if (!m.redo()) this.say("Nothing to redo."); this.afterEdit(true); return; }
    if (a === "delcols") {
      cols = this.selCols(n); var k = m.setKeep(cols, 0);
      this.say(k + " column(s) left out of the result (they can be restored)."); this.afterEdit(!this.showRm); return;
    }
    if (a === "restore") {
      cols = this.selCols(n); var k2 = m.setKeep(cols, 1);
      this.say(k2 + " column(s) restored."); this.afterEdit(false); return;
    }
    if (a === "restoreall") {
      var all = []; for (var c = 0; c < m.nc; c++) if (!m.keep[c]) all.push(c);
      this.say(m.setKeep(all, 1) + " column(s) restored."); this.afterEdit(true); return;
    }
    if (a === "delrows") {
      if (n.r1 - n.r0 + 1 >= m.rows.length) { this.say("You cannot delete every sequence."); return; }
      this.say(m.delRows(n.r0, n.r1) + " sequence(s) deleted (Ctrl+Z undoes)."); this.sel = null; this.afterEdit(true); return;
    }
    if (a === "gap") { this.say(m.setCells(n.r0, n.r1, this.selCols(n), "-") + " base(s) set to gap."); this.afterEdit(false); return; }
    if (a === "insgap") {
      var cnt = (n.mode === "rows") ? 1 : (n.c1 - n.c0 + 1);
      m.insertGap(n.r0, n.r1, this.cols[n.c0], cnt);
      this.say(cnt + " gap(s) inserted; the other sequences got gaps at the end."); this.rebuildCols(); this.afterEdit(true); return;
    }
    if (a === "shiftdel") {
      var p0 = this.cols[n.c0], p1 = this.cols[n.c1];
      m.delShift(n.r0, n.r1, p0, p1);
      this.say("Bases deleted and the rest shifted left (gaps added at the end)."); this.afterEdit(false); return;
    }
    if (a === "allgap") {
      var ag = m.allGapCols().filter(function (c) { return m.keep[c]; });
      this.say(m.setKeep(ag, 0) + " all-gap column(s) left out of the result."); this.afterEdit(true); return;
    }
    if (a === "setref") {
      var rn = this.norm(); if (!rn) { this.say("Click a sequence name first, then press Set reference."); return; }
      this.refName = m.names[rn.r0]; this.say("Reference: " + this.refName + ". Select bases and press = (or 'Copy from reference').");
      this.updateInfo(); this.schedule(); return;
    }
    if (a === "clearref") { this.refName = null; this.say("Reference cleared."); this.updateInfo(); this.schedule(); return; }
    if (a === "copyref" || a === "nref") {
      var ref = this.refIndex();
      if (ref < 0) { this.say("Set a reference first: click a sequence name and press 'Set reference'."); return; }
      var rr0, rr1, rcols;
      if (n) { rr0 = n.r0; rr1 = n.r1; rcols = this.selCols(n); }
      else if (a === "nref") { rr0 = 0; rr1 = m.rows.length - 1; rcols = this.cols.slice(); }
      else { this.say("Select the bases to copy first."); return; }
      var cnt2 = m.setCellsFromRef(rr0, rr1, rcols, ref, a === "nref");
      this.say(cnt2 + (a === "nref" ? " N replaced by the reference base" : " base(s) copied from the reference") + (n ? "" : " (whole alignment)") + "."); this.afterEdit(false); return;
    }
    if (a === "commit") { this.commit(); return; }
    if (a === "dl_all") { this.download(m.toFasta("all"), "alignment_edited_all_columns.fasta"); return; }
    if (a === "dl_kept") { this.download(m.toFasta("kept"), "alignment_edited_kept_columns.fasta"); return; }
    if (a === "dl_ungapped") { this.download(m.toFasta("ungapped"), "sequences_edited_ungapped.fasta"); return; }
  };
  Editor.prototype.afterEdit = function (relayout) {
    this.rebuildCols();
    var n = this.cols.length;
    if (this.sel) {      // keep the selection inside the (possibly changed) alignment
      var nr = this.model.rows.length;
      this.sel.ar = clamp(this.sel.ar, 0, nr - 1); this.sel.er = clamp(this.sel.er, 0, nr - 1);
      this.sel.ac = clamp(this.sel.ac, 0, Math.max(0, n - 1)); this.sel.ec = clamp(this.sel.ec, 0, Math.max(0, n - 1));
      if (!n) this.sel = null;
    }
    this.layout();
  };
  Editor.prototype.commit = function () {
    var st = this.model.state();
    if (window.Shiny && window.Shiny.setInputValue) window.Shiny.setInputValue(this.id + "_state", st, { priority: "event" });
    this.model.savedRev = this.model.rev;
    this.say("Edits applied: the other tabs and the downloads now use this alignment.");
  };
  Editor.prototype.download = function (text, name) {
    var blob = new Blob([text], { type: "text/plain" });
    var a = document.createElement("a");
    a.href = URL.createObjectURL(blob); a.download = name;
    document.body.appendChild(a); a.click();
    setTimeout(function () { URL.revokeObjectURL(a.href); a.remove(); }, 500);
  };

  Editor.prototype.updateInfo = function () {
    var m = this.model; if (!m || !this.infoEl) return;
    var nk = m.keptCount(), parts = [];
    if (this.msg) parts.push(this.msg);
    parts.push(m.rows.length + " sequences x " + m.nc + " columns (" + nk + " kept, " + (m.nc - nk) + " removed)");
    var n = this.norm();
    if (n) {
      var a = this.cols[n.c0] + 1, b = this.cols[n.c1] + 1;
      if (n.mode === "rows") parts.push("sequences " + (n.r0 + 1) + "-" + (n.r1 + 1));
      else if (n.mode === "cols") parts.push("columns " + a + "-" + b);
      else parts.push("rows " + (n.r0 + 1) + "-" + (n.r1 + 1) + ", columns " + a + "-" + b);
    }
    if (this.refName !== null) parts.push("reference: " + this.refName);
    if (!this.ro && m.rev !== m.savedRev) parts.push("NOT APPLIED YET: press \"Apply edits\"");
    this.infoEl.textContent = parts.join("   |   ");
  };

  /* ---- statistics (cached per edit) ---------------------------------------------------- */
  Editor.prototype.stats = function () {
    var m = this.model, c = this.statCache;
    if (c && c.rev === m.rev && c.nr === m.rows.length && c.nc === m.nc) return c;
    var nr = m.rows.length, nc = m.nc, occ = new Uint16Array(nc), A = new Uint16Array(nc), C = new Uint16Array(nc),
        G = new Uint16Array(nc), T = new Uint16Array(nc), rowPct = [], r, j, code, row, nb;
    for (r = 0; r < nr; r++) {
      row = m.rows[r]; nb = 0;
      for (j = 0; j < nc; j++) {
        code = row[j];
        if (code !== GAP) {
          occ[j]++; nb++;
          if (code === 65) A[j]++; else if (code === 67) C[j]++; else if (code === 71) G[j]++; else if (code === 84 || code === 85) T[j]++;
        }
      }
      rowPct.push(nc ? Math.round(100 * nb / nc) : 0);
    }
    var cons = new Float32Array(nc), occF = new Float32Array(nc), mx;
    for (j = 0; j < nc; j++) {
      mx = Math.max(A[j], C[j], G[j], T[j]);
      cons[j] = nr ? mx / nr : 0; occF[j] = nr ? occ[j] / nr : 0;
    }
    this.statCache = { rev: m.rev, nr: nr, nc: nc, cons: cons, occ: occF, rowPct: rowPct };
    return this.statCache;
  };

  /* ---- overview: profile of bases per column (mountains = complete, valleys = gaps) + one stripe per sequence ---- */
  Editor.prototype.ovGo = function (e) {
    if (!this.model || !this.ov) return;
    var rect = this.ov.getBoundingClientRect(), x = e.clientX - rect.left, y = e.clientY - rect.top;
    var n = this.cols.length, nr = this.model.rows.length; if (!n || !rect.width) return;
    var ci = clamp(x / rect.width, 0, 1) * n;
    this.scroll.scrollLeft = Math.max(0, ci * this.cw - (this.W - LW) / 2);
    if (y >= this.ovMapY && this.ovMapH) {
      var row = clamp((y - this.ovMapY) / this.ovMapH, 0, 1) * nr;
      this.scroll.scrollTop = Math.max(0, row * ROW_H - (this.H - this.head) / 2);
    }
    this.schedule();
  };
  Editor.prototype.buildOv = function (key) {
    var m = this.model, st = this.stats(), W = this.ovW, n = this.cols.length, nr = m.rows.length, cols = this.cols;
    var prof = new Float32Array(W), rem = new Float32Array(W), lo = new Int32Array(W), hi = new Int32Array(W);
    var px, ci, r, a, b, s, rm, col;
    for (px = 0; px < W; px++) {
      a = Math.floor(px * n / W); b = Math.max(a + 1, Math.floor((px + 1) * n / W));
      if (a > n - 1) a = n - 1; if (b > n) b = n;
      lo[px] = a; hi[px] = b; s = 0; rm = 0;
      for (ci = a; ci < b; ci++) { col = cols[ci]; s += st.occ[col]; if (!m.keep[col]) rm++; }
      prof[px] = b > a ? s / (b - a) : 0; rem[px] = b > a ? rm / (b - a) : 0;
    }
    var mc = this.ovMap || (this.ovMap = document.createElement("canvas"));
    mc.width = W; mc.height = Math.max(1, nr);
    var mctx = mc.getContext("2d"), img = mctx.createImageData(W, Math.max(1, nr)), d = img.data, row, cnt, tot, p, f, k;
    for (r = 0; r < nr; r++) {
      row = m.rows[r];
      for (px = 0; px < W; px++) {
        cnt = 0; tot = hi[px] - lo[px];
        for (ci = lo[px]; ci < hi[px]; ci++) if (row[cols[ci]] !== GAP) cnt++;
        p = tot ? cnt / tot : 0; f = rem[px] > 0.5 ? 0.35 : 1; k = (r * W + px) * 4;
        d[k] = 255 - (255 - 44) * p * f; d[k + 1] = 255 - (255 - 123) * p * f; d[k + 2] = 255 - (255 - 229) * p * f; d[k + 3] = 255;
      }
    }
    mctx.putImageData(img, 0, 0);
    this.ovCache = { key: key, prof: prof, rem: rem };
  };
  Editor.prototype.drawOv = function () {
    if (!this.ov || !this.model || !this.ovW) return;
    var m = this.model, n = this.cols.length; if (!n) return;
    var key = m.rev + "|" + this.ovW + "|" + n + "|" + m.rows.length + "|" + (this.showRm ? 1 : 0);
    if (!this.ovCache || this.ovCache.key !== key) this.buildOv(key);
    var ctx = this.ov.getContext("2d"), W = this.ovW, H = this.ovH, PH = this.ovPH, dpr = this.ovDpr || 1, c = this.ovCache, px, y;
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    ctx.fillStyle = "#f5f7fa"; ctx.fillRect(0, 0, W, H);
    ctx.beginPath(); ctx.moveTo(0, PH);
    for (px = 0; px < W; px++) ctx.lineTo(px + 0.5, PH - 2 - c.prof[px] * (PH - 10));
    ctx.lineTo(W, PH); ctx.closePath(); ctx.fillStyle = "#a9c3e4"; ctx.fill();
    ctx.strokeStyle = "#5f86b8"; ctx.lineWidth = 1; ctx.beginPath();
    for (px = 0; px < W; px++) { y = PH - 2 - c.prof[px] * (PH - 10); if (px === 0) ctx.moveTo(0.5, y); else ctx.lineTo(px + 0.5, y); }
    ctx.stroke();
    ctx.fillStyle = "#e67e22";
    for (px = 0; px < W; px++) if (c.rem[px] > 0.5) ctx.fillRect(px, 0, 1, 3);
    ctx.imageSmoothingEnabled = m.rows.length * 3 > this.ovMapH;
    ctx.drawImage(this.ovMap, 0, this.ovMapY, W, this.ovMapH);
    ctx.strokeStyle = "#dfe3e8"; ctx.lineWidth = 1; ctx.strokeRect(0.5, this.ovMapY - 0.5, W - 1, this.ovMapH + 1);
    var sx = this.scroll.scrollLeft, vis = Math.max(1, (this.W - LW) / this.cw), x0 = sx / this.cw / n * W, wv = Math.max(4, vis / n * W);
    ctx.fillStyle = "rgba(44,123,229,.12)"; ctx.fillRect(x0, 0, wv, H);
    ctx.strokeStyle = "#c0392b"; ctx.lineWidth = 1.5; ctx.strokeRect(x0 + 0.5, 0.5, wv - 1, H - 1);
    ctx.fillStyle = "#52606d"; ctx.font = "11px Consolas, monospace"; ctx.textAlign = "left"; ctx.textBaseline = "top";
    ctx.fillText("bases per column (top = every sequence has a base)   |   one stripe per sequence (white = gaps)", 6, 3);
  };


  /* ---- right-click menu ------------------------------------------------------------------------- */
  var MENU = null;
  function hideMenu() { if (MENU) MENU.style.display = "none"; }
  function ensureMenu() {
    if (MENU) return MENU;
    var st = document.createElement("style");
    st.textContent = ".msaed-menu{display:none;position:fixed;z-index:10000;background:#fff;border:1px solid #cfd6df;border-radius:8px;" +
      "box-shadow:0 6px 20px rgba(0,0,0,.22);padding:4px;min-width:210px;font:13px system-ui,Segoe UI,Arial,sans-serif}" +
      ".msaed-menu .mi{padding:6px 12px;cursor:pointer;border-radius:5px;color:#1f3c5f}.msaed-menu .mi:hover{background:#eaf2ff}" +
      ".msaed-menu .mi.danger{color:#c0392b}.msaed-menu .sep{height:1px;background:#e3e8ee;margin:4px 2px}";
    document.head.appendChild(st);
    MENU = document.createElement("div"); MENU.className = "msaed-menu";
    document.body.appendChild(MENU);
    document.addEventListener("mousedown", function (e) { if (MENU && !MENU.contains(e.target)) hideMenu(); }, true);
    window.addEventListener("scroll", hideMenu, true);
    window.addEventListener("resize", hideMenu);
    window.addEventListener("keydown", function (e) { if (e.key === "Escape") hideMenu(); }, true);
    return MENU;
  }
  var MENUS = {
    row: [["rename", "Rename..."], ["setref", "Set as reference"], ["-"],
          ["moveup", "Move up   (Alt+Up)"], ["movedown", "Move down   (Alt+Down)"], ["movetop", "Move to top"], ["movebottom", "Move to bottom"], ["-"],
          ["delrows", "Delete sequence(s)", 1]],
    col: [["delcols", "Remove column(s) from the result"], ["restore", "Restore column(s)"], ["-"],
          ["insgap", "Insert gap before the column(s)"], ["allgap", "Remove all-gap columns"]],
    cell: [["gap", "Set to gap (-)"], ["insgap", "Insert gap"], ["shiftdel", "Delete + shift left"], ["copyref", "Copy from reference"], ["-"],
           ["delcols", "Remove these columns from the result"]]
  };
  Editor.prototype.onMenu = function (e) {
    if (!this.model || this.ro) return;
    var nr = this.model.rows.length, nc = this.cols.length; if (!nr || !nc) return;
    var h = this.hit(e); if (h.type === "all") return;
    e.preventDefault(); ACTIVE = this; this.scroll.focus();
    var n = this.norm(), inside = false;
    if (h.type === "row") inside = n && n.mode !== "cols" && h.r >= n.r0 && h.r <= n.r1;
    else if (h.type === "col") inside = n && n.mode === "cols" && h.ci >= n.c0 && h.ci <= n.c1;
    else inside = n && n.mode === "cells" && h.r >= n.r0 && h.r <= n.r1 && h.ci >= n.c0 && h.ci <= n.c1;
    if (!inside) {
      if (h.type === "row") this.sel = { mode: "rows", ar: h.r, er: h.r, ac: 0, ec: nc - 1 };
      else if (h.type === "col") this.sel = { mode: "cols", ar: 0, er: nr - 1, ac: h.ci, ec: h.ci };
      else this.sel = { mode: "cells", ar: h.r, er: h.r, ac: h.ci, ec: h.ci };
      this.updateInfo(); this.schedule();
    }
    var self = this, menu = ensureMenu(), defs = MENUS[h.type];
    menu.innerHTML = "";
    defs.forEach(function (d) {
      var el = document.createElement("div");
      if (d[0] === "-") { el.className = "sep"; menu.appendChild(el); return; }
      el.className = "mi" + (d[2] ? " danger" : ""); el.textContent = d[1];
      el.addEventListener("click", function () { hideMenu(); self.act(d[0]); self.scroll.focus(); });
      menu.appendChild(el);
    });
    menu.style.display = "block";
    menu.style.left = Math.max(2, Math.min(e.clientX, window.innerWidth - menu.offsetWidth - 4)) + "px";
    menu.style.top = Math.max(2, Math.min(e.clientY, window.innerHeight - menu.offsetHeight - 4)) + "px";
  };

  /* ---- drawing ------------------------------------------------------------------------------ */
  Editor.prototype.draw = function () {
    var m = this.model; if (!m) return;
    this.sizeCanvas();
    var W = this.W, H = this.H; if (!W || !H) return;
    var HEAD = this.head, oy = this.bar, st = this.stats();
    var ctx = this.ctx, cw = this.cw, cols = this.cols, nr = m.rows.length;
    var sx = this.scroll.scrollLeft, sy = this.scroll.scrollTop, dpr = this.dpr || 1;
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    ctx.clearRect(0, 0, W, H);
    ctx.font = "12px Consolas, 'Courier New', monospace";
    ctx.textAlign = "center"; ctx.textBaseline = "middle";

    var c0 = Math.max(0, Math.floor(sx / cw)), c1 = Math.min(cols.length - 1, Math.floor((sx + W - LW) / cw));
    var r0 = Math.max(0, Math.floor(sy / ROW_H)), r1 = Math.min(nr - 1, Math.floor((sy + H - HEAD) / ROW_H));
    var ci, r, x, y, code, col, txt = cw >= 8;
    var cons = [], refIdx = this.refIndex(), refRow = refIdx >= 0 ? m.rows[refIdx] : null;
    for (ci = c0; ci <= c1; ci++) cons.push(m.consensusCode(cols[ci]));

    // sequence cells
    ctx.save(); ctx.beginPath(); ctx.rect(LW, HEAD, W - LW, H - HEAD); ctx.clip();
    for (r = r0; r <= r1; r++) {
      y = HEAD + r * ROW_H - sy;
      for (ci = c0; ci <= c1; ci++) {
        col = cols[ci]; code = m.rows[r][col]; x = LW + ci * cw - sx;
        ctx.globalAlpha = m.keep[col] ? 1 : 0.3;
        var same = this.dots && code !== GAP && r !== refIdx && code === (refRow ? refRow[col] : cons[ci - c0]);
        if (code !== GAP && !same) { ctx.fillStyle = FILL[code] || "#f4f0c8"; ctx.fillRect(x, y, cw, ROW_H - 1); }
        if (txt) {
          if (code === GAP) { ctx.fillStyle = "#9aa5b1"; ctx.fillText("-", x + cw / 2, y + ROW_H / 2); }
          else if (same) { ctx.fillStyle = "#9aa5b1"; ctx.fillText(".", x + cw / 2, y + ROW_H / 2); }
          else { ctx.fillStyle = "#1f2d3d"; ctx.fillText(String.fromCharCode(code), x + cw / 2, y + ROW_H / 2); }
        } else if (code === GAP) { ctx.fillStyle = "#e6e9ee"; ctx.fillRect(x, y + ROW_H / 2 - 1, cw, 2); }
      }
    }
    ctx.globalAlpha = 1; ctx.restore();

    // names panel
    ctx.fillStyle = "#fff"; ctx.fillRect(0, HEAD, LW, H - HEAD);
    var n = this.norm();
    ctx.textAlign = "left";
    for (r = r0; r <= r1; r++) {
      y = HEAD + r * ROW_H - sy;
      if (n && r >= n.r0 && r <= n.r1 && n.mode !== "cols") { ctx.fillStyle = "#dbe8fb"; ctx.fillRect(0, y, LW, ROW_H - 1); }
      var nm = m.names[r], maxc = Math.floor((LW - 14 - 30) / 7.2);
      if (nm.length > maxc) nm = nm.substring(0, maxc - 1) + "\u2026";
      ctx.fillStyle = r === refIdx ? "#a04000" : "#1f3c5f"; ctx.fillText(nm, 8, y + ROW_H / 2);
      ctx.textAlign = "right";
      if (r === refIdx) { ctx.fillStyle = "#a04000"; ctx.fillText("REF", LW - 6, y + ROW_H / 2); }
      else { ctx.fillStyle = "#9aa5b1"; ctx.fillText(st.rowPct[r] + "%", LW - 6, y + ROW_H / 2); }
      ctx.textAlign = "left";
    }
    ctx.strokeStyle = "#dfe3e8"; ctx.beginPath(); ctx.moveTo(LW - 0.5, 0); ctx.lineTo(LW - 0.5, H); ctx.stroke();

    // header: ruler + consensus
    ctx.fillStyle = "#fff"; ctx.fillRect(0, 0, W, HEAD);
    ctx.save(); ctx.beginPath(); ctx.rect(LW, 0, W - LW, HEAD); ctx.clip();
    var step = cw >= 14 ? 5 : (cw >= 8 ? 10 : 50);
    ctx.textAlign = "center";
    if (this.bar) {
      var bsrc = this.barMode === "occ" ? st.occ : st.cons, bh;
      ctx.fillStyle = "#f5f7fa"; ctx.fillRect(LW, 0, W - LW, this.bar);
      for (ci = c0; ci <= c1; ci++) {
        col = cols[ci]; x = LW + ci * cw - sx; bh = Math.round(bsrc[col] * (this.bar - 6));
        ctx.globalAlpha = m.keep[col] ? 1 : 0.3; ctx.fillStyle = "#8d99a6";
        ctx.fillRect(x, this.bar - 3 - bh, Math.max(1, cw - 1), bh);
      }
      ctx.globalAlpha = 1;
    }
    for (ci = c0; ci <= c1; ci++) {
      col = cols[ci]; x = LW + ci * cw - sx;
      if (!m.keep[col]) { ctx.fillStyle = "#e67e22"; ctx.fillRect(x, oy + RULER_H - 4, cw, 3); }
      if ((col + 1) % step === 0) {
        ctx.fillStyle = "#7b8794"; ctx.fillRect(x + cw / 2 - 0.5, oy + RULER_H - 7, 1, 4);
        ctx.fillStyle = "#52606d"; ctx.fillText(String(col + 1), x + cw / 2, oy + 8);
      }
      code = cons[ci - c0]; ctx.globalAlpha = m.keep[col] ? 1 : 0.3;
      if (code !== GAP) { ctx.fillStyle = FILL[code] || "#f4f0c8"; ctx.fillRect(x, oy + RULER_H, cw, CONS_H - 2); }
      if (txt) { ctx.fillStyle = code === GAP ? "#9aa5b1" : "#1f2d3d"; ctx.fillText(String.fromCharCode(code), x + cw / 2, oy + RULER_H + CONS_H / 2 - 1); }
      ctx.globalAlpha = 1;
    }
    ctx.restore();
    ctx.textAlign = "left"; ctx.fillStyle = "#1f3c5f"; ctx.font = "bold 12px Consolas, 'Courier New', monospace";
    ctx.fillText("Consensus", 8, oy + RULER_H + CONS_H / 2 - 1);
    ctx.fillStyle = "#7b8794"; ctx.font = "11px Consolas, monospace"; ctx.fillText("column", 8, oy + 8);
    if (this.bar) { ctx.fillText(this.barMode === "occ" ? "bases (not gaps)" : "conservation", 8, 14); }
    ctx.strokeStyle = "#dfe3e8"; ctx.beginPath(); ctx.moveTo(0, HEAD - 0.5); ctx.lineTo(W, HEAD - 0.5); ctx.stroke();

    // selection
    if (n) {
      ctx.save(); ctx.beginPath(); ctx.rect(LW, HEAD, W - LW, H - HEAD); ctx.clip();
      var sx0 = LW + n.c0 * cw - sx, sw = (n.c1 - n.c0 + 1) * cw, sy0 = HEAD + n.r0 * ROW_H - sy, sh = (n.r1 - n.r0 + 1) * ROW_H;
      ctx.fillStyle = "rgba(44,123,229,.16)"; ctx.fillRect(sx0, sy0, sw, sh);
      ctx.strokeStyle = "#2c7be5"; ctx.lineWidth = 1.5; ctx.strokeRect(sx0 + 0.5, sy0 + 0.5, sw - 1, sh - 1);
      ctx.restore();
      if (n.mode !== "rows") { ctx.fillStyle = "rgba(44,123,229,.22)"; ctx.fillRect(Math.max(LW, sx0), 0, Math.max(0, Math.min(sw, W - Math.max(LW, sx0))), HEAD - 1); }
    }
    // drop position while moving sequences
    if (this.drag === "rowmove" && this.mv && this.mv.moved && this.mv.at !== null) {
      var ly = HEAD + this.mv.at * ROW_H - sy;
      ctx.fillStyle = "#e67e22"; ctx.fillRect(0, ly - 1.5, W, 3);
    }
    this.drawOv();
  };

  /* ---- Shiny wiring --------------------------------------------------------------------------- */
  // Keys keep working when the focus has drifted away from the sequence area (after a click on a button, a tab...).
  window.addEventListener("keydown", function (e) {
    if (e._msaed || !ACTIVE || !ACTIVE.model) return;
    var a = ACTIVE; if (!a.root.offsetParent) return;                         // its tab is hidden
    var t = e.target, tag = (t && t.tagName ? t.tagName.toLowerCase() : "");
    var typing = tag === "textarea" || tag === "select" || (t && t.isContentEditable) ||
                 (tag === "input" && !/^(range|checkbox|button)$/.test(t.type || ""));
    if (typing) return;
    if (!(a.root.contains(t) || t === document.body || t === document.documentElement)) return;   // focus is on something else
    a.onKey(e);
  });

  function get(id) {
    if (REG[id]) return REG[id];
    var root = document.getElementById(id);
    if (!root) return null;
    REG[id] = new Editor(id, root);
    return REG[id];
  }
  function register() {
    window.Shiny.addCustomMessageHandler("msaed_clear", function (m) {
      var e = get(m.id); if (e) e.clear();
    });
    window.Shiny.addCustomMessageHandler("msaed_load", function (m) {
      var e = get(m.id); if (!e) return;
      try { e.setData(m); } catch (err) { e.fail(err); }
    });
  }
  (function wait() {
    if (window.Shiny && window.Shiny.addCustomMessageHandler) register();
    else setTimeout(wait, 30);
  })();

  window.SangerlyMSA = { get: get, Model: Model, _registry: REG };
})();
