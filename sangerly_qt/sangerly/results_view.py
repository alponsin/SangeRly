"""Document tab that shows the result of an assembly: summary tiles, sample list, colour-coded consensus, variants and notes."""
from pathlib import Path

from PySide6.QtCore import Qt
from PySide6.QtGui import QColor
from PySide6.QtWidgets import (QAbstractItemView, QApplication, QFileDialog, QFrame, QHBoxLayout, QHeaderView, QLabel, QPlainTextEdit,
                               QPushButton, QSplitter, QTableWidget, QTableWidgetItem, QTabWidget, QVBoxLayout, QWidget)

from . import theme
from .widgets import Chip, SequenceView, StatTile, base_legend, panel_head

OK_COL, BAD_COL = QColor(theme.OK), QColor(theme.BAD)


def fasta_wrap(name, seq, n=70):
    return f">{name}\n" + "\n".join(seq[i:i + n] for i in range(0, len(seq), n)) + "\n"


def _num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return 0.0


class ResultsView(QWidget):
    def __init__(self, title, request, result, parent=None):
        super().__init__(parent)
        self.title = title
        self.samples = result.get("samples", [])
        self.saved = result.get("saved", [])
        ok_s = [s for s in self.samples if s.get("ok")]
        v = QVBoxLayout(self)
        v.setContentsMargins(18, 16, 18, 14)
        v.setSpacing(12)

        head = QHBoxLayout()
        head.setSpacing(10)
        t = QLabel(title)
        t.setObjectName("h2")
        head.addWidget(t)
        head.addWidget(Chip("reference-guided" if request.get("mode") == "reference" else "de novo", "info"))
        head.addStretch(1)
        b1 = QPushButton("Save all consensus FASTA…")
        b1.clicked.connect(self._save_all)
        head.addWidget(b1)
        v.addLayout(head)

        tiles = QHBoxLayout()
        tiles.setSpacing(10)
        mean_len = (sum(_num(s.get("length")) for s in ok_s) / len(ok_s)) if ok_s else 0
        review = int(sum(_num(s.get("review")) for s in ok_s))
        for val, lab, col in [(len(self.samples), "Samples", theme.NAVY), (len(ok_s), "Assembled", theme.OK),
                              (len(self.samples) - len(ok_s), "Failed", theme.BAD if len(self.samples) > len(ok_s) else theme.MUTED),
                              (f"{mean_len:.0f} bp", "Mean length", theme.NAVY), (review, "Positions to review", theme.WARN if review else theme.MUTED)]:
            tiles.addWidget(StatTile(val, lab, col))
        v.addLayout(tiles)
        if self.saved:
            s = QLabel(f"{len(self.saved)} file(s) saved in {Path(self.saved[0]).parent}")
            s.setObjectName("hint")
            v.addWidget(s)

        sp = QSplitter(Qt.Horizontal)
        v.addWidget(sp, 1)
        left = QFrame()
        left.setObjectName("panel")
        ll = QVBoxLayout(left)
        ll.setContentsMargins(0, 0, 0, 0)
        ll.setSpacing(0)
        ll.addWidget(panel_head("Samples", Chip(f"{len(self.samples)}", "muted")))
        self.tbl = QTableWidget(len(self.samples), 5)
        self.tbl.setHorizontalHeaderLabels(["Sample", "Status", "Length", "N / IUPAC", "To review"])
        self.tbl.verticalHeader().setVisible(False)
        self.tbl.setAlternatingRowColors(True)
        self.tbl.setShowGrid(False)
        self.tbl.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.tbl.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.tbl.horizontalHeader().setSectionResizeMode(0, QHeaderView.Stretch)
        for c in (1, 2, 3, 4):
            self.tbl.horizontalHeader().setSectionResizeMode(c, QHeaderView.ResizeToContents)
        for i, s in enumerate(self.samples):
            ok = s.get("ok")
            vals = [s["sample"], "● Assembled" if ok else "● Failed", s.get("length", ""), s.get("ambiguous", ""), s.get("review", "")]
            for j, val in enumerate(vals):
                it = QTableWidgetItem(str(val))
                if j == 1:
                    it.setForeground(OK_COL if ok else BAD_COL)
                if j in (2, 3, 4):
                    it.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)
                self.tbl.setItem(i, j, it)
        self.tbl.itemSelectionChanged.connect(self._select)
        self.tbl.setFrameShape(QFrame.NoFrame)
        ll.addWidget(self.tbl, 1)
        sp.addWidget(left)

        right = QTabWidget()
        right.setObjectName("inner")
        w = QWidget()
        wl = QVBoxLayout(w)
        wl.setContentsMargins(12, 10, 12, 10)
        wl.setSpacing(8)
        bar = QHBoxLayout()
        bar.setSpacing(8)
        self.chips = QHBoxLayout()
        self.chips.setSpacing(6)
        bar.addLayout(self.chips)
        bar.addStretch(1)
        bc = QPushButton("Copy")
        bc.clicked.connect(self._copy)
        bs = QPushButton("Save FASTA…")
        bs.clicked.connect(self._save_one)
        bar.addWidget(bc)
        bar.addWidget(bs)
        wl.addLayout(bar)
        self.cons = SequenceView()
        wl.addWidget(self.cons, 1)
        wl.addWidget(base_legend())
        right.addTab(w, "Consensus")
        self.var = QTableWidget(0, 6)
        self.var.setHorizontalHeaderLabels(["Ref pos", "Reference", "Consensus", "Kind", "Supported", "Needs review"])
        self.var.verticalHeader().setVisible(False)
        self.var.setAlternatingRowColors(True)
        self.var.setShowGrid(False)
        self.var.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.var.horizontalHeader().setSectionResizeMode(3, QHeaderView.Stretch)
        right.addTab(self.var, "Variants")
        self.notes = QPlainTextEdit()
        self.notes.setReadOnly(True)
        right.addTab(self.notes, "Notes")
        sp.addWidget(right)
        sp.setSizes([420, 780])
        sp.setChildrenCollapsible(False)
        if self.samples:
            self.tbl.selectRow(0)

    def current(self):
        r = self.tbl.currentRow()
        return self.samples[r] if 0 <= r < len(self.samples) else None

    def _set_chips(self, items):
        while self.chips.count():
            it = self.chips.takeAt(0)
            if it.widget():
                it.widget().deleteLater()
        for text, kind in items:
            self.chips.addWidget(Chip(text, kind))

    def _select(self):
        s = self.current()
        if not s:
            return
        if not s.get("ok"):
            self.cons.set_sequence("")
            self.var.setRowCount(0)
            self.notes.setPlainText(s.get("msg", ""))
            self._set_chips([("failed", "bad")])
            return
        seq = s["consensus"]
        self.cons.set_sequence(seq)
        acgt = sum(seq.upper().count(b) for b in "ACGT")
        gc = (seq.upper().count("G") + seq.upper().count("C")) / acgt * 100 if acgt else 0
        amb = len(seq) - acgt
        self._set_chips([(f"{len(seq)} bp", "info"), (f"GC {gc:.1f}%", "muted"), (f"{amb} N / IUPAC", "warn" if amb else "ok"),
                         (f"{int(_num(s.get('review')))} to review", "warn" if _num(s.get("review")) else "ok")])
        vs = s.get("variants") or []
        self.var.setRowCount(len(vs))
        for i, x in enumerate(vs):
            for j, k in enumerate(["RefPos", "Reference", "Consensus", "Kind", "Supported", "Needs_review"]):
                val = x.get(k)
                self.var.setItem(i, j, QTableWidgetItem("" if val is None else ("yes" if val is True else "" if val is False else str(val))))
        notes = s.get("notes") or []
        if isinstance(notes, str):
            notes = [notes]
        self.notes.setPlainText("\n".join(notes + ["", "Reads: " + ", ".join(s.get("reads") or [])]))

    def _copy(self):
        s = self.current()
        if s and s.get("ok"):
            QApplication.clipboard().setText(fasta_wrap(s["sample"], s["consensus"]))

    def _save_one(self):
        s = self.current()
        if not s or not s.get("ok"):
            return
        p, _ = QFileDialog.getSaveFileName(self, "Save consensus", f"{s['sample']}_consensus.fasta", "FASTA (*.fasta *.fa)")
        if p:
            Path(p).write_text(fasta_wrap(s["sample"], s["consensus"]), encoding="utf-8")

    def _save_all(self):
        ok = [s for s in self.samples if s.get("ok")]
        if not ok:
            return
        p, _ = QFileDialog.getSaveFileName(self, "Save all consensus sequences", "all_consensus.fasta", "FASTA (*.fasta *.fa)")
        if p:
            Path(p).write_text("".join(fasta_wrap(s["sample"], s["consensus"]) for s in ok), encoding="utf-8")
