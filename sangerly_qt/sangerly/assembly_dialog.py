"""Assembly dialog (Base / Advanced tabs, numbers instead of sliders) with selectable mode cards and drag-and-drop of .ab1 files."""
from pathlib import Path

from PySide6.QtCore import QSettings, Qt
from PySide6.QtGui import QColor
from PySide6.QtGui import QGuiApplication
from PySide6.QtWidgets import (QScrollArea, QAbstractItemView, QButtonGroup, QCheckBox, QComboBox, QDialog, QDoubleSpinBox, QFileDialog, QFormLayout,
                               QFrame, QGridLayout, QGroupBox, QHBoxLayout, QHeaderView, QLabel, QLineEdit, QMessageBox, QPushButton,
                               QRadioButton, QSpinBox, QTableWidget, QTableWidgetItem, QTabWidget, QVBoxLayout, QWidget)

from . import theme
from .widgets import DropFrame

NAMING = [   # (key used by the R engine, text)
    ("name", "By file name  (the F / R tag is removed: S1_F + S1_R = sample S1)"),
    ("name_norun", "By file name, ignoring long run numbers  (1st BASE / Macrogen IDs)"),
    ("all", "All files are one sample"),
    ("one_per_file", "One sample per file  (every file is its own sample)"),
    ("custom", "One sample with a name I type"),
]


def spin(lo, hi, val, suffix="", step=1):
    s = QSpinBox()
    s.setRange(lo, hi)
    s.setValue(val)
    s.setSingleStep(step)
    if suffix:
        s.setSuffix(suffix)
    s.setAlignment(Qt.AlignRight)
    s.setMinimumWidth(110)
    s.setButtonSymbols(QSpinBox.UpDownArrows)
    return s


def dspin(lo, hi, val, dec=2, step=0.01, suffix=""):
    s = QDoubleSpinBox()
    s.setRange(lo, hi)
    s.setDecimals(dec)
    s.setSingleStep(step)
    s.setValue(val)
    if suffix:
        s.setSuffix(suffix)
    s.setAlignment(Qt.AlignRight)
    s.setMinimumWidth(110)
    return s


def expand_paths(paths, recursive=True):
    """Dropped / chosen paths -> list of .ab1 files (folders are searched)."""
    out = []
    for p in paths:
        q = Path(p)
        if q.is_dir():
            it = q.rglob("*") if recursive else q.glob("*")
            out += sorted(str(x) for x in it if x.suffix.lower() == ".ab1")
        elif q.suffix.lower() == ".ab1":
            out.append(str(q))
    return out


class OptionCard(QFrame):
    """A selectable card wrapping a radio button (the card is highlighted when selected)."""

    def __init__(self, radio, title, hint):
        super().__init__()
        self.setObjectName("optCard")
        self.radio = radio
        self.setCursor(Qt.PointingHandCursor)
        v = QVBoxLayout(self)
        v.setContentsMargins(14, 12, 14, 12)
        v.setSpacing(3)
        radio.setText(title)
        radio.setObjectName("optRadio")
        h = QLabel(hint)
        h.setObjectName("hint")
        h.setWordWrap(True)
        h.setContentsMargins(26, 0, 0, 0)
        v.addWidget(radio)
        v.addWidget(h)
        radio.toggled.connect(self._sync)
        self._sync(radio.isChecked())

    def _sync(self, on):
        self.setProperty("selected", bool(on))
        theme.repolish(self)

    def mousePressEvent(self, e):
        self.radio.setChecked(True)


def path_row(edit: QLineEdit, button_text="Browse…"):
    w = QWidget()
    h = QHBoxLayout(w)
    h.setContentsMargins(0, 0, 0, 0)
    h.setSpacing(6)
    b = QPushButton(button_text)
    h.addWidget(edit, 1)
    h.addWidget(b)
    return w, b


def _form_group(title):
    g = QGroupBox(title)
    f = QFormLayout(g)
    f.setFieldGrowthPolicy(QFormLayout.FieldsStayAtSizeHint)
    f.setHorizontalSpacing(14)
    f.setVerticalSpacing(8)
    return g, f


DIR_TEXT = {"F": ("Forward  →", "#1f6feb"), "R": ("←  Reverse", "#8b5cf6")}


class AssemblyDialog(QDialog):
    def __init__(self, bridge, parent=None, files=None):
        super().__init__(parent)
        self.bridge = bridge
        self.cfg = QSettings()
        self.setWindowTitle("Sanger Assembly")
        self.setMinimumSize(640, 420)
        scr = (parent.screen() if parent is not None else None) or QGuiApplication.primaryScreen()
        av = scr.availableGeometry() if scr else None
        self.resize(min(1000, av.width() - 60) if av else 1000, min(760, av.height() - 80) if av else 760)
        self.rows = []                       # [{name, path, sample, direction}]
        self.result_request = None

        root = QVBoxLayout(self)
        root.setContentsMargins(0, 0, 0, 0)
        root.setSpacing(0)

        head = QFrame()
        head.setObjectName("dlgHead")
        hl = QHBoxLayout(head)
        hl.setContentsMargins(22, 14, 22, 14)
        hl.setSpacing(14)
        ic = QLabel()
        ic.setPixmap(theme.icon("assembly", 44, "#2f6fed").pixmap(44, 44))
        hl.addWidget(ic)
        tv = QVBoxLayout()
        tv.setSpacing(1)
        t = QLabel("Sanger assembly")
        t.setObjectName("h2")
        s = QLabel("Trim, pair and assemble .ab1 reads into one consensus per sample.")
        s.setObjectName("hint")
        tv.addWidget(t)
        tv.addWidget(s)
        hl.addLayout(tv, 1)
        root.addWidget(head)

        tabs = QTabWidget()
        tabs.addTab(self._scrolled(self._base_tab()), "Base")
        tabs.addTab(self._scrolled(self._advanced_tab()), "Advanced")
        root.addWidget(tabs, 1)

        foot = QFrame()
        foot.setObjectName("dlgFoot")
        fl = QHBoxLayout(foot)
        fl.setContentsMargins(22, 12, 22, 12)
        fl.setSpacing(8)
        self.foot_info = QLabel("")
        self.foot_info.setObjectName("hint")
        fl.addWidget(self.foot_info, 1)
        defaults = QPushButton("Defaults")
        defaults.clicked.connect(self._defaults)
        cancel = QPushButton("Cancel")
        cancel.clicked.connect(self.reject)
        self.run_btn = QPushButton("  Run assembly  ")
        self.run_btn.setObjectName("primary")
        self.run_btn.setDefault(True)
        self.run_btn.clicked.connect(self._accept)
        fl.addWidget(defaults)
        fl.addWidget(cancel)
        fl.addWidget(self.run_btn)
        root.addWidget(foot)
        self._mode_changed()
        self._trim_changed()
        if files:
            self.add_paths(files)

    @staticmethod
    def _scrolled(w):
        sa = QScrollArea()
        sa.setWidgetResizable(True)
        sa.setFrameShape(QFrame.NoFrame)
        sa.setWidget(w)
        return sa

    # ---- Base tab ---------------------------------------------------------------------------------------------------
    def _base_tab(self):
        w = QWidget()
        v = QVBoxLayout(w)
        v.setContentsMargins(20, 16, 20, 12)
        v.setSpacing(12)

        row = QHBoxLayout()
        row.setSpacing(12)
        self.rb_denovo = QRadioButton()
        self.rb_ref = QRadioButton()
        self.mode_group = QButtonGroup(self)
        self.mode_group.addButton(self.rb_denovo)
        self.mode_group.addButton(self.rb_ref)
        self.rb_denovo.setChecked(True)
        row.addWidget(OptionCard(self.rb_denovo, "De novo", "One Forward + one Reverse read per sample, no reference needed."))
        row.addWidget(OptionCard(self.rb_ref, "Reference-guided", "Any number of reads per sample, mapped to a reference sequence."))
        v.addLayout(row)
        self.ref_edit = QLineEdit()
        self.ref_edit.setPlaceholderText("Input a file with the reference sequence(s) (FASTA)")
        rw, rb = path_row(self.ref_edit)
        rb.clicked.connect(self._pick_ref)
        self.ref_row = QWidget()
        rl = QFormLayout(self.ref_row)
        rl.setContentsMargins(0, 0, 0, 0)
        rl.addRow("Reference", rw)
        self.ref_acc = QLineEdit()
        self.ref_acc.setPlaceholderText("or an NCBI accession, for example NC_007605.1 (needs the internet)")
        rl.addRow("", self.ref_acc)
        v.addWidget(self.ref_row)
        self.rb_denovo.toggled.connect(self._mode_changed)

        g = DropFrame()
        g.dropped.connect(lambda ps: self.add_paths(expand_paths(ps, self.sub_chk.isChecked())))
        gl = QVBoxLayout(g)
        gl.setContentsMargins(14, 12, 14, 10)
        gl.setSpacing(8)
        bar = QHBoxLayout()
        bar.setSpacing(8)
        tt = QLabel("Input files  (.ab1)")
        tt.setObjectName("cardTitle")
        bar.addWidget(tt)
        bar.addStretch(1)
        for text, fn in [("Add files…", self._add_files), ("Add folder…", self._add_folder), ("Remove", self._remove), ("Clear", self._remove_all)]:
            b = QPushButton(text)
            b.clicked.connect(fn)
            bar.addWidget(b)
        gl.addLayout(bar)

        top = QHBoxLayout()
        top.setSpacing(10)
        top.addWidget(QLabel("Sample names"))
        self.naming = QComboBox()
        for k, t in NAMING:
            self.naming.addItem(t, k)
        self.naming.currentIndexChanged.connect(self._naming_changed)
        top.addWidget(self.naming, 1)
        self.custom_name = QLineEdit()
        self.custom_name.setPlaceholderText("Sample name")
        self.custom_name.setVisible(False)
        self.custom_name.editingFinished.connect(self._rescan)
        top.addWidget(self.custom_name)
        gl.addLayout(top)

        self.table = QTableWidget(0, 4)
        self.table.setHorizontalHeaderLabels(["File", "Sample", "Direction", "Folder"])
        hh = self.table.horizontalHeader()
        hh.setSectionResizeMode(0, QHeaderView.Stretch)
        hh.setSectionResizeMode(1, QHeaderView.ResizeToContents)
        hh.setSectionResizeMode(2, QHeaderView.ResizeToContents)
        hh.setSectionResizeMode(3, QHeaderView.Stretch)
        self.table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.table.setAlternatingRowColors(True)
        self.table.setShowGrid(False)
        self.table.verticalHeader().setVisible(False)
        self.table.setEditTriggers(QAbstractItemView.DoubleClicked)
        self.table.itemChanged.connect(self._item_changed)
        self.table.setMinimumHeight(170)
        gl.addWidget(self.table, 1)

        foot = QHBoxLayout()
        self.sub_chk = QCheckBox("Include subfolders when adding a folder")
        self.sub_chk.setChecked(True)
        foot.addWidget(self.sub_chk)
        foot.addStretch(1)
        self.info = QLabel("Drop .ab1 files or folders here, or use Add files / Add folder.")
        self.info.setObjectName("hint")
        foot.addWidget(self.info)
        gl.addLayout(foot)
        v.addWidget(g, 1)

        return w

    # ---- Advanced tab: every numeric setting is a number box -------------------------------------------
    def _advanced_tab(self):
        w = QWidget()
        grid = QGridLayout(w)
        grid.setContentsMargins(20, 16, 20, 12)
        grid.setHorizontalSpacing(14)
        grid.setVerticalSpacing(8)

        g, f = _form_group("Trimming")
        self.trim_method = QComboBox()
        self.trim_method.setMinimumWidth(250)
        self.trim_method.addItem("M2: sliding window", "M2")
        self.trim_method.addItem("M1: modified Mott", "M1")
        self.trim_method.currentIndexChanged.connect(self._trim_changed)
        self.q_trim = spin(0, 60, 20, "  Phred")
        self.win = spin(1, 200, 10, "  bases")
        self.mott = dspin(0.00001, 0.5, 0.0001, dec=5, step=0.0001)
        f.addRow("Trimming method", self.trim_method)
        f.addRow("Window mean Phred ≥  (M2)", self.q_trim)
        f.addRow("Window size  (M2)", self.win)
        f.addRow("Error cutoff  (M1)", self.mott)
        grid.addWidget(g, 0, 0)

        g, f = _form_group("Base quality")
        self.q_good = spin(0, 60, 30, "  Phred")
        self.q_poor = spin(0, 60, 20, "  Phred")
        f.addRow("Good base: Phred ≥", self.q_good)
        f.addRow("Poor base: Phred <", self.q_poor)
        grid.addWidget(g, 0, 1)

        g, f = _form_group("Forward / Reverse overlap")
        self.min_overlap = spin(1, 2000, 30, "  bp")
        self.min_ident = spin(10, 100, 75, " %")
        f.addRow("Minimum overlap", self.min_overlap)
        f.addRow("Minimum identity in the overlap", self.min_ident)
        grid.addWidget(g, 1, 0)

        g, f = _form_group("Double peaks and conflicts")
        self.ratio = dspin(0.01, 1.0, 0.33, dec=2, step=0.01)
        self.dq = spin(0, 60, 15, "  Phred")
        self.type4 = QComboBox()
        self.type4.setMinimumWidth(220)
        self.type4.addItem("N  (unresolved)", "N")
        self.type4.addItem("IUPAC ambiguity code", "iupac")
        self.type4.addItem("Follow the best-supported base", "F")
        f.addRow("Secondary / primary peak ratio ≥", self.ratio)
        f.addRow("Phred difference that favours one base", self.dq)
        f.addRow("Conflict that cannot be decided", self.type4)
        grid.addWidget(g, 1, 1)

        g, f = _form_group("Reference-guided assembly")
        self.map_mode = QComboBox()
        self.map_mode.setMinimumWidth(250)
        self.map_mode.addItem("Whole read (poor parts stay)", "whole")
        self.map_mode.addItem("Only the part that matches (local)", "local")
        self.map_min_id = spin(10, 100, 55, " %")
        self.map_min_len = spin(10, 500, 50, "  bp", 10)
        self.full_ref = QCheckBox("Show the whole reference region (positions without a read become N)")
        self.full_ref.setChecked(True)
        self.keep_ins = QCheckBox("Keep insertions and read ends beyond the reference as extra columns\n"
                                  "(the reference numbering does not change; I decide what to delete)")
        self.keep_ins.setChecked(True)
        f.addRow("How the reads are placed", self.map_mode)
        f.addRow("Minimum identity of a read", self.map_min_id)
        f.addRow("Minimum aligned bases of a read", self.map_min_len)
        f.addRow(self.full_ref)
        f.addRow(self.keep_ins)
        hint = QLabel("Reads below these limits are skipped and listed in the notes. Lower the identity (30-50 %) for messy or "
                      "divergent reads: they are then placed even where they match poorly.")
        hint.setObjectName("hint")
        hint.setWordWrap(True)
        f.addRow(hint)
        self.ref_group = g
        grid.addWidget(g, 2, 0)

        g, f = _form_group("Automatic deletion")
        self.auto_del = QCheckBox("Delete the bases below a Phred threshold while assembling")
        self.auto_q = spin(0, 60, 15, "  Phred")
        f.addRow(self.auto_del)
        f.addRow("Delete bases with Phred <", self.auto_q)
        self.auto_del.toggled.connect(self.auto_q.setEnabled)
        self.auto_q.setEnabled(False)
        grid.addWidget(g, 2, 1, Qt.AlignTop)
        grid.setRowStretch(3, 1)
        grid.setColumnStretch(0, 1)
        grid.setColumnStretch(1, 1)
        return w

    def _defaults(self):
        self.trim_method.setCurrentIndex(0)
        self.q_trim.setValue(20)
        self.win.setValue(10)
        self.mott.setValue(0.0001)
        self.q_good.setValue(30)
        self.q_poor.setValue(20)
        self.min_overlap.setValue(30)
        self.min_ident.setValue(75)
        self.ratio.setValue(0.33)
        self.dq.setValue(15)
        self.type4.setCurrentIndex(0)
        self.map_mode.setCurrentIndex(0)
        self.map_min_id.setValue(55)
        self.map_min_len.setValue(50)
        self.full_ref.setChecked(True)
        self.keep_ins.setChecked(True)
        self.auto_del.setChecked(False)
        self.auto_q.setValue(15)

    def params(self):
        return {"trim_method": self.trim_method.currentData(), "q_trim": self.q_trim.value(), "win": self.win.value(),
                "mott_cutoff": self.mott.value(), "q_good": self.q_good.value(), "q_poor": self.q_poor.value(),
                "ratio": self.ratio.value(), "dq_pref": self.dq.value(), "min_overlap": self.min_overlap.value(),
                "min_ident": self.min_ident.value() / 100.0, "type4": self.type4.currentData(),
                "map_mode": self.map_mode.currentData(), "map_min_id": self.map_min_id.value(), "map_min_len": self.map_min_len.value(),
                "full_ref": self.full_ref.isChecked(), "keep_ins": self.keep_ins.isChecked(),
                "auto_del": self.auto_del.isChecked(), "auto_q": self.auto_q.value()}

    # ---- behaviour --------------------------------------------------------------------------------------------------
    def _mode_changed(self):
        self.ref_row.setVisible(self.rb_ref.isChecked())
        self.ref_group.setEnabled(self.rb_ref.isChecked())
        self.ref_group.setToolTip("" if self.rb_ref.isChecked() else "These settings are used in the Reference-guided mode.")

    def _trim_changed(self):
        m2 = self.trim_method.currentData() == "M2"
        self.q_trim.setEnabled(m2)
        self.win.setEnabled(m2)
        self.mott.setEnabled(not m2)

    def _naming_changed(self):
        self.custom_name.setVisible(self.naming.currentData() == "custom")
        self._rescan()

    def _pick_ref(self):
        p, _ = QFileDialog.getOpenFileName(self, "Reference sequence", self.cfg.value("assembly/ref_dir", ""),
                                           "FASTA (*.fa *.fasta *.fas *.fna *.txt);;All files (*)")
        if p:
            self.ref_edit.setText(p)
            self.cfg.setValue("assembly/ref_dir", str(Path(p).parent))

    def _add_files(self):
        ps, _ = QFileDialog.getOpenFileNames(self, "Add .ab1 files", self.cfg.value("assembly/in_dir", ""), "AB1 files (*.ab1);;All files (*)")
        if ps:
            self.cfg.setValue("assembly/in_dir", str(Path(ps[0]).parent))
            self.add_paths(ps)

    def _add_folder(self):
        d = QFileDialog.getExistingDirectory(self, "Folder with .ab1 files", self.cfg.value("assembly/in_dir", ""))
        if not d:
            return
        self.cfg.setValue("assembly/in_dir", d)
        ps = expand_paths([d], self.sub_chk.isChecked())
        if not ps:
            QMessageBox.information(self, "Add folder", "No .ab1 files were found in that folder.")
            return
        self.add_paths(ps)

    def add_paths(self, paths):
        have = {r["path"] for r in self.rows}
        for p in paths:
            if p not in have:
                self.rows.append({"name": Path(p).name, "path": p, "sample": Path(p).stem, "direction": None})
                have.add(p)
        self._rescan()

    def _remove(self):
        sel = sorted({i.row() for i in self.table.selectedItems()}, reverse=True)
        for r in sel:
            del self.rows[r]
        self._fill()

    def _remove_all(self):
        self.rows = []
        self._fill()

    def _rescan(self):
        if not self.rows or not self.bridge.is_ready:
            self._fill()
            return
        self.bridge.call("scan", {"paths": [r["path"] for r in self.rows], "naming": self.naming.currentData(), "custom": self.custom_name.text()},
                         on_done=lambda ok, res: self._scanned(ok, res))

    def _scanned(self, ok, res):
        if ok:
            self.rows = [{"name": f["name"], "path": f["path"], "sample": f["sample"], "direction": f.get("direction")} for f in res["files"]]
        self._fill()

    def _dir_item(self, d):
        key = str(d or "").strip().upper()[:1]
        text, color = DIR_TEXT.get(key, ("?", theme.WARN))
        it = QTableWidgetItem(text)
        it.setForeground(QColor(color))
        it.setFlags(it.flags() & ~Qt.ItemIsEditable)
        return it

    def _fill(self):
        self.table.blockSignals(True)
        self.table.setRowCount(len(self.rows))
        for i, r in enumerate(self.rows):
            for j, val in enumerate([r["name"], r["sample"], None, str(Path(r["path"]).parent)]):
                if j == 2:
                    self.table.setItem(i, j, self._dir_item(r.get("direction")))
                    continue
                it = QTableWidgetItem(val)
                if j != 1:
                    it.setFlags(it.flags() & ~Qt.ItemIsEditable)
                self.table.setItem(i, j, it)
        self.table.blockSignals(False)
        self._update_info()

    def _update_info(self):
        ns = len({r["sample"] for r in self.rows})
        if self.rows:
            self.info.setText(f"{len(self.rows)} file(s) in {ns} sample(s) — double-click a Sample cell to rename it.")
            self.foot_info.setText("Files with the same sample name are assembled together.")
        else:
            self.info.setText("Drop .ab1 files or folders here, or use Add files / Add folder.")
            self.foot_info.setText("")

    def _item_changed(self, it):
        if it.column() == 1 and it.row() < len(self.rows):
            self.rows[it.row()]["sample"] = it.text().strip() or self.rows[it.row()]["sample"]
            self.rows[it.row()]["manual"] = True
            self._update_info()

    def _accept(self):
        if not self.rows:
            QMessageBox.warning(self, "Sanger Assembly", "Add at least one .ab1 file.")
            return
        if self.rb_ref.isChecked() and not Path(self.ref_edit.text().strip()).is_file() and not self.ref_acc.text().strip():
            QMessageBox.warning(self, "Sanger Assembly", "Choose the reference FASTA file, or type an NCBI accession.")
            return
        self.result_request = {
            "mode": "reference" if self.rb_ref.isChecked() else "denovo",
            "reference": self.ref_edit.text().strip() if self.rb_ref.isChecked() and Path(self.ref_edit.text().strip()).is_file() else "",
            "ref_acc": self.ref_acc.text().strip() if self.rb_ref.isChecked() else "",
            "files": [{"name": r["name"], "path": r["path"], "sample": r["sample"]} for r in self.rows],
            "params": self.params(), "out_dir": "", "save": []}      # saving is the last step: File > Save assembly results
        self.accept()
