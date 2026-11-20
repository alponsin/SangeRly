"""Dialogs for the tools that used to be stacked in one web page.

MsaDialog   Tools > Multiple sequence alignment > Align with ...   (input and output, method and its details, advanced options)
TrimDialog  Tools > Multiple sequence alignment > Trim the alignment ...
TreeDialog  Tools > Phylogeny > Build a tree ...                    (tree building method, tabs: model / bootstrapping / display)

A dialog only collects the choices. The main window sends them to the Shiny tool (sgQt.msaRun, sgQt.msaTrimRun, sgQt.apply), which
does the work with the same R code as before; the result then shows in a document tab.
"""
from pathlib import Path

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (QAbstractItemView, QCheckBox, QComboBox, QDialog, QDialogButtonBox, QDoubleSpinBox, QFileDialog, QFormLayout, QGroupBox,
                               QHBoxLayout, QLabel, QLineEdit, QListWidget, QMessageBox, QPlainTextEdit, QPushButton, QRadioButton,
                               QSpinBox, QStackedWidget, QTabWidget, QVBoxLayout, QWidget)

from . import theme

FASTA_FILTER = "FASTA files (*.fa *.fasta *.fas *.fna *.ffn *.txt);;All files (*)"
FASTA_EXT = (".fa", ".fasta", ".fas", ".fna", ".ffn", ".txt")

MSA_METHODS = [("Align with DECIPHER (AlignSeqs)", "decipher"), ("Align with Clustal Omega", "ClustalOmega"),
               ("Align with ClustalW", "ClustalW"), ("Align with MUSCLE", "Muscle"),
               ("Align to a reference (like ViralMSA)", "refguided"), ("I already have an alignment (trim only)", "prealigned")]
MSA_INFO = {
    "decipher": "Progressive alignment of DNA with the R package DECIPHER. A good default for barcodes and Sanger consensus sequences.",
    "ClustalOmega": "Progressive alignment that stays fast with many sequences.",
    "ClustalW": "The classic progressive method.",
    "Muscle": "Iterative refinement: accurate for small and medium sets of sequences.",
    "refguided": "Every sequence is aligned to one reference (like ViralMSA) and the numbering follows the reference. "
                 "minimap2 is used when it is installed, otherwise a built-in aligner. Sequences that cannot be mapped are left out and listed.",
    "prealigned": "The sequences already are an alignment (all the same length). Nothing is aligned: the sequences go straight to trimming.",
}
TREE_METHODS = [("Neighbor-Joining (distance matrix)", "nj"), ("Maximum Likelihood (IQ-TREE)", "ml"),
                ("Best substitution model only (ModelFinder)", "mf")]


def clean(text):
    return text.replace("…", "").replace("...", "").strip()


def spin(lo, hi, val, step=1, suffix=""):
    s = QSpinBox()
    s.setRange(lo, hi)
    s.setValue(val)
    s.setSingleStep(step)
    if suffix:
        s.setSuffix(suffix)
    s.setAlignment(Qt.AlignRight)
    s.setMinimumWidth(110)
    return s


def dspin(lo, hi, val, step=0.05, dec=2):
    s = QDoubleSpinBox()
    s.setRange(lo, hi)
    s.setDecimals(dec)
    s.setSingleStep(step)
    s.setValue(val)
    s.setAlignment(Qt.AlignRight)
    s.setMinimumWidth(110)
    return s


def combo(items, current=None):
    """items: [(text, data)]"""
    c = QComboBox()
    for t, d in items:
        c.addItem(t, d)
    if current is not None:
        i = c.findData(current)
        if i >= 0:
            c.setCurrentIndex(i)
    return c


class _Dialog(QDialog):
    def __init__(self, parent, title, ok_text, help_text=""):
        super().__init__(parent)
        self.setWindowTitle(title)
        self.setModal(True)
        self.setMinimumWidth(580)
        self.help_text = help_text
        self.root = QVBoxLayout(self)
        self.root.setContentsMargins(14, 10, 14, 12)
        self.root.setSpacing(8)
        self.buttons = QDialogButtonBox()
        self.ok = self.buttons.addButton(ok_text, QDialogButtonBox.AcceptRole)
        self.ok.setObjectName("primary")
        self.ok.setDefault(True)
        self.buttons.addButton(QDialogButtonBox.Cancel)
        self.buttons.addButton(QDialogButtonBox.Help)
        self.buttons.accepted.connect(self._accept)
        self.buttons.rejected.connect(self.reject)
        self.buttons.helpRequested.connect(self._help)

    def finish(self):
        """Put everything added so far into a scroll area, so a tall dialog never pushes the buttons off the screen."""
        from PySide6.QtGui import QGuiApplication
        from PySide6.QtWidgets import QFrame, QScrollArea, QVBoxLayout, QWidget
        inner = QWidget()
        il = QVBoxLayout(inner)
        il.setContentsMargins(0, 0, 0, 0)
        il.setSpacing(self.root.spacing())
        while self.root.count():
            it = self.root.takeAt(0)
            if it.widget() is not None:
                il.addWidget(it.widget())
            elif it.layout() is not None:
                il.addLayout(it.layout())
            else:
                il.addItem(it)
        sa = QScrollArea()
        sa.setWidgetResizable(True)
        sa.setFrameShape(QFrame.NoFrame)
        sa.setWidget(inner)
        self.root.addWidget(sa, 1)
        self.root.addWidget(self.buttons)
        scr = self.screen() or QGuiApplication.primaryScreen()
        if scr:
            self.setMaximumHeight(scr.availableGeometry().height() - 40)

    def _help(self):
        QMessageBox.information(self, self.windowTitle(), self.help_text or "No help for this window yet.")

    def _accept(self):
        if self.check():
            self.accept()

    def check(self):
        return True


# ======================================================================================================================================
class MsaDialog(_Dialog):
    """Input (files, sequences already in the MSA list, BLAST hits, Assembly consensus, pasted FASTA) + output + method + options."""

    def __init__(self, parent, code, state, cfg):
        super().__init__(parent, "Align", "Align",
                         "Input: the sequences already in the MSA list (for example sent from BLAST), FASTA files, the consensus "
                         "sequences of the Assembly and / or pasted FASTA are put together and aligned.\n\n"
                         "Output: optional. The alignment is also saved as a FASTA file there.\n\n"
                         "The alignment opens in the MSA Editor tab; Tools > Multiple sequence alignment > Trim the alignment trims it.")
        self.cfg = cfg
        st = state or {}
        self.known = state is not None
        self.pool = [str(x) for x in (st.get("pool") or [])]
        self.n_asm = int(st.get("asm") or 0)

        # ---- input and output -------------------------------------------------------------------------------------------------
        g = QGroupBox("Input and output")
        gl = QVBoxLayout(g)
        self.chk_keep = QCheckBox(f"Use the {len(self.pool)} sequence(s) already in the MSA list (BLAST hits, earlier work)" if self.known
                                  else "Use the sequences already in the MSA list (BLAST hits, earlier work)")
        self.chk_keep.setChecked(bool(self.pool) or not self.known)
        self.chk_keep.setEnabled(bool(self.pool) or not self.known)
        gl.addWidget(self.chk_keep)
        self.chk_asm = QCheckBox(f"Add the {self.n_asm} consensus sequence(s) of the Assembly" if self.known else "Add the consensus sequences of the Assembly")
        self.chk_asm.setEnabled(self.n_asm > 0 or not self.known)
        gl.addWidget(self.chk_asm)
        gl.addWidget(QLabel("FASTA files (one or many sequences per file)"))
        row = QHBoxLayout()
        self.files = QListWidget()
        self.files.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.files.setMinimumHeight(84)
        row.addWidget(self.files, 1)
        col = QVBoxLayout()
        for text, slot in (("Add files…", self._add_files), ("Add folder…", self._add_folder), ("Remove", self._remove), ("Clear", self.files.clear)):
            b = QPushButton(text)
            b.clicked.connect(slot)
            col.addWidget(b)
        col.addStretch(1)
        row.addLayout(col)
        gl.addLayout(row)
        self.chk_sub = QCheckBox("with subfolders (Add folder)")
        self.chk_sub.setChecked(True)
        gl.addWidget(self.chk_sub)
        gl.addWidget(QLabel("Or paste FASTA"))
        self.paste = QPlainTextEdit()
        self.paste.setPlaceholderText(">sample1\nACGT…\n>sample2\nACGT…")
        self.paste.setFixedHeight(66)
        gl.addWidget(self.paste)
        orow = QHBoxLayout()
        orow.addWidget(QLabel("Output file"))
        self.out = QLineEdit()
        self.out.setPlaceholderText("(optional) save the alignment as FASTA")
        orow.addWidget(self.out, 1)
        ob = QPushButton("Browse…")
        ob.clicked.connect(self._pick_out)
        orow.addWidget(ob)
        gl.addLayout(orow)
        self.root.addWidget(g)

        # ---- method ---------------------------------------------------------------------------------------------------------------
        mrow = QHBoxLayout()
        mrow.addWidget(QLabel("Method"))
        self.method = combo([(clean(t), c) for t, c in MSA_METHODS], code)
        mrow.addWidget(self.method, 1)
        self.root.addLayout(mrow)
        self.details = QLabel()
        self.details.setWordWrap(True)
        self.details.setObjectName("hint")
        self.details.setMinimumHeight(44)
        self.details.setStyleSheet(f"background:#fff;border:1px solid {theme.BORDER2};border-radius:6px;padding:6px 8px;")
        self.root.addWidget(self.details)

        # ---- advanced options -----------------------------------------------------------------------------------------------------
        self.adv = QGroupBox("Advanced options")
        al = QVBoxLayout(self.adv)
        self.chk_orient = QCheckBox("Reverse-complement sequences that fit better in the opposite orientation")
        self.chk_orient.setChecked(True)
        al.addWidget(self.chk_orient)
        self.ref_box = QWidget()
        rl = QFormLayout(self.ref_box)
        rl.setContentsMargins(0, 0, 0, 0)
        self.ref = QComboBox()
        self.ref.addItem("(the first sequence of the list)", "")
        for n in self.pool:
            self.ref.addItem(n, n)
        rl.addRow("Reference sequence", self.ref)
        frow = QHBoxLayout()
        self.ref_file = QLineEdit()
        self.ref_file.setPlaceholderText("(optional) a FASTA file with the reference")
        frow.addWidget(self.ref_file, 1)
        fb = QPushButton("Browse…")
        fb.clicked.connect(self._pick_ref)
        frow.addWidget(fb)
        rl.addRow("Reference file", frow)
        self.chk_ins = QCheckBox("Keep the insertions as extra columns (otherwise they are left out, like ViralMSA)")
        rl.addRow("", self.chk_ins)
        al.addWidget(self.ref_box)
        self.root.addWidget(self.adv)

        self.method.currentIndexChanged.connect(self._method_changed)
        self._method_changed()
        self.finish()
        self.resize(620, 640)

    # ---- helpers ------------------------------------------------------------------------------------------------------------------
    def code(self):
        return self.method.currentData()

    def _method_changed(self):
        c = self.code()
        self.setWindowTitle(clean(self.method.currentText()))
        self.details.setText("<b>Method details:</b> " + MSA_INFO.get(c, ""))
        self.chk_orient.setEnabled(c not in ("prealigned", "refguided"))
        self.ref_box.setVisible(c == "refguided")
        self.adv.setVisible(c != "prealigned")
        self.ok.setText("Load alignment" if c == "prealigned" else "Align")

    def _dir(self):
        return self.cfg.value("msa/in_dir", "")

    def _add_files(self):
        ps, _ = QFileDialog.getOpenFileNames(self, "Add FASTA files", self._dir(), FASTA_FILTER)
        for p in ps:
            self._add(p)
        if ps:
            self.cfg.setValue("msa/in_dir", str(Path(ps[0]).parent))

    def _add_folder(self):
        d = QFileDialog.getExistingDirectory(self, "Folder with FASTA files", self._dir())
        if not d:
            return
        self.cfg.setValue("msa/in_dir", d)
        it = Path(d).rglob("*") if self.chk_sub.isChecked() else Path(d).glob("*")
        found = sorted(str(x) for x in it if x.is_file() and x.suffix.lower() in FASTA_EXT)
        if not found:
            QMessageBox.information(self, "Add folder", "No FASTA files (.fa .fasta .fas .fna .ffn .txt) in that folder.")
        for p in found:
            self._add(p)

    def _add(self, p):
        if p not in [self.files.item(i).text() for i in range(self.files.count())]:
            self.files.addItem(p)

    def _remove(self):
        for it in self.files.selectedItems():
            self.files.takeItem(self.files.row(it))

    def _pick_out(self):
        p, _ = QFileDialog.getSaveFileName(self, "Save the alignment as", self._dir(), "FASTA (*.fasta *.fa *.fas);;All files (*)")
        if p:
            self.out.setText(p)

    def _pick_ref(self):
        p, _ = QFileDialog.getOpenFileName(self, "Reference sequence", self._dir(), FASTA_FILTER)
        if p:
            self.ref_file.setText(p)

    def file_list(self):
        return [self.files.item(i).text() for i in range(self.files.count())]

    def check(self):
        has = (self.chk_keep.isChecked() and self.chk_keep.isEnabled()) or (self.chk_asm.isChecked() and self.chk_asm.isEnabled()) \
            or self.files.count() > 0 or bool(self.paste.toPlainText().strip())
        if not has:
            QMessageBox.information(self, self.windowTitle(), "Give the sequences first: add FASTA files, paste FASTA, add the Assembly "
                                    "consensus or keep the sequences that are already in the MSA list.")
            return False
        return True

    def request(self):
        c = self.code()
        return {"method": c, "keep": bool(self.chk_keep.isChecked() and self.chk_keep.isEnabled()),
                "files": self.file_list(), "asm": bool(self.chk_asm.isChecked() and self.chk_asm.isEnabled()),
                "paste": self.paste.toPlainText().strip(), "out": self.out.text().strip(),
                "orient": bool(self.chk_orient.isChecked() and self.chk_orient.isEnabled()),
                "ref": self.ref.currentData() or "", "ref_file": self.ref_file.text().strip(), "ref_ins": self.chk_ins.isChecked()}


# ======================================================================================================================================
class TrimDialog(_Dialog):
    def __init__(self, parent, state, cfg):
        super().__init__(parent, "Trim the alignment", "Trim",
                         "Columns of the alignment that do not pass the rule are removed from the trimmed result. They stay in the Editor "
                         "(faded), so you can bring them back there.\n\nBuilt-in rules follow trimAl: gap threshold (-gt), conservation "
                         "threshold (-st), both, or ragged ends only. The trimAl program is used when it is installed.")
        st = state or {}
        n = int(st.get("nseq") or 0)
        lab = QLabel(f"Alignment in the MSA Editor: {n} sequences." if n else "Alignment in the MSA Editor.")
        lab.setObjectName("hint")
        self.root.addWidget(lab)
        self.engine = combo([("Built-in rules (trimAl style)", "native"), ("trimAl program (if installed)", "trimal")], "native")
        er = QHBoxLayout()
        er.addWidget(QLabel("Trimming engine"))
        er.addWidget(self.engine, 1)
        self.root.addLayout(er)
        self.stack = QStackedWidget()

        nat = QGroupBox("Built-in rules")
        f = QFormLayout(nat)
        self.mode = combo([("Gap threshold (-gt)", "gt"), ("Conservation threshold (-st)", "st"), ("Gap AND conservation", "both"),
                           ("Ragged ends only", "ends"), ("No trimming", "none")], "gt")
        f.addRow("Rule", self.mode)
        self.gt = dspin(0, 1, 0.5)
        self.st = dspin(0, 1, 0.5)
        self.ends = dspin(0, 1, 0.5)
        self.cons = spin(0, 100, 0, 5, " %")
        f.addRow("Minimum fraction of sequences with a base in the column", self.gt)
        f.addRow("Minimum conservation of the column (most common base)", self.st)
        f.addRow("Ends: first / last column where this fraction has a base", self.ends)
        f.addRow("Always keep at least this % of the columns (-cons)", self.cons)
        self.stack.addWidget(nat)

        tm = QGroupBox("trimAl")
        g = QFormLayout(tm)
        self.tflag = combo([("automated1 (recommended for ML trees)", "-automated1"), ("gappyout", "-gappyout"), ("strict", "-strict"),
                            ("strictplus", "-strictplus"), ("nogaps (remove all columns with gaps)", "-nogaps"),
                            ("noallgaps (remove only all-gap columns)", "-noallgaps"), ("custom (-gt / -st / -cons)", "custom")], "-automated1")
        g.addRow("trimAl mode", self.tflag)
        self.t_gt = dspin(0, 1, 0.5)
        self.t_st = dspin(0, 1, 0.0)
        self.t_cons = spin(0, 100, 0, 5, " %")
        g.addRow("-gt", self.t_gt)
        g.addRow("-st (0 = not used)", self.t_st)
        g.addRow("-cons", self.t_cons)
        self.stack.addWidget(tm)
        self.root.addWidget(self.stack)

        self.engine.currentIndexChanged.connect(lambda: self.stack.setCurrentIndex(0 if self.engine.currentData() == "native" else 1))
        self.mode.currentIndexChanged.connect(self._enable)
        self.tflag.currentIndexChanged.connect(self._enable)
        self._enable()
        self.finish()

    def _enable(self):
        m = self.mode.currentData()
        self.gt.setEnabled(m in ("gt", "both"))
        self.st.setEnabled(m in ("st", "both"))
        self.ends.setEnabled(m == "ends")
        self.cons.setEnabled(m != "none")
        c = self.tflag.currentData() == "custom"
        for w in (self.t_gt, self.t_st, self.t_cons):
            w.setEnabled(c)

    def request(self):
        return {"engine": self.engine.currentData(), "mode": self.mode.currentData(), "gt": self.gt.value(), "st": self.st.value(),
                "ends_cov": self.ends.value(), "cons": self.cons.value(), "tflag": self.tflag.currentData(),
                "t_gt": self.t_gt.value(), "t_st": self.t_st.value(), "t_cons": self.t_cons.value()}


# ======================================================================================================================================
class TreeDialog(_Dialog):
    """Build a tree: the tree building method on top, then tabs for the model, bootstrapping and display."""

    def __init__(self, parent, code, state, cfg):
        super().__init__(parent, "Build Phylogenetic Tree", "Build",
                         "Neighbor-Joining uses a distance matrix (R package ape). Maximum Likelihood runs IQ-TREE (ModelFinder picks the "
                         "substitution model unless you fix one). ModelFinder only reports the best model without building the final tree.\n\n"
                         "The alignment is the one of the MSA tab. The tree opens in its own tab; the display options can be changed there too.")
        self.cfg = cfg
        st = state or {}
        n = int(st.get("nseq") or 0)
        top = QFormLayout()
        self.method = combo(TREE_METHODS, code)
        top.addRow("Tree building method", self.method)
        self.src = combo([("MSA tab: kept columns (after trimming / editing)", "trim"), ("MSA tab: whole alignment", "all"),
                           ("An aligned FASTA file", "file")], "trim")
        top.addRow("Alignment", self.src)
        self.fa = QLineEdit()
        self.fa.setPlaceholderText("aligned FASTA file (at least 3 sequences, same length)")
        self.fa_btn = QPushButton("Browse…")
        self.fa_btn.clicked.connect(self._pick_fa)
        fr = QHBoxLayout()
        fr.addWidget(self.fa, 1)
        fr.addWidget(self.fa_btn)
        self.fa_row = QWidget()
        self.fa_row.setLayout(fr)
        fr.setContentsMargins(0, 0, 0, 0)
        top.addRow("Aligned FASTA", self.fa_row)
        self.src.currentIndexChanged.connect(self._src_changed)
        self.root.addLayout(top)
        if st and not n:  # no alignment in the MSA tab (an aligned FASTA file still works)
            warn = QLabel("There is no alignment yet. Align the sequences first (Tools > MSA and trimming > Multiple sequence alignment).")
            warn.setWordWrap(True)
            warn.setStyleSheet(f"color:{theme.WARN};")
            self.root.addWidget(warn)

        self.tabs = QTabWidget()
        self.tabs.setObjectName("inner")

        # tab 1: model ---------------------------------------------------------------------------------------------------------------
        self.s_model = QStackedWidget()
        p = QWidget()
        f = QFormLayout(p)
        self.nj_model = combo([("Auto: best-fit model from IQ-TREE ModelFinder", "auto"), ("JC69", "JC69"), ("K80 (Kimura 2-parameter)", "K80"),
                               ("F81", "F81"), ("K81", "K81"), ("F84", "F84"), ("TN93", "TN93"), ("p-distance (raw)", "raw"), ("LogDet", "logdet")], "K80")
        f.addRow("Distance matrix model", self.nj_model)
        self.nj_alpha = dspin(0, 100, 0.0, 0.1, 2)
        f.addRow("Gamma shape alpha (0 = equal rates)", self.nj_alpha)
        self.s_model.addWidget(p)
        p = QWidget()
        f = QFormLayout(p)
        self.iq_model = combo([("ModelFinder + tree (-m MFP)", "MFP"), ("Fixed model", "fixed")], "MFP")
        f.addRow("Substitution model", self.iq_model)
        self.iq_fixed = QLineEdit("GTR+G")
        f.addRow("Model", self.iq_fixed)
        self.iq_threads = QLineEdit("AUTO")
        f.addRow("Threads (AUTO or a number)", self.iq_threads)
        self.iq_seed = QLineEdit()
        self.iq_seed.setPlaceholderText("empty = random")
        f.addRow("Random seed", self.iq_seed)
        self.iq_extra = QLineEdit()
        self.iq_extra.setPlaceholderText("-mset JC,K80,HKY,GTR")
        f.addRow("Extra IQ-TREE arguments", self.iq_extra)
        self.s_model.addWidget(p)
        p = QWidget()
        f = QFormLayout(p)
        self.mf_threads = QLineEdit("AUTO")
        f.addRow("Threads (AUTO or a number)", self.mf_threads)
        self.s_model.addWidget(p)
        self.tabs.addTab(self.s_model, "Model")

        # tab 2: bootstrapping -------------------------------------------------------------------------------------------------------
        self.s_boot = QStackedWidget()
        p = QWidget()
        f = QFormLayout(p)
        self.nj_boot = spin(0, 10000, 100, 100)
        f.addRow("Bootstrap replicates (0 = none)", self.nj_boot)
        self.s_boot.addWidget(p)
        p = QWidget()
        f = QFormLayout(p)
        self.iq_boot = combo([("Ultrafast bootstrap (-B)", "ufboot"), ("Standard bootstrap (-b, slow)", "std"), ("None", "none")], "ufboot")
        f.addRow("Branch support", self.iq_boot)
        self.iq_nboot = spin(100, 100000, 1000, 100)
        f.addRow("Replicates (ultrafast: at least 1000)", self.iq_nboot)
        self.iq_alrt = QCheckBox("Also SH-aLRT (-alrt 1000)")
        f.addRow("", self.iq_alrt)
        self.s_boot.addWidget(p)
        p = QWidget()
        QVBoxLayout(p).addWidget(QLabel("Not used: ModelFinder only reports the best model."))
        self.s_boot.addWidget(p)
        self.tabs.addTab(self.s_boot, "Bootstrapping and Consensus tree")

        # tab 3: display -------------------------------------------------------------------------------------------------------------
        p = QWidget()
        f = QFormLayout(p)
        self.layout_ = combo([("Phylogram", "phylogram"), ("Cladogram", "cladogram"), ("Unrooted", "unrooted"), ("Fan", "fan"), ("Radial", "radial")], "phylogram")
        f.addRow("Layout", self.layout_)
        self.root_ = combo([("As computed", "none"), ("On the outgroup", "outgroup"), ("Midpoint (needs phangorn)", "midpoint")], "none")
        f.addRow("Root the drawn tree", self.root_)
        self.show_sup = QCheckBox("Show branch support")
        self.show_sup.setChecked(True)
        f.addRow("", self.show_sup)
        self.sup_dots = QCheckBox("Support as coloured dots (green ≥ 95, orange 70–94, red < 70)")
        f.addRow("", self.sup_dots)
        self.min_sup = spin(0, 100, 0, 5)
        f.addRow("Hide support below", self.min_sup)
        self.cex = dspin(0.4, 2.0, 0.9, 0.1, 1)
        f.addRow("Label size", self.cex)
        self.scale = QCheckBox("Scale bar (phylogram)")
        self.scale.setChecked(True)
        f.addRow("", self.scale)
        self.tabs.addTab(p, "Display Options")

        # tab 4: where the files are saved (IQ-TREE writes many files: log, report, treefile, ...) ---------------------------------
        p = QWidget()
        v = QVBoxLayout(p)
        f = QFormLayout()
        default_dir = self.cfg.value("tree/out_dir", "") or str(Path.home() / "Documents" / "SangeRly" / "trees")
        self.out_dir = QLineEdit(default_dir)
        self.out_dir.setPlaceholderText("empty = a temporary folder (deleted when the program closes)")
        b = QPushButton("Browse…")
        b.clicked.connect(self._pick_out)
        r = QHBoxLayout()
        r.setContentsMargins(0, 0, 0, 0)
        r.addWidget(self.out_dir, 1)
        r.addWidget(b)
        w = QWidget()
        w.setLayout(r)
        f.addRow("Save the files in this folder", w)
        self.out_name = QLineEdit()
        self.out_name.setPlaceholderText("empty = tree / iqtree")
        f.addRow("Name for the files", self.out_name)
        v.addLayout(f)
        note = QLabel("Every run gets its own sub-folder here (the name, the date and the time) with ALL the files of IQ-TREE: "
                      "the log, the report (.iqtree), the tree file with support values, the model table and the rest. "
                      "Afterwards, the section Tree > Save results saves figures (Newick, PNG, PDF), distance tables and heat maps.")
        note.setWordWrap(True)
        note.setStyleSheet(f"color:{theme.MUTED};" if hasattr(theme, "MUTED") else "")
        v.addWidget(note)
        v.addStretch(1)
        self.tabs.addTab(p, "Save results")
        self.root.addWidget(self.tabs, 1)

        self.iq_model.currentIndexChanged.connect(self._enable)
        self.method.currentIndexChanged.connect(self._method_changed)
        self._method_changed()
        self._src_changed()
        self.finish()
        self.resize(600, 520)

    def code(self):
        return self.method.currentData()

    def _pick_fa(self):
        p, _ = QFileDialog.getOpenFileName(self, "Aligned FASTA", "", FASTA_FILTER)
        if p:
            self.fa.setText(p)

    def _src_changed(self):
        self.fa_row.setVisible(self.src.currentData() == "file")

    def _pick_out(self):
        d = QFileDialog.getExistingDirectory(self, "Folder for the tree files", self.out_dir.text() or str(Path.home()))
        if d:
            self.out_dir.setText(d)

    def check(self):
        if self.src.currentData() == "file" and not Path(self.fa.text().strip()).is_file():
            QMessageBox.information(self, self.windowTitle(), "Choose the aligned FASTA file first.")
            return False
        d = self.out_dir.text().strip()
        if d and self.code() != "nj":
            try:
                Path(d).mkdir(parents=True, exist_ok=True)
            except OSError as e:
                QMessageBox.information(self, self.windowTitle(), f"This folder cannot be used:\n{d}\n\n{e}")
                return False
            self.cfg.setValue("tree/out_dir", d)
        return True

    def _enable(self):
        self.iq_fixed.setEnabled(self.iq_model.currentData() == "fixed")

    def _method_changed(self):
        i = {"nj": 0, "ml": 1, "mf": 2}[self.code()]
        self.s_model.setCurrentIndex(i)
        self.s_boot.setCurrentIndex(i)
        self.tabs.setTabEnabled(2, self.code() != "mf")
        self.tabs.setTabEnabled(3, self.code() != "nj")
        self.ok.setText("Find the model" if self.code() == "mf" else "Build")
        self._enable()

    def _fasta(self):
        return self.fa.text().strip() if self.src.currentData() == "file" else ""

    def request(self):
        """{set: {input id of the Shiny tree tool: value}, click: id of its button, open: tab to show}"""
        c = self.code()
        t = "tree-"
        s = {t + "src": self.src.currentData(), t + "layout": self.layout_.currentData(), t + "root": self.root_.currentData(),
             t + "show_sup": self.show_sup.isChecked(), t + "sup_dots": self.sup_dots.isChecked(), t + "min_sup": self.min_sup.value(),
             t + "cex": self.cex.value(), t + "scale": self.scale.isChecked()}
        if self.out_name.text().strip():
            s[t + "sv_name"] = self.out_name.text().strip()
        send = [] if c == "nj" else [{"name": t + "qt_outdir", "value": self.out_dir.text().strip().replace("\\", "/")}]
        if c == "nj":
            s.update({t + "nj_model": self.nj_model.currentData(), t + "nj_alpha": self.nj_alpha.value(), t + "nj_boot": self.nj_boot.value()})
            return {"set": s, "send": send, "click": t + "run_nj", "open": "NJ tree", "fasta": self._fasta()}
        if c == "ml":
            s.update({t + "iq_model": self.iq_model.currentData(), t + "iq_fixed": self.iq_fixed.text().strip() or "GTR+G",
                      t + "iq_boot": self.iq_boot.currentData(), t + "iq_nboot": self.iq_nboot.value(), t + "iq_alrt": self.iq_alrt.isChecked(),
                      t + "iq_threads": self.iq_threads.text().strip() or "AUTO", t + "iq_seed": self.iq_seed.text().strip(),
                      t + "iq_extra": self.iq_extra.text().strip()})
            return {"set": s, "send": send, "click": t + "run_iq", "open": "IQ-TREE tree", "fasta": self._fasta()}
        s.update({t + "iq_model": "MF", t + "iq_threads": self.mf_threads.text().strip() or "AUTO"})
        return {"set": s, "send": send, "click": t + "run_iq", "open": "ModelFinder", "fasta": self._fasta()}


# ======================================================================================================================================
class StructureFields:
    """Settings of Structure (K clusters): method, model, K range, runs per K, seed."""

    def __init__(self, form):
        self.method = combo([("Maximum likelihood, EM (fast)", "em"), ("Bayesian MCMC, Gibbs sampler (slower)", "mcmc")], "em")
        self.model = combo([("Admixture: every sequence can have several ancestries", "admix"),
                            ("No admixture: every sequence in one cluster", "mix")], "admix")
        self.kmin = spin(1, 50, 1)
        self.kmax = spin(2, 50, 5)
        self.rep = spin(1, 100, 5)
        self.seed = QLineEdit()
        self.seed.setPlaceholderText("empty = random")
        self.rows = [("Method", self.method), ("Model", self.model), ("K from", self.kmin), ("K to", self.kmax),
                     ("Runs per K (the best likelihood is shown)", self.rep), ("Seed", self.seed)]
        for lab, w in self.rows:
            form.addRow(lab, w)

    def show_rows(self, form, on):
        for _lab, w in self.rows:
            if hasattr(form, "setRowVisible"):
                form.setRowVisible(w, on)
            else:
                w.setVisible(on)

    def values(self, t):
        v = {t + "st_method": self.method.currentData(), t + "st_model": self.model.currentData(), t + "st_kmin": self.kmin.value(),
             t + "st_kmax": max(self.kmax.value(), self.kmin.value() + 1), t + "st_rep": self.rep.value()}
        if self.seed.text().strip():
            v[t + "st_seed"] = self.seed.text().strip()
        return v


class PopDialog(_Dialog):
    """Population structure with tabs: Alignment | Populations and groups | Run. 'By hand' is long (100 sequences), so it opens its own tab."""

    def __init__(self, parent, state, cfg):
        super().__init__(parent, "Population structure", "Run",
                         "Alignment: the one of the MSA tab, or an aligned FASTA file (at least 3 sequences of the same length).\n\n"
                         "Populations and groups: from the sequence names (a rule), from a CSV / Excel table (Sequence, Population and "
                         "optionally Group) or by hand. By hand is long, so it opens the tab 'Populations and groups' where you select "
                         "sequences and name the sets; run the analysis afterwards (Tools > Population structure > Run analysis).")
        self.cfg = cfg
        self.setMinimumWidth(620)
        self.tabs = QTabWidget()
        self.tabs.setObjectName("inner")

        # ---- tab 1: alignment ------------------------------------------------------------------------------------------------------
        p = QWidget()
        f1 = QFormLayout(p)
        self.src = combo([("MSA tab: kept columns (after trimming / editing)", "trim"), ("MSA tab: whole alignment", "all"),
                          ("An aligned FASTA file", "file")], "trim")
        f1.addRow("Alignment", self.src)
        self.fa = QLineEdit()
        self.fa.setPlaceholderText("aligned FASTA file (at least 3 sequences, same length)")
        self.fa_row = self._file_row(self.fa, FASTA_FILTER)
        f1.addRow("Aligned FASTA", self.fa_row)
        self._f1 = f1
        self.tabs.addTab(p, "Alignment")

        # ---- tab 2: populations and groups -----------------------------------------------------------------------------------------
        p = QWidget()
        f2 = QFormLayout(p)
        self.pop_src = combo([("From the sequence names (a rule)", "names"), ("From a CSV / Excel file", "file"),
                              ("By hand (select sequences, name the set)", "manual")], "names")
        f2.addRow("Where the populations come from", self.pop_src)
        self.rule = combo([("Text before a separator", "prefix"), ("First characters", "chars"), ("Regular expression (first group)", "regex"),
                           ("One population", "all")], "prefix")
        f2.addRow("Population taken from the sequence name", self.rule)
        self.sep = QLineEdit("_")
        self.nfield = spin(1, 20, 1)
        self.nchars = spin(1, 50, 3)
        self.regex = QLineEdit("^([A-Za-z]+)")
        f2.addRow("Separator", self.sep)
        f2.addRow("Number of fields", self.nfield)
        f2.addRow("Number of characters", self.nchars)
        f2.addRow("Expression", self.regex)
        self.csv = QLineEdit()
        self.csv.setPlaceholderText("columns: Sequence, Population, (Group)")
        self.csv_row = self._file_row(self.csv, "Tables (*.csv *.tsv *.txt *.xlsx *.xls);;All files (*)")
        f2.addRow("CSV / Excel file", self.csv_row)
        self.hand_note = QLabel("By hand opens the tab \u201cPopulations and groups\u201d: select sequences in the table, name the set and press "
                                "\u2018Update population\u2019. Run the analysis afterwards (Tools > Population structure > Run analysis).")
        self.hand_note.setWordWrap(True)
        self.hand_note.setObjectName("hint")
        f2.addRow(self.hand_note)
        self._f2 = f2
        self.tabs.addTab(p, "Populations and groups")

        # ---- tab 3: run -------------------------------------------------------------------------------------------------------------
        p = QWidget()
        f3 = QFormLayout(p)
        self.run = combo([("Haplotypes, diversity, differentiation", "div"), ("AMOVA and pairwise FST", "amova"),
                          ("Structure (K clusters)", "struct"), ("Open the data only (no analysis yet)", "none")], "div")
        f3.addRow("Run", self.run)
        self.nsim = spin(0, 100000, 1000, 500)
        self.nperm = spin(100, 100000, 1000, 500)
        f3.addRow("Simulated samples for Fu's Fs (0 = none)", self.nsim)
        f3.addRow("Permutations (AMOVA)", self.nperm)
        self.st = StructureFields(f3)
        self._f3 = f3
        self.tabs.addTab(p, "Run")
        self.root.addWidget(self.tabs, 1)

        for c in (self.src, self.pop_src, self.rule, self.run):
            c.currentIndexChanged.connect(self._vis)
        self._vis()
        self.finish()
        self.resize(660, 420)

    def _file_row(self, edit, flt):
        b = QPushButton("Browse\u2026")
        b.clicked.connect(lambda: self._pick(edit, flt))
        r = QHBoxLayout()
        r.setContentsMargins(0, 0, 0, 0)
        r.addWidget(edit, 1)
        r.addWidget(b)
        w = QWidget()
        w.setLayout(r)
        return w

    def _pick(self, edit, flt):
        p, _ = QFileDialog.getOpenFileName(self, "Open", self.cfg.value("pop/dir", ""), flt)
        if p:
            edit.setText(p)
            self.cfg.setValue("pop/dir", str(Path(p).parent))

    @staticmethod
    def _row(form, widget, on):
        if hasattr(form, "setRowVisible"):
            form.setRowVisible(widget, on)
        else:
            widget.setVisible(on)

    def _vis(self):
        ps = self.pop_src.currentData()
        r = self.rule.currentData()
        self._row(self._f1, self.fa_row, self.src.currentData() == "file")
        self._row(self._f2, self.rule, ps == "names")
        self._row(self._f2, self.sep, ps == "names" and r == "prefix")
        self._row(self._f2, self.nfield, ps == "names" and r == "prefix")
        self._row(self._f2, self.nchars, ps == "names" and r == "chars")
        self._row(self._f2, self.regex, ps == "names" and r == "regex")
        self._row(self._f2, self.csv_row, ps == "file")
        self._row(self._f2, self.hand_note, ps == "manual")
        self.tabs.setTabEnabled(2, ps != "manual")          # by hand: run later, after the sets are made
        k = self.run.currentData()
        self._row(self._f3, self.nsim, k == "div")
        self._row(self._f3, self.nperm, k == "amova")
        self.st.show_rows(self._f3, k == "struct")
        self.ok.setText("Open" if (ps == "manual" or k == "none") else "Run")

    def check(self):
        if self.src.currentData() == "file" and not Path(self.fa.text().strip()).is_file():
            self.tabs.setCurrentIndex(0)
            QMessageBox.information(self, self.windowTitle(), "Choose the aligned FASTA file first.")
            return False
        if self.pop_src.currentData() == "file" and not Path(self.csv.text().strip()).is_file():
            self.tabs.setCurrentIndex(1)
            QMessageBox.information(self, self.windowTitle(), "Choose the CSV / Excel file with the populations first.")
            return False
        return True

    def request(self):
        t = "popgen-"
        ps = self.pop_src.currentData()
        s = {t + "src": self.src.currentData(), t + "pop_src": ps, t + "nsim": self.nsim.value(), t + "nperm": self.nperm.value()}
        send, clicks = [], []
        if self.src.currentData() == "file":
            send.append({"name": t + "qt_fasta", "value": self.fa.text().strip(), "wait": 100})
        if ps == "names":
            s.update({t + "pop_method": self.rule.currentData(), t + "pop_sep": self.sep.text(), t + "pop_nfield": self.nfield.value(),
                      t + "pop_nchars": self.nchars.value(), t + "pop_regex": self.regex.text()})
            clicks.append({"id": t + "apply_rule", "wait": 900})
        elif ps == "file":
            send.append({"name": t + "qt_popfile", "value": self.csv.text().strip(), "wait": 1200})
        run = "none" if ps == "manual" else self.run.currentData()
        if run == "struct":
            s.update(self.st.values(t))
        for c in clicks:
            c["ready"] = True                      # wait until the sequence table of the populations card is filled
            c["wait"] = 300
        if run != "none":
            clicks.append({"id": t + {"div": "run_div", "amova": "run_amova", "struct": "run_struct"}[run], "wait": 2500, "ready": True})
        opened = {"div": "Haplotype", "amova": "AMOVA", "struct": "Structure"}.get(run, "view:pops")
        return {"set": s, "send": send, "clicks": clicks, "open": opened}


class PopRunDialog(_Dialog):
    """Population structure, last step: what to run (the populations come from the tab 'Populations and groups')."""

    def __init__(self, parent, state, cfg):
        super().__init__(parent, "Population structure: run", "Run",
                         "Set the populations and groups first (tab 'Populations and groups'). The other settings (gaps, distance method, "
                         "...) are on the left of the result tabs.")
        f = QFormLayout()
        self.run = combo([("Haplotypes, diversity, differentiation", "div"), ("AMOVA and pairwise FST", "amova"),
                          ("Structure (K clusters)", "struct")], "div")
        f.addRow("Run", self.run)
        self.nsim = spin(0, 100000, 1000, 500)
        self.nperm = spin(100, 100000, 1000, 500)
        f.addRow("Simulated samples for Fu's Fs (0 = none)", self.nsim)
        f.addRow("Permutations (AMOVA)", self.nperm)
        self.st = StructureFields(f)
        self._f = f
        self.root.addLayout(f)
        self.run.currentIndexChanged.connect(self._en)
        self._en()
        self.finish()

    def _en(self):
        k = self.run.currentData()
        for w, on in ((self.nsim, k == "div"), (self.nperm, k == "amova")):
            if hasattr(self._f, "setRowVisible"):
                self._f.setRowVisible(w, on)
            else:
                w.setVisible(on)
        self.st.show_rows(self._f, k == "struct")

    def request(self):
        t = "popgen-"
        r = self.run.currentData()
        s = {t + "nsim": self.nsim.value(), t + "nperm": self.nperm.value()}
        if r == "struct":
            s.update(self.st.values(t))
        return {"set": s, "send": [], "clicks": [{"id": t + {"div": "run_div", "amova": "run_amova", "struct": "run_struct"}[r], "wait": 500}],
                "open": {"div": "Haplotype", "amova": "AMOVA", "struct": "Structure"}[r]}
