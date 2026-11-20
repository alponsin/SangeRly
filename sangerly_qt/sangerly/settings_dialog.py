"""Settings: where Rscript and the SangeRly R code are. Saved with QSettings; the engine is restarted afterwards."""
from pathlib import Path

from PySide6.QtCore import QSettings
from PySide6.QtWidgets import (QDialog, QFileDialog, QFrame, QGridLayout, QHBoxLayout, QLabel, QLineEdit, QPushButton, QVBoxLayout)

from . import theme
from .rbridge import detect_project_root, detect_rscript, find_project_root, find_rscript, is_project_root


class SettingsDialog(QDialog):
    def __init__(self, bridge, parent=None):
        super().__init__(parent)
        self.bridge = bridge
        self.cfg = QSettings()
        self.setWindowTitle("Settings")
        self.setMinimumWidth(700)
        root = QVBoxLayout(self)
        root.setContentsMargins(0, 0, 0, 0)
        root.setSpacing(0)

        head = QFrame()
        head.setObjectName("dlgHead")
        hl = QHBoxLayout(head)
        hl.setContentsMargins(22, 16, 22, 16)
        hl.setSpacing(14)
        ic = QLabel()
        ic.setPixmap(theme.icon("gear", 44, "#475569").pixmap(44, 44))
        hl.addWidget(ic)
        tv = QVBoxLayout()
        tv.setSpacing(2)
        t = QLabel("R engine")
        t.setObjectName("h2")
        s = QLabel("SangeRly does its calculations in R. Tell it where R and the SangeRly R code are.")
        s.setObjectName("hint")
        tv.addWidget(t)
        tv.addWidget(s)
        hl.addLayout(tv, 1)
        root.addWidget(head)

        body = QGridLayout()
        body.setContentsMargins(22, 18, 22, 10)
        body.setHorizontalSpacing(10)
        body.setVerticalSpacing(6)
        self.rs_edit = QLineEdit(self.cfg.value("engine/rscript", "") or bridge.rscript or find_rscript())
        self.rs_edit.setPlaceholderText("Rscript.exe  (e.g. C:\\Program Files\\R\\R-4.4.1\\bin\\Rscript.exe)")
        self.rs_state = QLabel()
        self.root_edit = QLineEdit(self.cfg.value("engine/root", "") or bridge.root or find_project_root())
        self.root_edit.setPlaceholderText("Folder that contains app.R and R\\assembly_core.R")
        self.root_state = QLabel()
        for r, (lab, ed, st, fn) in enumerate([("Rscript", self.rs_edit, self.rs_state, self._pick_rscript),
                                                ("SangeRly folder", self.root_edit, self.root_state, self._pick_root)]):
            body.addWidget(QLabel(lab), r * 3, 0, 1, 2)
            body.addWidget(ed, r * 3 + 1, 0)
            b = QPushButton("Browse…")
            b.clicked.connect(fn)
            body.addWidget(b, r * 3 + 1, 1)
            body.addWidget(st, r * 3 + 2, 0, 1, 2)
        body.setColumnStretch(0, 1)
        root.addLayout(body)
        self.rs_edit.textChanged.connect(self._check)
        self.root_edit.textChanged.connect(self._check)

        foot = QFrame()
        foot.setObjectName("dlgFoot")
        fl = QHBoxLayout(foot)
        fl.setContentsMargins(22, 12, 22, 12)
        auto = QPushButton("Auto-detect")
        auto.clicked.connect(self._auto)
        fl.addWidget(auto)
        fl.addStretch(1)
        cancel = QPushButton("Cancel")
        cancel.clicked.connect(self.reject)
        self.save_btn = QPushButton("Save and restart engine")
        self.save_btn.setObjectName("primary")
        self.save_btn.clicked.connect(self._save)
        fl.addWidget(cancel)
        fl.addWidget(self.save_btn)
        root.addWidget(foot)
        self._check()

    def _mark(self, label, ok, good, bad):
        c = theme.OK if ok else theme.BAD
        label.setText(f"<span style='color:{c}'>{'✓ ' + good if ok else '✗ ' + bad}</span>")

    def _check(self):
        p = self.rs_edit.text().strip()
        self._mark(self.rs_state, bool(p) and Path(p).is_file(), "Rscript found", "Rscript not found")
        r = self.root_edit.text().strip()
        self._mark(self.root_state, is_project_root(r), "SangeRly R code found", "R\\assembly_core.R not found in this folder")

    def _pick_rscript(self):
        p, _ = QFileDialog.getOpenFileName(self, "Choose Rscript", self.rs_edit.text() or "", "Rscript (Rscript.exe Rscript);;All files (*)")
        if p:
            self.rs_edit.setText(p)

    def _pick_root(self):
        p = QFileDialog.getExistingDirectory(self, "Choose the SangeRly folder", self.root_edit.text() or "")
        if p:
            self.root_edit.setText(p)

    def _auto(self):
        self.rs_edit.setText(detect_rscript() or self.rs_edit.text())
        self.root_edit.setText(detect_project_root() or self.root_edit.text())

    def _save(self):
        self.cfg.setValue("engine/rscript", self.rs_edit.text().strip())
        self.cfg.setValue("engine/root", self.root_edit.text().strip())
        self.accept()
