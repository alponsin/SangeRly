"""Pages shown as tabs: the Start page (hero + workflow cards + recent work) and a page for modules that still live in classic SangeRly."""
from pathlib import Path

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (QFrame, QGridLayout, QHBoxLayout, QLabel, QListWidget, QListWidgetItem, QPushButton, QScrollArea,
                               QVBoxLayout, QWidget)

from . import APP_NAME, theme
from .modules import MODULES
from .widgets import ChromatogramBanner, ModuleCard, section_title


class StartPage(QScrollArea):
    def __init__(self, win):
        super().__init__()
        self.win = win
        self.setWidgetResizable(True)
        self.setFrameShape(QFrame.NoFrame)
        page = QWidget()
        page.setObjectName("page")
        self.setWidget(page)
        v = QVBoxLayout(page)
        v.setContentsMargins(32, 26, 32, 28)
        v.setSpacing(22)
        v.addWidget(self._hero())
        v.addWidget(section_title("Workflow", "Open any step at any time. A green dot on the toolbar marks a step that already has results."))
        v.addLayout(self._cards())
        v.addLayout(self._bottom())
        v.addStretch(1)

    def _hero(self):
        hero = ChromatogramBanner()
        h = QHBoxLayout(hero)
        h.setContentsMargins(34, 26, 34, 26)
        h.setSpacing(26)
        logo = QLabel()
        pm = theme.logo_pixmap(104)
        logo.setPixmap(pm.scaled(104, 104, Qt.KeepAspectRatio, Qt.SmoothTransformation))
        logo.setFixedSize(104, 104)
        h.addWidget(logo, 0, Qt.AlignVCenter)
        t = QVBoxLayout()
        t.setSpacing(6)
        title = QLabel(f"Welcome to {APP_NAME}")
        title.setObjectName("heroTitle")
        sub = QLabel("Sanger assembly, alignment, phylogeny and population genetics — a desktop app powered by the R engine.")
        sub.setObjectName("heroSub")
        sub.setWordWrap(True)
        sub.setMaximumWidth(520)
        t.addStretch(1)
        t.addWidget(title)
        t.addWidget(sub)
        row = QHBoxLayout()
        row.setSpacing(10)
        row.setContentsMargins(0, 10, 0, 0)
        for text, name, slot in [("Start with Assembly", "heroPrimary", lambda: self.win.open_module("assembly")),
                                 ("Open .ab1 files", "heroGhost", self.win.open_ab1),
                                 ("Example data", "heroGhost", self.win.open_examples)]:
            b = QPushButton(text)
            b.setObjectName(name)
            b.setCursor(Qt.PointingHandCursor)
            b.clicked.connect(slot)
            row.addWidget(b)
        row.addStretch(1)
        t.addLayout(row)
        t.addStretch(1)
        h.addLayout(t, 1)
        h.addStretch(1)
        return hero

    def _cards(self):
        g = QGridLayout()
        g.setSpacing(16)
        n = 0
        for m in MODULES:
            ico = theme.logo_icon(self.win.logo_path(m), 48) or theme.icon(m.glyph, 48, m.color)
            c = ModuleCard(ico, m.title, m.summary, "Open", "info")
            c.clicked.connect(lambda k=m.key: self.win.open_module(k))
            g.addWidget(c, n // 4, n % 4)
            n += 1
        for col in range(4):
            g.setColumnStretch(col, 1)
        return g

    def _bottom(self):
        row = QHBoxLayout()
        row.setSpacing(16)
        box = QFrame()
        box.setObjectName("panel")
        bl = QVBoxLayout(box)
        bl.setContentsMargins(16, 14, 16, 14)
        bl.addWidget(section_title("Recent assemblies"))
        self.recent = QListWidget()
        self.recent.setFrameShape(QFrame.NoFrame)
        self.recent.setMinimumHeight(130)
        self.recent.itemDoubleClicked.connect(self._open_recent)
        bl.addWidget(self.recent, 1)
        row.addWidget(box, 3)

        box = QFrame()
        box.setObjectName("panel")
        bl = QVBoxLayout(box)
        bl.setContentsMargins(16, 14, 16, 14)
        bl.addWidget(section_title("Get started"))
        for text, slot in [("Open the example data folder", self.win.open_examples),
                           ("Take the tutorial", lambda: self.win.open_module("tutorial")),
                           ("R engine and folders…", self.win.settings_dialog)]:
            b = QPushButton(text)
            b.setObjectName("link")
            b.setCursor(Qt.PointingHandCursor)
            b.clicked.connect(slot)
            bl.addWidget(b)
        bl.addStretch(1)
        tip = QLabel("Tip: drop .ab1 files anywhere in this window to start an assembly.")
        tip.setObjectName("hint")
        tip.setWordWrap(True)
        bl.addWidget(tip)
        row.addWidget(box, 2)
        return row

    def _open_recent(self, item):
        p = item.data(Qt.UserRole)
        if p:
            self.win.open_assembly([p])

    def refresh(self, recents):
        self.recent.clear()
        if not recents:
            it = QListWidgetItem("Nothing yet — your latest assemblies will appear here.")
            it.setFlags(Qt.NoItemFlags)
            self.recent.addItem(it)
            return
        for p in recents[:12]:
            it = QListWidgetItem(theme.icon("assembly", 20, theme.ACCENT, flat=True), Path(p).name)
            it.setToolTip(p)
            it.setData(Qt.UserRole, p)
            self.recent.addItem(it)
