"""Main window: brand bar with menus and tool state, ribbon with the 8 workflow steps (+ More), project panel, the Start page and the
Workspace, which shows the real Shiny tools of SangeRly (same R code, same options) with a new look. Drop .ab1 files anywhere to start an assembly."""
import html
import time
from datetime import datetime
from pathlib import Path

import json

from PySide6.QtCore import QEvent, QSettings, Qt, QUrl
from PySide6.QtGui import QAction, QCursor, QDesktopServices, QKeySequence
from PySide6.QtWidgets import (QApplication, QDockWidget, QFileDialog, QFrame, QHBoxLayout, QLabel, QLineEdit, QMainWindow, QMenu,
                               QDialog, QMessageBox, QPushButton, QSplitter, QStackedWidget, QStyle, QStatusBar, QTabBar, QTextEdit, QToolButton,
                               QTreeWidget, QTreeWidgetItem, QVBoxLayout, QWidget)

from . import APP_NAME, __version__, theme
from .navarea import NavArea
from .assembly_dialog import AssemblyDialog, expand_paths
from .dialogs import MSA_METHODS, TREE_METHODS, MsaDialog, TreeDialog, TrimDialog, clean
from .modules import BY_KEY, DIST_ANALYSES, EXTRAS, MODULES, SECTIONS, SUB_LABELS, TAB_TITLES, sub_id
from .pages import StartPage
from .rbridge import app_base, find_project_root
from .settings_dialog import SettingsDialog
from .webhost import Workspace
from .widgets import DockTitle, RibbonButton

NEED_JS = 16          # www/sangerly_embed.js of the R folder must be at least this version (sgQt.status().ver)


class MainWindow(QMainWindow):
    def __init__(self, bridge, host):
        super().__init__()
        self.bridge = bridge
        self.host = host
        self.cfg = QSettings()
        self._step_btns = {}
        self._tree_sig = None
        self._tools_state = ("starting", "Tools: starting…", "")
        self.msa_state = None            # what the R tool tells: names in the MSA list, is there an alignment, consensus count
        self._nav_guard = 0.0            # time of the last navigation we started (the page status is stale for a moment after it)
        self._ver_warned = False
        self.setWindowTitle(APP_NAME)
        self.setWindowIcon(theme.app_icon())
        self.resize(1400, 880)
        self.work = Workspace(host, self.log_msg)
        self.work.status.connect(self._status)
        self.work.problem.connect(self._patch_problem)
        self.work.command.connect(self._web_command)
        self.work.connected.connect(self._page_connected)
        self._build_central()
        self._build_docks()
        self._build_brand_bar()
        self._build_status()
        bridge.log.connect(self.log_msg)
        bridge.ready.connect(self._r_ready)
        bridge.failed.connect(lambda why: self.log_msg("[R] " + why))
        bridge.stopped.connect(lambda code: self.log_msg(f"[R] stopped (code {code})"))
        bridge.busy.connect(self._busy)
        host.log.connect(self.log_msg)
        host.starting.connect(lambda: self._pill("starting", "Tools: starting…", "Starting the Shiny tools"))
        host.ready.connect(lambda _u: self._pill("ready", "Tools ready", self._ready_tip()))
        host.failed.connect(lambda why: self._pill("error", "Tools: not available", why + "\nClick to open the settings"))
        self.start.refresh(self.cfg.value("recent", [], type=list))
        g = self.cfg.value("win/geometry")
        if g:
            self.restoreGeometry(g)
        QApplication.instance().installEventFilter(self)       # .ab1 files dropped on the web page are caught here

    # ---- helpers -------------------------------------------------------------------------------------------------------------
    def logo_path(self, m):
        root = find_project_root()
        return str(Path(root) / "www" / "img" / m.logo) if root and m.logo else ""

    def _step_icon(self, m, size=38):
        return theme.logo_icon(self.logo_path(m), size) or theme.icon(m.glyph, size, m.color)

    # ---- layout ----------------------------------------------------------------------------------------------------------------
    def _build_central(self):
        box = QWidget()
        v = QVBoxLayout(box)
        v.setContentsMargins(0, 0, 0, 0)
        v.setSpacing(0)
        v.addWidget(self._build_ribbon())
        # document tabs: every tool or section you open gets its own tab; all tabs show the same Shiny
        # session, so the data (assembly, alignment, trees) are shared between them
        self.start = StartPage(self)
        self.doc_bar = QTabBar()
        self.doc_bar.setDocumentMode(True)
        self.doc_bar.setTabsClosable(True)
        self.doc_bar.setMovable(False)
        self.doc_bar.setExpanding(False)
        self.doc_bar.setUsesScrollButtons(True)
        self.doc_bar.setElideMode(Qt.ElideRight)
        self.doc_bar.addTab(theme.icon("home", 18, theme.MUTED, flat=True), "Start page")
        self.doc_bar.setTabData(0, {"key": "home"})
        self.doc_bar.setTabButton(0, QTabBar.ButtonPosition.RightSide, None)
        v.addWidget(self.doc_bar)
        line = QFrame()
        line.setFixedHeight(1)
        line.setStyleSheet(f"background:{theme.BORDER};")
        v.addWidget(line)
        self.stack = QStackedWidget()
        self.stack.addWidget(self.start)
        self.stack.addWidget(self.work)
        v.addWidget(self.stack, 1)
        self.doc_bar.setContextMenuPolicy(Qt.CustomContextMenu)
        self.doc_bar.customContextMenuRequested.connect(self._doc_menu)
        self.doc_bar.currentChanged.connect(self._doc_changed)
        self.doc_bar.tabCloseRequested.connect(self._doc_close)
        self.setCentralWidget(box)

    # ---- document tabs -------------------------------------------------------------------------------------------------------------
    def _doc_icon(self, key):
        m = BY_KEY.get(key)
        return self._step_icon(m, 18) if m else theme.icon("open", 18, theme.ACCENT, flat=True)

    def _doc_find(self, key, sid):
        for i in range(1, self.doc_bar.count()):
            d = self.doc_bar.tabData(i) or {}
            if d.get("key") == key and (d.get("sub") or None) == (sid or None):
                return i
        return -1

    def _doc_title(self, key, sid):
        base = TAB_TITLES.get(key, key)
        lab = SUB_LABELS.get((key, sid)) if sid else None
        return f"{base} – {lab}" if lab else base

    def _doc_changed(self, i):
        if i <= 0:
            self.stack.setCurrentWidget(self.start)
            return
        self.stack.setCurrentWidget(self.work)
        self._navigate(i)

    def _navigate(self, i):
        d = self.doc_bar.tabData(i) or {}
        if d.get("key"):
            self._nav_guard = time.time()
            self.work.goto(d["key"], d.get("nav"))
            if d.get("file"):
                self.work.call("openPaths", [d["file"]])      # shows this file again (it is already open: only selected)

    def _doc_close(self, i):
        if i > 0:
            self.doc_bar.removeTab(i)
            if self.doc_bar.count() <= 1:
                self.stack.setCurrentWidget(self.start)

    def close_other_tabs(self, keep=None):
        """Close every document tab except one (the current one unless `keep` is given)."""
        keep = self.doc_bar.currentIndex() if keep is None or isinstance(keep, bool) else keep
        if keep <= 0:
            return
        for i in range(self.doc_bar.count() - 1, 0, -1):
            if i != keep:
                self.doc_bar.removeTab(i)

    def _doc_menu(self, pos):
        i = self.doc_bar.tabAt(pos)
        if i <= 0:
            return
        m = QMenu(self)
        m.addAction("Close tab", lambda: self._doc_close(i))
        m.addAction("Close all other tabs", lambda: self.close_other_tabs(i))
        m.addAction("Close all tabs", self.close_all_tabs)
        m.exec(self.doc_bar.mapToGlobal(pos))

    def close_all_tabs(self):
        for i in range(self.doc_bar.count() - 1, 0, -1):
            self.doc_bar.removeTab(i)
        self.show_start()

    def _build_ribbon(self):
        fr = QFrame()
        fr.setObjectName("toolbar")
        h = QHBoxLayout(fr)
        h.setContentsMargins(14, 6, 14, 4)
        h.setSpacing(0)

        def sep():
            s = QFrame()
            s.setObjectName("ribbonSep")
            s.setFixedSize(1, 58)
            h.addSpacing(10)
            h.addWidget(s, 0, Qt.AlignTop)
            h.addSpacing(10)

        def group(caption, buttons):
            box = QWidget()
            v = QVBoxLayout(box)
            v.setContentsMargins(0, 0, 0, 0)
            v.setSpacing(1)
            row = QHBoxLayout()
            row.setSpacing(2)
            for b in buttons:
                row.addWidget(b)
            v.addLayout(row)
            cap = QLabel(caption)
            cap.setObjectName("groupCaption")
            cap.setAlignment(Qt.AlignCenter)
            v.addWidget(cap)
            h.addWidget(box)

        def btn(ico, text, tip, slot, step=None):
            b = RibbonButton(ico, text, tip, step, None)
            if slot:
                b.clicked.connect(slot)
            return b

        group("FILE", [btn(theme.icon("home", 38, "#64748b"), "START", "Start page (open files: File > Open files, Ctrl+O)", self.show_start)])
        sep()
        steps = []
        for m in MODULES:
            b = btn(self._step_icon(m), m.short, m.title, None)
            menu = QMenu(b)
            self._fill_tool_menu(menu, m)
            b.setMenu(menu)
            b.setPopupMode(QToolButton.InstantPopup)
            self._step_btns[m.key] = b
            steps.append(b)
        group("WORKFLOW", steps)
        sep()
        more = btn(theme.icon("folder", 38, "#334155"), "MORE", "Open and view files, project and report, GenBank export, tutorial, help", None)
        menu = QMenu(more)
        for key, title, tip in EXTRAS:
            a = QAction(title, self)
            a.setToolTip(tip)
            a.triggered.connect(lambda _=False, k=key: self.open_module(k))
            menu.addAction(a)
        more.setMenu(menu)
        more.setPopupMode(QToolButton.InstantPopup)
        group("TOOLS", [more])
        h.addStretch(1)
        sep()
        group("APPLICATION", [btn(theme.icon("gear", 38, "#475569"), "SETTINGS", "R engine and folders", self.settings_dialog)])
        return fr

    def _build_brand_bar(self):
        bar = QFrame()
        bar.setObjectName("topBar")
        bar.setFixedHeight(46)
        h = QHBoxLayout(bar)
        h.setContentsMargins(14, 0, 14, 0)
        h.setSpacing(10)
        logo = QLabel()
        logo.setPixmap(theme.logo_pixmap(28).scaled(28, 28, Qt.KeepAspectRatio, Qt.SmoothTransformation))
        logo.setFixedSize(28, 28)
        h.addWidget(logo)
        name = QLabel(APP_NAME)
        name.setObjectName("brandName")
        h.addWidget(name)
        h.addSpacing(14)
        from PySide6.QtWidgets import QMenuBar
        self.mb = QMenuBar()
        self.mb.setObjectName("brandMenu")
        self.mb.setNativeMenuBar(False)
        self._build_menus(self.mb)
        h.addWidget(self.mb)
        h.addStretch(1)
        self.pill = QPushButton("●  Tools: starting…")
        self.pill.setObjectName("enginePill")
        self.pill.setCursor(Qt.PointingHandCursor)
        self.pill.setProperty("state", "starting")
        self.pill.clicked.connect(self.settings_dialog)
        h.addWidget(self.pill)
        self.setMenuWidget(bar)

    def _build_docks(self):
        box = QWidget()
        v = QVBoxLayout(box)
        v.setContentsMargins(8, 8, 8, 8)
        v.setSpacing(8)
        self.proj_filter = QLineEdit()
        self.proj_filter.setPlaceholderText("Filter project…")
        self.proj_filter.setClearButtonEnabled(True)
        self.proj_filter.textChanged.connect(self._filter_project)
        v.addWidget(self.proj_filter)
        self.proj = QTreeWidget()
        self.proj.setHeaderHidden(True)
        self.proj.setFrameShape(QFrame.NoFrame)
        self.proj.itemDoubleClicked.connect(self._proj_open)
        # Navigation area (like the one of a sequence workbench): the folders of the person, with import, rename, delete
        roots0 = [Path(find_project_root() or ".") / "Data_example"] if find_project_root() else []
        self.nav = NavArea(self.cfg, lambda k: (self._step_icon(BY_KEY[k], 16) if k in BY_KEY else None), roots0)
        self.nav.open_files.connect(self.nav_open_files)
        self.nav.open_session.connect(self.nav_open_session)
        nav_box = QWidget()
        nv = QVBoxLayout(nav_box)
        nv.setContentsMargins(0, 4, 0, 0)
        nv.setSpacing(4)
        cap = QLabel("NAVIGATION AREA")
        cap.setObjectName("panelCaption")
        cap.setStyleSheet("font-size:10px;font-weight:700;letter-spacing:1px;color:#64748b;padding:2px 4px;")
        nv.addWidget(cap)
        nv.addWidget(self.nav, 1)
        wf_box = QWidget()
        wv = QVBoxLayout(wf_box)
        wv.setContentsMargins(0, 0, 0, 0)
        wv.setSpacing(4)
        wv.addWidget(self.proj, 1)
        split = QSplitter(Qt.Vertical)
        split.setChildrenCollapsible(False)
        split.addWidget(wf_box)
        split.addWidget(nav_box)
        split.setStretchFactor(0, 2)
        split.setStretchFactor(1, 3)
        split.setSizes([230, 340])
        v.addWidget(split, 1)
        self._nav_data = None
        self._fill_project(set(), [], None)
        self.dock_proj = QDockWidget("Project", self)
        self.dock_proj.setWidget(box)
        self.dock_proj.setTitleBarWidget(DockTitle("Project", self.dock_proj))
        self.dock_proj.setFeatures(QDockWidget.DockWidgetClosable | QDockWidget.DockWidgetMovable)
        self.addDockWidget(Qt.LeftDockWidgetArea, self.dock_proj)

        self.logw = QTextEdit()
        self.logw.setObjectName("log")
        self.logw.setReadOnly(True)
        self.logw.document().setMaximumBlockCount(5000)
        self.logw.setFont(theme.mono_font(9))
        self.dock_log = QDockWidget("Log", self)
        self.dock_log.setWidget(self.logw)
        self.dock_log.setTitleBarWidget(DockTitle("Log", self.dock_log))
        self.dock_log.setFeatures(QDockWidget.DockWidgetClosable | QDockWidget.DockWidgetMovable)
        self.addDockWidget(Qt.BottomDockWidgetArea, self.dock_log)
        self.dock_proj.setMinimumWidth(220)
        self.logw.setMinimumHeight(70)
        self.resizeDocks([self.dock_proj], [240], Qt.Horizontal)
        self.resizeDocks([self.dock_log], [110], Qt.Vertical)

    def _act(self, text, slot=None, kind=None, shortcut=None):
        a = QAction(text, self)
        if kind:
            a.setIcon(theme.icon(kind, 32, theme.ACCENT, flat=True))
        if shortcut:
            a.setShortcut(QKeySequence(shortcut))
        if slot:
            a.triggered.connect(slot)
        return a

    def _fill_tool_menu(self, menu, mod):
        """Drop-down of one tool : open it, its dialogs, then every section of the tool."""
        a = self._act(f"Open {mod.title}", lambda _=False, k=mod.key: self.open_module(k))
        f = a.font()
        f.setBold(True)
        a.setFont(f)
        menu.addAction(a)
        if mod.key == "assembly":
            menu.addAction(self._act("Sanger assembly…", lambda: self.open_assembly(), "new"))
        if mod.key == "msa":
            # Tools > Multiple sequence alignment > Align with ... opens a dialog (input, method, options)
            sub = menu.addMenu("Multiple sequence alignment")
            for text, code in MSA_METHODS:
                sub.addAction(self._act(text + "…", lambda _=False, c=code: self.open_msa_dialog(c)))
            menu.addAction(self._act("Trim the alignment…", self.open_trim_dialog))
        if mod.key == "barcode":
            # Tools > Genetic distances > Compute ... chooses the analysis; the page then asks for the groups and the settings
            sub = menu.addMenu("Distance")
            for text, code in DIST_ANALYSES:
                sub.addAction(self._act(text + "…", lambda _=False, c=code: self.open_distance(c)))
        if mod.key == "tree":
            sub = menu.addMenu("Build a tree")
            for text, code in TREE_METHODS:
                sub.addAction(self._act(text + "…", lambda _=False, c=code: self.open_tree_dialog(c)))
        menu.addSeparator()
        for e in SECTIONS.get(mod.key, []):
            if e == "-":
                menu.addSeparator()
            else:
                menu.addAction(self._act(e[0], lambda _=False, k=mod.key, sub=e[1]: self.open_module(k, sub)))

    def _build_menus(self, mb):
        m = mb.addMenu("&File")
        m.addAction(self._act("New assembly…", lambda: self.open_assembly(), "new", "Ctrl+N"))
        m.addAction(self._act("Open files…", self.file_open, "open", "Ctrl+O"))
        m.addAction(self._act("Open folder…", self.file_open_folder))
        self.m_recent = m.addMenu("Recent assemblies")
        self._fill_recent_menu()
        m.addSeparator()
        mi = m.addMenu("Import")
        mi.addAction(self._act("Sequences from files into the open alignment…", self.import_files))
        mi.addAction(self._act("Sequences from the clipboard", self.import_clipboard))
        mi.addAction(self._act("Type or paste sequences…", lambda: self.open_cmd("add_seq")))
        mi.addSeparator()
        mi.addAction(self._act("Assemble .ab1 files…", self.open_ab1))
        m.addSeparator()
        m.addAction(self._act("Save", self.file_save, "save", "Ctrl+S"))
        m.addAction(self._act("Save as…", lambda: self.open_cmd("save_as"), None, "Ctrl+Shift+S"))
        m.addAction(self._act("Undo the last edit", lambda: self.open_cmd("undo")))
        m.addAction(self._act("Send to the MSA tab", lambda: self.open_cmd("to_msa")))
        m.addSeparator()
        m.addAction(self._act("Rename the file…", lambda: self.open_cmd("rename")))
        m.addAction(self._act("Delete the file from disk…", lambda: self.open_cmd("del_disk")))
        m.addSeparator()
        m.addAction(self._act("Save assembly results…", lambda: self.open_module("assembly", "Consensus and export")))
        m.addAction(self._act("Save or open the project…", lambda: self.open_module("project")))
        m.addAction(self._act("GenBank export…", lambda: self.open_module("genbank")))
        m.addSeparator()
        m.addAction(self._act("Open the example data folder", self.open_examples))
        m.addSeparator()
        m.addAction(self._act("Close tab", lambda: self._doc_close(self.doc_bar.currentIndex()), None, "Ctrl+W"))
        m.addAction(self._act("Close all other tabs", self.close_other_tabs))
        m.addAction(self._act("Close all tabs", self.close_all_tabs, None, "Ctrl+Shift+W"))
        m.addSeparator()
        m.addAction(self._act("Exit", self.close))

        m = mb.addMenu("&Actions")
        m.addAction(self._act("Sanger assembly…", lambda: self.open_assembly(), "assembly"))
        m.addSeparator()
        m.addAction(self._act("BLAST the consensus", lambda: self.open_module("blast")))
        m.addAction(self._act("Align sequences…", lambda: self.open_msa_dialog("ClustalOmega")))
        m.addAction(self._act("Construct a Neighbor-Joining tree…", lambda: self.open_tree_dialog("nj")))
        m.addAction(self._act("Construct a Maximum-Likelihood tree (IQ-TREE)…", lambda: self.open_tree_dialog("ml")))
        m.addAction(self._act("Find the best substitution model (ModelFinder)…", lambda: self.open_tree_dialog("mf")))
        m.addAction(self._act("Design primers", lambda: self.open_module("primer", "design")))

        m = mb.addMenu("&Tools")
        for mod in MODULES:
            sub = m.addMenu(self._step_icon(mod, 20), mod.title)
            self._fill_tool_menu(sub, mod)
        m.addSeparator()
        for key, title, _tip in EXTRAS:
            m.addAction(self._act(title, lambda _=False, k=key: self.open_module(k)))

        m = mb.addMenu("&Settings")
        m.addAction(self._act("Application settings…", self.settings_dialog, "gear"))
        m.addAction(self._act("Restart the tools", self.restart_tools))

        m = mb.addMenu("&Window")
        m.addAction(self._act("Start page", self.show_start))
        m.addAction(self._act("Close all tabs", self.close_all_tabs))
        self.act_tabs = QAction("Show the tabs inside the tools", self)
        self.act_tabs.setCheckable(True)
        self.act_tabs.setChecked(self.cfg.value("view/inner_tabs", False, type=bool))
        self.act_tabs.toggled.connect(self._inner_tabs_toggled)
        m.addAction(self.act_tabs)
        self.act_classic = QAction("Classic panels (input and run on the page, as in Shiny)", self)
        self.act_classic.setCheckable(True)
        self.act_classic.setChecked(self.cfg.value("view/classic", False, type=bool))
        self.act_classic.toggled.connect(self._classic_toggled)
        m.addAction(self.act_classic)
        m.addSeparator()
        m.addAction(self.dock_proj.toggleViewAction())
        m.addAction(self.dock_log.toggleViewAction())
        m.addSeparator()
        m.addAction(self._act("Reload the workspace", lambda: self.work.view.reload()))
        m.addAction(self._act("Open the workspace in the browser", self.open_in_browser))

        m = mb.addMenu("&Help")
        m.addAction(self._act("Tutorial", lambda: self.open_module("tutorial")))
        m.addAction(self._act("Help", lambda: self.open_module("help"), None, "F1"))
        m.addAction(self._act("About and citation", lambda: self.open_module("about")))
        m.addAction(self._act("About this window", self.about))

    def _build_status(self):
        sb = QStatusBar()
        sb.setSizeGripEnabled(False)
        self.setStatusBar(sb)
        self.task_lbl = QLabel("No active tasks")
        sb.addWidget(self.task_lbl, 1)
        logo = QLabel()
        pm = theme.logo_pixmap(18)
        if not pm.isNull():
            logo.setPixmap(pm.scaled(18, 18, Qt.KeepAspectRatio, Qt.SmoothTransformation))
        logo.setFixedSize(20, 20)
        logo.setStyleSheet("background: transparent; border: none;")
        sb.addPermanentWidget(logo)
        sb.addPermanentWidget(QLabel(f"{APP_NAME} {__version__}"))

    # ---- state -----------------------------------------------------------------------------------------------------------------
    def log_msg(self, s):
        low = s.lower()
        col = theme.BAD if any(w in low for w in ("error", "not found", "failed", "stopped")) else theme.WARN if "warn" in low else theme.MUTED if s.startswith("[task]") else theme.TEXT
        ts = datetime.now().strftime("%H:%M:%S")
        self.logw.append(f"<span style='color:#98a2b3'>{ts}</span>&nbsp; <span style='color:{col}'>{html.escape(s).replace(chr(10), '<br>')}</span>")

    def _pill(self, state, text, tip=""):
        self._tools_state = (state, text, tip)
        self.pill.setText("●  " + text)
        self.pill.setToolTip(tip or "Click to open the settings")
        self.pill.setProperty("state", state)
        theme.repolish(self.pill)

    def _ready_tip(self):
        return (self._r_version + "\n" if getattr(self, "_r_version", "") else "") + "Click to open the settings"

    def _r_ready(self, msg):
        self._r_version = msg.get("r", "")
        self.log_msg(f"[R] ready ({self._r_version}, pid {msg.get('pid')})")
        if self._tools_state[0] == "ready":
            self._pill("ready", "Tools ready", self._ready_tip())

    def _web_command(self, cmd):
        """Buttons inside the web page (sgQt.cmd) that need a Qt window."""
        if cmd == "assembly":
            self.open_assembly()
        elif cmd == "assembly-ref":
            self.open_assembly(mode="reference")

    def _patch_problem(self, msg):
        self._pill("error", "R code needs the Qt patch", msg)
        if not getattr(self, "_patch_warned", False):
            self._patch_warned = True
            QMessageBox.warning(self, "SangeRly", msg)

    def _busy(self, b):
        self.task_lbl.setText("Working…" if b else "No active tasks")

    def settings_dialog(self):
        dlg = SettingsDialog(self.bridge, self)
        if dlg.exec() == SettingsDialog.Accepted:
            self.bridge.restart()
            self.host.restart()

    def restart_tools(self):
        self.bridge.restart()
        self.host.restart()

    # ---- status of the Shiny workspace (every second) -----------------------------------------------------------------------
    def _status(self, st):
        done = set(st.get("done") or [])
        for k, b in self._step_btns.items():
            b.set_done(k in done)
        if "msa" in st:
            self.msa_state = st.get("msa")
        ver = st.get("ver")
        if st.get("ready") and (ver or 0) < NEED_JS and not self._ver_warned:
            self._ver_warned = True
            self._patch_problem(self.work._patch_msg())
        # the page moved by itself (for example "Send to MSA" in BLAST): follow it with a tab
        tab, sub = st.get("tab"), st.get("sub")
        cur = self.doc_bar.currentIndex()
        if cur > 0 and tab and tab != "home" and time.time() - self._nav_guard > 2.5:
            d = self.doc_bar.tabData(cur) or {}
            if tab != d.get("key"):
                self.open_module(tab)
            elif sub and d.get("sub") and sub != d.get("sub") and (tab, sub) in SUB_LABELS:
                self.open_module(tab, sub)
        nav = st.get("nav")
        nav_sig = json.dumps(nav, sort_keys=True) if nav else None
        if nav_sig != getattr(self, "_nav_sig", None):
            self._nav_sig = nav_sig
            self.nav.set_session(nav)
        samples = st.get("samples") or []
        sig = (tuple(sorted(done)), tuple(samples))
        if sig != self._tree_sig:
            self._tree_sig = sig
            self._fill_project(done, samples, st.get("cur"))

    def _classic_toggled(self, on):
        self.cfg.setValue("view/classic", bool(on))
        self.work.js(f"window.sgQt && sgQt.classic && sgQt.classic({'true' if on else 'false'});")

    def _page_connected(self):
        self.work.js(f"window.sgQt && sgQt.classic && sgQt.classic({'true' if self.act_classic.isChecked() else 'false'});")
        self.work.js(f"window.sgQt && sgQt.showTabs && sgQt.showTabs({'true' if self.act_tabs.isChecked() else 'false'});")

    def _inner_tabs_toggled(self, on):
        self.cfg.setValue("view/inner_tabs", bool(on))
        self.work.js(f"window.sgQt && sgQt.showTabs && sgQt.showTabs({'true' if on else 'false'});")

    def _fill_project(self, done, samples, cur):
        self.proj.clear()
        root = QTreeWidgetItem(["Workflow"])
        root.setIcon(0, theme.icon("assembly", 20, theme.ACCENT, flat=True))
        self.proj.addTopLevelItem(root)
        for m in MODULES:
            it = QTreeWidgetItem([("✓ " if m.key in done else "") + m.title])
            it.setIcon(0, self._step_icon(m, 20))
            it.setData(0, Qt.UserRole, ("step", m.key))
            root.addChild(it)
            if m.key == "assembly":
                for s in samples:
                    c = QTreeWidgetItem([s])
                    c.setData(0, Qt.UserRole, ("sample", s))
                    if s == cur:
                        f = c.font(0)
                        f.setBold(True)
                        c.setFont(0, f)
                    it.addChild(c)
                it.setExpanded(bool(samples))
        root.setExpanded(True)
        self._filter_project(self.proj_filter.text() if hasattr(self, "proj_filter") else "")

    def _proj_open(self, item, _col):
        kind = item.data(0, Qt.UserRole)
        if not kind:
            return
        if kind[0] == "sample":
            self.open_module("assembly")
            self.work.js(f"window.sgQt && sgQt.pickSample({__import__('json').dumps(kind[1])});")
        else:
            self.open_module(kind[1])

    def nav_open_session(self, key, sub, name, kind):
        if kind == "sample":
            self.open_module("assembly")
            self.work.js(f"window.sgQt && sgQt.pickSample({json.dumps(name)});")
        else:
            self.open_module(key, sub or None)

    def nav_open_files(self, paths, action):
        """Open files of the Navigation area in the tool that fits them."""
        paths = [str(p) for p in paths]
        if action == "default":
            # what the "Open and view files" tool can show (ab1, FASTA, FASTQ, GenBank, Newick, tables, text): open it there;
            # every other file (Excel, pictures ...) goes to the application of the computer
            known = {".ab1", ".abi", ".fa", ".fasta", ".fas", ".fna", ".fsa", ".seq", ".txt", ".fastq", ".fq", ".gb", ".gbk", ".gbff",
                     ".nwk", ".tre", ".tree", ".treefile", ".csv", ".tsv"}
            ok = [p for p in paths if Path(p).suffix.lower() in known]
            if not ok:
                for p in paths:
                    QDesktopServices.openUrl(QUrl.fromLocalFile(p))
                return
            self.open_file_tabs(ok)
            return
        if action == "assembly":
            self.open_assembly(paths)
        elif action == "msa":
            self.open_msa_dialog("ClustalOmega", files=paths)
        elif action == "distances":
            self.open_module("barcode")
            self.work.call("alnPath", paths[0])

    def _filter_project(self, text):
        if hasattr(self, "nav"):
            self.nav.filter(text)
        t = (text or "").strip().lower()
        root = self.proj.topLevelItem(0)
        if root is None:
            return
        for i in range(root.childCount()):
            a = root.child(i)
            a_hit = not t or t in a.text(0).lower()
            any_child = False
            for j in range(a.childCount()):
                c = a.child(j)
                hit = a_hit or t in c.text(0).lower()
                c.setHidden(not hit)
                any_child = any_child or hit
            a.setHidden(not (a_hit or any_child))

    # ---- actions ---------------------------------------------------------------------------------------------------------------
    def show_start(self):
        self.doc_bar.setCurrentIndex(0)

    def open_module(self, key, sub=None):
        """Open a tool (and a section of it) in a document tab; an already open tab is selected instead of opening a second one."""
        if key == "home":
            self.show_start()
            return
        sid = sub_id(sub)
        i = self._doc_find(key, sid)
        if i < 0:
            i = self.doc_bar.addTab(self._doc_icon(key), self._doc_title(key, sid))
            self.doc_bar.setTabData(i, {"key": key, "sub": sid, "nav": sub})
            self.doc_bar.setTabToolTip(i, self._doc_title(key, sid))
        if self.doc_bar.currentIndex() == i:
            self.stack.setCurrentWidget(self.work)
            self._navigate(i)
        else:
            self.doc_bar.setCurrentIndex(i)

    # ---- File menu: open, import and save the files of the Open tool ----------------------------------------------------------
    FILE_FILTER = ("Sequences and data (*.ab1 *.abi *.fa *.fasta *.fas *.fna *.fsa *.seq *.txt *.fastq *.fq *.gb *.gbk *.gbff "
                   "*.nwk *.tre *.tree *.treefile *.csv *.tsv);;All files (*)")

    def _last_dir(self):
        return self.cfg.value("files/last_dir", str(Path.home()), type=str)

    def _remember_dir(self, p):
        self.cfg.setValue("files/last_dir", str(Path(p).parent if Path(p).is_file() else p))

    def _in_open_tool(self):
        i = self.doc_bar.currentIndex()
        return i > 0 and (self.doc_bar.tabData(i) or {}).get("key") == "open"

    def file_open(self):
        paths, _ = QFileDialog.getOpenFileNames(self, "Open files", self._last_dir(), self.FILE_FILTER)
        if paths:
            self._remember_dir(paths[0])
            self.nav_open_files(paths, "default")

    def file_open_folder(self):
        d = QFileDialog.getExistingDirectory(self, "Open a folder", self._last_dir())
        if not d:
            return
        self._remember_dir(d)
        known = {".ab1", ".abi", ".fa", ".fasta", ".fas", ".fna", ".fsa", ".seq", ".fastq", ".fq", ".gb", ".gbk", ".gbff",
                 ".nwk", ".tre", ".tree", ".treefile", ".csv", ".tsv"}
        files = sorted(str(p) for p in Path(d).iterdir() if p.is_file() and p.suffix.lower() in known)
        if files:
            self.open_file_tabs(files)
        else:
            self.statusBar().showMessage("No sequence or data files in that folder.", 5000)

    def open_cmd(self, cmd):
        """Press a button of the Open tool (save, save_as, undo, to_msa, rename, del_disk, add_seq ...) for the file that is shown."""
        if not self._in_open_tool():
            self.statusBar().showMessage("Open a file first (File > Open files…).", 4000)
            return
        self.work.call("openCmd", cmd)

    def file_save(self):
        if self._in_open_tool():
            self.work.call("openCmd", "save")
        else:
            self.open_module("assembly", "Consensus and export")

    def import_files(self):
        paths, _ = QFileDialog.getOpenFileNames(self, "Import sequences", self._last_dir(),
                                                "Sequences (*.fa *.fasta *.fas *.fna *.fsa *.seq *.txt *.fastq *.fq *.gb *.gbk *.gbff);;All files (*)")
        if not paths:
            return
        self._remember_dir(paths[0])
        if not self._in_open_tool():
            self.open_module("open")
        self.work.call("importSeqs", {"paths": paths})

    def import_clipboard(self):
        txt = QApplication.clipboard().text()
        if not txt.strip():
            self.statusBar().showMessage("The clipboard is empty.", 4000)
            return
        if not self._in_open_tool():
            self.open_module("open")
        self.work.call("importSeqs", {"text": txt})

    def open_file_tabs(self, paths):
        """Every file opened from the Navigation area gets its own document tab (all stay open); the tab shows it in the Open tool."""
        self.work.call("openPaths", paths)               # the Open tool adds them all
        first = None
        for p in paths:
            sid = "file:" + p
            i = self._doc_find("open", sid)
            if i < 0:
                i = self.doc_bar.addTab(self._doc_icon("open"), Path(p).name)
                self.doc_bar.setTabData(i, {"key": "open", "sub": sid, "nav": None, "file": p})
                self.doc_bar.setTabToolTip(i, p)
            first = i if first is None else first
        if first is not None:
            if self.doc_bar.currentIndex() == first:
                self._navigate(first)
            else:
                self.doc_bar.setCurrentIndex(first)

    def open_distance(self, code):
        """Genetic distances page, with the analysis chosen (the page then asks for the taxa and groups and the settings)."""
        self.open_module("barcode")
        self.work.goto("barcode", "analysis:" + code)

    # ---- dialogs ------------------------------------------------------------------------------------------------------------
    def open_msa_dialog(self, code="ClustalOmega", files=None):
        dlg = MsaDialog(self, code, self.msa_state, self.cfg)
        for f in files or []:
            dlg._add(f)
        if dlg.exec() != QDialog.Accepted:
            return
        req = dlg.request()
        self.log_msg(f"[task] MSA: {clean(dlg.method.currentText())}, {len(req['files'])} file(s)"
                     + (", sequences of the list" if req["keep"] else "") + (", Assembly consensus" if req["asm"] else "")
                     + (", pasted FASTA" if req["paste"] else ""))
        self.open_module("msa", "Editor")
        self.work.call("msaRun", req)

    def open_trim_dialog(self):
        st = self.msa_state or {}
        if self.msa_state is not None and not st.get("aln"):
            if QMessageBox.question(self, "Trim the alignment", "There is no alignment yet. Align the sequences first?\n\n"
                                    "(Choose 'I already have an alignment' if your sequences are already aligned.)") == QMessageBox.Yes:
                self.open_msa_dialog("ClustalOmega")
            return
        dlg = TrimDialog(self, self.msa_state, self.cfg)
        if dlg.exec() != QDialog.Accepted:
            return
        self.log_msg("[task] MSA: trimming the alignment")
        self.open_module("msa", "Trimmed result")
        self.work.call("msaTrimRun", dlg.request())

    def open_tree_dialog(self, code="nj"):
        dlg = TreeDialog(self, code, self.msa_state, self.cfg)
        if dlg.exec() != QDialog.Accepted:
            return
        req = dlg.request()
        self.log_msg(f"[task] Tree: {dlg.method.currentText()}")
        self.open_module("tree", req["open"])
        self.work.call("apply", {"set": req["set"], "send": req.get("send", []), "click": req["click"], "fasta": req.get("fasta", "")})

    def open_in_browser(self):
        if self.host.url and self.host.is_ready:
            QDesktopServices.openUrl(QUrl(self.host.url + "/"))
        else:
            QMessageBox.information(self, "Workspace", "The tools are not ready yet.")

    def open_ab1(self):
        ps, _ = QFileDialog.getOpenFileNames(self, "Open .ab1 files", self.cfg.value("assembly/in_dir", ""), "AB1 files (*.ab1);;All files (*)")
        if ps:
            self.open_assembly(ps)

    def open_examples(self):
        for c in [Path(find_project_root() or ".") / "Data_example", app_base() / "Data_example", Path.home() / "Documents" / "Sangerly" / "Data_example"]:
            if c.exists():
                QDesktopServices.openUrl(QUrl.fromLocalFile(str(c)))
                return
        QMessageBox.information(self, "Example data", "The Data_example folder was not found.")

    def open_assembly(self, files=None, mode=None):
        files = files if isinstance(files, list) else None
        dlg = AssemblyDialog(self.bridge, self, files)
        if mode == "reference":
            dlg.rb_ref.setChecked(True)
        if dlg.exec() != AssemblyDialog.Accepted or not dlg.result_request:
            return
        self.run_assembly(dlg.result_request)

    def run_assembly(self, req):
        """The Assembly tab of the workspace does the work (same code as its ASSEMBLE button), then saves the chosen files."""
        self.log_msg(f"[task] Assembly: {len(req['files'])} file(s), mode {req['mode']}" +
                     (f", results go to {req['out_dir']}" if req.get("out_dir") else ""))
        self.open_module("assembly")
        self.work.assemble(req)
        rec = self.cfg.value("recent", [], type=list)
        for f in req["files"][:1]:
            if f["path"] in rec:
                rec.remove(f["path"])
            rec.insert(0, f["path"])
        self.cfg.setValue("recent", rec[:20])
        self.start.refresh(rec)
        self._fill_recent_menu()

    def _fill_recent_menu(self):
        self.m_recent.clear()
        for p in self.cfg.value("recent", [], type=list)[:10]:
            a = QAction(Path(p).name, self)
            a.setToolTip(p)
            a.triggered.connect(lambda _=False, p=p: self.open_assembly([p]))
            self.m_recent.addAction(a)

    def about(self):
        QMessageBox.about(self, f"About {APP_NAME}", f"<b>{APP_NAME}</b> {__version__}<br>Desktop window (Qt / Python) that runs the SangeRly tools "
                          "(R / Shiny) inside it. The tools and their options are those of the Shiny program.")

    # ---- drag and drop: .ab1 files / folders dropped anywhere start an assembly ---------------------------------------------------
    def eventFilter(self, obj, ev):
        if ev.type() in (QEvent.DragEnter, QEvent.DragMove, QEvent.Drop) and isinstance(obj, QWidget) and obj.window() is self:
            md = ev.mimeData()
            if md is not None and md.hasUrls():
                local = [u.toLocalFile() for u in md.urls() if u.isLocalFile()]
                if local:
                    if ev.type() == QEvent.Drop:
                        ps = expand_paths(local, True)
                        if ps:
                            ev.acceptProposedAction()
                            self.open_assembly(ps)
                            return True
                    else:
                        ev.acceptProposedAction()
                        return True
        return super().eventFilter(obj, ev)

    def closeEvent(self, e):
        self.cfg.setValue("win/geometry", self.saveGeometry())
        self.host.stop()
        self.bridge.stop()
        super().closeEvent(e)
