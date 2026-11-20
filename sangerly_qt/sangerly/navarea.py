"""Navigation area: the folders of the person on this computer, like the Navigation Area of a sequence workbench.

The person adds folders (locations) and sees their files in a tree. Files can be imported (copied) into a folder, renamed, deleted (to the
recycle bin), dragged in from the file manager, and opened in the tool that fits them. The last item, "Current session", lists what the
open tools hold now (assembly consensus, alignment, trees ...), read-only."""
import os
import shutil
from pathlib import Path

from PySide6.QtCore import QFile, QFileInfo, QSettings, Qt, QUrl, Signal
from PySide6.QtGui import QDesktopServices
from PySide6.QtWidgets import (QAbstractItemView, QApplication, QFileDialog, QFileIconProvider, QFrame, QInputDialog, QLabel, QMenu, QMessageBox,
                               QStyle, QToolButton, QTreeWidget, QTreeWidgetItem, QVBoxLayout, QHBoxLayout, QWidget)

KIND_AB1 = {".ab1", ".abi"}
KIND_SEQ = {".fa", ".fasta", ".fas", ".fna", ".ffn", ".aln", ".fst", ".fsa"}
KIND_TREE = {".nwk", ".newick", ".tre", ".tree", ".treefile", ".contree"}
ROLE = Qt.UserRole
MAX_FILTER_NODES = 6000


def kind_of(path):
    e = Path(path).suffix.lower()
    return "ab1" if e in KIND_AB1 else "seq" if e in KIND_SEQ else "tree" if e in KIND_TREE else "other"


class NavArea(QWidget):
    open_files = Signal(list, str)            # paths, action: "default" | "assembly" | "msa" | "distances"
    open_session = Signal(str, str, str, str)  # module key, tab title, name, kind

    def __init__(self, cfg: QSettings, icon_for, default_roots=(), parent=None):
        super().__init__(parent)
        self.cfg = cfg
        self.icon_for = icon_for              # module key -> QIcon (or None)
        self.provider = QFileIconProvider()   # the icon the operating system shows for a file (Excel, image, tree viewer ...)
        st = self.style()
        self.ic_dir, self.ic_file = st.standardIcon(QStyle.SP_DirIcon), st.standardIcon(QStyle.SP_FileIcon)
        roots = self.cfg.value("nav/roots", None)
        if roots is None:
            roots = [str(p) for p in default_roots if Path(p).is_dir()]
            self.cfg.setValue("nav/roots", roots)
        self.roots = [r for r in (roots if isinstance(roots, list) else [roots]) if r]
        self._session = None
        v = QVBoxLayout(self)
        v.setContentsMargins(0, 0, 0, 0)
        v.setSpacing(4)
        bar = QHBoxLayout()
        bar.setContentsMargins(2, 0, 2, 0)
        bar.setSpacing(2)
        for icon, tip, slot, txt in (
                (QStyle.SP_DirOpenIcon, "Add a folder to the navigation area", self.add_location, ""),
                (QStyle.SP_FileDialogNewFolder, "New folder (inside the selected folder)", self.new_folder, ""),
                (QStyle.SP_ArrowDown, "Import files into the selected folder", self.import_files, " Import"),
                (QStyle.SP_BrowserReload, "Refresh", self.refresh, "")):
            b = QToolButton()
            b.setIcon(st.standardIcon(icon))
            b.setText(txt)
            b.setToolButtonStyle(Qt.ToolButtonTextBesideIcon if txt else Qt.ToolButtonIconOnly)
            b.setToolTip(tip)
            b.setAutoRaise(True)
            b.clicked.connect(slot)
            bar.addWidget(b)
        bar.addStretch(1)
        v.addLayout(bar)
        self.tree = QTreeWidget()
        self.tree.setHeaderHidden(True)
        self.tree.setFrameShape(QFrame.NoFrame)
        self.tree.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.tree.setAcceptDrops(True)
        self.tree.setDragDropMode(QAbstractItemView.DropOnly)
        self.tree.setContextMenuPolicy(Qt.CustomContextMenu)
        self.tree.customContextMenuRequested.connect(self._menu)
        self.tree.itemExpanded.connect(self._expanded)
        self.tree.itemDoubleClicked.connect(self._double)
        self.tree.dragEnterEvent = self._drag_enter
        self.tree.dragMoveEvent = self._drag_enter
        self.tree.dropEvent = self._drop
        v.addWidget(self.tree, 1)
        self.hint = QLabel("Add a folder with the button above, or drop files here.")
        self.hint.setWordWrap(True)
        self.hint.setStyleSheet("color:#64748b;font-size:11px;padding:2px 6px;")
        v.addWidget(self.hint)
        self.refresh()

    # ---- building the tree ----------------------------------------------------------------------------------------------------
    def refresh(self):
        open_paths = self._expanded_paths()
        self.tree.clear()
        for r in self.roots:
            it = self._dir_item(r, root=True)
            self.tree.addTopLevelItem(it)
            if r in open_paths:
                it.setExpanded(True)
        self._add_session_item()
        self._restore(open_paths)
        self.hint.setVisible(not self.roots)

    def _expanded_paths(self):
        out = set()

        def walk(it):
            d = it.data(0, ROLE)
            if d and d[0] == "dir" and it.isExpanded():
                out.add(d[1])
            for i in range(it.childCount()):
                walk(it.child(i))
        for i in range(self.tree.topLevelItemCount()):
            walk(self.tree.topLevelItem(i))
        return out

    def _restore(self, paths):
        def walk(it):
            d = it.data(0, ROLE)
            if d and d[0] == "dir" and d[1] in paths:
                it.setExpanded(True)
            for i in range(it.childCount()):
                walk(it.child(i))
        for i in range(self.tree.topLevelItemCount()):
            walk(self.tree.topLevelItem(i))

    def _dir_item(self, path, root=False):
        it = QTreeWidgetItem([Path(path).name or path])
        it.setIcon(0, self.ic_dir)
        it.setData(0, ROLE, ("dir", path, root))
        it.setToolTip(0, path)
        if root:
            f = it.font(0)
            f.setBold(True)
            it.setFont(0, f)
        if not os.path.isdir(path):
            it.setForeground(0, Qt.gray)
            it.setToolTip(0, path + "  (not found)")
        else:
            it.addChild(QTreeWidgetItem(["…"]))          # placeholder: the children are read when the folder is opened
            it.setChildIndicatorPolicy(QTreeWidgetItem.ShowIndicator)
        return it

    def _file_item(self, path):
        k = kind_of(path)
        it = QTreeWidgetItem([Path(path).name])
        m = {"ab1": "assembly", "seq": "msa"}.get(k)
        ic = self.icon_for(m) if m else None
        if ic is None:
            ic = self.provider.icon(QFileInfo(path))
        it.setIcon(0, ic if not ic.isNull() else self.ic_file)
        it.setData(0, ROLE, ("file", path, k))
        try:
            it.setToolTip(0, f"{path}\n{os.path.getsize(path) / 1024:.1f} KB")
        except OSError:
            it.setToolTip(0, path)
        return it

    def _load(self, it):
        d = it.data(0, ROLE)
        if not d or d[0] != "dir":
            return
        if it.childCount() == 1 and it.child(0).data(0, ROLE) is None:
            it.takeChild(0)
        elif it.childCount():
            return
        try:
            ents = sorted(os.scandir(d[1]), key=lambda e: (not e.is_dir(), e.name.lower()))
        except OSError:
            return
        for e in ents:
            if e.name.startswith("."):
                continue
            it.addChild(self._dir_item(e.path) if e.is_dir() else self._file_item(e.path))

    def _expanded(self, it):
        self._load(it)

    # ---- session data (read-only) -------------------------------------------------------------------------------------------------
    def set_session(self, nav):
        self._session = nav
        self._add_session_item(replace=True)

    def _add_session_item(self, replace=False):
        if replace:
            for i in range(self.tree.topLevelItemCount()):
                d = self.tree.topLevelItem(i).data(0, ROLE)
                if d and d[0] == "session":
                    was = self.tree.topLevelItem(i).isExpanded()
                    self.tree.takeTopLevelItem(i)
                    break
            else:
                was = False
        else:
            was = False
        root = QTreeWidgetItem(["Current session"])
        root.setIcon(0, self.style().standardIcon(QStyle.SP_ComputerIcon))
        root.setData(0, ROLE, ("session",))
        f = root.font(0)
        f.setBold(True)
        root.setFont(0, f)
        root.setToolTip(0, "What the open tools hold now (read-only)")
        for g in (self._session or {}).get("groups") or []:
            items = g.get("items") or []
            fo = QTreeWidgetItem([g.get("title", "")])
            fo.setIcon(0, self.ic_dir)
            fo.setData(0, ROLE, ("sfolder", g.get("go"), g.get("sub") or "", g.get("title"), ""))
            if not items:
                fo.setForeground(0, Qt.gray)
            ic = self.icon_for(g.get("go"))
            for x in items:
                c = QTreeWidgetItem([x.get("name", "")])
                if ic:
                    c.setIcon(0, ic)
                c.setToolTip(0, x.get("name", "") + (f"  ({x['info']})" if x.get("info") else ""))
                c.setData(0, ROLE, ("sitem", x.get("go"), x.get("sub") or "", x.get("name"), x.get("kind") or ""))
                fo.addChild(c)
            root.addChild(fo)
            fo.setExpanded(bool(items) and g.get("key") != "seqs")
        self.tree.addTopLevelItem(root)
        root.setExpanded(was or not replace)

    # ---- the files: folders, import, new folder -------------------------------------------------------------------------------------
    def selected_dir(self):
        """The folder that new things go into: the selected folder, the folder of the selected file, else the first location."""
        for it in self.tree.selectedItems():
            d = it.data(0, ROLE)
            if d and d[0] == "dir":
                return d[1]
            if d and d[0] == "file":
                return str(Path(d[1]).parent)
        return self.roots[0] if self.roots else None

    def add_location(self):
        p = QFileDialog.getExistingDirectory(self, "Add a folder to the navigation area", self.cfg.value("nav/last", str(Path.home())))
        if not p:
            return
        p = os.path.normpath(p)
        self.cfg.setValue("nav/last", p)
        if p not in self.roots:
            self.roots.append(p)
            self.cfg.setValue("nav/roots", self.roots)
        self.refresh()

    def remove_location(self, path):
        self.roots = [r for r in self.roots if r != path]
        self.cfg.setValue("nav/roots", self.roots)
        self.refresh()

    def new_folder(self):
        base = self.selected_dir()
        if not base:
            QMessageBox.information(self, "New folder", "Add a folder to the navigation area first.")
            return
        name, ok = QInputDialog.getText(self, "New folder", f"Name of the new folder in\n{base}")
        if ok and name.strip():
            try:
                os.makedirs(os.path.join(base, name.strip()), exist_ok=False)
            except OSError as e:
                QMessageBox.warning(self, "New folder", str(e))
            self._reload_dir(base)

    def import_files(self, base=None):
        base = base if isinstance(base, str) else self.selected_dir()
        if not base:
            QMessageBox.information(self, "Import", "Add a folder to the navigation area first, then import files into it.")
            return
        ps, _ = QFileDialog.getOpenFileNames(self, f"Import files into {Path(base).name}", self.cfg.value("nav/import", str(Path.home())),
                                             "Sequence and data files (*.ab1 *.abi *.fa *.fasta *.fas *.fna *.aln *.nwk *.tre *.tree *.csv *.tsv *.xlsx *.txt);;All files (*)")
        if ps:
            self.cfg.setValue("nav/import", str(Path(ps[0]).parent))
            self._copy_in(ps, base)

    def import_folder(self, base):
        p = QFileDialog.getExistingDirectory(self, "Import a folder", self.cfg.value("nav/import", str(Path.home())))
        if p:
            self._copy_in([p], base)

    def _unique(self, base, name):
        t = os.path.join(base, name)
        if not os.path.exists(t):
            return t
        stem, ext = os.path.splitext(name)
        i = 1
        while os.path.exists(os.path.join(base, f"{stem}_{i}{ext}")):
            i += 1
        return os.path.join(base, f"{stem}_{i}{ext}")

    def _copy_in(self, paths, base):
        n, err = 0, []
        for p in paths:
            try:
                t = self._unique(base, os.path.basename(p.rstrip("/\\")))
                if os.path.isdir(p):
                    shutil.copytree(p, t)
                else:
                    shutil.copy2(p, t)
                n += 1
            except (OSError, shutil.Error) as e:
                err.append(f"{os.path.basename(p)}: {e}")
        if err:
            QMessageBox.warning(self, "Import", "\n".join(err[:6]))
        self._reload_dir(base)
        return n

    def _reload_dir(self, path):
        """Read a folder again and keep it open."""
        def find(it):
            d = it.data(0, ROLE)
            if d and d[0] == "dir" and d[1] == path:
                return it
            for i in range(it.childCount()):
                r = find(it.child(i))
                if r:
                    return r
        for i in range(self.tree.topLevelItemCount()):
            it = find(self.tree.topLevelItem(i))
            if it:
                while it.childCount():
                    it.takeChild(0)
                self._load(it)
                it.setExpanded(True)
                return
        self.refresh()

    # ---- rename, delete, open ----------------------------------------------------------------------------------------------------------
    def rename(self, path):
        old = os.path.basename(path)
        name, ok = QInputDialog.getText(self, "Rename", "New name:", text=old)
        name = name.strip()
        if not ok or not name or name == old:
            return
        t = os.path.join(os.path.dirname(path), name)
        if os.path.exists(t):
            QMessageBox.warning(self, "Rename", f"'{name}' already exists in this folder.")
            return
        try:
            os.rename(path, t)
        except OSError as e:
            QMessageBox.warning(self, "Rename", str(e))
            return
        if path in self.roots:
            self.roots = [t if r == path else r for r in self.roots]
            self.cfg.setValue("nav/roots", self.roots)
            self.refresh()
        else:
            self._reload_dir(os.path.dirname(path))

    def delete(self, paths):
        names = "\n".join(os.path.basename(p) for p in paths[:8]) + ("\n…" if len(paths) > 8 else "")
        if QMessageBox.question(self, "Delete", f"Move {len(paths)} item(s) to the recycle bin?\n\n{names}") != QMessageBox.Yes:
            return
        parents = set()
        for p in paths:
            ok = QFile.moveToTrash(p)
            if not ok:
                if QMessageBox.question(self, "Delete", f"The recycle bin is not available for\n{p}\n\nDelete it permanently?") == QMessageBox.Yes:
                    try:
                        shutil.rmtree(p) if os.path.isdir(p) else os.remove(p)
                        ok = True
                    except OSError as e:
                        QMessageBox.warning(self, "Delete", str(e))
            if ok:
                parents.add(os.path.dirname(p))
                if p in self.roots:
                    self.roots = [r for r in self.roots if r != p]
                    self.cfg.setValue("nav/roots", self.roots)
        self.refresh()

    def _double(self, it, _c):
        d = it.data(0, ROLE)
        if not d:
            return
        if d[0] == "file":
            self.open_files.emit([d[1]], "default")
        elif d[0] in ("sitem", "sfolder") and d[1]:
            self.open_session.emit(d[1], d[2], d[3] or "", d[4] or "")

    def _paths(self):
        out = []
        for it in self.tree.selectedItems():
            d = it.data(0, ROLE)
            if d and d[0] in ("file", "dir"):
                out.append(d[1])
        return out

    def _menu(self, pos):
        it = self.tree.itemAt(pos)
        d = it.data(0, ROLE) if it else None
        m = QMenu(self)
        if not d:
            m.addAction("Add a folder…", self.add_location)
            m.exec(self.tree.viewport().mapToGlobal(pos))
            return
        if d[0] in ("sitem", "sfolder"):
            if d[1]:
                m.addAction("Open", lambda: self.open_session.emit(d[1], d[2], d[3] or "", d[4] or ""))
            if d[0] == "sitem":
                m.addAction("Copy the name", lambda: QApplication.clipboard().setText(d[3] or ""))
            m.exec(self.tree.viewport().mapToGlobal(pos))
            return
        if d[0] == "session":
            return
        if not it.isSelected():
            self.tree.clearSelection()
            it.setSelected(True)
        sel = [x.data(0, ROLE) for x in self.tree.selectedItems()]
        files = [x[1] for x in sel if x and x[0] == "file"]
        kinds = {kind_of(f) for f in files}
        if d[0] == "file":
            m.addAction("Open", lambda: self.open_files.emit(files, "default"))
            m.addSeparator()
            if kinds == {"ab1"}:
                m.addAction(f"Assemble {len(files)} .ab1 file(s)…", lambda: self.open_files.emit(files, "assembly"))
            if kinds == {"seq"}:
                m.addAction("Align in MSA…", lambda: self.open_files.emit(files, "msa"))
                if len(files) == 1:
                    m.addAction("Open in Genetic distances (aligned file)", lambda: self.open_files.emit(files, "distances"))
            m.addAction("Open with the system application", lambda: [QDesktopServices.openUrl(QUrl.fromLocalFile(f)) for f in files])
            m.addSeparator()
        else:
            ab1 = self._ab1_in(d[1])
            if ab1:
                m.addAction(f"Assemble the {len(ab1)} .ab1 file(s) of this folder…", lambda: self.open_files.emit(ab1, "assembly"))
                m.addSeparator()
            m.addAction("Import files here…", lambda: self.import_files(d[1]))
            m.addAction("Import a folder here…", lambda: self.import_folder(d[1]))
            m.addAction("New folder", lambda: (self.tree.clearSelection(), it.setSelected(True), self.new_folder()))
            m.addSeparator()
        if len(sel) == 1:
            m.addAction("Rename…", lambda: self.rename(d[1]))
        m.addAction("Delete…", lambda: self.delete(self._paths()))
        m.addSeparator()
        m.addAction("Show in the file manager", lambda: QDesktopServices.openUrl(QUrl.fromLocalFile(d[1] if d[0] == "dir" else os.path.dirname(d[1]))))
        m.addAction("Copy the path", lambda: QApplication.clipboard().setText(d[1]))
        if d[0] == "dir":
            m.addAction("Refresh", lambda: self._reload_dir(d[1]))
            if d[2]:
                m.addAction("Remove from the navigation area", lambda: self.remove_location(d[1]))
        m.exec(self.tree.viewport().mapToGlobal(pos))

    @staticmethod
    def _ab1_in(folder):
        try:
            return sorted(os.path.join(folder, f) for f in os.listdir(folder) if Path(f).suffix.lower() in KIND_AB1)
        except OSError:
            return []

    # ---- drag in from the file manager ----------------------------------------------------------------------------------------------
    def _drag_enter(self, ev):
        if ev.mimeData().hasUrls():
            ev.acceptProposedAction()
        else:
            ev.ignore()

    def _drop(self, ev):
        urls = [u.toLocalFile() for u in ev.mimeData().urls() if u.isLocalFile()]
        if not urls:
            return
        it = self.tree.itemAt(ev.position().toPoint())
        d = it.data(0, ROLE) if it else None
        base = d[1] if d and d[0] == "dir" else os.path.dirname(d[1]) if d and d[0] == "file" else (self.roots[0] if self.roots else None)
        if not base:
            if QMessageBox.question(self, "Navigation area", "Add the dropped folder to the navigation area?") == QMessageBox.Yes:
                for u in urls:
                    if os.path.isdir(u) and u not in self.roots:
                        self.roots.append(u)
                self.cfg.setValue("nav/roots", self.roots)
                self.refresh()
            ev.acceptProposedAction()
            return
        self._copy_in(urls, base)
        ev.acceptProposedAction()

    # ---- search ----------------------------------------------------------------------------------------------------------------------------------
    def filter(self, text):
        t = (text or "").strip().lower()
        budget = [MAX_FILTER_NODES]

        def walk(it, depth):
            d = it.data(0, ROLE)
            if d and d[0] == "dir" and t and depth < 8 and budget[0] > 0:
                self._load(it)
            name_hit = not t or t in it.text(0).lower()
            any_child = False
            for i in range(it.childCount()):
                c = it.child(i)
                budget[0] -= 1
                any_child = walk(c, depth + 1) or any_child
            hit = name_hit or any_child
            it.setHidden(not hit)
            if t and any_child:
                it.setExpanded(True)
            return hit
        for i in range(self.tree.topLevelItemCount()):
            top = self.tree.topLevelItem(i)
            d = top.data(0, ROLE)
            if d and d[0] == "session":
                walk(top, 0)
            else:
                walk(top, 0)
