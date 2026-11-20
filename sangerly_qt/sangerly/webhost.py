"""The Shiny workspace inside the Qt window.

ShinyHost  starts the Shiny app of SangeRly (engine/shiny_host.R) on a free local port and tells when it answers.
Workspace  a QWebEngineView that shows that app with ?embed=1 (compact look, no web menu bar). The ribbon of the Qt window opens the tabs
           of the app (sgQt.go), the "New assembly" dialog sends its files and settings to the Assembly tab (sgQt.assemble),
           downloads of the app are saved with a normal Save dialog, and links to other sites open in the system browser.
"""
import json
import threading
import time
import urllib.request
from pathlib import Path

from PySide6.QtCore import QObject, QProcess, QSettings, Qt, QTimer, QUrl, Signal
from PySide6.QtGui import QDesktopServices
from PySide6.QtWidgets import QFileDialog, QLabel, QProgressBar, QStackedLayout, QVBoxLayout, QWidget
from PySide6.QtWebEngineCore import QWebEnginePage, QWebEngineProfile, QWebEngineSettings
from PySide6.QtWebEngineWidgets import QWebEngineView

from .rbridge import _free_port, app_base, find_project_root, find_rscript, has_embed_support


class ShinyHost(QObject):
    log = Signal(str)
    starting = Signal()
    ready = Signal(str)            # base url, e.g. http://127.0.0.1:41234
    failed = Signal(str)
    _answered = Signal(bool, str)  # from the polling thread

    def __init__(self, bridge, parent=None):
        super().__init__(parent)
        self.bridge = bridge
        self.proc = QProcess(self)
        self.proc.setProcessChannelMode(QProcess.MergedChannels)
        self.proc.readyReadStandardOutput.connect(self._on_out)
        self.proc.finished.connect(self._on_finished)
        self.url = ""
        self.is_ready = False
        self._stop = False
        self._answered.connect(self._on_answered)

    def start(self, rscript="", root=""):
        rscript = rscript or self.bridge.rscript or find_rscript()
        root = root or self.bridge.root or find_project_root()
        host = app_base() / "engine" / "shiny_host.R"
        if not rscript or not root or not (Path(root) / "app.R").exists() or not host.exists():
            why = ("Rscript was not found." if not rscript else
                   "The SangeRly folder (with app.R) was not found. Choose it in Settings." if not root or not (Path(root) / "app.R").exists() else
                   f"{host} not found.")
            self.log.emit("[shiny] " + why)
            self.failed.emit(why)
            return False
        env, portable = self.bridge._env()
        self.proc.setProcessEnvironment(env)
        port = _free_port()
        self.url = f"http://127.0.0.1:{port}"
        self.is_ready = False
        self._stop = False
        args = (["--vanilla"] if portable else ["--no-save", "--no-restore", "--no-init-file"]) + [str(host), root, str(port)]
        self.proc.setWorkingDirectory(root)
        self.log.emit(f"[shiny] starting the tools on {self.url}")
        self.starting.emit()
        self.proc.start(rscript, args)
        threading.Thread(target=self._poll, args=(self.url,), daemon=True).start()
        return True

    def _poll(self, url):
        """The first request builds the page (theme, UI), so it can take a while; ask again until it answers 200."""
        t0 = time.time()
        while not self._stop and time.time() - t0 < 600:
            try:
                with urllib.request.urlopen(url + "/", timeout=90) as r:
                    if r.status == 200:
                        self._answered.emit(True, url)
                        return
            except Exception:
                time.sleep(0.6)
        if not self._stop:
            self._answered.emit(False, url)

    def _on_answered(self, ok, url):
        if url != self.url:
            return
        if ok:
            self.is_ready = True
            self.log.emit("[shiny] tools ready")
            self.ready.emit(url)
        else:
            self.failed.emit("The Shiny tools did not start within 10 minutes. See the Log and %LOCALAPPDATA%\\Sangerly\\session.log.")

    def _on_out(self):
        txt = bytes(self.proc.readAllStandardOutput()).decode("utf-8", "replace").rstrip()
        if txt:
            self.log.emit("[shiny] " + txt)

    def _on_finished(self, code, _status):
        if not self._stop:
            self.is_ready = False
            self.log.emit(f"[shiny] stopped (code {code})")
            self.failed.emit(f"The Shiny tools stopped (code {code}). See the Log.")

    def stop(self):
        self._stop = True
        self.is_ready = False
        if self.proc.state() != QProcess.NotRunning:
            self.proc.terminate()
            if not self.proc.waitForFinished(1500):
                self.proc.kill()
                self.proc.waitForFinished(1000)

    def restart(self):
        self.stop()
        return self.start()


class _Page(QWebEnginePage):
    """Only the local app is shown here; other addresses (and links that open a new window) go to the system browser."""

    def __init__(self, profile, parent, log, command=None):
        super().__init__(profile, parent)
        self._log = log
        self._command = command

    def acceptNavigationRequest(self, url, nav_type, is_main_frame):
        if url.scheme() in ("http", "https") and url.host() not in ("127.0.0.1", "localhost"):
            QDesktopServices.openUrl(url)
            return False
        return super().acceptNavigationRequest(url, nav_type, is_main_frame)

    def createWindow(self, _type):
        p = QWebEnginePage(self.profile(), self)
        p.urlChanged.connect(lambda u: (QDesktopServices.openUrl(u), p.deleteLater()))
        return p

    def javaScriptConsoleMessage(self, level, msg, line, src):
        if msg.startswith("SGQT:"):                 # the page asks the window to do something (sangerly_embed.js: console.info("SGQT:...")) 
            if self._command:
                self._command(msg[5:].strip())
            return
        if level == QWebEnginePage.ErrorMessageLevel:
            self._log(f"[web] {msg}  ({Path(src).name}:{line})")


class Workspace(QWidget):
    status = Signal(dict)          # {"ready", "done": [...], "tab", "saved"} every second while the page is open
    loaded = Signal()
    connected = Signal()           # the Shiny session of the page is connected (sgQt can be used)
    command = Signal(str)          # request from the page, e.g. "assembly" (open the assembly window)
    problem = Signal(str)          # the R code cannot be controlled from this window (see NO_PATCH)

    def __init__(self, host: ShinyHost, log, parent=None):
        super().__init__(parent)
        self.host = host
        self.log = log
        self.cfg = QSettings()
        self.want = "assembly"
        self.want_sub = None
        self._queue = []                       # JavaScript waiting for the page to be connected
        self._connected = False
        self._saved_seen = None
        self._loaded_url = ""
        self._loaded_at = 0.0
        self._no_bridge = False

        self.profile = QWebEngineProfile("SangeRlyWeb", self)
        self.profile.setPersistentCookiesPolicy(QWebEngineProfile.NoPersistentCookies)
        self.profile.downloadRequested.connect(self._download)
        self.view = QWebEngineView(self)
        self.page = _Page(self.profile, self.view, log, self.command.emit)
        self.view.setPage(self.page)
        s = self.page.settings()
        s.setAttribute(QWebEngineSettings.LocalContentCanAccessRemoteUrls, False)
        s.setAttribute(QWebEngineSettings.ScrollAnimatorEnabled, True)
        s.setAttribute(QWebEngineSettings.FullScreenSupportEnabled, False)
        self.view.setAcceptDrops(False)

        # overlay shown while the tools start
        self.overlay = QWidget()
        self.overlay.setObjectName("loadingPage")
        ov = QVBoxLayout(self.overlay)
        ov.addStretch(2)
        self.ov_title = QLabel("Starting the SangeRly tools…")
        self.ov_title.setObjectName("h1")
        self.ov_title.setAlignment(Qt.AlignCenter)
        self.ov_sub = QLabel("The first start can take a minute while R loads its packages.")
        self.ov_sub.setObjectName("hint")
        self.ov_sub.setAlignment(Qt.AlignCenter)
        bar = QProgressBar()
        bar.setRange(0, 0)
        bar.setTextVisible(False)
        bar.setMaximumWidth(360)
        bar.setFixedHeight(6)
        self.ov_bar = bar
        ov.addWidget(self.ov_title)
        ov.addWidget(self.ov_sub)
        ov.addWidget(bar, 0, Qt.AlignCenter)
        ov.addStretch(3)

        self.stack = QStackedLayout(self)
        self.stack.setContentsMargins(0, 0, 0, 0)
        self.stack.addWidget(self.overlay)
        self.stack.addWidget(self.view)
        self.stack.setCurrentWidget(self.overlay)

        self.view.loadFinished.connect(self._load_finished)
        self.timer = QTimer(self)
        self.timer.setInterval(1000)
        self.timer.timeout.connect(self._poll_status)
        self.timer.start()
        host.ready.connect(self._host_ready)
        host.failed.connect(self._host_failed)
        host.starting.connect(self._host_starting)

    # ---- host state ----------------------------------------------------------------------------------------------------
    def _host_starting(self):
        self._connected = False
        self.stack.setCurrentWidget(self.overlay)
        self.ov_title.setText("Starting the SangeRly tools…")
        self.ov_sub.setText("The first start can take a minute while R loads its packages.")
        self.ov_bar.setVisible(True)

    def _host_failed(self, why):
        self._connected = False
        self.stack.setCurrentWidget(self.overlay)
        self.ov_title.setText("The tools could not start")
        self.ov_sub.setText(why + "\nCheck Settings (Rscript and the SangeRly folder), then restart.")
        self.ov_bar.setVisible(False)

    def _host_ready(self, url):
        self._loaded_url = url + "/?embed=1"
        self.view.load(QUrl(self._loaded_url))

    def _load_finished(self, ok):
        if not ok:
            self.log("[web] the page did not load; retrying")
            QTimer.singleShot(1500, lambda: self.view.load(QUrl(self._loaded_url)) if self._loaded_url else None)
            return
        self.stack.setCurrentWidget(self.view)
        self._loaded_at = time.time()
        self._no_bridge = False
        self.loaded.emit()

    # ---- talking to the page -----------------------------------------------------------------------------------------------
    def js(self, code, cb=None):
        if cb:
            self.page.runJavaScript(code, cb)
        else:
            self.page.runJavaScript(code)

    def goto(self, tab, sub=None):
        """Open a tool; `sub` is the title of a tab inside it (or a list of titles for nested tabs)."""
        self.want = tab
        self.want_sub = sub
        if self._connected:
            self.js(f"window.sgQt && sgQt.go({json.dumps(tab)}, {json.dumps(sub)});")

    def assemble(self, req):
        """Send the result of the New assembly dialog to the Assembly tab; it runs there with the same code as ASSEMBLE."""
        code = f"window.sgQt ? sgQt.assemble({json.dumps(req, ensure_ascii=False)}) : false;"
        if self._no_bridge:
            self.problem.emit(self._patch_msg())
            return
        if self._connected:
            self.js(code)
        else:
            self._queue.append(code)
            self.log("[task] the tools are still starting; the assembly begins as soon as they are ready")

    def call(self, fn, obj=None):
        """Call sgQt.<fn>(obj) in the page (the dialogs use this); waits until the page is connected."""
        code = f"window.sgQt ? sgQt.{fn}({json.dumps(obj, ensure_ascii=False)}) : false;"
        if self._no_bridge:
            self.problem.emit(self._patch_msg())
        elif self._connected:
            self.js(code)
        else:
            self._queue.append(code)

    def _patch_msg(self):
        root = self.host.bridge.root or find_project_root()
        return (f"The R code in\n{root}\ndoes not match this window (it is older), so the dialogs and tabs cannot work with it "
                "(and the page keeps the old look).\n\nCopy everything inside the folder sangerly_r of this package (app.R, R, www) into that "
                "folder, overwrite the files, and restart SangeRly.")

    def _poll_status(self):
        if self.stack.currentWidget() is not self.view:
            return

        def got(res):
            if not res and self._loaded_at and time.time() - self._loaded_at > 8 and not self._no_bridge:
                # the page is open, but sangerly_embed.js did not run: the R code is not patched
                self._no_bridge = True
                self.log("[error] " + self._patch_msg().replace("\n\n", "  ").replace("\n", " "))
                self.problem.emit(self._patch_msg())
            try:
                st = json.loads(res) if res else {}
            except (TypeError, ValueError):
                return
            now = bool(st.get("ready"))
            if now and not self._connected:
                self._connected = True
                self.connected.emit()
                self.goto(self.want, self.want_sub)
                for c in self._queue:
                    self.js(c)
                self._queue.clear()
            elif not now:
                self._connected = False
            sv = st.get("saved")
            if sv and sv != self._saved_seen:
                self._saved_seen = sv
                self.log(f"[task] {sv.get('n', 0)} file(s) saved in {sv.get('dir', '')}")
            self.status.emit(st)

        self.page.runJavaScript("window.sgQt ? sgQt.status() : ''", got)

    # ---- downloads of the app (FASTA, CSV, figures, reports) ------------------------------------------------------------
    def _download(self, d):
        name = d.downloadFileName() or "download"
        start = self.cfg.value("download/dir", str(Path.home() / "Documents" / "SangeRly"))
        Path(start).mkdir(parents=True, exist_ok=True)
        path, _ = QFileDialog.getSaveFileName(self.window(), "Save file", str(Path(start) / name))
        if not path:
            d.cancel()
            return
        self.cfg.setValue("download/dir", str(Path(path).parent))
        d.setDownloadDirectory(str(Path(path).parent))
        d.setDownloadFileName(Path(path).name)
        d.isFinishedChanged.connect(lambda: self.log(f"[task] saved {path}") if d.isFinished() else None)
        d.accept()
