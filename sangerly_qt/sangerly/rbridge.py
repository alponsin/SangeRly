"""Talks to the R engine (engine/engine.R): one long-running Rscript process, one JSON line per request/answer.
Also finds Rscript / the SangeRly R code (Settings > env variables > portable R > PATH > usual install folders) and can start the classic Shiny app."""
import json
import os
import shutil
import socket
import sys
import tempfile
from pathlib import Path

from PySide6.QtCore import QObject, QProcess, QProcessEnvironment, QSettings, QUrl, Signal
from PySide6.QtGui import QDesktopServices


def app_base() -> Path:
    """Folder that holds the program: next to SangeRly.exe when frozen, the repository folder otherwise."""
    if getattr(sys, "frozen", False):
        return Path(sys.executable).resolve().parent
    return Path(__file__).resolve().parent.parent


def _setting(key: str) -> str:
    v = QSettings().value(key, "")
    return str(v) if v else ""


def common_rscripts():
    """Rscript locations of a normal R installation (Windows, macOS, Linux), newest first."""
    out = []
    if os.name == "nt":
        roots = {os.environ.get(k) for k in ("ProgramFiles", "ProgramFiles(x86)", "LOCALAPPDATA")} - {None}
        for r in roots:
            for sub in ("R", "Programs/R"):
                for d in sorted((Path(r) / sub).glob("R-*"), reverse=True):
                    for rel in ("bin/Rscript.exe", "bin/x64/Rscript.exe"):
                        if (d / rel).exists():
                            out.append(str(d / rel))
        try:
            import winreg
            for hive in (winreg.HKEY_LOCAL_MACHINE, winreg.HKEY_CURRENT_USER):
                try:
                    with winreg.OpenKey(hive, r"SOFTWARE\R-core\R") as k:
                        d = Path(winreg.QueryValueEx(k, "InstallPath")[0])
                        for rel in ("bin/Rscript.exe", "bin/x64/Rscript.exe"):
                            if (d / rel).exists():
                                out.append(str(d / rel))
                except OSError:
                    pass
        except ImportError:
            pass
    else:
        out += ["/usr/bin/Rscript", "/usr/local/bin/Rscript", "/opt/homebrew/bin/Rscript", "/Library/Frameworks/R.framework/Resources/bin/Rscript"]
    return [p for p in out if Path(p).exists()]


def detect_rscript() -> str:
    """Automatic search, without the saved Settings."""
    env = os.environ.get("SANGERLY_RSCRIPT")
    if env and Path(env).exists():
        return env
    exe = "Rscript.exe" if os.name == "nt" else "Rscript"
    cand = app_base() / "R" / "bin" / exe          # portable R shipped with the installer
    if cand.exists():
        return str(cand)
    w = shutil.which("Rscript")
    if w:
        return w
    c = common_rscripts()
    return c[0] if c else ""


def find_rscript() -> str:
    s = _setting("engine/rscript")
    if s and Path(s).exists():
        return s
    return detect_rscript()


def is_project_root(p) -> bool:
    return bool(p) and (Path(p) / "R" / "assembly_core.R").exists()


def has_embed_support(p) -> bool:
    """True when the R code carries the Qt patch (www/sangerly_embed.js + app.R loads it).
    Without it the Workspace keeps the plain Shiny look and the New assembly dialog cannot reach the Assembly tab."""
    try:
        root = Path(p)
        return ((root / "www" / "sangerly_embed.js").exists() and (root / "www" / "sangerly_embed.css").exists()
                and "sangerly_embed" in (root / "app.R").read_text(encoding="utf-8", errors="replace")
                and "qt_cmd" in (root / "R" / "assembly_module.R").read_text(encoding="utf-8", errors="replace"))
    except OSError:
        return False


def _root_candidates():
    env = os.environ.get("SANGERLY_ROOT")
    base = app_base()
    return ([Path(env)] if env else []) + [base / "app", base / "sangerly_r", base.parent / "sangerly_r", base.parent, base]


def detect_project_root() -> str:
    """First valid R-code folder; a folder that already has the Qt patch wins over one that does not."""
    valid = [c for c in _root_candidates() if is_project_root(c)]
    for c in valid:
        if has_embed_support(c):
            return str(c)
    return str(valid[0]) if valid else ""


def find_project_root() -> str:
    """Folder that contains R/assembly_core.R (the R code of SangeRly)."""
    s = _setting("engine/root")
    if is_project_root(s):
        if has_embed_support(s):
            return s
        d = detect_project_root()                 # the saved folder is not patched: use a patched one when there is one
        return d if d and has_embed_support(d) else s
    return detect_project_root()


def _free_port() -> int:
    with socket.socket() as s:
        s.bind(("127.0.0.1", 0))
        return s.getsockname()[1]


class RBridge(QObject):
    log = Signal(str)                   # text from R (stderr) and bridge messages
    starting = Signal()                 # engine process is being started
    ready = Signal(dict)                # engine started: {"r": "R version ...", "pid": ...}
    failed = Signal(str)                # engine could not be started (reason)
    stopped = Signal(int)               # engine ended (exit code)
    busy = Signal(bool)                 # at least one request is waiting for an answer

    def __init__(self, parent=None):
        super().__init__(parent)
        self.proc = QProcess(self)
        self.proc.setProcessChannelMode(QProcess.SeparateChannels)
        self.proc.readyReadStandardOutput.connect(self._on_out)
        self.proc.readyReadStandardError.connect(self._on_err)
        self.proc.finished.connect(lambda code, _s: self.stopped.emit(code))
        self.proc.errorOccurred.connect(lambda e: self.log.emit(f"[R] process error: {e}"))
        self.classic = None             # QProcess of the classic Shiny app, when started
        self._classic_url = ""
        self._buf = b""
        self._next = 1
        self._pending = {}              # id -> (on_done, on_progress)
        self.is_ready = False
        self.rscript = ""
        self.root = ""

    # ---- environment ------------------------------------------------------------------------------------------------
    def _env(self):
        """With the portable R the library folder is fixed; with a normal R installation the user's own library is kept."""
        env = QProcessEnvironment.systemEnvironment()
        lib = app_base() / "R" / "library"
        if lib.exists():
            env.insert("R_LIBS", str(lib))
            env.insert("R_LIBS_USER", "")
        return env, lib.exists()

    # ---- start / stop -----------------------------------------------------------------------------------------------
    def start(self, rscript: str = "", root: str = "") -> bool:
        self.rscript = rscript or find_rscript()
        self.root = root or find_project_root()
        if not self.rscript:
            msg = "Rscript was not found. Install R, or choose Rscript in Settings."
            self.log.emit("[R] " + msg)
            self.failed.emit(msg)
            return False
        if not self.root:
            msg = "The SangeRly R code was not found (a folder with R/assembly_core.R). Choose it in Settings."
            self.log.emit("[R] " + msg)
            self.failed.emit(msg)
            return False
        engine = app_base() / "engine" / "engine.R"
        if not engine.exists():
            msg = f"engine.R not found: {engine}"
            self.log.emit(f"[R] {msg}")
            self.failed.emit(msg)
            return False
        env, portable = self._env()
        self.proc.setProcessEnvironment(env)
        self.log.emit(f"[R] starting: {self.rscript}")
        self.log.emit(f"[R] R code: {self.root}")
        # --vanilla would hide the user's own R library, so it is only used with the portable R that ships its packages
        # --no-init-file: ignore any .Rprofile (e.g. one that calls renv/activate.R); the working folder is a neutral temp folder
        args = (["--vanilla"] if portable else ["--no-save", "--no-restore", "--no-init-file"]) + [str(engine), self.root]
        self.proc.setWorkingDirectory(tempfile.gettempdir())
        self.starting.emit()
        self.proc.start(self.rscript, args)
        return True

    def stop(self):
        if self.proc.state() != QProcess.NotRunning:
            self.proc.write(b'{"action":"quit"}\n')
            if not self.proc.waitForFinished(1500):
                self.proc.kill()
                self.proc.waitForFinished(1000)
        self.is_ready = False
        if self.classic is not None and self.classic.state() != QProcess.NotRunning:
            self.classic.kill()

    def restart(self) -> bool:
        """Stop the engine and start it again with the current Settings."""
        self.stop()
        self._buf = b""
        self._pending.clear()
        self.busy.emit(False)
        return self.start()

    # ---- classic Shiny app ---------------------------------------------------------------------------------------------
    def launch_classic(self) -> bool:
        """Start the classic SangeRly (Shiny) in the browser; it uses the same R and the same project folder."""
        if self.classic is not None and self.classic.state() != QProcess.NotRunning and self._classic_url:
            QDesktopServices.openUrl(QUrl(self._classic_url))
            return True
        rs = self.rscript or find_rscript()
        root = self.root or find_project_root()
        if not rs or not root or not (Path(root) / "app.R").exists():
            self.log.emit("[classic] Rscript or the SangeRly folder (with app.R) was not found. Check Settings.")
            return False
        port = _free_port()
        self._classic_url = f"http://127.0.0.1:{port}"
        self.classic = QProcess(self)
        env, _portable = self._env()
        self.classic.setProcessEnvironment(env)
        self.classic.setProcessChannelMode(QProcess.MergedChannels)
        self.classic.readyReadStandardOutput.connect(lambda: self.log.emit("[classic] " + bytes(self.classic.readAllStandardOutput()).decode("utf-8", "replace").rstrip()))
        code = f"shiny::runApp({json.dumps(Path(root).as_posix())}, port={port}, host='127.0.0.1', launch.browser=TRUE)"
        self.log.emit(f"[classic] starting the classic SangeRly on {self._classic_url}")
        self.classic.setWorkingDirectory(tempfile.gettempdir())
        self.classic.start(rs, ["--no-save", "--no-restore", "--no-init-file", "-e", code])
        return True

    # ---- requests ---------------------------------------------------------------------------------------------------
    def call(self, action: str, args: dict | None = None, on_done=None, on_progress=None) -> int:
        """Send a request; on_done(ok: bool, result_or_error) is called when R answers."""
        rid = self._next
        self._next += 1
        self._pending[rid] = (on_done, on_progress)
        line = json.dumps({"id": rid, "action": action, "args": args or {}}, ensure_ascii=False) + "\n"
        self.proc.write(line.encode("utf-8"))
        self.busy.emit(True)
        return rid

    # ---- reading ----------------------------------------------------------------------------------------------------
    def _on_err(self):
        txt = bytes(self.proc.readAllStandardError()).decode("utf-8", "replace").rstrip()
        if txt:
            self.log.emit("[R] " + txt)

    def _on_out(self):
        self._buf += bytes(self.proc.readAllStandardOutput())
        while b"\n" in self._buf:
            raw, self._buf = self._buf.split(b"\n", 1)
            raw = raw.strip()
            if not raw:
                continue
            try:
                msg = json.loads(raw.decode("utf-8", "replace"))
            except json.JSONDecodeError:
                self.log.emit("[R] " + raw.decode("utf-8", "replace"))
                continue
            ev = msg.get("event")
            if ev == "ready":
                self.is_ready = True
                self.ready.emit(msg)
            elif ev == "progress":
                cb = self._pending.get(msg.get("id"), (None, None))[1]
                if cb:
                    cb(msg)
            elif "id" in msg:
                on_done = self._pending.pop(msg["id"], (None, None))[0]
                ok = bool(msg.get("ok"))
                if on_done:
                    on_done(ok, msg.get("result") if ok else msg.get("error", "Unknown error"))
                if not self._pending:
                    self.busy.emit(False)
