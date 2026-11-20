"""Reusable widgets: hero banner with a chromatogram, toolbar button, module card, chips, stat tiles, colour-coded sequence viewer, drop frame."""
import html
import math
import random

from PySide6.QtCore import QRectF, QSize, Qt, Signal
from PySide6.QtGui import QBrush, QColor, QFont, QLinearGradient, QPainter, QPainterPath, QPen
from PySide6.QtWidgets import QFrame, QHBoxLayout, QLabel, QTextEdit, QToolButton, QVBoxLayout, QWidget

from . import theme

CHIP_COLORS = {
    "ok": ("#dcfce7", "#166534"),
    "warn": ("#fef3c7", "#92400e"),
    "bad": ("#fee2e2", "#991b1b"),
    "info": (theme.ACCENT_SOFT, theme.ACCENT_DARK),
    "muted": ("#f1f5f9", theme.MUTED),
}


class Chip(QLabel):
    """Small coloured pill: Chip("Ready", "ok")."""

    def __init__(self, text, kind="muted"):
        super().__init__(text)
        bg, fg = CHIP_COLORS.get(kind, CHIP_COLORS["muted"])
        self.setStyleSheet(f"background:{bg}; color:{fg}; border-radius:9px; padding:2px 10px; font-size:8pt; font-weight:600;")


def base_legend():
    """One-line colour key of the sequence viewer (A / C / G / T / N-IUPAC), built from the same colours as the viewer."""
    w = QWidget()
    h = QHBoxLayout(w)
    h.setContentsMargins(0, 0, 0, 0)
    h.setSpacing(6)
    for key, text in [("A", "A"), ("C", "C"), ("G", "G"), ("T", "T"), ("N", "N / IUPAC")]:
        bg, fg = theme.BASE_COLORS.get(key, theme.BASE_OTHER)
        l = QLabel(text)
        l.setStyleSheet(f"background:{bg}; color:{fg}; border-radius:4px; padding:1px 8px; font-size:8pt; font-weight:700;")
        h.addWidget(l)
    h.addStretch(1)
    return w


def section_title(text, hint=""):
    w = QWidget()
    v = QVBoxLayout(w)
    v.setContentsMargins(0, 0, 0, 0)
    v.setSpacing(2)
    t = QLabel(text)
    t.setObjectName("h2")
    v.addWidget(t)
    if hint:
        h = QLabel(hint)
        h.setObjectName("hint")
        h.setWordWrap(True)
        v.addWidget(h)
    return w


def panel_head(title, right=None):
    """Header strip of a panel (grey bar with a bold title, optional widget on the right)."""
    f = QFrame()
    f.setObjectName("panelHead")
    h = QHBoxLayout(f)
    h.setContentsMargins(14, 8, 12, 8)
    h.setSpacing(8)
    t = QLabel(title)
    t.setObjectName("panelTitle")
    h.addWidget(t)
    h.addStretch(1)
    if right is not None:
        h.addWidget(right)
    return f


class StatTile(QFrame):
    def __init__(self, value, label, color=theme.NAVY):
        super().__init__()
        self.setObjectName("stat")
        v = QVBoxLayout(self)
        v.setContentsMargins(14, 10, 14, 10)
        v.setSpacing(0)
        a = QLabel(str(value))
        a.setObjectName("statValue")
        a.setStyleSheet(f"color:{color};")
        b = QLabel(f"<span style='color:{color}'>●</span>&nbsp; {html.escape(label.upper())}")
        b.setObjectName("statLabel")
        v.addWidget(a)
        v.addWidget(b)


class ModuleCard(QFrame):
    """Clickable card of the workflow grid (soft shadow, accent border on hover)."""
    clicked = Signal()

    def __init__(self, ico, title, desc, chip_text, chip_kind, step=None):
        super().__init__()
        self.setObjectName("card")
        self.setCursor(Qt.PointingHandCursor)
        self.setMinimumSize(210, 138)
        self._fx = theme.shadow(self, 14, 2, 22)
        h = QHBoxLayout(self)
        h.setContentsMargins(16, 14, 16, 14)
        h.setSpacing(14)
        i = QLabel()
        i.setPixmap(ico.pixmap(QSize(46, 46)))
        i.setFixedSize(46, 46)
        h.addWidget(i, 0, Qt.AlignTop)
        v = QVBoxLayout()
        v.setSpacing(4)
        head = QHBoxLayout()
        head.setSpacing(6)
        t = QLabel(title)
        t.setObjectName("cardTitle")
        head.addWidget(t, 1)
        if step:
            n = QLabel(str(step))
            n.setFixedSize(18, 18)
            n.setAlignment(Qt.AlignCenter)
            n.setStyleSheet(f"background:{theme.HOVER}; color:{theme.MUTED}; border-radius:9px; font-size:8pt; font-weight:700;")
            head.addWidget(n)
        v.addLayout(head)
        d = QLabel(desc)
        d.setObjectName("hint")
        d.setWordWrap(True)
        v.addWidget(d, 1)
        v.addWidget(Chip(chip_text, chip_kind), 0, Qt.AlignLeft)
        h.addLayout(v, 1)
        for c in self.findChildren(QWidget):          # the card itself receives hover / click, not the labels on it
            c.setAttribute(Qt.WA_TransparentForMouseEvents, True)

    def enterEvent(self, e):
        self._fx.setBlurRadius(24)
        self._fx.setOffset(0, 5)
        super().enterEvent(e)

    def leaveEvent(self, e):
        self._fx.setBlurRadius(14)
        self._fx.setOffset(0, 2)
        super().leaveEvent(e)

    def mouseReleaseEvent(self, e):
        if e.button() == Qt.LeftButton and self.rect().contains(e.position().toPoint()):
            self.clicked.emit()
        super().mouseReleaseEvent(e)


class RibbonButton(QToolButton):
    """Icon button of the toolbar, with a step number (top-left) and a status dot (top-right, green = the step has results)."""

    def __init__(self, ico, text, tip, step=None, ready=None):
        super().__init__()
        self.setObjectName("ribbonBtn")
        self._step, self._ready = step, ready
        self.setIcon(ico)
        self.setIconSize(QSize(30, 30))
        self.setText(text)
        self.setToolTip(tip)
        self.setToolButtonStyle(Qt.ToolButtonTextUnderIcon)
        self.setCursor(Qt.PointingHandCursor)
        self.setFixedSize(84, 62)

    def set_done(self, done: bool):
        """Green dot = the step already has results in the workspace."""
        if (self._ready is True) != bool(done):
            self._ready = True if done else None
            self.update()

    def paintEvent(self, e):
        super().paintEvent(e)
        p = QPainter(self)
        p.setRenderHint(QPainter.Antialiasing)
        if self._step:
            r = QRectF(7, 5, 14, 14)
            p.setPen(Qt.NoPen)
            p.setBrush(QColor(theme.HOVER))
            p.drawEllipse(r)
            f = QFont(self.font())
            f.setPointSizeF(7)
            f.setBold(True)
            p.setFont(f)
            p.setPen(QColor(theme.MUTED))
            p.drawText(r, Qt.AlignCenter, str(self._step))
        if self._ready is not None:
            p.setPen(Qt.NoPen)
            p.setBrush(QColor(theme.OK if self._ready else theme.WARN))
            p.drawEllipse(QRectF(self.width() - 15, 7, 8, 8))
        p.end()


class ChromatogramBanner(QFrame):
    """Hero background: blue gradient with a synthetic four-colour chromatogram on the right."""
    TRACE = {"A": "#4ade80", "C": "#60a5fa", "G": "#fbbf24", "T": "#f87171"}

    def __init__(self):
        super().__init__()
        self.setMinimumHeight(190)
        rnd = random.Random(11)
        n = 16
        self._calls = [(rnd.choice("ACGT"), (i + 0.5) / n + rnd.uniform(-0.01, 0.01), rnd.uniform(0.5, 1.0)) for i in range(n)]
        self._noise = [(rnd.choice("ACGT"), rnd.uniform(0.05, 0.95), rnd.uniform(0.06, 0.16)) for _ in range(14)]

    def _level(self, base, t):
        y = 0.0
        for b, c, h in self._calls:
            if b == base:
                y += h * math.exp(-((t - c) / 0.014) ** 2)
        for b, c, h in self._noise:
            if b == base:
                y += h * math.exp(-((t - c) / 0.02) ** 2)
        return y

    def paintEvent(self, ev):
        r = self.rect()
        p = QPainter(self)
        p.setRenderHint(QPainter.Antialiasing)
        clip = QPainterPath()
        clip.addRoundedRect(QRectF(r), 14, 14)
        p.setClipPath(clip)
        g = QLinearGradient(0, 0, r.width(), r.height())
        g.setColorAt(0, QColor("#0f1f3d"))
        g.setColorAt(0.55, QColor("#16399a"))
        g.setColorAt(1, QColor("#2563eb"))
        p.fillRect(r, QBrush(g))
        p.setPen(QPen(QColor(255, 255, 255, 12), 1))
        for x in range(0, r.width(), 36):
            p.drawLine(x, 0, x, r.height())
        for y in range(0, r.height(), 36):
            p.drawLine(0, y, r.width(), y)
        w, h = r.width(), r.height()
        x0 = w * 0.44
        span = max(1.0, w - x0 - 16)
        base_y, amp = h * 0.88, h * 0.60
        for base, color in self.TRACE.items():
            pen = QPen(QColor(color), 2.2)
            pen.setCapStyle(Qt.RoundCap)
            p.setPen(pen)
            p.setBrush(Qt.NoBrush)
            path = QPainterPath()
            x = x0
            first = True
            while x <= w - 8:
                t = (x - x0) / span
                env = min(1.0, t * 7.0, (1.0 - t) * 14.0 + 0.2)
                yy = base_y - amp * self._level(base, t) * env
                if first:
                    path.moveTo(x, yy)
                    first = False
                else:
                    path.lineTo(x, yy)
                x += 3
            p.drawPath(path)
        f = QFont()
        f.setPointSizeF(9)
        f.setBold(True)
        p.setFont(f)
        for b, c, hgt in self._calls:
            p.setPen(QColor(self.TRACE[b]))
            p.drawText(QRectF(x0 + c * span - 8, 12, 16, 16), Qt.AlignCenter, b)
        p.end()


class DockTitle(QWidget):
    def __init__(self, text, dock):
        super().__init__()
        h = QHBoxLayout(self)
        h.setContentsMargins(12, 7, 6, 7)
        l = QLabel(text.upper())
        l.setObjectName("dockTitle")
        b = QToolButton()
        b.setObjectName("dockClose")
        b.setText("×")
        b.setCursor(Qt.PointingHandCursor)
        b.clicked.connect(dock.close)
        h.addWidget(l, 1)
        h.addWidget(b)


class DropFrame(QFrame):
    """A frame that accepts files / folders dropped from the file manager."""
    dropped = Signal(list)

    def __init__(self):
        super().__init__()
        self.setObjectName("dropFrame")
        self.setAcceptDrops(True)

    def _hl(self, on):
        self.setProperty("drag", bool(on))
        theme.repolish(self)

    def dragEnterEvent(self, e):
        if e.mimeData().hasUrls():
            e.acceptProposedAction()
            self._hl(True)

    def dragLeaveEvent(self, e):
        self._hl(False)

    def dropEvent(self, e):
        self._hl(False)
        ps = [u.toLocalFile() for u in e.mimeData().urls() if u.isLocalFile()]
        if ps:
            self.dropped.emit(ps)
            e.acceptProposedAction()


class SequenceView(QTextEdit):
    """Read-only nucleotide viewer: colour-coded bases in blocks of 10, 60 per line, with a column ruler and position numbers
    ."""

    def __init__(self):
        super().__init__()
        self.setObjectName("seq")
        self.setReadOnly(True)
        self.setLineWrapMode(QTextEdit.NoWrap)
        self.setFont(theme.mono_font(10))

    def set_sequence(self, seq, per_line=60, block=10):
        if not seq:
            self.clear()
            return
        spans = {}
        for b in set(seq.upper()):
            bg, fg = theme.BASE_COLORS.get(b, theme.BASE_OTHER)
            spans[b] = f"<span style='background-color:{bg}; color:{fg}; font-weight:600;'>{html.escape(b)}</span>"
        faint = theme.FAINT
        out = ["<pre style='margin:0; line-height:135%;'>"]
        ruler = " ".join(f"{i + block:>{block}}" for i in range(0, per_line, block))
        out.append(f"<span style='color:{faint};'>{'':>8}{ruler}</span>\n")
        for start in range(0, len(seq), per_line):
            chunk = seq[start:start + per_line].upper()
            parts = ["".join(spans[b] for b in chunk[i:i + block]) for i in range(0, len(chunk), block)]
            out.append(f"<span style='color:{faint};'>{start + 1:>6}  </span>" + " ".join(parts) + "\n")
        out.append("</pre>")
        self.setHtml("".join(out))
