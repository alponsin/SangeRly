"""Look of SangeRly (v3): clean light chrome.

White top bar (brand, menus, engine state), grouped toolbar with captions, project tree on the left, document tabs with an
accent underline, soft cards and panels, slate neutrals with one blue accent. All colours are tokens below; the style sheet
uses them as @TOKEN@. Everything small (icons, check marks, arrows) is drawn in code, so no extra image files are needed
except the logo in sangerly/assets/.
"""
import math
import re
import tempfile
from pathlib import Path

from PySide6.QtCore import QPointF, QRectF, Qt
from PySide6.QtGui import (QBrush, QColor, QFont, QIcon, QLinearGradient, QPainter, QPainterPath, QPalette, QPen, QPixmap,
                           QPolygonF)

# ---- palette ----------------------------------------------------------------------------------------------------------
NAVY = "#0f1f3d"          # headings
NAVY2 = "#1d3a73"
ACCENT = "#2563eb"
ACCENT_DARK = "#1d4ed8"
ACCENT_SOFT = "#e8f0ff"
BG = "#f4f6fa"            # work area behind panels
PANEL = "#ffffff"
SUBTLE = "#f8fafc"        # panel headers, table headers, alternate rows
HOVER = "#eef2f7"
BORDER = "#e2e8f0"
BORDER2 = "#cbd5e1"       # inputs, hover borders
TEXT = "#0f172a"
MUTED = "#64748b"
FAINT = "#94a3b8"
OK = "#16a34a"
WARN = "#d97706"
BAD = "#dc2626"

# nucleotide colours (background, text) used by the sequence viewer
BASE_COLORS = {
    "A": ("#bbf7d0", "#14532d"),
    "C": ("#bfdbfe", "#1e3a8a"),
    "G": ("#fde68a", "#78350f"),
    "T": ("#fecaca", "#7f1d1d"),
    "-": ("#ffffff", "#94a3b8"),
}
BASE_OTHER = ("#e2e8f0", "#475569")          # N and IUPAC ambiguity codes

ASSETS = Path(__file__).resolve().parent / "assets"
UI_FAMILIES = ["Segoe UI Variable Text", "Segoe UI", "Inter", "SF Pro Text", "Noto Sans", "DejaVu Sans"]


def repolish(w):
    """Re-apply the style sheet after a dynamic property (state, selected, drag) changed."""
    w.style().unpolish(w)
    w.style().polish(w)
    w.update()


def mono_font(size=10):
    f = QFont()
    f.setFamilies(["Cascadia Mono", "Consolas", "JetBrains Mono", "DejaVu Sans Mono", "Menlo", "Courier New", "monospace"])
    f.setStyleHint(QFont.Monospace)
    f.setPointSize(size)
    return f


# ---- logo -------------------------------------------------------------------------------------------------------------
def logo_pixmap(size=64) -> QPixmap:
    pm = QPixmap(str(ASSETS / "logo.png"))
    if pm.isNull():                       # logo file missing: draw a simple tile instead
        return icon("assembly", size, ACCENT).pixmap(size, size)
    return pm.scaled(size * 2, size * 2, Qt.KeepAspectRatio, Qt.SmoothTransformation)


def app_icon() -> QIcon:
    for name in ("logo.ico", "logo.png"):
        p = ASSETS / name
        if p.exists():
            ic = QIcon(str(p))
            if not ic.isNull():
                return ic
    return icon("assembly", 64, ACCENT)


# ---- small images used by the style sheet (check mark, arrows ...) drawn once, saved as PNG (+@2x) --------------------------
def _make_assets():
    d = Path(tempfile.gettempdir()) / "sangerly_qt_assets"
    d.mkdir(exist_ok=True)

    def png(name, w, h, draw):
        for sc, suffix in ((1, ""), (2, "@2x")):
            pm = QPixmap(w * sc, h * sc)
            pm.fill(Qt.transparent)
            p = QPainter(pm)
            p.setRenderHint(QPainter.Antialiasing)
            p.scale(sc, sc)
            draw(p)
            p.end()
            pm.save(str(d / f"{name}{suffix}.png"))
        return (d / f"{name}.png").as_posix()

    def line(p, color, width, pts):
        p.setPen(QPen(QColor(color), width, Qt.SolidLine, Qt.RoundCap, Qt.RoundJoin))
        p.setBrush(Qt.NoBrush)
        p.drawPolyline(QPolygonF([QPointF(x, y) for x, y in pts]))

    def dot(p):
        p.setPen(QPen(QColor(ACCENT), 1.2))
        p.setBrush(QColor("white"))
        p.drawEllipse(QRectF(1, 1, 14, 14))
        p.setPen(Qt.NoPen)
        p.setBrush(QColor(ACCENT))
        p.drawEllipse(QRectF(4.5, 4.5, 7, 7))

    return {
        "CHECK": png("check", 16, 16, lambda p: line(p, "white", 2, [(3.5, 8.5), (6.5, 11.5), (12.5, 4.5)])),
        "UP": png("up", 10, 10, lambda p: line(p, MUTED, 1.6, [(2, 7), (5, 3.5), (8, 7)])),
        "DOWN": png("down", 10, 10, lambda p: line(p, MUTED, 1.6, [(2, 3.5), (5, 7), (8, 3.5)])),
        "CLOSE": png("close", 12, 12, lambda p: (line(p, MUTED, 1.5, [(3, 3), (9, 9)]), line(p, MUTED, 1.5, [(9, 3), (3, 9)]))),
        "DOT": png("dot", 16, 16, dot),
    }


QSS = """
QMainWindow, QDialog { background: @BG@; }
QWidget#page { background: @BG@; }
QScrollArea { border: none; background: @BG@; }
QToolTip { background: @NAVY@; color: white; border: none; padding: 5px 9px; }
QWidget { color: @TEXT@; }

/* ---- top bar: logo, menus, engine state --------------------------------------------------------------------------- */
QFrame#topBar { background: @PANEL@; border: none; border-bottom: 1px solid @BORDER@; }
QLabel#brandName { color: @NAVY@; font-size: 12pt; font-weight: 700; background: transparent; }
QMenuBar#brandMenu { background: transparent; border: none; }
QMenuBar#brandMenu::item { background: transparent; color: @TEXT@; padding: 5px 11px; border-radius: 6px; margin: 0 1px; }
QMenuBar#brandMenu::item:selected { background: @HOVER@; }
QMenuBar#brandMenu::item:pressed { background: @ACCENT_SOFT@; color: @ACCENT_DARK@; }
QMenu { background: @PANEL@; border: 1px solid @BORDER@; padding: 6px; }
QMenu::item { padding: 7px 30px 7px 12px; border-radius: 5px; color: @TEXT@; }
QMenu::item:selected { background: @ACCENT_SOFT@; color: @ACCENT_DARK@; }
QMenu::item:disabled { color: @FAINT@; }
QMenu::separator { height: 1px; background: @BORDER@; margin: 5px 8px; }
QMenu::icon { padding-left: 6px; }
QPushButton#enginePill { border-radius: 11px; padding: 2px 12px; font-size: 8.5pt; font-weight: 600; border: 1px solid transparent; min-height: 18px; }
QPushButton#enginePill[state="starting"] { background: #fef3c7; color: #92400e; border-color: #fde68a; }
QPushButton#enginePill[state="ready"] { background: #dcfce7; color: #166534; border-color: #bbf7d0; }
QPushButton#enginePill[state="error"] { background: #fee2e2; color: #991b1b; border-color: #fecaca; }
QPushButton#enginePill:hover { border-color: @BORDER2@; }

/* ---- toolbar (grouped module buttons, with captions) --------------------------------------------------------------- */
QFrame#toolbar { background: @PANEL@; border: none; border-bottom: 1px solid @BORDER@; }
QLabel#groupCaption { color: @FAINT@; font-size: 7pt; font-weight: 700; background: transparent; }
QToolButton#ribbonBtn { border: 1px solid transparent; border-radius: 8px; padding: 5px 2px 2px 2px; font-size: 7.5pt; font-weight: 600; color: @TEXT@; background: transparent; }
QToolButton#ribbonBtn:hover { background: @HOVER@; border-color: @BORDER@; }
QToolButton#ribbonBtn:pressed { background: @ACCENT_SOFT@; border-color: #bfd3fb; }
QFrame#ribbonSep { background: @BORDER@; border: none; }

/* ---- docks, tabs, splitters ---------------------------------------------------------------------------------------- */
QLabel#dockTitle { font-size: 7.5pt; font-weight: 700; color: @MUTED@; background: transparent; }
QToolButton#dockClose { border: none; border-radius: 4px; font-size: 11pt; color: @MUTED@; padding: 0 6px; background: transparent; }
QToolButton#dockClose:hover { background: @HOVER@; color: @TEXT@; }
QDockWidget > QWidget { background: @PANEL@; }
QTabWidget::pane { border: none; border-top: 1px solid @BORDER@; background: @BG@; }
QTabBar { background: transparent; }
QTabBar::tab { background: transparent; padding: 9px 14px 8px 14px; margin-right: 2px; color: @MUTED@; border-bottom: 2px solid transparent; }
QTabBar::tab:selected { color: @TEXT@; border-bottom: 2px solid @ACCENT@; font-weight: 600; }
QTabBar::tab:hover:!selected { color: @TEXT@; background: @HOVER@; border-top-left-radius: 6px; border-top-right-radius: 6px; }
QTabBar::close-button { image: url(@CLOSE@); subcontrol-position: right; margin: 2px; border-radius: 4px; }
QTabBar::close-button:hover { background: #fee2e2; }
QTabWidget#inner::pane { border: 1px solid @BORDER@; background: @PANEL@; border-radius: 10px; top: -1px; }
QTabWidget#inner QTabBar::tab { padding: 8px 16px 7px 16px; }
QTabWidget#inner QTabBar::tab:selected { background: transparent; }
QTabWidget#inner QTableView, QTabWidget#inner QPlainTextEdit, QTabWidget#inner QTextEdit#seq { border: none; border-radius: 0; background: @PANEL@; }
QSplitter::handle { background: transparent; }
QSplitter::handle:horizontal { width: 8px; }
QSplitter::handle:vertical { height: 8px; }
QSplitter::handle:hover { background: @ACCENT_SOFT@; }
QMainWindow::separator { background: @BG@; width: 5px; height: 5px; }
QMainWindow::separator:hover { background: @ACCENT_SOFT@; }

/* ---- text ---------------------------------------------------------------------------------------------------------- */
QLabel { background: transparent; }
QLabel#h1 { font-size: 20pt; font-weight: 600; color: @NAVY@; }
QWidget#loadingPage { background: @BG@; }
QLabel#h2 { font-size: 13pt; font-weight: 700; color: @NAVY@; }
QLabel#muted { color: @MUTED@; }
QLabel#hint { color: @MUTED@; font-size: 9pt; }
QLabel#cardTitle { font-size: 10.5pt; font-weight: 700; color: @NAVY@; }
QLabel#panelTitle { font-size: 9.5pt; font-weight: 700; color: @NAVY@; }
QLabel#heroTitle { color: white; font-size: 24pt; font-weight: 600; }
QLabel#heroSub { color: #d6e2ff; font-size: 11pt; }
QLabel#statValue { font-size: 18pt; font-weight: 700; }
QLabel#statLabel { color: @MUTED@; font-size: 8pt; font-weight: 700; }

/* ---- cards, panels ------------------------------------------------------------------------------------------------- */
QFrame#card { background: @PANEL@; border: 1px solid @BORDER@; border-radius: 12px; }
QFrame#card:hover { border-color: @ACCENT@; }
QFrame#panel { background: @PANEL@; border: 1px solid @BORDER@; border-radius: 12px; }
QFrame#panel QTableView, QFrame#panel QListView, QFrame#panel QTreeView { border: none; border-radius: 0; background: transparent; }
QFrame#panelHead { background: @SUBTLE@; border: none; border-bottom: 1px solid @BORDER@; border-top-left-radius: 12px; border-top-right-radius: 12px; }
QFrame#stat { background: @PANEL@; border: 1px solid @BORDER@; border-radius: 10px; }
QFrame#dlgHead { background: @PANEL@; border: none; border-bottom: 1px solid @BORDER@; }
QFrame#dlgFoot { background: @PANEL@; border: none; border-top: 1px solid @BORDER@; }
QFrame#optCard { background: @PANEL@; border: 1px solid @BORDER@; border-radius: 10px; }
QFrame#optCard:hover { border-color: #9bbcf5; }
QFrame#optCard[selected="true"] { border: 2px solid @ACCENT@; background: #f5f9ff; }
QFrame#dropFrame { background: @PANEL@; border: 1px solid @BORDER@; border-radius: 12px; }
QFrame#dropFrame[drag="true"] { border: 2px dashed @ACCENT@; background: @ACCENT_SOFT@; }
QGroupBox { background: @PANEL@; border: 1px solid @BORDER@; border-radius: 10px; margin-top: 14px; padding: 16px 12px 10px 12px; }
QGroupBox::title { subcontrol-origin: margin; left: 14px; padding: 0 6px; color: @NAVY@; font-weight: 700; }

/* ---- buttons ------------------------------------------------------------------------------------------------------- */
QPushButton { background: @PANEL@; border: 1px solid @BORDER2@; border-radius: 7px; padding: 5px 15px; min-height: 18px; }
QPushButton:hover { border-color: @ACCENT@; background: @ACCENT_SOFT@; color: @ACCENT_DARK@; }
QPushButton:pressed { background: #d8e5fd; }
QPushButton:disabled { color: #a3acb8; background: #f1f4f8; border-color: @BORDER@; }
QPushButton#primary { background: @ACCENT@; border-color: @ACCENT@; color: white; font-weight: 700; }
QPushButton#primary:hover { background: @ACCENT_DARK@; border-color: @ACCENT_DARK@; color: white; }
QPushButton#primary:disabled { background: #a9c6f7; border-color: #a9c6f7; color: white; }
QPushButton#heroPrimary { background: white; border: 1px solid white; color: @NAVY@; font-weight: 700; padding: 8px 20px; border-radius: 8px; }
QPushButton#heroPrimary:hover { background: #e8f0ff; color: @NAVY@; }
QPushButton#heroGhost { background: rgba(255, 255, 255, 28); border: 1px solid rgba(255, 255, 255, 110); color: white; padding: 8px 18px; border-radius: 8px; }
QPushButton#heroGhost:hover { background: rgba(255, 255, 255, 60); color: white; }
QPushButton#link { background: transparent; border: none; color: @ACCENT_DARK@; text-align: left; padding: 5px 2px; font-weight: 600; }
QPushButton#link:hover { color: @ACCENT@; background: transparent; text-decoration: underline; }
QMessageBox QPushButton { min-width: 84px; }

/* ---- inputs -------------------------------------------------------------------------------------------------------- */
QLineEdit, QSpinBox, QDoubleSpinBox, QComboBox, QPlainTextEdit, QTextEdit {
  background: @PANEL@; border: 1px solid @BORDER2@; border-radius: 7px; padding: 4px 8px; selection-background-color: @ACCENT@; selection-color: white; }
QLineEdit:hover, QSpinBox:hover, QDoubleSpinBox:hover, QComboBox:hover { border-color: #94a3b8; }
QLineEdit:focus, QSpinBox:focus, QDoubleSpinBox:focus, QComboBox:focus { border: 1px solid @ACCENT@; }
QLineEdit:disabled, QSpinBox:disabled, QDoubleSpinBox:disabled, QComboBox:disabled { background: #f1f4f8; color: #a3acb8; }
QSpinBox, QDoubleSpinBox { padding-right: 4px; min-height: 20px; }
QComboBox { min-height: 20px; padding-right: 22px; }
QComboBox QAbstractItemView { background: @PANEL@; border: 1px solid @BORDER@; selection-background-color: @ACCENT_SOFT@; selection-color: @ACCENT_DARK@; outline: none; }
QComboBox::drop-down { subcontrol-origin: padding; subcontrol-position: center right; width: 22px; border: none; }
QComboBox::down-arrow { image: url(@DOWN@); width: 10px; height: 10px; }
QSpinBox::up-button, QDoubleSpinBox::up-button { subcontrol-origin: border; subcontrol-position: top right; width: 20px; border-left: 1px solid @BORDER@; border-top-right-radius: 6px; }
QSpinBox::down-button, QDoubleSpinBox::down-button { subcontrol-origin: border; subcontrol-position: bottom right; width: 20px; border-left: 1px solid @BORDER@; border-bottom-right-radius: 6px; }
QSpinBox::up-button:hover, QSpinBox::down-button:hover, QDoubleSpinBox::up-button:hover, QDoubleSpinBox::down-button:hover { background: @ACCENT_SOFT@; }
QSpinBox::up-arrow, QDoubleSpinBox::up-arrow { image: url(@UP@); width: 10px; height: 10px; }
QSpinBox::down-arrow, QDoubleSpinBox::down-arrow { image: url(@DOWN@); width: 10px; height: 10px; }
QCheckBox, QRadioButton { spacing: 8px; background: transparent; }
QCheckBox::indicator, QRadioButton::indicator { width: 16px; height: 16px; }
QCheckBox::indicator { border: 1px solid #aeb7c4; border-radius: 4px; background: white; }
QCheckBox::indicator:hover { border-color: @ACCENT@; }
QCheckBox::indicator:checked { background: @ACCENT@; border-color: @ACCENT@; image: url(@CHECK@); }
QRadioButton::indicator { border: 1px solid #aeb7c4; border-radius: 8px; background: white; }
QRadioButton::indicator:checked { border: none; image: url(@DOT@); }
QRadioButton#optRadio { font-weight: 700; color: @NAVY@; }

/* ---- lists, tables, trees ------------------------------------------------------------------------------------------ */
QTableView, QTreeView, QListView { background: @PANEL@; alternate-background-color: @SUBTLE@; border: 1px solid @BORDER@; border-radius: 8px;
  gridline-color: transparent; selection-background-color: @ACCENT_SOFT@; selection-color: @TEXT@; outline: none; }
QTreeView, QListView { padding: 4px; }
QTableView::item { padding: 5px 6px; border: none; }
QTreeView::item, QListView::item { padding: 5px 4px; border-radius: 5px; }
QTreeView::item:hover, QListView::item:hover { background: @HOVER@; }
QTreeView::item:selected, QListView::item:selected, QTableView::item:selected { background: @ACCENT_SOFT@; color: @ACCENT_DARK@; }
QHeaderView::section { background: @SUBTLE@; border: none; border-right: 1px solid @BORDER@; border-bottom: 1px solid @BORDER@; padding: 7px 10px; font-weight: 600; font-size: 8.5pt; color: @MUTED@; }
QTableCornerButton::section { background: @SUBTLE@; border: none; }
QTextEdit#seq, QTextEdit#log { background: @PANEL@; border: 1px solid @BORDER@; border-radius: 8px; padding: 8px; }
QTextEdit#log { border: none; border-radius: 0; padding: 6px 10px; }

/* ---- status bar, progress, scroll bars ----------------------------------------------------------------------------- */
QStatusBar { background: @PANEL@; border-top: 1px solid @BORDER@; color: @MUTED@; }
QStatusBar QLabel { color: @MUTED@; font-size: 8.5pt; }
QStatusBar::item { border: none; }
QProgressBar { border: none; border-radius: 4px; background: #e2e8f0; text-align: center; max-height: 8px; color: transparent; }
QProgressBar::chunk { background: @ACCENT@; border-radius: 4px; }
QScrollBar:vertical { background: transparent; width: 12px; margin: 0; }
QScrollBar::handle:vertical { background: @BORDER2@; border-radius: 5px; min-height: 28px; margin: 2px; }
QScrollBar::handle:vertical:hover { background: @FAINT@; }
QScrollBar:horizontal { background: transparent; height: 12px; }
QScrollBar::handle:horizontal { background: @BORDER2@; border-radius: 5px; min-width: 28px; margin: 2px; }
QScrollBar::handle:horizontal:hover { background: @FAINT@; }
QScrollBar::add-line, QScrollBar::sub-line { width: 0; height: 0; }
QScrollBar::add-page, QScrollBar::sub-page { background: transparent; }
"""

_TOKENS = {"NAVY2": NAVY2, "NAVY": NAVY, "ACCENT_DARK": ACCENT_DARK, "ACCENT_SOFT": ACCENT_SOFT, "ACCENT": ACCENT, "BG": BG,
           "PANEL": PANEL, "SUBTLE": SUBTLE, "HOVER": HOVER, "BORDER2": BORDER2, "BORDER": BORDER, "TEXT": TEXT, "MUTED": MUTED, "FAINT": FAINT}


def build_qss(assets):
    tokens = dict(_TOKENS)
    tokens.update(assets)
    return re.sub(r"@([A-Z0-9_]+)@", lambda m: tokens[m.group(1)], QSS)      # an unknown token raises KeyError right away


def _palette():
    """Explicit light palette, so the window looks the same when the operating system is in dark mode."""
    pal = QPalette()
    for role, col in [(QPalette.Window, BG), (QPalette.WindowText, TEXT), (QPalette.Base, PANEL), (QPalette.AlternateBase, SUBTLE),
                      (QPalette.Text, TEXT), (QPalette.Button, PANEL), (QPalette.ButtonText, TEXT), (QPalette.Highlight, ACCENT),
                      (QPalette.HighlightedText, "#ffffff"), (QPalette.ToolTipBase, NAVY), (QPalette.ToolTipText, "#ffffff"),
                      (QPalette.PlaceholderText, FAINT), (QPalette.Link, ACCENT_DARK)]:
        pal.setColor(role, QColor(col))
    for role in (QPalette.WindowText, QPalette.Text, QPalette.ButtonText):
        pal.setColor(QPalette.Disabled, role, QColor("#a3acb8"))
    return pal


def apply(app):
    app.setStyle("Fusion")
    f = QFont()
    f.setFamilies(UI_FAMILIES)
    f.setPointSizeF(9.5)
    app.setFont(f)
    app.setPalette(_palette())
    app.setStyleSheet(build_qss(_make_assets()))


def shadow(widget, blur=18, dy=2, alpha=26):
    """Soft drop shadow under a card."""
    from PySide6.QtWidgets import QGraphicsDropShadowEffect
    e = QGraphicsDropShadowEffect(widget)
    e.setBlurRadius(blur)
    e.setOffset(0, dy)
    e.setColor(QColor(15, 23, 42, alpha))
    widget.setGraphicsEffect(e)
    return e


# ---- icons: a gradient tile with a white pictogram -----------------------------------------------------------------------
def _draw(p: QPainter, kind: str, s: float):
    pen = QPen(QColor("white"), s * 0.07, Qt.SolidLine, Qt.RoundCap, Qt.RoundJoin)
    p.setPen(pen)
    p.setBrush(Qt.NoBrush)
    c = s / 2
    if kind == "assembly":            # staggered reads over a consensus line
        for i, (x0, x1) in enumerate([(.22, .62), (.38, .78), (.30, .70)]):
            y = s * (.30 + .14 * i)
            p.drawLine(QPointF(s * x0, y), QPointF(s * x1, y))
        p.drawLine(QPointF(s * .2, s * .74), QPointF(s * .8, s * .74))
    elif kind == "align":             # aligned bars
        for i, w in enumerate([.56, .42, .50, .34]):
            y = s * (.28 + .15 * i)
            p.drawLine(QPointF(s * .22, y), QPointF(s * (.22 + w), y))
    elif kind == "tree":
        p.drawLine(QPointF(s * .25, s * .5), QPointF(s * .45, s * .5))
        p.drawLine(QPointF(s * .45, s * .3), QPointF(s * .45, s * .7))
        p.drawLine(QPointF(s * .45, s * .3), QPointF(s * .75, s * .3))
        p.drawLine(QPointF(s * .45, s * .7), QPointF(s * .6, s * .7))
        p.drawLine(QPointF(s * .6, s * .6), QPointF(s * .6, s * .8))
        p.drawLine(QPointF(s * .6, s * .6), QPointF(s * .75, s * .6))
        p.drawLine(QPointF(s * .6, s * .8), QPointF(s * .75, s * .8))
    elif kind == "popgen":            # pie
        p.drawEllipse(QRectF(s * .24, s * .24, s * .52, s * .52))
        p.drawLine(QPointF(c, c), QPointF(c, s * .24))
        p.drawLine(QPointF(c, c), QPointF(s * .74, s * .6))
    elif kind == "blast":             # magnifier
        p.drawEllipse(QRectF(s * .24, s * .22, s * .38, s * .38))
        p.drawLine(QPointF(s * .56, s * .54), QPointF(s * .76, s * .76))
    elif kind == "primer":            # two opposing arrows
        p.drawLine(QPointF(s * .22, s * .38), QPointF(s * .70, s * .38)); p.drawLine(QPointF(s * .62, s * .30), QPointF(s * .70, s * .38)); p.drawLine(QPointF(s * .62, s * .46), QPointF(s * .70, s * .38))
        p.drawLine(QPointF(s * .30, s * .64), QPointF(s * .78, s * .64)); p.drawLine(QPointF(s * .38, s * .56), QPointF(s * .30, s * .64)); p.drawLine(QPointF(s * .38, s * .72), QPointF(s * .30, s * .64))
    elif kind == "barcode":
        for i in range(7):
            x = s * (.22 + .08 * i)
            p.drawLine(QPointF(x, s * .28), QPointF(x, s * .72))
    elif kind == "open":
        p.drawRoundedRect(QRectF(s * .22, s * .34, s * .56, s * .38), s * .05, s * .05)
        p.drawLine(QPointF(s * .22, s * .34), QPointF(s * .36, s * .26)); p.drawLine(QPointF(s * .36, s * .26), QPointF(s * .50, s * .34))
    elif kind == "save":
        p.drawRoundedRect(QRectF(s * .26, s * .26, s * .48, s * .48), s * .05, s * .05)
        p.drawRect(QRectF(s * .36, s * .26, s * .26, s * .14)); p.drawRect(QRectF(s * .36, s * .54, s * .28, s * .20))
    elif kind == "new":
        p.drawLine(QPointF(c, s * .26), QPointF(c, s * .74)); p.drawLine(QPointF(s * .26, c), QPointF(s * .74, c))
    elif kind == "folder":
        p.drawRoundedRect(QRectF(s * .22, s * .36, s * .56, s * .36), s * .04, s * .04)
        p.drawLine(QPointF(s * .22, s * .36), QPointF(s * .30, s * .28)); p.drawLine(QPointF(s * .30, s * .28), QPointF(s * .44, s * .28)); p.drawLine(QPointF(s * .44, s * .28), QPointF(s * .50, s * .36))
    elif kind == "run":
        path = QPainterPath(); path.moveTo(s * .38, s * .28); path.lineTo(s * .38, s * .72); path.lineTo(s * .74, c); path.closeSubpath()
        p.setBrush(QBrush(QColor("white"))); p.drawPath(path)
    elif kind == "home":
        p.drawPolyline(QPolygonF([QPointF(s * .22, s * .50), QPointF(c, s * .26), QPointF(s * .78, s * .50)]))
        p.drawPolyline(QPolygonF([QPointF(s * .30, s * .46), QPointF(s * .30, s * .74), QPointF(s * .70, s * .74), QPointF(s * .70, s * .46)]))
    elif kind == "classic":           # a small window
        p.drawRoundedRect(QRectF(s * .22, s * .26, s * .56, s * .48), s * .06, s * .06)
        p.drawLine(QPointF(s * .22, s * .40), QPointF(s * .78, s * .40))
        p.drawLine(QPointF(s * .32, s * .54), QPointF(s * .60, s * .54)); p.drawLine(QPointF(s * .32, s * .64), QPointF(s * .50, s * .64))
    elif kind == "gear":
        p.drawEllipse(QRectF(s * .34, s * .34, s * .32, s * .32))
        for k in range(8):
            a = k * math.pi / 4
            p.drawLine(QPointF(c + math.cos(a) * s * .21, c + math.sin(a) * s * .21), QPointF(c + math.cos(a) * s * .30, c + math.sin(a) * s * .30))


_cache = {}


def icon(kind: str, size: int = 48, color: str = ACCENT, flat: bool = False) -> QIcon:
    """kind: assembly, align, tree, popgen, blast, primer, barcode, open, save, new, folder, run, home, classic, gear.
    flat=True gives the bare pictogram in `color` (for menus); otherwise a rounded gradient tile in `color`."""
    key = (kind, size, color, flat)
    if key in _cache:
        return _cache[key]
    s = size * 2
    pm = QPixmap(s, s)
    pm.fill(Qt.transparent)
    p = QPainter(pm)
    p.setRenderHint(QPainter.Antialiasing)
    if flat:
        tmp = QPixmap(s, s)
        tmp.fill(Qt.transparent)
        q = QPainter(tmp)
        q.setRenderHint(QPainter.Antialiasing)
        _draw(q, kind, s)
        q.setCompositionMode(QPainter.CompositionMode_SourceIn)       # tint the white pictogram
        q.fillRect(tmp.rect(), QColor(color))
        q.end()
        p.drawPixmap(0, 0, tmp)
    else:
        g = QLinearGradient(0, 0, 0, s)
        g.setColorAt(0, QColor(color).lighter(120))
        g.setColorAt(1, QColor(color).darker(112))
        p.setPen(Qt.NoPen)
        p.setBrush(QBrush(g))
        p.drawRoundedRect(QRectF(s * .04, s * .04, s * .92, s * .92), s * .26, s * .26)
        _draw(p, kind, s)
    p.end()
    pm.setDevicePixelRatio(2.0)
    ic = QIcon(pm)
    _cache[key] = ic
    return ic


def logo_icon(path, size: int = 38):
    """Icon from a picture file (the logos in www/img of the Shiny program); None when the file is missing."""
    try:
        if path and Path(path).exists():
            pm = QPixmap(str(path))
            if not pm.isNull():
                return QIcon(pm.scaled(size * 2, size * 2, Qt.KeepAspectRatio, Qt.SmoothTransformation))
    except Exception:
        pass
    return None
