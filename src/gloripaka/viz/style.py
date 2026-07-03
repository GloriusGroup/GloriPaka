"""A single, central figure style so every manuscript figure looks identical.

All plots go through :class:`JACSStyle`, which owns the project colour palette, the
JACS two-column geometry, and the font/spine/dpi conventions. Nothing else in the
codebase should call ``plt.rcParams`` or pick colours ad hoc - that is what keeps
Figure 1 and Figure 7 visually comparable.

Palette source: ``modeling_strategyB/visualisation/color_palette.md``.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.figure import Figure


# --- brand palette (hex) ---------------------------------------------------- #
MAIN_COLORS: list[str] = [
    "#002554",  # 0 deep navy      (primary / reference)
    "#005689",  # 1 blue
    "#006e89",  # 2 teal-blue
    "#4a8882",  # 3 muted teal
    "#6baf9e",  # 4 green
    "#b9decc",  # 5 pale green
    "#ffcb05",  # 6 yellow         (accent / highlight)
    "#fedf63",  # 7 pale yellow
    "#585c7d",  # 8 slate          (neutral / secondary)
    "#d0b8be",  # 9 dusty rose     (contrast)
]

# 20-step sequential ramp (cream -> deep teal) for heatmaps / continuous encodings.
SEQUENTIAL_HEX: list[str] = [
    "#fffff5", "#edf7ea", "#dcefe0", "#cbe7d6", "#bbdfcd", "#acd6c4", "#9dcdbc",
    "#90c4b4", "#82bbad", "#76b2a6", "#6aa9a0", "#5e9f99", "#539694", "#488d8e",
    "#3e8389", "#337a84", "#29717f", "#1e677b", "#125e77", "#005573",
]

# ACS/JACS column geometry (inches).
SINGLE_COLUMN_IN: float = 3.33
DOUBLE_COLUMN_IN: float = 7.0


class JACSStyle:
    """Owns the manuscript's visual identity and hands out styled figures.

    Typical use::

        from gloripaka.viz import STYLE
        fig, ax = STYLE.figure("single")
        ax.scatter(x, y, color=STYLE.color(0))
        STYLE.save(fig, "e2_scatter")
    """

    def __init__(self) -> None:
        self.main = MAIN_COLORS
        self.cmap = LinearSegmentedColormap.from_list("gdp_seq", SEQUENTIAL_HEX, N=256)
        self.cmap_r = self.cmap.reversed()
        self._font = self._resolve_font()
        self._apply_rc()

    # -- fonts ---------------------------------------------------------------- #
    # JACS house style is Arial. We vendor Liberation Sans (metric-IDENTICAL to
    # Arial, Apache/GPL-licensed) under assets/fonts so every figure renders in an
    # Arial-equivalent face regardless of what the host machine has installed.
    _PREFERRED = ("Arial", "Liberation Sans", "Arimo", "Aptos", "Helvetica")

    @classmethod
    def _bundled_font_dir(cls) -> Path:
        from gloripaka.config import PATHS

        return PATHS.workspace / "assets" / "fonts"

    @classmethod
    def _resolve_font(cls) -> str:
        """Register bundled fonts, then return the best available Arial-family face.

        Falls back to DejaVu Sans only if nothing Arial-like is present (record it in
        the journal if that ever happens, since it changes figure typography).
        """
        font_dir = cls._bundled_font_dir()
        if font_dir.exists():
            for ttf in font_dir.glob("*.ttf"):
                try:
                    mpl.font_manager.fontManager.addfont(str(ttf))
                except Exception:
                    continue
        available = {f.name for f in mpl.font_manager.fontManager.ttflist}
        for name in cls._PREFERRED:
            if name in available:
                return name
        return "DejaVu Sans"

    @property
    def font(self) -> str:
        return self._font

    def _apply_rc(self) -> None:
        plt.rcParams.update(
            {
                "font.family": "sans-serif",
                "font.sans-serif": [self._font, "Aptos", "Helvetica", "DejaVu Sans"],
                "font.size": 8,
                "axes.titlesize": 8,
                "axes.labelsize": 8,
                "xtick.labelsize": 7,
                "ytick.labelsize": 7,
                "legend.fontsize": 7,
                "axes.linewidth": 0.8,
                "axes.spines.top": False,
                "axes.spines.right": False,
                "axes.grid": False,
                "xtick.direction": "out",
                "ytick.direction": "out",
                "xtick.major.width": 0.8,
                "ytick.major.width": 0.8,
                "lines.linewidth": 1.2,
                "lines.markersize": 4,
                "figure.dpi": 150,
                "savefig.dpi": 400,
                "savefig.bbox": "tight",
                "savefig.pad_inches": 0.02,
                "pdf.fonttype": 42,  # embed TrueType (editable text in Illustrator)
                "ps.fonttype": 42,
            }
        )

    # -- colours -------------------------------------------------------------- #
    def color(self, i: int) -> str:
        """Cycle through the main palette by index."""
        return self.main[i % len(self.main)]

    def discrete(self, n: int) -> list[str]:
        """``n`` evenly-spaced colours sampled from the sequential ramp."""
        if n <= 1:
            return [SEQUENTIAL_HEX[len(SEQUENTIAL_HEX) // 2]]
        idx = [round(k * (len(SEQUENTIAL_HEX) - 1) / (n - 1)) for k in range(n)]
        return [SEQUENTIAL_HEX[i] for i in idx]

    # -- figures -------------------------------------------------------------- #
    def figure(
        self,
        width: str = "single",
        height: float | None = None,
        nrows: int = 1,
        ncols: int = 1,
        **kwargs,
    ) -> tuple[Figure, object]:
        """Create a correctly-sized figure + axes for one/two-column placement.

        Args:
            width: ``"single"`` (3.33 in) or ``"double"`` (7.0 in).
            height: figure height in inches; defaults to a pleasant aspect per width.
        """
        w = {"single": SINGLE_COLUMN_IN, "double": DOUBLE_COLUMN_IN}[width]
        if height is None:
            height = 2.6 if width == "single" else 3.2
        fig, axes = plt.subplots(nrows, ncols, figsize=(w, height), **kwargs)
        return fig, axes

    def save(self, fig: Figure, name: str, formats: tuple[str, ...] = ("pdf", "png")) -> list[Path]:
        """Save a figure into the workspace ``figures/`` dir in each format."""
        from gloripaka.config import PATHS

        PATHS.figures_dir.mkdir(parents=True, exist_ok=True)
        out: list[Path] = []
        for fmt in formats:
            p = PATHS.figures_dir / f"{name}.{fmt}"
            fig.savefig(p, format=fmt)
            out.append(p)
        return out


# Import-time singleton: importing the module applies the rcParams once.
STYLE = JACSStyle()
