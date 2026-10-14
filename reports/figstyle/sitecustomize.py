"""House figure style for the Overcharged reports (user, 2026-10-13): every label, title, legend entry and annotation starts with a capital letter,
and all fonts are larger. Loaded automatically by Python when this directory is on PYTHONPATH:

    PYTHONPATH=~/Physics/Overcharged/reports/figstyle python3 some_figure_script.py
    (or: reports/figstyle/run.sh some_figure_script.py [args])

Capitalisation only touches the first word (after an optional panel prefix "(a) "), and only when that word is an ordinary lower-case word of two or more
letters. Symbols and variables are left alone: "k = 6", "s = 0.01", "q ≈ 0.03", "λ", "$\\lambda=2$", "cos τ", numbers. "rms" becomes "RMS".
Font sizes: every numeric size (rc defaults and explicit fontsize=...) is multiplied by FIGSTYLE_SCALE (default 1.3).
"""
import os, re

SCALE = float(os.environ.get("FIGSTYLE_SCALE", "1.3"))
KEEP = {"cos", "sin", "tan", "exp", "log", "ln", "arccos", "arcsin", "arctan", "eps", "vs", "max", "min", "det", "tr", "mod", "e", "x", "y", "z"}
MAP = {"rms": "RMS", "mpcd": "MPCD"}
_PREFIX = re.compile(r"^(\s*\([a-z]\d?\)\s*)?(.*)$", re.S)
_WORD = re.compile(r"^([A-Za-z]+)")


def cap(s):
    if not isinstance(s, str) or not s:
        return s
    m = _PREFIX.match(s)
    pre, rest = m.group(1) or "", m.group(2)
    w = _WORD.match(rest)
    if not w:
        return s
    word = w.group(1)
    low = word.lower()
    if low in MAP and word == low:
        return pre + MAP[low] + rest[len(word):]
    if len(word) < 2 or low in KEEP or not word.islower():
        return s
    # a word directly followed by "_" or "^" or a digit is an identifier (e.g. R_eq, k02): leave it
    nxt = rest[len(word):len(word) + 1]
    if nxt in ("_", "^") or nxt.isdigit():
        return s
    return pre + word[0].upper() + rest[1:]


def _install():
    try:
        import matplotlib
        from matplotlib import axis as maxis, axes as maxes, figure as mfig, legend as mleg, font_manager as fm
    except Exception:
        return

    # 1. capitalised labels, titles, suptitles, annotations
    _set_label_text = maxis.Axis.set_label_text
    def set_label_text(self, label, *a, **k):
        return _set_label_text(self, cap(label), *a, **k)
    maxis.Axis.set_label_text = set_label_text

    _set_title = maxes.Axes.set_title
    def set_title(self, label, *a, **k):
        return _set_title(self, cap(label), *a, **k)
    maxes.Axes.set_title = set_title

    for cls in (mfig.Figure, getattr(mfig, "SubFigure", mfig.Figure)):
        _sup = cls.suptitle
        def suptitle(self, t, *a, _sup=_sup, **k):
            return _sup(self, cap(t), *a, **k)
        cls.suptitle = suptitle

    _text = maxes.Axes.text
    def text(self, x, y, s, *a, **k):
        return _text(self, x, y, cap(s), *a, **k)
    maxes.Axes.text = text

    _ftext = mfig.Figure.text
    def ftext(self, x, y, s, *a, **k):
        return _ftext(self, x, y, cap(s), *a, **k)
    mfig.Figure.text = ftext

    _leg_init = mleg.Legend.__init__
    def leg_init(self, *a, **k):
        _leg_init(self, *a, **k)
        for t in self.get_texts():
            t.set_text(cap(t.get_text()))
        tt = self.get_title()
        if tt.get_text():
            tt.set_text(cap(tt.get_text()))
    mleg.Legend.__init__ = leg_init

    # 2. larger fonts: scale every size set on a FontProperties (rc defaults, fontsize=..., legend prop, ticks)
    _set_size = fm.FontProperties.set_size
    def set_size(self, size):
        _set_size(self, size)
        try:
            self._size = float(self._size) * SCALE
        except Exception:
            pass
    fm.FontProperties.set_size = set_size
    fm.FontProperties.set_fontsize = set_size

    # 3. vector copies for journals: every savefig("x.png") also writes x.pdf (TrueType fonts embedded). Heavy artists (collections with many
    #    elements: dense scatter, director/quiver segments, streamplots, meshes) are rasterised in the PDF only, at >= 300 dpi, to keep files small.
    #    Disable with FIGSTYLE_PDF=0; threshold FIGSTYLE_RASTER_N (default 1500 elements).
    matplotlib.rcParams["pdf.fonttype"] = 42
    matplotlib.rcParams["ps.fonttype"] = 42
    from matplotlib import collections as mcoll
    NMAX = int(os.environ.get("FIGSTYLE_RASTER_N", "1500"))

    def _heavy(art):
        if isinstance(art, mcoll.QuadMesh):
            return True
        if isinstance(art, mcoll.Collection):
            try:
                n = max(len(art.get_offsets()), len(art.get_paths()))
            except Exception:
                return False
            return n > NMAX
        return False

    _savefig = mfig.Figure.savefig
    def savefig(self, fname, *a, **k):
        out = _savefig(self, fname, *a, **k)
        if os.environ.get("FIGSTYLE_PDF", "1") != "0" and isinstance(fname, (str, os.PathLike)) and str(fname).lower().endswith(".png"):
            changed = []
            for ax in self.get_axes():
                for art in ax.get_children():
                    if _heavy(art) and not art.get_rasterized():
                        art.set_rasterized(True); changed.append(art)
            k2 = dict(k); k2["dpi"] = max(float(k.get("dpi") or 0), 300.0); k2.pop("format", None)
            try:
                _savefig(self, str(fname)[:-4] + ".pdf", *a, **k2)
            finally:
                for art in changed:
                    art.set_rasterized(False)
        return out
    mfig.Figure.savefig = savefig

    # matplotlib < 3.7 bug: rasterised images with few colours are written as palette images with bit depth < 8, but the PNG predictor
    # parameters omit BitsPerComponent, so poppler/evince report "insufficient data for an image" (fixed upstream in 3.7). Patch the method.
    try:
        import inspect, textwrap
        from matplotlib.backends import backend_pdf as bpdf
        src = textwrap.dedent(inspect.getsource(bpdf.PdfFile._writeImg))
        old = "png = {'Predictor': 10, 'Colors': color_channels, 'Columns': width}"
        if old in src and "'BitsPerComponent': obj['BitsPerComponent']" not in src:
            src = src.replace(old, "png = {'Predictor': 10, 'Colors': color_channels, 'Columns': width, 'BitsPerComponent': obj['BitsPerComponent']}")
            ns = {}
            exec(compile(src, bpdf.__file__, "exec"), vars(bpdf), ns)
            bpdf.PdfFile._writeImg = ns["_writeImg"]
    except Exception:
        pass


_install()
