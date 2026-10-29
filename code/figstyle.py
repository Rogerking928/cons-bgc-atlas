"""本專案共用的繪圖設定（腸球菌生合成基因簇圖譜）。

沿用嚴重氣喘那篇圖表本的做法：DejaVu Sans、細軸線、無格線，
類別（物種、基因簇家族）**同時**用顏色與線紋區分，灰階列印仍分得開。

尺寸一律以**印出來的寬度**下去畫：正文單欄 180 mm = 7.09 in
（見全域筆記「圖要照印出來的寬度畫」）。存檔 700 dpi PNG ＋ 向量 PDF。
"""
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(ROOT, "out")

WIDTH_IN = 7.09          # 180 mm
DPI = 700

NAVY = "#31485d"         # S. epidermidis
RED = "#a4303f"          # S. lugdunensis
AMBER = "#b8791f"        # S. capitis
GREEN = "#4f7d5b"        # S. haemolyticus
PURPLE = "#6b5b95"       # S. hominis
TEAL = "#2f7d8c"         # S. saprophyticus
BROWN = "#7a5c3e"        # S. simulans
PLUM = "#9c6b8a"         # S. warneri
GREY = "#8a939c"
FAINT = "#c9d1d8"
INK = "#1f2933"

# ---------------------------------------------------------------- palette
# 一個顏色一個意思，整份文件通用。讀者在第一張圖學到「深藍＝S. epidermidis」，
# 之後就不能在別張圖看到深藍代表別的東西。灰階列印用線紋分辨。
SPECIES = {"epidermidis": NAVY, "lugdunensis": RED, "capitis": AMBER,
           "haemolyticus": GREEN, "hominis": PURPLE, "saprophyticus": TEAL,
           "simulans": BROWN, "warneri": PLUM, "other": GREY}
SPECIES_HATCH = {"epidermidis": "", "lugdunensis": "///", "capitis": "...",
                 "haemolyticus": "\\\\", "hominis": "xxx",
                 "saprophyticus": "---", "simulans": "+++",
                 "warneri": "ooo", "other": ""}
SPECIES_ORDER = ["epidermidis", "haemolyticus", "capitis", "lugdunensis",
                 "hominis", "saprophyticus", "simulans", "warneri", "other"]
REPLICON = {"chromosome": "#6f8ca3", "plasmid": AMBER}
REPLICON_HATCH = {"chromosome": "", "plasmid": "\\\\"}
# 方法對照（模型 vs 基準、crude vs 分層）一律灰階＋形狀，不用顏色
METHOD = {"main": "#2b333b", "other": "#9aa3ab", "faint": "#c9d1d8"}
METHOD_MARKER = {"main": "o", "other": "s"}
# 基因角色的顏色（座位圖用）
ROLE = {"structural": "#c4322f", "maturation": "#31485d", "export": "#b8791f",
        "immunity": "#4f7d5b", "regulation": "#6b5b95", "mobility": "#8a939c",
        "resistance": "#000000", "other": "#c9d1d8"}
KIND = {"NRPS/PKS": NAVY, "hybrid": PURPLE, "RiPP": AMBER}


def sp_color(sp):
    return SPECIES.get(sp, GREY)


def sp_hatch(sp):
    return SPECIES_HATCH.get(sp, "")


HATCH_NC = ""            # 第一類：實心
HATCH_C = "///"          # 第二類：斜線，灰階下仍可辨

plt.rcParams.update({
    "font.family": ["DejaVu Sans"],
    "font.size": 9,
    "axes.titlesize": 9.5,
    "axes.labelsize": 9,
    "legend.fontsize": 8.5,
    "xtick.labelsize": 8.5,
    "ytick.labelsize": 8.5,
    "axes.unicode_minus": False,
    "axes.edgecolor": "#5b6b7a",
    "axes.labelcolor": INK,
    "text.color": INK,
    "xtick.color": "#5b6b7a",
    "ytick.color": "#5b6b7a",
    "figure.facecolor": "white",
    "hatch.linewidth": 0.6,
})


def bare(ax, keep=("bottom", "left")):
    for side, sp in ax.spines.items():
        sp.set_visible(side in keep)
    ax.grid(False)
    return ax


def panel(ax, letter, dx=-34, dy=8):
    """粗體 panel 字母，用 offset points 定位。

    先前寫成 fraction × 72 的混合寫法，字母會落在標題文字上。offset points
    不隨軸寬變動，換圖尺寸也不會再撞。
    """
    ax.annotate(letter, xy=(0, 1), xycoords="axes fraction",
                xytext=(dx, dy), textcoords="offset points",
                fontsize=11, fontweight="bold", va="bottom", ha="left")


def italicise(fig):
    """基因符號與學名改斜體，字彙統一用 code/italics.py。

    逐一在畫的地方包一層會漏掉刻度標籤與圖例，所以跟 .docx 一樣做成出場前的
    一次性掃描：整段都該斜體的直接設 fontstyle（字型不變，比 mathtext 好看），
    只有夾雜的才退回 mathtext。

    刻度標籤要**連 formatter 一起改**：`set_xticklabels()` 留下的是一個
    formatter（3.10 是 FuncFormatter，不是 FixedFormatter），每次 draw 都會
    拿它的輸出重新寫回標籤，所以只改 Text 物件的話，下一次 draw（check() 自己
    就會 draw 一次）就把字洗回去、只剩 fontstyle 留著。Figure 4 的
    F001-ausA 就是這樣整條印成正體的。這裡不去猜 formatter 是哪一種，直接把
    原本那個包起來，讓它的輸出再過一次同樣的轉換。
    """
    import matplotlib.text as mtext
    from matplotlib.ticker import Formatter
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import italics

    def convert(raw):
        """(新字串, 底層字型是否斜體)；None 表示這段不必動。

        整段都該斜體的直接設 fontstyle（字型不變，比 mathtext 好看）；夾雜的
        交給 italics.mpl()，它會自己決定底層字型是正體還是斜體——基因符號帶
        括號時要用斜體當底，mathtext 不肯把括號畫成斜的。
        """
        if not raw or "$" in raw:
            return None
        parts = italics.split(raw)
        if not any(it for _, it in parts):
            return None
        if all(it for _, it in parts):
            return raw, True
        return italics.mpl(raw)

    class Delegating(Formatter):
        """原本的 formatter 外面包一層斜體轉換。

        不能只把它塞進 FuncFormatter：軸呼叫的是**外層**的 set_locs()，
        ScalarFormatter 沒拿到刻度就不知道要印什麼，數字軸會整排變空白
        （Figure 1 的 0/50/100 就這樣不見過）。這裡把每一個呼叫都轉交回去。
        """

        def __init__(self, inner):
            self.inner = inner

        def __call__(self, x, pos=None):
            got = convert(self.inner(x, pos))
            return got[0] if got else self.inner(x, pos)

        def set_locs(self, locs):
            self.inner.set_locs(locs)

        def get_offset(self):
            return self.inner.get_offset()

        def set_axis(self, axis):
            self.inner.set_axis(axis)
            super().set_axis(axis)

    n = 0
    for ax in fig.axes:
        for axis in (ax.xaxis, ax.yaxis):
            inner = axis.get_major_formatter()
            if isinstance(inner, Delegating):
                continue
            try:
                labels = inner.format_ticks(list(axis.get_majorticklocs()))
            except Exception:
                continue
            # 數字軸沒有東西要斜體，就不要動它
            if not any(convert(lab) for lab in labels):
                continue
            axis.set_major_formatter(Delegating(inner))
    for t in fig.findobj(mtext.Text):
        got = convert(t.get_text())
        if not got:
            continue
        n += 1
        label, base_italic = got
        if base_italic:
            t.set_fontstyle("italic")
        if label != t.get_text():
            t.set_text(label)
    return n


def save(fig, stem):
    """同時存 700 dpi PNG（放進 Word）與向量 PDF（投稿用）。"""
    os.makedirs(OUT, exist_ok=True)
    n_it = italicise(fig)
    check(fig, stem)
    png = os.path.join(OUT, f"{stem}.png")
    pdf = os.path.join(OUT, f"{stem}.pdf")
    fig.savefig(png, dpi=DPI, bbox_inches="tight", pad_inches=0.02)
    fig.savefig(pdf, bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)
    print(f"  fig -> {os.path.basename(png)} + .pdf"
          + (f"  ({n_it} italic)" if n_it else ""), flush=True)
    return png


_MATH = __import__("re").compile(r"\$\\(?:it|mathrm)\{([^}]*)\}\$")


def _ink_box(t, bb, r):
    """The box a mathtext label actually inks.

    matplotlib measures `$\\it{ausA}$` four pixels taller than the same word in
    plain text: a math box carries the ascent and descent the formula layout
    reserves, and none of it is ink. Two tick labels one line apart are then
    reported as overlapping by a pixel although the reader sees a clear gap.
    The label is re-measured with the markup stripped and the box shrunk to
    that size about its centre, so the gate keeps judging ink.
    """
    import matplotlib.text as mtext
    from matplotlib.transforms import Bbox
    raw = t.get_text()
    if "$" not in raw:
        return bb
    plain = _MATH.sub(r"\1", raw).replace("$", "")
    probe = mtext.Text(0, 0, plain, fontproperties=t.get_fontproperties(),
                       rotation=t.get_rotation())
    probe.set_figure(t.get_figure())
    try:
        pb = probe.get_window_extent(renderer=r)
    except Exception:
        return bb
    w, h = min(bb.width, pb.width), min(bb.height, pb.height)
    cx, cy = (bb.x0 + bb.x1) / 2, (bb.y0 + bb.y1) / 2
    return Bbox.from_extents(cx - w / 2, cy - h / 2, cx + w / 2, cy + h / 2)


def check(fig, stem=""):
    """Report overlapping or out-of-canvas text before the figure is saved.

    Eyeballing a 700 dpi panel does not catch a label sitting on top of a title
    (see the global note on checking overlap programmatically), so every figure
    goes through this on the way out.
    """
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    # A tick label's .axes is not set, so asking the text which axes it belongs
    # to drops every one of them and the gate passes a figure whose tick labels
    # sit on each other. The owner is looked up from the axes instead.
    owner = {}
    for ax_ in fig.axes:
        for t_ in list(ax_.get_xticklabels()) + list(ax_.get_yticklabels()):
            owner[id(t_)] = ax_
    items = []
    for t in fig.findobj(plt.Text):
        if not t.get_visible() or not (t.get_text() or "").strip():
            continue
        try:
            bb = t.get_window_extent(renderer=r)
        except Exception:
            continue
        bb = _ink_box(t, bb, r)
        if bb.width <= 1 or bb.height <= 1:
            continue
        ax = owner.get(id(t)) or getattr(t, "axes", None)
        if ax is None and t not in list(fig.texts):
            # an orphan tick label kept for a tick outside the view limits; it
            # belongs to no axes and is never drawn
            continue
        if id(t) in owner:
            ab = ax.get_window_extent(renderer=r)
            if (bb.x1 < ab.x0 - 60 or bb.x0 > ab.x1 + 60
                    or bb.y1 < ab.y0 - 60 or bb.y0 > ab.y1 + 60):
                continue
        items.append((t, bb))
    # text outside the canvas is fine: the figures are saved with a tight
    # bounding box, which grows to include it
    problems = []
    for i in range(len(items)):
        for j in range(i + 1, len(items)):
            a, b = items[i][1], items[j][1]
            dx = min(a.x1, b.x1) - max(a.x0, b.x0)
            dy = min(a.y1, b.y1) - max(a.y0, b.y0)
            if dx > 1 and dy > 1 and dx * dy > 12:
                problems.append(f"overlap {dx*dy:6.0f}px2: "
                                f"{items[i][0].get_text()[:28]!r} x {items[j][0].get_text()[:28]!r}")
    # text on top of a drawn shape is the other half of the problem, and the
    # text-against-text pass above is blind to it: a label centred between two
    # bars can sit squarely on one of them (2026-09-20, Figure 6F). Bars are
    # compared against the labels of their own axes; a label the author put
    # inside a bar on purpose is written in white, so white text is exempt.
    import matplotlib.patches as _mpatches
    from matplotlib.colors import to_rgba as _to_rgba
    for ax_ in fig.axes:
        bars = [q for q in ax_.patches
                if isinstance(q, _mpatches.Rectangle) and q.get_visible()
                and _to_rgba(q.get_facecolor())[3] > 0.05]
        if not bars:
            continue
        for t, bb in items:
            if (owner.get(id(t)) or getattr(t, "axes", None)) is not ax_:
                continue
            if _to_rgba(t.get_color())[:3] == (1.0, 1.0, 1.0):
                continue
            area = max(1.0, bb.width * bb.height)
            for q in bars:
                qb = q.get_window_extent(renderer=r)
                dx = min(bb.x1, qb.x1) - max(bb.x0, qb.x0)
                dy = min(bb.y1, qb.y1) - max(bb.y0, qb.y0)
                # a label drawn inside a cell on purpose sits wholly within it
                # (Figure S2 writes the amino-acid code in every coloured box);
                # the accident is a label that straddles an edge
                if (bb.x0 >= qb.x0 - 1 and bb.x1 <= qb.x1 + 1
                        and bb.y0 >= qb.y0 - 1 and bb.y1 <= qb.y1 + 1):
                    continue
                if dx > 1 and dy > 1 and dx * dy > 12 and dx * dy > 0.15 * area:
                    problems.append(f"text on a bar {dx*dy:6.0f}px2: "
                                    f"{t.get_text()[:34]!r}")
                    break

    # 文字壓在「線或散點」上，是第三種漏網的：Figure 4B 的計數落在 100% 的
    # 圓點上、Figure 5A 的圖例被 PR 曲線穿過，兩者文字沒撞到文字、也沒撞到
    # 長條，所以上面兩關都放行。比對的是該軸自己畫的資料點座標。
    for ax_ in fig.axes:
        leg_ = ax_.get_legend()
        skip = {id(h) for h in (leg_.get_lines() if leg_ else [])}
        curves = []
        for ln in ax_.get_lines():
            if not ln.get_visible() or id(ln) in skip:
                continue
            xy = ln.get_xydata()
            # 兩點的線就是 axvline／axhline 這種參考線，貼著它標字是刻意的
            if len(xy) <= 2 and (ln.get_marker() in (None, "", "None")):
                continue
            curves.append(ln.get_transform().transform(xy))
        for col in ax_.collections:
            off = col.get_offsets()
            if off is not None and len(off):
                curves.append(col.get_offset_transform().transform(off))
        if not curves:
            continue
        for t, bb in items:
            if (owner.get(id(t)) or getattr(t, "axes", None)) is not ax_:
                continue
            if _to_rgba(t.get_color())[:3] == (1.0, 1.0, 1.0):
                continue
            for pts in curves:
                inside = ((pts[:, 0] > bb.x0 + 1) & (pts[:, 0] < bb.x1 - 1)
                          & (pts[:, 1] > bb.y0 + 1) & (pts[:, 1] < bb.y1 - 1))
                if inside.sum() >= 2:
                    problems.append("text on plotted data: "
                                    f"{t.get_text()[:34]!r}")
                    break

    # 圖例伸出自己的格子本身沒問題（擺在軸下方就一定會寬一點），
    # 真正的毛病是它伸進**別的** panel 裡（Figure 6B 的三欄圖例壓到 C）。
    boxes = [(ax_, ax_.get_window_extent(renderer=r)) for ax_ in fig.axes]
    for ax_ in fig.axes:
        leg_ = ax_.get_legend()
        if leg_ is None or not leg_.get_visible():
            continue
        try:
            lb = leg_.get_window_extent(renderer=r)
        except Exception:
            continue
        for other, ob in boxes:
            if other is ax_:
                continue
            dx = min(lb.x1, ob.x1) - max(lb.x0, ob.x0)
            dy = min(lb.y1, ob.y1) - max(lb.y0, ob.y0)
            if dx > 1 and dy > 1:
                problems.append(f"legend reaches into another panel "
                                f"({dx*dy:.0f}px2)")
                break

        # 圖例壓在自己的軸線上（Figure 5C：錨點寫 upper center + (0.5, 1.06)，
        # 等於把圖例的上緣釘在軸框**內**，比格子寬的部分就落在 y 軸線上，
        # 紅色的 _ 標記被軸線穿過）。標記不是文字，前面每一關都看不到它。
        for side, sp in ax_.spines.items():
            if not sp.get_visible():
                continue
            sb = sp.get_window_extent(renderer=r)
            dx = min(lb.x1, sb.x1) - max(lb.x0, sb.x0)
            dy = min(lb.y1, sb.y1) - max(lb.y0, sb.y0)
            if dx > 1 and dy > 1 and dx * dy > 20:
                problems.append(f"legend sits on the {side} axis line "
                                f"({dx*dy:.0f}px2)")
                break

    if problems:
        print(f"  !! {stem or 'figure'}: {len(problems)} layout problems", flush=True)
        for p in problems[:12]:
            print("     -", p, flush=True)
    return problems
