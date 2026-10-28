"""Step 8. header.png (1200 x 630 social card) and the static figures for the paper (release/paper/f_*.png, 200 dpi,
white) and blog (release/blog/b_*.png, 1520 px wide, light background). Maps are drawn from the published page images
(data/world_model.png, data/bi_model.png) and the station list; numbers from data/results.json.
Fonts: FONT_DIR holding FrauncesVar.ttf, SourceSerif4-Regular.ttf, IBMPlexMono-Regular.ttf, IBMPlexMono-Medium.ttf.
Run: FONT_DIR=/path/to/fonts python3 code/08_figures.py
"""
import json, math, os, sys
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from PIL import Image

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.join(HERE, "..")
sys.path.insert(0, HERE)
from common import CLASSES, BREAKS, CLASS_RANGES, GROUPS  # noqa

FD = os.environ.get("FONT_DIR", "/tmp/claude-0/-home-claude/c232619d-81e0-53f3-9de6-41e973630753/scratchpad/fonts")
for f in ["FrauncesVar.ttf", "IBMPlexMono-Regular.ttf", "IBMPlexMono-Medium.ttf", "SourceSerif4-Regular.ttf"]:
    fm.fontManager.addfont(os.path.join(FD, f))
DISP = fm.FontProperties(fname=os.path.join(FD, "FrauncesVar.ttf"), weight=900)
plt.rcParams["font.family"] = "IBM Plex Mono"
R = json.load(open(os.path.join(ROOT, "data", "results.json")))
GM = json.load(open(os.path.join(ROOT, "data", "grids.json")))
PO = os.path.join(ROOT, "release", "paper"); BO = os.path.join(ROOT, "release", "blog")
os.makedirs(PO, exist_ok=True); os.makedirs(BO, exist_ok=True)
FG, MUT, RULE, DEEP, HOT, COOL = "#10232E", "#5A6B74", "#D3DEE3", "#0E3F66", "#B5562A", "#1F7FAE"
BLOGBG = "#F4F7F8"
CC = [c[4] for c in CLASSES]
CRGB = np.array([[int(c[i:i + 2], 16) for i in (1, 3, 5)] for c in CC], float)
LAND, NOEST = np.array([227, 225, 218.]), np.array([214, 212, 206.])
LMAX = GM["code_lmax"]


def cls(v):
    return np.digitize(v, BREAKS)


def rgb(name, sea=(255, 255, 255), model=True, hill=1.0):
    a = np.asarray(Image.open(os.path.join(ROOT, "data", f"{name}_model.png"))).astype(float)
    code, hs, b = a[..., 0], a[..., 1] / 255, a[..., 2]
    v = 10 ** ((code - 1) / 249 * LMAX)
    out = np.empty(code.shape + (3,)); out[:] = sea
    out[b == 1] = NOEST if model else LAND; out[b >= 2] = LAND
    if model:
        m = b >= 2
        out[m] = CRGB[cls(v[m])]
        lo = b == 3                       # low confidence (coarse geology): faded halfway to land and hatched
        out[lo] = out[lo] * .5 + LAND * .5
        yy, xx = np.indices(b.shape)
        out[lo & ((xx + yy) % 6 == 0)] *= .72
    f = np.clip(1 - hill * (.77 - hs) * 1.25, .35, 1.25)
    land = b > 0
    out[land] = out[land] * f[land][:, None]
    return np.clip(out, 0, 255).astype(np.uint8)


def merc(lat):
    return np.log(np.tan(np.pi / 4 + np.radians(lat) / 2))


def to_px(name, lat, lon):
    g = GM[name]
    y = (merc(g["lat_top"]) - merc(np.asarray(lat))) / (merc(g["lat_top"]) - merc(g["lat_bot"])) * g["h"]
    x = (np.asarray(lon) - g["lon0"]) / (g["lon1"] - g["lon0"]) * g["w"]
    return x, y


st = pd.read_csv(os.path.join(ROOT, "data", "stations.csv"))
sm = pd.read_csv(os.path.join(ROOT, "data", "station_model.csv"))
st = st.merge(sm, on="id", how="left")
st["c"] = cls(st["median"].values)
MARK = {"River": "o", "Lake": "s", "Reservoir": "s", "Groundwater": "^", "Wetland": "D"}


def style(ax, grid="y"):
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    ax.spines["left"].set_color(MUT); ax.spines["bottom"].set_color(MUT)
    ax.tick_params(colors=FG, length=0, labelsize=11)
    if grid: ax.grid(axis=grid, color=RULE, lw=.8); ax.set_axisbelow(True)


def save(fig, name, blog):
    if blog: fig.savefig(os.path.join(BO, f"b_{name}.png"), dpi=160, facecolor=BLOGBG)
    else: fig.savefig(os.path.join(PO, f"f_{name}.png"), dpi=200, facecolor="white")
    plt.close(fig)


def bgc(blog):
    plt.rcParams["axes.facecolor"] = BLOGBG if blog else "white"
    return BLOGBG if blog else "white"


KEY = [Patch(fc=c, ec="#00000022", label=f"{l} ({r})") for c, l, r in zip(CC, [c[1] for c in CLASSES], CLASS_RANGES)]
SHAPES = [Line2D([], [], marker="o", ls="", mfc="#86C1D9", mec=FG, ms=7, label="River"), Line2D([], [], marker="s", ls="", mfc="#86C1D9", mec=FG, ms=7, label="Lake or reservoir"),
          Line2D([], [], marker="^", ls="", mfc="#86C1D9", mec=FG, ms=8, label="Groundwater")]


def world_window(img, lat_top=78, lat_bot=-56):
    g = GM["world"]
    _, y0 = to_px("world", lat_top, 0); _, y1 = to_px("world", lat_bot, 0)
    return img[int(y0):int(y1)], int(y0)


def plot_stations(ax, name, df, y0=0, size=10, edge=.25):
    for t, mk in MARK.items():
        d = df[df["wtype"] == t].sort_values("median")
        x, y = to_px(name, d["lat"].values, d["lon"].values)
        ax.scatter(x, y - y0, s=size * (1.25 if mk == "^" else 1), c=[CC[i] for i in d["c"]], marker=mk, edgecolors=FG, linewidths=edge, zorder=3)


def fig_world_measured(blog):
    bg = BLOGBG if blog else "white"
    img, y0 = world_window(rgb("world", sea=[int(bg[i:i + 2], 16) for i in (1, 3, 5)] if blog else (255, 255, 255), model=False))
    fig = plt.figure(figsize=(9.5, 6.3) if blog else (10, 6.4), facecolor=bgc(blog))
    ax = fig.add_axes([0, .16, 1, .84]); ax.imshow(img, interpolation="lanczos"); ax.axis("off")
    plot_stations(ax, "world", st, y0, size=6 if blog else 5, edge=.15)
    ax.set_xlim(0, img.shape[1]); ax.set_ylim(img.shape[0], 0)
    fig.legend(handles=KEY, loc="lower left", ncol=3, frameon=False, fontsize=10, bbox_to_anchor=(0.01, 0.0), handlelength=1.2, columnspacing=1.2)
    save(fig, "world_measured", blog)


def fig_world_model(blog):
    bg = BLOGBG if blog else "white"
    img, y0 = world_window(rgb("world", sea=[int(bg[i:i + 2], 16) for i in (1, 3, 5)] if blog else (255, 255, 255)))
    fig = plt.figure(figsize=(9.5, 6.3) if blog else (10, 6.4), facecolor=bgc(blog))
    ax = fig.add_axes([0, .16, 1, .84]); ax.imshow(img, interpolation="lanczos"); ax.axis("off")
    fig.legend(handles=KEY + [Patch(fc="#C9DCE4", hatch="////", ec="#86A9BB", label="low confidence (coarse geology map)"), Patch(fc="#D6D4CE", label="no estimate (under 200 mm of rain, or unmapped)")], loc="lower left", ncol=3, frameon=False, fontsize=9.5,
               bbox_to_anchor=(0.01, 0.0), handlelength=1.2, columnspacing=1.0)
    save(fig, "world_model", blog)


def fig_britain(blog):
    bg = BLOGBG if blog else "white"
    img = rgb("bi", sea=[int(bg[i:i + 2], 16) for i in (1, 3, 5)] if blog else (255, 255, 255), hill=1.1)
    fig = plt.figure(figsize=(9.5, 10.6) if blog else (8.2, 9.4), facecolor=bgc(blog))
    ax = fig.add_axes([0, .1, .74, .9]); ax.imshow(img, interpolation="lanczos"); ax.axis("off")
    ie = st[st["country"] == "Ireland"]
    plot_stations(ax, "bi", ie, 0, size=22, edge=.5)
    for t in R["bi_towns"]:
        if t["name"] in ("London", "Cambridge", "Norwich", "Birmingham", "Manchester", "Glasgow", "Edinburgh", "Cardiff", "Dublin", "Galway", "Cork", "Belfast", "Aberdeen", "Hull", "Exeter", "Brighton", "Inverness"):
            pass
    import zipfile, io
    pl = pd.read_csv(os.path.join(ROOT, "data", "places.csv"), keep_default_na=False, na_values=[""])
    names = ["London", "Cambridge", "Norwich", "Birmingham", "Manchester", "Leeds", "Glasgow", "Edinburgh", "Cardiff", "Dublin", "Galway", "Cork", "Belfast", "Aberdeen", "Exeter", "Brighton", "Inverness", "Newcastle upon Tyne"]
    for n in names:
        g = pl[(pl["name"] == n) & pl["cc"].isin(["GB", "IE"])].sort_values("pop", ascending=False)
        if not len(g): continue
        r = g.iloc[0]; x, y = to_px("bi", r["lat"], r["lon"])
        ax.plot(x, y, "o", ms=3.5, color=FG, zorder=4)
        dy = -14 if n == "Edinburgh" else 3
        ax.text(x + 8, y + dy, n.replace(" upon Tyne", ""), fontsize=10, color=FG, zorder=4,
                path_effects=[__import__("matplotlib.patheffects", fromlist=["x"]).withStroke(linewidth=2.5, foreground=bg)])
    lg = fig.add_axes([.70, .3, .29, .5]); lg.axis("off")
    lg.legend(handles=KEY[::-1], loc="upper left", frameon=False, fontsize=9.2, handlelength=1.2, title="Modelled river water\nmg/l as CaCO3", title_fontsize=10.5, alignment="left")
    lg2 = fig.add_axes([.72, .12, .27, .2]); lg2.axis("off")
    lg2.legend(handles=SHAPES, loc="upper left", frameon=False, fontsize=10, title="Measured, Ireland", title_fontsize=10.5, alignment="left")
    save(fig, "britain", blog)


def fig_rocks(blog):
    fig, axs = plt.subplots(1, 2, figsize=(9.5, 4.4) if blog else (10, 4.2), facecolor=bgc(blog), sharex=True)
    for ax, rows, col, t in zip(axs, (R["by_crystalline"], R["by_carbonate"]), ("#C1606A", "#2B6F9E"), ("Share of crystalline rock nearby", "Share of carbonate rock nearby")):
        for i, r in enumerate(rows):
            y = len(rows) - 1 - i
            ax.plot([r["q25"], r["q75"]], [y, y], lw=8, color=col, alpha=.3, solid_capstyle="round")
            ax.plot(r["median"], y, "o", ms=10, color=col)
            ax.text(r["median"], y + .28, f"{r['median']:.0f}", ha="center", fontsize=10.5, color=FG)
        ax.set_yticks(range(len(rows))[::-1]); ax.set_yticklabels([f"{r['lo'] * 100:.0f} to {r['hi'] * 100:.0f}%\n{r['n']:,} stations" for r in rows], fontsize=10)
        ax.set_xscale("log"); ax.set_xlim(8, 700); ax.set_xticks([10, 20, 50, 100, 200, 500]); ax.set_xticklabels(["10", "20", "50", "100", "200", "500"])
        style(ax, "x"); ax.set_title(t, fontsize=12, color=FG, loc="left"); ax.set_xlabel("median station hardness, mg/l (log scale)", fontsize=10, color=MUT)
        ax.set_ylim(-.6, len(rows) - .4)
    fig.tight_layout()
    save(fig, "rocks", blog)


def fig_check(blog):
    d = st[st["pred_cv"].notna()]
    fig, ax = plt.subplots(figsize=(7.6, 7.2) if blog else (6.6, 6.2), facecolor=bgc(blog))
    lo, hi = 2, 3000
    xs = np.array([lo, hi])
    ax.fill_between(xs, xs / 2, xs * 2, color=COOL, alpha=.08, lw=0)
    for b in BREAKS:
        ax.axvline(b, color=RULE, lw=.8, zorder=0); ax.axhline(b, color=RULE, lw=.8, zorder=0)
    ax.scatter(d["pred_cv"].clip(lo, hi), d["median"].clip(lo, hi), s=4, c=[CC[i] for i in d["c"]], alpha=.7, linewidths=.25, edgecolors="#5A6B7466")
    ax.plot(xs, xs, color=FG, lw=1)
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
    tk = [5, 10, 50, 100, 300, 1000]
    ax.set_xticks(tk); ax.set_xticklabels([f"{t:,}" for t in tk]); ax.set_yticks(tk); ax.set_yticklabels([f"{t:,}" for t in tk])
    style(ax, None); ax.set_xlabel("predicted with its 5 degree square hidden, mg/l", fontsize=11, color=MUT); ax.set_ylabel("measured station median, mg/l", fontsize=11, color=MUT)
    m = R["model"]["world_cv"]["rocks_rain_and_type"]
    ax.text(.03, .97, f"{m['n']:,} stations\nright class {m['class_exact'] * 100:.0f}%\nwithin one class {m['class_within1'] * 100:.0f}%\nshaded: within a factor of 2", transform=ax.transAxes, va="top", fontsize=10.5, color=FG)
    fig.tight_layout()
    save(fig, "check", blog)


def fig_coef(blog):
    M = R["model"]; G = {g[0]: g for g in GROUPS}
    lab = {"carb": "Carbonate (limestone, chalk, dolomite)", "evap": "Evaporites (gypsum, salt)", "clas": "Sandstone, mudstone, shale",
           "unc": "Loose young sediments", "sedu": "Sedimentary, type not stated", "volc": "Volcanic (basalt, tuff)", "crys": "Crystalline (granite, gneiss)"}
    w = sorted(M["world_typical"], key=lambda k: M["world_typical"][k])
    bl = {"carb": "Carbonate and evaporite", "clas": "Sandstone, mudstone, other sediments", "volc": "Volcanic", "crys": "Crystalline"}
    bcol = {"carb": "#2B6F9E", "clas": "#C9A66B", "volc": "#8E6C8A", "crys": "#C1606A"}
    b = sorted(M["bi_typical"], key=lambda k: M["bi_typical"][k])
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(9.5, 6.6) if blog else (9, 6.2), facecolor=bgc(blog), gridspec_kw={"height_ratios": [7, 4]}, sharex=True)
    for ax, keys, vals, labs, cols, t in ((a1, w, M["world_typical"], lab, {k: G[k][2] for k in w}, "World model, seven rock groups"),
                                         (a2, b, M["bi_typical"], bl, bcol, "Britain and Ireland model, four groups")):
        ys = range(len(keys))
        ax.barh(list(ys), [vals[k] for k in keys], color=[cols[k] for k in keys], height=.62)
        for y, k in zip(ys, keys):
            ax.text(vals[k] * 1.06, y, f"{vals[k]:.0f}", va="center", fontsize=10.5, color=FG)
        ax.set_yticks(list(ys)); ax.set_yticklabels([labs[k] for k in keys], fontsize=10.5)
        ax.set_xscale("log"); ax.set_xlim(5, 1000); style(ax, "x"); ax.set_title(t, fontsize=11.5, color=FG, loc="left")
    a2.set_xticks([5, 10, 20, 50, 100, 200, 500]); a2.set_xticklabels(["5", "10", "20", "50", "100", "200", "500"])
    a2.set_xlabel("typical hardness of river water from that rock alone, mg/l\nat 1,000 mm of rain a year (log scale)", fontsize=10, color=MUT)
    fig.tight_layout()
    save(fig, "coef", blog)


def fig_countries(blog):
    c = pd.DataFrame(R["countries_measured"]); c = c[c["stations"] >= 20].sort_values("median")
    fig, ax = plt.subplots(figsize=(9.5, 8.2) if blog else (8.6, 7.6), facecolor=bgc(blog))
    ax.barh(range(len(c)), c["median"], color=[CC[i] for i in cls(c["median"].values)], height=.7, edgecolor="#00000022")
    for y, (_, r) in enumerate(c.iterrows()):
        ax.text(r["median"] + 4, y, f"{r['median']:.0f}  ({r['stations']:,})", va="center", fontsize=9.5, color=FG)
    ax.set_yticks(range(len(c))); ax.set_yticklabels(c["country"], fontsize=10.5)
    for b in BREAKS: ax.axvline(b, color=RULE, lw=.8, zorder=0)
    style(ax, None); ax.set_xlim(0, 400); ax.set_xlabel("median of station medians, mg/l (stations in brackets)", fontsize=10.5, color=MUT)
    ax.set_ylim(-.6, len(c) - .4)
    fig.tight_layout()
    save(fig, "countries", blog)


def header():
    fig = plt.figure(figsize=(12, 6.3), dpi=100, facecolor=BLOGBG)
    img, y0 = world_window(rgb("world", sea=(244, 247, 248)), 74, -50)
    ax = fig.add_axes([0.40, 0.05, 0.6, 0.9]); ax.imshow(img, interpolation="lanczos"); ax.axis("off")
    ax.set_xlim(img.shape[1] * .02, img.shape[1] * .98)
    fig.text(0.045, 0.86, "DATA STUDY  ·  WATER", family="IBM Plex Mono", fontsize=12, color=MUT)
    fig.text(0.04, 0.82, "Hard", fontproperties=DISP, fontsize=70, color=DEEP, va="top")
    fig.text(0.04, 0.635, "water", fontproperties=DISP, fontsize=70, color=COOL, va="top")
    o = R["overall"]
    fig.text(0.045, 0.43, f"The world's water hardness from\n{o['n']:,} river, lake and groundwater\nstations, and an estimate from the\nrocks and rainfall where none exist.",
             family="Source Serif 4", fontsize=15.5, color=FG, va="top", linespacing=1.3)
    for i, c in enumerate(CC):
        fig.patches.append(plt.Rectangle((0.045 + i * 0.045, 0.115), 0.043, 0.025, transform=fig.transFigure, color=c, ec="#00000022"))
    fig.text(0.045, 0.075, "soft", family="IBM Plex Mono", fontsize=10, color=MUT); fig.text(0.312, 0.075, "very hard", family="IBM Plex Mono", fontsize=10, color=MUT, ha="right")
    fig.text(0.045, 0.03, "Dr Tristan Hotham · Independent researcher", family="IBM Plex Mono", fontsize=11, color=MUT)
    fig.savefig(os.path.join(ROOT, "header.png"), dpi=100, facecolor=BLOGBG)
    plt.close(fig)


if __name__ == "__main__":
    header()
    for blog in (False, True):
        fig_world_measured(blog); fig_world_model(blog); fig_britain(blog); fig_rocks(blog); fig_check(blog); fig_coef(blog); fig_countries(blog)
    print(sorted(os.listdir(PO)), sorted(os.listdir(BO)))
