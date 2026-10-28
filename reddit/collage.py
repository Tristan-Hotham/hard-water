"""r/dataisbeautiful collage: reddit/hard_water.png (2400 x 3000), pale water blue and white.
Maps are drawn from the published page images (../data/world_model.png, ../data/bi_model.png) and ../data/stations.csv;
every number comes from ../data/results.json. Fonts: FONT_DIR (Fraunces, Source Serif 4, IBM Plex Mono).
Run from anywhere: FONT_DIR=/path/to/fonts python3 reddit/collage.py
"""
import json, os, sys
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
from matplotlib.patches import Patch, PathPatch
from matplotlib.path import Path
import matplotlib.patheffects as pe
from PIL import Image

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.join(HERE, "..")
sys.path.insert(0, os.path.join(ROOT, "code"))
from common import CLASSES, BREAKS, CLASS_RANGES  # noqa

FD = os.environ.get("FONT_DIR", "/tmp/claude-0/-home-claude/c232619d-81e0-53f3-9de6-41e973630753/scratchpad/fonts")
for f in ["FrauncesVar.ttf", "IBMPlexMono-Regular.ttf", "IBMPlexMono-Medium.ttf", "SourceSerif4-Regular.ttf", "SourceSerif4-Semibold.ttf"]:
    fm.fontManager.addfont(os.path.join(FD, f))
DISP = fm.FontProperties(fname=os.path.join(FD, "FrauncesVar.ttf"), weight=900)
SER, MONO = "Source Serif 4", "IBM Plex Mono"
BG, FG, MUT, DEEP, COOL, HOT = "#EEF5F8", "#10232E", "#4E5F68", "#0E3F66", "#1F7FAE", "#B5562A"
BGRGB = (238, 245, 248)
R = json.load(open(os.path.join(ROOT, "data", "results.json")))
GM = json.load(open(os.path.join(ROOT, "data", "grids.json")))
CC = [c[4] for c in CLASSES]
CRGB = np.array([[int(c[i:i + 2], 16) for i in (1, 3, 5)] for c in CC], float)
S, O, M = R["stories"], R["overall"], R["model"]


def rgb(name, hill=1.0):
    a = np.asarray(Image.open(os.path.join(ROOT, "data", f"{name}_model.png"))).astype(float)
    code, hs, b = a[..., 0], a[..., 1] / 255, a[..., 2]
    v = 10 ** ((code - 1) / 249 * GM["code_lmax"])
    out = np.empty(code.shape + (3,)); out[:] = BGRGB
    out[b == 1] = (214, 213, 208); out[b >= 2] = CRGB[np.digitize(v[b >= 2], BREAKS)]
    lo = b == 3                           # low confidence (coarse geology): faded and hatched
    out[lo] = out[lo] * .5 + np.array([227, 225, 218.]) * .5
    yy, xx = np.indices(b.shape)
    out[lo & ((xx + yy) % 6 == 0)] *= .72
    f = np.clip(1 - hill * (.77 - hs) * 1.25, .35, 1.25)
    out[b > 0] *= f[b > 0][:, None]
    return np.clip(out, 0, 255).astype(np.uint8)


def merc(lat):
    return np.log(np.tan(np.pi / 4 + np.radians(lat) / 2))


def to_px(name, lat, lon):
    g = GM[name]
    y = (merc(g["lat_top"]) - merc(np.asarray(lat))) / (merc(g["lat_top"]) - merc(g["lat_bot"])) * g["h"]
    x = (np.asarray(lon) - g["lon0"]) / (g["lon1"] - g["lon0"]) * g["w"]
    return x, y


fig = plt.figure(figsize=(12, 15), dpi=200, facecolor=BG)
fig.text(0.045, 0.977, "Why is some water hard?", fontproperties=DISP, fontsize=50, color=DEEP, va="top")
fig.text(0.045, 0.928, "Look at the rocks beneath it", fontproperties=DISP, fontsize=50, color=COOL, va="top")
fig.text(0.045, 0.872, f"{O['n']:,} river, lake and groundwater stations in {R['stations']['countries']} countries (dots), over a simple estimate\n"
         "from the rocks within 10 km and the rainfall (shading). Soft over granite, hard over limestone and chalk.",
         family=SER, fontsize=15, color=FG, va="top", linespacing=1.2)

# hero: world map, stations over the modelled layer
img = rgb("world")
_, y0 = to_px("world", 73, 0); _, y1 = to_px("world", -46, 0); y0, y1 = int(y0), int(y1)
x0, x1 = int(GM["world"]["w"] * 0.02), int(GM["world"]["w"] * 0.985)
crop = img[y0:y1, x0:x1]
hfrac = 0.95 * 12 * crop.shape[0] / crop.shape[1] / 15
top = 0.826
ax = fig.add_axes([0.025, top - hfrac, 0.95, hfrac]); ax.imshow(crop, interpolation="lanczos"); ax.axis("off")
st = pd.read_csv(os.path.join(ROOT, "data", "stations.csv"))
st["c"] = np.digitize(st["median"].values, BREAKS)
for t, mk in {"River": "o", "Lake": "s", "Reservoir": "s", "Groundwater": "^"}.items():
    d = st[st["wtype"] == t].sort_values("median")
    x, y = to_px("world", d["lat"].values, d["lon"].values)
    ax.scatter(x - x0, y - y0, s=5.5, marker=mk, c=[CC[i] for i in d["c"]], edgecolors=FG, linewidths=.18, zorder=3)
ax.set_xlim(0, crop.shape[1]); ax.set_ylim(crop.shape[0], 0)
# call-outs on the map
def call(lat, lon, text, dx, dy, ha="left"):
    x, y = to_px("world", lat, lon); x -= x0; y -= y0
    ax.annotate(text, (x, y), (x + dx, y + dy), fontsize=11.5, family=MONO, color=FG, ha=ha, va="center",
                arrowprops=dict(arrowstyle="-", color=FG, lw=.9, shrinkA=2, shrinkB=4),
                path_effects=[pe.withStroke(linewidth=3.5, foreground=BG)], zorder=5)
call(62, 15, f"Nordic stations\nmedian {S['scandinavia']['median']:.0f} mg/l", -330, -95, "right")
call(20.5, -89, f"Yucatán limestone\nmedian {S['yucatan']['median']:.0f} mg/l", 140, 120)
call(49.6, 1.8, f"Chalk, northern France\nmedian {S['northern_france_chalk']['median']:.0f} mg/l", -170, 150, "right")
call(58, -95, f"Canadian Shield\nmedian {S['canada']['shield_median']:.0f} mg/l", -150, 160, "right")
key = [Patch(fc=c, ec="#00000022", label=f"{l.lower()} {r}") for c, l, r in zip(CC, [c[1] for c in CLASSES], CLASS_RANGES)]
fig.legend(handles=key, loc="upper center", bbox_to_anchor=(0.5, top - hfrac - 0.004), ncol=6, frameon=False, handlelength=1.0, columnspacing=.9,
           prop={"family": MONO, "size": 9.6})
fig.text(0.5, top - hfrac - 0.026, "mg/l as CaCO3 (= ppm).  Circles rivers, squares lakes, triangles wells.  Hatched: low confidence (coarse geology).  Grey: too dry",
         family=MONO, fontsize=9.6, color=MUT, ha="center", va="top")
ylow = top - hfrac - 0.05
print("map bottom", top - hfrac)

# big numbers, two by two on the left
cv = M["world_cv"]["rocks_rain_and_type"]
nums = [
    (f"{R['by_crystalline'][3]['median']:.0f} mg/l", "#C1606A", f"median station where granite,\ngneiss and other crystalline\nrock cover over 80% nearby\n(against {R['by_crystalline'][0]['median']:.0f} under 20%)"),
    (f"{O['by_type']['Groundwater']['median']:.0f} mg/l", DEEP, f"median groundwater station,\nharder than rivers ({O['by_type']['River']['median']:.0f})\nand lakes ({O['by_type']['Lake']['median']:.0f}): it sits\nin the rock for longer"),
    (f"{cv['class_within1'] * 100:.0f}%", COOL, f"of stations the rock model\nputs within one class, tested\non places it never saw ({cv['class_exact'] * 100:.0f}%\nexactly right)"),
    (f"{S['yucatan']['model_map_median']:.0f}", HOT, f"mg/l: the model's guess for\nthe Yucatán (measured {S['yucatan']['median']:.0f}). The\nworld map calls that limestone\njust \"Eocene sedimentary\""),
]
for i, (big, c, txt) in enumerate(nums):
    x = 0.045 + (i % 2) * 0.27; y = ylow - (i // 2) * 0.112
    fig.text(x, y, big, fontproperties=DISP, fontsize=31, color=c, va="top")
    fig.text(x, y - 0.037, txt, family=SER, fontsize=12, color=FG, va="top", linespacing=1.1)

# drops: six places, filled with their class colour
dt = ylow - 0.252
fig.text(0.045, dt + 0.006, "Measured medians, mg/l", family=MONO, fontsize=11, color=MUT, va="bottom", weight="medium")
drops = [(S["canada"]["shield_median"], "Canadian\nShield"), (S["scandinavia"]["by_country"]["Norway"], "Norway"),
         (S["ireland"]["low_carb_median"], "Ireland,\nno limestone"), (S["ireland"]["high_carb_median"], "Ireland,\nlimestone"),
         (S["northern_france_chalk"]["median"], "N. France\nchalk"), (S["yucatan"]["median"], "Yucatán")]
_t = np.linspace(0, 2 * np.pi, 120)
dp = Path(np.column_stack([0.62 * np.sin(_t) * np.sin(_t / 2) ** 1.3, 0.62 * np.cos(_t) - 0.05]), closed=True)   # a drop, point at the top
dax = fig.add_axes([0.035, dt - 0.094, 0.53, 0.094]); dax.set_xlim(-.6, 5.6); dax.set_ylim(-1.55, 1.15); dax.axis("off")
for i, (v, lab) in enumerate(drops):
    c = CC[int(np.digitize(v, BREAKS))]
    p = dp.transformed(matplotlib.transforms.Affine2D().scale(.78, .78).translate(i, .2))
    dax.add_patch(PathPatch(p, fc=c, ec=FG, lw=1.2))
    dax.text(i, .08, f"{v:.0f}", ha="center", va="center", fontproperties=DISP, fontsize=15, color="#FFFFFF" if v > 150 else FG)
    dax.text(i, -.38, lab, ha="center", va="top", family=MONO, fontsize=9.5, color=FG, linespacing=1.05)

# Britain and Ireland inset
bimg = rgb("bi", hill=1.1)
gb = GM["bi"]
_, by0 = to_px("bi", 59.0, 0); _, by1 = to_px("bi", 50.0, 0)
bcrop = bimg[int(by0):int(by1), int(gb["w"] * 0.06):int(gb["w"] * 1.0)]
bh = ylow - 0.085 - 0.096; bw = bh * 15 / 12 * bcrop.shape[1] / bcrop.shape[0]
bax = fig.add_axes([0.975 - bw, 0.096, bw, bh]); bax.imshow(bcrop, interpolation="lanczos"); bax.axis("off")
ie = st[st["country"] == "Ireland"]
for t, mk in {"River": "o", "Lake": "s", "Reservoir": "s", "Groundwater": "^"}.items():
    d = ie[ie["wtype"] == t]
    x, y = to_px("bi", d["lat"].values, d["lon"].values)
    bax.scatter(x - int(gb["w"] * 0.06), y - int(by0), s=9, marker=mk, c=[CC[i] for i in d["c"]], edgecolors=FG, linewidths=.3, zorder=3)
bax.set_xlim(0, bcrop.shape[1]); bax.set_ylim(bcrop.shape[0], 0)
bt = R["bi_tees_exe"]; bc = M["bi_cv"]["rocks_type_and_world_rain_effect"]
fig.text(0.61, ylow + 0.0, "Britain and Ireland at 1 km", fontproperties=DISP, fontsize=21, color=DEEP, va="top")
fig.text(0.61, ylow - 0.027, f"from bedrock maps; no British stations in the\nopen archive, so calibrated on {R['stories']['ireland']['n']} Irish ones\n"
         f"({bc['class_within1'] * 100:.0f}% within one class). Chalk lights up:\n{R['bi_counties_hardest'][0]['name']} {R['bi_counties_hardest'][0]['median']:.0f} mg/l,"
         f" Highland {R['bi_counties_softest'][1]['median']:.0f}.", family=SER, fontsize=11.5, color=FG, va="top", linespacing=1.12)

# footer
fig.text(0.045, 0.083, "Monitoring stations are not taps, and the shading is an estimate, tested by spatial cross-validation on 5 degree squares.\n"
         "A window of rock is not a catchment: the model misses London, whose water comes from chalk and rivers far upstream.",
         family=MONO, fontsize=9.4, color=MUT, va="top", linespacing=1.25)
fig.text(0.045, 0.048, "Sources: UNEP GEMS/Water GFQA v3 open subset (CC BY 4.0); Macrostrat map tiles (CC BY 4.0); BGS Geology 625k (OGL v3);\n"
         "GSI and GSNI 1:500k bedrock (CC BY 4.0, OGL); CRU CL 2.0 rainfall; Natural Earth; AWS Terrain Tiles.\n"
         "Analysis and chart: Tristan Hotham, with Python and matplotlib.  Interactive: tristan-hotham.github.io/hard-water",
         family=MONO, fontsize=9.4, color=FG, va="top", linespacing=1.25)
out = os.path.join(HERE, "hard_water.png")
fig.savefig(out, dpi=200, facecolor=BG)
print(out)
