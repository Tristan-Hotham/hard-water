"""Step 6. Compact files for the interactive page.

Writes
  data/world_model.png    Web Mercator image of the world estimate (85 N to 60 S): R = hardness code, G = hill shade,
                          B = 2 estimate, 3 low-confidence estimate (coarse geology), 1 land without one (too dry or unmapped), 0 sea
  data/bi_model.png       the same for the sharper Britain and Ireland estimate
  data/page_stations.json station list for the map (compact arrays)
  data/places.json        GeoNames places with modelled values, for search and map names
  data/borders.json       simplified country borders (Natural Earth 1:50m)
  data/grids.json         where the two images sit on the map
Hardness code: 0 none, else 1 + round(log10(mg/l) / log10(3000) * 249), so mg/l = 10^((code - 1) / 249 * log10(3000)).
Run: python3 code/06_page_data.py
"""
import base64, gzip, io, json, math, os, sys
import numpy as np
import pandas as pd
from PIL import Image

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import *  # noqa

LMAX = math.log10(3000)


def code(v):
    c = np.where(np.isfinite(v), 1 + np.round(np.log10(np.clip(np.nan_to_num(v, nan=1), 1, 3000)) / LMAX * 249), 0)
    return c.astype(np.uint8)


def merc_y(lat):
    la = np.radians(np.clip(lat, -85, 85))
    return np.log(np.tan(np.pi / 4 + la / 2))


def inv_merc(y):
    return np.degrees(2 * np.arctan(np.exp(y)) - np.pi / 2)


def decode_tile(b64):
    a = np.asarray(Image.open(io.BytesIO(base64.b64decode(b64))).convert("RGB")).astype(np.float32)
    return a[..., 0] * 256 + a[..., 1] + a[..., 2] / 256 - 32768


def hillshade(z, cell_m, ex=1.6):
    """0..1 shade (1 = flat, lower = in shadow) from a north-west light."""
    gy, gx = np.gradient(z, cell_m)
    slope = np.arctan(np.hypot(gx, gy) * ex)
    aspect = np.arctan2(-gx, gy)
    az, alt = np.radians(315), np.radians(45)
    sh = np.sin(alt) * np.cos(slope) + np.cos(alt) * np.sin(slope) * np.cos(az - aspect)
    flat = np.sin(alt)
    return np.clip(sh / flat, 0, 1.3) / 1.3


T = json.load(gzip.open(TERRARIUM))


def mosaic(z, x0, x1, y0, y1):
    out = np.zeros(((y1 - y0 + 1) * 256, (x1 - x0 + 1) * 256), np.float32)
    for x in range(x0, x1 + 1):
        for y in range(y0, y1 + 1):
            k = f"{z}/{x}/{y}"
            if k in T["tiles"]:
                out[(y - y0) * 256:(y - y0 + 1) * 256, (x - x0) * 256:(x - x0 + 1) * 256] = decode_tile(T["tiles"][k])
    return out


GM = {}
# ---------------------------------------------------------------------------
# world image: 3072 px wide, rows from 85 N to 60 S in Web Mercator
W = 3072
yT, yB = merc_y(85), merc_y(-60)
H = int(round(W * (yT - yB) / (2 * math.pi)))
ys = yT - (np.arange(H) + 0.5) / H * (yT - yB)
lat = inv_merc(ys)
lon = -180 + (np.arange(W) + 0.5) / W * 360
ZP = np.load(os.path.join(WORK, "world_pred.npz")); PW = ZP["pred"].astype(np.float32)
ZG = np.load(os.path.join(WORK, "world_geo.npz")); LAND = ZG["land"]
ri = ((90 - lat) / 0.1).astype(int).clip(0, 1799); ci = ((lon + 180) / 0.1).astype(int).clip(0, 3599)
pv = PW[np.ix_(ri, ci)]; ld = LAND[np.ix_(ri, ci)]
R_ = code(pv)
cf = ZP["conf"].astype(np.float32)[np.ix_(ri, ci)]
CONF_MIN = json.load(open(os.path.join(WORK, "r_model.json")))["conf_min"]
B_ = np.where(np.isfinite(pv), np.where(np.nan_to_num(cf, nan=0) >= CONF_MIN, 2, 3), np.where(ld, 1, 0)).astype(np.uint8)
# hill shade from the zoom-4 terrain (4096 px wide world), sampled to the same pixels
dem = mosaic(4, 0, 15, 0, 15)
n4 = 4096
gy4 = (1 - np.log(np.tan(np.pi / 4 + np.radians(np.clip(lat, -85, 85)) / 2)) / math.pi) / 2 * n4
gx4 = (lon + 180) / 360 * n4
hs = hillshade(np.maximum(dem, 0), 2 * math.pi * 6371008.8 / n4 * 0.6, ex=9)
G_ = (hs[np.ix_(gy4.astype(int).clip(0, n4 - 1), gx4.astype(int).clip(0, n4 - 1))] * 255).astype(np.uint8)
G_[B_ == 0] = 255
Image.fromarray(np.dstack([R_, G_, B_])).save(os.path.join(DATA, "world_model.png"), optimize=True)
GM["world"] = {"w": W, "h": H, "lat_top": 85, "lat_bot": -60, "lon0": -180, "lon1": 180}
del dem, hs

# ---------------------------------------------------------------------------
# Britain and Ireland image
ZB = np.load(os.path.join(WORK, "bi_pred.npz")); PB = ZB["pred"].astype(np.float32)
LON0, LON1, LAT0, LAT1, DY, DX = ZB["box"]
WB = 1100
yT, yB = merc_y(LAT1), merc_y(LAT0)
HB = int(round(WB * (yT - yB) / math.radians(LON1 - LON0)))
latb = inv_merc(yT - (np.arange(HB) + 0.5) / HB * (yT - yB))
lonb = LON0 + (np.arange(WB) + 0.5) / WB * (LON1 - LON0)
rb = ((LAT1 - latb) / DY).astype(int).clip(0, PB.shape[0] - 1); cb = ((lonb - LON0) / DX).astype(int).clip(0, PB.shape[1] - 1)
pbv = PB[np.ix_(rb, cb)]
RB_ = code(pbv)
# land from the 0.1 degree land mask is too coarse here: land = has an estimate
BB_ = np.where(np.isfinite(pbv), 2, 0).astype(np.uint8)
dem7 = mosaic(7, 58, 69, 34, 45)
n7 = 256 * 2 ** 7
gy7 = (1 - np.log(np.tan(np.pi / 4 + np.radians(latb) / 2)) / math.pi) / 2 * n7 - 34 * 256
gx7 = (lonb + 180) / 360 * n7 - 58 * 256
hs7 = hillshade(np.maximum(dem7, 0), 2 * math.pi * 6371008.8 / n7 * 0.6, ex=4)
GB_ = (hs7[np.ix_(gy7.astype(int).clip(0, hs7.shape[0] - 1), gx7.astype(int).clip(0, hs7.shape[1] - 1))] * 255).astype(np.uint8)
GB_[BB_ == 0] = 255
Image.fromarray(np.dstack([RB_, GB_, BB_])).save(os.path.join(DATA, "bi_model.png"), optimize=True)
GM["bi"] = {"w": WB, "h": HB, "lat_top": float(LAT1), "lat_bot": float(LAT0), "lon0": float(LON0), "lon1": float(LON1)}
GM["code_lmax"] = LMAX
json.dump(GM, open(os.path.join(DATA, "grids.json"), "w"))
del dem7, hs7

# ---------------------------------------------------------------------------
# stations
st = pd.read_csv(os.path.join(DATA, "stations.csv"))
sm = pd.read_csv(os.path.join(DATA, "station_model.csv"))
d = st.merge(sm, on="id", how="left")
TYPES = ["River", "Lake", "Reservoir", "Groundwater", "Wetland"]
countries = sorted(d["country"].unique())
rows = []
for r in d.itertuples():
    rows.append([r.id, round(r.lat, 4), round(r.lon, 4), TYPES.index(r.wtype) if r.wtype in TYPES else 0, r.median, int(r.n),
                 int(str(r.first)[:4]), int(str(r.last)[:4]),
                 None if pd.isna(r.pred) else round(r.pred), None if pd.isna(r.pred_cv) else round(r.pred_cv),
                 None if pd.isna(r.pred_bi) else round(r.pred_bi), r.sname if isinstance(r.sname, str) else "",
                 r.water_body if isinstance(r.water_body, str) else "", countries.index(r.country)])
json.dump({"types": TYPES, "countries": countries,
           "cols": ["id", "lat", "lon", "type", "median", "samples", "first", "last", "model", "model_cv", "model_bi", "name", "water_body", "country"],
           "rows": rows}, open(os.path.join(DATA, "page_stations.json"), "w"), separators=(",", ":"), ensure_ascii=False)

# places
pl = pd.read_csv(os.path.join(DATA, "places.csv"), keep_default_na=False, na_values=[""])
pl = pl.sort_values("pop", ascending=False)
prow = [[r.name, r.cc, round(r.lat, 3), round(r.lon, 3), int(r.pop), None if pd.isna(r.world) else int(r.world), None if pd.isna(r.bi) else int(r.bi)] for r in pl.itertuples()]
json.dump(prow, open(os.path.join(DATA, "places.json"), "w"), separators=(",", ":"), ensure_ascii=False)

# borders: Natural Earth 1:50m country outlines, rounded and thinned
ne = json.load(open(os.path.join(NE, "ne_50m_admin_0_countries.geojson")))
lines = []
for f in ne["features"]:
    g = f["geometry"]
    parts = [g["coordinates"]] if g["type"] == "Polygon" else g["coordinates"]
    for rings in parts:
        r = rings[0]
        kept = []
        for x, y in r:
            q = [round(y, 2), round(x, 2)]
            if not kept or abs(q[0] - kept[-1][0]) + abs(q[1] - kept[-1][1]) >= 0.18:
                kept.append(q)
        if len(kept) >= 4:
            lines.append([c for q in kept for c in q])
json.dump(lines, open(os.path.join(DATA, "borders.json"), "w"), separators=(",", ":"))
for f in ["world_model.png", "bi_model.png", "page_stations.json", "places.json", "borders.json"]:
    print(f, os.path.getsize(os.path.join(DATA, f)))
