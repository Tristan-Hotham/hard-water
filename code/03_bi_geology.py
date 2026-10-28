"""Step 3. Britain and Ireland rock groups at about 1 km.

Great Britain and the Isle of Man: BGS Geology 625k bedrock (rock type text rcs_d; named igneous
intrusions count as crystalline rock). The island of Ireland: GSI / GSNI 1:500k bedrock (unit names),
used for all of Ireland so that the calibration stations and the map share one source there.
Superficial deposits (till, alluvium and so on) are not used; see the caveats.

Polygons are drawn on a fine grid (0.0025 degrees of latitude by 0.00375 of longitude) and summed into
the analysis grid of 0.01 by 0.015 degrees (about 1.1 by 1.0 km at 54 degrees north).

Writes
  data/units_bi.csv        one row per map unit with its group shares
  work/bi_geo.npz          group shares (7 x rows x cols), mapped fraction, source (1 BGS, 2 GSI)
  work/r_bi_geo.json
Run: python3 code/03_bi_geology.py
"""
import csv, glob, gzip, json, os, sys, zipfile
import numpy as np
from PIL import Image, ImageDraw
import shapefile
from pyproj import Transformer

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import *  # noqa

LON0, LON1, LAT0, LAT1 = -11.0, 2.0, 49.8, 61.0
DY, DX = 0.01, 0.015
SUP = 4
NYc, NXc = int(round((LAT1 - LAT0) / DY)), int(round((LON1 - LON0) / DX))
NYf, NXf = NYc * SUP, NXc * SUP


def topx(lon, lat):
    return (np.asarray(lon) - LON0) / DX * SUP, (LAT1 - np.asarray(lat)) / DY * SUP


units, ukey = [], {}


def add_unit(key, **kw):
    if key not in ukey:
        kw["id"] = len(units) + 1
        units.append(kw); ukey[key] = kw["id"]
    return ukey[key]


def ring_area(r):
    a = np.asarray(r)
    return abs(np.dot(a[:-1, 0], a[1:, 1]) - np.dot(a[1:, 0], a[:-1, 1])) / 2 if len(a) > 2 else 0


def paint(polys):
    im = Image.new("I", (NXf, NYf), 0)
    d = ImageDraw.Draw(im)
    polys.sort(key=lambda p: -ring_area(p[1][0]))
    for v, rings in polys:
        if len(rings[0]) >= 3:
            d.polygon(rings[0], fill=v)
        for h in rings[1:]:
            if len(h) >= 3:
                d.polygon(h, fill=0)
    return np.asarray(im, dtype=np.int32)


# ---------------------------------------------------------------------------
# BGS (Great Britain, Isle of Man; Northern Ireland is taken from GSNI below)
B = json.load(gzip.open(BGS))
bgs_polys = []
for f in B["collections"]["bgsgeology625kbedrock"]:
    p = f["properties"]
    t = p["rcs_d"]
    intr = p["lex"].startswith("UI") or p["rank"] == "INTRUSION"
    if intr and not any(w in t for w in ("LAVA", "TUFF")):
        sh = {"crys": 1.0}
    elif t.startswith("SEDIMENTARY ROCK CYCLES"):
        sh = {"clas": 0.75, "carb": 0.25}   # Carboniferous cyclic sequences with limestone beds (Yoredale, Strathclyde, Clackmannan types)
    else:
        sh = shares_from_text(t)
    uid = add_unit("BGS:" + t + (":intr" if intr else ""), src="BGS", text=clean(t.capitalize()) + (" (igneous intrusion)" if intr else ""), shares=sh)
    g = f["geometry"]
    parts = [g["coordinates"]] if g["type"] == "Polygon" else g["coordinates"] if g["type"] == "MultiPolygon" else []
    for rings in parts:
        rr = []
        for r in rings:
            a = np.asarray(r)
            x, y = topx(a[:, 0], a[:, 1])
            rr.append(list(zip(x.tolist(), y.tolist())))
        bgs_polys.append((uid, rr))
print("BGS polygons", len(bgs_polys))

# ---------------------------------------------------------------------------
# GSI / GSNI (island of Ireland)
GSI_OVERRIDE = {65: {"carb": 0.7, "clas": 0.3},   # "Calp": dark-grey argillaceous and cherty limestone and shale
                62: {"carb": 1.0}, 67: {"carb": 0.9, "clas": 0.1}, 77: {"carb": 0.85, "clas": 0.15}}
tmp = os.path.join(WORK, "gsi")
os.makedirs(tmp, exist_ok=True)
with zipfile.ZipFile(GSI_ZIP) as zf:
    for n in zf.namelist():
        if "_MS." in n and n.split(".")[-1] in ("shp", "shx", "dbf", "prj", "cpg", "CPG"):
            zf.extract(n, tmp)
shp = glob.glob(os.path.join(tmp, "**", "*_MS.shp"), recursive=True)[0]
itm = Transformer.from_crs("EPSG:2157", "EPSG:4326", always_xy=True)
sr = shapefile.Reader(shp)
gsi_polys = []
for sh_, rec in zip(sr.shapes(), sr.records()):
    un = rec["UNITNAME"]
    num = int(un.split(",")[0]); nm = un.split(",", 1)[1].strip()
    lith = nm.split(";", 1)[1] if ";" in nm else nm
    if num in GSI_OVERRIDE:
        shs = GSI_OVERRIDE[num]
    elif num <= 11:
        shs = {"volc": 1.0} if num == 10 else {"crys": 1.0}
    else:
        shs = shares_from_text(lith)
    uid = add_unit("GSI:%d" % num, src="GSI", text=clean(nm.replace(" - ", " to ")), shares=shs)
    pts = np.asarray(sh_.points)
    parts = list(sh_.parts) + [len(pts)]
    cur = None
    for a, b in zip(parts[:-1], parts[1:]):
        lon, lat = itm.transform(pts[a:b, 0], pts[a:b, 1])
        x, y = topx(lon, lat)
        r = list(zip(x.tolist(), y.tolist()))
        arr = np.asarray(r)
        s = np.dot(arr[:-1, 0], arr[1:, 1]) - np.dot(arr[1:, 0], arr[:-1, 1])
        if s > 0 or cur is None:   # clockwise exterior in y-down pixel space
            if cur:
                gsi_polys.append((uid, cur))
            cur = [r]
        else:
            cur.append(r)
    if cur:
        gsi_polys.append((uid, cur))
print("GSI polygons", len(gsi_polys))
for f in glob.glob(os.path.join(tmp, "**", "*"), recursive=True)[::-1]:
    (os.remove if os.path.isfile(f) else os.rmdir)(f)
os.rmdir(tmp) if os.path.isdir(tmp) else None

gsi = paint(gsi_polys)
bgs = paint(bgs_polys)
fine = np.where(gsi > 0, gsi, bgs)
src_f = np.where(gsi > 0, 2, np.where(bgs > 0, 1, 0)).astype(np.uint8)
del gsi, bgs

NG = len(GKEYS)
SH = np.zeros((len(units) + 1, NG), np.float32)
for u in units:
    for k, v in (u["shares"] or {}).items():
        SH[u["id"], GKEYS.index(k)] = v
f4 = SH[fine]                                    # (NYf, NXf, NG)
blk = f4.reshape(NYc, SUP, NXc, SUP, NG).sum((1, 3))
cnt = (fine > 0).reshape(NYc, SUP, NXc, SUP).sum((1, 3))
tot = blk.sum(-1)
shares = np.where(tot[..., None] > 0, blk / np.maximum(tot[..., None], 1e-9), 0).transpose(2, 0, 1)
mapped = cnt / SUP ** 2
src = np.maximum.reduce([src_f[i::SUP, j::SUP] for i in range(SUP) for j in range(SUP)])
np.savez_compressed(os.path.join(WORK, "bi_geo.npz"), shares=(shares * 250).round().astype(np.uint8),
                    mapped=(mapped * 250).round().astype(np.uint8), src=src, box=np.array([LON0, LON1, LAT0, LAT1, DY, DX]))
with open(os.path.join(DATA, "units_bi.csv"), "w", newline="", encoding="utf-8") as fh:
    w = csv.writer(fh)
    w.writerow(["id", "source", "rock_description"] + ["share_" + k for k in GKEYS])
    for u in units:
        s = u["shares"] or {}
        w.writerow([u["id"], u["src"], u["text"]] + [round(s.get(k, 0), 3) for k in GKEYS])
R = {"units": len(units), "units_bgs": sum(u["src"] == "BGS" for u in units), "units_gsi": sum(u["src"] == "GSI" for u in units),
     "grid": [NYc, NXc], "box": [LON0, LON1, LAT0, LAT1], "cell_deg": [DY, DX]}
print(R)
json.dump(R, open(os.path.join(WORK, "r_bi_geo.json"), "w"), indent=1)
