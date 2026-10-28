"""Step 2. World rock groups from the Macrostrat carto tiles at zoom 4.

Each zoom-4 vector tile is drawn at 512 by 512 pixels (about 4.9 km at the equator). Every distinct map
unit gets shares of the seven rock groups from its lithology text (code/common.py). The pixels are then
summed into a 0.1 degree latitude and longitude grid: for each cell, the mapped share of each group.
Submarine units and anything outside Natural Earth 1:10m land are left out.

Writes
  data/units_world.csv     one row per Macrostrat unit: name, lithology, age, group shares, mapped pixels
  work/world_geo.npz       group shares (7 x 1800 x 3600, 0.1 degree), mapped fraction and land mask
  work/r_world_geo.json    numbers for results.json
Run: python3 code/02_world_geology.py
"""
import base64, csv, gzip, json, math, os, re, sys
import numpy as np
from PIL import Image, ImageDraw
import mapbox_vector_tile as mvt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import *  # noqa

TZ = 4; RES = 512; EXT = 4096; NT = 2 ** TZ; S = RES / EXT
M = json.load(gzip.open(MACRO))
units, ukey = [], {}
# The zoom-4 tiles stitch together six regional maps. Two are very coarse world compilations (a few hundred
# units, or 21 for Africa) with huge polygons and straight clip edges where they meet the detailed maps.
COARSE_SOURCES = {"Generalized Geology of the World", "Africa"}
GRID = np.zeros((NT * RES, NT * RES), np.uint16)


def ring_area(r):
    a = np.asarray(r)
    return abs(np.dot(a[:-1, 0], a[1:, 1]) - np.dot(a[1:, 0], a[:-1, 1])) / 2 if len(a) > 2 else 0


for key, b64 in M["tiles"].items():
    z, tx, ty = map(int, key.split("/"))
    t = mvt.decode(base64.b64decode(b64), default_options={"y_coord_down": True})
    if "units" not in t:
        continue
    polys = []
    for f in t["units"]["features"]:
        p = f["properties"]
        k = (p.get("name"), p.get("lith"), p.get("age"), p.get("best_age_top"), p.get("best_age_bottom"))
        if k not in ukey:
            ukey[k] = len(units) + 1
            nm = p.get("name") or ""
            sub = "submarine" in nm.lower() or "ocean" in nm.lower() or "sea floor" in nm.lower()
            sh = None if sub else shares_from_text(p.get("lith"), nm)
            src = p.get("ref_name") or ""
            # "detailed": from a regional map (not one of the two coarse world compilations) and with a lithology
            # that says what the rock is (less than half "sedimentary, type not stated")
            # units named only after whole tectonic plates or seas (generalised legend units of the Europe map) count as coarse
            platey = bool(re.search(r"\bplates?\b|\bsea\b|ocean", nm, re.I))
            det = int(bool(sh) and src not in COARSE_SOURCES and sh.get("sedu", 0) < 0.5 and not platey)
            units.append({"id": ukey[k], "name": clean(nm), "lith": clean(p.get("lith") or ""), "age": clean(p.get("age") or ""),
                          "submarine": int(sub), "shares": sh, "px": 0, "source": clean(src), "detailed": det})
        g = f["geometry"]
        parts = [g["coordinates"]] if g["type"] == "Polygon" else g["coordinates"] if g["type"] == "MultiPolygon" else []
        for rings in parts:
            polys.append((ukey[k], [[(x * S, y * S) for x, y in r] for r in rings]))
    polys.sort(key=lambda q: -ring_area(q[1][0]))
    im = Image.new("I", (RES, RES), 0)
    d = ImageDraw.Draw(im)
    for v, rings in polys:
        if len(rings[0]) >= 3:
            d.polygon(rings[0], fill=v)
        for h in rings[1:]:
            if len(h) >= 3:
                d.polygon(h, fill=0)
    GRID[ty * RES:(ty + 1) * RES, tx * RES:(tx + 1) * RES] = np.asarray(im, dtype=np.int32)
print("units", len(units))

# ---------------------------------------------------------------------------
# Natural Earth land at 0.1 degrees (lat/lon grid, rows from 90 N)
NX, NY, D = 3600, 1800, 0.1
land_img = Image.new("L", (NX, NY), 0)
dl = ImageDraw.Draw(land_img)
for fn, val in (("ne_10m_land.geojson", 1),):
    path = os.path.join(NE, fn)
    if not os.path.exists(path):
        print("missing", path); continue
    for f in json.load(open(path))["features"]:
        g = f["geometry"]
        parts = [g["coordinates"]] if g["type"] == "Polygon" else g["coordinates"]
        for rings in parts:
            for i, r in enumerate(rings):
                pts = [((x + 180) / D, (90 - y) / D) for x, y in r]
                if len(pts) >= 3:
                    dl.polygon(pts, fill=val if i == 0 else (0 if val else 1))
LAND = np.asarray(land_img) > 0

# ---------------------------------------------------------------------------
# sum unit pixels into the 0.1 degree grid
NU = len(units) + 1
NG = len(GKEYS)
SH = np.zeros((NU, NG), np.float32)
for u in units:
    if u["shares"]:
        for k, v in u["shares"].items():
            SH[u["id"], GKEYS.index(k)] = v
mapped_u = SH.sum(1) > 0
DET = np.zeros(NU, np.float32)
for u in units:
    DET[u["id"]] = u["detailed"]
N = 256 * 2 ** (TZ + 1)
xs = np.arange(NT * RES) + 0.5
lon = xs / N * 360 - 180
ci = np.clip((lon + 180) / D, 0, NX - 1).astype(np.int32)
acc = np.zeros((NG + 2, NY * NX), np.float64)
for y0 in range(0, NT * RES, 512):
    blk = GRID[y0:y0 + 512].astype(np.int64)
    yy = np.arange(y0, y0 + blk.shape[0]) + 0.5
    lat = np.degrees(np.arctan(np.sinh(math.pi * (1 - 2 * yy / N))))
    ri = np.clip((90 - lat) / D, 0, NY - 1).astype(np.int32)
    # weight each mercator pixel by its true area relative to the cell (cos^2 lat for web mercator pixels)
    wrow = np.cos(np.radians(lat)) ** 2
    idx = (ri[:, None] * NX + ci[None, :]).ravel()
    w = np.repeat(wrow, blk.shape[1])
    u = blk.ravel()
    land = LAND.ravel()[idx]
    ok = (u > 0) & mapped_u[u] & land
    acc[NG] += np.bincount(idx, weights=w * land, minlength=NY * NX)
    acc[NG + 1] += np.bincount(idx[ok], weights=(w * DET[u])[ok], minlength=NY * NX)
    for j in range(NG):
        acc[j] += np.bincount(idx[ok], weights=(w * SH[u, j])[ok], minlength=NY * NX)
    cnt = np.bincount(u, minlength=NU)
    for i in np.nonzero(cnt)[0]:
        if i:
            units[i - 1]["px"] += int(cnt[i])
tot = acc[:NG].sum(0)
mapped = np.where(acc[NG] > 0, tot / np.maximum(acc[NG], 1e-12), 0).reshape(NY, NX)
shares = np.where(tot > 0, acc[:NG] / np.maximum(tot, 1e-12), 0).reshape(NG, NY, NX).astype(np.float32)
del GRID
detailed = np.where(tot > 0, acc[NG + 1] / np.maximum(tot, 1e-12), 0).reshape(NY, NX)
del acc
np.savez_compressed(os.path.join(WORK, "world_geo.npz"), shares=(shares * 250).round().astype(np.uint8), detailed=(detailed * 250).round().astype(np.uint8),
                    mapped=(np.clip(mapped, 0, 1) * 250).round().astype(np.uint8), land=LAND, wsum=(tot.reshape(NY, NX)).astype(np.float32))

# unit table
with open(os.path.join(DATA, "units_world.csv"), "w", newline="", encoding="utf-8") as fh:
    w = csv.writer(fh)
    w.writerow(["id", "name", "lithology", "age", "source_map", "detailed", "submarine", "pixels_z5"] + ["share_" + k for k in GKEYS])
    for u in units:
        s = u["shares"] or {}
        w.writerow([u["id"], u["name"], u["lith"], u["age"], u["source"], u["detailed"], u["submarine"], u["px"]] + [round(s.get(k, 0), 3) for k in GKEYS])

cell_area = (111.32 * D) ** 2 * np.cos(np.radians(90 - (np.arange(NY) + 0.5) * D))[:, None] * np.ones((1, NX))
landA = (cell_area * LAND)
ant = (90 - (np.arange(NY) + 0.5) * D) < -60
landA[ant] = 0
srcs = {}
for u in units:
    srcs.setdefault(u["source"], [0, 0]); srcs[u["source"]][0] += 1; srcs[u["source"]][1] += u["px"]
R = {"sources": {k: {"units": v[0], "pixels": v[1]} for k, v in srcs.items()}, "coarse_sources": sorted(COARSE_SOURCES),
     "detailed_share_of_mapped_land": None, "units": len(units), "units_classified": int(sum(1 for u in units if u["shares"])),
     "units_submarine": int(sum(u["submarine"] for u in units)),
     "land_km2": float(landA.sum()), "mapped_share_of_land": float((landA * mapped).sum() / landA.sum()),
     "group_share_of_mapped_land": {k: float((landA * mapped * shares[j]).sum() / (landA * mapped).sum()) for j, k in enumerate(GKEYS)}}
R["detailed_share_of_mapped_land"] = float((landA * mapped * detailed).sum() / (landA * mapped).sum())
print(R)
json.dump(R, open(os.path.join(WORK, "r_world_geo.json"), "w"), indent=1)
