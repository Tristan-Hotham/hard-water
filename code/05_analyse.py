"""Step 5. All the statistics used anywhere, written to data/results.json, plus the tables.

Writes
  data/results.json            every number used in the page, paper, blog and collage
  data/countries_measured.csv  country league table from the stations
  data/countries_modelled.csv  modelled river-baseline hardness by country (Natural Earth admin 0)
  data/bi_areas.csv            Britain and Ireland: modelled class shares by nation and county
  data/places.csv              GeoNames places with their modelled hardness (world and, where covered, the 1 km estimate)
Run: python3 code/05_analyse.py
"""
import csv, io, json, math, os, sys, zipfile
import numpy as np
import pandas as pd
from PIL import Image, ImageDraw

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import *  # noqa

R = {}
RS = json.load(open(os.path.join(WORK, "r_stations.json")))
RW = json.load(open(os.path.join(WORK, "r_world_geo.json")))
RB = json.load(open(os.path.join(WORK, "r_bi_geo.json")))
RM = json.load(open(os.path.join(WORK, "r_model.json")))
R["stations"] = {k: v for k, v in RS.items() if k not in ("country_medians",)}
R["world_geo"] = RW; R["bi_geo"] = RB; R["model"] = RM
R["classes"] = [{"key": c[0], "label": c[1], "range": r, "colour": c[4]} for c, r in zip(CLASSES, CLASS_RANGES)]
R["units"] = UNITS

st = pd.read_csv(os.path.join(DATA, "stations.csv"))
sm = pd.read_csv(os.path.join(DATA, "station_model.csv"))
d = st.merge(sm, on="id", how="left")
d["cls"] = d["median"].apply(hclass)

# ---------------------------------------------------------------------------
# overall
R["overall"] = {"median": float(d["median"].median()), "n": int(len(d)),
                "share_by_class": [float((d["cls"] == i).mean()) for i in range(6)],
                "by_type": {t: {"n": int(len(g)), "median": float(g["median"].median())} for t, g in d.groupby("wtype")},
                "samples_median_per_station": float(d["n"].median()), "stations_10plus_samples": int((d["n"] >= 10).sum()),
                "first_year_median": int(pd.to_datetime(d["first"]).dt.year.median()), "last_year_median": int(pd.to_datetime(d["last"]).dt.year.median()),
                "stations_with_data_since_2010": int((pd.to_datetime(d["last"]).dt.year >= 2010).sum()),
                "max_station": d.loc[d["median"].idxmax(), ["sname", "country", "median", "wtype"]].to_dict(),
                "min_station": d.loc[d["median"].idxmin(), ["sname", "country", "median", "wtype"]].to_dict()}

# country league (measured)
rows = []
for c, g in d.groupby("country"):
    rows.append({"country": c, "stations": int(len(g)), "samples": int(g["n"].sum()), "median": round(float(g["median"].median()), 1),
                 "river_median": round(float(g[g["wtype"] == "River"]["median"].median()), 1) if (g["wtype"] == "River").any() else None,
                 "groundwater_median": round(float(g[g["wtype"] == "Groundwater"]["median"].median()), 1) if (g["wtype"] == "Groundwater").any() else None,
                 "share_soft": round(float((g["cls"] <= 1).mean()), 3), "share_hard": round(float((g["cls"] >= 4).mean()), 3),
                 "first": str(g["first"].min())[:4], "last": str(g["last"].max())[:4]})
cm = pd.DataFrame(rows).sort_values("median", ascending=False)
cm.to_csv(os.path.join(DATA, "countries_measured.csv"), index=False)
R["countries_measured"] = cm.to_dict("records")
big = cm[cm["stations"] >= 20]
R["countries_measured_n20"] = int(len(big))
R["softest_country"] = big.iloc[-1].to_dict(); R["hardest_country"] = big.iloc[0].to_dict()

# class shares by type of water
R["type_shares"] = []
for lab, ts in [("Rivers", ["River"]), ("Lakes and reservoirs", ["Lake", "Reservoir"]), ("Groundwater", ["Groundwater"]), ("All stations", None)]:
    g = d if ts is None else d[d["wtype"].isin(ts)]
    R["type_shares"].append({"label": lab, "n": int(len(g)), "shares": [round(float((g["cls"] == i).mean()), 4) for i in range(6)]})

# carbonate share against measured hardness (world windows)
bins = [0, 0.05, 0.15, 0.3, 0.5, 1.01]
carb = d["s_carb"] + d["s_evap"]
R["by_carbonate"] = []
for lo, hi in zip(bins[:-1], bins[1:]):
    g = d[(carb >= lo) & (carb < hi)]
    R["by_carbonate"].append({"lo": lo, "hi": min(hi, 1), "n": int(len(g)), "median": float(g["median"].median()),
                              "river": float(g[g["wtype"] == "River"]["median"].median()), "q25": float(g["median"].quantile(.25)), "q75": float(g["median"].quantile(.75))})
crys = d["s_crys"]
R["by_crystalline"] = []
for lo, hi in zip([0, 0.2, 0.5, 0.8], [0.2, 0.5, 0.8, 1.01]):
    g = d[(crys >= lo) & (crys < hi)]
    R["by_crystalline"].append({"lo": lo, "hi": min(hi, 1), "n": int(len(g)), "median": float(g["median"].median()),
                                "q25": float(g["median"].quantile(.25)), "q75": float(g["median"].quantile(.75))})
R["corr_carb_spearman"] = float(pd.Series(carb).corr(np.log10(d["median"]), method="spearman"))
R["corr_crys_spearman"] = float(pd.Series(crys).corr(np.log10(d["median"]), method="spearman"))
R["corr_rain_spearman"] = float(d["rain_mm"].corr(np.log10(d["median"]), method="spearman"))

# residuals: where the model is most wrong (by country, spatial CV)
d["cv_resid"] = np.log10(d["median"]) - np.log10(d["pred_cv"])
res = []
for c, g in d[d["pred_cv"].notna()].groupby("country"):
    if len(g) >= 20:
        res.append({"country": c, "n": int(len(g)), "median_resid_factor": float(10 ** g["cv_resid"].median())})
R["resid_by_country"] = sorted(res, key=lambda r: r["median_resid_factor"])


# ---------------------------------------------------------------------------
# stories
def polys(g):
    """Polygon ring lists of a GeoJSON geometry (Polygon, MultiPolygon or GeometryCollection)."""
    if g is None:
        return []
    if g["type"] == "Polygon":
        return [g["coordinates"]]
    if g["type"] == "MultiPolygon":
        return g["coordinates"]
    if g["type"] == "GeometryCollection":
        return [p for gg in g["geometries"] for p in polys(gg)]
    return []


def box(df, lat0, lat1, lon0, lon1):
    return df[df["lat"].between(lat0, lat1) & df["lon"].between(lon0, lon1)]


S = {}
yu = box(d, 17.8, 21.7, -92.0, -86.5)
S["yucatan"] = {"box": [17.8, 21.7, -92.0, -86.5], "n": int(len(yu)), "median": float(yu["median"].median()),
                "groundwater_n": int((yu["wtype"] == "Groundwater").sum()), "share_vhard": float((yu["cls"] == 5).mean()),
                "share_hard_or_vhard": float((yu["cls"] >= 4).mean()),
                "pred_cv_median": float(yu["pred_cv"].median()), "carb_share_median": float((yu["s_carb"] + yu["s_evap"]).median()),
                "sedu_share_median": float(yu["s_sedu"].median())}
mx = d[d["country"] == "Mexico"]
S["mexico"] = {"n": int(len(mx)), "median": float(mx["median"].median()), "groundwater_n": int((mx["wtype"] == "Groundwater").sum())}
nf = box(d[d["country"] == "France"], 48.6, 50.6, 0.0, 3.6)
S["northern_france_chalk"] = {"box": [48.6, 50.6, 0.0, 3.6], "n": int(len(nf)), "median": float(nf["median"].median()),
                              "river_median": float(nf[nf["wtype"] == "River"]["median"].median()),
                              "share_hard_or_vhard": float((nf["cls"] >= 4).mean())}
sc = d[d["country"].isin(["Norway", "Sweden", "Finland"])]
S["scandinavia"] = {"n": int(len(sc)), "median": float(sc["median"].median()), "share_soft": float((sc["cls"] == 0).mean()),
                    "by_country": {c: float(g["median"].median()) for c, g in sc.groupby("country")}}
ca = d[d["country"] == "Canada"]
shield = ca[ca["s_crys"] >= 0.6]; rest = ca[ca["s_crys"] < 0.6]
S["canada"] = {"n": int(len(ca)), "median": float(ca["median"].median()), "shield_n": int(len(shield)), "shield_median": float(shield["median"].median()),
               "other_n": int(len(rest)), "other_median": float(rest["median"].median())}
ie = d[d["country"] == "Ireland"]
S["ireland"] = {"n": int(len(ie)), "median": float(ie["median"].median()),
                "by_type": {t: {"n": int(len(g)), "median": float(g["median"].median())} for t, g in ie.groupby("wtype")},
                "high_carb_n": int((ie["b_carb"] >= 0.5).sum()), "high_carb_median": float(ie[ie["b_carb"] >= 0.5]["median"].median()),
                "low_carb_n": int((ie["b_carb"] < 0.1).sum()), "low_carb_median": float(ie[ie["b_carb"] < 0.1]["median"].median())}
it = d[d["country"] == "Italy"]
S["italy"] = {"n": int(len(it)), "median": float(it["median"].median())}
inn = d[d["country"] == "India"]
S["india"] = {"n": int(len(inn)), "median": float(inn["median"].median()), "groundwater_n": int((inn["wtype"] == "Groundwater").sum()),
              "groundwater_median": float(inn[inn["wtype"] == "Groundwater"]["median"].median())}

# modelled values for south China karst (no stations in the open archive)
ZP = np.load(os.path.join(WORK, "world_pred.npz")); PW = ZP["pred"].astype(np.float32); CW = ZP["carb"].astype(np.float32); CONF = ZP["conf"].astype(np.float32)
CONF_MIN = RM["conf_min"]
ZG = np.load(os.path.join(WORK, "world_geo.npz")); SHW = ZG["shares"].astype(np.float32) / 250
NY, NX = PW.shape
lat_rows = 90 - (np.arange(NY) + 0.5) * 0.1
AREA = np.cos(np.radians(lat_rows))[:, None] * np.ones((1, NX))


def cellbox(lat0, lat1, lon0, lon1):
    r0, r1 = int((90 - lat1) / 0.1), int((90 - lat0) / 0.1); c0, c1 = int((lon0 + 180) / 0.1), int((lon1 + 180) / 0.1)
    return slice(r0, r1), slice(c0, c1)


rs, cs = cellbox(22, 27, 104, 112)
pv = PW[rs, cs]; ok = np.isfinite(pv)
S["south_china"] = {"box": [22, 27, 104, 112], "stations": int(len(box(d, 22, 27, 104, 112))), "pred_median": float(np.median(pv[ok])),
                    "low_conf_share": float((np.nan_to_num(CONF[rs, cs][ok], nan=0) < CONF_MIN).mean()),
                    "carb_share_mapped": float(np.nanmean(SHW[0][rs, cs][ok] + SHW[1][rs, cs][ok])),
                    "sedu_share_mapped": float(np.nanmean(SHW[GKEYS.index("sedu")][rs, cs][ok]))}
rs, cs = cellbox(17.8, 21.7, -92.0, -86.5)
pv = PW[rs, cs]; ok = np.isfinite(pv)
S["yucatan"]["model_map_median"] = float(np.median(pv[ok]))
S["yucatan"]["low_conf_share"] = float((np.nan_to_num(CONF[rs, cs][ok], nan=0) < CONF_MIN).mean())
R["stories"] = S

# ---------------------------------------------------------------------------
# modelled countries (Natural Earth admin 0 at 0.1 degrees)
ne = json.load(open(os.path.join(NE, "ne_10m_admin_0_countries.geojson")))
cid = Image.new("I", (NX, NY), 0); dr = ImageDraw.Draw(cid)
names = [None]
for f in ne["features"]:
    p = f["properties"]; names.append((p.get("NAME") or p.get("ADMIN"), p.get("ADM0_A3")))
    for rings in polys(f["geometry"]):
        for i, r in enumerate(rings):
            pts = [((x + 180) / 0.1, (90 - y) / 0.1) for x, y in r]
            if len(pts) >= 3:
                dr.polygon(pts, fill=len(names) - 1 if i == 0 else 0)
CID = np.asarray(cid, dtype=np.int32)
mrows = []
for k in range(1, len(names)):
    m = CID == k
    tot = AREA[m].sum()
    if tot == 0:
        continue
    km2 = float(tot * (111.32 * 0.1) ** 2)
    okm = m & np.isfinite(PW)
    cov = float(AREA[okm].sum() / tot)
    if km2 < 20000 or cov < 0.5:
        continue
    v = PW[okm]; w = AREA[okm]
    o = np.argsort(v); cw = np.cumsum(w[o]) / w.sum()
    med = float(v[o][np.searchsorted(cw, 0.5)])
    cl = np.digitize(v, BREAKS)
    lowc = float(AREA[okm & (np.nan_to_num(CONF, nan=0) < CONF_MIN)].sum() / AREA[okm].sum())
    mrows.append({"country": names[k][0], "a3": names[k][1], "land_km2": round(km2), "coverage": round(cov, 3), "modelled_median": round(med, 1), "share_low_confidence": round(lowc, 3),
                  "share_soft": round(float(w[cl <= 1].sum() / w.sum()), 3), "share_hard": round(float(w[cl >= 4].sum() / w.sum()), 3)})
mm = pd.DataFrame(mrows).sort_values("modelled_median", ascending=False)
mm.to_csv(os.path.join(DATA, "countries_modelled.csv"), index=False)
R["countries_modelled_n"] = int(len(mm))
R["countries_modelled"] = mm.to_dict("records")
# measured against modelled for countries with both
mmm = cm.merge(mm, on="country", how="inner")
mmm = mmm[mmm["stations"] >= 20]
R["measured_vs_modelled_countries"] = {"n": int(len(mmm)), "spearman": float(mmm["median"].corr(mmm["modelled_median"], method="spearman"))}

# ---------------------------------------------------------------------------
# Britain and Ireland at 1 km
ZB = np.load(os.path.join(WORK, "bi_pred.npz")); PB = ZB["pred"].astype(np.float32); BOX = ZB["box"]
LON0, LON1, LAT0, LAT1, DY, DX = BOX
nyb, nxb = PB.shape
latb = LAT1 - (np.arange(nyb) + 0.5) * DY
AB = np.cos(np.radians(latb))[:, None] * np.ones((1, nxb)) * (111.32 * DY) * (111.32 * DX)
cty = json.load(open("/home/claude/hist/how-roman-is-your-road/data/counties.geojson"))
im = Image.new("I", (nxb, nyb), 0); dr = ImageDraw.Draw(im)
cnames = [None]
for f in cty["features"]:
    p = f["properties"]; cnames.append((p["county"], p["nation"]))
    for rings in polys(f["geometry"]):
        for i, r in enumerate(rings):
            pts = [((x - LON0) / DX, (LAT1 - y) / DY) for x, y in r]
            if len(pts) >= 3:
                dr.polygon(pts, fill=len(cnames) - 1 if i == 0 else 0)
CB = np.asarray(im, dtype=np.int32)
# island of Ireland from Natural Earth (Ireland, and the UK part west of 5.4 W and south of 55.4 N for Northern Ireland)
imi = Image.new("L", (nxb, nyb), 0); dri = ImageDraw.Draw(imi)
for f in ne["features"]:
    if f["properties"].get("ADM0_A3") not in ("IRL", "GBR"):
        continue
    a3 = f["properties"]["ADM0_A3"]
    for rings in polys(f["geometry"]):
        r = rings[0]
        xs = [x for x, _ in r]; ys = [y for _, y in r]
        if a3 == "GBR" and not (max(xs) < -5.3 and min(ys) > 53.9 and max(ys) < 55.5):
            continue
        dri.polygon([((x - LON0) / DX, (LAT1 - y) / DY) for x, y in r], fill=1 if a3 == "IRL" else 2)
IRL = np.asarray(imi)


def summarise(mask):
    v = PB[mask]; w = AB[mask]; ok = np.isfinite(v); v, w = v[ok], w[ok]
    if w.sum() == 0:
        return None
    o = np.argsort(v); cwv = np.cumsum(w[o]) / w.sum()
    cl = np.digitize(v, BREAKS)
    return {"km2": round(float(w.sum())), "median": round(float(v[o][np.searchsorted(cwv, .5)]), 1),
            "shares": [round(float(w[cl == i].sum() / w.sum()), 4) for i in range(6)]}


areas = []
for nat in ["England", "Wales", "Scotland"]:
    ks = [i for i, c in enumerate(cnames) if c and c[1] == nat]
    s_ = summarise(np.isin(CB, ks)); s_.update(name=nat, kind="nation"); areas.append(s_)
for code, nm in [(2, "Northern Ireland"), (1, "Republic of Ireland")]:
    s_ = summarise(IRL == code); s_.update(name=nm, kind="nation"); areas.append(s_)
for i, c in enumerate(cnames):
    if c:
        s_ = summarise(CB == i)
        if s_:
            s_.update(name=c[0], kind=c[1]); areas.append(s_)
pd.DataFrame([{"name": a["name"], "kind": a["kind"], "km2": a["km2"], "modelled_median": a["median"], **{"share_" + CLASSES[i][0]: a["shares"][i] for i in range(6)}} for a in areas]).to_csv(os.path.join(DATA, "bi_areas.csv"), index=False)
R["bi_nations"] = [a for a in areas if a["kind"] == "nation"]
cts = sorted([a for a in areas if a["kind"] != "nation"], key=lambda a: -a["median"])
R["bi_counties_n"] = len(cts)
R["bi_counties_hardest"] = cts[:8]; R["bi_counties_softest"] = cts[-8:][::-1]
R["bi_counties"] = [{"name": a["name"], "nation": a["kind"], "median": a["median"], "hard_share": round(a["shares"][4] + a["shares"][5], 3)} for a in cts]
gbm = np.isin(CB, [i for i, c in enumerate(cnames) if c])
R["bi_gb"] = summarise(gbm); R["bi_ireland_island"] = summarise(IRL > 0)
# rainfall over Great Britain's 1 km cells (CRU CL 2.0, half-degree), for the note on the rainfall effect
import pandas as _pd
_rn = _pd.read_csv(RAIN)
_RA = np.full((360, 720), np.nan)
_RA[((90 - _rn["lat"]) / 0.5).astype(int).clip(0, 359), ((_rn["lon"] + 180) / 0.5).astype(int).clip(0, 719)] = _rn["rain_mm"]
_lonb = LON0 + (np.arange(nxb) + 0.5) * DX
_rr = ((90 - latb) / 0.5).astype(int).clip(0, 359)[:, None] * np.ones((1, nxb), int); _cc = np.ones((nyb, 1), int) * ((_lonb + 180) / 0.5).astype(int).clip(0, 719)[None, :]
_v = _RA[_rr, _cc][gbm & np.isfinite(PB)]
_v = _v[np.isfinite(_v)]
R["bi_gb_rain"] = {"p05": float(np.percentile(_v, 5)), "p50": float(np.percentile(_v, 50)), "p95": float(np.percentile(_v, 95))}
# south and east against north and west of a line from the Severn to the Tees (Tees-Exe line, roughly)
lonb = LON0 + (np.arange(nxb) + 0.5) * DX
LON2, LAT2 = np.meshgrid(lonb, latb)
# Tees-Exe line: from about (-3.5, 50.6) to (-1.2, 54.6); east of it is "lowland" Britain
side = (LON2 - (-3.5)) * (54.6 - 50.6) - (LAT2 - 50.6) * (-1.2 - (-3.5))
R["bi_tees_exe"] = {"line": [[-3.5, 50.6], [-1.2, 54.6]], "south_east": summarise(gbm & (side > 0)), "north_west": summarise(gbm & (side <= 0))}

# ---------------------------------------------------------------------------
# places (GeoNames cities15000) with modelled values
with zipfile.ZipFile(GEONAMES) as z:
    txt = z.read("cities15000.txt").decode("utf-8")
gn = pd.read_csv(io.StringIO(txt), sep="\t", header=None, quoting=3, low_memory=False, keep_default_na=False,
                 usecols=[1, 4, 5, 8, 14], names=["name", "lat", "lon", "cc", "pop"])
gn = gn[gn["pop"] >= 15000].reset_index(drop=True)
r_ = ((90 - gn["lat"]) / 0.1).astype(int).clip(0, NY - 1); c_ = ((gn["lon"] + 180) / 0.1).astype(int).clip(0, NX - 1)
gn["world"] = PW[r_, c_]
gn["world_low_conf"] = (np.nan_to_num(CONF[r_, c_], nan=0) < CONF_MIN).astype(int)
# nearest valid cell within 3 cells for coastal places
miss = gn["world"].isna()
for k in np.nonzero(miss.values)[0]:
    rr, cc = r_[k], c_[k]
    win = PW[max(0, rr - 3):rr + 4, max(0, cc - 3):cc + 4]
    if np.isfinite(win).any():
        gn.loc[k, "world"] = float(np.nanmedian(win))
inbi = gn["lat"].between(LAT0, LAT1) & gn["lon"].between(LON0, LON1)
rb_ = ((LAT1 - gn["lat"]) / DY).astype(int).clip(0, nyb - 1); cb_ = ((gn["lon"] - LON0) / DX).astype(int).clip(0, nxb - 1)
gn["bi"] = np.where(inbi, PB[rb_, cb_], np.nan)
for k in np.nonzero((inbi & gn["bi"].isna()).values)[0]:
    rr, cc = rb_[k], cb_[k]
    win = PB[max(0, rr - 4):rr + 5, max(0, cc - 4):cc + 5]
    if np.isfinite(win).any():
        gn.loc[k, "bi"] = float(np.nanmedian(win))
gn["world"] = gn["world"].round(0); gn["bi"] = gn["bi"].round(0)
gn["name"] = gn["name"].map(clean)
gn.to_csv(os.path.join(DATA, "places.csv"), index=False)
R["places_n"] = int(len(gn)); R["places_with_value"] = int(gn["world"].notna().sum())
TOWNS = ["London", "Birmingham", "Manchester", "Leeds", "Cambridge", "Brighton", "Norwich", "Bristol", "Cardiff", "Swansea", "Glasgow", "Edinburgh",
         "Aberdeen", "Inverness", "Belfast", "Dublin", "Galway", "Cork", "Exeter", "Plymouth", "Hull", "Newcastle upon Tyne", "Liverpool", "Lancaster", "Oxford"]
R["bi_towns"] = []
for t in TOWNS:
    g = gn[(gn["name"] == t) & gn["cc"].isin(["GB", "IE"])].sort_values("pop", ascending=False)
    if len(g):
        r = g.iloc[0]
        R["bi_towns"].append({"name": t, "bi": None if pd.isna(r["bi"]) else float(r["bi"]), "world": None if pd.isna(r["world"]) else float(r["world"])})
json.dump(R, open(os.path.join(DATA, "results.json"), "w"), indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o))
print(json.dumps({k: R[k] for k in ["overall", "stories", "bi_nations", "bi_tees_exe", "bi_towns", "measured_vs_modelled_countries"]}, indent=0, default=str)[:6000])
