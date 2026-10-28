"""Step 4. The model: expected hardness from the rocks around a place and its rainfall.

A deliberately simple, transparent model fitted to the station medians:

    hardness = (sum over rock groups k of h_k * share_k) * (rain / 1000 mm)^c * water-type factor

share_k is the share of the mapped land in a square window of plus or minus R km around the place (R chosen
by cross-validation) that is rock group k. The waters mix by concentration, so h_k is the typical
hardness of river water draining only group k under 1,000 mm of rain a year. The offsets let groundwater
and lake water differ from rivers; the maps show the river baseline. Each station is weighted by one
over the number of stations in its one-degree square, so dense national networks do not dominate.
h_k is bounded between 5 and about 3,000 mg/l (no rock gives water softer than 5 mg/l in this model). The fit minimises squared errors of log10 hardness.

Accuracy is measured honestly by spatial cross-validation: stations are grouped into 5 by 5 degree
blocks and whole blocks are held out (10 folds), so a station is never predicted from its neighbours.
A harder test leaves out one country at a time. The window radius is chosen by the same cross-validation
among radii up to a cap (50 km for the world, 10 km for Britain and Ireland) so the maps keep local detail;
the scores for every radius tried are kept in results.json.

The world model uses the Macrostrat rock groups at 0.1 degrees. A second, sharper model for Britain and
Ireland uses the 1 km BGS and GSI grid with four merged groups, calibrated on the Irish stations (the
open archive has no stations in Great Britain) and then applied to the whole of both islands.

Each estimate also carries a confidence: the share of the mapped land in the same window that comes from a detailed
map unit (a regional map, not one of the two coarse world compilations, and a lithology that says what the rock is).

Writes
  data/station_model.csv   per station: window rock shares, confidence, rainfall, predicted hardness, residual
  data/model_coefficients.csv
  work/world_pred.npz, work/bi_pred.npz   predicted river-baseline hardness grids (temporary)
  work/r_model.json
Run: python3 code/04_model.py
"""
import json, math, os, sys
import numpy as np
import pandas as pd
from scipy.optimize import lsq_linear
from scipy import ndimage

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import *  # noqa

R = {}
st = pd.read_csv(os.path.join(DATA, "stations.csv"))
st = st[st["lat"].between(-90, 90) & st["lon"].between(-180, 180)].reset_index(drop=True)

# ---------------------------------------------------------------------------
# rainfall (CRU CL 2.0 annual, 0.5 degrees), filled to the coast by nearest cell, at 0.1 degrees
rn = pd.read_csv(RAIN)
RA = np.full((360, 720), np.nan)
RA[((90 - rn["lat"]) / 0.5).astype(int).clip(0, 359), ((rn["lon"] + 180) / 0.5).astype(int).clip(0, 719)] = rn["rain_mm"]
idx = ndimage.distance_transform_edt(np.isnan(RA), return_distances=False, return_indices=True)
RA = RA[idx[0], idx[1]]
RAIN01 = np.repeat(np.repeat(RA, 5, 0), 5, 1)  # 1800 x 3600


def rain_at(lat, lon):
    r = ((90 - np.asarray(lat)) / 0.5).astype(int).clip(0, 359); c = ((np.asarray(lon) + 180) / 0.5).astype(int).clip(0, 719)
    return RA[r, c]


# ---------------------------------------------------------------------------
def box_smooth(arrs, lat_rows, dlat_deg, dlon_deg, Rkm, wrap):
    """Sum of each array over a window of +/- Rkm (box in km, so wider in degrees towards the poles)."""
    ny, nx = arrs[0].shape
    out = []
    hv = max(0, int(round(Rkm / (111.32 * dlat_deg))))
    hl = np.minimum(nx // 2 - 1, np.round(Rkm / (111.32 * dlon_deg * np.maximum(np.cos(np.radians(lat_rows)), 0.02)))).astype(int)
    for a in arrs:
        b = np.empty_like(a, dtype=np.float64)
        for i in range(ny):
            h = hl[i]
            row = a[i].astype(np.float64)
            if wrap:
                p = np.concatenate([row[-h - 1:] if h + 1 <= nx else row, row, row[:h]]) if h > 0 else np.concatenate([[0.0], row])
                cs = np.cumsum(p)
                if h > 0:
                    b[i] = cs[2 * h + 1:2 * h + 1 + nx] - cs[:nx]
                else:
                    b[i] = row
            else:
                cs = np.concatenate([[0.0], np.cumsum(row)])
                lo = np.clip(np.arange(nx) - h, 0, nx); hi = np.clip(np.arange(nx) + h + 1, 0, nx)
                b[i] = cs[hi] - cs[lo]
        cs = np.concatenate([np.zeros((1, nx)), np.cumsum(b, 0)])
        lo = np.clip(np.arange(ny) - hv, 0, ny); hi = np.clip(np.arange(ny) + hv + 1, 0, ny)
        out.append(cs[hi] - cs[lo])
    return out


def fit(X, y, w, nb):
    """Weighted bounded least squares on log10 hardness.

    Mixing form: hardness = (sum_k 10^a_k * share_k) * 10^(c . z), so the waters from each rock group mix by
    concentration; 10^a_k is the typical hardness from rock group k alone, bounded between 5 and about 3,000 mg/l.
    """
    from scipy.optimize import least_squares
    sw = np.sqrt(w)
    p0 = np.r_[np.full(nb, 2.0), np.zeros(X.shape[1] - nb)]
    lb = np.r_[np.full(nb, math.log10(5)), np.full(X.shape[1] - nb, -3.0)]
    ub = np.r_[np.full(nb, 3.5), np.full(X.shape[1] - nb, 3.0)]
    r = least_squares(lambda p: (predict(X, p, nb) - y) * sw, p0, bounds=(lb, ub), method="trf", x_scale=1.0)
    return r.x


def predict(X, p, nb):
    """log10 hardness from shares (first nb columns) and other terms."""
    mix = np.maximum(X[:, :nb] @ (10 ** p[:nb]), 1e-6)
    return np.log10(mix) + X[:, nb:] @ p[nb:]


def classes(v):
    return np.digitize(v, BREAKS)


def metrics(y, yh, w=None):
    """y, yh in log10 mg/l."""
    e = yh - y
    ss = ((y - y.mean()) ** 2).sum()
    co = classes(10 ** y); cp = classes(10 ** yh)
    return {"n": int(len(y)), "r2": float(1 - (e ** 2).sum() / ss), "rmse_log10": float(np.sqrt((e ** 2).mean())),
            "median_factor": float(10 ** np.median(np.abs(e))), "within_factor_2": float((np.abs(e) <= math.log10(2)).mean()),
            "class_exact": float((co == cp).mean()), "class_within1": float((np.abs(co - cp) <= 1).mean()),
            "spearman": float(pd.Series(y).corr(pd.Series(yh), method="spearman"))}


def design(df, gcols, use_rain=True, use_type=True):
    cols = [df[c].values for c in gcols]
    names = list(gcols)
    if use_rain:
        cols.append(np.log10(np.maximum(df["rain_mm"].values, 20) / 1000)); names.append("log10_rain_per_1000mm")
    if use_type:
        cols.append((df["wtype"] == "Groundwater").astype(float).values); names.append("groundwater")
        cols.append(df["wtype"].isin(["Lake", "Reservoir"]).astype(float).values); names.append("lake_or_reservoir")
    return np.column_stack(cols), names


def lrain(df):
    return np.log10(np.maximum(df["rain_mm"].values, 20) / 1000)


def cv(df, gcols, groups, use_rain=True, use_type=True, k=10, seed=1, fixed_rain=None):
    """Spatially blocked cross-validated predictions (log10). fixed_rain: use this rain coefficient instead of fitting one."""
    if fixed_rain is not None:
        return cv(df.assign(_off=fixed_rain * lrain(df)), gcols, groups, False, use_type, k, seed, None) if False else _cv_fixed(df, gcols, groups, use_type, k, seed, fixed_rain)
    X, _ = design(df, gcols, use_rain, use_type)
    y = np.log10(df["median"].values); w = df["w"].values
    ug = np.unique(groups)
    rng = np.random.default_rng(seed); perm = rng.permutation(ug)
    fold_of = {g: i % k for i, g in enumerate(perm)}
    f = np.array([fold_of[g] for g in groups])
    yh = np.empty_like(y)
    for i in range(k):
        te = f == i
        if not te.any():
            continue
        b = fit(X[~te], y[~te], w[~te], len(gcols))
        yh[te] = predict(X[te], b, len(gcols))
    return yh


def _cv_fixed(df, gcols, groups, use_type, k, seed, c):
    off = c * lrain(df)
    X, _ = design(df, gcols, False, use_type)
    y = np.log10(df["median"].values) - off; w = df["w"].values
    ug = np.unique(groups); rng = np.random.default_rng(seed); perm = rng.permutation(ug)
    fold_of = {g: i % k for i, g in enumerate(perm)}
    f = np.array([fold_of[g] for g in groups])
    yh = np.empty_like(y)
    for i in range(k):
        te = f == i
        if te.any():
            yh[te] = predict(X[te], fit(X[~te], y[~te], w[~te], len(gcols)), len(gcols))
    return yh + off


# ===========================================================================
# WORLD
Z = np.load(os.path.join(WORK, "world_geo.npz"))
DETW = Z["detailed"].astype(np.float32) / 250
SHW = Z["shares"].astype(np.float32) / 250; MAPW = Z["mapped"].astype(np.float32) / 250; LAND = Z["land"]
NY, NX = MAPW.shape
lat_rows = 90 - (np.arange(NY) + 0.5) * 0.1
area = np.cos(np.radians(lat_rows))[:, None] * np.ones((1, NX))
wmap = MAPW * area * LAND                     # mapped land weight per cell
ri = ((90 - st["lat"].values) / 0.1).astype(int).clip(0, NY - 1); ci = ((st["lon"].values + 180) / 0.1).astype(int).clip(0, NX - 1)
st["rain_mm"] = rain_at(st["lat"], st["lon"])
cell1 = (np.floor(st["lat"]).astype(int).astype(str) + "_" + np.floor(st["lon"]).astype(int).astype(str))
st["w"] = 1.0 / cell1.map(cell1.value_counts())
blocks5 = (np.floor(st["lat"] / 5).astype(int).astype(str) + "_" + np.floor(st["lon"] / 5).astype(int).astype(str)).values
G = len(GKEYS)

RADII = [10, 25, 50, 100, 200]
cache = {}
radius_scores = []
for Rk in RADII:
    sm = box_smooth([wmap] + [SHW[j] * wmap for j in range(G)], lat_rows, 0.1, 0.1, Rk, True)
    den = sm[0]
    sh = np.stack([np.where(den > 0, s / np.maximum(den, 1e-12), np.nan) for s in sm[1:]])
    cache[Rk] = (sh, den)
    d = st.copy()
    for j, k in enumerate(GKEYS):
        d["s_" + k] = sh[j, ri, ci]
    ok = d["s_carb"].notna().values
    d = d[ok].reset_index(drop=True)
    yh = cv(d, ["s_" + k for k in GKEYS], blocks5[ok])
    m = metrics(np.log10(d["median"].values), yh)
    radius_scores.append({"radius_km": Rk, **m})
    print("R", Rk, {k: round(v, 3) for k, v in m.items()})
RMAX_W = 50   # the map should keep detail at about the scale of a county: choose among windows up to +/- 50 km
best = max([r for r in radius_scores if r["radius_km"] <= RMAX_W], key=lambda r: r["r2"])["radius_km"]
R["world_radius_cap_km"] = RMAX_W
R["world_radius_scores"] = radius_scores; R["world_radius_km"] = best
sh, den = cache[best]
# confidence: the share of the mapped land in the same window whose map unit is "detailed" (a regional map, and a
# lithology that says what the rock is). Below CONF_MIN the estimate rests mainly on coarse world compilations.
CONF_MIN = 0.5
_cs = box_smooth([wmap, DETW * wmap], lat_rows, 0.1, 0.1, best, True)
conf = np.where(_cs[0] > 0, _cs[1] / np.maximum(_cs[0], 1e-12), np.nan)
del _cs
R["conf_min"] = CONF_MIN
for Rk in list(cache):
    if Rk != best:
        del cache[Rk]
for j, k in enumerate(GKEYS):
    st["s_" + k] = sh[j, ri, ci]
ok = st["s_carb"].notna().values
R["stations_without_geology"] = int((~ok).sum())
d = st[ok].reset_index(drop=True); bl = blocks5[ok]
y = np.log10(d["median"].values)
gcols = ["s_" + k for k in GKEYS]

# comparisons under the same spatial cross-validation
variants = {}
variants["baseline_global_median"] = metrics(y, np.full_like(y, np.average(y, weights=d["w"])))
variants["rain_and_type_only"] = metrics(y, cv(d.assign(one=1.0), ["one"], bl))
variants["rocks_only"] = metrics(y, cv(d, gcols, bl, use_rain=False, use_type=False))
variants["rocks_and_type"] = metrics(y, cv(d, gcols, bl, use_rain=False))
yh_full = cv(d, gcols, bl)
variants["rocks_rain_and_type"] = metrics(y, yh_full)
cc = d["country"].values
yh_loco = np.empty_like(y)
X, names = design(d, gcols)
for c in np.unique(cc):
    te = cc == c
    b = fit(X[~te], y[~te], d["w"].values[~te], G)
    yh_loco[te] = predict(X[te], b, G)
variants["rocks_rain_and_type_leave_country_out"] = metrics(y, yh_loco)
R["world_cv"] = variants
for k, v in variants.items():
    print(k, {a: round(b, 3) for a, b in v.items()})

# by water type and by country, spatial CV predictions
d["pred_cv"] = 10 ** yh_full
d["conf"] = conf[ri[ok], ci[ok]]
hc = d["conf"].values >= CONF_MIN
R["world_cv_by_confidence"] = {"high": metrics(y[hc], yh_full[hc]), "low": metrics(y[~hc], yh_full[~hc])}
R["world_cv_by_type"] = {t: metrics(y[d["wtype"].values == t], yh_full[d["wtype"].values == t]) for t in ["River", "Groundwater", "Lake", "Reservoir"]}
byc = []
for c, g in d.groupby("country"):
    if len(g) >= 10:
        byc.append({"country": c, "n": int(len(g)), "measured_median": float(g["median"].median()), "cv_pred_median": float(g["pred_cv"].median())})
R["world_cv_country_medians"] = byc
lm = np.log10([r["measured_median"] for r in byc]); lp = np.log10([r["cv_pred_median"] for r in byc])
R["world_cv_country_r2"] = float(1 - ((lm - lp) ** 2).sum() / ((lm - lm.mean()) ** 2).sum())
R["world_cv_country_spearman"] = float(pd.Series(lm).corr(pd.Series(lp), method="spearman"))

# final world fit
bw = fit(X, y, d["w"].values, G)
R["world_coef"] = {n: float(v) for n, v in zip(names, bw)}
R["world_typical"] = {k: float(10 ** bw[j]) for j, k in enumerate(GKEYS)}
CW = float(bw[G])   # the world rain effect, reused (fixed) in the Britain and Ireland model
print("coef", {n: round(v, 3) for n, v in zip(names, bw)})
d["pred"] = 10 ** predict(X, bw, G)
d["resid_log10"] = y - predict(X, bw, G)
d["pred_river"] = 10 ** predict(X[:, :G + 1], bw[:G + 1], G)

# world prediction grid (river baseline)
lr = np.log10(np.maximum(RAIN01, 20) / 1000)
pw = np.log10(np.maximum(np.einsum("kij,k->ij", np.nan_to_num(sh), 10 ** bw[:G]), 1e-6)) + bw[G] * lr
DRY = 200   # mm a year: drier land has little surface water to measure and no stations to learn from; left blank
valid = (den > 0) & LAND & (lat_rows[:, None] > -60) & (MAPW > 0) & (RAIN01 >= DRY)
R["dry_mask_mm"] = DRY
R["station_rain_p01"] = float(np.percentile(d["rain_mm"], 1)); R["stations_below_dry"] = int((d["rain_mm"] < DRY).sum())
R["land_share_too_dry"] = float(np.where(LAND & (lat_rows[:, None] > -60) & (RAIN01 < DRY), area, 0).sum() / np.where(LAND & (lat_rows[:, None] > -60), area, 0).sum())
predw = np.where(valid, 10 ** pw, np.nan).astype(np.float32)
np.savez_compressed(os.path.join(WORK, "world_pred.npz"), pred=predw.astype(np.float16), carb=np.where(valid, sh[0] + sh[1], np.nan).astype(np.float16),
                    conf=np.where(valid, conf, np.nan).astype(np.float16))
A = area * 1.0
cls = np.digitize(np.nan_to_num(predw, nan=-1), BREAKS)
wv = np.where(valid, A, 0)
R["world_pred_class_share"] = [float(wv[(cls == i)].sum() / wv.sum()) for i in range(6)]
R["world_pred_low_conf_share"] = float(wv[np.nan_to_num(conf, nan=0) < CONF_MIN].sum() / wv.sum())
hv = np.where(valid & (np.nan_to_num(conf, nan=0) >= CONF_MIN), A, 0)
R["world_pred_class_share_high_conf"] = [float(hv[(cls == i)].sum() / hv.sum()) for i in range(6)]
R["world_pred_land_share_with_prediction"] = float(wv.sum() / np.where(LAND & (lat_rows[:, None] > -60), A, 0).sum())
del cache

# ===========================================================================
# BRITAIN AND IRELAND (1 km)
ZB = np.load(os.path.join(WORK, "bi_geo.npz"))
SB = ZB["shares"].astype(np.float32) / 250; MB = ZB["mapped"].astype(np.float32) / 250; SRC = ZB["src"]
LON0, LON1, LAT0, LAT1, DY, DX = ZB["box"]
nyb, nxb = MB.shape
latb = LAT1 - (np.arange(nyb) + 0.5) * DY
lonb = LON0 + (np.arange(nxb) + 0.5) * DX
BGRP = [("carb", ["carb", "evap"]), ("clas", ["clas", "unc", "sedu"]), ("volc", ["volc"]), ("crys", ["crys"])]
SB4 = np.stack([sum(SB[GKEYS.index(k)] for k in ks) for _, ks in BGRP])
wb = MB * np.cos(np.radians(latb))[:, None]
inb = st["lat"].between(LAT0, LAT1) & st["lon"].between(LON0, LON1)
rb = ((LAT1 - st["lat"]) / DY).astype(int).clip(0, nyb - 1); cb = ((st["lon"] - LON0) / DX).astype(int).clip(0, nxb - 1)
st["bi_src"] = np.where(inb, SRC[rb, cb], 0)
ie = st[(st["country"] == "Ireland") & inb].copy()
R["ireland_stations"] = int(len(ie))
gb4 = ["b_" + k for k, _ in BGRP]
bi_scores = []; bcache = {}
blocks_ie_all = (np.floor(ie["lat"] / 0.5).astype(int).astype(str) + "_" + np.floor(ie["lon"] / 0.5).astype(int).astype(str)).values
for Rk in [1, 2, 5, 10, 20]:
    smb = box_smooth([wb] + [SB4[j] * wb for j in range(4)], latb, DY, DX, Rk, False)
    shb = np.stack([np.where(smb[0] > 0, s / np.maximum(smb[0], 1e-12), np.nan) for s in smb[1:]])
    bcache[Rk] = (shb, smb[0])
    e = ie.copy()
    for j, (k, _) in enumerate(BGRP):
        e["b_" + k] = shb[j, rb[e.index], cb[e.index]]
    okk = e["b_carb"].notna().values
    e = e[okk]; e["w"] = 1.0
    yh = cv(e, gb4, blocks_ie_all[okk], fixed_rain=CW)
    m = metrics(np.log10(e["median"].values), yh)
    bi_scores.append({"radius_km": Rk, **m})
    print("BI R", Rk, {k: round(v, 3) for k, v in m.items()})
RMAX_B = 10   # likewise for the 1 km map: windows up to +/- 10 km
bbest = max([r for r in bi_scores if r["radius_km"] <= RMAX_B], key=lambda r: r["r2"])["radius_km"]
R["bi_radius_cap_km"] = RMAX_B
R["bi_radius_scores"] = bi_scores; R["bi_radius_km"] = bbest
shb, denb = bcache[bbest]; del bcache
for j, (k, _) in enumerate(BGRP):
    ie["b_" + k] = shb[j, rb[ie.index], cb[ie.index]]
ie = ie[ie["b_carb"].notna()].copy(); ie["w"] = 1.0
bl_ie = (np.floor(ie["lat"] / 0.5).astype(int).astype(str) + "_" + np.floor(ie["lon"] / 0.5).astype(int).astype(str)).values
yi = np.log10(ie["median"].values)
yh_ie = cv(ie, gb4, bl_ie, fixed_rain=CW)
R["bi_cv"] = {"rocks_type_and_world_rain_effect": metrics(yi, yh_ie), "baseline_median": metrics(yi, np.full_like(yi, np.median(yi))),
              "rocks_only": metrics(yi, cv(ie, gb4, bl_ie, use_rain=False, use_type=False)),
              "rocks_and_type": metrics(yi, cv(ie, gb4, bl_ie, use_rain=False)),
              "rocks_type_and_free_rain_effect": metrics(yi, cv(ie, gb4, bl_ie))}
Xf, nf = design(ie, gb4)
R["bi_free_rain_coef"] = float(fit(Xf, yi, ie["w"].values, 4)[4])
# the world model at the same Irish stations, for comparison
wm = d.set_index("id").reindex(ie["id"])
okw = wm["pred_cv"].notna().values
R["bi_cv"]["world_model_at_irish_stations"] = metrics(yi[okw], np.log10(wm["pred_cv"].values[okw]))
for k, v in R["bi_cv"].items():
    print("BI", k, {a: round(b, 3) for a, b in v.items()})
XB0, nb0 = design(ie, gb4, use_rain=False)
b0 = fit(XB0, yi - CW * lrain(ie), ie["w"].values, 4)
bb = np.r_[b0[:4], CW, b0[4:]]
XB, nb_names = design(ie, gb4)
R["bi_coef"] = {n: float(v) for n, v in zip(nb_names, bb)}
R["bi_typical"] = {k: float(10 ** bb[j]) for j, (k, _) in enumerate(BGRP)}
print("BI coef", R["bi_coef"])
ie["pred_bi"] = 10 ** predict(XB, bb, 4); ie["resid_bi_log10"] = yi - predict(XB, bb, 4)
RAINB = rain_at(latb[:, None] * np.ones((1, nxb)), lonb[None, :] * np.ones((nyb, 1)))
pb = np.log10(np.maximum(np.einsum("kij,k->ij", np.nan_to_num(shb), 10 ** bb[:4]), 1e-6)) + bb[4] * np.log10(np.maximum(RAINB, 20) / 1000)
validb = (denb > 0) & (MB > 0)
predb = np.where(validb, 10 ** pb, np.nan).astype(np.float32)
np.savez_compressed(os.path.join(WORK, "bi_pred.npz"), pred=predb.astype(np.float16), carb=np.where(validb, shb[0], np.nan).astype(np.float16),
                    box=ZB["box"])

# ---------------------------------------------------------------------------
# outputs
out = st[["id"]].copy()
out = out.merge(d[["id"] + gcols + ["rain_mm", "conf", "pred", "pred_river", "pred_cv", "resid_log10"]], on="id", how="left")
out["conf"] = out["conf"].round(3)
out = out.merge(ie[["id"] + gb4 + ["pred_bi", "resid_bi_log10"]], on="id", how="left")
for c in gcols + gb4:
    out[c] = out[c].round(3)
for c in ["pred", "pred_river", "pred_cv", "pred_bi"]:
    out[c] = out[c].round(1)
out["rain_mm"] = out["rain_mm"].round(0)
out["resid_log10"] = out["resid_log10"].round(3); out["resid_bi_log10"] = out["resid_bi_log10"].round(3)
out.to_csv(os.path.join(DATA, "station_model.csv"), index=False)
rows = [["world", n, round(v, 4), round(10 ** v, 1) if j < G else ""] for j, (n, v) in enumerate(zip(names, bw))]
rows += [["britain_ireland", n, round(v, 4), round(10 ** v, 1) if j < 4 else ""] for j, (n, v) in enumerate(zip(nb_names, bb))]
pd.DataFrame(rows, columns=["model", "term", "coefficient_log10", "typical_mg_l_if_all_this_rock"]).to_csv(os.path.join(DATA, "model_coefficients.csv"), index=False)
R["bi_groups"] = {k: ks for k, ks in BGRP}
R["n_model_stations"] = int(len(d))
json.dump(R, open(os.path.join(WORK, "r_model.json"), "w"), indent=1)
