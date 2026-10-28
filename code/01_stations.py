"""Step 1. Station hardness from the GEMS/Water archive (GFQA v3 open subset).

For each sample (station and date) the hardness in mg/l as CaCO3 is the total hardness (H-T) where it was
measured; otherwise it is computed from calcium and magnesium measured on the same day
(2.497 Ca + 4.118 Mg, preferring dissolved over total values). Values below 0.5 mg/l and above
10,000 mg/l are dropped as implausible for fresh water. Each station gets its median, sample count
and date range, and an annual median series for the page.

Writes
  data/stations.csv          one row per station with a hardness median
  data/station_series.json   compact annual medians per station (for the click charts)
  work/samples.parquet-like  work/samples.csv.gz (temporary, every usable sample)
  work/r_stations.json       numbers for results.json
Run: python3 code/01_stations.py
"""
import gzip, io, json, os, sys
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import *  # noqa

G = json.load(gzip.open(GEMSTAT))
F = G["files"]
R = {"source": G["source"], "retrieved": G["retrieved"][:10]}


def rd(name):
    return pd.read_csv(io.StringIO(F[name]), low_memory=False, encoding_errors="replace")


st = rd("GEMStat_station_metadata.csv")
st = st.rename(columns={"GEMS Station Number": "id", "Country Name": "country", "Water Type": "wtype", "Latitude": "lat", "Longitude": "lon",
                        "Station Identifier": "sname", "Water Body Name": "water_body", "Main Basin": "basin"})
R["stations_in_metadata"] = int(len(st))
st = st.drop_duplicates("id", keep="first")   # a few stations are listed twice with the same details
R["stations_unique_ids"] = int(len(st))

H = rd("Hardness.csv"); CA = rd("Calcium.csv"); MG = rd("Magnesium.csv")
COLS = {"GEMS Station Number": "id", "Sample Date": "date", "Parameter Code": "code", "Value": "v", "Unit": "unit",
        "Data Quality": "q", "Value Flags": "flag", "License Information": "lic", "Depth": "depth"}
for d in (H, CA, MG):
    d.rename(columns=COLS, inplace=True)
R["rows"] = {"hardness": int(len(H)), "calcium": int(len(CA)), "magnesium": int(len(MG))}
R["licence_field_filled"] = int(sum(d["lic"].notna().sum() for d in (H, CA, MG)))
R["hardness_codes"] = {k: int(v) for k, v in H["code"].value_counts().items()}

# total hardness
ht = H[(H["code"] == "H-T") & (H["unit"] == "mg/l")].copy()
ht["src"] = "H-T"

# calcium and magnesium on the same station and date (dissolved preferred, then total)


def pick(d, pref):
    d = d[d["unit"] == "mg/l"].copy()
    d["rank"] = d["code"].map({c: i for i, c in enumerate(pref)})
    d = d[d["rank"].notna()].sort_values("rank")
    return d.groupby(["id", "date"], as_index=False).agg(v=("v", "median"), code=("code", "first"), q=("q", "first"))


ca = pick(CA, ["Ca-Dis", "Ca-Tot", "Ca-Ext"]); mg = pick(MG, ["Mg-Dis", "Mg-Tot", "Mg-Ext"])
cm = ca.merge(mg, on=["id", "date"], suffixes=("_ca", "_mg"))
cm = cm[(cm["v_ca"] >= 0) & (cm["v_mg"] >= 0)]
cm["v"] = CA_F * cm["v_ca"] + MG_F * cm["v_mg"]
R["camg_pairs"] = int(len(cm))

# cross-check H-T against Ca and Mg where both exist on the same day
htd = ht.groupby(["id", "date"], as_index=False)["v"].median()
chk = htd.merge(cm[["id", "date", "v"]], on=["id", "date"], suffixes=("_ht", "_cm"))
chk = chk[(chk["v_ht"] > 0) & (chk["v_cm"] > 0)]
chk["ratio"] = chk["v_ht"] / chk["v_cm"]
chk = chk.merge(st[["id", "country"]], on="id")
bycountry = chk.groupby("country")["ratio"].agg(["median", "count"]).reset_index()
R["crosscheck"] = {"pairs": int(len(chk)), "median_ratio": float(chk["ratio"].median()),
                   "share_within_20pct": float(((chk["ratio"] > 0.8) & (chk["ratio"] < 1.25)).mean()),
                   "by_country": [{"country": r["country"], "median_ratio": float(f"{r['median']:.4g}"), "n": int(r["count"])} for _, r in bycountry.iterrows()]}
print("H-T vs Ca+Mg ratio by country:\n", bycountry.to_string())

# Countries whose H-T values disagree with Ca + Mg by more than a factor of 1.25 (different units, for example
# meq/l) lose their H-T values; their hardness comes from Ca + Mg only.
bad = bycountry[((bycountry["median"] < 0.8) | (bycountry["median"] > 1.25)) & (bycountry["count"] >= 50)]["country"].tolist()
R["ht_dropped_countries"] = bad
badids = set(st[st["country"].isin(bad)]["id"])
R["ht_rows_dropped_units"] = int(ht["id"].isin(badids).sum())
ht = ht[~ht["id"].isin(badids)]
print("H-T dropped for", bad)

# one value per station and date: H-T if measured, else Ca + Mg
s1 = ht.groupby(["id", "date"], as_index=False).agg(v=("v", "median"), q=("q", "first"), flag=("flag", "first"))
s1["src"] = "ht"
s2 = cm[["id", "date", "v", "q_ca"]].rename(columns={"q_ca": "q"})
s2["src"] = "camg"; s2["flag"] = np.nan
s2 = s2.merge(s1[["id", "date"]], on=["id", "date"], how="left", indicator=True)
s2 = s2[s2["_merge"] == "left_only"].drop(columns="_merge")
S = pd.concat([s1, s2], ignore_index=True)
n0 = len(S)
S = S[(S["v"] >= 0.5) & (S["v"] <= 10000)]
R["samples_dropped_implausible"] = int(n0 - len(S))
S = S.merge(st[["id", "lat", "lon"]], on="id")
S = S[S["lat"].notna() & S["lon"].notna() & ~((S["lat"] == 0) & (S["lon"] == 0))]
S["year"] = S["date"].str[:4].astype(int)
R["samples"] = int(len(S)); R["samples_from_ht"] = int((S["src"] == "ht").sum()); R["samples_from_camg"] = int((S["src"] == "camg").sum())
R["quality"] = {k: int(v) for k, v in S["q"].fillna("Unknown").value_counts().items()}
R["below_detection_flags"] = int((S["flag"] == "<").sum())
R["years"] = [int(S["year"].min()), int(S["year"].max())]
S[["id", "date", "v", "src", "q"]].to_csv(os.path.join(WORK, "samples.csv.gz"), index=False, compression="gzip")

# station summaries
g = S.groupby("id")
sm = pd.DataFrame({"n": g.size(), "median": g["v"].median(), "p25": g["v"].quantile(0.25), "p75": g["v"].quantile(0.75),
                   "first": g["date"].min(), "last": g["date"].max(), "share_camg": g["src"].apply(lambda x: (x == "camg").mean())}).reset_index()
recent = S[S["year"] >= 2000].groupby("id")["v"].median().rename("median_2000on")
sm = sm.merge(recent, on="id", how="left")
sm = sm.merge(st[["id", "sname", "country", "wtype", "lat", "lon", "water_body", "basin"]], on="id", how="left")
sm["wtype"] = sm["wtype"].fillna("Unknown").str.replace(" station", "", regex=False)
sm["class"] = sm["median"].apply(hclass)
for c in ("sname", "water_body", "basin", "country"):
    sm[c] = sm[c].fillna("").astype(str).map(clean)
sm["country"] = sm["country"].replace({"Netherlands (-the )": "Netherlands", "United States of America (the)": "United States",
                                       "Russian Federation (the)": "Russia", "Korea (the Republic of)": "South Korea",
                                       "Czechia": "Czechia", "Bolivia (Plurinational State of)": "Bolivia",
                                       "Iran (Islamic Republic of)": "Iran", "Tanzania, United Republic of": "Tanzania",
                                       "Moldova (the Republic of)": "Moldova", "Philippines (the)": "Philippines",
                                       "Dominican Republic (the)": "Dominican Republic", "Congo (the Democratic Republic of the)": "DR Congo",
                                       "Niger (the)": "Niger", "Macedonia (the former Yugoslav Republic of)": "North Macedonia", "Sudan (the)": "Sudan", "Gambia (the)": "Gambia",
                                       "Viet Nam": "Vietnam", "Lao People's Democratic Republic (the)": "Laos",
                                       "Taiwan (Province of China)": "Taiwan", "Venezuela (Bolivarian Republic of)": "Venezuela"})
sm = sm.sort_values("id").reset_index(drop=True)
sm["median"] = sm["median"].round(1); sm["p25"] = sm["p25"].round(1); sm["p75"] = sm["p75"].round(1)
sm["median_2000on"] = sm["median_2000on"].round(1); sm["share_camg"] = sm["share_camg"].round(3)
sm.to_csv(os.path.join(DATA, "stations.csv"), index=False)
R["stations"] = int(len(sm))
R["countries"] = int(sm["country"].nunique())
R["by_wtype"] = {k: int(v) for k, v in sm["wtype"].value_counts().items()}
R["class_counts"] = [int((sm["class"] == i).sum()) for i in range(6)]
R["median_all"] = float(sm["median"].median())
print(R["stations"], "stations", R["countries"], "countries", R["samples"], "samples")

# annual medians for the click charts: {id: [[year, median, n], ...]}
ys = S.groupby(["id", "year"])["v"].agg(["median", "size"]).reset_index()
ser = {}
for sid, d in ys.groupby("id"):
    ser[sid] = [[int(r.year), round(float(r["median"]), 1 if r["median"] < 100 else 0), int(r["size"])] for _, r in d.iterrows()]
json.dump(ser, open(os.path.join(DATA, "station_series.json"), "w"), separators=(",", ":"))
R["station_years"] = int(len(ys))
R["country_medians"] = {c: [round(float(d["median"].median()), 1), int(len(d))] for c, d in sm.groupby("country")}
print(sorted(R["country_medians"].items(), key=lambda x: x[1][0]))
json.dump(R, open(os.path.join(WORK, "r_stations.json"), "w"), indent=1)
