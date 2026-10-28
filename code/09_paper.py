"""Step 9. Build the paper PDF: fill code/paper_template.tex with numbers and tables from data/results.json,
write release/paper/Hotham_Hard_Water.tex and compile it with xelatex (figures from step 8), then copy the PDF to release/.
Run: python3 code/09_paper.py
"""
import json, math, os, re, shutil, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.join(HERE, "..")
sys.path.insert(0, HERE)
R = json.load(open(os.path.join(ROOT, "data", "results.json")))
OUT = os.path.join(ROOT, "release", "paper"); NAME = "Hotham_Hard_Water"
M = R["model"]; WC = M["world_coef"]
rc = {r["country"]: r["median_resid_factor"] for r in R["resid_by_country"]}
xc = {c["country"]: c["median_ratio"] for c in R["stations"]["crosscheck"]["by_country"]}
rs = {r["radius_km"]: r["r2"] for r in M["world_radius_scores"]}
se, nw = R["bi_tees_exe"]["south_east"], R["bi_tees_exe"]["north_west"]
ws = M["world_pred_class_share"]
cm = R["countries_modelled"]
exact = M["world_cv"]["rocks_rain_and_type"]["class_exact"]
R["d"] = D = {
    "typ_crys": M["world_typical"]["crys"], "typ_carb": M["world_typical"]["carb"], "typ_unc": M["world_typical"]["unc"], "typ_clas": M["world_typical"]["clas"],
    "typ_volc": M["world_typical"]["volc"], "typ_evap": M["world_typical"]["evap"],
    "bi_crys": M["bi_typical"]["crys"], "bi_clas": M["bi_typical"]["clas"], "bi_volc": M["bi_typical"]["volc"], "bi_carb": M["bi_typical"]["carb"],
    "rain_drop_pct": (1 - 2 ** WC["log10_rain_per_1000mm"]) * 100, "rain_c": WC["log10_rain_per_1000mm"],
    "gw_factor": 10 ** WC["groundwater"], "lake_factor": 10 ** WC["lake_or_reservoir"],
    "nor_resid": rc["Norway"], "esp_resid": rc["Spain"], "cyp_resid": rc["Cyprus"],
    "bg_factor": 1 / xc["Bulgaria"], "ee_factor": 1 / xc["Estonia"],
    "r10": rs[10], "r25": rs[25], "r50": rs[50], "r100": rs[100], "r200": rs[200],
    "free_rain_factor": (R["bi_gb_rain"]["p95"] / R["bi_gb_rain"]["p05"]) ** abs(M["bi_free_rain_coef"]),
    "se_hard": se["shares"][4] + se["shares"][5], "nw_soft": nw["shares"][0],
    "london_bi": {t["name"]: t for t in R["bi_towns"]}["London"]["bi"],
    "china_sedu": R["stories"]["south_china"]["sedu_share_mapped"], "sedu_share": R["world_geo"]["group_share_of_mapped_land"]["sedu"],
    "wshare_soft": ws[0] + ws[1], "wshare_hard": ws[4] + ws[5],
    "cm_min": cm[-1]["country"], "cm_min_v": cm[-1]["modelled_median"], "cm_max": cm[0]["country"], "cm_max_v": cm[0]["modelled_median"],
    "detailed_share": R["world_geo"]["detailed_share_of_mapped_land"],
    "hc_soft": M["world_pred_class_share_high_conf"][0] + M["world_pred_class_share_high_conf"][1],
    "hc_hard": M["world_pred_class_share_high_conf"][4] + M["world_pred_class_share_high_conf"][5],
    "exact_tenths": {2: "two", 3: "three", 4: "four"}[round(exact * 10)],
}
# facts stated in words
assert D["free_rain_factor"] > 10
assert abs(rs[100] - rs[10]) < 0.02 and rs[200] > rs[10]                      # "about the same at 100 km, a little better at 200 km"
assert M["world_cv"]["rocks_rain_and_type"]["class_within1"] > 0.6 and M["world_cv"]["rocks_rain_and_type"]["class_within1"] < 0.7  # "two times in three"
assert M["bi_cv"]["rocks_type_and_world_rain_effect"]["r2"] < M["bi_cv"]["rocks_type_and_free_rain_effect"]["r2"]  # "costs some accuracy"
assert M["world_cv_by_confidence"]["low"]["r2"] < 0
assert {c["a3"]: c["share_low_confidence"] for c in R["countries_modelled"]}["MNG"] == 1
assert R["by_crystalline"][3]["lo"] == 0.8 and R["by_carbonate"][4]["lo"] == 0.5


def get(path):
    v = R
    for k in path.split("."):
        v = v[int(k)] if isinstance(v, list) else v[k]
    return v


tex = lambda s: str(s).replace("&", "\\&").replace("%", "\\%").replace("_", "\\_").replace("#", "\\#")


def fmt(v, f):
    if f == "int": return f"{int(round(v)):,}"
    if f == "pct": return f"{v * 100:.0f}\\%"
    if f == "pct1": return f"{v * 100:.1f}\\%"
    if f == "f1": return f"{v:.1f}"
    if f == "f2": return f"{v:.2f}".replace("-", "$-$")
    return tex(str(v))


def sub(m):
    return fmt(get(m.group(1)), m.group(3))


# tables
c = [r for r in R["countries_measured"] if r["stations"] >= 20]
c = sorted(c, key=lambda r: r["median"])
soft, hard = c[:6], c[::-1][:6]
TAB_COUNTRIES = "\n".join(f"{tex(a['country'])} & {a['stations']:,} & {a['median']:.0f} && {tex(b['country'])} & {b['stations']:,} & {b['median']:.0f}\\\\" for a, b in zip(soft, hard))
CVROWS = [("Global median for every station", M["world_cv"]["baseline_global_median"]), ("Rainfall and type of water only", M["world_cv"]["rain_and_type_only"]),
          ("Rocks only", M["world_cv"]["rocks_only"]), ("Rocks and type of water", M["world_cv"]["rocks_and_type"]),
          ("Rocks, rainfall and type of water", M["world_cv"]["rocks_rain_and_type"]), ("The same, leaving out whole countries", M["world_cv"]["rocks_rain_and_type_leave_country_out"]),
          ("Britain and Ireland model, Irish stations", M["bi_cv"]["rocks_type_and_world_rain_effect"]), ("World model at the same Irish stations", M["bi_cv"]["world_model_at_irish_stations"])]
TAB_CV = "\n".join(f"{n} & {r['n']:,} & {r['r2']:.2f} & {r['class_exact'] * 100:.0f}\\% & {r['class_within1'] * 100:.0f}\\% & {r['median_factor']:.1f}\\\\".replace("-0.", "$-$0.") for n, r in CVROWS)
TAB_BI = "\n".join(f"{n['name']} & {n['median']:.0f} & " + " & ".join(f"{x * 100:.0f}\\%" for x in n["shares"]) + "\\\\" for n in R["bi_nations"])

if __name__ == "__main__":
    t = open(os.path.join(HERE, "paper_template.tex"), encoding="utf-8").read()
    t = t.replace("{{TAB_COUNTRIES}}", TAB_COUNTRIES).replace("{{TAB_CV}}", TAB_CV).replace("{{TAB_BI}}", TAB_BI)
    t = re.sub(r"\{\{([^}|]+)(\|(\w+))?\}\}", sub, t)
    assert "{{" not in t
    body = t.split("\\begin{document}")[1]
    assert "--" not in body and "\u2014" not in body and "\u2013" not in body
    os.makedirs(OUT, exist_ok=True)
    open(os.path.join(OUT, NAME + ".tex"), "w", encoding="utf-8").write(t)
    for _ in range(2):
        r = subprocess.run(["xelatex", "-interaction=nonstopmode", NAME + ".tex"], cwd=OUT, capture_output=True, text=True)
    if r.returncode:
        print(r.stdout[-3000:]); raise SystemExit("xelatex failed")
    for ext in (".aux", ".log", ".out"):
        p = os.path.join(OUT, NAME + ext)
        if os.path.exists(p): os.remove(p)
    shutil.move(os.path.join(OUT, NAME + ".pdf"), os.path.join(ROOT, "release", NAME + ".pdf"))
    txt = subprocess.run(["pdftotext", os.path.join(ROOT, "release", NAME + ".pdf"), "-"], capture_output=True, text=True).stdout
    print("dashes:", txt.count("\u2014"), txt.count("\u2013"), txt.count("--"), "words:", len(txt.split()))
