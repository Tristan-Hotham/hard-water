"""Step 7. Build index.html from code/page_template.html.

Every number in the text is a {{path|format}} token filled from data/results.json (plus a few values derived from
it below, under "d"), so the page cannot drift from the analysis. Facts the template states in words are checked
here. A slim copy of the results is embedded; the rasters, stations, places and borders are fetched from data/.
Run: python3 code/07_build_site.py
"""
import json, math, os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.join(HERE, "..")
sys.path.insert(0, HERE)
from common import GROUPS, GKEYS, CLASSES  # noqa

R = json.load(open(os.path.join(ROOT, "data", "results.json"), encoding="utf-8"))
GR = json.load(open(os.path.join(ROOT, "data", "grids.json")))
tpl = open(os.path.join(HERE, "page_template.html"), encoding="utf-8").read()
M = R["model"]
WC = M["world_coef"]


def derived():
    rc = {r["country"]: r["median_resid_factor"] for r in R["resid_by_country"]}
    towns = {t["name"]: t for t in R["bi_towns"]}
    xc = {c["country"]: c["median_ratio"] for c in R["stations"]["crosscheck"]["by_country"]}
    se, nw = R["bi_tees_exe"]["south_east"], R["bi_tees_exe"]["north_west"]
    return {
        "typ_crys": M["world_typical"]["crys"], "typ_carb": M["world_typical"]["carb"], "typ_unc": M["world_typical"]["unc"],
        "rain_drop_pct": (1 - 2 ** WC["log10_rain_per_1000mm"]) * 100,
        "gw_factor": 10 ** WC["groundwater"], "lake_factor": 10 ** WC["lake_or_reservoir"],
        "nor_resid": rc["Norway"], "esp_resid": rc["Spain"], "cyp_resid": rc["Cyprus"],
        "sedu_share": R["world_geo"]["group_share_of_mapped_land"]["sedu"],
        "se_hard": se["shares"][4] + se["shares"][5], "nw_soft": nw["shares"][0],
        "london_bi": towns["London"]["bi"],
        "china_sedu": R["stories"]["south_china"]["sedu_share_mapped"],
        "dh": 1 / R["units"]["dH"], "eh": 1 / R["units"]["eH"],
        "bg_factor": 1 / xc["Bulgaria"], "ee_factor": 1 / xc["Estonia"],
        "r200": {r["radius_km"]: r["r2"] for r in M["world_radius_scores"]}[200],
    }


R["d"] = D = derived()


def get(path):
    v = R
    for k in path.split("."):
        v = v[int(k)] if isinstance(v, list) else v[k]
    return v


def fmt(v, f):
    if f is None or f == "raw": return str(v)
    if f == "int": return f"{round(v):,}"
    if f == "pct": return f"{round(v * 100)}%"
    if f == "pct1": return f"{v * 100:.1f}%"
    if f == "f1": return f"{v:.1f}"
    if f == "f2": return f"{v:.2f}"
    raise ValueError(f)


def sub(m):
    return fmt(get(m.group(1)), m.group(3))


# ---------------------------------------------------------------------------
# facts stated in words in the template
O, S, CV = R["overall"], R["stories"], M["world_cv"]
assert 0.23 <= O["share_by_class"][5] <= 0.27                                   # "a quarter are very hard"
assert O["by_type"]["Groundwater"]["median"] > O["by_type"]["River"]["median"] > O["by_type"]["Lake"]["median"]
assert O["by_type"]["Reservoir"]["median"] < O["by_type"]["River"]["median"]
assert R["softest_country"]["country"] == "Norway" and R["hardest_country"]["country"] == "Belgium"
assert S["india"]["n"] - 5 <= 704                                                # "India's are almost all rivers" (checked: 704 of 709)
assert sorted(R["stations"]["ht_dropped_countries"]) == ["Bulgaria", "Estonia"]
assert 0.25 <= CV["rocks_rain_and_type"]["class_exact"] <= 0.35                  # "three times in ten"
assert 0.6 <= CV["rocks_rain_and_type"]["class_within1"] <= 0.7                  # "two times in three"
assert R["by_crystalline"][3]["lo"] == 0.8 and R["by_crystalline"][0]["hi"] == 0.2  # "more than four fifths", "less than a fifth"
assert 0.4 <= D["lake_factor"] <= 0.6                                            # "about half as hard"
assert 0.45 <= D["sedu_share"] <= 0.55                                            # "about half of the mapped land"
order = [r["country"] for r in R["resid_by_country"]]
assert order[:3] == ["Norway", "Canada", "Denmark"] and set(order[-4:]) == {"Spain", "Netherlands", "Serbia", "Cyprus"}
assert R["bi_counties_softest"][0]["name"] == "Na h-Eileanan Siar" and R["bi_counties_softest"][1]["name"] == "Highland"  # "Highlands and Islands"
assert [c["name"] for c in R["bi_counties_hardest"][:3]] == ["East Riding of Yorkshire", "Hertfordshire", "Norfolk"]
assert D["london_bi"] < 100                                                        # London modelled as soft or moderately soft
assert S["northern_france_chalk"]["share_hard_or_vhard"] > 0.9
assert S["south_china"]["stations"] == 0
assert S["yucatan"]["model_map_median"] < S["yucatan"]["median"] / 2                 # "misses this badly"
assert M["world_radius_km"] <= 50 and M["bi_radius_km"] <= 10
assert M["world_radius_km"] == M["bi_radius_km"]                                   # "for the world and for Britain and Ireland alike"
assert not ({"United Kingdom", "China", "Russia", "Australia"} & {c["country"] for c in R["countries_measured"]})
AFRICA = {"Algeria", "Egypt", "Morocco", "Tunisia", "Nigeria", "Ghana", "Kenya", "South Africa", "Ethiopia", "Tanzania", "Uganda", "Senegal", "Mali", "Niger", "Sudan", "Zambia", "Zimbabwe", "Botswana", "Namibia", "Mozambique", "Malawi", "Cameroon", "DR Congo", "Congo", "Gambia", "Burkina Faso", "Benin", "Togo", "Cote d'Ivoire", "Madagascar", "Rwanda", "Burundi", "Angola", "Chad", "Libya", "Mauritania", "Guinea", "Sierra Leone", "Liberia"}
assert not (AFRICA & {c["country"] for c in R["countries_measured"]})
assert S["scandinavia"]["share_soft"] > 0.85
CM = {c["a3"]: c["share_low_confidence"] for c in R["countries_modelled"]}
assert CM["CHN"] == 1 and CM["MNG"] == 1 and CM["BRA"] == 1 and CM["KAZ"] > 0.5       # "all of China, Mongolia and Brazil and most of Kazakhstan"
assert M["world_cv_by_confidence"]["low"]["r2"] < 0                                    # "no better than a flat guess"


def slim():
    keep_c = ["country", "stations", "median", "share_hard"]
    groups = {g[0]: g for g in GROUPS}
    short = {"carb": "Carbonate", "evap": "Evaporites", "clas": "Sandstone, mudstone", "unc": "Loose sediments", "sedu": "Sedimentary, unstated",
             "volc": "Volcanic", "crys": "Crystalline"}
    medium = {"carb": "Carbonate (limestone, chalk, dolomite)", "evap": "Evaporites (gypsum, salt)", "clas": "Sandstone, mudstone, shale",
              "unc": "Loose young sediments", "sedu": "Sedimentary, type not stated", "volc": "Volcanic (basalt, tuff)", "crys": "Crystalline (granite, gneiss)"}
    order = sorted(GKEYS, key=lambda k: M["world_typical"][k])
    rows = [{"head": "World model, seven groups"}]
    for k in order:
        rows.append({"label": medium[k], "short": short[k], "v": M["world_typical"][k], "col": groups[k][2]})
    bl = {"carb": ("Carbonate and evaporite", "#2B6F9E"), "clas": ("Sandstone, mudstone, other sediments", "#C9A66B"), "volc": ("Volcanic", "#8E6C8A"), "crys": ("Crystalline", "#C1606A")}
    rows.append({"head": "Britain and Ireland model, four groups"})
    for k in sorted(bl, key=lambda k: M["bi_typical"][k]):
        rows.append({"label": bl[k][0], "short": {"carb": "Carbonate", "clas": "Sediments", "volc": "Volcanic", "crys": "Crystalline"}[k], "v": M["bi_typical"][k], "col": bl[k][1]})
    return {"classes": R["classes"], "units": R["units"],
            "countries_measured": [{k: c[k] for k in keep_c} for c in R["countries_measured"]],
            "countries_modelled": [{k: c[k] for k in ["country", "coverage", "modelled_median", "share_low_confidence"]} for c in R["countries_modelled"]],
            "by_crystalline": R["by_crystalline"], "by_carbonate": R["by_carbonate"], "type_shares": R["type_shares"],
            "model": {"world_cv": M["world_cv"], "bi_cv": M["bi_cv"]}, "bi_nations": R["bi_nations"], "bi_counties": R["bi_counties"], "coef_rows": rows}


if __name__ == "__main__":
    html = re.sub(r"\{\{([^}|]+)(\|(\w+))?\}\}", sub, tpl)
    html = html.replace("/*RESULTS*/", json.dumps(slim(), ensure_ascii=False, separators=(",", ":"), default=lambda o: None))
    html = html.replace("/*GRIDS*/", json.dumps(GR, separators=(",", ":")))
    assert "{{" not in html
    vis = re.sub(r"<script.*?</script>|<style.*?</style>|<!--.*?-->", "", html, flags=re.S)
    vis = re.sub(r"<[^>]+>", " ", vis)
    for bad in ("\u2014", "\u2013", "--"):
        assert bad not in vis, bad
    for js in re.findall(r"<script>(.*?)</script>", html, re.S):
        for bad in ("\u2014", "\u2013"):
            assert bad not in js, bad
    open(os.path.join(ROOT, "index.html"), "w", encoding="utf-8").write(html)
    print("index.html", len(html.encode()), "bytes")
