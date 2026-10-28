"""Step 10. Blog post (release/blog/post.html from code/blog_template.html), blog and Zenodo meta files,
and the Zenodo folder (paper, key CSVs, results.json and a code-and-data zip).
Image placeholders {{IMG|file|alt|caption}} are left for the uploader (numbers inside them are filled).
Run: python3 code/10_release.py
"""
import importlib.util, json, os, re, shutil, zipfile

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.join(HERE, "..")
spec = importlib.util.spec_from_file_location("b", os.path.join(HERE, "07_build_site.py")); b = importlib.util.module_from_spec(spec); spec.loader.exec_module(b)
R = b.R
t = open(os.path.join(HERE, "blog_template.html"), encoding="utf-8").read()
t = re.sub(r"\{\{(?!IMG\|)([^}|]+)(\|(\w+))?\}\}", b.sub, t)
assert "{{" not in t.replace("{{IMG|", "")
vis = re.sub(r"<!--.*?-->", "", t, flags=re.S)
for bad in ("\u2014", "\u2013", "--"):
    assert bad not in vis, bad
os.makedirs(os.path.join(ROOT, "release", "blog"), exist_ok=True)
open(os.path.join(ROOT, "release", "blog", "post.html"), "w", encoding="utf-8").write(t)
for m in re.findall(r"\{\{IMG\|([^|]+)\|", t):
    assert os.path.exists(os.path.join(ROOT, "release", "blog", m)), m
O, S = R["overall"], R["stories"]
meta = dict(
    title="Hard water: the world's water hardness and the rocks behind it",
    slug="hard-water",
    excerpt=(f"A world map of water hardness from {O['n']:,} river, lake and groundwater stations, and a simple model that guesses it from the rocks: "
             f"Scandinavia's median is {S['scandinavia']['median']:.0f} mg/l, the Yucatán's {S['yucatan']['median']:.0f}."),
    seo_description="Why is water soft over granite and hard over chalk? A world map of measured hardness, a rock-based estimate and a sharper map of Britain and Ireland.",
    tags=["water", "hard water", "geology", "limestone", "maps", "Britain", "interactive", "data"],
    header_alt=("Title card reading 'Hard water' in deep and mid blue on a pale water-blue background, beside a world map shaded in six blues from pale (soft) "
                "to deep (very hard), with deserts in grey. Text: The world's water hardness from "
                f"{O['n']:,} river, lake and groundwater stations, and an estimate from the rocks and rainfall where none exist. A six-step colour bar runs from soft to very hard."))
assert len(meta["seo_description"]) < 160, len(meta["seo_description"])
for v in meta.values():
    for bad in ("\u2014", "\u2013", "--"):
        assert bad not in str(v)
json.dump(meta, open(os.path.join(ROOT, "release", "blog", "meta.json"), "w"), indent=2, ensure_ascii=False)

NAME = "Hotham_Hard_Water"
pt = open(os.path.join(ROOT, "release", "paper", NAME + ".tex"), encoding="utf-8").read()
ab = pt.split("\\begin{center}\\textbf{Abstract}\\end{center}")[1].split("\\textbf{Keywords:}")[0].strip()
ab = ab.replace("\\%", "%").replace("\\&", "&").replace("$-$", "−")
assert "\\" not in ab and "$" not in ab, ab
Z = os.path.join(ROOT, "release", "zenodo")
if os.path.isdir(Z):
    shutil.rmtree(Z)
os.makedirs(Z)
zmeta = dict(
    title="Hard Water: Measured and Modelled Hardness of the World's Rivers, Lakes and Groundwater, and the Rocks Behind It",
    description_html=f"<p>{ab}</p><p>Interactive version: https://tristan-hotham.github.io/hard-water/. Data and code: https://github.com/Tristan-Hotham/hard-water. Blog post: to follow.</p>"
                     "<p>Sources: UNEP GEMS/Water Global Freshwater Quality Archive v3, open subset (CC BY 4.0); Macrostrat map tiles (CC BY 4.0); BGS Geology 625k "
                     "(Open Government Licence v3.0); Geological Survey Ireland and Geological Survey of Northern Ireland 1:500k bedrock (CC BY 4.0, Open Government Licence); "
                     "CRU CL 2.0 rainfall; GeoNames (CC BY 4.0); Natural Earth (public domain); ONS boundaries (Open Government Licence v3.0); AWS Terrain Tiles.</p>",
    keywords=["water hardness", "calcium carbonate", "GEMS/Water", "geology", "lithology", "limestone", "chalk", "spatial cross-validation", "Britain and Ireland", "open data"],
    licence="cc-by-4.0",
    related=[
        dict(identifier="https://github.com/Tristan-Hotham/hard-water", relation="issupplementedby", resource_type="software"),
        dict(identifier="https://tristan-hotham.github.io/hard-water/", relation="issupplementedby", resource_type="other"),
        dict(identifier="https://doi.org/10.5281/zenodo.18459694", relation="isderivedfrom", resource_type="dataset"),
        dict(identifier="https://macrostrat.org", relation="isderivedfrom", resource_type="dataset"),
        dict(identifier="https://www.bgs.ac.uk", relation="isderivedfrom", resource_type="dataset"),
        dict(identifier="https://www.gsi.ie", relation="isderivedfrom", resource_type="dataset"),
        dict(identifier="https://crudata.uea.ac.uk/cru/data/hrg/tmc/", relation="isderivedfrom", resource_type="dataset"),
        dict(identifier="https://doi.org/10.5281/zenodo.23732580", relation="isderivedfrom", resource_type="dataset"),
        dict(identifier="https://download.geonames.org/export/dump/", relation="isderivedfrom", resource_type="dataset"),
        dict(identifier="https://www.naturalearthdata.com", relation="isderivedfrom", resource_type="dataset"),
        dict(identifier="https://registry.opendata.aws/terrain-tiles/", relation="isderivedfrom", resource_type="dataset"),
    ])
json.dump(zmeta, open(os.path.join(Z, "meta.json"), "w"), indent=2, ensure_ascii=False)
shutil.copy(os.path.join(ROOT, "release", NAME + ".pdf"), Z)
for f in ["stations.csv", "station_model.csv", "countries_measured.csv", "countries_modelled.csv", "bi_areas.csv", "model_coefficients.csv", "units_world.csv", "units_bi.csv", "results.json"]:
    shutil.copy(os.path.join(ROOT, "data", f), Z)
with zipfile.ZipFile(os.path.join(Z, "hard-water_code_and_data.zip"), "w", zipfile.ZIP_DEFLATED) as zf:
    for sub in ("code", "data"):
        for f in sorted(os.listdir(os.path.join(ROOT, sub))):
            p = os.path.join(ROOT, sub, f)
            if os.path.isfile(p) and not f.endswith(".pyc"):
                zf.write(p, f"hard-water/{sub}/{f}")
    for f in ("README.md", "LICENSE"):
        zf.write(os.path.join(ROOT, f), f"hard-water/{f}")
words = len(re.sub(r"<[^>]+>", " ", re.sub(r"\{\{IMG\|.*?\}\}", "", re.sub(r"<!--.*?-->", "", t, flags=re.S))).split())
print("release files written; post words:", words)
