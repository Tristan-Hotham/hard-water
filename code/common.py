"""Shared settings for Hard water: paths, hardness classes and units, and the rules that turn a
rock description into shares of seven rock groups.

Inputs are read in place (nothing large is copied into the project):
  UNEP GEMS/Water GFQA v3 open subset (CC BY 4.0)        gemstat_hardness.json.gz
  Macrostrat carto vector tiles at zoom 4 (CC BY 4.0)    macrostrat_carto_z4.json.gz
  BGS Geology 625k bedrock (Open Government Licence v3)  bgs625k.json.gz
  GSI / GSNI 1:500k bedrock (CC BY 4.0 / OGL v3)         IE_GSI_GSNI_..._500K_IE32_ITM.zip
  CRU CL 2.0 annual rainfall on a 0.5 degree grid, as derived and published by the sun-and-rain project
  Natural Earth 1:10m (public domain), GeoNames cities15000 (CC BY 4.0), OS Open Names (OGL v3)
"""
import math, os, re

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA = os.path.join(ROOT, "data")
WORK = os.path.join(ROOT, "work")          # temporary, not published
os.makedirs(WORK, exist_ok=True)

D15 = os.environ.get("D15", "/home/claude/hist/_data15")
D14 = os.environ.get("D14", "/home/claude/hist/_data14")
D13 = os.environ.get("D13", "/home/claude/hist/_data13")
GEMSTAT = os.path.join(D15, "gemstat_hardness.json.gz")
MACRO = os.path.join(D14, "macrostrat_carto_z4.json.gz")
BGS = os.path.join(D14, "bgs625k.json.gz")
GSI_ZIP = os.path.join(D14, "IE_GSI_GSNI_Bedrock_Geology_Datasets_500K_IE32_ITM.zip")
RAIN = os.environ.get("RAIN", "/home/claude/hist/sun-and-rain/data/grid_annual_30min.csv")
NE = os.environ.get("NE", "/home/claude/hist/_data3/ne/geojson")
GEONAMES = os.path.join(D13, "cities15000.zip")
TERRARIUM = os.path.join(D13, "terrarium_tiles.json.gz")
OSNAMES = os.environ.get("OSNAMES", "/home/claude/hist/osnames_gb.tsv")

# ---------------------------------------------------------------------------
# Hardness classes (mg/l as CaCO3, the same as ppm), as on the familiar England and Wales maps.
CLASSES = [  # key, label, lower bound (inclusive), upper bound (exclusive), colour (light to deep blue, my own)
    ("soft", "Soft", 0, 50.5, "#E6F2F7"),
    ("msoft", "Moderately soft", 50.5, 100.5, "#BCDDEA"),
    ("shard", "Slightly hard", 100.5, 150.5, "#86C1D9"),
    ("mhard", "Moderately hard", 150.5, 200.5, "#4D9CC2"),
    ("hard", "Hard", 200.5, 300.5, "#22709E"),
    ("vhard", "Very hard", 300.5, 1e9, "#0E3F66"),
]
CLASS_RANGES = ["0 to 50", "51 to 100", "101 to 150", "151 to 200", "201 to 300", "over 300"]
BREAKS = [50.5, 100.5, 150.5, 200.5, 300.5]


def hclass(v):
    """Index 0 to 5 of the class for a hardness in mg/l CaCO3."""
    for i, b in enumerate(BREAKS):
        if v < b:
            return i
    return 5


# Unit conversions from mg/l CaCO3.
UNITS = {"ppm": 1.0, "dH": 1 / 17.848, "eH": 1 / 14.254, "fH": 1 / 10.0}
CA_F, MG_F = 2.497, 4.118   # mg/l CaCO3 per mg/l of Ca and of Mg

# ---------------------------------------------------------------------------
# Rock groups for the model. Each map unit gets shares of these groups from its description.
GROUPS = [  # key, label, colour
    ("carb", "Carbonate rocks (limestone, chalk, dolomite, marble, marl)", "#2B6F9E"),
    ("evap", "Evaporites (gypsum, anhydrite, salt)", "#7A4F9A"),
    ("clas", "Sandstone, mudstone and other hard sedimentary rock", "#C9A66B"),
    ("unc", "Loose young sediments (alluvium, sand, gravel, till, loess)", "#E3D3A8"),
    ("sedu", "Sedimentary rock, type not stated on the map", "#B9B2A0"),
    ("volc", "Volcanic rocks (basalt, andesite, tuff)", "#8E6C8A"),
    ("crys", "Granite, gneiss, schist, slate and other crystalline rock", "#C1606A"),
]
GKEYS = [g[0] for g in GROUPS]

WORDS = [  # (regex on the lower-case text, group, weight of a carbonate word meaning only partly carbonate)
    (r"\bevaporit|\bgypsum|\banhydrit|\bhalite|\brock salt|\bsalt\b|\bsalina", "evap"),
    (r"\bmarl", "marl"),
    (r"limestone|dolomit|dolostone|\bchalk|marble|carbonate|calcarenite|travertine|tufa\b|calcrete|\bcalc", "carb"),
    (r"basalt|andesit|rhyolit|dacit|\btuff|volcan|\blava|pyroclast|trachyt|ignimbrit|tholeiit|latite|phonolit|extrusive|agglomerat|spilite", "volc"),
    (r"granit|gneiss|schist|gabbr|diorit|tonalit|monzonit|syenit|plutonic|intrusive|quartzite|migmatit|phyllite|slate|amphibol|"
     r"serpentin|ultramafic|anorthosit|charnockit|granulit|metamorph|crystalline|norite|diabase|dolerite|felsite|porphyr|greenstone|"
     r"felsic|mafic|igneous|metasediment|metavolcan|eclogite|komatiite|hornfels|mylonit|orthogneiss|paragneiss|appinite|granophyre|"
     r"pelit|psammit|granitoid|ultramafit|mafite|calcsilicate", "crys"),
    (r"alluvi|\bsand\b|gravel|\bsilt\b|\bclay\b|unconsolidated|quaternary|\bloess|\btill\b|glacial|moraine|colluvi|fluvial deposit|"
     r"lacustrine|aeolian|eolian|dune|\bpeat|surficial|regolith|laterite|bauxite|lithomarge", "unc"),
    (r"sandstone|shale|mudstone|siltstone|claystone|conglomerat|greywacke|graywacke|wacke|arkose|chert|\bcoal|lignite|turbidit|"
     r"argillite|argillaceous|breccia|flysch|molasse|redbed|red bed|siliciclastic|clastic|iron.formation|diamictite|tillite|"
     r"arenite|grit|melange|glauconit|radiolarit|phosphorit", "clas"),
    (r"sedimentary|sediments?\b|supracrustal", "sedu"),
]
WRX = [(re.compile(p), g) for p, g in WORDS]


def word_group(w):
    """Group of one rock word or short phrase, or None. 'marl' and partial carbonate words split."""
    w = w.lower().strip()
    if not w:
        return None
    if "calcsilicate" in w:
        return {"crys": 1.0}
    # qualifiers that make a carbonate word only a minor part
    for rx, g in WRX:
        if rx.search(w):
            if g == "carb" and re.search(r"calcareous|argillaceous limestone|cherty limestone|limy", w) and not re.search(r"^limestone", w):
                return {"carb": 0.5, "clas": 0.5}
            if g == "carb" and re.search(r"^(meta)?(volcanic|basalt|schist|gneiss)", w):
                break
            if g == "marl":
                return {"carb": 0.5, "clas": 0.5}
            if g == "carb" and "marble" in w and re.search(r"schist|gneiss|psammit|pelit", w):
                return {"carb": 0.3, "crys": 0.7}
            return {g: 1.0}
    return None


def split_list(t):
    parts = re.split(r",|;|/|\band\b|&|\bor\b|\bwith\b|\+", t)
    return [p.strip() for p in parts if p.strip()]


def shares_from_text(lith, name="", minor_text=None):
    """Shares of the six groups (a dict summing to 1) from a lithology text, or None if unknown.

    'Major:{a,b}, Minor{c,d}' (Macrostrat): majors share 0.75 equally, minors 0.25 equally.
    Otherwise the listed rock words share equally; words such as 'thin', 'minor' or 'subsidiary'
    give that word a quarter weight. The unit name is used only when the text says nothing useful.
    """
    t = (lith or "").lower()
    out = {}

    def add(words, total):
        got = []
        for w in words:
            wt = 0.25 if re.search(r"\b(minor|thin|subsidiary|rare|local|some|interbedded)\b", w) else 1.0
            g = word_group(w)
            if g:
                got.append((g, wt))
        s = sum(wt for _, wt in got)
        for g, wt in got:
            for k, v in g.items():
                out[k] = out.get(k, 0) + total * v * wt / s
        return bool(got)

    m = re.search(r"major:\{([^}]*)\}", t)
    if m:
        mi = re.search(r"minor:?\{([^}]*)\}", t)
        okM = add(split_list(m.group(1)), 0.75 if mi else 1.0)
        if mi:
            okm = add(split_list(mi.group(1)), 0.25)
            if not okM and okm:
                pass
    else:
        # "A with subordinate B and C" or "A with minor B": everything after the qualifier is minor
        mm = re.search(r"\bwith (subordinate|minor)\b|\[subequal/subordinate\]|\[subordinate\]|\bsubordinate\b", t)
        if mm:
            head, tail = t[:mm.start()], t[mm.end():]
            okh = add(split_list(head), 0.75)
            okt = add(split_list(tail), 0.25)
        else:
            add(split_list(t), 1.0)
    if not out:
        n = (name or "").lower()
        if "carbonate" in n:
            out = {"carb": 1.0}
        elif re.search(r"plutonic|intrusive|metamorphic|crystalline", n):
            out = {"crys": 1.0}
        elif re.search(r"volcanic|extrusive", n):
            out = {"volc": 1.0}
        elif re.search(r"unconsolidated|quaternary|alluvi", n):
            out = {"unc": 1.0}
        elif re.search(r"clastic|siliciclastic", n):
            out = {"clas": 1.0}
        elif re.search(r"sedimentary", n):
            out = {"sedu": 1.0}
        else:
            g = word_group(n) if n else None
            out = g or {}
    s = sum(out.values())
    if s <= 0:
        return None
    return {k: v / s for k, v in out.items()}


def clean(s):
    """No dashes in visible text."""
    s = (s or "").replace("�", "?")   # characters already lost in the source bundle
    s = s.replace("\u2014", ", ").replace("\u2013", " to ").replace("-" * 2, " to ")
    return re.sub(r"\s+", " ", s).strip()


def haversine_km(lat1, lon1, lat2, lon2):
    import numpy as np
    p1, p2 = np.radians(lat1), np.radians(lat2)
    a = np.sin((p2 - p1) / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(np.radians(lon2 - lon1) / 2) ** 2
    return 2 * 6371.0088 * np.arcsin(np.sqrt(a))
