# Hard water: the world's water hardness and the rocks behind it

A world map of water hardness from the open UNEP GEMS/Water archive of river, lake, reservoir and groundwater monitoring, with a simple, cross-validated estimate of hardness from the rocks and rainfall wherever nobody has measured, and a sharper 1 km estimate for Britain and Ireland. The six classes are those of the familiar England and Wales maps: soft (0 to 50 mg/l as CaCO3), moderately soft (51 to 100), slightly hard (101 to 150), moderately hard (151 to 200), hard (201 to 300) and very hard (over 300).

**Interactive version:** https://tristan-hotham.github.io/hard-water/

**Paper and archived release:** https://doi.org/10.5281/zenodo.23735824

A data study by Dr Tristan Hotham (Independent researcher, ORCID 0000-0003-4931-9349).

**This is not a tap water map.** The measurements are of rivers, lakes and groundwater at monitoring stations, and the modelled layer is an estimate from geology and rainfall. Tap water depends on where a supplier draws its water, how it blends sources and whether it softens it. No single open dataset of tap water hardness by supply zone exists for the UK; each water company publishes its own lookup.

## Headline results

- 16,682 stations in 34 countries with 611,858 usable samples (1906 to 2023). The median station is 202 mg/l; 26% are very hard and 12% soft. Groundwater is hardest (median 276 mg/l, 3,861 stations), then rivers (199 mg/l, 10,515), lakes (139 mg/l) and reservoirs (111 mg/l).
- Norway, Sweden and Finland: 571 stations, median 14 mg/l, 90% soft. The Yucatán Peninsula: 187 stations (131 wells), median 514 mg/l, 92% very hard. Northern France's chalk: 141 stations, median 301 mg/l.
- Crystalline rock is the clearest signal: stations with more than four fifths crystalline rock within 10 km have a median of 34 mg/l, against 220 mg/l with less than a fifth.
- The model (rock group shares within 10 km, annual rainfall, type of water) under spatial cross-validation by 5 degree squares: 29% of stations in the right class, 64% within one class (38% for a constant guess), typical error a factor of 1.8, R² 0.24 on log hardness. Leaving out whole countries: 61% within one class, R² 0.10.
- Fitted typical hardness of river water from one rock group alone at 1,000 mm of rain: crystalline 22 mg/l, sandstone and mudstone 67, volcanic 72, carbonate 211, loose young sediments 257, evaporites 512. Each doubling of rainfall makes water about 21% softer; groundwater is about twice as hard as river water nearby.
- About half (50%) of the world's mapped land is described by the world geology only as sedimentary rock of a given age, which is the model's main weakness (it predicts about 117 mg/l on the Yucatán).
- Confidence: where at least half the land in a station's window comes from detailed map units, the model gets 68% within one class (R² 0.36, 10,363 stations); where it rests mainly on coarse geology, 56% (R² -0.12, 4,430 stations). 71% of the modelled land is marked low confidence and shown faded and hatched, including all of China, Mongolia and Brazil.
- Britain and Ireland at 1 km, calibrated on 388 Irish stations (the open archive has none in Great Britain): 82% of Irish stations within one class, R² 0.60. Modelled medians: England 75 mg/l, Wales 38, Scotland 31, Northern Ireland 100, Republic of Ireland 105; south and east of the Tees to Exe line 92 mg/l, north and west 35. The estimate misses London (modelled 40 mg/l), whose water comes from the chalk and from rivers draining far upstream.

Every number is in `data/results.json`, which the page, paper, blog and collage are built from.

## Files

| Path | What it is |
|---|---|
| `index.html` | The interactive page: map with measured stations, modelled layer and the 1 km Britain and Ireland layer, a model check view, units (mg/l, German, Clark and French degrees), place search, station charts over time, charts and tables. |
| `header.png` | Social card (1200 by 630). |
| `data/stations.csv` | Every station with a hardness figure: country, type of water, position, median, quartiles, number of samples, first and last date, share computed from calcium and magnesium, median since 2000, class. |
| `data/station_series.json` | Yearly medians and sample counts per station (for the page's charts). |
| `data/station_model.csv` | Per station: rock group shares within the window, rainfall, model prediction, cross-validated prediction and residual; for Irish stations the 1 km model too. |
| `data/model_coefficients.csv` | Both models' fitted parameters. |
| `data/units_world.csv`, `data/units_bi.csv` | Every Macrostrat, BGS and GSI map unit with its rock description and the rock group shares given to it (world: also its source map and whether it counts as detailed). |
| `data/countries_measured.csv`, `data/countries_modelled.csv` | Country tables from the stations and from the modelled layer. |
| `data/bi_areas.csv` | Britain and Ireland: modelled class shares by nation and by county, unitary area or council area. |
| `data/places.csv`, `data/places.json` | GeoNames places over 15,000 people with their modelled values (world and 1 km). |
| `data/world_model.png`, `data/bi_model.png`, `data/grids.json` | The modelled layers as Web Mercator images: R = hardness code (mg/l = 10^((code − 1) / 249 × log10 3000)), G = hill shade, B = 2 estimate, 3 low-confidence estimate, 1 land without one, 0 sea (the Britain and Ireland image has no low-confidence class). |
| `data/page_stations.json`, `data/borders.json` | Compact station list and Natural Earth borders for the page. |
| `data/results.json` | Every number used anywhere. |
| `code/` | The pipeline, run in order (below). |
| `reddit/` | Collage (`collage.py`, PNG) and post text. |
| `release/` | Paper PDF and source, blog post and figures, Zenodo files. |

## Method in brief

- **Stations.** Each sample is the measured total hardness, or 2.497 × Ca + 4.118 × Mg (mg/l) measured on the same day, preferring dissolved values. Bulgaria and Estonia report total hardness in other units (about 50 and 2,491 times smaller than their calcium and magnesium imply) and use calcium and magnesium only. Values below 0.5 or above 10,000 mg/l are dropped. Each station's figure is the median of its samples. The licence field is empty on every row of the open subset, which is released as a whole under CC BY 4.0.
- **Rocks.** Every map unit gets shares of seven rock groups (carbonate, evaporite, sandstone and mudstone, loose young sediments, sedimentary not stated, volcanic, crystalline) from its lithology text by fixed word rules (`code/common.py`). World: Macrostrat zoom-4 tiles summed to 0.1 degrees. Britain: BGS 625k bedrock; island of Ireland: GSI and GSNI 1:500k bedrock; 0.01 by 0.015 degree grid (about 1 km). Superficial deposits are not used.
- **Model.** H = (Σ h_k s_k) × (rain / 1000 mm)^c × f_type, with s_k the share of group k in a square window reaching R km either side; fitted by least squares on log10 hardness, stations weighted by one over the number in their one degree square, h_k bounded 5 to 3,000 mg/l. R chosen by spatial cross-validation among windows up to 50 km (world) and 10 km (Britain and Ireland); both chose 10 km. Windows of 100 and 200 km scored about the same or a little better (R² up to 0.27) but blur the map. Land with under 200 mm of rain a year is left blank. Britain and Ireland: four merged groups, fitted to Irish stations with the world rainfall exponent held fixed.
- **Confidence.** The zoom-4 tiles stitch together six regional maps. Two (the Generalized Geology of the World and Africa) are very coarse, with country-sized units and straight clip edges where they meet the detailed maps; the Europe map also has plate-scale legend units east of the Volga. A unit counts as detailed if it is from a regional map, not named only after a tectonic plate or sea, and less than half "sedimentary, type not stated". An estimate is low confidence when under half of the mapped land in its window is detailed. Straight seams between source maps remain visible on the modelled layer, for example along about 50 and 60 degrees north in Asia and Russia.
- **Validation.** Ten folds of 5 by 5 degree squares (half-degree squares in Ireland), and leave-one-country-out.

## Known limits

- Monitoring stations are not taps; supply zones blend and transfer water and may soften it.
- The open archive covers 34 countries unevenly, none in Great Britain, China, Russia, Africa or Australia; country figures reflect where stations are. Station medians pool all years and seasons.
- The world geology is coarse and often names only the age of the rock; the word rules will misread some units. A square window is not a catchment. The model knows nothing of soils, evaporation, land use, pollution or treatment.
- The Britain and Ireland estimate is calibrated on Irish stations only and could not be checked against any British measurement.
- Some station names lost accented letters in the copy of the archive used here; they show as "?".

## How to rerun

Inputs are read in place: `gemstat_hardness.json.gz` (folder `D15`), `macrostrat_carto_z4.json.gz`, `bgs625k.json.gz` and the GSI zip (`D14`), `terrarium_tiles.json.gz` and `cities15000.zip` (`D13`), Natural Earth GeoJSON (`NE`), the half-degree CRU CL 2.0 annual grid `grid_annual_30min.csv` from https://github.com/Tristan-Hotham/sun-and-rain (`RAIN`) and the county boundaries `counties.geojson` from https://github.com/Tristan-Hotham/how-roman-is-your-road.

```
python3 code/01_stations.py
python3 code/02_world_geology.py
python3 code/03_bi_geology.py
python3 code/04_model.py
python3 code/05_analyse.py
python3 code/06_page_data.py
python3 code/07_build_site.py
FONT_DIR=path/to/fonts python3 code/08_figures.py
python3 code/09_paper.py
python3 code/10_release.py
FONT_DIR=path/to/fonts python3 reddit/collage.py
```

Needs Python 3 with numpy, scipy, pandas, Pillow, matplotlib, mapbox_vector_tile, pyshp and pyproj, and xelatex with TeX Gyre Pagella and Inter for the paper. Temporary grids go in `work/` (not published).

## Sources and licences

- UNEP GEMS/Water, Global Freshwater Quality Archive v3 (GFQA_v3), open subset, https://doi.org/10.5281/zenodo.18459694 (CC BY 4.0).
- Macrostrat geological map tiles, https://macrostrat.org (CC BY 4.0); Peters, S.E., Husson, J.M. and Czaplewski, J. (2018) Geochemistry, Geophysics, Geosystems 19. Compiled from regional maps including the Geologic Map of North America (Garrity and Soller 2009, USGS DS 424), the Generalized Geology of the World (Chorlton 2007, GSC Open File 5529) and the International Geological Map of Europe (Asch 2003).
- BGS Geology 625k (DiGMapGB-625), British Geological Survey (Open Government Licence v3.0).
- Geological Survey Ireland and Geological Survey of Northern Ireland 1:500k bedrock geology (CC BY 4.0; Open Government Licence).
- CRU CL 2.0, Climatic Research Unit, University of East Anglia: New, M., Lister, D., Hulme, M. and Makin, I. (2002) Climate Research 21: 1 to 25. No licence is stated with the files; only a derived half-degree annual grid is used.
- GeoNames cities15000 (CC BY 4.0); Natural Earth 1:10m and 1:50m (public domain); ONS local authority boundaries (Open Government Licence v3.0); AWS Terrain Tiles (Mapzen and Tilezen; ETOPO1, GMTED2010, SRTM and others) for hill shading.
- Map library: Leaflet 1.9.4. Base map tiles © Esri when they load.

Derived data: CC BY 4.0 with the credits above. Text and figures: CC BY 4.0. Code: MIT (see `LICENSE`).

## Citation

Hotham, T. (2026) Hard Water: Measured and Modelled Hardness of the World's Rivers, Lakes and Groundwater, and the Rocks Behind It. Zenodo. https://doi.org/10.5281/zenodo.23735824

I wrote the text, with drafting and editing assistance from Claude (Anthropic), which also helped prepare the analysis code. I reviewed and checked all text and results.
