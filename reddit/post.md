## Title

[OC] How hard is the world's water? 16,682 river, lake and groundwater stations in 34 countries, and a simple model that guesses hardness from the rocks: soft over granite (median 34 mg/l), hard over limestone and chalk (Yucatán 514 mg/l)

## Body (optional)

I took every station with a hardness figure in the open UN GEMS/Water archive and coloured each by its median, in the six classes used on UK water company maps. Norway, Sweden and Finland: median 14 mg/l. The limestone of the Yucatán: 514 mg/l. Then I fitted a simple model from the rocks within 10 km and the rainfall. Tested on places it never saw, it gets 64% of stations within one class. The archive has no stations in Great Britain, so for Britain and Ireland I used 1 km bedrock maps calibrated on 388 Irish stations (82% within one class): the chalk belt appears from geology alone.

Caveat: these are rivers, lakes and wells, not taps, and the shading is an estimate. It misses London, whose water comes from far upstream.

Interactive version: https://tristan-hotham.github.io/hard-water/

## Source comment

Source: UNEP GEMS/Water, Global Freshwater Quality Archive v3, open subset (Zenodo 10.5281/zenodo.18459694, CC BY 4.0): hardness, calcium and magnesium, 611,858 samples, station medians; Macrostrat geological map tiles at zoom 4 (CC BY 4.0); BGS Geology 625k bedrock (Open Government Licence v3.0); Geological Survey Ireland and Geological Survey of Northern Ireland 1:500k bedrock (CC BY 4.0; Open Government Licence); CRU CL 2.0 annual rainfall 1961 to 1990 (New et al. 2002); Natural Earth (public domain); hill shading from AWS Terrain Tiles.

Tools: Python (pandas, NumPy, SciPy for the model fit, mapbox_vector_tile, pyshp and pyproj for the geology, Pillow, matplotlib).

Interactive: https://tristan-hotham.github.io/hard-water/
Blog post: https://blog.tristanhotham.com/2026/10/28/hard-water/
Paper: https://doi.org/10.5281/zenodo.23735824
