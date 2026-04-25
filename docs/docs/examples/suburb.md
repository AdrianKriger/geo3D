---
layout: default
title: Suburb
parent: Examples
nav_order: 4
---

# Suburb
{: .no_toc }

---

[osm3DwStock_param.json](https://github.com/AdrianKriger/geo3D/blob/main/suburb/osm3DwStock_param.json) defines the settings to produce a CityJSON 3D City Model of a high-density urban suburb with an approximate population of 17 000 (3 700 buildings). 

`runtime: 7:09:59` with a 25-m elevation model and `1:57:28` with a 5-m elevation model.

[osm_LoD1_3DCityModel-walkthrough.ipynb](https://github.com/AdrianKriger/geo3D/blob/main/suburb/osm_LoD1_3DCityModel-walkthrough.ipynb) will prduce a LoD1 3D City Model.  
Parse the model through [CityJSONspatialDataScience.ipynb](https://github.com/AdrianKriger/geo3D/blob/main/suburb/CityJSONspatialDataScience.ipynb) to estimate population, BVPC and produce of a dynamic-colorized `.html`.

Due to the volume of the data this processing option will typically require an 'osm.pbf'.  
If internet access is a challange in your region **No Internet** options are possible.  

Example data (raster DEM, `osm.pbf` and `.mbtiles`) is available in the `./data` folder.
