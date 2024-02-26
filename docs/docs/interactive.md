---
layout: default
title: Interactive Visualization
nav_order: 4
---

# Interactive Visualization
{: .no_toc }


The example `iframe` below is the product of iether [interactiveOnly](https://github.com/AdrianKriger/geo3D/blob/main/interactiveOnly.ipynb) or [osm_LoD1_3DCityModel-walkthrough](https://github.com/AdrianKriger/geo3D/blob/main/osm_LoD1_3DCityModel-walkthrough.ipynb). Building stock is differentiated through color. A school, housing, retail, healthcare and community focused facilities are easily identified while the tooltips highlight the underlying data. Additional features unique to an aoi can also be included. Here farmland, streams, recreational spaces and bus rapid transit routes have been added *- you are thus limited only through data and your imagination*. <!-- {: .fs-6 .fw-300 } -->

**To navigate on a laptop without a mouse**:

- `trackpad left-click drag-left` and `-right`;
- `Ctrl left-click drag-up`, `-down`, `-left` and `-right` to rotate and so-on and
- `+` next to Backspace zoom-in and `-` next to `+` zoom-out.


<iframe src="{{site.baseurl | prepend: site.url}}/img/interactiveOnly.html" style="width: 800px; height: 400px; border: 0px"></iframe>

The visualisation above employs the default [Carto Dark Matter](https://github.com/CartoDB/basemap-styles) basemap. [Pydeck](https://deckgl.readthedocs.io/en/latest/index.html) supports a number of [map_styles](https://deckgl.readthedocs.io/en/latest/deck.html) including the extensive [mapbox gallery](https://www.mapbox.com/gallery/) and [Maptiler](https://www.maptiler.com/) urls (e.g.: `https://api.maptiler.com/maps/{style}/style.json?key={your API key}`).

