---
layout: default
title: Processing Strategies
parent: geo3D
nav_order: 2
---

# Processing Strategies
<!-- {: .no_toc } -->
&nbsp;

<p align="center"><b>There are two procesing strategies</b></p>

<!--| [Village/Campus](https://github.com/AdrianKriger/osm_LoD1_3DCityModel/tree/main/village_campus) | [District](https://github.com/AdrianKriger/osm_LoD1_3DCityModel/tree/main/districts) *[on hold]*|
| :-----: | :-----: |
| [village/campus]((https://github.com/AdrianKriger/osm_LoD1_3DCityModel/tree/main/village_campus)) is designed for extremely focused analysis at a **neighbourhood** level. These are areas <br /> with a population of no more than 10 000| For **larger** areas with populations of more <br /> than 10 000;  one or many suburbs, census wards or tracts; please execute [districts](https://github.com/AdrianKriger/osm_LoD1_3DCityModel/tree/main/districts).|
| [village/campus]((https://github.com/AdrianKriger/osm_LoD1_3DCityModel/tree/main/village_campus)) harvests [osm contributions](https://www.openstreetmap.org/about) via [overpass-turbo](https://wiki.openstreetmap.org/wiki/Overpass_turbo) in [GeoJSON](https://geojson.org/) format| With more substantial volumes of data;<br />[districts]((https://github.com/AdrianKriger/osm_LoD1_3DCityModel/tree/main/districts)) extracts the necessary building outlines from the [osm.pbf format](https://wiki.openstreetmap.org/wiki/PBF_Format) (Protocolbuffer Binary Format) <!--via [Pyrosm](https://pyrosm.readthedocs.io/en/latest/)|-->

<table>
  <tr>
    <th align="center"><a href="https://github.com/AdrianKriger/geo3D/blob/main/osm_LoD1_3DCityModel-walkthrough.ipynb">osm_LoD1_3DCityModel </a></th>
    <th align="center"><a href="https://github.com/AdrianKriger/geo3D/blob/main/interactiveOnly.ipynb">InteractiveOnly</a> </th>
  </tr>
  <tr>
    <td align="center"> If you need a topologically correct LoD1 3D City Model <br> please choose <br> <a href="https://github.com/AdrianKriger/geo3D/blob/main/osm_LoD1_3DCityModel walkthrough.ipynb">osm_LoD1_3DCityModel-walkthrough</a> <br> followed by <br> <a href="https://github.com/AdrianKriger/geo3D/blob/main/CityJSONspatialDataScience.ipynb">CityJSONspatialDataScience.ipynb </a> </td>
    <td align="center"> Please choose <a href="https://github.com/AdrianKriger/geo3D/blob/main/interactiveOnly.ipynb">InteractiveOnly</a> <br> if you do not need a LoD1 3D City Model.     </td>
  </tr>
  <tr>
    <td align="center"> <strong>Product:</strong> An [Open Geospatial Consortium (OGC)](https://www.ogc.org/standard/citygml/) standard LoD1 3D model adhering to the [International Standards Organization (ISO19107)](https://www.iso.org/standard/66175.html) rules. <br><br> These siimulation-ready models provide insights for urban planning, resource allocation, and community development with such use cases as estimating noise propogation, energy demand and simulating wind comfort factor. </td>
    <td align="center"> <strong>Product:</strong> A 3D HTML-based visualization for user interaction, navigation, and sharing, promoting community engagement and understanding </td>
 </tr>
</table>


<!--  Table of contents
{: .no_toc .text-delta }

<!-- |<td colspan=3><b>The reason for this is</b></td> -->
<!-- ||<b>The reason for this is</b>|| -->
