---
layout: default
title: What do I need?
parent: geo3D
nav_order: 3
---

# What do I need?
{: .no_toc }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc} 

---

## Ingredients
 
Depending of the processing strategy the needs are slightly different. 

<!--| [osm_LoD1_3DCityModel](https://github.com/AdrianKriger/geo3D/blob/main/osm_LoD1_3DCityModel-walkthrough.ipynb) | [interactiveOnly](https://github.com/AdrianKriger/geo3D/blob/main/interactiveOnly.ipynb) |
| :-----: | :-----: |
|With the [osm_LoD1_3DCityModel](https://github.com/AdrianKriger/geo3D/blob/main/osm_LoD1_3DCityModel-walkthrough.ipynb) a DEM nothing <br />more than the raster DEM is necessary. [osm_LoD1_3DModel](https://github.com/AdrianKriger/osm_LoD1_3DCityModel) will call [overpass-turbo](https://wiki.openstreetmap.org/wiki/Overpass_turbo) for the [osm contributions](https://www.openstreetmap.org/about)| [osm_LoD1_3DModel](https://github.com/AdrianKriger/osm_LoD1_3DCityModel) will access [osm contributions](https://www.openstreetmap.org/about) through [Pyrosm](https://pyrosm.readthedocs.io/en/latest/index.html).<br /><br />Due to the substantial amounts of data in the osm.pbf extract; [districts]((https://github.com/AdrianKriger/osm_LoD1_3DCityModel/tree/main/districts)) requires a bit more.<br /><br />[osm_LoD1_3DModel](https://github.com/AdrianKriger/osm_LoD1_3DCityModel) uses [osmconvert](https://wiki.openstreetmap.org/wiki/Osmconvert) to make the osm.pbf more manageable. <br />It does this through selecting only the data from a specific area. <br /><br />[osm.poly](https://wiki.openstreetmap.org/wiki/Osmosis/Polygon_Filter_File_Format) files, that cover [various regions around the world](https://github.com/JamesChevalier/cities), are available for this very purpose.|
|raster DEM | raster DEM, [osmconvert](https://wiki.openstreetmap.org/wiki/Osmconvert) and a [osm.poly](https://wiki.openstreetmap.org/wiki/Osmosis/Polygon_Filter_File_Format) file for a [region of choice](https://github.com/JamesChevalier/cities)|-->

<table>
  <tr>
    <th align="center"><a href="https://github.com/AdrianKriger/geo3D/blob/main/osm_LoD1_3DCityModel-walkthrough.ipynb">osm_LoD1_3DCityModel </a> </th>
    <th align="center"><a href="https://github.com/AdrianKriger/geo3D/blob/main/interactiveOnly.ipynb">InteractiveOnly </a> </th>
  </tr>
  <tr>
    <td align="center"> With the <a href="https://github.com/AdrianKriger/geo3D/blob/main/osm_LoD1_3DCityModel-walkthrough.ipynb">osm_LoD1_3DCityModel</a> <br /> a raster DEM is necessary. <a href="https://github.com/AdrianKriger/osm_LoD1_3DCityModel">osm_LoD1_3DModel</a> will call <a href="https://wiki.openstreetmap.org/wiki/Overpass_turbo">overpass-turbo</a> for the <a href="https://www.openstreetmap.org/about">osm contributions</a></td>
    <td style="background: repeating-linear-gradient(-45deg, transparent, transparent 5px, rgba(0,0,0,0.2) 5px, rgba(0,0,0,0.2) 10px); text-align: center; border: 1px solid black;"> <a href="https://github.com/AdrianKriger/osm_LoD1_3DCityModel">osm_LoD1_3DModel</a> will access <a href="https://www.openstreetmap.org/about">osm contributions</a> through <a href="https://pyrosm.readthedocs.io/en/latest/index.html">Pyrosm</a><br /><br />Due to the substantial amounts of data in the osm.pbf extract; <a href="https://github.com/AdrianKriger/osm_LoD1_3DCityModel/tree/main/districts">Districts</a> requires a bit more.<br /><br /> <a href="https://github.com/AdrianKriger/osm_LoD1_3DCityModel">osm_LoD1_3DModel</a> uses <a href="https://wiki.openstreetmap.org/wiki/Osmconvert">osmconvert</a> to make the osm.pbf more manageable. <br />It does this through selecting only the data from a specific area. <br /><br /> <a href="https://wiki.openstreetmap.org/wiki/Osmosis/Polygon_Filter_File_Format">osm.poly</a> files, that cover <a href="https://github.com/JamesChevalier/cities">various regions around the world</a>, are available for this very purpose.</td>
  </tr>
  <tr>
    <td align="center">raster DEM</td>
    <tdalign="center"> nothing is necessary </td>
  </tr>
</table>

## Folder Structure

The recommended folder structure is:

```
project
│   osm_LoD1_3DCityModel.ipynb
│   param.json
|   interactive.ipynb
|
└───raster
│   │   dem.tif
|
└───data
│   │   
│      
└───result

```
