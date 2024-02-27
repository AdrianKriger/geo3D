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

<table>
  <tr>
    <th align="center"><a href="https://github.com/AdrianKriger/geo3D/blob/main/osm_LoD1_3DCityModel-walkthrough.ipynb">osm_LoD1_3DCityModel </a> </th>
    <th align="center"><a href="https://github.com/AdrianKriger/geo3D/blob/main/interactiveOnly.ipynb">InteractiveOnly </a> </th>
  </tr>
  <tr>
    <td align="center"> This workflow creates a LoD1 3D City Model <br>(buildings and terrain) from <br> <a href="https://www.openstreetmap.org/about">OpenStreetMap (osm) contributions</a> <br> with elevation from a Digital Elevation Model (DEM)<sup>*</sup> </td>
    <td align="center"> This strategy creates a basic 3D visualization from <br> <a href="https://www.openstreetmap.org/about">OpenStreetMap (osm) contributions</a> <br> </td>
  </tr>
  <tr>
    <td align="center"> raster DEM </td>
    <td align="center"> nothing is necessary </td>
  </tr>
</table>

<sup>* ***for the purposes of [geo3D](https://github.com/AdrianKriger/geo3D) a DEM is a bare earth raster grid; free of man-made and natural features.***

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
