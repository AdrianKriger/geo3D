---
layout: default
title: What does it do?
parent: geo3D
nav_order: 1
---

# What does it do?
{: .no_toc }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## Extrusion

An extremely well documented method of producing 3D Models is through extrusion. With extrusion; 2D features are *lifted*{: style="color: red; opacity: 0.80" } from an existing surface creating a volumetric 3D object.  <code><b>geo3D</b></code> inferes the height with which to *lift*{: style="color: red; opacity: 0.80" } 2D features from osm contributions. 

The osm tag `building:level` is taken as a [proxy for the height of a building](https://wiki.openstreetmap.org/wiki/Key:building:levels). The calculation is simply `building:level * 2.8 + 1.3`. If a structure does not have a `building:level` tag no LoD1 model is created.
 &nbsp; &nbsp;
 <figure><center>
  <img src="{{site.baseurl | prepend: site.url}}/img/extrusion_tuDelft.png" alt="alt text" width="650" height="350">
  <!-- <figcaption>Fig.1 - <code><b>The osm_LoD1_3DCityModel</b></code> process. <span style="color:blue;opacity:0.8;"><em>--image TUDelft</em></span>.</figcaption> -->
  <figcaption>Fig.1 - The <code><b>osm_LoD1_3DCityModel</b></code> process. <em>--image adapted from</em><cite><a href="https://github.com/tudelft3d/3dfier"> 3dfier</a></cite> by the<cite><a href="https://3d.bk.tudelft.nl/"> 3D geoinformation research group at TU Delft</a></cite>.</figcaption>
</center></figure>
<!--
<p align="center">
  <img src="{{site.baseurl | prepend: site.url}}/img/extrusion_tuDelft.png" alt="alt text" width="650" height="350">
 </p> 
<p align="center"> 
    Fig 1. The osm_LoD1_3DCityModel process. <span style="color:blue"><em>--image TUDelft</em></span>.
</p> -->

Fig 1 illustrates the process where the osm *proxy `building:level` height*  is added to the raster DEM to create a 3D topologically connected surface ~ containing 2D polygons as 3D objects.

## `geo3D` products

<!--### Trianglated MultiSurfaces

MultiSurface outputs are the walls and rooves of buildings, along with the terrain, as a collection of connected triangles. This surface is created in the [Wavefront OBJ](https://en.wikipedia.org/wiki/Wavefront_.obj_file) format. An accompanying [material file](https://en.wikipedia.org/wiki/Wavefront_.obj_file#Material_template_library) (.mtl) to associate objects with a respective color is [available](https://github.com/AdrianKriger/osm_LoD1_3DCityModel/blob/main/village_campus/result/osm_LoD1_3DCityModel.mtl). 

<figure><center>
  <img src="{{site.baseurl | prepend: site.url}}/img/objects_horizontal_view_multisurface_tuDelft.png" alt="alt text" width="650" height="350">
  <figcaption>Fig.2 - illustrates a horizontal view of the 2.75D surface with the exterior of all features together. <em>--image </em><cite><a href="https://github.com/tudelft3d/3dfier"> 3dfier</a></cite> by the<cite><a href="https://3d.bk.tudelft.nl/"> 3D geoinformation research group at TU Delft</a></cite>.</figcaption>
</center></figure>   
-->
### Solids

Solid objects are walls, floors and rooves stored as rectangular faces. In the [CityJSON](https://www.cityjson.org/) format these are [Building](https://www.cityjson.org/specs/1.0.1/#building) [City Ojects](https://www.cityjson.org/specs/1.0.1/#cityjson-object) separate from the [TINRelief](https://www.cityjson.org/specs/1.0.1/#tinrelief).

<figure><center>
  <img src="{{site.baseurl | prepend: site.url}}/img/objects_horizontal_view_solid_tuDelft.png" alt="alt text" width="650" height="350">
  <figcaption>Fig.2 - solid Building CityObjects's connected to the terrain. <em>--image </em><cite><a href="https://github.com/tudelft3d/3dfier"> 3dfier</a></cite> by the<cite><a href="https://3d.bk.tudelft.nl/"> 3D geoinformation research group at TU Delft</a></cite>.</figcaption>
</center></figure>

### Interactive Visualisation with Spatial Data Science

A dynamic visualisation and spatial analysis is possible. [Interactive Visualization](https://adriankriger.github.io/geo3D/docs/interactive/) and [Spatial Data Science](https://adriankriger.github.io/geo3D/docs/spatial/) discusses this further.

## Is it useful?

The LoD1 City Model, while basic, offers many advantages over 2D datasets. These may be used for shadow analyses, line of sight predictions, advanced flood simulation, or more advanced quantitative evaluations such as estimating wind comfort factor and simulating noise propogation.

As the  [CityJSONspatialDataScience.ipynb](https://github.com/AdrianKriger/geo3D/blob/main/CityJSONspatialDataScience.ipynb) illustrates population estimation and the calculation of [Building Volume per Capita (Ghosh, T.; et. al.)](https://www.frontiersin.org/articles/10.3389/frsc.2020.00037/full) are also possible.

With the coming revolution in air traffic control, to accommodate newer forms of air services (delivery drones and urban air mobility), an accurate *digital* representation of the built environment will become crucial. A 3D City Model is one component for the effective air space management of the future.

Challenges do exist. Of primary concern are errors in the source data that propagate to the generated 3D model. Care must be taken to ensure the quality of both the vector building outlines and raster DEM.
