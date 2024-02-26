---
layout: default
title: Spatial Data Science
nav_order: 5
---

# Spatial Data Science
{: .no_toc }

The [Jupyter](https://eis n.wikipedia.org/wiki/Project_Jupyter#Jupyter_Notebook) environment allows for extensive customization and deep analysis through *spatial data science*.

[geo3D](https://github.com/AdrianKriger/geo3D/tree/main) illustrates a example of population estimation and the calculation of [Building Volume per Capita (Ghosh, T.; et. al.)](https://www.frontiersin.org/articles/10.3389/frsc.2020.00037/full).

While the prefered process would proceed: 

<figure><center>
  <img src="{{site.baseurl | prepend: site.url}}/img/flow1.png" style="width: 800px; height: 300px; border: 0px">
</center></figure> 

[osm_LoD1_3DCityModel](https://github.com/AdrianKriger/geo3D/blob/main/osm_LoD1_3DCityModel-walkthrough.ipynb)) -> [Spatial Data Science](https://github.com/AdrianKriger/geo3D/blob/main/CityJSONspatialDataScience.ipynb) -> Interactive visualization; an alternate does exist. 

[interactiveOnly.ipynb](https://github.com/AdrianKriger/geo3D/blob/main/interactiveOnly.ipynb) will create a basic 3D model visualisation followed by population estimation and the calculation of [Building Volume per Capita (Ghosh, T.; et. al.)](https://www.frontiersin.org/articles/10.3389/frsc.2020.00037/full).

<figure><center>
  <img src="{{site.baseurl | prepend: site.url}}/img/flow2.png" style="width: 500px; height: 250px; border: 0px">
</center></figure> 

Please consider your needs before executing the solution. We do not want to burden the [OpenStreetMap](https://www.openstreetmap.org/#map=5/-28.676/24.677) server with repeat calls for data. 
