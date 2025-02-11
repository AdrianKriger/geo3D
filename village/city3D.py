# -*- coding: utf-8 -*-
# env/osm3D_v4
#########################
# helper functions to create LoD1 3D City Model from volunteered public data (OpenStreetMap) with elevation via a raster DEM.

# author: arkriger - 2023-2025
# github: https://github.com/AdrianKriger/geo3D

# script credit:
#    - building height from osm building:level: https://github.com/ualsg/hdb3d-code/blob/master/hdb2d.py - Filip Biljecki <filip@nus.edu.sg>
#    - extruder: https://github.com/cityjson/misc-example-code/blob/master/extruder/extruder.py - Hugo Ledoux <h.ledoux@tudelft.nl>

# additional thanks:
#    - cityjson community: https://github.com/cityjson
#########################

import json
import fiona
import copy

import numpy as np
import pandas as pd

import shapely.geometry as sg
from shapely.geometry import Point, LineString, Polygon, MultiPolygon, LinearRing, shape, mapping
from shapely.ops import snap
from shapely.ops import transform

import pyproj

from openlocationcode import openlocationcode as olc

from cjio import cityjson, geom_help

dps = 3

# #-- calculate building height and write to geojson
# def process_geometry(geometry):
#     """Ensure valid polygon geometry for buildings."""
#     if geometry.geom_type == 'LineString':
#         return Polygon(geometry)
#     elif geometry.geom_type == 'MultiPolygon':
#         return Polygon(list(geometry.geoms)[0])
#     return geometry

# def extract_address(row):
#     # """Extract and format address components from OSM data."""
#     # address_keys = [
#     #     'addr:housename', 'addr:flats', 'addr:housenumber', 'addr:street',
#     #     'addr:suburb', 'addr:postcode', 'addr:city', 'addr:province'
#     # ]
#     # return " ".join([row[col] for col in address_keys if col in columns and row[col] is not None])

#     """Extract and format address components from the 'tags' dictionary in OSM data."""
#     address_keys = [
#         'addr:housename', 'addr:flats', 'addr:housenumber', 'addr:street',
#         'addr:suburb', 'addr:postcode', 'addr:city', 'addr:province'
#     ]
    
#     address_parts = [row.get('tags', {}).get(key) for key in address_keys if row.get('tags', {}).get(key) is not None]
    
#     return " ".join(address_parts) if address_parts else None  # Return None if empty

# def calculate_building_heights(row, storeyheight=2.8):
#     """Compute ground, building, and roof heights based on building type."""
#     ground_height = round(row.get("mean"), 2)  # Default 0 to prevent errors
#     levels = float(row.get('tags', {}).get('building:levels', 1)) * storeyheight
    
    
#     if row['building'] == 'bridge':
#         #min_height = float(row.get('min_height', 0)) if 'min_height' in row else levels
#         min_height = ( 
#             float(row.get('tags', {}).get('min_height')) if row.get('tags', {}).get('min_height') is not None else float(row.get('tags', {}).get('building:min_level', 0)) * storeyheight)

#         return {
#             'building:levels': row.get('tags', {}).get('building:levels'),
#             'ground_height': ground_height,
#             'bottom_bridge_height': round(min_height + ground_height, 2),
#             'building_height': round(levels, 2),
#             'roof_height': round(levels + ground_height, 2)

#         }
#     elif row.get('building') == 'roof':
#         return {
#             'building:levels': row.get('tags', {}).get('building:levels'),
#             'ground_height': ground_height,
#             'bottom_roof_height': round(levels + ground_height, 2),
#             'roof_height': round(levels + ground_height + 1.5, 2)
#         }
#     else:
#         return {
#             'building:levels': row.get('tags', {}).get('building:levels'),
#             'ground_height': ground_height,
#             'building_height': round(levels + 1.3, 2),
#             'roof_height': round(levels + 1.3 + ground_height, 2)
#         }

# def write_geojson(ts, jparams, epsg):
#     """Process buildings and write results to GeoJSON."""
#     storeyheight = 2.8
#     #columns = ts.columns
#     footprints = {"type": "FeatureCollection", "features": []}
    
#     for _, row in ts.iterrows():
#         if row.geometry.geom_type == 'LineString' and len(row.geometry.coords) < 3:
#             continue  # Skip invalid geometries
        
#         f = {"type": "Feature", "properties": {"osm_id": row.id}}
#         f["properties"]["address"] = extract_address(row)
        
#         # Store 'building' attribute directly
#         f["properties"]["building"] = row.get("building")

#         # Harvest attributes from 'tags' dictionary
#         for key in [
#             'building:use', 'building:levels', 'building:flats', 'building:units',
#             'beds', 'rooms', 'residential', 'amenity', 'social_facility'
#         ]:
#             value = row.get('tags', {}).get(key)  # Extract safely from 'tags' dictionary
#             if value is not None:
#                 f["properties"][key] = value  # Store only if not None
        
#         osm_shape = process_geometry(row["geometry"])
#         f["geometry"] = mapping(osm_shape)
#         f["properties"]["footprint"] = mapping(osm_shape)
        
#         # Compute plus_code
#         p = osm_shape.representative_point()
#         f["properties"]["plus_code"] = olc.encode(p.y, p.x, 11)
        
#         # Compute height attributes
#         height_attributes = calculate_building_heights(row, storeyheight) 
#         for key, value in height_attributes.items(): 
#             if key not in f["properties"]: 
#                 f["properties"][key] = value  # Add new values only
        
#         footprints['features'].append(f)
    
#     with open(jparams['osm_bldings'], 'w') as outfile:
#         json.dump(footprints, outfile, indent=2)



# import json
# import openlocationcode as olc
# from shapely.geometry import shape, Polygon, mapping

def process_geometry(geometry):
    """Ensure valid polygon geometry for buildings."""
    if geometry.geom_type == 'LineString':
        return Polygon(geometry)
    elif geometry.geom_type == 'MultiPolygon':
        return Polygon(list(geometry.geoms)[0])
    return geometry

def extract_address(row):
    """Extract and format address components from the 'tags' dictionary in OSM data."""
    address_keys = [
        'addr:housename', 'addr:flats', 'addr:housenumber', 'addr:street',
        'addr:suburb', 'addr:postcode', 'addr:city', 'addr:province'
    ]
    address_parts = [row.get('tags', {}).get(key) for key in address_keys if row.get('tags', {}).get(key) is not None]
    return " ".join(address_parts) if address_parts else None  # Return None if empty

def calculate_building_heights(row, storeyheight=2.8):
    """Compute ground, building, and roof heights based on building type."""
    ground_height = round(row.get("mean", 0), 2)  
    levels = float(row.get('tags', {}).get('building:levels', 1)) * storeyheight

    if row.get('building') == 'bridge':
        min_height = (
            float(row.get('tags', {}).get('min_height')) 
            if row.get('tags', {}).get('min_height') is not None 
            else float(row.get('tags', {}).get('building:min_level', 0)) * storeyheight
        )
        return {
            'ground_height': ground_height,
            'bottom_bridge_height': round(min_height + ground_height, 2),
            'building_height': round(levels, 2),
            'roof_height': round(levels + ground_height, 2)
        }
    elif row.get('building') == 'roof':
        return {
            'ground_height': ground_height,
            'bottom_roof_height': round(levels + ground_height, 2),
            'roof_height': round(levels + ground_height + 1.5, 2)
        }
    else:
        return {
            'ground_height': ground_height,
            'building_height': round(levels + 1.3, 2),
            'roof_height': round(levels + 1.3 + ground_height, 2)
        }

def write_geojson(ts, jparams, epsg):
    """Process buildings and write results to GeoJSON (handles both standard and district processing)."""
    storeyheight = 2.8
    footprints = {"type": "FeatureCollection", "features": []}

    for _, row in ts.iterrows():
        if row.geometry.geom_type == 'LineString' and len(row.geometry.coords) < 3:
            continue  # Skip invalid geometries
        
        #- Ensure only valid buildings with 'building:levels' are processed 
        if row.get('type') == 'node' or row.get('tags') is None or 'building:levels' not in row.get('tags'): 
            continue  # Skip rows that don't meet the condition
        
        f = {"type": "Feature", "properties": {"osm_id": row.id}}
        f["properties"]["address"] = extract_address(row)
        f["properties"]["building"] = row.get("building")

        # Harvest attributes from 'tags' dictionary
        for key in [
            'building:use', 'building:levels', 'building:flats', 'building:units',
            'beds', 'rooms', 'residential', 'amenity', 'social_facility'
        ]:
            value = row.get('tags', {}).get(key)  
            if value is not None:
                f["properties"][key] = value 
        
        osm_shape = process_geometry(row["geometry"])
        f["geometry"] = mapping(osm_shape)
        f["properties"]["footprint"] = mapping(osm_shape)
        
        # Compute plus_code
        p = osm_shape.representative_point()
        f["properties"]["plus_code"] = olc.encode(p.y, p.x, 11)

        # Compute height attributes
        height_attributes = calculate_building_heights(row, storeyheight) 
        for key, value in height_attributes.items(): 
            if key not in f["properties"]: 
                f["properties"][key] = value  # Add new values only
        
        footprints['features'].append(f)

    # Store the data as GeoJSON
    with open(jparams['osm_bldings'], 'w') as outfile:
        json.dump(footprints, outfile, indent=2)
    print(f"GeoJSON saved successfully to {jparams['osm_bldings']}")



def getBldVertices(dis, gt_forward, rb): #
    """
    retrieve vertices from building footprints ~ without duplicates 
    - these vertices already have a z attribute
    """  
    all_coords = []
    min_zbld = []
    dps = 3
    segs = set()
    
    for ids, row in dis.iterrows():
        oring = list(row.geometry.exterior.coords)
        
        if row.geometry.exterior.is_ccw == False:
            #-- to get proper orientation of the normals
            oring.reverse()
        
        coords_rounded = [(round(x, dps), round(y, dps), round(float(rasterQuery2(x, y, gt_forward, rb)), 2)) for x, y in oring]
        all_coords.extend(coords_rounded)
        zbld = [z for x, y, z in coords_rounded]
        min_zbld.append(min(zbld))
        
        segs.update({(x1, y1, x2, y2) if (x1 < x2) else (x2, y2, x1, y1) for (x1, y1, z1), (x2, y2, z2) in zip(coords_rounded[:-1], coords_rounded[1:])})
        
        for interior in row.geometry.interiors:
            iring = list(interior.coords)
            
            if interior.is_ccw == True:
                #-- to get proper orientation of the normals
                iring.reverse() 
            
            coords_rounded = [(round(x, dps), round(y, dps), round(float(rasterQuery2(x, y, gt_forward, rb)), 2)) for x, y in iring]
            all_coords.extend(coords_rounded)
            
            segs.update({(x1, y1, x2, y2) if (x1 < x2) else (x2, y2, x1, y1) for (x1, y1, z1), (x2, y2, z2) in zip(coords_rounded[:-1], coords_rounded[1:])})
    
    c = pd.DataFrame.from_dict({"coords": list(segs)}).groupby("coords").size().reset_index(name="count")
    
    ac = pd.DataFrame(all_coords, 
                      columns=["x", "y", "z"]).sort_values(by="z", ascending=False).drop_duplicates(subset=["x", "y"]).reset_index(drop=True)
        
    return ac, c, min_zbld

def rasterQuery2(mx, my, gt_forward, rb):
    
    px = int((mx - gt_forward[0]) / gt_forward[1])
    py = int((my - gt_forward[3]) / gt_forward[5])

    intval = rb.ReadAsArray(px, py, 1, 1)

    return intval[0][0]

def getAOIVertices(aoi, gt_forward, rb): 
    """
    retrieve vertices from aoi ~ without duplicates 
    - these vertices are assigned a z attribute
    """   
    aoi_coords = []
    dps = 3
    segs = set()
    
    for ids, row in aoi.iterrows():
        oring = list(row.geometry.exterior.coords)
        
        if row.geometry.exterior.is_ccw == False:
            #-- to get proper orientation of the normals
            oring.reverse()
        
        coords_rounded = [(round(x, dps), round(y, dps), round(float(rasterQuery2(x, y, gt_forward, rb)), 2)) for x, y in oring]
        aoi_coords.extend(coords_rounded)
        
        segs.update({(x1, y1, x2, y2) if (x1 < x2) else (x2, y2, x1, y1) for (x1, y1, z1), (x2, y2, z2) in zip(coords_rounded[:-1], coords_rounded[1:])})
        
        for interior in row.geometry.interiors:
            iring = list(interior.coords)
            
            if interior.is_ccw == True:
                #-- to get proper orientation of the normals
                iring.reverse() 
            
            coords_rounded = [(round(x, dps), round(y, dps), round(float(rasterQuery2(x, y, gt_forward, rb)), 2)) for x, y in iring]
            aoi_coords.extend(coords_rounded)
            
            segs.update({(x1, y1, x2, y2) if (x1 < x2) else (x2, y2, x1, y1) for (x1, y1, z1), (x2, y2, z2) in zip(coords_rounded[:-1], coords_rounded[1:])})
    
    ca = pd.DataFrame.from_dict({"coords": list(segs)}).groupby("coords").size().reset_index(name="count")
    
    acoi = pd.DataFrame(aoi_coords, 
                      columns=["x", "y", "z"]).sort_values(by="z", ascending=False).drop_duplicates(subset=["x", "y"]).reset_index(drop=True)
    
    return acoi, ca

def concatCoords(gdf, ac):
    df2 = pd.concat([gdf, ac])
    
    return df2

def createSgmts(ac, c, gdf, idx):
    """
    create a segment list for Triangle
    - indices of vertices [from, to]
    """
    
    l = len(gdf) #- 1
    lr = 0
    idx01 = []
    
    for i, row in c.iterrows():
        frx, fry = row.coords[0], row.coords[1]
        tox, toy = row.coords[2], row.coords[3]

        [index_f] = (ac[(ac['x'] == frx) & (ac['y'] == fry)].index.values)
        [index_t] = (ac[(ac['x'] == tox) & (ac['y'] == toy)].index.values)
        idx.append([l + index_f, l + index_t])
        idx01.append([lr + index_f, lr + index_t])
    
    return idx, idx01


# -- create CityJSON
def doVcBndGeomRd(lsgeom, lsattributes, extent, minz, maxz, TerrainT, pts, acoi, jparams, min_zbld, result): 
    
    #-- create the JSON data structure for the City Model
    cm = {}
    cm["type"] = "CityJSON"
    cm["version"] = "1.1"
    #cm["transform"] = {
        #"scale": [0.0, 0.0, 0.0],
        #"translate": [1.0, 1.0, 1.0]
    #},
    cm["CityObjects"] = {}
    cm["vertices"] = []
    #-- Metadata is added manually
    cm["metadata"] = {
    "title": jparams['cjsn_title'],
    "referenceDate": jparams['cjsn_referenceDate'],
    #"dataSource": jparams['cjsn_source'],
    #"geographicLocation": jparams['cjsn_Locatn'],
    "referenceSystem": jparams['cjsn_referenceSystem'],
    "geographicalExtent": [
        extent[0],
        extent[1],
        minz ,
        extent[1],
        extent[1],
        maxz
      ],
    "datasetPointOfContact": {
        "contactName": jparams['cjsn_contactName'],
        "emailAddress": jparams['cjsn_emailAddress'],
        "contactType": jparams['cjsn_contactType'],
        "website": jparams['cjsn_website']
        },
    "+metadata-extended": {
        "lineage":
            [{"featureIDs": ["TINRelief"],
             "source": [
                 {
                     "description": jparams['cjsn_+meta-description'],
                     "sourceSpatialResolution": jparams['cjsn_+meta-sourceSpatialResolution'],
                     "sourceReferenceSystem": jparams['cjsn_+meta-sourceReferenceSystem'],
                     "sourceCitation":jparams['cjsn_+meta-sourceCitation'],
                     }],
             "processStep": {
                 "description" : "Processing of raster DEM using osm_LoD1_3DCityModel workflow",
                 "processor": {
                     "contactName": jparams['cjsn_contactName'],
                     "contactType": jparams['cjsn_contactType'],
                     "website": jparams['cjsn_website']
                     }
                 }
            },
            {"featureIDs": ["Building", "Road"],
             "source": [
                 {
                     "description": "OpenStreetMap contributors",
                     "sourceReferenceSystem": "urn:ogc:def:crs:EPSG:4326",
                     "sourceCitation": "https://www.openstreetmap.org",
                 }],
             "processStep": {
                 "description" : "Processing of building vector contributions using osm_LoD1_3DCityModel workflow",
                 "processor": {
                     "contactName": jparams['cjsn_contactName'],
                     "contactType": jparams['cjsn_contactType'],
                     "website": "https://github.com/AdrianKriger/osm_LoD1_3DCityModel"
                     }
                 }
            }]
        }
    #"metadataStandard": jparams['metaStan'],
    #"metadataStandardVersion": jparams['metaStanV']
    }
      ##-- do terrain
    add_terrain_v(pts, cm)
    grd = {}
    grd['type'] = 'TINRelief'
    grd['geometry'] = [] #-- a cityobject can have >1 
      #-- the geometry
    g = {} 
    g['type'] = 'CompositeSurface'
    g['lod'] = 1
    g['boundaries'] = []
    allsurfaces = [] #-- list of surfaces
    add_terrain_b(TerrainT, allsurfaces)
    g['boundaries'] = allsurfaces
      #-- add the geom 
    grd['geometry'].append(g)
      #-- insert the terrain as one new city object
    cm['CityObjects']['terrain01'] = grd
    

    count = 0
     #-- then buildings
    for (i, geom) in enumerate(lsgeom):

        poly = list(result[lsattributes[i]['osm_id']].values())
        footprint = geom
        footprint = sg.polygon.orient(footprint, 1)

        #-- one building
        oneb = {}
        oneb['type'] = 'Building'
        oneb['attributes'] = {}
        for (k, v)in list(lsattributes[i].items()):
            if v is None:
                del lsattributes[i][k]
        for a in lsattributes[i]:
            oneb['attributes'][a] = lsattributes[i][a]                   
        oneb['geometry'] = [] #-- a cityobject can have > 1
        
        #-- the geometry
        g = {} 
        g['type'] = 'Solid'
        g['lod'] = 1
        allsurfaces = [] #-- list of surfaces forming the shell of the solid
        #-- exterior ring of each footprint
        oring = list(footprint.exterior.coords)
        oring.pop() #-- remove last point since first==last
        
        if footprint.exterior.is_ccw == False:
            #-- to get proper orientation of the normals
            oring.reverse()
        
        if lsattributes[i]['building'] == 'bridge':
            edges = [[ele for ele in sub if ele <= lsattributes[i]['roof_height']] for sub in poly]
            extrude_walls(oring, lsattributes[i]['roof_height'], lsattributes[i]['bottom_bridge_height'], 
                          allsurfaces, cm, edges)
            count = count + 1

        if lsattributes[i]['building'] == 'roof':
            edges = [[ele for ele in sub if ele <= lsattributes[i]['roof_height']] for sub in poly]
            extrude_walls(oring, lsattributes[i]['roof_height'], lsattributes[i]['bottom_roof_height'], 
                          allsurfaces, cm, edges)
            count = count + 1

        if lsattributes[i]['building'] != 'bridge' and lsattributes[i]['building'] != 'roof':
            new_edges = [[ele for ele in sub if ele <= lsattributes[i]['roof_height']] for sub in poly]
            new_edges = [[min_zbld[i-count]] + sub_list for sub_list in new_edges]
            extrude_walls(oring, lsattributes[i]['roof_height'], min_zbld[i-count], 
                          allsurfaces, cm, new_edges)
       
        #-- interior rings of each footprint
        irings = []
        interiors = list(footprint.interiors)
        for each in interiors:
            iring = list(each.coords)
            iring.pop() #-- remove last point since first==last
            
            if each.is_ccw == True:
                #-- to get proper orientation of the normals
                iring.reverse() 
            
            irings.append(iring)
            extrude_int_walls(iring, lsattributes[i]['roof_height'], min_zbld[i-count], allsurfaces, cm)
            
        #-- top-bottom surfaces
        if lsattributes[i]['building'] == 'bridge':
            extrude_roof_ground(oring, irings, lsattributes[i]['roof_height'], 
                                False, allsurfaces, cm)
            extrude_roof_ground(oring, irings, lsattributes[i]['bottom_bridge_height'], 
                                True, allsurfaces, cm)
        if lsattributes[i]['building'] == 'roof':
            extrude_roof_ground(oring, irings, lsattributes[i]['roof_height'], 
                                False, allsurfaces, cm)
            extrude_roof_ground(oring, irings, lsattributes[i]['bottom_roof_height'], 
                                True, allsurfaces, cm)
        if lsattributes[i]['building'] != 'bridge' and lsattributes[i]['building'] != 'roof':
            extrude_roof_ground(oring, irings, lsattributes[i]['roof_height'], 
                            False, allsurfaces, cm)
            extrude_roof_ground(oring, irings, min_zbld[i-count], True, allsurfaces, cm)

        #-- add the extruded geometry to the geometry
        g['boundaries'] = []
        g['boundaries'].append(allsurfaces)
        
        #-- add the geom to the building 
        oneb['geometry'].append(g)
        #-- insert the building as one new city object
        cm['CityObjects'][lsattributes[i]['osm_id']] = oneb

    return cm

def add_terrain_v(pts, cm):
    for p in pts:
        cm['vertices'].append([p[0], p[1], p[2]])
    
def add_terrain_b(Terr, allsurfaces):
    for i in Terr:
        allsurfaces.append([[i[0], i[1], i[2]]]) 
        
def extrude_roof_ground(orng, irngs, height, reverse, allsurfaces, cm):
    oring = copy.deepcopy(orng)
    irings = copy.deepcopy(irngs)
    #irings2 = []
    if reverse == True:
        oring.reverse()
        for each in irings:
            each.reverse()
    for (i, pt) in enumerate(oring):
        cm['vertices'].append([round(pt[0], dps), round(pt[1], dps), height])
        oring[i] = (len(cm['vertices']) - 1)
    for (i, iring) in enumerate(irings):
        for (j, pt) in enumerate(iring):
            cm['vertices'].append([round(pt[0], dps), round(pt[1], dps), height])
            irings[i][j] = (len(cm['vertices']) - 1)
    output = []
    output.append(oring)
    for each in irings:
        output.append(each)
    allsurfaces.append(output)
    
def extrude_walls(ring, height, ground, allsurfaces, cm, edges):  
    #-- each edge become a wall, ie a rectangle
    for (j, v) in enumerate(ring[:-1]):
        #- if iether the left or right vertex has more than 2 heights [grnd and roof] incident:
        if len(edges[j]) > 2 or len(edges[j+1]) > 2:
            cm['vertices'].append([round(ring[j][0], dps), round(ring[j][1], dps), edges[j][0]])
            cm['vertices'].append([round(ring[j+1][0], dps), round(ring[j+1][1], dps), edges[j+1][0]])
            c = 0
            #- traverse up [grnd-roof]:
            for i, o in enumerate(edges[j+1][1:]):
                cm['vertices'].append([round(ring[j+1][0], dps), round(ring[j+1][1], dps), o])
                c = c + 1
            #- traverse down [roof-grnd]:
            for i in edges[j][::-1][:-1]:
                cm['vertices'].append([round(ring[j][0], dps), round(ring[j][1], dps), i])
                c = c + 1
            t = len(cm['vertices'])
            c = c + 2
            b = c
            l = []
            for i in range(c):
                l.append(t-b)
                b = b - 1 
            allsurfaces.append([l])

        #- if iether the left and right vertex has only 2 heights [grnd and roof] incident: 
        if len(edges[j]) == 2 and len(edges[j+1]) == 2:
            cm['vertices'].append([round(ring[j][0], dps),   round(ring[j][1], dps),   edges[j][0]])
            cm['vertices'].append([round(ring[j+1][0], dps), round(ring[j+1][1], dps), edges[j+1][0]])
            cm['vertices'].append([round(ring[j+1][0], dps), round(ring[j+1][1], dps), edges[j+1][1]])
            cm['vertices'].append([round(ring[j][0], dps),   round(ring[j][1], dps),   edges[j][1]])
            t = len(cm['vertices'])
            allsurfaces.append([[t-4, t-3, t-2, t-1]])
    
    #- last edge polygon
    if len(edges[-1]) == 2 and len(edges[0]) == 2:
        cm['vertices'].append([round(ring[-1][0], dps),  round(ring[-1][1], dps), edges[-1][0]]) 
        cm['vertices'].append([round(ring[0][0], dps), round(ring[0][1], dps), edges[0][0]])
        cm['vertices'].append([round(ring[0][0], dps),  round(ring[0][1], dps),  edges[0][1]])
        cm['vertices'].append([round(ring[-1][0], dps), round(ring[-1][1], dps), edges[-1][1]])
        t = len(cm['vertices'])
        allsurfaces.append([[t-4, t-3, t-2, t-1]])
        
    #- last edge polygon   
    if len(edges[-1]) > 2 or len(edges[0]) > 2:
        c = 0
        cm['vertices'].append([round(ring[-1][0], dps),   round(ring[-1][1], dps),   edges[-1][0]])
        cm['vertices'].append([round(ring[0][0], dps), round(ring[0][1], dps), edges[0][0]])
        for i, o in enumerate(edges[0][1:]):
            cm['vertices'].append([round(ring[0][0], dps), round(ring[0][1], dps), o])
            c = c + 1
        for i in edges[-1][::-1][:-1]:
            cm['vertices'].append([round(ring[-1][0], dps),   round(ring[-1][1], dps),   i])
            c = c + 1
        t = len(cm['vertices'])
        c = c + 2
        b = c
        l = []
        for i in range(c): 
            l.append(t-b)
            b = b - 1 
        allsurfaces.append([l])
               
def extrude_int_walls(ring, height, ground, allsurfaces, cm):
    #-- each edge become a wall, ie a rectangle
    for (j, v) in enumerate(ring[:-1]):
        l = []
        cm['vertices'].append([round(ring[j][0], dps),   round(ring[j][1], dps),   ground])
        #values.append(0)
        cm['vertices'].append([round(ring[j+1][0], dps), round(ring[j+1][1], dps), ground])
        #values.append(0)
        cm['vertices'].append([round(ring[j+1][0], dps), round(ring[j+1][1], dps), height])
        cm['vertices'].append([round(ring[j][0], dps), round(ring[j][1], dps), height])
        t = len(cm['vertices'])
        allsurfaces.append([[t-4, t-3, t-2, t-1]])    
    #-- last-first edge
    #l = []
    cm['vertices'].append([round(ring[-1][0], dps), round(ring[-1][1], dps), ground])
    #values.append(0)
    cm['vertices'].append([round(ring[0][0], dps), round(ring[0][1], dps), ground])
    cm['vertices'].append([round(ring[0][0], dps), round(ring[0][1], dps), height])
    #values.append(0)
    cm['vertices'].append([round(ring[-1][0], dps), round(ring[-1][1], dps), height])
    t = len(cm['vertices'])
    allsurfaces.append([[t-4, t-3, t-2, t-1]])
    
def output_cityjson(extent, minz, maxz, TerrainT, pts, jparams, min_zbld, acoi, result):
    """
    basic function to produce LoD1 City Model
    - buildings and terrain
    """
     ##- open buildings ---fiona object
    c = fiona.open(jparams['osm_bldings'])
    lsgeom = [] #-- list of the geometries
    lsattributes = [] #-- list of the attributes
    for each in c:
        lsgeom.append(shape(each['geometry'])) #-- geom are casted to Fiona's 
        lsattributes.append(each['properties'])
               
    #- 3D Model
    cm = doVcBndGeomRd(lsgeom, lsattributes, extent, minz, maxz, TerrainT, pts, acoi, jparams, min_zbld, result)    
    
    json_str = json.dumps(cm)#, indent=2)
    fout = open(jparams['cjsn_out'], "w")                 
    fout.write(json_str)  
    ##- close fiona object
    c.close() 
    #clean cityjson
    cm = cityjson.load(jparams['cjsn_out'])               
    cityjson.save(cm, jparams['cjsn_solid'])   
    
    # #- 2.75D surface
    # cm = doVcBndGeomOBJ(lsgeom, lsattributes, extent, minz, maxz, TerrainT, pts, acoi, jparams, min_zbld, result)    
    
    # json_str = json.dumps(cm)#, indent=2)
    # fout = open(jparams['cjsn_out'], "w")                 
    # fout.write(json_str)  
    # ##- close fiona object
    # c.close() 
    # #clean 2.75D cityjson
    # cm = cityjson.load(jparams['cjsn_out'])               
    # cityjson.save(cm, jparams['cjsn_2_75D'])  
    # #- write obj
    # write275obj(jparams)
    

# def calc_Bldheight(gj):

#     storeyheight = 2.8
    
#     #-- iterate through the list of buildings and create GeoJSON features rich in attributes
#     footprints = {
#         "type": "FeatureCollection",
#         "features": []
#         }
    
#     for i in gj['features']:
#         f = {
#         "type" : "Feature"
#         }
#         # at a minimum we only want building:levels tagged
#         if i['properties']['type'] != 'node' and 'tags' in i['properties'] \
#         and 'building:levels' in i['properties']['tags'] is not None:
            
#             f["properties"] = {}
#             for p in i["properties"]:             
#             #-- store OSM attributes and prefix them with osm_
#                 f['properties']['osm_id'] = i['properties']['id']
#                 f["properties"]["osm_%s" % p] = i["properties"][p]
#                 if 'amenities' in i['properties']:
#                     f['properties']['osm_tags']['amenity'] = i['properties']['amenities']
#                 osm_shape = shape(i["geometry"])
#             #-- a few buildings are not polygons, rather linestrings. This converts them to polygons
#             #-- rare, but if not done it breaks the code later
#             #if osm_shape.type == 'LineString':
#             if osm_shape.geom_type == 'LineString':
#                 osm_shape = Polygon(osm_shape)
#             #-- and multipolygons must be accounted for
#             #elif osm_shape.type == 'MultiPolygon':
#             elif osm_shape.geom_type == 'MultiPolygon':
#                 #osm_shape = Polygon(osm_shape.geoms[0])
#                 polys = list(osm_shape.geoms) 
#                 #for poly in list(osm_shape.geoms):
#                 for poly in polys:
#                     osm_shape = Polygon(poly)#[0])
#             #-- convert the shapely object to geojson
#             f["geometry"] = mapping(osm_shape)
#             f["properties"]['footprint'] = mapping(osm_shape)
    
#             #-- calculate the height and store it as an attribute
#             f["properties"]['height'] = float(i["properties"]['tags']['building:levels']) * storeyheight + 1.3 
#             #-- plus code
#             p = osm_shape.representative_point()
#             f["properties"]["plus_code"] = olc.encode(p.y, p.x, 11)
                
#             footprints['features'].append(f)
    
    
#     #-- store the data as GeoJSON
#     with open('./data/fp_j.geojson', 'w') as outfile:
#         json.dump(footprints, outfile)        


def calc_Bldheight(data, is_geojson=True):
    """Calculate building height and write to GeoJSON from either a GeoJSON dictionary or a DataFrame."""
    
    storeyheight = 2.8  # Default storey height assumption
    footprints = {"type": "FeatureCollection", "features": []}

    # Decide whether data comes from GeoJSON (`gj`) or a DataFrame (`df2`)
    iterable = data["features"] if is_geojson else data.iterrows()

    for item in iterable:
        f = {"type": "Feature", "properties": {}}

        # Extract row based on input type
        row = item if is_geojson else item[1]
        properties = row["properties"] if is_geojson else row.to_dict()
        tags = properties.get("tags", {}) if is_geojson else row  # Tags dictionary for GeoJSON

        # Skip nodes (only process buildings with 'building:levels')
        if properties.get("type") == "node" or "building:levels" not in tags:
            continue

        # Store OSM attributes
        f["properties"]["osm_id"] = properties.get("id") if is_geojson else row["id"]

        # Extract address components
        address_keys = [
            'name', 'addr:housename', 'addr:flats', 'addr:housenumber', 
            'addr:street', 'addr:suburb', 'addr:postcode', 'addr:city', 'addr:province'
        ]
        address_parts = [tags.get(k) for k in address_keys if tags.get(k) is not None]
        f["properties"]["address"] = " ".join(address_parts) if address_parts else None

        # # Extract building attributes
        # for key in [
        #     'building', 'building:use', 'building:levels', 'building:flats', 
        #     'building:units', 'beds', 'rooms', 'residential', 'amenity', 'social_facility'
        # ]:
        #     value = tags.get(key)
        #     if value is not None:
        #         f["properties"][f"osm_{key}"] = value
        
        # Harvest attributes from 'tags' dictionary
        for key in [
            'building', 'building:use', 'building:levels', 'building:flats', 'building:units',
            'beds', 'rooms', 'residential', 'amenity', 'social_facility'
        ]:
            value = tags.get(key)
            if value is not None:
                f["properties"][key] = value  # Store only if not None
        

        # Convert geometry to a valid polygon
        osm_shape = shape(row["geometry"]) if is_geojson else row["geometry"]
        if osm_shape.geom_type == 'LineString':
            osm_shape = Polygon(osm_shape)
        elif osm_shape.geom_type == 'MultiPolygon':
            polys = list(osm_shape.geoms)
            osm_shape = Polygon(polys[0])  # Take the first polygon

        f["geometry"] = mapping(osm_shape)
        f["properties"]["footprint"] = mapping(osm_shape)

        # Compute Plus Code
        p = osm_shape.representative_point()
        f["properties"]["plus_code"] = olc.encode(p.y, p.x, 11)

        # Compute building height
        levels = float(tags.get('building:levels', 1))  # Default to 1 level if missing
        f["properties"]['building_height'] = round(levels * storeyheight + 1.3, 2)

        footprints["features"].append(f)

    # Store data as GeoJSON
    with open('./data/fp_j.geojson', 'w') as outfile:
        json.dump(footprints, outfile, indent=2)
               