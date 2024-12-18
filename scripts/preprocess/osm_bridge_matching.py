#!/usr/bin/env python
# coding: utf-8
import sys
import os
import pandas as pd
import ast
import geopandas as gpd
import shapely
from utils import *
from tqdm import tqdm
tqdm.pandas()

def assign_lanes_to_roads(x,rc_column,fw_column,lane_width,shoulder_width,collapsed_roads=False):
    width = x["averageWidth"]
    if (x[rc_column] == "Motorway") and (x[fw_column] == "Collapsed Dual Carriageway"):
        width = width - 2*shoulder_width
    elif x[rc_column] == "Motorway":
        width = width - shoulder_width 

    lanes = np.floor(1.0*width/lane_width).astype(int)
    if collapsed_roads is True:
        lanes = np.floor(1.0*width/lane_width).astype(int)
        if (lanes > 1) and (lanes % 2 == 0):
            lanes += 0
        else:
            lanes += 1

        if x[rc_column] == "Motorway":
            lanes = max(lanes,4)
    else:
        if (x[fw_column] == "Single Carriageway") and (lanes % 2 != 0):
            lanes = lanes - 1

    
    if lanes == 0:
        lanes += 1
    
    return lanes

def main(config):
    incoming_data_path = config['paths']['incoming_data']
    data_path = config['paths']['data']
    epsg_meters = 27700
    buffer_distance = 6.0
    lane_width = 3.65 # Standard design width in meters of a single lane in UK
    shoulder_width = 3.3 # Apply to Motorways
    save_file = False
    step = False
    if step is True:
        # This section was written to get the original data and extract major roads 
        osm_roads = gpd.read_parquet(
                                os.path.join(
                                    incoming_data_path,
                                    "roads",
                                    "edges.gpq"
                                    )
                                )
        
        osm_roads = osm_roads[
                                (
                                    osm_roads["from_iso_a3"] == "GBR"
                                ) & (
                                    osm_roads["to_iso_a3"] == "GBR"
                                    )
                                ]

        osm_roads.to_parquet(
                                os.path.join(
                                        incoming_data_path,
                                        "roads",
                                        "GB_osm_roads.gpq"
                                    )
                                )
        osm_roads.to_file(
                            os.path.join(
                                        incoming_data_path,
                                        "roads",
                                        "GB_osm_roads.gpkg"
                                    ),driver="GPKG"
                            )

    gb_roadlinks = gpd.read_parquet(
                            os.path.join(
                                data_path,
                                "networks",
                                "GB_RoadLink.geoparquet"
                                )
                            )
    gb_roadlinks = gb_roadlinks.to_crs(epsg=epsg_meters)
    gb_roadlinks = gb_roadlinks[
                                [
                                    "gml_id","roadClassification",
                                    "formOfWay","averageWidth",
                                    "minimumWidth","geometry"
                                ]
                                ]
    osm_bridges = gpd.read_parquet(
                                os.path.join(
                                        incoming_data_path,
                                        "roads",
                                        "GB_osm_roads.gpq"
                                    )
                                )
    osm_bridges = osm_bridges[
                            (
                                osm_bridges["asset_type"] == "road_bridge"
                            ) & (
                                osm_bridges["tag_highway"] != "tertiary"
                                )
                            ]
    osm_bridges = osm_bridges.to_crs(epsg=epsg_meters)
    osm_bridges["length"] = osm_bridges.geometry.length
    osm_bridges["length_uom"] = "m"
    road_classes = ["Motorway","Others"]
    reduced_bridge_matches = []
    for r_c in road_classes:
        if r_c == "Motorway":
            bridges = osm_bridges[osm_bridges["tag_highway"] == "motorway"]
            gb_links = gb_roadlinks[gb_roadlinks["roadClassification"] == "Motorway"]
        else:
            bridges = osm_bridges[osm_bridges["tag_highway"] != "motorway"]
            gb_links = gb_roadlinks[gb_roadlinks["roadClassification"] != "Motorway"]
        
        bridges["geometry"] = bridges.geometry.buffer(buffer_distance)
        bridge_matches = gpd.sjoin(
                                bridges,gb_links,
                                how="inner",
                                predicate="intersects",
                                ).reset_index()
        
        bridge_matches.rename(columns={"geometry":"bridge_geometry"},inplace=True)
        bridge_matches = pd.merge(bridge_matches,gb_links[["gml_id","geometry"]],how="left",on=["gml_id"])
        bridge_matches["os_length_m"] = bridge_matches.geometry.length
        bridge_matches["intersect_length_m"
            ] = bridge_matches.progress_apply(
                            lambda x:x["bridge_geometry"].intersection(x["geometry"]).length,
                            axis=1)

        bridge_matches["fraction_length"
            ] = 1.0*bridge_matches["intersect_length_m"]/bridge_matches["os_length_m"]

        bridge_matches = bridge_matches.sort_values(by=["fraction_length"],ascending=False)
        # bridge_matches = bridge_matches.drop_duplicates(subset=["gml_id"],keep="first")
        bridge_matches.drop("bridge_geometry",axis=1,inplace=True)
        bridge_ids = list(set(bridge_matches["id"].values.tolist()))

        reduced_bridge_matches.append(bridge_matches[
                            (
                                bridge_matches["fraction_length"] > 0.60
                            ) | (
                                bridge_matches["intersect_length_m"] > bridge_matches["length"]
                            ) | (
                                bridge_matches["intersect_length_m"] > 2*buffer_distance  
                            )
                            ])

    reduced_bridge_matches = pd.concat(reduced_bridge_matches,axis=0,ignore_index=True)
    bridge_attributes = reduced_bridge_matches.groupby(
                            "id"
                            ).agg(
                                    {
                                        "averageWidth":"mean",
                                        "minimumWidth":"mean",
                                        "roadClassification":list,
                                        "formOfWay":list
                                    }
                            ).reset_index()
    bridge_attributes["roadClassification"] = bridge_attributes["roadClassification"].progress_apply(lambda x:','.join(list(set(x))))
    bridge_attributes["formOfWay"] = bridge_attributes["formOfWay"].progress_apply(lambda x:','.join(list(set(x))))
    reduced_bridge_ids = list(set(reduced_bridge_matches["id"].values.tolist()))
    bridges = osm_bridges[osm_bridges["id"].isin(reduced_bridge_ids)]
    bridges = pd.merge(bridges,bridge_attributes,how="left",on=["id"])
    bridges["lanes"] = bridges.progress_apply(
                                    lambda x: assign_lanes_to_roads(
                                        x,"roadClassification","formOfWay",
                                        lane_width,shoulder_width
                                        ),axis=1
                                    )


    bridges.to_parquet(
                     os.path.join(
                                data_path,
                                "networks",
                                "GB_osm_bridges.geoparquet"
                            )
                    )

    bridges.to_file(
                     os.path.join(
                                data_path,
                                "networks",
                                "GB_osm_bridges.gpkg"
                            ),driver="GPKG"
                    )

    reduced_bridge_matches = reduced_bridge_matches[["id","gml_id"]]
    gb_roads_mapping = pd.read_csv(
                                    os.path.join(
                                            data_path,
                                            "networks",
                                            "OSOpenRoadLookUpTable_major_roads.csv"
                                        )
                                    )
    reduced_bridge_matches = pd.merge(
                            reduced_bridge_matches,
                            gb_roads_mapping,
                            how="left",right_on=["Roadlink_Id"],left_on=["gml_id"])
    reduced_bridge_matches.drop("gml_id",axis=1,inplace=True)
    reduced_bridge_matches.rename(columns={"id":"osm_bridge_id"},inplace=True)
    reduced_bridge_matches.to_csv(
                                    os.path.join(
                                            data_path,
                                            "networks",
                                            "OSM_OSMasterMap_OSOpenRoadLookUpTable_bridges.csv"
                                        ),index=False
                                    )

    
if __name__ == '__main__':
    CONFIG = load_config() 
    main(CONFIG)