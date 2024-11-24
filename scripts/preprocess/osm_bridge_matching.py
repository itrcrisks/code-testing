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

def main(config):
    incoming_data_path = config['paths']['incoming_data']
    data_path = config['paths']['data']
    epsg_meters = 27700
    buffer_distance = 4.0
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
    gb_roadlinks = gb_roadlinks[["gml_id","geometry"]]
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
    bridges = osm_bridges.copy()
    bridges["geometry"] = bridges.geometry.buffer(buffer_distance)
    bridge_matches = gpd.sjoin(
                            bridges,gb_roadlinks,
                            how="inner",
                            predicate="intersects",
                            ).reset_index()
    
    bridge_matches.rename(columns={"geometry":"bridge_geometry"},inplace=True)
    bridge_matches = pd.merge(bridge_matches,gb_roadlinks,how="left",on=["gml_id"])
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

    reduced_bridge_matches = bridge_matches[
                        (
                            bridge_matches["fraction_length"] > 0.60
                        ) | (
                            bridge_matches["intersect_length_m"] > bridge_matches["length"]
                        ) | (
                            bridge_matches["intersect_length_m"] > 2*buffer_distance  
                        )
                        ]
    reduced_bridge_ids = list(set(reduced_bridge_matches["id"].values.tolist()))
    bridges = osm_bridges[osm_bridges["id"].isin(reduced_bridge_ids)]
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