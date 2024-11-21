#!/usr/bin/env python
# coding: utf-8
import sys
import os
import pandas as pd
import geopandas as gpd
from utils import *


def main(config):
    data_path = config['paths']['data']
    epsg_meters = 27700
    gb_osm_edges = gpd.read_parquet(
                            os.path.join(
                                data_path,
                                "miraca",
                                "inputs",
                                "incoming_data",
                                "europe-latest_filter-road-tertiary",
                                "edges.gpq"
                                )
                            )
    gb_osm_nodes = gpd.read_parquet(
                            os.path.join(
                                data_path,
                                "miraca",
                                "inputs",
                                "incoming_data",
                                "europe-latest_filter-road-tertiary",
                                "nodes.gpq"
                                )
                            )
    gb_osm_nodes = gb_osm_nodes[gb_osm_nodes["iso_a3"] == "GBR"]
    gb_osm_nodes = gb_osm_nodes.to_crs(epsg=epsg_meters)

    gb_osm_edges = gb_osm_edges[
                            (
                                gb_osm_edges["from_iso_a3"] == "GBR"
                            ) & (
                                gb_osm_edges["to_iso_a3"] == "GBR"
                            )]
    gb_osm_edges.drop(["id","from_id","to_id"],axis=1,inplace=True)
    gb_osm_edges = gb_osm_edges.to_crs(epsg=epsg_meters)
    gb_osm_edges["length_m"] = gb_osm_edges.geometry.length
    
    os_nodes = gpd.read_parquet(
                            os.path.join(
                                data_path,
                                "DAFNI_NIRD",
                                "processed_data",
                                "networks",
                                "road_node_file.geoparquet"
                                )
                            )                                   
    os_nodes.rename(columns={"id":"os_id"},inplace=True)
    os_nodes = os_nodes.to_crs(epsg=epsg_meters)

    os_nodes = ckdnearest(os_nodes,gb_osm_nodes)
    os_nodes.to_parquet(
                            os.path.join(
                                data_path,
                                "DAFNI_NIRD",
                                "processed_data",
                                "networks",
                                "os_osm_proximity.geoparquet"
                                )
                            )

    # network = create_network(gb_osm_edges,os_nodes,"road")
    # edges = network.edges
    # nodes = network.nodes
    # edges, nodes = components(edges,nodes,"id")
    # del os_nodes,gb_osm_edges
    # gpd.GeoDataFrame(
    #                 nodes,geometry="geometry",crs=f"EPSG:{epsg_meters}"
    #                 ).to_parquet(
    #                         os.path.join(
    #                             data_path,
    #                             "DAFNI_NIRD",
    #                             "processed_data",
    #                             "networks",
    #                             "gb_osm_os_nodes.gpq"
    #                             )
    #                     )
    # nodes.rename(columns={"id":"nid"},inplace=True)
    # os_roads = gpd.read_parquet(
    #                         os.path.join(
    #                             data_path,
    #                             "DAFNI_NIRD",
    #                             "processed_data",
    #                             "networks",
    #                             "road_link_file.geoparquet"
    #                             )
    #                         )
    # os_roads = os_roads[
    #                     [
    #                         "start_node","end_node",
    #                         "id","fictitious",
    #                         "road_classification","road_function",
    #                         "form_of_way","road_classification_number",
    #                         "road_structure",
    #                         "average_toll_cost",
    #                         "urban"
    #                     ]
    #                     ]
    # os_roads = pd.merge(
    #                     os_roads,nodes[["os_id","nid"]],
    #                     how="left",
    #                     left_on=["start_node"],
    #                     right_on=["os_id"])
    # os_roads.rename(columns={"nid":"origin_id"},inplace=True)
    # os_roads.drop(["start_node","os_id"],axis=1,inplace=True)

    # os_roads = pd.merge(
    #                     os_roads,nodes[["os_id","nid"]],
    #                     how="left",
    #                     left_on=["end_node"],
    #                     right_on=["os_id"])
    # os_roads.rename(columns={"nid":"destination_id"},inplace=True)
    # os_roads.drop(["end_node","os_id"],axis=1,inplace=True)

    # osm_graph = create_igraph_from_dataframe(
    #                     edges[
    #                         ["from_id","to_id","id","length_m"]
    #                         ]
    #                     )
    # os_roads = network_od_paths_assembly(os_roads, osm_graph,
    #                             "length_m","id")
    # os_roads.drop(["origin_id","destination_id","length_m"],axis=1,inplace=True)
    # os_roads.rename(columns={"id","os_id"},inplace=True)

    # osm_matches = [] 
    # cols = [c for c in os_roads.columns.values.tolist() if c != "edge_path"]
    # for r in os_roads.itertuples():
    #     rl = []
    #     ep = r.edge_path
    #     rl.append(ep)
    #     for c in cols:
    #         rl.append([getattr(r,c)]*len(ep))

    #     osm_matches += list(zip(*rl))

    # osm_matches = pd.DataFrame(osm_matches,columns=["id"] + cols)
    # osm_matches.to_csv("osm_matches.csv",index=False)
    # edges = pd.merge(edges,osm_matches,how="left",on=["id"])

    # gpd.GeoDataFrame(
    #                 edges,geometry="geometry",crs=f"EPSG:{epsg_meters}"
    #                 ).to_parquet(
    #                         os.path.join(
    #                             data_path,
    #                             "DAFNI_NIRD",
    #                             "processed_data",
    #                             "networks",
    #                             "gb_osm_os_edges.gpq"
    #                             )
    #                     )
    
if __name__ == '__main__':
    CONFIG = load_config() 
    main(CONFIG)