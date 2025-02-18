#!/usr/bin/env python
# coding: utf-8
# Merge the UNRA national roads dataset with an OSM Roads dataset 
import sys
import os
import pandas as pd
import igraph as ig
import snkit
import geopandas as gpd
from uganda_utils import *
from tqdm import tqdm
tqdm.pandas()

def main():
    incoming_data_path = "folder/path/to/UNRA roads/and/OSM roads files" # Probably best to have both in same folder
    output_data_path = "folder/path/to/final/results/file"
    epsg_meters = 32735 # Local projection system for Uganda

    # Fill in the names of the columns in the UNRA dataset, whose properties you want to migrate to OSM 
    unra_attributes = ["col1","col2"]
    # Read the UNRA roads file. Assumption below is that it is GPKG
    unra_roads = gpd.read_file(
                     os.path.join(
                                incoming_data_path,
                                "UNRA_data_file.gpkg"
                            )
                    )
    # Convert to the local projection system for accurate distance calculations
    unra_roads = unra_roads.to_crs(epsg=epsg_meters)
    # First convert the UNRA roads file to a network with edges and nodes
    # This is done in order to identify the end points of roads, 
    # that will be matched to OSM nodes  
    network = create_network(unra_roads,id_prefix="unra_id")
    # Get the edges and nodes file as GeoDataFrames
    unra_edges = gpd.GeoDataFrame(network.edges,geometry="geometry",crs=f"EPSG:{epsg_meters}")
    unra_nodes = gpd.GeoDataFrame(network.nodes,geometry="geometry",crs=f"EPSG:{epsg_meters}")
    # Will have to rename ID columns because they most probably have the same name as OSM ID columns
    unra_edges.rename(columns={"id":"unra_id","from_id":"from_unra_id","to_id":"to_unra_id"},inplace=True)
    unra_nodes.rename(columns={"id":"nid"},inplace=True)
    # If you want, then save the intermediate result to see how it looks
    save_file = False
    if save_file is True:
        unra_edges.to_file(os.path.join(output_data_path,"unra_roads.gpkg"),layer="edges")
        unra_nodes.to_file(os.path.join(output_data_path,"unra_roads.gpkg"),layer="nodes")

    # Read the OSM edges and nodes files created with Open-GIRA
    osm_id_column = "id" # Name of the ID column in OSM. It is generally same in edges and nodes files
    from_id_column = "from_id"
    to_id_column = "to_id"
    # Select the OSM road classes that best match UNRA roads. 
    # This will make sure the matches are done with the right kind of road
    osm_road_classes = ["motorway","primary","secondary"]
    osm_edges = gpd.read_parquet(
                                    os.path.join(
                                            incoming_data_path,
                                            "uganda_osm_edges.geoparquet")
                                    )
    osm_nodes = gpd.read_parquet(
                                    os.path.join(
                                            incoming_data_path,
                                            "uganda_osm_nodes.geoparquet")
                                    )
    edges = osm_edges[osm_edges["tag_highway"].isin(osm_road_classes)]
    selected_node_ids = list(
                            set(
                                edges[from_id_column].values.tolist() + edges[to_id_column].values.tolist()
                                )
                            )
    del edges
    nodes = osm_nodes[osm_nodes[osm_id_column].isin(selected_node_ids)]
    # Convert to the local projection system for accurate distance calculations
    nodes = nodes.to_crs(epsg=epsg_meters)

    # Find the nearest OSM nodes to the UNRA nodes
    gdf_p = gpd.sjoin_nearest(
                            unra_nodes,nodes,how="left",distance_col="distance_m"
                            )
    # Assume some distance threshold, within which we assume the UNRA and OSM points are the same
    # This is done to avoid chopping the OSM roads into small line segments
    distance_threshold = 5.0 # Value in meters
    gdf_p["new_id"] = np.where(
                            gdf_p["distance_m"] < distance_threshold,
                            gdf_p[osm_id_column],
                            gdf_p["nid"]
                            )
    # Will have to modify the OSM network by adding in the new nodes
    network = create_network(osm_edges,nodes=gdf_p[["unra_id","new_id","geometry"]]) # This might take time
    edges = gpd.GeoDataFrame(network.edges,geometry="geometry",crs=f"EPSG:{epsg_meters}")
    edges["length_m"] = edges.geometry.length

    # Add the new matched OSM id to the UNRA edges 
    unra_edges = pd.merge(unra_edges,gdf_p[["nid","new_id"]],how="left",left_on=["from_id"],right_on=["nid"])
    unra_edges.rename(columns={"new_id":"from_new_id"},inplace=True)
    unra_edges.drop(["nid"],axis=1,inplace=True)
    unra_edges = pd.merge(unra_edges,gdf_p[["nid","new_id"]],how="left",left_on=["to_id"],right_on=["nid"])
    unra_edges.rename(columns={"new_id":"to_new_id"},inplace=True)
    unra_edges.drop(["nid"],axis=1,inplace=True)

    # Find the shortest path between the UNRA points of the OSM network and assign the new properties 
    graph = create_igraph_from_dataframe(edges[["from_id","to_id","id","length_m"]])
    roads_with_attributes = []
    for row in unra_edges.itertuples():
        path = graph.get_shortest_paths(
                                        row.from_new_id, 
                                        row.to_new_id, 
                                        weights="length_m", output="epath")[0]

        connected_roads = []
        connected_attributes = []
        if path:
            for n in path:
                connected_roads.append(graph.es[n]["id"])
            connected_attributes.append(connected_roads)
            for unra_at in unra_attributes:
                attr_names = [getattr(row,unra_at)]*len(connected_roads)
                    connected_attributes.append(attr_names)

            roads_with_attributes += list(zip(*connected_attributes))

    roads_with_attributes = pd.DataFrame(roads_with_attributes,columns=["id"] + unra_attributes)
                                 
    # Could be possible that a road might have multiple attributes
    roads_with_attributes = roads_with_attributes.groupby("id").agg(
                                                            dict([(u,list) for u in unra_attributes])
                                                            ).reset_index()
    
    

    # This will covert values from a list to a string seperated by / 
    # Example: ['a','b','c'] becomes 'a/b/c' 
    for u in unra_attributes:
        roads_with_attributes[u] = roads_with_attributes[u].progress_apply(lambda a: "/".join(list(set(a))))
                         
    
    edges = pd.merge(edges,roads_with_attributes,how="left",on=["id"])
    edges = gpd.GeoDataFrame(edges,geometry="geometry",crs=f"EPSG:{epsg_meters}")
    nodes = gpd.GeoDataFrame(network.nodes,geometry="geometry",crs=f"EPSG:{epsg_meters}")

    edges,nodes = components(edges,nodes) # This tells us about the connectivity of our network
    nodes = add_node_degree(edges,nodes) # This tells us about degree of nodes

    edges.to_file(os.path.join(output_data_path,"uganda_unra_osm_edges.geoparquet"))
    nodes.to_file(os.path.join(output_data_path,"uganda_unra_osm_nodes.geoparquet"))
    
if __name__ == '__main__':
    main()