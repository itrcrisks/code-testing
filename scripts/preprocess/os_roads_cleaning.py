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

        if x[fw_column] == "Motorway":
            lanes = max(lanes,4)
    else:
        if (x[fw_column] == "Single Carriageway") and (lanes % 2 != 0):
            lanes = lanes - 1

    
    if lanes == 0:
        lanes += 1
    
    return lanes

def main(config):
    data_path = config['paths']['data']
    epsg_meters = 27700
    lane_width = 3.65 # Standard design width in meters of a single lane in UK
    shoulder_width = 3.3 # Apply to Motorways
    gb_roadlinks_columns = [
                                'gml_id','roadClassification',
                                'formOfWay','averageWidth',
                                'minimumWidth','startGradeSeparation',
                                'endGradeSeparation',
                                'max_z','min_z','mean_z'
                            ]

    save_file = False
    gb_os_roads = gpd.read_parquet(
                            os.path.join(
                                data_path,
                                "networks",
                                "road_link_file.geoparquet"
                                )
                            )
    step = False
    if step is True:
        # This section was written to get the original data and extract major roads 
        gb_roads_mapping = pd.read_csv(
                                os.path.join(
                                    data_path,
                                    "networks",
                                    "OSOpenRoadLookUpTable_2023_09.csv"
                                )
                                )
        gb_roads_mapping = gb_roads_mapping[
                            gb_roads_mapping["OSOpenRoads_RoadLinkIdentifier"
                            ].isin(gb_os_roads["id"].values.tolist())
                            ]
        if save_file is True:
            gb_roads_mapping.to_csv(
                                    os.path.join(
                                            data_path,
                                            "networks",
                                            "OSOpenRoadLookUpTable_major_roads.csv"
                                        ),
                                    index=False)
        gb_roadlinks = gpd.read_parquet(
                                os.path.join(
                                    data_path,
                                    "networks",
                                    "GB_Highways.pq"
                                    )
                                )
        
        gb_roadlinks = gb_roadlinks[
                            (
                                gb_roadlinks["gml_id"
                                    ].isin(gb_roads_mapping["Roadlink_Id"].values.tolist())
                            ) | (
                                gb_roadlinks["roadClassification"].isin(["Motorway","A Road","B Road"])
                            )
                            ]

        gb_roadlinks.to_parquet(
                                os.path.join(
                                        data_path,
                                        "networks",
                                        "GB_RoadLink_OSOpenroads.gpq"
                                    )
                                )

    save_file = False
    if save_file is True:
        gb_roadlinks = gpd.read_parquet(
                                os.path.join(
                                        data_path,
                                        "networks",
                                        "GB_RoadLink_OSOpenroads.gpq"
                                    )
                                )
        gb_roadlinks["averageWidth"] = gb_roadlinks["averageWidth"].fillna(2*lane_width)
        gb_roadlinks["minimumWidth"] = np.where(
                                                gb_roadlinks["minimumWidth"].isna(),
                                                gb_roadlinks["averageWidth"],
                                                gb_roadlinks["minimumWidth"]
                                                )
        gb_roadlinks["max_z"] = gb_roadlinks['geometry'].progress_apply(lambda geom: np.max([coord[2] for coord in geom.coords]))
        gb_roadlinks["min_z"] = gb_roadlinks['geometry'].progress_apply(lambda geom: np.min([coord[2] for coord in geom.coords]))
        gb_roadlinks["mean_z"] = gb_roadlinks['geometry'].progress_apply(lambda geom: np.mean([coord[2] for coord in geom.coords]))
        gb_roadlinks["geometry"] = gb_roadlinks["geometry"].progress_apply(lambda geom: shapely.force_2d(geom))
        gb_roadlinks["lanes"] = gb_roadlinks.progress_apply(
                                    lambda x: assign_lanes_to_roads(
                                        x,"roadClassification","formOfWay",
                                        lane_width,shoulder_width
                                        ),axis=1
                                    )
        gb_roadlinks.to_parquet(
                                os.path.join(
                                        data_path,
                                        "networks",
                                        "GB_RoadLink_OSOpenroads.gpq"
                                    )
                                )
        gb_roadlinks = gb_roadlinks[gb_roadlinks["roadClassification"].isin(["Motorway","A Road","B Road"])]
        # This was just to save the file to visualise. Can ignore during re-run

        gb_roadlinks[gb_roadlinks_columns + ["lanes","geometry"]
            ].to_file(
                        os.path.join(
                        data_path,
                        "networks",
                        "GB_RoadLink_major.gpkg"
                    ),driver="GPKG")

    os_roads_mapping = False
    if os_roads_mapping is True:
        gb_roadlinks = gpd.read_parquet(
                                os.path.join(
                                        data_path,
                                        "networks",
                                        "GB_RoadLink_OSOpenroads.gpq"
                                    )
                                )
        gb_roads_mapping = pd.read_csv(
                                    os.path.join(
                                            data_path,
                                            "networks",
                                            "OSOpenRoadLookUpTable_major_roads.csv"
                                        )
                                    )
        gb_roads_mapping = pd.merge(
                                gb_roads_mapping,
                                gb_roadlinks[gb_roadlinks_columns],
                                how="left",left_on=["Roadlink_Id"],right_on=["gml_id"])
        gb_roads_mapping = gb_roads_mapping.groupby(
                            "OSOpenRoads_RoadLinkIdentifier"
                            ).agg(
                                    {
                                        "averageWidth":"mean",
                                        "minimumWidth":"mean",
                                        "max_z":"max",
                                        "min_z":"min",
                                        "mean_z":"mean"
                                    }
                            ).reset_index()
        if save_file is True:
            gb_roads_mapping.to_csv(
                                        os.path.join(
                                                data_path,
                                                "networks",
                                                "OSOpenRoadLookUpTable_widths.csv"
                                            ),index=False
                                        )
        gb_os_roads = pd.merge(
                            gb_os_roads,
                            gb_roads_mapping,
                            how="left",left_on=["id"],
                            right_on=["OSOpenRoads_RoadLinkIdentifier"]
                            )
        # Assign widths to the non-matched roads
        # Some roads are unmatched between the OS roads and the OS Mastermap highways
        # Widths for these unmatched are estimated by simply assigning the widths of next roads
        nan_roads = gb_os_roads[gb_os_roads["averageWidth"].isna()]
        non_nan_roads = gb_os_roads[~gb_os_roads["averageWidth"].isna()]
        df = pd.DataFrame()
        while len(nan_roads.index) > 0:
            if len(df.index) == len(nan_roads.index):
                nan_roads["averageWidth"] = nan_roads["minimumWidth"] = 2*lane_width
            else:
                nan_roads.drop(["averageWidth","minimumWidth"],axis=1,inplace=True)
                from_nodes = non_nan_roads[["from_id","averageWidth","minimumWidth"]]
                from_nodes.rename(columns={"from_id":"nid"},inplace=True)
                to_nodes = non_nan_roads[["to_id","averageWidth","minimumWidth"]]
                to_nodes.rename(columns={"to_id":"nid"},inplace=True)

                nodes = pd.concat(
                                    [
                                        from_nodes,
                                        to_nodes
                                    ],
                                axis=0,ignore_index=True)
                nodes = nodes.groupby(
                                    "nid"
                                    ).agg(
                                            {
                                                "averageWidth":"mean",
                                                "minimumWidth":"min",
                                            }
                                    ).reset_index()

                nan_roads = pd.merge(nan_roads,nodes,how="left",left_on=["from_id"],right_on=["nid"])
                nan_roads.rename(columns={"averageWidth":"from_width","minimumWidth":"from_width"},inplace=True)
                nan_roads.drop(["nid"],axis=1,inplace=True)

                nan_roads = pd.merge(nan_roads,nodes,how="left",left_on=["to_id"],right_on=["nid"])
                nan_roads.rename(columns={"averageWidth":"to_width","minimumWidth":"to_width"},inplace=True)
                nan_roads.drop(["nid"],axis=1,inplace=True)

                nan_roads["averageWidth"] = nan_roads[["from_width","to_width"]].mean(axis=1,skipna=True)
                nan_roads["minimumWidth"] = nan_roads[["from_width","to_width"]].min(axis=1,skipna=True)

                nan_roads.drop(["from_width","to_width"],axis=1,inplace=True) 
            non_nan_roads = pd.concat(
                                [
                                    non_nan_roads,
                                    nan_roads[~nan_roads["averageWidth"].isna()]
                                ],axis=0,ignore_index=True)
            df = nan_roads.copy()
            nan_roads = nan_roads[nan_roads["averageWidth"].isna()]

        gb_os_roads = non_nan_roads.copy()
        for c in ["averageWidth","minimumWidth"]:
            gb_os_roads[c] = np.where(
                                        gb_os_roads["form_of_way"] == "Collapsed Dual Carriageway",
                                        2*gb_os_roads[c],
                                        gb_os_roads[c]
                                    )
        gb_os_roads["lanes"] = gb_os_roads.progress_apply(
                                    lambda x: assign_lanes_to_roads(
                                        x,"road_classification","form_of_way",
                                        lane_width,shoulder_width,collapsed_roads=True
                                        ),axis=1
                                    )

        gb_os_roads.to_parquet(
                                os.path.join(
                                    data_path,
                                    "networks",
                                    "GB_road_link_file.geoparquet"
                                    )
                                )
        gb_os_roads.to_file(
                                os.path.join(
                                    data_path,
                                    "networks",
                                    "GB_road_link_file.gpkg"
                                    ),
                                driver="GPKG"
                                )
    
    create_topology = False
    if create_topology is True:
        gb_roadlinks = gpd.read_parquet(
                                os.path.join(
                                        data_path,
                                        "networks",
                                        "GB_RoadLink_OSOpenroads.gpq"
                                    )
                                )
        gb_roadlinks = gb_roadlinks[gb_roadlinks["roadClassification"].isin(["Motorway","A Road","B Road"])]
        gb_roadlinks = gb_roadlinks.to_crs(epsg=epsg_meters)
        network = create_network(gb_roadlinks,id_prefix="road")
        edges = network.edges
        nodes = network.nodes
        edges, nodes = components(edges,nodes,"id")

        nodes = gpd.GeoDataFrame(
                        nodes,geometry="geometry",crs=f"EPSG:{epsg_meters}"
                        )
        nodes = add_node_degree(edges,nodes)
        nodes.to_parquet(
                            os.path.join(
                                data_path,
                                "networks",
                                "GB_RoadNode.geoparquet"
                                )
                        )
        edges = gpd.GeoDataFrame(
                        edges,geometry="geometry",crs=f"EPSG:{epsg_meters}"
                        )
        edges.to_parquet(
                            os.path.join(
                                data_path,
                                "networks",
                                "GB_RoadLink.geoparquet"
                                )
                        )

        nodes.to_file(
                            os.path.join(
                                data_path,
                                "networks",
                                "GB_Roads.gpkg"
                                ),layer="nodes",driver="GPKG"
                        )
        edges[
                gb_roadlinks_columns + ["id","from_id","to_id",
                "lanes","component","geometry"]
            ].to_file(
                        os.path.join(
                            data_path,
                            "networks",
                            "GB_Roads.gpkg"
                            ),layer="edges",driver="GPKG"
                    )

    node_seperation = True
    if node_seperation is True:
        edges = gpd.read_parquet(
                            os.path.join(
                                data_path,
                                "networks",
                                "GB_RoadLink.geoparquet"
                                )
                        )
        print (edges)
        nodes = gpd.read_parquet(
                            os.path.join(
                                data_path,
                                "networks",
                                "GB_RoadNode.geoparquet"
                                )
                        )
        print (nodes)
        
        from_nodes = edges[["from_id","startGradeSeparation"]]
        from_nodes.rename(columns={"from_id":"id","startGradeSeparation":"grade_separation"},inplace=True)

        to_nodes = edges[["to_id","endGradeSeparation"]]
        to_nodes.rename(columns={"to_id":"id","endGradeSeparation":"grade_separation"},inplace=True)

        from_to_nodes = pd.concat([from_nodes,to_nodes],axis=0,ignore_index=True)
        # from_to_nodes = from_to_nodes.groupby(["id"])["grade_separation"].apply(list).reset_index()
        # from_to_nodes["grade_separation"] = from_to_nodes["grade_separation"].progress_apply(lambda x: list(set(x)))
        # from_to_nodes["grade_difference"] = from_to_nodes["grade_separation"].progress_apply(lambda x: 1 if len(x) > 1 else 0)

        save_file = True
        if save_file is True:
            # nodes = pd.merge(nodes,from_to_nodes,how="left",on=["id"])
            # nodes["grade_separation"] = nodes["grade_separation"].astype(str)
            from_to_nodes["nid"] = from_to_nodes.progress_apply(lambda x: f"{x['id']}_{x['grade_separation']}",axis=1)
            from_to_nodes = from_to_nodes.drop_duplicates(subset=["nid"],keep="first")
            from_to_nodes = pd.merge(from_to_nodes[["id","nid"]],nodes,how="left",on=["id"])
            nodes = from_to_nodes.copy()
            nodes.drop(["id","degree","component"],axis=1,inplace=True)
            nodes.rename(columns={"nid":"id"},inplace=True)
            
            edges["from_id"] = edges.progress_apply(lambda x: f"{x['from_id']}_{x['startGradeSeparation']}",axis=1)
            edges["to_id"] = edges.progress_apply(lambda x: f"{x['to_id']}_{x['endGradeSeparation']}",axis=1)

            edges, nodes = components(edges,nodes,"id")

            nodes = gpd.GeoDataFrame(
                            nodes,geometry="geometry",crs=f"EPSG:{epsg_meters}"
                            )

            nodes = add_node_degree(edges,nodes)

            edges.to_parquet(
                        os.path.join(
                            data_path,
                            "networks",
                            "GB_RoadLink.geoparquet"
                            )
                    )
            print (edges)
            nodes.to_parquet(
                            os.path.join(
                                data_path,
                                "networks",
                                "GB_RoadNode.geoparquet"
                                )
                        )
            print (nodes)
            save_gpkg = False
            if save_gpkg is True:
                edges[
                        gb_roadlinks_columns + [
                            "id","from_id","to_id",
                            "lanes","component","geometry"
                            ]
                    ].to_file(
                                os.path.join(
                                    data_path,
                                    "networks",
                                    "GB_Roads.gpkg"
                                    ),layer="edges",driver="GPKG"
                            )

                nodes.to_file(
                                    os.path.join(
                                        data_path,
                                        "networks",
                                        "GB_Roads.gpkg"
                                        ),layer="nodes",driver="GPKG"
                                )

if __name__ == '__main__':
    CONFIG = load_config() 
    main(CONFIG)