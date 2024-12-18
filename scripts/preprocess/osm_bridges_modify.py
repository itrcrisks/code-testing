#!/usr/bin/env python
# coding: utf-8
import sys
import os
import pandas as pd
import igraph as ig
import snkit
import shapely
from collections import defaultdict
import ast
import geopandas as gpd
import shapely
from utils import *
from tqdm import tqdm
tqdm.pandas()

def snap_bridges_to_roads(bridges,bridge_nodes,edges,distance_threshold=20):
    road_classes = ["motorway","secondary","primary","trunk"]
    br_nodes = []
    for r_c in road_classes:
        osm_bridges = bridges[bridges["tag_highway"] == r_c]
        if r_c == "motorway":
            gb_links = edges[edges["road_classification"] == "Motorway"] 
        elif r_c == "secondary":
            gb_links = edges[edges["road_classification"] == "B Road"]
        else:
            gb_links = edges[edges["road_classification"] == "A Road"]

        osm_nodes_ids = list(
                            set(
                                osm_bridges.from_id.values.tolist() + osm_bridges.to_id.values.tolist()
                                )
                            )
        osm_nodes = bridge_nodes[bridge_nodes["id"].isin(osm_nodes_ids)] 
        gb_links = gb_links.reset_index()
        network = snkit.Network(osm_nodes,gb_links)
        print("* Done with network creation")
        network = snkit.network.snap_nodes(network)
        print ('* Done with snapping nodes to edges')

        br_nodes.append(network.nodes)
    
    br_nodes = gpd.GeoDataFrame(
                    pd.concat(br_nodes,axis=0,ignore_index=True),
                    geometry="geometry",crs=bridge_nodes.crs)
    br_nodes = br_nodes.drop_duplicates(subset=["id"],keep="first")
    return br_nodes

def main(config):
    incoming_data_path = config['paths']['incoming_data']
    data_path = config['paths']['data']
    epsg_meters = 27700
    buffer_distance = 6.0
    lane_width = 3.65 # Standard design width in meters of a single lane in UK
    shoulder_width = 3.3 # Apply to Motorways
    save_file = False

    bridges = gpd.read_parquet(
                     os.path.join(
                                data_path,
                                "networks",
                                "GB_osm_bridges.geoparquet"
                            )
                    )
    bridge_nodes = gpd.read_parquet(
                                    os.path.join(
                                            data_path,
                                            "networks",
                                            "GB_osm_bridges_nodes.geoparquet")
                                    )
    bridge_nodes["form_of_road_node"] = "bridge_end_point"
    bridge_nodes = bridge_nodes[["id","form_of_road_node","geometry"]]
    bridge_nodes = bridge_nodes.to_crs(epsg=epsg_meters)
    edges = gpd.read_parquet(
                                os.path.join(
                                    data_path,
                                    "networks",
                                    "GB_road_link_file.geoparquet")
                                )
    edges = edges.to_crs(epsg=epsg_meters)
    nodes = gpd.read_parquet(
                            os.path.join(
                                data_path,
                                "networks",
                                "road_node_file.geoparquet")
                            )
    nodes = nodes[["id","form_of_road_node","geometry"]]

    step = False
    if step is True:
        add_bridges = []
        snapped_nodes = []
        unique_ids = []
        gdf_p = bridges[["from_id","to_id"]]
        for nt in ["from_id","to_id"]:
            gdf_p = gpd.GeoDataFrame(
                            pd.merge(
                                gdf_p,
                                bridge_nodes[["id","geometry"]],
                                how="left",left_on=[nt],
                                right_on=["id"]),
                            geometry="geometry",crs = bridge_nodes.crs
                                )
            gdf_p.drop("id",axis=1,inplace=True)
            gdf_p = gpd.sjoin_nearest(
                            gdf_p, edges[["id","geometry"]]
                            ).merge(edges[["geometry"]], left_on="index_right", right_index=True)
            gdf_p[f"{nt}_distance"] = gdf_p.progress_apply(
                                                lambda r: r["geometry_x"].distance(r["geometry_y"]), axis=1)
            gdf_p.rename(columns={"id":f"{nt}_edge"},inplace=True)
            gdf_p.drop(["index_right","geometry_x","geometry_y"],axis=1,inplace=True)

        gdf_p = pd.merge(gdf_p,bridges[["id","length","from_id","to_id"]],how="left",on=["from_id","to_id"])
        unique_snap = gdf_p[gdf_p["from_id_edge"] == gdf_p["to_id_edge"]]
        unique_snap = unique_snap.drop_duplicates(subset=["id"],keep="first")
        nonunique_snap = gdf_p[~gdf_p["id"].isin(unique_snap.id.values.tolist())]
        unique_snap = unique_snap[
                                    (
                                        unique_snap["from_id_distance"] < 30
                                    ) & (
                                        unique_snap["to_id_distance"] < 30
                                    )
                                ]
        nonunique_snap = nonunique_snap[
                                    (
                                        nonunique_snap["from_id_distance"] < 30
                                    ) & (
                                        nonunique_snap["to_id_distance"] < 30
                                    )
                                ]
        nonunique_snap = nonunique_snap.sort_values(
                                    by=["from_id_distance","to_id_distance"],
                                    ascending=[True,True]
                                    )
        nonunique_snap = nonunique_snap.drop_duplicates(subset=["id"],keep="first")
        add_bridges.append(unique_snap)

        unique_ids += list(set(unique_snap.from_id.values.tolist() + unique_snap.to_id.values.tolist()))
        gdf_n = bridge_nodes[bridge_nodes["id"].isin(unique_ids)]
        network = snkit.Network(gdf_n, edges)
        network = snkit.network.snap_nodes(network)
        snapped_nodes.append(network.nodes)

        iteration = 0
        snap_dict = defaultdict(list)
        remaining_edges = bridges["id"].values.tolist()
        while (len(nonunique_snap.index) > 0) and (len(set(remaining_edges) - set(nonunique_snap["id"].values.tolist())) > 0):
            remaining_edges = nonunique_snap["id"].values.tolist()
            unadded_bridges = []
            nonunique_snap["iteration"] = iteration
            nonunique_snap["from_id_remain_fixed"
                ] = np.where(
                        nonunique_snap["from_id"].isin(unique_ids),
                        1,
                        0
                    )
            nonunique_snap["to_id_remain_fixed"
                ] = np.where(
                        nonunique_snap["to_id"].isin(unique_ids),
                        1,
                        0
                    )
            print (iteration)
            print (nonunique_snap)

            nonunique_snap = nonunique_snap[
                                        ~(
                                            (
                                                nonunique_snap["from_id_remain_fixed"] == 1
                                            ) & (
                                                nonunique_snap["to_id_remain_fixed"] == 1
                                            )
                                        )]
            if len(nonunique_snap.index) > 0:
                # nonunique_snap = nonunique_snap.head(100)
                checked_ids = []
                for row in nonunique_snap.itertuples():
                    gdf_n = bridge_nodes[bridge_nodes["id"].isin([row.from_id,row.to_id])]
                    # print (gdf_n)
                    from_edge = new_from_edge = row.from_id_edge
                    to_edge = row.to_id_edge
                    pairs = [(from_edge,to_edge)]
                    if iteration > 0:
                        if row.from_id_remain_fixed == 0:
                            gdf_p = bridge_nodes[bridge_nodes["id"] == row.from_id]
                            if row.from_id not in checked_ids:
                                snap_dict[row.from_id].append(from_edge)
                                checked_ids.append(row.from_id)

                            gdf_p = gpd.sjoin_nearest(
                                gdf_p, edges[~edges["id"].isin(snap_dict[row.from_id])][["id","geometry"]]
                                ).reset_index()
                            new_from_edge = gdf_p["id_right"].values[0]
                            pairs.append((new_from_edge,to_edge))

                        if row.to_id_remain_fixed == 0:
                            gdf_p = bridge_nodes[bridge_nodes["id"] == row.to_id]
                            if row.to_id not in checked_ids:
                                snap_dict[row.to_id].append(to_edge)
                                checked_ids.append(row.from_id)
                            gdf_p = gpd.sjoin_nearest(
                                gdf_p, edges[~edges["id"].isin(snap_dict[row.to_id])][["id","geometry"]]
                                ).reset_index()
                            new_to_edge = gdf_p["id_right"].values[0]
                            pairs.append((from_edge,new_to_edge))
                            pairs.append((new_from_edge,new_to_edge))

                        pairs = list(set(pairs[1:]))


                        # mask = (nonunique_snap["id"] == row.id)
                        # nonunique_snap.loc[mask, "from_id_edge"] = from_edge
                        # nonunique_snap.loc[mask, "to_id_edge"] = to_edge

                    counter = 0
                    for idx,(from_edge,to_edge) in enumerate(pairs):
                        if counter == 1:
                            break
                        mask = (nonunique_snap["id"] == row.id)
                        nonunique_snap.loc[mask, "from_id_edge"] = from_edge
                        nonunique_snap.loc[mask, "to_id_edge"] = to_edge
                        gdf_e = edges[edges["id"].isin([from_edge,to_edge])].reset_index()
                        # print (gdf_e)
                        network = snkit.Network(gdf_n, gdf_e)
                        network = snkit.network.snap_nodes(network)
                        # if from_edge == to_edge:
                        #     snapped_nodes.append(network.nodes[network.nodes.id.isin([row.from_id,row.to_id])])
                        #     unique_ids += [row.from_id,row.to_id]
                        # else:
                        network = snkit.network.split_edges_at_nodes(network,tolerance=1.0e-3)
                        n = network.nodes
                        # e = gpd.GeoDataFrame(
                        #                 pd.concat(
                        #                         [
                        #                             network.edges,
                        #                             edges[~edges["id"].isin([from_edge,to_edge])]
                        #                         ],
                        #                     axis=0,ignore_index=True),
                        #                 geometry="geometry",crs=gdf_e.crs)
                        e = gpd.GeoDataFrame(network.edges,geometry="geometry",crs=gdf_e.crs)
                        n = gpd.GeoDataFrame(
                                        pd.concat(
                                                [
                                                    n,nodes[nodes["id"].isin(
                                                    list(
                                                        set(
                                                            gdf_e.from_id.values.tolist() + gdf_e.to_id.values.tolist()
                                                            )
                                                        )
                                                    )]],
                                                axis=0,ignore_index=True),
                                        geometry="geometry",crs=nodes.crs)
                        # n = gpd.GeoDataFrame(
                        #                 pd.concat(
                        #                         [
                        #                             n,nodes
                        #                         ],
                        #                         axis=0,ignore_index=True),
                        #                 geometry="geometry",crs=nodes.crs)
                        n = snkit.network.drop_duplicate_geometries(n).reset_index()
                        network = snkit.Network(n, e)
                        network = snkit.network.add_topology(network, id_col="id")
                        e = gpd.GeoDataFrame(
                                        pd.concat(
                                                [
                                                    network.edges,
                                                    edges[~edges["id"].isin([from_edge,to_edge])]
                                                ],
                                            axis=0,ignore_index=True),
                                        geometry="geometry",crs=gdf_e.crs)
                        e["length_m"] = e.geometry.length
                        graph = create_igraph_from_dataframe(e[["from_id","to_id","id","length_m"]])
                        graph_vs = [x['name'] for x in graph.vs]
                        if (row.from_id in graph_vs) and (row.to_id in graph_vs):
                            path = graph.get_shortest_paths(
                                                            row.from_id, 
                                                            row.to_id, 
                                                            weights="length_m", output="epath")[0]
                            if len(path) > 0:
                                path_length = 0
                                for n in path:
                                    path_length += graph.es[n]["length_m"]

                                diff = abs(row.length - path_length) 
                                if diff <= 20 or diff/row.length < 0.3:
                                    snapped_nodes.append(network.nodes[network.nodes.id.isin([row.from_id,row.to_id])])
                                    unique_ids += [row.from_id,row.to_id]
                                    counter = 1
                                elif idx == (len(pairs) - 1):
                                    unadded_bridges.append(row)
                            elif idx == (len(pairs) - 1):
                                unadded_bridges.append(row)

                if len(unadded_bridges) > 0:
                    nonunique_snap = pd.DataFrame(unadded_bridges)
                else:
                    nonunique_snap = pd.DataFrame()
                unique_ids = list(set(unique_ids))
                iteration += 1
        
        snapped_nodes = gpd.GeoDataFrame(
                                pd.concat(
                                    snapped_nodes,axis=0,ignore_index=True),
                                geometry="geometry",crs=bridge_nodes.crs)
        snapped_nodes = snapped_nodes[["id","form_of_road_node","geometry"]]
        snapped_nodes = snkit.network.drop_duplicate_geometries(snapped_nodes).reset_index()
        print (snapped_nodes) 
        snapped_nodes.to_parquet(os.path.join(data_path,"networks","GB_bridge_snapped_nodes.gpq")) 
        if len(nonunique_snap.index) > 0:
            nonmatched_bridges = bridges[bridges["id"].isin(nonunique_snap.id.values.tolist())]
            nonmatched_bridges.to_parquet(os.path.join(data_path,"networks","GB_bridge_unmatched.gpq")) 

        
    step = False
    if step is True:
        b_nodes = gpd.read_parquet(os.path.join(data_path,"networks","GB_bridge_snapped_nodes.gpq"))
        b_nodes.to_file(os.path.join(data_path,"networks","GB_bridge_snapped_nodes.gpkg"),driver="GPKG")

        nonmatched_bridges = gpd.read_parquet(os.path.join(data_path,"networks","GB_bridge_unmatched.gpq"))
        nonmatched_bridges.to_file(os.path.join(data_path,"networks","GB_bridge_unmatched.gpkg"),driver="GPKG")  


    step = False
    if step is True:
        """
        Test of snap and split
        """
        bids = ["europe-latest_28_277319","europe-latest_28_277320"]
        eids = ["DE2AB266-6ADD-475B-96B0-8C9046E8B49F"]
        bridge_nodes = bridge_nodes[bridge_nodes["id"].isin(bids)]
        edges = edges[edges["id"].isin(eids)].reset_index()
        print (bridge_nodes)
        print (edges)
        network = snkit.Network(bridge_nodes, edges)
        print("* Done with network creation")
        network = snkit.network.snap_nodes(network)
        print ('* Done with snapping nodes to edges')
        network = snkit.network.split_edges_at_nodes(network,tolerance=1.0e-3)
        print ('* Done with splitting edges at nodes')

        print (network.edges[["id","geometry"]])
        print (network.nodes)

    step = False
    if step is True:
        new_bridge_nodes = gpd.read_parquet(os.path.join(data_path,"networks","GB_bridge_snapped_nodes.gpq"))
        new_bridge_nodes = new_bridge_nodes.drop_duplicates(subset=["id"],keep="first")
        br_nodes = new_bridge_nodes.copy()
        br_nodes["geometry"] = br_nodes.progress_apply(lambda x:x.geometry.buffer(2.0),axis=1)
        bridge_copy = br_nodes.copy()
        bridge_copy.rename(columns={"id":"b_id"},inplace=True)

        matches = gpd.sjoin(
                            br_nodes[["id","geometry"]],
                            bridge_copy[["b_id","geometry"]],
                            how='inner',
                            predicate='intersects'
                            )
        # print (matches)

        matches = matches[matches["id"] != matches["b_id"]]
        # print (matches)
        graph = create_igraph_from_dataframe(matches[["id","b_id"]],simple=True)
        df = graph.get_edge_dataframe()
        df_vert = graph.get_vertex_dataframe()
        df['source'].replace(df_vert['name'], inplace=True)
        df['target'].replace(df_vert['name'], inplace=True)
        df = df.reset_index()
        # print (df)

        source_targets = list(zip(df.source.values.tolist(),df.target.values.tolist()))
        bridge_pairs = []
        for row in df.itertuples():
            s = row.source
            t = row.target

            s_df = bridges[(bridges.from_id == s) | (bridges.to_id == s)]
            s_ends = [i for i in s_df.from_id.values.tolist() + s_df.to_id.values.tolist() if i != s]

            t_df = bridges[(bridges.from_id == t) | (bridges.to_id == t)]
            t_ends = [i for i in t_df.from_id.values.tolist() + t_df.to_id.values.tolist() if i != t]

            for se in s_ends:
                for te in t_ends:
                    if (se,te) in source_targets or (te,se) in source_targets:
                        edg_s = bridges[
                                        (
                                            (bridges.from_id == s) & (bridges.to_id == se)
                                        ) | (
                                            (bridges.from_id == se) & (bridges.to_id == s)
                                        )].id.values.tolist()
                        edg_t = bridges[
                                        (
                                            (bridges.from_id == t) & (bridges.to_id == te)
                                        ) | (
                                            (bridges.from_id == te) & (bridges.to_id == t)
                                        )].id.values.tolist()
                        bridge_pairs += list(zip(edg_s,edg_t))

        bridge_pairs = pd.DataFrame(bridge_pairs,columns=["id","b_id"])
        # print (bridge_pairs)

        graph = create_igraph_from_dataframe(bridge_pairs[["id","b_id"]],simple=True)
        df = graph.get_edge_dataframe()
        df_vert = graph.get_vertex_dataframe()
        df['source'].replace(df_vert['name'], inplace=True)
        df['target'].replace(df_vert['name'], inplace=True)
        df = df.reset_index()
        # print (df)

        df.rename(columns={"source":"id"},inplace=True)
        df = df[df["id"] != df["target"]]
        exclude_bridges = list(set(df["target"].values.tolist()))
        print (len(exclude_bridges))
        # df = gpd.GeoDataFrame(
        #                     pd.merge(df,bridges[["id","geometry"]],how="left",on=["id"]),
        #                     geometry="geometry",
        #                     crs=bridges.crs
        #                     )
        # df.to_file(os.path.join(data_path,"networks","common_bridges.gpkg"),driver="GPKG")

        e_df = gpd.GeoDataFrame(bridges[bridges["id"].isin(exclude_bridges)],geometry="geometry",crs=bridges.crs)
        print (e_df)
        e_df.to_file(
                    os.path.join(
                        data_path,"networks","exclude_bridges.gpkg"),
                    driver="GPKG")

        include_bridges = bridges[~bridges["id"].isin(exclude_bridges)]
        print (include_bridges)
        bridge_ids = list(set(include_bridges.from_id.values.tolist() + include_bridges.to_id.values.tolist()))

        new_bridge_nodes = new_bridge_nodes[new_bridge_nodes["id"].isin(bridge_ids)]
        print (new_bridge_nodes)
        
        network = snkit.Network(new_bridge_nodes, edges)
        network = snkit.network.split_edges_at_nodes(network,tolerance=1.0e-3)
        print ('* Done with splitting edges at nodes')
        assert set(network.edges.geometry.type.values) == {"LineString"}

        edges = gpd.GeoDataFrame(network.edges,geometry="geometry",crs=edges.crs)
        bridge_nodes = network.nodes
        nodes = gpd.GeoDataFrame(
                        pd.concat(
                                [bridge_nodes,nodes],axis=0,ignore_index=True),
                        geometry="geometry",crs=f"EPSG:{epsg_meters}")
        nodes = snkit.network.drop_duplicate_geometries(nodes).reset_index()
        network = snkit.Network(nodes,edges)
        network = snkit.network.add_topology(network, id_col="id")
        edges = network.edges
        edges["id"] = edges.progress_apply(lambda x:f"roade_{x.name}",axis=1)
        edges = edges.set_crs(epsg=epsg_meters)
        edges.to_parquet(os.path.join(data_path,"networks","GB_roads_edges_test.gpq"))

        nodes = network.nodes
        nodes = nodes.set_crs(epsg=epsg_meters)
        node_ids = list(set(edges.from_id.values.tolist() + edges.to_id.values.tolist()))
        nodes = nodes[nodes["id"].isin(node_ids)] 
        nodes.to_parquet(os.path.join(data_path,"networks","GB_roads_nodes_test.gpq"))

    step = False
    if step is True:
        edges = gpd.read_parquet(os.path.join(data_path,"networks","GB_roads_edges_test.gpq"))
        edges.to_file(os.path.join(data_path,"networks","GB_roads_test.gpkg"),layer="edges",driver="GPKG")
        nodes = gpd.read_parquet(os.path.join(data_path,"networks","GB_roads_nodes_test.gpq"))
        nodes.to_file(os.path.join(data_path,"networks","GB_roads_test.gpkg"),layer="nodes",driver="GPKG")

    step = False
    if step is True:
        edges = gpd.read_parquet(os.path.join(data_path,"networks","GB_roads_edges_test.gpq"))
        nodes = gpd.read_parquet(os.path.join(data_path,"networks","GB_roads_nodes_test.gpq"))
        edges["length_m"] = edges.geometry.length
        print (edges)
        print (nodes)

        graph = create_igraph_from_dataframe(edges[["from_id","to_id","id","length_m"]])
        node_ids = list(set(edges.from_id.values.tolist() + edges.to_id.values.tolist()))
        print (len(node_ids))
        roads_with_bridges = []
        no_of_bridges = len(bridges.index)
        for row in bridges.itertuples():
            bid = row.id
            source = row.from_id
            target = row.to_id
            bridge_length = row.length
            if source in node_ids and target in node_ids:
                path = graph.get_shortest_paths(
                                            source, 
                                            target, 
                                            weights="length_m", output="epath")[0]
                edge_path = []
                path_length = 0
                if path:
                    for n in path:
                        edge_path.append(graph.es[n]["id"])
                        path_length += graph.es[n]["length_m"]
                
                roads_with_bridges.append((bid,bridge_length,path_length,1,edge_path))
            print (f"* Done with bridge {row.Index + 1} out of {no_of_bridges}")

        roads_with_bridges = pd.DataFrame(
                                        roads_with_bridges,
                                        columns=[
                                                    "id","bridge_length",
                                                    "route_length",
                                                    "has_bridge",
                                                    "road_edges"
                                                ]
                                            )

        roads_with_bridges.to_parquet(os.path.join(data_path,"networks","GB_bridge_matches.pq"))
        roads_with_bridges[
                    [
                        "id",
                        "bridge_length",
                        "route_length"]
                    ].to_csv(os.path.join(data_path,"networks","GB_bridge_matches.csv"),index=False)

    step = True
    if step is True:
        m_bridges = pd.read_parquet(os.path.join(data_path,"networks","GB_bridge_matches.pq"))
        m_bridges["length_diff"] = abs(m_bridges["route_length"] - m_bridges["bridge_length"])
        m_bridges["length_ratio"] = m_bridges["length_diff"]/m_bridges["bridge_length"]

        m_bridges = m_bridges[(m_bridges["length_diff"] <= 20) | (m_bridges["length_ratio"] <= 0.3)]
        print (m_bridges[["id","road_edges"]])
        m_bridges = gpd.GeoDataFrame(
                            pd.merge(
                                m_bridges,bridges[["id","tag_highway","geometry"]],
                                how="left",on=["id"]),
                            geometry="geometry",crs=bridges.crs)

        bridge_edges = m_bridges["road_edges"].values.tolist()
        bridge_edges = list(set([item for sublist in bridge_edges for item in sublist]))

        edges,nodes = components(edges,nodes)

        edges = gpd.read_parquet(os.path.join(data_path,"networks","GB_roads_edges_test.gpq"))
        nodes = gpd.read_parquet(os.path.join(data_path,"networks","GB_roads_nodes_test.gpq"))
        edges["road_bridge"] = np.where(
                                        edges["id"].isin(bridge_edges),
                                        "bridge",
                                        "no"    
                                )
        edges["length"] = edges.geometry.length
        edges["length"] = edges["length"].astype(int)
        print (edges[["id","length"]])
        print (nodes)

        edges,nodes = components(edges,nodes)

        if save_file is True:
            m_bridges[
                        [
                            "id",
                            "tag_highway",
                            "bridge_length",
                            "route_length",
                            "length_diff",
                            "length_ratio",
                            "geometry"]
                        ].to_file(os.path.join(
                                            data_path,"networks","matched_bridges.gpkg"),
                                        driver="GPKG")
            nodes[["id","form_of_road_node","geometry"]].to_parquet(
                        os.path.join(
                            data_path,"networks","GB_road_nodes_with_bridges.gpq")
                        )
            edges[
                    [
                        'from_id','to_id','id', 'fictitious', 'road_classification', 'road_function', 
                        'form_of_way', 'road_classification_number', 
                        'name_1', 'name_1_lang', 'name_2', 'name_2_lang', 
                        'road_structure','road_bridge', 'length', 'length_uom', 
                        'loop', 'primary_route', 'trunk_road', 
                        'e_id','average_toll_cost', 'urban', 
                        'averageWidth', 'minimumWidth', 
                        'max_z', 'min_z', 'mean_z', 'lanes','geometry'
                    ]
                ].to_parquet(
                        os.path.join(
                            data_path,"networks","GB_road_links_with_bridges.gpq")
                        )




    
if __name__ == '__main__':
    CONFIG = load_config() 
    main(CONFIG)