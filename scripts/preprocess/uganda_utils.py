"""Functions for preprocessing road data
    WILL MODIFY LATER
"""
import logging
import warnings
import sys
import os
import json
import snkit
import numpy as np
import pandas as pd
import igraph as ig
import networkx
import geopandas as gpd


def components(edges,nodes,
                node_id_column="id",edge_id_column="id",
                from_node_column="from_id",to_node_column="to_id"):
    G = networkx.Graph()
    G.add_nodes_from(
        (getattr(n, node_id_column), {"geometry": n.geometry}) for n in nodes.itertuples()
    )
    G.add_edges_from(
        (getattr(e,from_node_column), getattr(e,to_node_column), 
            {edge_id_column: getattr(e,edge_id_column), "geometry": e.geometry})
        for e in edges.itertuples()
    )
    components = networkx.connected_components(G)
    for num, c in enumerate(components):
        print(f"Component {num} has {len(c)} nodes")
        edges.loc[(edges[from_node_column].isin(c) | edges[to_node_column].isin(c)), "component"] = num
        nodes.loc[nodes[node_id_column].isin(c), "component"] = num

    return edges, nodes

def add_node_degree(edges_dataframe,nodes_dataframe):
    degree_df = edges_dataframe[["from_id","to_id"]].stack().value_counts().rename_axis('id').reset_index(name='degree')
    nodes_dataframe = pd.merge(nodes_dataframe,degree_df,how="left",on=["id"])

    nodes_crs = nodes_dataframe.crs
    return gpd.GeoDataFrame(nodes_dataframe,geometry="geometry",crs=nodes_crs)


def create_network(
    edges: gpd.GeoDataFrame,
    nodes: gpd.GeoDataFrame = None,
    id_prefix: str = ""
    ) -> snkit.network.Network:
    """
    Create snkit network from edges and (optional) nodes and clean the result.

    Arguments:
        edges (gpd.GeoDataFrame): Expected to contain geometry column of linestrings
        nodes (gpd.GeoDataFrame): Optional nodes to include. snkit will try to snap to edges

    Returns:
        snkit.network.Network: Built network
    """

    logging.info("Starting network creation")

    # drop edges with no geometry
    empty_idx = edges.geometry.apply(lambda e: e is None or e.is_empty)
    if empty_idx.sum():
        empty_edges = edges[empty_idx]
        logging.info(f"Found {len(empty_edges)} empty edges.")
        logging.info(empty_edges)
        edges = edges[~empty_idx].copy()

    logging.info("Creating network")
    network = snkit.Network(nodes, edges)

    logging.info("Splitting multilines")
    network = snkit.network.split_multilinestrings(network)

    if nodes is not None and not nodes.empty:
        # check we have only point nodes
        assert set(network.nodes.geometry.type.values) == {"Point"}

        logging.info("Dropping duplicate geometries")
        # silence shapely.ops.split ShapelyDeprecationWarning regarding:
        # shapley.ops.split failure on split by empty geometry collection
        # this is currently caught by snkit.network.split_edge_at_points,
        # but won't be for shapely==2.0
        warnings.filterwarnings(
            "ignore",
            message=(
                ".*GeometryTypeError will derive from ShapelyError "
                "and not TypeError or ValueError in Shapely 2.0*"
            )
        )
        network.nodes = snkit.network.drop_duplicate_geometries(network.nodes)

        logging.info("Snapping nodes to edges")
        network = snkit.network.snap_nodes(network)

    logging.info("Adding endpoints")
    network = snkit.network.add_endpoints(network)

    # logging.info("Splitting edges at nodes")
    network = snkit.network.split_edges_at_nodes(network)

    # check we have only linestrings
    assert set(network.edges.geometry.type.values) == {"LineString"}

    logging.info("Renaming nodes and edges")
    network = snkit.network.add_ids(network, edge_prefix=id_prefix, node_prefix=id_prefix)

    logging.info("Creating network topology")
    network = snkit.network.add_topology(network, id_col="id")

    return network

def create_igraph_from_dataframe(graph_dataframe, directed=False, simple=False):
    graph = ig.Graph.TupleList(
        graph_dataframe.itertuples(index=False),
        edge_attrs=list(graph_dataframe.columns)[2:],
        directed=directed
    )
    if simple:
        graph.simplify()

    es, vs, simple = graph.es, graph.vs, graph.is_simple()
    d = "directed" if directed else "undirected"
    s = "simple" if simple else "multi"

    return graph
