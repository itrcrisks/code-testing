#!/usr/bin/env python
# coding: utf-8
import sys
import os
import pandas as pd
import geopandas as gpd
from utils import *


def main(config):
    incoming_data_path = config['paths']['incoming_data']
    processed_data_path = config['paths']['data']
    
    port_edges = gpd.read_file(os.path.join(incoming_data_path,
                                            "ports",
                                            "africa_ports.gpkg"),
                                        layer="edges")
    port_nodes = gpd.read_file(os.path.join(incoming_data_path,
                                            "ports",
                                            "africa_ports.gpkg"),
                                        layer="nodes")

    print (port_nodes)
    print (port_edges)


if __name__ == '__main__':
    CONFIG = load_config()
    main(CONFIG)