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
    folder_name = os.path.join(
                            data_path,
                            "DAFNI_NIRD",
                            "incoming_data",
                            "MasterMap Highways Network_roads_5731659"
                            )
    file_names = ["FerryLink","FerryNode","FerryTerminal","Road","RoadJunction","RoadLink","RoadNode","Street"]
    for fn in file_names:
        dfs = []
        for root, dirs, files in os.walk(folder_name):
            for file in files:
                if file.startswith(f"Highways_Roads_{fn}") and file.endswith(".gz"):
                    dfs.append(pd.read_csv(os.path.join(folder_name,file),compression='gzip'))
        
        dfs = pd.concat(dfs,axis=0,ignore_index=True)

        dfs.to_parquet(os.path.join(folder_name,f"GB_{fn}.pq"))

if __name__ == '__main__':
    CONFIG = load_config() 
    main(CONFIG)