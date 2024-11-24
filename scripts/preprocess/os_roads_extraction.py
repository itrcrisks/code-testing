# Script to extract GML data from OS Mastermap 
# coding: utf-8
import sys
import os
import gzip
import pandas as pd
import geopandas as gpd
from utils import *

def get_chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

def fix_columns():

def main(config):
    data_path = config['paths']['data']
    epsg_meters = 27700
    folder_names = ["MasterMap Highways Network_roads_5731659","MasterMap Highways Network_rami_5731064"]
    file_start_strings ["Highways_Roads","Highways_Rrami"]
    file_names = [
                    [
                        "FerryLink","FerryNode","FerryTerminal",
                        "Road","RoadJunction","RoadLink",
                        "RoadNode","Street"
                    ],
                    [
                        "FerryLink","FerryNode","FerryTerminal",
                        "Hazard","Maintenance","Reinstatement",
                        "AccessRestriction","Dedication",
                        "RestrictionnForVehicles","Road",
                        "RoadJunction","RoadLink","RoadNode",
                        "Street","SpecialDesigation","Structure",
                        "TurnRestriction"
                    ]
                ]
    fix_cols = ["roadName","alternateName"]
    chunks = 20
    for idx,(fld,fst,fnames) in zip(folder_names,file_start_strings,file_names):
        folder_name = os.path.join(
                                data_path,
                                "DAFNI_NIRD",
                                "incoming_data",
                                f"{fld}"
                                )
        for fn in fnames:
            all_names = []
            dfs = []
            for root, dirs, files in os.walk(folder_name):
                for file in files:
                    if file.startswith(f"{fst}_{fn}_") and file.endswith(".gz"):
                        all_names.append(file)

            if len(all_names) > 20:
                all_names = get_chunks(all_names, chunks)
            else:
                for file in all_names:
                    with gzip.open(os.path.join(folder_name,file), 'rb') as f:
                        df = gpd.read_file(f, driver='GML')
                        dfs.append(df)
    
                dfs = pd.concat(dfs,axis=0,ignore_index=True)
                    
            for fx_col in fix_cols:
                if fx_col in dfs.columns.values.tolist():
                    dfs[fx_col] = dfs[fx_col].astype(str)
            
            dfs.to_parquet(os.path.join(folder_name,f"GB_{fn}.pq"))

if __name__ == '__main__':
    CONFIG = load_config() 
    main(CONFIG)