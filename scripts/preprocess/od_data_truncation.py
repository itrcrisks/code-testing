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

def truncate_by_threshold(flow_dataframe, flow_column='flux', threshold=.99):
    print(f"Truncating paths with threshold {threshold * 100:.0f}%.")
    flows_sorted = flow_dataframe.reset_index(drop=True).sort_values(by=flow_column, ascending=False)
    fluxes_sorted = flows_sorted[flow_column]
    total_flux = fluxes_sorted.sum()
    flux_percentiles = fluxes_sorted.cumsum() / total_flux
    excess = flux_percentiles[flux_percentiles >= threshold]
    cutoff = excess.idxmin()
    keep = flux_percentiles[flux_percentiles <= threshold].index
    flows_truncated = flows_sorted.loc[keep, :]
    print(f"Number of paths before: {len(flows_sorted):,.0f}.")
    print(f"Number of paths after: {len(flows_truncated):,.0f}.")
    return flows_truncated

def main(config):
    incoming_data_path = config['paths']['incoming_data']
    data_path = config['paths']['data']
    
    od_df = pd.read_csv(os.path.join(data_path,"census_datasets","od_gb_oa_2021.csv"))
    od_df = od_df[od_df["Car21"] > 0]
    
    od_df = od_df.sort_values(by=["Car21"],ascending=False)
    od_df["car_total"] = od_df.groupby(["Area of usual residence"])["Car21"].transform("sum")
    od_df["car_cumsum"] = od_df.groupby(["Area of usual residence"])["Car21"].transform("cumsum")
    od_df["fraction"] = od_df["car_cumsum"]/od_df["car_total"]
    od_df["count"] = od_df.groupby(["Area of usual residence"])["Area of usual residence"].transform("count")
    print (od_df)
    before = od_df["Car21"].sum()
    print (before)
    od_df = od_df[(od_df["fraction"] <= 0.90) | (od_df["count"] == 1.0)]
    print (od_df)
    after = od_df["Car21"].sum()
    print (after)
    print (before-after)
    print (1.0*after/before)


    # origins = od_df["Area of usual residence"].values.tolist()
    # new_od_df = []
    # for o in origins:
    #     df = od_df[od_df["Area of usual residence"] == o]
    #     if len(df.index) == 1:
    #         new_od_df.append(df)
    #     else:
    #         df = truncate_by_threshold(df, flow_column='Car21', threshold=.90)
    #         new_od_df.append(df)

    # new_od_df = pd.concat(new_od_df,axis=0,inplace=True)
    # print (new_od_df)


    
if __name__ == '__main__':
    CONFIG = load_config() 
    main(CONFIG)