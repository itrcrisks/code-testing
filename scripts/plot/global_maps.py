"""Road network risks and adaptation maps
"""
import os
import pandas as pd
import geopandas as gpd
from shapely import wkt
from ast import literal_eval
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from map_plotting_utils import *
from tqdm import tqdm
tqdm.pandas()


def main(config):
    data_path = config['paths']['data']
    output_path = config['paths']['results']
    figure_path = config['paths']['figures']

    figures = os.path.join(figure_path,"figures")
    if os.path.exists(figures) == False:
        os.mkdir(figures)

    """Identify countries from the SDG index
    """
    sdr_columns = ["Country Code ISO3","Country","Goal 9 Dash"]
    rename_columns = ["ISO_A3_EH","Country","SDG_9_level"]
    sdr_level_colors = ["green","yellow","orange","red"]
    sdr_level_progress = ["SDG achieved",
                            "Challenges remain",
                            "Significant challenges remain",
                            "Major challenges remain"]
    original_colors = sdr_level_colors + ["grey"]
    new_colors = ["#000000","#525252","#969696","#d9d9d9","#f0f0f0"]
    progress_labels = sdr_level_progress + ["No data"]
    infra_sectors = ["Transport","Energy","Water"]
    # infra_colors = ["#fd8d3c","#41ab5d","#3690c0"]
    infra_colors = {
                    "Transport":"#fdae61",
                    "Energy":"#ffffbf",
                    "Water": "#4575b4",
                    "Waste":"#abd9e9",
                    "Extractives": "#d73027"
                    }

    admin_gdf = gpd.read_file(os.path.join(
                    data_path,
                    'ne_10m_admin_0_countries',
                    'ne_10m_admin_0_countries.shp'),
                    encoding="utf-8")
    admin_gdf = admin_gdf[["ISO_A3_EH","ADMIN","geometry"]]

    sdg_title = "SDG 9 'infrastructure' index score"
    sdr_df = pd.read_excel(
                    os.path.join(data_path,
                            "Raw data","SDR2023-data_RAW.xlsx"),
                    sheet_name="SDR 2023 Data")
    sdr_df = sdr_df[sdr_columns]
    sdr_df.columns = rename_columns
    sdr_df["identify_country"] = sdr_df.progress_apply(lambda x:0 if x["ISO_A3_EH"][0] == "_" else 1, axis=1)
    sdr_df["sdr_progress"] = ''
    for idx,(lc,lp) in enumerate(zip(sdr_level_colors,sdr_level_progress)):
        sdr_df.loc[sdr_df["SDG_9_level"] == lc,"sdr_progress"] = lp

    sdr_df = sdr_df[sdr_df["identify_country"] == 1]
    admin_gdf = gpd.read_file(os.path.join(
                    data_path,
                    'ne_10m_admin_0_countries',
                    'ne_10m_admin_0_countries.shp'),
                    encoding="utf-8")
    admin_gdf = admin_gdf[["ISO_A3_EH","ADMIN","geometry"]]
    sdr_df = pd.merge(admin_gdf,sdr_df,how="left",on=["ISO_A3_EH"])
    sdr_df["identify_country"] = sdr_df["identify_country"].fillna(0)
    sdr_df["sdr_progress"] = sdr_df["sdr_progress"].fillna("No data")
    sdr_df["SDG_9_level"] = sdr_df["SDG_9_level"].fillna("grey")

    """
    Clean the data on the conflict typology linked to infrastructure classes
    """
    infra_typology_df = pd.read_excel(
                            os.path.join(
                                data_path,"EJATLAS_look-up.xlsx"),
                            sheet_name="Sheet1")
    chars = [" ","(",")",",",".","/","--","--"]
    infra_typology_df["type"] = infra_typology_df["type"].str.lower()
    for ch in chars:
        infra_typology_df["type"] = infra_typology_df["type"].str.replace(ch,"-")

    infra_typology_df["type"] = infra_typology_df["type"].str.strip("-")
    infra_typology_df["UDHR_type_new"] = infra_typology_df["UDHR_type_new"].fillna("Other")

    """
    EJATLAS data
    """
    ejatlas_columns = ["id","country","type","geometry"]
    ejatlas_df = pd.read_csv(os.path.join(data_path,
                            "Raw data",
                            "ejatlas_org_combined.csv"))[ejatlas_columns]
    ejatlas_df["type"] = ejatlas_df["type"].fillna("['Unkown']")
    ejatlas_df["type"] = ejatlas_df.progress_apply(lambda x: literal_eval(x["type"]),axis=1)
    ej_df = []
    for row in ejatlas_df.itertuples():
        id_v = row.id
        t_v = row.type
        ej_df += list(zip([id_v]*len(t_v),t_v))

    ej_df = pd.DataFrame(ej_df,columns=["id","type"])
    ej_df = pd.merge(ej_df,ejatlas_df[["id","country","geometry"]],how="left",on=["id"])
    ej_df = pd.merge(ej_df,infra_typology_df[["type","UDHR_type_new"]],how="left",on=["type"])
    ej_df["UDHR_type_new"] = ej_df["UDHR_type_new"].fillna("Unkown")        
    ej_df['geometry'] = ej_df['geometry'].apply(wkt.loads)
    ej_df = gpd.GeoDataFrame(ej_df,geometry="geometry", crs='epsg:4326')

    """
    Testing global map plot
    """
    
    ax_proj = get_projection(epsg=4326)
    fig, ax_plots = plt.subplots(1,1,
                    subplot_kw={'projection': ax_proj},
                    figsize=(12,6),
                    dpi=500)
    ax = plot_global_basemap(ax_plots)
    plt.tight_layout()
    save_fig(os.path.join(figures,"global_map_template.png"))
    plt.close()

    """Replicate the SDG color map
    """
    ax_proj = get_projection(epsg=4326)
    fig, ax_plots = plt.subplots(1,1,
                    subplot_kw={'projection': ax_proj},
                    figsize=(12,6),
                    dpi=500)
    ax = plot_global_basemap(ax_plots)
    legend_handles = []
    for idx,(c,l) in enumerate(zip(original_colors,progress_labels)):
        ax = plot_polygon_assets(
                        ax,
                        sdr_df[sdr_df["SDG_9_level"] == c],
                        l,
                        facecolor=c)
        legend_handles.append(mpatches.Patch(color=c,
                                        label=l))
    plt.legend(
                handles=legend_handles,
                fontsize=12,
                title=sdg_title,
                title_fontproperties={'weight':'bold'},
                loc='lower left',frameon=False)
    plt.tight_layout()
    save_fig(os.path.join(figures,"sdg_colors_original.png"))
    plt.close()

    """Create the SDG color map as per new color scheme for paper
    """
    ax_proj = get_projection(epsg=4326)
    fig, ax_plots = plt.subplots(1,1,
                    subplot_kw={'projection': ax_proj},
                    figsize=(12,6),
                    dpi=500)
    ax = plot_global_basemap(ax_plots)
    legend_handles = []
    for idx,(c,nc,l) in enumerate(zip(original_colors,new_colors,progress_labels)):
        ax = plot_polygon_assets(
                        ax,
                        sdr_df[sdr_df["SDG_9_level"] == c],
                        l,
                        facecolor=nc)
        legend_handles.append(mpatches.Patch(color=nc,
                                        label=l))
    plt.legend(
                handles=legend_handles,
                fontsize=12,
                title=sdg_title,
                title_fontproperties={'weight':'bold'},
                loc='lower left',frameon=False)
    plt.tight_layout()
    save_fig(os.path.join(figures,"sdg_colors_new.png"))
    plt.close()

    """Create the SDG color map as per new color scheme for paper
        And now add the points with the categorisations
    """
    infra_title = "\n$\\bf{Infrastructure \, related \, conflicts}$"
    sdg_title = "\n$\\bf{SDG \, 9 \, 'infrastructure' \, index \, score}$"
    ax_proj = get_projection(epsg=4326)
    fig, ax_plots = plt.subplots(1,1,
                    subplot_kw={'projection': ax_proj},
                    figsize=(12,6),
                    dpi=500)
    ax = plot_global_basemap(ax_plots,include_scalebar=False)
    legend_handles = []
    legend_handles.append(plt.plot([],[],
                                    color="none",
                                    label="\n$\\bf{Infrastructure \, related \, conflicts}$")[0])
    # for idx,(infra,color) in enumerate(zip(infra_sectors,infra_colors)):
    for infra,color in infra_colors.items():
        nodes = ej_df[ej_df["UDHR_type_new"] == infra]
        ax = plot_point_assets(ax,nodes,infra,colors=color,size=9,marker='.')
        legend_handles.append(plt.plot([],[],
                                    marker=".", 
                                    ms=9, 
                                    ls="",
                                    color=color,
                                    label=infra)[0])

    
    legend_handles.append(mpatches.Patch(color="none",
                                        label="\n$\\bf{SDG \, 9 \, 'infrastructure' \, index \, score}$"))
    for idx,(c,nc,l) in enumerate(zip(original_colors,new_colors,progress_labels)):
        ax = plot_polygon_assets(
                        ax,
                        sdr_df[sdr_df["SDG_9_level"] == c],
                        l,
                        facecolor=nc,edgecolor=nc)
        legend_handles.append(mpatches.Patch(color=nc,
                                        label=l))
    leg = ax.legend(
            handles=legend_handles, 
            fontsize=9, 
            loc='lower left',
            frameon=False)

    ## Move titles to the left 
    for item, label in zip(leg.legend_handles, leg.texts):
        if label._text  in [infra_title,sdg_title]:
            width=item.get_window_extent(fig.canvas.get_renderer()).width
            label.set_ha('left')
            label.set_position((-4*width,0))

    plt.tight_layout()
    save_fig(os.path.join(figures,"sdg_infra_conflicts.png"))
    plt.close()
    

if __name__ == '__main__':
    CONFIG = load_config()
    main(CONFIG)
