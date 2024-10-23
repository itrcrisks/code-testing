"""Road network risks and adaptation maps
"""
import os
import pandas as pd
import geopandas as gpd
from scipy import integrate
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from map_plotting_utils import *
from tqdm import tqdm
tqdm.pandas()

def risks(dataframe,index_columns,probabilities,
            expected_risk_column,
            flood_protection_period=0,flood_protection_name=None):
    
    """
    Organise the dataframe to pivot with respect to index columns
    Find the expected risks
    """
    if flood_protection_name is None and flood_protection_period == 0:
        # When there is no flood protection at all
        expected_risk_column = f"{expected_risk_column}_undefended"
        probability_columns = [str(p) for p in probabilities]
        
    elif flood_protection_period > 0:
        if flood_protection_name is None:
            expected_risk_column = f"{expected_risk_column}_{flood_protection_period}_year_protection"
        else:
            expected_risk_column = f"{expected_risk_column}_{flood_protection_name}"
        
        probabilities = [pr for pr in probabilities if pr <= 1.0/flood_protection_period]
        probability_columns = [str(p) for p in probabilities]
    else:
        # When there is no flood protection at all
        expected_risk_column = f"{expected_risk_column}_{flood_protection_name}"
        probability_columns = [str(p) for p in probabilities]
        
    dataframe.columns = dataframe.columns.astype(str)
    dataframe[expected_risk_column] = list(integrate.trapz(dataframe[probability_columns].to_numpy(),
                                            np.array([probabilities*len(dataframe.index)]).reshape(dataframe[probability_columns].shape)))

    # dataframe[expected_risk_column] = list(integrate.simpson(dataframe[probability_columns].to_numpy(),
    #                                         np.array([probabilities*len(dataframe.index)]).reshape(dataframe[probability_columns].shape)))
    
    # dataframe = dataframe[index_columns + [expected_risk_column]].set_index(index_cols)
    return dataframe[index_columns + [expected_risk_column]].set_index(index_columns)

def risk_estimations(hazard_data_details,hazard_index_columns,hazard_dataframe,asset_id,damage_type,flood_protection_name):
    hazard_data_details = hazard_data_details.set_index(hazard_index_columns)
    haz_index_vals = list(set(hazard_data_details.index.values.tolist()))
    expected_damages = []
    for hc in haz_index_vals:
        haz_df = hazard_data_details[hazard_data_details.index == hc]
        haz_cols, haz_rps = map(list,list(zip(*sorted(
                                    list(zip(haz_df.key.values.tolist(),
                                    haz_df.rp.values.tolist()
                                    )),key=lambda x:x[-1],reverse=True))))
        
        # haz_cols = [f"{c}_mod" for c in haz_cols]
        haz_prob = [1.0/rp for rp in haz_rps]
        damages = hazard_dataframe[[asset_id,flood_protection_name] + haz_cols] 
        damages.columns = [asset_id,flood_protection_name] + haz_prob
        if min(haz_prob) > 0:
            damages[0] = damages[min(haz_prob)]
            haz_prob = [0] + haz_prob
        if max(haz_prob) < 1:
            # damages[1] = damages[max(haz_prob)]
            damages[1] = 0
            haz_prob += [1]
        # hz_st = '_'.join([str(h_c) for h_c in hc])
        # damage_col = f"{damage_type}_{hz_st}"
        damage_col = damage_type
        damages = damages[[asset_id,flood_protection_name] + haz_prob]
        expected_damage_df = risks(damages,[asset_id, flood_protection_name],haz_prob,
                                damage_col,
                                flood_protection_name=flood_protection_name 
                                )
        # print (expected_damage_df)
        expected_damages.append(expected_damage_df)
        del expected_damage_df

    expected_damages = pd.concat(expected_damages,axis=1)
    return expected_damages.reset_index() 

def main(config):
    data_path = config['paths']['data']
    output_path = config['paths']['results']
    figure_path = config['paths']['figures']

    figures = os.path.join(figure_path,"p3_data","figures")

    epsg_meters = 27700
    clip_layer = False
    if clip_layer is True:
        gb_boundary = gpd.read_file(
                        os.path.join(
                                data_path,
                                "flood_areas",
                                "gb_boundary",
                                "GB_Boundary.shp"
                                )
                            )
        gb_boundary = gb_boundary.to_crs(epsg=epsg_meters)
        hydro_basins = gpd.read_file(os.path.join(
                                data_path,
                                "flood_areas",
                                "eu",
                                "hybas_eu_lev08_v1c",
                                "hybas_eu_lev08_v1c.shp"
                                )
                            )
        hydro_basins = hydro_basins.to_crs(epsg=epsg_meters)

        gb_basins = gpd.clip(hydro_basins,gb_boundary)
        gb_basins.to_file(os.path.join(
                                data_path,
                                "flood_areas",
                                "gb_boundary",
                                "gb_basins.gkpg"
                                ),driver="GPKG")
    gb_basins = gpd.read_file(os.path.join(
                                    data_path,
                                    "flood_areas",
                                    "gb_boundary",
                                    "gb_basins.gkpg"
                                    ))

    return_periods = [20,100,200]
    asset_types = ["Elect","Station","LX","Track","Total damage cost","passengers_delay_cost"]
    asset_labels = ["Traction substations","Rail stations",
                    "Level Crossings","Track length (km)",
                    "Direct damage costs","Passengers delay costs"]
    asset_bbox = [(1.26,0.89),(1.18,0.93),(1.21,0.86),(1.23,0.86),(1.25,0.85),(1.28,0.85)]
    all_colors = ["#ffffcc","#ffeda0","#fed976","#feb24c","#fd8d3c","#fc4e2a","#e31a1c","#bd0026","#800026"]
    all_colors = ["#fed976","#feb24c","#fd8d3c","#fc4e2a","#e31a1c","#bd0026","#800026"]
    no_exposure_color = "#d9d9d9"
    exp_df = []
    for rp in return_periods:
        exposures = pd.read_excel(
                        os.path.join(
                            data_path,
                            "p3_data",
                            "For_Figure_7_analysis_plots.xlsx"),
                        sheet_name=f"all_catchment_EAD_EAEL_{rp}",
                        header=[0,1])
        cols = ["HYBAS_ID"] + [f"{c[0]}_{c[1]}_{rp}" for c in exposures.columns.values.tolist()[1:]]
        exposures.columns = cols
        exp_df.append(exposures.set_index("HYBAS_ID"))
    
    exp_df = pd.concat(exp_df,axis=1)
    exp_df = exp_df.reset_index()
    gb_exposures = pd.merge(gb_basins,exp_df,how="left",on=["HYBAS_ID"]).fillna(0)
    # print (gb_exposures.columns.values.tolist())

    ead_eael_df = []
    cost_cols = ["Total damage cost","passengers_delay_cost"]
    risk_cols = ["EAD","EAEL"]
    for idx,(col,rcol) in enumerate(zip(cost_cols,risk_cols)):
        vals = []
        for rp in return_periods:
            vals.append((rp,f"Cost_{col}_{rp}"))
            if col == "passengers_delay_cost":
                gb_exposures[f"Cost_{col}_{rp}"] = gb_exposures[f"Cost_{col}_{rp}"]/60.0

        haz_det_df = pd.DataFrame(vals,columns=["rp","key"])
        haz_det_df["hazard"] = "flood"
        gb_exposures["protection"] = "river"
        risk_df = risk_estimations(
                        haz_det_df,
                        ["hazard"],gb_exposures,"HYBAS_ID",rcol,"protection")
        ead_eael_df.append(risk_df.set_index("HYBAS_ID"))
    
    ead_eael_df = pd.concat(ead_eael_df,axis=1).reset_index()
    print (ead_eael_df)
    ead_eael_df.to_csv("EAD_EAEL.csv",index=False)
    ead_eael_df["Risk"] = ead_eael_df["EAD_protection"] + ead_eael_df["EAEL_protection"]
    gb_exposures = pd.merge(gb_exposures,ead_eael_df,how="left",on=["HYBAS_ID"])

    plot_maps = True
    if plot_maps is True:
        # gb_basins = gpd.read_file(os.path.join(
        #                             data_path,
        #                             "flood_areas",
        #                             "gb_boundary",
        #                             "gb_basins.gkpg"
        #                             ))

        # return_periods = [20,100,200]
        # asset_types = ["Elect","Station","LX","Track","Total damage cost","passengers_delay_cost"]
        # asset_labels = ["Traction substations","Rail stations",
        #                 "Level Crossings","Track length (km)",
        #                 "Direct damage costs","Passengers delay costs"]
        # all_colors = ["#ffffcc","#ffeda0","#fed976","#feb24c","#fd8d3c","#fc4e2a","#e31a1c","#bd0026","#800026"]
        # all_colors = ["#fed976","#feb24c","#fd8d3c","#fc4e2a","#e31a1c","#bd0026","#800026"]
        # no_exposure_color = "#d9d9d9"
        # exp_df = []
        # for rp in return_periods:
        #     exposures = pd.read_excel(
        #                     os.path.join(
        #                         data_path,
        #                         "p3_data",
        #                         "For_Figure_7_analysis_plots.xlsx"),
        #                     sheet_name=f"all_catchment_EAD_EAEL_{rp}",
        #                     header=[0,1])
        #     cols = ["HYBAS_ID"] + [f"{c[0]}_{c[1]}_{rp}" for c in exposures.columns.values.tolist()[1:]]
        #     exposures.columns = cols
        #     exp_df.append(exposures.set_index("HYBAS_ID"))
        
        # exp_df = pd.concat(exp_df,axis=1)
        # exp_df = exp_df.reset_index()
        # gb_exposures = pd.merge(gb_basins,exp_df,how="left",on=["HYBAS_ID"]).fillna(0)

        ax_proj = ccrs.epsg(epsg_meters)

        figure_texts = ['(a) ','(b) ','(c) ']
        for idx,(asset,asset_label,ab) in enumerate(zip(asset_types,asset_labels,asset_bbox)):
            if asset in ["Total damage cost","passengers_delay_cost"]:
                columns = [f"Cost_{asset}_{rp}" for rp in return_periods]
                gb_exposures[columns] = gb_exposures[columns]/1.0e3
            else:
                columns = [f"failed assets_{asset}_{rp}" for rp in return_periods]
            # if asset == "passengers_delay_cost":
            #     gb_exposures[columns] = gb_exposures[columns]/60.0
            vals = []
            for c in columns:
                vals += gb_exposures[gb_exposures[c]>0][c].tolist()
            
            fig, ax_plots = plt.subplots(1,3,
                    subplot_kw={'projection': ax_proj},
                    figsize=(17,8),
                    dpi=500)
            ax_plots = ax_plots.flatten()
            for j, (rp,c) in enumerate(zip(return_periods,columns)):
                ax = ax_plots[j]
                total = gb_exposures[c].sum()
                legend_handles = []
                if (asset not in ["Total damage cost","passengers_delay_cost"]) and (max(vals) < len(all_colors)):
                    vals_by_range = sorted(list(set(vals)))
                    for v in vals_by_range:
                        ax = plot_polygon_assets(ax,
                                gb_exposures[gb_exposures[c]==v],"Exposed",
                                zorder=2 + int(v),
                                facecolor=all_colors[int(v)],
                                projection=ax_proj)
                        legend_handles.append(mpatches.Patch(color=all_colors[int(v)],
                                                label=v))

                else:
                    if asset not in ["Total damage cost","passengers_delay_cost"]:
                        vals_by_range = generate_weight_bins(vals,n_steps=7)
                    else:
                        vals_by_range = generate_weight_bins(vals,n_steps=7,interpolation='fisher-jenks')
                    for (i, ((nmin, nmax), width)) in enumerate(vals_by_range.items()):
                        ax = plot_polygon_assets(ax,
                                gb_exposures[(gb_exposures[c] >= nmin) & (gb_exposures[c] < nmax+1)],"Exposed",
                                zorder=2 + i + 1,
                                facecolor=all_colors[i],
                                projection=ax_proj)
                        legend_handles.append(mpatches.Patch(color=all_colors[i],
                                                label=f"{int(nmin)} - {int(nmax)}"))
                ax = plot_polygon_assets(ax,
                            gb_exposures[gb_exposures[c]==0],"No exposure",
                            zorder=2,
                            facecolor=no_exposure_color,
                            projection=ax_proj)
                legend_handles.append(mpatches.Patch(color=no_exposure_color,
                                                label="No exposure"))
                ax.text(
                        0.05,
                        0.95,
                        f"{figure_texts[j]} {rp}-year Return period",
                        horizontalalignment='left',
                        transform=ax.transAxes,
                        size=12,
                        weight='bold')
                if asset == "Track":
                    st = f"Total Track length flooded = {round(total,0)} km"
                elif asset in ["Total damage cost"]:
                    st = f"Total {asset_label} \n = £ {round(1e-3*total,0)} million"
                elif asset in ["passengers_delay_cost"]:
                    st = f"Total {asset_label} \n = £ {round(1e-3*total,2)} million"
                else:
                    st = f"Total {asset_label} flooded = {round(total,0)}"
                ax.text(
                        0.10,
                        0.85,
                        st,
                        horizontalalignment='left',
                        transform=ax.transAxes,
                        size=12,
                        weight='bold')
            title = f"Flooded {asset_label}"
            fontsize = 12
            bbox_to_anchor=ab
            if asset in ["Total damage cost","passengers_delay_cost"]:
                title = f"Flooded {asset_label} \n (£ '000)"
            plt.legend(
                    handles=legend_handles,
                    fontsize=fontsize,
                    title=title,
                    title_fontproperties={'weight':'bold'},
                    bbox_to_anchor=bbox_to_anchor,loc="center",frameon=False)
            plt.tight_layout()
            save_fig(os.path.join(figures,f"{asset.replace(' ','_')}_exposures.png"))


    # Plot of AEP
    plot_aep = False
    if plot_aep is True:
        mpl.style.use('ggplot')
        mpl.rcParams['font.size'] = 10.
        mpl.rcParams['font.family'] = 'tahoma'
        mpl.rcParams['axes.labelsize'] = 10.
        mpl.rcParams['xtick.labelsize'] = 12.
        mpl.rcParams['ytick.labelsize'] = 12.
        figure_texts = ['(a)','(b)']
        aep_df = pd.read_excel(
                    os.path.join(
                        data_path,
                        "p3_data",
                        "For Figure 8 - 9_plot.xlsx"),
                    sheet_name="Sheet1")
        y_column = "Joint AEP"
        x_columns = ["no. of trains delayed","no. of passengers delayed"]
        x_labels = ["No. of trains delayed","No. of passengers disrupted (million)"]
        divisor = [1,1.0e6]
        colors = ["#252525","#000000"]
        ytick_locs = [0.05,2.5e-3,2.5e-4,1.25e-5,1.25e-6,1.25e-7,1.25e-8,2.5e-9,2.5e-10,3e-11]
        ytick_lables = ["0.05","2.5e-3","2.5e-4","1.25e-5","1.25e-6","1.25e-7","1.25e-8","2.5e-9","2.5e-10","3e-11"]
        fig, ax_plots = plt.subplots(1,2,
                    sharey=True,
                    figsize=(16,8),
                    dpi=500)
        ax_plots = ax_plots.flatten()
        for j, (x,xl,c,d) in enumerate(zip(x_columns,x_labels,colors,divisor)):
            ax = ax_plots[j]
            ax.scatter(
                    aep_df[x]/d,
                    aep_df[y_column],
                    marker='o',
                    s=20,
                    color=c,
                    )
            ax.set_xlabel(f"{figure_texts[j]} {xl}",fontsize=16,fontweight='bold')
        
        ax_plots[0].set_ylabel("Joint annual exceedance probability (AEP)",fontsize=16,fontweight='bold')
        plt.yscale('log')
        plt.yticks(ytick_locs,ytick_lables)
        plt.tight_layout()
        save_fig(os.path.join(figures,"joint_disruptions.png"))

    plot_ead_eael = False
    if plot_ead_eael is True:
        figure_texts = ['(a) ','(b) ','(c) ']
        asset_types = columns = ["EAD_protection","EAEL_protection","Risk"]
        asset_labels = ["EAD", "EAEL", "EAL"]
        gb_exposures[columns] = gb_exposures[columns]/1.0e3
        vals = []
        for c in columns:
            vals += gb_exposures[gb_exposures[c]>0][c].tolist()
        
        vals_by_range = generate_weight_bins(vals,n_steps=7,interpolation='fisher-jenks')

        ax_proj = ccrs.epsg(epsg_meters)
        fig, ax_plots = plt.subplots(1,3,
                subplot_kw={'projection': ax_proj},
                figsize=(17,8),
                dpi=500)
        ax_plots = ax_plots.flatten()
        legend_handles = []
        for (i, ((nmin, nmax), width)) in enumerate(vals_by_range.items()):
            legend_handles.append(mpatches.Patch(color=all_colors[i],
                                    label=f"{int(nmin)} - {int(nmax)}"))
            
        legend_handles.append(mpatches.Patch(color=no_exposure_color,
                                        label="No exposure"))
        for idx,(c,asset_label) in enumerate(zip(asset_types,asset_labels)):
            ax = ax_plots[idx]
            total = gb_exposures[c].sum()
            for (i, ((nmin, nmax), width)) in enumerate(vals_by_range.items()):
                ax = plot_polygon_assets(ax,
                        gb_exposures[(gb_exposures[c] >= nmin) & (gb_exposures[c] < nmax+1)],"Exposed",
                        zorder=2 + i + 1,
                        facecolor=all_colors[i],
                        projection=ax_proj)
                ax = plot_polygon_assets(ax,
                            gb_exposures[gb_exposures[c]==0],"No exposure",
                            zorder=2,
                            facecolor=no_exposure_color,
                            projection=ax_proj)
                ax.text(
                        0.05,
                        0.95,
                        f"{figure_texts[idx]} {asset_label}",
                        horizontalalignment='left',
                        transform=ax.transAxes,
                        size=12,
                        weight='bold')

                st = f"Total {asset_label} = £ {round(1e-3*total,2)} million"
                ax.text(
                        0.10,
                        0.85,
                        st,
                        horizontalalignment='left',
                        transform=ax.transAxes,
                        size=12,
                        weight='bold')
        title = f"Flooded Expected Annual Losses \n (£ '000)" 
        fontsize = 12
        bbox_to_anchor=(1.3, 0.85)
        plt.legend(
                handles=legend_handles,
                fontsize=fontsize,
                title=title,
                title_fontproperties={'weight':'bold'},
                bbox_to_anchor=bbox_to_anchor,loc="center",frameon=False)
        plt.tight_layout()
        save_fig(os.path.join(figures,f"EAD_EAEL_risks.png"))



        




    

if __name__ == '__main__':
    CONFIG = load_config()
    main(CONFIG)
