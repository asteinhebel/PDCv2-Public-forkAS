#----------------------------------------------------------------------------------
#-- Company: GRAMS
#-- Designer: Gabriel Lessard
#--
#-- Create Date: 2025-03-11
#-- Description:
#--
#-- Dependencies:
#-- Revision:
#-- Revision 1.0 - File Created
#-- Additional Comments:
#----------------------------------------------------------------------------------


import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.ticker as mticker
from matplotlib.ticker import FormatStrFormatter
import seaborn as sns
from matplotlib.colors import *
import math
import matplotlib.lines as mlines
import argparse
from pdcv2_modules.pixMap import xy_map
from matplotlib import colors, cm



def plot_results_df(df: pd.DataFrame, criteria="SPAD_TCR (cps)", head_id=None):
    fig, axes = plt.subplots(2, 4, figsize=(14, 8), width_ratios=[1, 1, 1.5, 1.5])

    palette="deep"

    #sns.set(rc={'axes.yscale': 'log'})

    # Plot scatter per PDC
    ax = sns.scatterplot(data=df,
                         x="SPAD_idx",
                         y=criteria,
                         hue="iPDC",
                         palette=palette,
                         ax=axes[0,0],
                         s=10,
                         legend=False)
    ax.set_yscale("log")



    # Plot sorted population
    ax_pop = sns.ecdfplot(data=df, y=criteria, hue="iPDC", palette=palette, ax=axes[0,1], legend=False)
    #labels = []
    for iPDC, df_pdc in df.groupby("iPDC"):
        ax_pop.axhline(y=df_pdc[criteria].median(), color=sns.color_palette(palette)[iPDC], linestyle="--")
        ax_pop.axhline(y=df_pdc[criteria].mean(), color=sns.color_palette(palette)[iPDC], linestyle="-."  )
        #labels.append(f'PDC{iPDC}, med.: {int(df_pdc[criteria].median())}, mean: {int(df_pdc[criteria].mean())}')
    ax_pop.set_yscale("log")
    #ax_pop.legend(labels=labels)

    # Plot Histogram
    gs = axes[1, 0].get_gridspec()
    for axi in axes[1, :2]:
        axi.remove()
    ax_merged = fig.add_subplot(gs[1, :2])
    df_histo = df[df[criteria].between(10, df[criteria].mean())]
    ax_merged = sns.histplot(data=df_histo, x=criteria, palette=palette, hue="iPDC",
                             ax=ax_merged, kde=True, bins=40, multiple="dodge", log_scale=False)
    ax_merged.set_xscale("log")
    ax_merged.set_yscale("log")
    ax_merged.set_ylim(ymin=1)


    # Plot heatmaps
    heatmap_loc = [(0,2), (1,2) ,(1,3), (0,3)]

    df[['X', 'Y']] = df['SPAD_idx'].apply(xy_map).apply(pd.Series)
    #df["TCR_log"] =  np.log10(df[criteria])

    for iPDC, df_pdc in df.groupby("iPDC"):
        ax = axes[heatmap_loc[iPDC]]
        ax.set_title(f"PDC{iPDC}")
        #df_pdc = df_pdc[df_pdc[criteria] < 500_000]
        df_hm = df_pdc.pivot_table(index='Y', columns='X', values=criteria)
        
        #ax = sns.heatmap(data=df_hm, square=True, clip_on=False, xticklabels=10,yticklabels=10,
        #                 robust=False, cmap="viridis",  ax=ax, cbar_kws={"fraction":0.045}, vmax=math.log10(100_000), norm=LogNorm())
        im =ax.imshow(df_hm, norm=colors.LogNorm(vmin=df_pdc[criteria].min()+1, vmax=200_000)) # vmax=df_pdc[criteria].max()
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        ax.invert_yaxis()

    items = []
    for iPDC, df_pdc in df.groupby("iPDC"):
        label= f'PDC{iPDC}, med.: {(df_pdc[criteria].median()):.3f}, mean: {(df_pdc[criteria].mean()):.3f}'
        items.append(mlines.Line2D([], [], color=sns.color_palette(palette)[iPDC], marker='o', ls='', label=label))
    # etc etc


    axes[0,0].legend(handles=items, bbox_to_anchor=(0., 1.02, 3.2, .102), loc=3,
               ncol=2,  borderaxespad=0, )
    #fig.legend(handles=items, ncols=4, loc="lower right")
    plt.subplots_adjust(left=0.05, bottom=0.08, right=0.96, top=0.92, wspace=0.3, hspace=0.2)
    #fig.tight_layout()
    if head_id is not None:
        fig.text(x=0.70, y=0.95, s=head_id, fontsize=20)

    return fig



def main():
    # Create an ArgumentParser object
    parser = argparse.ArgumentParser(description="Process a file.")

    # Add an argument for the file, specifying its type as argparse.FileType('r')
    parser.add_argument('--file', type=argparse.FileType('r'), help='The input file to process')

    # Parse the command-line arguments
    args = parser.parse_args()

    # Access the file object
    if args.file:
        # Get the file name using the .name attribute
        file_name = args.file.name
        print(file_name)
        head_id = file_name.split('_')[3]
        if head_id[0] != "H":
            head_id = None
        df = csv_to_df(args.file.name)
        fig = plot_results_df(df, criteria="SPAD_TCR (cps)", head_id=head_id)
    plt.show()
    return



def csv_to_df(filepath:str):
    data = pd.read_csv(filepath, sep=';', header=0)
    num_pdc = math.floor(len(data.keys()) / 4)


    pdcs = [int(k.split("SPAD_idx")[1]) for k in data.keys() if "SPAD_idx" in k]
    data = data[data[f"SPAD_idx{min(pdcs)}"] < 4096]
    print(data)

    df = pd.DataFrame()

    for iPDC in pdcs:
        df_pdc = pd.DataFrame()
        df_pdc["SPAD_idx"] = data[f"SPAD_idx{iPDC}"]
        df_pdc["SPAD_TCR (cps)"] = data[f"SPAD_TCR{iPDC}"]
        df_pdc["iPDC"] = iPDC
        df = pd.concat([df, df_pdc], ignore_index=True)
        print(f"PDC{iPDC} has {(df_pdc['SPAD_TCR (cps)'] == 0).sum()} SPADs with no TCR value")


    spad_active_area= 42**2-(28**2-math.pi*14**2) #* 1e-6clear
    spad_active_area = 78*78*0.26 * 1e-6
    df["SPAD_TCR (cps/mm2 active)"] = df["SPAD_TCR (cps)"]/spad_active_area
    # SPAD total area in um2 (78um pitch)
    spad_total_area = 78*78 * 1e-6
    df["SPAD_TCR (cps/mm2 total)"] = df["SPAD_TCR (cps)"]/spad_total_area

    return df

if __name__ == "__main__":
    main()
