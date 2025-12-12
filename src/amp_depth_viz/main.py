import argparse
import pandas as pd
from bokeh.plotting import figure, show
from bokeh.io import export_png
from bokeh.plotting import figure, output_file, save
from bokeh.layouts import gridplot
import numpy as np
from bokeh.models import ColumnDataSource
import os
import sys

class _colors:
    """Some colours that someone thought were nice."""
    cerulean = "#0084A9"
    not_black = "#001A21"
    feldgrau = "#455556"
    dim_gray = "#666666"
    light_cornflower_blue = "#90C5E7"
    dark_gray = "#B5AEA7"
    isabelline = "#F0EFED"
    medium_spring_bud = "#B8E986"
    cinnabar = "#EF4134"
    sandstorm = "#F5CC49"
    fandango = "#A53F96"
    green = "#17BB75"
    verdigris = "#54B8B1"

def main():

    parser=argparse.ArgumentParser(description="create wf-artic like amplicon coverage plots using bokeh")
    parser.add_argument("coveragebed", help="a full sample bed files contains all genome depth information from mosdepth")
    parser.add_argument("--output", default="amplicon_coverage.html",help="html output file")
    parser.add_argument("--xlim", default=30000, type=int, help="max of the genome position")
    parser.add_argument("--ylim", default=800, type=int, help="the expected max depth")
    parser.add_argument("--threshold", default=20, type=int, help="Max number of job to run together, default as 20")
    parser.add_argument("--ncols", default=3, type=int, help="number of columns to grid the plots in the html pages, default as 3")
    args = parser.parse_args()
    

    ## check the inputs
    bed_path = args.coveragebed

    if not os.path.exists(bed_path):
        print(f"{bed_path} does not exist, please double check your input")
        sys.exit()
    #read the bed file with pandas
    df = pd.DataFrame()
    cols = ["chrome","start","end", "depth","pool","sample"]
    try:
        df = pd.read_csv(bed_path, sep="\t", header=None,names=cols)
        #print(df.iloc[0,0])
        if (df.iloc[0,0] == "chrome"):
            df = pd.read_csv(bed_path, sep="\t")
    except:
        print(f"Error when loading the bed file {bed_path}, please check the input format")
        sys.exit()
    df["pos"] = (df["start"] + df["end"])/2

    Colors = _colors()
    plot_list = []
    for sample in list(df["sample"].unique()):
        df_sub = df[df["sample"] == sample]
        df_sum = df_sub[["start","end","depth", "pos"]].groupby(by=["start","end","pos"]).sum()
        mean_depth = df_sum.depth.mean()
        pass_ratio = 100 * (df_sum.depth >= args.threshold).sum() / len(df_sum.depth)
        title="{}: {:.0f}X, {:.1f}% > {}X".format(
                    'test', mean_depth, pass_ratio, args.threshold)
        source = ColumnDataSource(data=dict(
                depth_pool_1 = df_sub[df_sub.pool == 1]["depth"],
                depth_pool_2 = df_sub[df_sub.pool == 2]["depth"],
                position = df_sub.pos.unique(),
            ))
        TOOLTIPS = [
            ("Position", "@position"),
            ("Pool-1", "@depth_pool_1"),
            ("Pool-2", "@depth_pool_2"),
            ]
        p = figure(width=400, height=250, tooltips=TOOLTIPS, title=title,
                    x_axis_label='position', y_axis_label='depth',
                    x_range=(0,args.xlim), y_range=(0,args.ylim))
        p.line(x='position', y='depth_pool_1', color=Colors.dark_gray, source=source)
        p.varea(x='position', y1=0, y2='depth_pool_1', source=source,
                                fill_color=Colors.dark_gray, alpha=0.7,
                                muted_color=Colors.dark_gray, muted_alpha=0.2)
        p.line(x='position', y='depth_pool_2',source=source, color=Colors.verdigris)
        p.varea(x='position', y1=0, y2='depth_pool_2', source=source,
                                fill_color=Colors.verdigris, alpha=0.7,
                                muted_color=Colors.verdigris, muted_alpha=0.2)
        plot_list.append(p)
    
    grid_plot= gridplot(plot_list, ncols=args.ncols)
    output_file(filename=args.output, title="Amplicon Coverage")
    save(grid_plot)