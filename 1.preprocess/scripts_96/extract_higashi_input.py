# Description: Transform barcode name to id(0,1,2,3,...), you should provide a barcode.csv file
# Author: genger
# Data: 2022-07-15-Friday

import pandas as pd
import seaborn as sns
from tqdm import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

import os
import sys
import getopt

## change the path 
path = "./"
os.chdir(path)

## usage
def usage():
    print("")
    print("Transform barcode name to id(0,1,2,3,...).")
    print("uage: python %s -option <argument>" %sys.argv[0])
    print(" -h/--help ")
    print(" --input_stat=<STRING> input stat file.")
    print(" --all=<STRING> all data file.")
    print(" --out=<STRING> out file for higashi.")
    print(" --fig=<STRING> output figure file.")

## deal with options
try: 
    opts, args = getopt.getopt(sys.argv[1:], "h", ["help","input_stat=","all=","out=","fig="])
except getopt.GetoptError:
    print("ERROR: Get option error.\nYou can contact the author through wechat 13958598285.")
    usage()
    sys.exit(2)

for opt, val in opts:
    if opt in ("-h", "--help"):
        usage()
        sys.exit(1)
    else:
        if opt in ("--input_stat"):
            input_stat_url = val
        if opt in ("--all"):
            rawdata = val
        if opt in ("--out"):
            out_url = val
        if opt in ("--fig"):
            fig_url = val

try:
    print("Your input_stat is:",input_stat_url)
    print("Your all is:",rawdata)
    print("Your out is:",out_url)
    print("Your fig is:",fig_url)
except:
    print("ERROR: Missing option.")
    usage()
    sys.exit(2)

raw_stat=pd.read_csv(input_stat_url)

all_df=pd.read_csv(rawdata,sep="\t",dtype={"barcode": str, "chrom1": str, "pos1": int, "chrom2": str, "pos2": int, "count": int})
out_df=all_df.loc[all_df["barcode"].isin(raw_stat["pixel"].values)]
out_df=out_df.loc[(out_df["chrom1"]==out_df["chrom2"]) & (abs(out_df["pos1"]-out_df["pos2"])>=10000)]

bc_dic={v:i for i,v in enumerate(raw_stat["pixel"].values)}
out_df["barcode"]=out_df["barcode"].apply(lambda x: bc_dic[x])
out_df.columns=["cell_id","chrom1","pos1","chrom2","pos2","count"]
out_df.to_csv(out_url,sep="\t",index=None)

plt.figure(figsize=(8,8))
normlizer=mpl.colors.Normalize(vmin=np.min(raw_stat["cis_more_10kb"]), vmax=np.quantile(raw_stat["cis_more_10kb"],0.99))
sns.scatterplot(data=raw_stat,x="iB",y="iA",marker="s",hue="cis_more_10kb",hue_norm=normlizer,palette="rainbow",alpha=0.9)
# plt.gca().invert_yaxis()
# plt.gca().invert_xaxis()
plt.title("higashi_more10kb_cis_contact heatmap")
#plt.axis("off")
plt.savefig(fig_url,dpi=500,bbox_inches="tight")



