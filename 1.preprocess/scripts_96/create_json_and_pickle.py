# Description: Create .json .pickle for higashi
# Author: genger
# Data: 2022-04-30-Saturday

import pickle
import pandas as pd

import os
import sys
import getopt

## change the path 
path = "./"
os.chdir(path)

## usage
def usage():
    print("")
    print("Create .json .pickle for higashi.\nMost sets are default, if you want to change them please contact genger@wechat:13958598285")
    print("uage: python %s -option <argument>" %sys.argv[0])
    print(" -h/--help ")
    print(" --projectName=<STRING> your project name.")
    print(" --chromNum=<STRING> set chrom number(default without chrXYM).")
    print(" --higashiInput=<STRING> data.txt for higashi input.")
    print(" --tempDir=<STRING> higashi temp dir.")
    print(" --chromsizeDir=<STRING> chrom size file.")
    print(" --cytoDir=<STRING> cytoband file.")
    print(" --pickleDir=<STRING> pickle file.")
    print(" --jsonDir=<STRING> json file.")
    print(" --input_stat=<STRING> input stat file.")

## deal with options
try: 
    opts, args = getopt.getopt(sys.argv[1:], "h", ["help","projectName=","chromNum=","higashiInput=","tempDir=","chromsizeDir=","cytoDir=","pickleDir=","jsonDir=","input_stat="])
except getopt.GetoptError:
    print("ERROR: Get option error.\nYou can contact the author through wechat 13958598285.")
    usage()
    sys.exit(2)

for opt, val in opts:
    if opt in ("-h", "--help"):
        usage()
        sys.exit(1)
    else:
        if opt in ("--projectName"):
            config_name = val
        if opt in ("--chromNum"):
            chrom_num = int(val)
        if opt in ("--higashiInput"):
            data_dir = val
        if opt in ("--tempDir"):
            temp_dir = val
        if opt in ("--chromsizeDir"):
            genome_reference_path = val
        if opt in ("--cytoDir"):
            cytoband_path = val
        if opt in ("--pickleDir"):
            pickle_path = val
        if opt in ("--jsonDir"):
            json_path = val
        if opt in ("--input_stat"):
            input_stat_url = val

try:
    print("Your projectName is:",config_name)
    print("Your chromNum is:",chrom_num)
    print("Your higashiInput is:",data_dir)
    print("Your tempDir is:",temp_dir)
    print("Your chromsizeDir is:",genome_reference_path)
    print("Your cytoDir is:",cytoband_path)
    print("Your pickleDir is:",pickle_path)
    print("Your jsonDir is:",json_path)
    print("Your input_stat is:",input_stat_url)
except:
    print("ERROR: Missing option.")
    usage()
    sys.exit(2)

input_format="higashi_v1"
structured="false"
chrom_list=["\"chr"+str(i)+"\"" for i in range(1,chrom_num+1)];chrom_list="["+",".join(chrom_list)+"]"
resolution=1000000
resolution_cell=1000000
local_transfer_range=1
dimensions=64
loss_mode="zinb"
rank_thres=1
embedding_epoch=80
no_nbr_epoch=80
with_nbr_epoch=60
embedding_name=config_name+"_"+loss_mode
impute_list=["\"chr"+str(i)+"\"" for i in range(1,chrom_num+1)];impute_list="["+",".join(impute_list)+"]"
minimum_distance=1000000
maximum_distance=-1
neighbor_num=5
cpu_num_torch=-1
cpu_num=1
gpu_num=0
UMAP_params=20

with open(json_path,"w") as w:
    w.write("{\n")
    w.write("\"config_name\": \"%s\",\n"%(config_name))
    w.write("\"data_dir\": \"%s\",\n"%(data_dir))
    w.write("\"input_format\": \"%s\",\n"%(input_format))
    w.write("\"structured\": \"%s\",\n"%(structured))
    w.write("\"temp_dir\": \"%s\",\n"%(temp_dir))
    w.write("\"genome_reference_path\": \"%s\",\n"%(genome_reference_path))
    w.write("\"cytoband_path\": \"%s\",\n"%(cytoband_path))
    w.write("\"chrom_list\": %s,\n"%(chrom_list))
    w.write("\"resolution\": %s,\n"%(resolution))
    w.write("\"resolution_cell\": %s,\n"%(resolution_cell))
    w.write("\"local_transfer_range\": %s,\n"%(local_transfer_range))
    w.write("\"dimensions\": %s,\n"%(dimensions))
    w.write("\"loss_mode\": \"%s\",\n"%(loss_mode))
    w.write("\"rank_thres\": %s,\n"%(rank_thres))
    w.write("\"embedding_epoch\": %s,\n"%(embedding_epoch))
    w.write("\"no_nbr_epoch\": %s,\n"%(no_nbr_epoch))
    w.write("\"with_nbr_epoch\": %s,\n"%(with_nbr_epoch))
    w.write("\"embedding_name\": \"%s\",\n"%(embedding_name))
    w.write("\"impute_list\": %s,\n"%(impute_list))
    w.write("\"minimum_distance\": %s,\n"%(minimum_distance))
    w.write("\"maximum_distance\": %s,\n"%(maximum_distance))
    w.write("\"neighbor_num\": %s,\n"%(neighbor_num))
    w.write("\"cpu_num_torch\": %s,\n"%(cpu_num_torch))
    w.write("\"cpu_num\": %s,\n"%(cpu_num))
    w.write("\"gpu_num\": %s,\n"%(gpu_num))
    w.write("\"UMAP_params\": {\"n_neighbors\": %s}\n"%(UMAP_params))
    w.write("}\n")

raw_stat=pd.read_csv(input_stat_url)
if "pixel" not in raw_stat.columns.values:
    output_label_file = open(pickle_path, "wb")
    label_info = {
      "cell_id": [i for i in range(raw_stat.shape[0])],
      "cell_name": raw_stat["cell_id"].tolist(),
      "project": raw_stat["sample"].tolist(),
      "cis_10kb_contact": raw_stat["cis_more_10kb"].tolist()
    }
    pickle.dump(label_info, output_label_file)
    output_label_file.close()
else:
    output_label_file = open(pickle_path, "wb")
    label_info = {
      "cell_id": [i for i in range(raw_stat.shape[0])],
      "cell_name": raw_stat["pixel"].tolist(),
      "nB": raw_stat["iB"].tolist(),
      "nA": raw_stat["iA"].tolist(),
      "bB": raw_stat["bc_B"].tolist(),
      "bA": raw_stat["bc_A"].tolist(),
      "cis_10kb_contact": raw_stat["cis_more_10kb"].tolist()
    }
    pickle.dump(label_info, output_label_file)
    output_label_file.close()



