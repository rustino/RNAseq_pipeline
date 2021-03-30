#!/usr/bin/env python3

import glob
import os
import re
#import seaborn as sns
#import matplotlib.pyplot as plt
#import matplotlib.patches as mpatches
#from matplotlib_venn import venn2, venn2_circles
import pandas as pd
#import numpy as np
import pathlib
import ntpath
import argparse
from os import path

# 1 - conda activate RNAseq_mac_2


# 2 - python /Users/miguel/Box/My_scripts/Deseq2_plus_rnk.py --folder /Users/miguel/Box/Marc/Mongolia/RNAseq/Analysis/Deseq2_groups_2/


parser = argparse.ArgumentParser(description='Generate Decoder for DESeq2 and analysis')
parser.add_argument('--folder', metavar='folder', type=str, help='folder to output stuff')
parser.add_argument('--Qfold', metavar='Qorts_folder', type=str, help='Qorts_folder')
parser.add_argument('--fext', metavar='file extension', type=str, help='file extension of Qorts file')
parser.add_argument('--folder_script', metavar='file extension', type=str, help='file extension of Qorts file')

#parser.add_argument('--model', metavar='model for DEseq', type=str, help='Model for DEseq2')

args = parser.parse_args()
folder_script = str(args.folder_script)
folder = str(args.folder)
QoRTs = str(args.Qfold)
QoRTs_ext = str(args.fext)
#Model = str(args.model)
#folder = '/Users/miguel/Box/Marc/Mongolia/RNAseq/Analysis/Deseq2_groups/'

files_deseq = glob.glob(f'{folder}*/decoder.deseq2*')

name_list = []

for file in files_deseq:
    name = re.findall(fr'(?<={folder}).*(?=/)', file)[0]
    name_list.append(name)

for name in name_list:
    if path.exists(f"{folder}{name}/decoder.{name}.txt"):
        try:
            os.system(f'Rscript {folder_script}Deseq2.R {name} {folder} {QoRTs} {QoRTs_ext}')
            os.system(f'python {folder_script}rnk_file_generation.py --file {folder}{name}/res_total_{name}_afternorm.txt')
        except:
            pass
    else:
        print(f'{folder}{name}/decoder.{name}.txt decoder not present')