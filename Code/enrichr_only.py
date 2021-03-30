#!/usr/bin/env python3


# 1 - conda activate gseapy

# 2 - python /Users/miguel/Box/My_scripts/enrichr_only.py --folder /Users/miguel/Box/Marc/Mongolia/RNAseq/Analysis/Deseq2_groups_2/


import gseapy as gp
from gseapy.plot import gseaplot
from gseapy.plot import barplot, dotplot
import pandas as pd
import glob
import os
import re
#import seaborn as sns
#import matplotlib.pyplot as plt
#import matplotlib.patches as mpatches
#from matplotlib_venn import venn2, venn2_circles
import numpy as np
import pathlib
import ntpath
import argparse
import pathlib

parser = argparse.ArgumentParser(description='Generate enrichr')
parser.add_argument('--folder', metavar='folder', type=str, help='folder to output stuff')
args = parser.parse_args()

#main_folder = '/Users/miguel/Box/Marc/Mongolia/RNAseq/Analysis/Deseq2_groups/'

main_folder = str(args.folder)

files_gsea = glob.glob(f'{main_folder}/*/decoder.deseq2*')

name_list = []

for file in files_gsea:
    name = re.findall(fr'(?<={main_folder}).*(?=/)', file)[0]
    name_list.append(name) 

#name_list = name_list[:1]

## enrichr

gene_sets_enrichr = ['KEGG_2019_Human', 'GO_Biological_Process_2018', 'GO_Molecular_Function_2018', 'OMIM_Disease', 'ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X', 'BioCarta_2016','VirusMINT','WikiPathways_2019_Human','Virus_Perturbations_from_GEO_down', 'Virus_Perturbations_from_GEO_up']

enrichfolder = ['padj_01_fold_1', 'padj_001_fold_1','padj_005_fold_0','padj_005_fold_1']

for enrichr_file in name_list:
    try:
        for outfolder in enrichfolder:
            pathlib.Path(f'{main_folder}{enrichr_file}/enrichr/{outfolder}/').mkdir(exist_ok=True)
            gene_list = pd.read_csv(f'{main_folder}{enrichr_file}/enrichr/res_total_{enrichr_file}_afternorm_{outfolder}.txt', header=None, sep="\t")
            glist = gene_list.squeeze().str.strip().tolist()
            enr = gp.enrichr(gene_list=glist, cutoff=0.1 , verbose=True, description=f'{outfolder}',gene_sets= gene_sets_enrichr,outdir=f'{main_folder}{enrichr_file}/enrichr/{outfolder}/',format='png')
            enr.results.to_csv(f'{main_folder}{enrichr_file}/enrichr/{outfolder}/{outfolder}.tsv', sep="\t", header=True, index=False, mode='w')
    except:
        pass