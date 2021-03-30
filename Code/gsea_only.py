#!/usr/bin/env python3

# SET UP GSEA.JAR FOLDER (LINE 31)

# 1 - conda activate gseapy

# 2 - python /Users/miguel/Box/My_scripts/gseaonly.py --folder /Users/miguel/Box/Marc/Mongolia/RNAseq/Analysis/Deseq2_groups/


#import gseapy as gp
#from gseapy.plot import gseaplot
#from gseapy.plot import barplot, dotplot
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

#parser = argparse.ArgumentParser(description='Get rnk file for GSEA from DESeq2 output')
#parser.add_argument('--file', metavar='file', type=str, help='file from DESeq2')
#args = parser.parse_args()

gsea_folder = '/PATHTO/gsea-3.0.jar'

parser = argparse.ArgumentParser(description='Generate gsea')
parser.add_argument('--folder', metavar='folder', type=str, help='folder to output stuff')
args = parser.parse_args()


#main_folder = '/Users/miguel/Box/Marc/Mongolia/RNAseq/Analysis/Deseq2_groups/'

main_folder = str(args.folder)

files_gsea = glob.glob(f'{main_folder}/*/*/*_afternorm_gsea_ranked.rnk')


name_list = []

for file in files_gsea:
    name = re.findall(r'(?<=gsea/res_total_deseq2_).*(?=_afternorm_gsea_ranked.rnk)', file)[0]
    name_list.append(name)

gmts = {'Hallmarks': 'gseaftp.broadinstitute.org://pub/gsea/gene_sets_final/h.all.v6.2.symbols.gmt',
                'KEGG': 'gseaftp.broadinstitute.org://pub/gsea/gene_sets_final/c2.cp.kegg.v6.2.symbols.gmt',
                'GO_Biological_Process': 'gseaftp.broadinstitute.org://pub/gsea/gene_sets_final/c5.bp.v6.2.symbols.gmt',
                'GO_Molecular_Function': 'gseaftp.broadinstitute.org://pub/gsea/gene_sets_final/c5.mf.v6.2.symbols.gmt',
                'Curated_Gene_Sets': 'gseaftp.broadinstitute.org://pub/gsea/gene_sets_final/c2.cgp.v6.2.symbols.gmt'
                }

for gsea_file in name_list:
    try:
        for name, gset in gmts.items():
        	total_name = f'{main_folder}deseq2_{gsea_file}/gsea/res_total_deseq2_{gsea_file}_afternorm_gsea_ranked.rnk'
        	out_folder = f'{main_folder}deseq2_{gsea_file}/gsea/'
        	set_dir = f'{out_folder}/{name}'
        	os.makedirs(set_dir, exist_ok=True)
        	os.system(f'java -Xmx5000m -cp {gsea_folder} xtools.gsea.GseaPreranked -gmx {gset} -norm meandiv -nperm 1000 -rnk "{total_name}" -scoring_scheme weighted -rpt_label {gsea_file}_wald -create_svgs false -make_sets true -plot_top_x 100 -rnd_seed timestamp -set_max 1000 -set_min 10 -zip_report false -out {set_dir} -gui false')       
    except:
        pass
