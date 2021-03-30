#!/usr/bin/env python

# Insert a gene_conversion table (line 25)
# Insert a master table with the description of the samples (line 26)

import glob
import os
import re
import pathlib
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import matplotlib
from scipy.stats import ttest_ind
import scipy.stats as stats
import scipy
import itertools
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import pathlib
import ntpath
import argparse

conversion_path = '/pathto/ENG_gene_conversion_nodup.txt'
master_table_path = '/pathto/Master_Table_RNAseq_mutation.csv'


parser = argparse.ArgumentParser(description='get higest SD and dedup')
parser.add_argument('--folder', metavar='folder', type=str, help='folder to output stuff')
#parser.add_argument('--Qfold', metavar='Qorts_folder', type=str, help='Qorts_folder')
#parser.add_argument('--fext', metavar='file extension', type=str, help='file extension of Qorts file')

args = parser.parse_args()

folder_output = str(args.folder)

Gene_conversion = pd.read_csv(f'{conversion_path}', sep='\t')

Master_table_RNAseq = pd.read_csv(f"{master_table_path}",       
                              sep=',')
Master_table_RNAseq = Master_table_RNAseq.where(Master_table_RNAseq['Delete'] == 0)


VST_files = glob.glob(f'{folder_output}/*/*_afternorm_vst.txt')

name_list = []

for file in VST_files:
    name = re.findall(fr'(?<={folder_output}).*.(?=/)', file)[0]
    name_list.append(name)

numbers_sd = [150,300,500,1000,2500,5000,10000,20000]
for name in name_list:
    try:
        for number_sd in numbers_sd:
            folder = f'{folder_output}{name}/table_reduced_sd/'
            pathlib.Path(f'{folder}').mkdir(exist_ok=True)
            file_VST = f'{folder_output}{name}/{name}_afternorm_vst.txt'
            Table_RNA_VST = pd.read_csv(file_VST, sep='\t')
            SD_RPKM = Table_RNA_VST.std(axis=1, skipna=None, level=None, ddof=1, numeric_only=None)
            SD_RPKM_DF = SD_RPKM.to_frame()
            SD_RPKM_DF.columns = ["SD"]
            #number_sd = -1000
            Table_RNA_VST_SD =  Table_RNA_VST.join(SD_RPKM_DF, 
                        how='outer').sort_values(by=['SD'], ascending = False).drop(axis=1,
                        labels="SD").iloc[:number_sd,:].transpose()
            Table_RNA_VST_SD.columns = Table_RNA_VST_SD.iloc[0]
            Table_RNA_VST_SD = Table_RNA_VST_SD.iloc[1:,:].reset_index()
            Table_RNA_VST_SD = Table_RNA_VST_SD.rename(columns={'index':'Sample_code_RNAseq'})
            Table_RNA_VST_SD_2 = Table_RNA_VST_SD.merge(Master_table_RNAseq, on='Sample_code_RNAseq', how = 'left')
            Table_RNA_VST_SD = Table_RNA_VST_SD.T
            Table_RNA_VST_SD.columns = Table_RNA_VST_SD.iloc[0]
            del Table_RNA_VST_SD.index.name
            Table_RNA_VST_SD = Table_RNA_VST_SD.iloc[1:,:]
            Table_RNA_VST_SD['ENSG']  = Table_RNA_VST_SD.index.str.split('.', 1).str[0]
            Table_RNA_VST_SD = Table_RNA_VST_SD.reset_index().drop('index', axis=1)
            Table_RNA_VST_SD = Table_RNA_VST_SD.merge(Gene_conversion, left_on='ENSG', 
                      right_on='#hg19.ensGene.name2').drop(columns=['#hg19.ensGene.name2']
                                        ).rename(columns={'hg19.ensemblToGeneName.value':'gene symbol',
                                                         'ENSG':'gene description'})
            cols_to_move = ['gene symbol','gene description']
            Table_RNA_VST_SD = Table_RNA_VST_SD[ cols_to_move + [ col for col in Table_RNA_VST_SD.columns if col not in cols_to_move ] ]
            #Table_RNA_VST_SD = Table_RNA_VST_SD[ ['gene description'] + [ col for col in Table_RNA_VST_SD.columns if col != 'gene description' ] ]
            #Table_RNA_VST_SD = Table_RNA_VST_SD[ ['gene symbol'] + [ col for col in Table_RNA_VST_SD.columns if col != 'gene symbol' ] ]
            number_samples = Table_RNA_VST_SD.shape[1] - 2
            number_genes = Table_RNA_VST_SD.shape[0]
            path_of_table = f'{folder}{name}_vst_{number_sd}_sd'
            Table_RNA_VST_SD.to_csv(path_of_table + str('.gct'), 
                                          sep="\t", header=True, index=False, index_label = False, mode='w')
            with open(path_of_table + str('.gct'),"r") as original:
                data = original.read()
            with open(path_of_table + str('.gct'),"w") as modified:
                modified.write(f'#1.2\n{number_genes}\t{number_samples}\n' + data)
            cols = Table_RNA_VST_SD.columns.drop(['gene symbol','gene description'])
            Table_RNA_VST_SD[cols] = Table_RNA_VST_SD[cols].apply(pd.to_numeric, errors='coerce')
            Table_RNA_VST_SD = Table_RNA_VST_SD.groupby(['gene symbol'], as_index=False).median()
            Table_RNA_VST_SD['gene description'] = Table_RNA_VST_SD['gene symbol']
            cols_to_move = ['gene symbol','gene description']
            Table_RNA_VST_SD = Table_RNA_VST_SD[ cols_to_move + [ col for col in Table_RNA_VST_SD.columns if col not in cols_to_move ] ]
            number_genes = Table_RNA_VST_SD.shape[0]
            Table_RNA_VST_SD.to_csv(path_of_table + str ('_nodup') + str('.gct'), 
                                    sep="\t", header=True, index=False, index_label = False, mode='w')
            with open(path_of_table + str ('_nodup') + str('.gct'),"r") as original:
                data = original.read()
            with open(path_of_table + str ('_nodup') + str('.gct'),"w") as modified:
                modified.write(f'#1.2\n{number_genes}\t{number_samples}\n' + data)
    except:
        pass