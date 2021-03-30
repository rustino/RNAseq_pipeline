#!/usr/bin/env python


import glob
import os
import re
from os import path
import os.path
import pandas as pd
import pathlib
import itertools
import ntpath
import argparse

# Insert a master table with the description of the samples (line 24)

# Customize from line 56. No way of doing pipeline for this, too many particular issues.

# 1 - conda activate RNAseq_mac_2


# 2 - python /Users/miguel/Box/My_scripts/01_Decoder.py --folder /Users/miguel/Box/Marc/Mongolia/RNAseq/Analysis/Deseq2_groups_2/ --Qfold /Users/miguel/Box/Marc/Mongolia/RNAseq/QoRTs_filtered/ --fext QC.filtered.geneCounts.formatted.for.DESeq.txt.gz

master_table_path = '/pathto/Master_Table_RNAseq_mutation.csv'

parser = argparse.ArgumentParser(description='Generate Decoder for DESeq2')
parser.add_argument('--folder', metavar='folder', type=str, help='folder to output stuff')
parser.add_argument('--Qfold', metavar='Qorts_folder', type=str, help='Qorts_folder')
parser.add_argument('--fext', metavar='file extension', type=str, help='file extension of Qorts file')

args = parser.parse_args()

folder_output = str(args.folder)
QoRTs = str(args.Qfold)
QoRTs_ext = str(args.fext)
#folder_output = '/Users/miguel/Box/Marc/Mongolia/RNAseq/Analysis/Deseq2_groups_filter/'

pathlib.Path(f'{folder_output}').mkdir(exist_ok=True)


Master_table_RNAseq = pd.read_csv(f'{master_table_path}',       
                              sep=',')

#QoRTs = '/Users/miguel/Box/Marc/Mongolia/RNAseq/QoRTs_filtered/'
#QoRTs_ext = 'QC.filtered.geneCounts.formatted.for.DESeq.txt.gz'
#QoRTs_files = glob.glob(f'{QoRTs}*/QC.geneCounts.formatted.for.DESeq.txt.gz')
QoRTs_files = glob.glob(f'{QoRTs}/*/{QoRTs_ext}')


name_list = []

for file in QoRTs_files:
    name = re.findall(fr'(?<={QoRTs}/).*.(?=/)', file)[0]
    name_list.append(name)

name_list = [x for x in name_list if x not in ['JML248-MGL358-NORMAL', 'JML249-MGL359-TUMOR', 
                                  'JML123-MGL39-TUMOR','JML166-MGL245-TUMOR','JML137-MGL105-TUMOR',
                                  'JML219-MGL154-NORMAL','JML094-NY1663-TUMOR','JML085-NY1645-TUMOR',
                                  'JML219-MGL153-TUMOR']]

def HBV_collapsing(row):
   if row['Etiology_master'] == 'HBV/HCV/HDV' :
      return 'HBV_HDV'
   if row['Etiology_master'] == 'HCV' :
      return 'HCV' 
   if row['Etiology_master'] == 'HBV/HDV' :
      return 'HBV_HDV'  
   if row['Etiology_master'] == 'Non-infected' :
      return 'Non_infected'  
   if row['Etiology_master'] == 'HBV' :
      return 'HBV'  

Master_table_RNAseq['Etiology_reduced'] = Master_table_RNAseq.apply (HBV_collapsing, axis=1)

Master_table_RNAseq = Master_table_RNAseq.where(Master_table_RNAseq['Delete'] == 0
                                               ).dropna(subset=['Delete'])
Master_table_RNAseq = Master_table_RNAseq.where(Master_table_RNAseq['Sample_code_RNAseq'] != 'MGL358'
                                               ).dropna(subset=['Sample_code_RNAseq'])
Master_table_RNAseq = Master_table_RNAseq.where(Master_table_RNAseq['Sample_code_RNAseq'] != 'MGL359'
                                               ).dropna(subset=['Sample_code_RNAseq'])
Master_table_RNAseq = Master_table_RNAseq.where(Master_table_RNAseq['Sample_code_RNAseq'] != 'MGL39'
                                               ).dropna(subset=['Sample_code_RNAseq'])
Master_table_RNAseq = Master_table_RNAseq.where(Master_table_RNAseq['Sample_code_RNAseq'] != 'MGL105'
                                               ).dropna(subset=['Sample_code_RNAseq'])
Master_table_RNAseq = Master_table_RNAseq.where(Master_table_RNAseq['Sample_code_RNAseq'] != 'MGL245'
                                               ).dropna(subset=['Sample_code_RNAseq'])
Master_table_RNAseq = Master_table_RNAseq.where(Master_table_RNAseq['Sample_code_RNAseq'] != 'MGL154'
                                               ).dropna(subset=['Sample_code_RNAseq'])
Master_table_RNAseq = Master_table_RNAseq.where(Master_table_RNAseq['Sample_code_RNAseq'] != 'NY1663'
                                               ).dropna(subset=['Sample_code_RNAseq'])
Master_table_RNAseq = Master_table_RNAseq.where(Master_table_RNAseq['Sample_code_RNAseq'] != 'NY1645'
                                               ).dropna(subset=['Sample_code_RNAseq'])
Master_table_RNAseq = Master_table_RNAseq.where(Master_table_RNAseq['Sample_code_RNAseq'] != 'MGL153'
                                               ).dropna(subset=['Sample_code_RNAseq'])
Master_table_RNAseq = Master_table_RNAseq.replace({'/': '_'}, regex=True)
Master_table_RNAseq = Master_table_RNAseq.replace({'-': '_'}, regex=True)
Master_table_RNAseq['sample.ID'] = Master_table_RNAseq['JML']+'-'+Master_table_RNAseq['Sample_code_RNAseq']+'-'+Master_table_RNAseq['Type_batch']

columns_reduce = ['JML','Sample_code','Sample_code_RNAseq', 'Etiology_reduced', 
                   'Library','Etiology_master','Dataset','Type','Dataset2','sample.ID',
                   'PCDH7']
columns_to_loop = ['PCDH7']
Tumor_control = ['TUMOR']
Table_to_do = Master_table_RNAseq[columns_reduce]
Table_to_do = Table_to_do.where(Table_to_do['Dataset'] == 'Mongolia').dropna()
fisher_exact_results = {}
results_df_fisher_exact_all = pd.DataFrame()
for kind in Tumor_control:
    Table_to_do_2 = Table_to_do.where(Table_to_do['Type'] == kind).dropna()
    for var in columns_to_loop:
        conds = list(Table_to_do[var].unique())
        conds_nona = [x for x in conds if str(x) != 'nan']
        for subset in itertools.combinations(conds_nona, 2):
            condition1 = subset[0]
            condition2 = subset[1]
            output_name = f'deseq2_{var}_comp_{condition1}_vs_{condition2}_Mongolia_{kind}'
            pathlib.Path(f'{folder_output}{output_name}').mkdir(exist_ok=True)
            output_total= f'{folder_output}{output_name}/decoder.{output_name}.txt'
            Table_to_do_2.where((Table_to_do_2 [var] == condition1
                    ) | (Table_to_do_2 [var] == condition2)
                    ).dropna(subset=[var]).to_csv(output_total, 
                        sep="\t", header=True, index=False, mode='w')

columns_reduce = ['JML','Sample_code','Sample_code_RNAseq', 'Etiology_reduced', 
                   'Library','Etiology_master','Dataset','Type','Dataset2','sample.ID',
                   'PCDH7']
columns_to_loop = ['PCDH7']
Tumor_control = ['TUMOR']
Table_to_do = Master_table_RNAseq[columns_reduce]
Table_to_do = Table_to_do.where(Table_to_do['Dataset'] != 'Mongolia').dropna()
fisher_exact_results = {}
results_df_fisher_exact_all = pd.DataFrame()
for kind in Tumor_control:
    Table_to_do_2 = Table_to_do.where(Table_to_do['Type'] == kind).dropna()
    for var in columns_to_loop:
        conds = list(Table_to_do[var].unique())
        conds_nona = [x for x in conds if str(x) != 'nan']
        for subset in itertools.combinations(conds_nona, 2):
            condition1 = subset[0]
            condition2 = subset[1]
            output_name = f'deseq2_{var}_comp_{condition1}_vs_{condition2}_BMS_{kind}'
            pathlib.Path(f'{folder_output}{output_name}').mkdir(exist_ok=True)
            output_total= f'{folder_output}{output_name}/decoder.{output_name}.txt'
            Table_to_do_2.where((Table_to_do_2 [var] == condition1
                    ) | (Table_to_do_2 [var] == condition2)
                    ).dropna(subset=[var]).to_csv(output_total, 
                        sep="\t", header=True, index=False, mode='w')

#columns_reduce = ['JML','Sample_code','Sample_code_RNAseq', 'Etiology_reduced',  'Etiology_HDV',
#                   'Library','Etiology_master','Dataset','Type','Dataset2','sample.ID','Cirrhosis']
#columns_to_loop = ['Etiology_reduced','Etiology_master','Dataset2','Dataset','Cirrhosis']
#Tumor_control = ['TUMOR']
#Table_to_do = Master_table_RNAseq[columns_reduce]
#Table_to_do = Table_to_do.where(Table_to_do['Dataset'] == 'Mongolia').dropna()
#fisher_exact_results = {}
#results_df_fisher_exact_all = pd.DataFrame()
#for kind in Tumor_control:
#    Table_to_do_2 = Table_to_do.where(Table_to_do['Type'] == kind).dropna()
#    for var in columns_to_loop:
#        conds = list(Table_to_do[var].unique())
#        conds_nona = [x for x in conds if str(x) != 'nan']
#        for subset in itertools.combinations(conds_nona, 2):
#            condition1 = subset[0]
#            condition2 = subset[1]
#            output_name = f'deseq2_{var}_comp_{condition1}_vs_{condition2}_Mon_{kind}'
#            pathlib.Path(f'{folder_output}{output_name}').mkdir(exist_ok=True)
#            output_total= f'{folder_output}{output_name}/decoder.{output_name}.txt'
#            Table_to_do_2.where((Table_to_do_2 [var] == condition1
#                    ) | (Table_to_do_2 [var] == condition2)
#                    ).dropna(subset=[var]).to_csv(output_total, 
#                        sep="\t", header=True, index=False, mode='w')