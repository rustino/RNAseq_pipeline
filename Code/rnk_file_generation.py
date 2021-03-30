#!/usr/bin/env python
# coding: utf-8

# Insert Gene Conversion Table at line 44

# example use: 

# 1 - conda activate RNAseq_mac_2

# 2 - python /Users/miguel/Box/My_scripts/rnk_file_generation.py --file /Users/miguel/Box/Marc/Mongolia/RNAseq/Analysis/Deseq2_groups/deseq2_Library_comp_Non_stranded_vs_Stranded_TUMOR/res_total_deseq2_Library_comp_Non_stranded_vs_Stranded_TUMOR_afternorm.txt

import glob
import os
import re
from os import path
import os.path
import pandas as pd
import ntpath
import argparse
import pathlib


# In[ ]:


parser = argparse.ArgumentParser(description='Get rnk file for GSEA from DESeq2 output')
parser.add_argument('--file', metavar='file', type=str, help='file from DESeq2')
args = parser.parse_args()


# In[ ]:


input_table = str(args.file)
path = os.path.dirname(input_table) + str("/")
file = re.match(r'(?<=).*(?=.txt)', os.path.basename(input_table))[0]

pathlib.Path(f'{path}/enrichr/').mkdir(exist_ok=True)
pathlib.Path(f'{path}/gsea/').mkdir(exist_ok=True)

# In[ ]:


Gene_conversion = pd.read_csv('/Users/MiguelT/Box/Genomes/ENG_gene_conversion_nodup.txt',       
    sep='\t')


# In[ ]:


res_table_RNAseq = pd.read_csv(input_table,  sep='\t').rename(columns={'Unnamed: 0':'ENSG'})
res_table_RNAseq['ENSG']  = res_table_RNAseq['ENSG'].str.split('.', 1).str[0]
res_table_RNAseq = res_table_RNAseq.merge(Gene_conversion, left_on='ENSG', 
          right_on='#hg19.ensGene.name2').drop(columns=['#hg19.ensGene.name2']
                            ).rename(columns={'hg19.ensemblToGeneName.value':'Gene'})
cols_2 = ['ENSG', 'Gene', 'baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'padj']
res_table_RNAseq = res_table_RNAseq[cols_2]

try:
	res_table_RNAseq.where((res_table_RNAseq['padj'] <= 0.05) & 
                       (res_table_RNAseq['log2FoldChange'] >= abs(1))
                      ).dropna(subset=['padj','Gene']).loc[:,'Gene'].to_csv(f'{path}/enrichr/{file}_padj_005_fold_1.txt', 
                      sep=",", header=False, index=False, mode='w')
except:
	pass

try:
	res_table_RNAseq.where((res_table_RNAseq['padj'] <= 0.05)
                      ).dropna(subset=['padj','Gene']).loc[:,'Gene'].to_csv(f'{path}/enrichr/{file}_padj_005_fold_0.txt', 
                      sep=",", header=False, index=False, mode='w')
except:
	pass
try:
	res_table_RNAseq.where((res_table_RNAseq['padj'] <= 0.1) & 
                       (res_table_RNAseq['log2FoldChange'] >= abs(1))
                      ).dropna(subset=['padj','Gene']).loc[:,'Gene'].to_csv(f'{path}/enrichr/{file}_padj_01_fold_1.txt', 
                      sep=",", header=False, index=False, mode='w')
except:
	pass
try:
	res_table_RNAseq.where((res_table_RNAseq['padj'] <= 0.01) & 
                       (res_table_RNAseq['log2FoldChange'] >= abs(1))
                      ).dropna(subset=['padj','Gene']).loc[:,'Gene'].to_csv(f'{path}/enrichr/{file}_padj_001_fold_1.txt', 
                      sep=",", header=False, index=False, mode='w')
except:
	pass
try:
	res_table_RNAseq.where((res_table_RNAseq['padj'] <= 0.001) & 
                       (res_table_RNAseq['log2FoldChange'] >= abs(1))
                      ).dropna(subset=['padj','Gene']).loc[:,'Gene'].to_csv(f'{path}/enrichr/{file}_padj_001_fold_1.txt', 
                      sep=",", header=False, index=False, mode='w')
except:
	pass

res_table_RNAseq_2 = res_table_RNAseq.dropna(subset=['padj']).drop(columns=['log2FoldChange'
    ,'pvalue','padj','lfcSE','ENSG'])
res_table_RNAseq_2['quartile'] = pd.qcut(res_table_RNAseq_2['baseMean'], 4, ['q1','q2','q3','q4'])
res_table_RNAseq_2['quartile'] = res_table_RNAseq_2.quartile.astype(str)
res_table_RNAseq_2 = res_table_RNAseq_2.where(res_table_RNAseq_2['quartile'] != 'q4'
                                         ).drop(columns=['quartile','baseMean']
                                         ).dropna().sort_values(by='stat', ascending = False)


# In[ ]:



path_file = f'{path}{file}_gene_names.txt'


# In[ ]:


file_2 = re.match(r'(?<=).*(?=.txt)', os.path.basename(input_table))[0]
path_file_2 = f'{path}/gsea/{file_2}_gsea_ranked.rnk'


# In[ ]:


res_table_RNAseq.to_csv(path_file, sep="\t", header=True, index=False, mode="w")  
res_table_RNAseq_2.to_csv(path_file_2, sep="\t", header=False, index=False, mode="w")  

