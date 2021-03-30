#!/usr/bin/env python

# Modify particular conditions at line 77 and 85

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

# conda activate RNAseq_mac_2

# python /Users/miguel/Box/My_scripts/PCA.py --folder /Users/miguel/Box/Marc/Mongolia/RNAseq/Analysis/Deseq2_groups_filter/

master_table_path = '/pathto/Master_Table_RNAseq_mutation.csv'


font = {'fontname':'Arial'}
sns.set_style("darkgrid", {"axes.facecolor": ".9"})
sns.set_context("talk")

parser = argparse.ArgumentParser(description='Generate PCA')
parser.add_argument('--folder', metavar='folder', type=str, help='PCA')
args = parser.parse_args()

folder_output = str(args.folder)

#folder_output = '/Users/MiguelT/Box/Marc/Mongolia/RNAseq/Analysis/Deseq2_1_batch_filter_10_01_2020_vst/'

Master_table_RNAseq = pd.read_csv(f"{master_table_path}",       
                              sep=',')

Master_table_RNAseq = Master_table_RNAseq.where(Master_table_RNAseq['Delete'] == 0)

VST_files = glob.glob(f'{folder_output}/*/*_afternorm_vst.txt')

name_list = []

for file in VST_files:
    name = re.findall(fr'(?<={folder_output}).*.(?=/)', file)[0]
    name_list.append(name)

#name_list = name_list[52:]

for name in name_list:
    try:
        folder = f'{folder_output}{name}/PCA/'
        pathlib.Path(f'{folder}').mkdir(exist_ok=True)
        file_VST = f'{folder_output}{name}/{name}_afternorm_vst.txt'
        Table_RNA_VST = pd.read_csv(file_VST, sep='\t')
        SD_RPKM = Table_RNA_VST.std(axis=1, skipna=None, level=None, ddof=1, numeric_only=None)
        SD_RPKM_DF = SD_RPKM.to_frame()
        SD_RPKM_DF.columns = ["SD"]
        number_sd = -1000
        Table_RNA_VST_SD =  Table_RNA_VST.join(SD_RPKM_DF, 
    		how='outer').sort_values(by=['SD'], ascending = 'False').drop(axis=1,
    		labels="SD").iloc[number_sd:,:].transpose()
        Table_RNA_VST_SD.columns = Table_RNA_VST_SD.iloc[0]
        Table_RNA_VST_SD = Table_RNA_VST_SD.iloc[1:,:].reset_index()
        Table_RNA_VST_SD = Table_RNA_VST_SD.rename(columns={'index':'Sample_code_RNAseq'})
        Table_RNA_VST_SD_2 = Table_RNA_VST_SD.merge(Master_table_RNAseq, on='Sample_code_RNAseq', how = 'left')
        subset_table = 'Sample_code_RNAseq'
        features = list(Table_RNA_VST_SD_2.columns.values)
        unwanted = {'index','Dataset','Dataset2','Sample_code','Sample_code_RNAseq','JML','File_name',
        'Type_batch','Type','Delete','Batch','Library','Age','Age_65','Gender','Etiology','Etiology_CAT',
        'Etiology_Mongolia','Etiology_master','Etiology_HBV','Etiology_HDV','Etiology_HCV','Etiology_NonInf',
        'Nodules','Size','Size_CAT5','BCLC','BCLC_CAT','AFP','AFP_CAT400','Fibrosis','Fibrosis_CAT',
        'Fibrosis_CAT_2','Cirrhosis','PCDH7','Order', 'TMB', 'TMB_4', 'TMB_3_5'}
        features = [e for e in features if e not in unwanted]
        DF_toplot_noname = Table_RNA_VST_SD_2.iloc[:,1:]
        DF_toplot = Table_RNA_VST_SD_2
        to_colour_plot = ['Dataset', 'Dataset2',  'Type',  'Batch', 'Library',
           'Age_65', 'Gender', 'Etiology_CAT',
           'Etiology_master', 'Etiology_HBV', 'Etiology_HCV', 
           'Nodules', 'AFP_CAT400', 'Etiology_HDV','PCDH7','Fibrosis_CAT',
        'Fibrosis_CAT_2']
        List_components = ['Comp_1', 'Comp_2', 'Comp_3', 'Comp_4']
        for Feature in to_colour_plot:
            for subset in itertools.combinations(List_components, 2):
                try:
                    Compon1 = subset[0]
                    Compon2 = subset[1]
                    name_master = f'{Feature}_{Compon1}_{Compon2}_PCA'
                    x = DF_toplot_noname.dropna(axis = 'index', subset=[Feature])
                    x = x.loc[:, features].values
                    y = DF_toplot_noname.loc[:,[Feature]].dropna().reset_index().drop('index', axis=1)
                    x = StandardScaler().fit_transform(x)
                    pca = PCA(n_components=4)
                    principalComponents = pca.fit_transform(x)
                    principalDf = pd.DataFrame(data = principalComponents
                                 , columns = List_components)
                    
                    finalDf = pd.concat([principalDf, pd.DataFrame(y)], axis = 1)
                    fig = plt.figure(figsize = (8,8))
                    ax = fig.add_subplot(1,1,1) 
                    ax.set_xlabel(Compon1, fontsize = 15)
                    ax.set_ylabel(Compon2, fontsize = 15)
                    ax.set_title(name_master, fontsize = 20)
                    targets = DF_toplot_noname[Feature].unique()
                    targets = [x for x in targets if str(x) != 'nan']
                    colors = ['r', 'b', 'g', 'y', 'orange']
                    for target, color in zip(targets,colors):
                        names = list(DF_toplot.loc[(DF_toplot[Feature] == target)].Sample_code_RNAseq.values)
                        indicesToKeep = finalDf[Feature] == target
                        z1 = finalDf.loc[indicesToKeep, Compon1].reset_index().drop('index', axis=1)
                        z2 = finalDf.loc[indicesToKeep, Compon2].reset_index().drop('index', axis=1)
                        ax.scatter(z1 , z2
                                   , c = color
                                   , s = 10)
                        for i, txt in enumerate(names):
                            ax.annotate(txt, (z1.loc[i], z2.loc[i]), fontsize = 'xx-small',
                                        horizontalalignment = 'center', color = color , rotation = 0, **font)
                    ax.legend(targets)
                    fig.savefig(f'{folder}{name_master}.pdf')
                    plt.close(fig)
                except:
                    pass
    except:
        pass
