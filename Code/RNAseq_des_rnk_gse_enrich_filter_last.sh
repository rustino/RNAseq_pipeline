#!/bin/bash

NewVar=$(echo '/Users/MiguelT/Box/Marc/Mongolia/RNAseq/Analysis')
NewVar2=$(echo '/Users/MiguelT/Box/Marc/Mongolia/RNAseq/QoRTs')
NewVar3=$(echo '/QC.geneCounts.formatted.for.DESeq.txt.gz')
NewVar4=$(echo '/Users/MiguelT/Box/My_scripts/Mongolia_RNAseq')
NewVar5=$(echo 'Deseq2_1_batch_filter_10_13_2020_vst')
NewVar6=$(echo '/Users/MiguelT/Box/My_scripts/Mongolia_RNAseq/Deseq2_1_batch_filter_10_13_2020_vst/')


mkdir -p $NewVar/$NewVar5/

source activate RNAseq_mac_2

python $NewVar4/$NewVar5/01_Decoder_filtering_2.py --folder $NewVar/$NewVar5/ --Qfold $NewVar2 --fext $NewVar3

echo Done 01_Decoder

conda deactivate 

source activate Deseq2_4

python $NewVar4/$NewVar5/Deseq2_plus_rnk.py --folder $NewVar/$NewVar5/ --Qfold $NewVar2 --fext $NewVar3  --folder_script $NewVar6

echo Done Deseq2_plus_rnk

conda deactivate 

source activate RNAseq_mac_2

python $NewVar4/$NewVar5/reduceSD.py --folder $NewVar/$NewVar5/

echo done reduceSD

python $NewVar4/$NewVar5/PCA.py --folder $NewVar/$NewVar5/

echo Done PCA

conda deactivate 

source activate gseapy

python $NewVar4/$NewVar5/enrichr_only.py --folder $NewVar/$NewVar5/

echo Done enrichr_only

python $NewVar4/$NewVar5/gsea_only.py --folder $NewVar/$NewVar5/

echo Done gseaonly

conda deactivate

echo Done all the stuff!