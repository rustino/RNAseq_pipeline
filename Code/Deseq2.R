#!/usr/bin/env Rscript

# Modify acordingly:
# design
# vst or rld normalization

library(RUVSeq)
library(DESeq2)
library(stringr)
library(qdapRegex)


args = commandArgs(trailingOnly=TRUE)

#name_ofcomparison = 'deseq2_CTNNB1_comp_WT_vs_Mutated_Mongolia_TUMOR'
name_ofcomparison = args[1]

#master_folder = '/Users/MiguelT/Box/Marc/Mongolia/RNAseq/Analysis/Deseq2_1_batch_filter_04_07_2020_rld/'

master_folder = args[2]

Input_table <- paste(master_folder,name_ofcomparison,"/","decoder.",name_ofcomparison,".txt",sep="")


name_purified <- qdapRegex::ex_between(name_ofcomparison, "deseq2_", "_comp")[[1]]

reference <- qdapRegex::ex_between(name_ofcomparison, "_comp_", "_vs_")[[1]]

decoder.bySample <- read.table(Input_table,header=T,stringsAsFactors=F)
#directory <- "/Users/MiguelT/Box/Marc/Mongolia/RNAseq/QoRTs_filtered"
directory = args[3]
#counts_deseq <- '/QC.filtered.geneCounts.formatted.for.DESeq.txt.gz'
counts_deseq = args[4]
#model = args[5]
#sampleFiles <- paste0(decoder.bySample$sample.ID,"/QC.geneCounts.formatted.for.DESeq.txt.gz")
sampleFiles <- paste0(decoder.bySample$sample.ID,counts_deseq)
#sampleCondition <- decoder.bySample$group.ID
#sampleCondition_2 <- decoder.bySample$name_purified
sampleCondition_2 <- decoder.bySample[, which(names(decoder.bySample) %in% name_purified)]
sampleName <- decoder.bySample$Sample_code_RNAseq
Library <- decoder.bySample$Library
sampleTable <- data.frame(sampleName = sampleName,
                          fileName = sampleFiles,
                          condition = sampleCondition_2,
                          Library = Library)

#Model for batch correction
#dds <-  DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
#                                   directory = directory,
#                                   design = ~ Library + condition + Library:condition)

#Regular
#dds <-  DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
#                                   directory = directory,
#                                   design = ~ Library + condition)


#Simple Model for batch correction
dds <-  DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                   directory = directory,
                                   design = ~ condition)

dds$condition <- relevel(dds$condition, ref = reference)

dds <- DESeq(dds)
res <- results(dds)
#vst=rlog(dds)
top <- as.data.frame(res)[order(as.data.frame(res)$padj, decreasing=FALSE), ]
original_counts = as.data.frame(assay(dds))
filter <- apply(original_counts, 1, function(x) length(x[x>5])>=2)
filtered <- original_counts[filter,]
genes <- rownames(filtered)[grep("^ENS", rownames(filtered))]
x <- as.factor(Library)
set <- newSeqExpressionSet(as.matrix(original_counts),
                           phenoData = data.frame(x, row.names=colnames(original_counts)))

output_file_5 <- paste(master_folder,name_ofcomparison,"/",name_ofcomparison,"_raw_counts.txt",sep="")

write.table(as.data.frame(counts(set)), file=output_file_5, sep="\t", quote=FALSE, row.names=TRUE, col.names=NA)

Input_figure <- paste(master_folder,name_ofcomparison,"/",name_ofcomparison,"_beforenorm_rle.png",sep="")
png(Input_figure)
try(plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors()[x]))
dev.off()

Input_figure <- paste(master_folder,name_ofcomparison,"/",name_ofcomparison,"_beforenorm_pca.png",sep="")
png(Input_figure)
plotPCA(set, col=colors()[x], cex=1.2)
dev.off()


empirical <- rownames(set)[which(!(rownames(set) %in% rownames(top)[1:10000]))]

set2 <- RUVg(set, empirical, k=1)


output_file_4 <- paste(master_folder,name_ofcomparison,"/",name_ofcomparison,"_normalized_counts.txt",sep="")

write.table(as.data.frame(normCounts(set2)), file=output_file_4, sep="\t", quote=FALSE, row.names=TRUE, col.names=NA)

Input_figure <- paste(master_folder,name_ofcomparison,"/",name_ofcomparison,"_afternorm_rle.png",sep="")
png(Input_figure)
plotRLE(set2, outline=FALSE, ylim=c(-4, 4), col=colors()[x])
dev.off()

Input_figure <- paste(master_folder,name_ofcomparison,"/",name_ofcomparison,"_afternorm_pca.png",sep="")
png(Input_figure)
try(plot(plotPCA(set2, col=colors()[x], cex=1.2)))
dev.off()


pData(set2)$condition <- sampleCondition_2
pData(set2)$Library <- Library


#dds_2 <- DESeqDataSetFromMatrix(countData = normCounts(set2),
#                                colData = pData(set2),
#                                design = ~ condition)

dds_2 <- DESeqDataSetFromMatrix(countData = normCounts(set2),
                                colData = pData(set2),
                                design = ~ condition)

#dds_2 <- DESeqDataSetFromMatrix(countData = normCounts(set2),
#                                colData = pData(set2),
#                                design = ~ Library + condition + Library:condition)

dds_2 <- DESeq(dds_2)

res_2 <- results(dds_2)

vst_2=vst(dds_2)

vst_counts_2=as.data.frame(assay(vst_2))

output_file <- paste(master_folder,name_ofcomparison,"/",name_ofcomparison,"_afternorm_vst.txt",sep="")

write.table(as.data.frame(vst_counts_2), file=output_file, 
            sep="\t", quote=FALSE, row.names=TRUE, col.names=NA)

output_file_2 <- paste(master_folder,name_ofcomparison,"/",'res_total_',name_ofcomparison,"_afternorm.txt",sep="")

write.table(as.data.frame(res_2), file=output_file_2, 
            sep="\t", quote=FALSE, row.names=TRUE, col.names=NA)

output_file_3 <- paste(master_folder,name_ofcomparison,"/",name_ofcomparison,"_afternorm_vst_4dec.txt",sep="")

write.table(as.data.frame(format(vst_counts_2, digits=5)), file=output_file_3, 
            sep="\t", quote=FALSE, row.names=TRUE, col.names=NA)
