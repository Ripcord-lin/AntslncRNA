# cur_dir = dirname(rstudioapi::getSourceEditorContext()$path)
# setwd(cur_dir)
library(tximport)
library(DESeq2)
library(dplyr)
# raed sample information
tp = read.csv("aech_raw_sampleinfo.csv")
rownames(tp) = tp$Run
exist_sample <- list.files("./salmon_results", full.names = FALSE)
tp <- tp[exist_sample,]

tp$Age = as.factor(tp$Age) 
tp$caste.x = as.factor(tp$caste.x)
tp$path = paste("./salmon_results",tp$Run,"quant.sf",sep = "/")
info.aech = tp
rm(tp)
# make count matrix
tx2 = read.table("./tx2gene")
txi = tximport(info.aech[["path"]],tx2gene = tx2,type = "salmon")
colnames(txi$abundance) <- info.aech$Run
colnames(txi$counts) <- info.aech$Run
colnames(txi$length) <- info.aech$Run
# DESeq
Aech_dev_dds = DESeqDataSetFromTximport(txi,colData = info.aech,design = ~caste.x)
data.aech = assay(vst(Aech_dev_dds))
abundance <- txi$abundance
counts <- txi$counts
colnames(abundance) <- info.aech$Run

geneinfo = read.csv(header = T,"./Aech_gene_info_LncAndProteinsFromChromosomes.csv")
used.geneinfo = geneinfo[which(geneinfo$gene %in% rownames(data.aech)),]

info.aech$Group <- paste0(info.aech$caste.x, "_", info.aech$Age)
group_levels <- unique(info.aech$Group)

library(preprocessCore)
abundance.log = log2(abundance + 1)
abundance.quantile = normalize.quantiles(as.matrix(abundance.log))
rownames(abundance.quantile) = rownames(abundance)
colnames(abundance.quantile) = colnames(abundance)
abundance.quantile <- as.data.frame(abundance.quantile)
abundance.quantile[which(abundance.quantile < 0)] = 0

print("data.aech: vst normalized")
print("abundance: tpm")
print("counts: raw counts")
print("abundance.quantile:log2, normalize.quantiles")
print("info.aech: sample info")


