# cur_dir = dirname(rstudioapi::getSourceEditorContext()$path)
# setwd(cur_dir)
library(tximport)
library(DESeq2)
library(dplyr)
# raed sample information
tp = read.csv("embryo.1st_RNAseqData_Info.csv")
rownames(tp) = tp$Run
tp$Age = as.factor(tp$Age) 
tp$path = paste("mpha_salmon_embryo_1st_results",tp$Run,"quant.sf",sep = "/")
info.embryo.1st = tp
rm(tp)
# make count matrix
tx2 = read.table("./tx2gene")
txi = tximport(info.embryo.1st[["path"]],tx2gene = tx2,type = "salmon")
colnames(txi$abundance) <- info.embryo.1st$Run
colnames(txi$counts) <- info.embryo.1st$Run
colnames(txi$length) <- info.embryo.1st$Run
Mpha_dev_dds = DESeqDataSetFromTximport(txi,colData = info.embryo.1st,design = ~caste.x)
data.mpha = assay(vst(Mpha_dev_dds))
abundance <- txi$abundance
counts <- txi$counts
colnames(abundance) <- info.embryo.1st$Run

geneinfo = read.csv(header = T,"./Mpha_gene_info_LncAndProteinsFromChromosomes.csv")
used.geneinfo = geneinfo[which(geneinfo$gene %in% rownames(abundance)),]
rm(geneinfo,tx2)

library(preprocessCore)
abundance.log = log2(abundance + 1)
abundance.quantile = normalize.quantiles(as.matrix(abundance.log))
rownames(abundance.quantile) = rownames(abundance)
colnames(abundance.quantile) = colnames(abundance)
abundance.quantile <- as.data.frame(abundance.quantile)
abundance.quantile[which(abundance.quantile < 0)] = 0

print("all.data was filtered by cor > 0.75")
print("abundance : tpm")
print("data.mpha: vst normalized")
print("counts: raw counts")
print("info.embryo.1st: data.info")
print("abundance.quantile:log2, normalize.quantiles")




