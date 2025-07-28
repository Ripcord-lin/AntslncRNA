cur_dir = dirname(rstudioapi::getSourceEditorContext()$path)
setwd(cur_dir)
source("./raw_data_process/mpha.dev.RNAseq/ini.R",chdir = T)
source("cscore_shared_function.R")
library(FactoMineR)
library(dplyr)
library(ggplot2)
library(DESeq2)

lnc.used.geneinfo = used.geneinfo %>% filter(identity == "lncRNA")
lnc <- lnc.used.geneinfo$gene %>% unlist() %>% unique()
tpdata = data.mpha[lnc,]

info.mpha$Age = info.mpha$Age %>% as.character() %>% as.factor()
Mpha_dev_dds = DESeqDataSetFromTximport(txi,colData = info.mpha,
                                        design = ~caste.x)
Mpha_dev_dds = Mpha_dev_dds[rownames(tpdata),] # choose lnc only
Mpha_dev_dds = DESeq(Mpha_dev_dds,test = "LRT",reduced = ~1)
# summary(res.dev)
res.dev <- results(Mpha_dev_dds)
res.dev.signif.RNZ = res.dev %>% as.data.frame()
res.dev.signif.RNZ = res.dev.signif.RNZ %>% filter(padj < 0.0001)
tpgene = res.dev.signif.RNZ %>% rownames()

# calculation
t.res = list()
info.mpha$Age = ordered(info.mpha$Age,levels = c("2nd","3rd",
                                                 "Pre.Pupa","Pupa.Young","Pupa.Old","Imago"))
for (tp_age in levels(info.mpha$Age)){
  t.res[[tp_age]] = calc.t(choose.genes = tpgene,
                           ages = tp_age,
                           x="Gyne",y="Worker",
                           tpinfo=info.mpha,
                           tpdata = data.mpha)
}
t.res = as.data.frame(t.res)
rownames(t.res) = tpgene
c.res = calc.pc(t.res,6,last.t = 5,cor.method = "p")
c.table = neat.merge(t.res,c.res)
rm(t.res,c.res)

c.table$bias = "n"
c.table[which(c.table$c > 1),"bias"] = "gyne"
c.table[which(c.table$c < -1),"bias"] = "worker"

print(paste0("All lncRNAs: ", nrow(tpdata),"; DE lncs: ", length(tpgene), "/", nrow(tpdata), ";",
             " Queen canalized: ", nrow(c.table[which(c.table$c > 1),]), "/", length(tpgene),";",
             " Worker canalized: ", nrow(c.table[which(c.table$c < -1),]), "/", length(tpgene)))
write.csv(c.table,file = "lncRNA_canalization.csv")
