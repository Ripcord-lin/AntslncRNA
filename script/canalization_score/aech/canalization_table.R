cur_dir = dirname(rstudioapi::getSourceEditorContext()$path)
setwd(cur_dir)
source("../../raw_data_process/aech.dev.RNAseq/ini.R",chdir = T)
source("cscore_shared_function.R")
library(FactoMineR)
library(dplyr)
library(ggplot2)
library(DESeq2)

lnc.used.geneinfo = used.geneinfo %>% filter(identity == "lncRNA")
lnc <- lnc.used.geneinfo$gene %>% unlist() %>% unique()
tpdata = data.aech[lnc,]

info.aech <- info.aech[info.aech$Age != '1st',]
info.aech <- info.aech[info.aech$Age != '2nd',]
info.aech <- info.aech[info.aech$Age != '4th',]
info.aech <- info.aech[info.aech$caste.x != 'Worker.Large',]
info.aech <- info.aech[info.aech$caste.x != 'Worker.Medium',]
tpdata = tpdata[,info.aech$Run]

info.aech$Age <- factor(info.aech$Age, levels = c("3rd","Pre.Pupa",
                                                  "Pupa.Young","Pupa.Old","Imago"))

filtered_txi <- list(
  counts = txi$counts[rownames(tpdata), info.aech$Run],
  abundance = txi$abundance[rownames(tpdata), info.aech$Run],
  length = txi$length[rownames(tpdata), info.aech$Run],
  countsFromAbundance = txi$countsFromAbundance
)

info.aech$Age = info.aech$Age %>% as.character() %>% as.factor()
Mpha_dev_dds = DESeqDataSetFromTximport(filtered_txi,colData = info.aech,
                                        design = ~caste.x)
Mpha_dev_dds = Mpha_dev_dds[rownames(tpdata),] # choose lnc only
Mpha_dev_dds = DESeq(Mpha_dev_dds,test = "LRT",reduced = ~1)
res.dev <- results(Mpha_dev_dds)
res.dev.signif.RNZ = res.dev %>% as.data.frame()
res.dev.signif.RNZ = res.dev.signif.RNZ %>% filter(padj < 0.0001)
tpgene = res.dev.signif.RNZ %>% rownames()

t.res = list()
info.aech$Age = ordered(info.aech$Age,levels = c("3rd",
                                                 "Pre.Pupa","Pupa.Young","Pupa.Old","Imago"))
for (tp_age in levels(info.aech$Age)){
  t.res[[tp_age]] = calc.t(choose.genes = tpgene,
                           ages = tp_age,
                           x="Gyne",y="Worker.Small",
                           tpinfo=info.aech,
                           tpdata = data.aech)
}
t.res = as.data.frame(t.res)
rownames(t.res) = tpgene
c.res = calc.pc(t.res,6,last.t = 5,cor.method = "p")
c.table = neat.merge(t.res,c.res)
c.table["Aech_lnc12873",]
c.table["Aech_lnc08121",]


c.table$bias = "n"
c.table[which(c.table$c > 1),"bias"] = "gyne"
c.table[which(c.table$c < -1),"bias"] = "worker.small"

print(paste0("All lncRNAs: ", nrow(tpdata),"; DE lncs: ", length(tpgene), "/", nrow(tpdata), ";",
             " gyne canalized: ", nrow(c.table[which(c.table$c > 1),]), "/", length(tpgene),";",
             " Worker canalized: ", nrow(c.table[which(c.table$c < -1),]), "/", length(tpgene)))
# write.csv(c.table,file = "lncRNA_canalization.csv")
