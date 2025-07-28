# cur_dir = dirname(rstudioapi::getSourceEditorContext()$path)
# setwd(cur_dir)

library(dplyr)
library(tximport)
# function
quick.qn = function(x,log2.transform){
  library(preprocessCore)
  tp = normalize.quantiles(x %>% as.matrix())
  rownames(tp) = rownames(x)
  colnames(tp) = colnames(x)
  x = tp
  x = as.data.frame(x)
  if (log2.transform){
    x = log2(x + 1)
  }
  return(x)
}

# tissue data metadata
neat.get.tissue.info <- function(raw.info.table, salmon.res.dir){
  tp = read.csv(raw.info.table)
  rownames(tp) = tp$V1
  tp$path = paste(salmon.res.dir,tp$V1,"quant.sf",sep = "/")
  tp = tp %>% select(V2.y,organ,path)
  colnames(tp)[1] <- "ID"
  return(tp)
}
info.tissue <- neat.get.tissue.info(raw.info.table = "info_gyne_tissue.csv",
                                    salmon.res.dir = "tissue_salmon_res")

tx2 = read.table("./tx2gene")
txi.tissue = tximport(info.tissue$path,tx2gene = tx2,type = "salmon")
data.tissue.tpm = txi.tissue$abundance
colnames(data.tissue.tpm) = info.tissue$ID
data.tissue.counts = txi.tissue$counts
colnames(data.tissue.counts) = info.tissue$ID
data.tissue.qn.log2 = quick.qn(data.tissue.tpm,log2.transform = T)

# gene metadata
info.gene = read.csv("./Mpha_gene_info_LncAndProteinsFromChromosomes.csv")
used.geneinfo = info.gene[which(info.gene$gene %in% rownames(data.tissue.counts)),]
# rm(list = ls()[-which(ls() %in% c("data.tissue.qn.log2","info.tissue","info.gene","quick.qn"))])

