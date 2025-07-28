cur_dir = dirname(rstudioapi::getSourceEditorContext()$path)
setwd(cur_dir)
library(EDASeq)
source("../raw_data_process/mpha.JH.RNAseq/ini.R",chdir = T)

lnc = used.geneinfo %>% filter(identity == "lncRNA") %>% select(gene) %>% unlist() %>% unique()
tpdata = data.mpha[lnc,]



used.sample.info = read.csv("used.sample.info.csv")
rownames(used.sample.info) = used.sample.info$run
used.sample.info$age = as.factor(used.sample.info$age) 
used.sample.info$group1 = as.factor(used.sample.info$group1)

gyne6D_sample <- rownames(used.sample.info[which(used.sample.info$group1 == "Control.Gyne" & used.sample.info$age == "3rd.6D"),])
gyne6D_matrix = tpdata[,gyne6D_sample]
colnames <- colnames(gyne6D_matrix)
new_colnames <- gsub("6D", "5D", colnames)
colnames(gyne6D_matrix) <- new_colnames
tpdata_add_newcol <- cbind(tpdata,gyne6D_matrix)
new_colnames <- gsub("6D", "10D", colnames)
colnames(gyne6D_matrix) <- new_colnames
tpdata_add_newcol <- cbind(tpdata_add_newcol,gyne6D_matrix)

used.sample.info <- used.sample.info[!(rownames(used.sample.info) %in% gyne6D_sample),]


age.levels <- factor(c("3rd.24h","3rd.5D","3rd.10D","Pre.Pupa","Young.Pupa"))

gene_stage_gyne_expression_df <- data.frame(matrix(NA, nrow = length(rownames(tpdata)),
                                                   ncol = length(age.levels)))
rownames(gene_stage_gyne_expression_df) <- rownames(tpdata)
colnames(gene_stage_gyne_expression_df) <- age.levels
# worker
gene_stage_worker_expression_df <- data.frame(matrix(NA, nrow = length(rownames(tpdata)),
                                                     ncol = length(age.levels)))
rownames(gene_stage_worker_expression_df) <- rownames(tpdata)
colnames(gene_stage_worker_expression_df) <- age.levels

for (period in unique(age.levels)) {
  caste = "Gyne"
  current_samples <- rownames(used.sample.info[which(used.sample.info$age == period & used.sample.info$caste == caste & used.sample.info$treatment == 'Control'),])
  period_expression <- as.matrix(tpdata_add_newcol[, current_samples])
  gene_stage_gyne_expression_df[,period] <- apply(period_expression, 1, mean)
  
  caste = "Worker"
  current_samples <- rownames(used.sample.info[which(used.sample.info$age == period & used.sample.info$caste == caste & used.sample.info$treatment == 'Control'),])
  period_expression <- as.matrix(tpdata_add_newcol[, current_samples])
  period_expression <- as.matrix(tpdata_add_newcol[, current_samples])
  gene_stage_worker_expression_df[,period] <- apply(period_expression, 1, mean)
  
}
library(EDASeq)

dist_df <- data.frame(matrix(NA, nrow = length(rownames(used.sample.info)),
                             ncol = 4))
rownames(dist_df) <- rownames(used.sample.info)
colnames(dist_df) <- c("dist2worker","dist2gyne","distworker2gyne","cts_score")

for (sample in rownames(used.sample.info)) {
  group = used.sample.info[sample,"group1"]
  age = as.character(used.sample.info[sample,'age'])
  dist2worker = dist(t(cbind(gene_stage_worker_expression_df[,age],tpdata_add_newcol[,sample])),method = 'euclidean')
  dist2gyne = dist(t(cbind(gene_stage_gyne_expression_df[,age],tpdata_add_newcol[,sample])),method = 'euclidean')
  distworker2gyne = dist(t(cbind(gene_stage_worker_expression_df[,age],gene_stage_gyne_expression_df[,age])),method = 'euclidean')
  cts_score = (dist2worker - dist2gyne)/distworker2gyne
  dist_df[sample, ] = c(dist2worker,dist2gyne, distworker2gyne, cts_score)
  
  
}
