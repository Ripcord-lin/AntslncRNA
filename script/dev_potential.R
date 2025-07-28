cur_dir = dirname(rstudioapi::getSourceEditorContext()$path)
setwd(cur_dir)
source("./raw_data_process/mpha.dev.RNAseq/ini.R",chdir = T)

# lncRNAs
lnc = used.geneinfo %>% filter(identity == "lncRNA") %>% select(gene) %>% unlist() %>% unique()
tpdata = data.mpha[lnc,]
# or protein
# pro = used.geneinfo %>% filter(identity == "protein") %>% select(gene) %>% unlist() %>% unique()
# tpdata = data.mpha[pro,]
info.mpha$Age <- factor(info.mpha$Age, levels = c("2nd","3rd",
                                                  "Pre.Pupa","Pupa.Young","Pupa.Old","Imago"))
library(sva)
iptdt = data.frame(sampleID=character(), x.coord = numeric(), 
                   Body_length = numeric(), Age = character(), Caste = character())
set.seed(1234) # boot strap is random
for (i in 1:5){
  current_age = levels(info.mpha$Age)[i]
  next_age = levels(info.mpha$Age)[i+1]
  n_gene = 6000
  lpdata = tpdata[order(rowSums(tpdata),decreasing = T),][1:n_gene,]
  n_bootstrap = 100
  current_age_samples = info.mpha %>% filter(Age == current_age)
  a = rownames(current_age_samples) %>% 
    sample((n_bootstrap - nrow(current_age_samples)),replace = T)
  a = c(rownames(current_age_samples),a) 
  current_age_data = lpdata[,a]
  next_age_samples = info.mpha %>% filter(Age == next_age)
  b = rownames(next_age_samples) %>% 
    sample((n_bootstrap - nrow(next_age_samples)),replace = T)
  b = c(rownames(next_age_samples),b)
  next_age_data = lpdata[,b]
  edata = cbind(current_age_data,next_age_data)
  age_batch = rep(c(1,2),c(n_bootstrap,n_bootstrap))
  combat_data = ComBat(dat = edata,batch = age_batch,mean.only = F)
  current_age_data_combat = combat_data[,1:100]
  next_age_data_combat = combat_data[,101:200]
  a = next_age_samples %>% filter(caste.x == "Worker")
  next_avg_w = next_age_data_combat[,which(colnames(next_age_data_combat) %in% rownames(a))] %>% rowMedians()
  a = next_age_samples %>% filter(caste.x == "Gyne")
  next_avg_g = next_age_data_combat[,which(colnames(next_age_data_combat) %in% rownames(a))] %>% rowMedians()
  cor.method = "p"
  dist_w.g = dist(t(cbind(next_avg_g,next_avg_w)),,method = "man")
  dist_w.g 
  x.coord = c()
  a = current_age_samples %>% filter(caste.x == "Worker")
  for (sp in rownames(a)){
    sp_exp = current_age_data_combat[,sp]
    cor.w = dist(t(cbind(sp_exp,next_avg_w)),method = "man")
    cor.g = dist(t(cbind(sp_exp,next_avg_g)),method = "man")
    x.coord = c(x.coord,(cor.w - cor.g)/dist_w.g)
  }
  b = current_age_samples %>% filter(caste.x == "Gyne")
  for (sp in rownames(b)){
    sp_exp = current_age_data_combat[,sp]
    cor.w = dist(t(cbind(sp_exp,next_avg_w)),method = "man")
    cor.g = dist(t(cbind(sp_exp,next_avg_g)),method = "man")
    x.coord = c(x.coord,(cor.w - cor.g)/dist_w.g)
  }
  iptdt = rbind(iptdt,
                data.frame(sampleID = c(rownames(a),rownames(b)),
                           x.coord = x.coord,
                           Body_length = c(a$body_length.x,b$body_length.x),
                           Age = current_age,
                           Caste = c(rep("Worker",nrow(a)),rep("Gyne",nrow(b)))
                           )
                )
}
iptdt$Age = ordered(iptdt$Age,levels = c("2nd","3rd",
                                            "Pre.Pupa","Pupa.Young","Pupa.Old"))

