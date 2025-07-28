cur_dir = dirname(rstudioapi::getSourceEditorContext()$path)
setwd(cur_dir)
library(tidyverse)
library(stringr)
library(KEGGREST)
library(AnnotationForge)
library(clusterProfiler)
library(dplyr)
library(jsonlite)
library(purrr)
library(RCurl)

emapper <- read.table("Mpha.emapper.annotations.tsv", header=TRUE, sep = "\t",quote = "")
emapper[emapper==""]<-NA

gene_info <- emapper %>% dplyr::select(GID = query, GENENAME = Preferred_name) %>% na.omit() 
gos <- emapper %>% dplyr::select(query, GOs) %>% na.omit()
gene2go = data.frame(GID = character(),
                     GO = character(),
                     EVIDENCE = character())
gos_list <- function(x){
  the_gos <- str_split(x[2], ",", simplify = FALSE)[[1]] 
  df_temp <- data.frame(GID = rep(x[1], length(the_gos)), 
                        GO = the_gos, 
                        EVIDENCE = rep("IEA", length(the_gos))) 
  return(df_temp) 
} 
gene2gol <- apply(as.matrix(gos),1,gos_list) 
gene2gol_df <- do.call(rbind.data.frame, gene2gol) 
gene2go <- gene2gol_df 
gene2go$GO[gene2go$GO=="-"]<-NA 
gene2go<-na.omit(gene2go)

gene2ko <- emapper %>% dplyr::select(GID = query, Ko = KEGG_ko) 
gene2ko$Ko[gene2ko$Ko=="-"]<-NA 
gene2ko<-na.omit(gene2ko) 
gene2kol <- apply(as.matrix(gene2ko),1,gos_list) 
gene2kol_df <- do.call(rbind.data.frame, gene2kol) 
gene2ko <- gene2kol_df[,1:2] 
colnames(gene2ko) <- c("GID","Ko") 
gene2ko$Ko <- gsub("ko:","",gene2ko$Ko)


update_kegg <- function(json = "ko00001.json") { 
  pathway2name <- tibble(Pathway = character(), Name = character()) 
  ko2pathway <- tibble(Ko = character(), Pathway = character()) 
  kegg <- fromJSON(json) 
  for (a in seq_along(kegg[["children"]][["children"]])) { 
    A <- kegg[["children"]][["name"]][[a]] 
    for (b in seq_along(kegg[["children"]][["children"]][[a]][["children"]])) { 
      B <- kegg[["children"]][["children"]][[a]][["name"]][[b]] 
      for (c in seq_along(kegg[["children"]][["children"]][[a]][["children"]][[b]][["children"]])) { 
        pathway_info <- kegg[["children"]][["children"]][[a]][["children"]][[b]][["name"]][[c]] 
        pathway_id <- str_match(pathway_info, "ko[0-9]{5}")[1] 
        pathway_name <- str_replace(pathway_info, " \\[PATH:ko[0-9]{5}\\]", "") %>% str_replace("[0-9]{5} ", "") 
        pathway2name <- rbind(pathway2name, tibble(Pathway = pathway_id, Name = pathway_name))
        kos_info <- kegg[["children"]][["children"]][[a]][["children"]][[b]][["children"]][[c]][["name"]] 
        kos <- str_match(kos_info, "K[0-9]*")[,1] 
        ko2pathway <- rbind(ko2pathway, tibble(Ko = kos, Pathway = rep(pathway_id, length(kos))))}}} 
  save(pathway2name, ko2pathway, file = "kegg_info.RData")} 
update_kegg() 
load(file = "kegg_info.RData")

gene2pathway <- gene2ko %>% left_join(ko2pathway, by = "Ko") %>% dplyr::select(GID, Pathway) %>% na.omit()

tax_id = "307658"
genus = "Monomorium"
species = "pharaonisncbi"
gene2go <- unique(gene2go) 
gene2go <- gene2go[!duplicated(gene2go),] 
gene2ko <- gene2ko[!duplicated(gene2ko),] 
gene2pathway <- gene2pathway[!duplicated(gene2pathway),]
gene_info <- gene_info[!duplicated(gene_info),]

makeOrgPackage(gene_info=gene_info,
               go=gene2go,
               # ko=gene2ko,
               maintainer='',
               author='',
               # pathway=gene2pathway,
               version="0.0.1",
               outputDir = "./build", #输出库的位置
               tax_id=tax_id,
               genus=genus,
               species=species,
               goTable="go")

install.packages("./build/org.Mpharaonisncbi.eg.db",repos = NULL, type="source") 
library("org.Mpharaonisncbi.eg.db")
columns(org.Mpharaonisncbi.eg.db) 
keys(org.Mpharaonisncbi.eg.db) 
pathway2name$Name <- gsub(" \\[BR:ko[0-9]{5}\\]", "",pathway2name$Name) 
pathway2name<- na.omit(pathway2name) 
pathway2gene <-gene2pathway[, c("Pathway","GID")] 
write.table(pathway2name,file = "./pathway2name", sep = "\t", quote = F, row.names = F) 
write.table(pathway2gene,file = "./pathway2gene", sep = "\t", quote = F, row.names = F)