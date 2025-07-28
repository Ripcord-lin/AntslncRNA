vis.tissue.genes = function(gene = "Mpha_lnc00484",display = "v",title = "title="){
  library(ggplot2)
  library(ggstatsplot)
  rownames(info.tissue) = info.tissue$ID
  iptdt = data.frame(Exp = data.tissue.qn.log2[gene,] %>% as.numeric(),
                     ID = colnames(data.tissue.qn.log2),
                     Organ = info.tissue[colnames(data.tissue.qn.log2),"organ"])
  if (display == "v"){
    ggplot(data = iptdt,
           aes(x=Organ,y=Exp,color=Organ))+
      geom_violin(alpha=.2,aes(fill = Organ))+
      geom_point(alpha=.5)+
      ylab("Exp: log2(TPM)")+
      theme_classic()+
      ggtitle(title)
  } else if (display == "b") {
    ggplot(data = iptdt,
           aes(x=Organ,y=Exp,color=Organ))+
      geom_boxplot(width=.5,
                   #color="black",fill="white",
                   outlier.shape = NA)+
      geom_point(alpha=.5)+
      ylab("Exp: log2(TPM)")+
      theme_classic()+
      ggtitle(title)
  } else if (display == "s"){
    ggbetweenstats(data = iptdt,
                   x=Organ,
                   y=Exp)
  } else {
    cat("Please specify 'display'.\nChoose from: v (violin plot), b (boxplot), s (show stats)")
  }
}

pca.enrichment <- function(){
  gene.set.sum = colSums(data.tissue[gene.id,]) # use abundance
  bkgd.sum = colSums(data.tissue)
  enrich.score = gene.set.sum/bkgd.sum
  # PCA
  pca.tissue <- PCA(t(data.tissue.TPM.QN),graph = F)
  dim1 = 1
  dim2 = 2
  var1 = pca.tissue$eig[dim1,2] %>% round(2)
  var2 = pca.tissue$eig[dim2,2] %>% round(2)
  iptdt = cbind(pca.tissue$ind$coord,info.tissue)
  iptdt$enrich = enrich.score
  ggplot(iptdt,aes_string(x=paste0("Dim.",dim1),y=paste0("Dim.",dim2),
                          color="enrich"))+
    stat_ellipse(aes(group=organ),color="grey",alpha=.5,level = 0.95)+
    geom_point(alpha=1)+
    scale_colour_gradientn(colours = c("yellow", "red","darkred"),
                           name = "Enrichment Score")+
    xlab(paste0("Dim.",dim1,": ", var1,"%"))+
    ylab(paste0("Dim.",dim2,": ", var2,"%"))+
    theme_classic()
}