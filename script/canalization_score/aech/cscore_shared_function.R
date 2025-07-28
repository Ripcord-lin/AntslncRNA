# 该函数 calc.t 通过以下步骤计算两个组之间基因表达的标准化差异：
# 
# 筛选出符合指定年龄和组别的样本。
# 计算每个基因在两个组中的方差和中位数。
# 计算并返回两个组之间基因表达的标准化差异。
# 标准化差异的计算考虑了差异的方向和大小，以及两个组内的表达方差。这种方法可以在基因表达数据分析中用于识别具有显著表达差异的基因
calc.t <- function(choose.genes,
                   ages,
                   x,y,
                   tpinfo,tpdata){
  x.id = tpinfo %>% filter(Age %in% ages) %>% filter(caste.x == x) %>% rownames()
  y.id = tpinfo %>% filter(Age %in% ages) %>% filter(caste.x == y) %>% rownames()
  a.v <- rowVars(tpdata[choose.genes,x.id],useNames = FALSE)
  a.m <- rowMedians(tpdata[choose.genes,x.id],useNames = FALSE)
  b.v <- rowVars(tpdata[choose.genes,y.id],useNames = FALSE)
  b.m <- rowMedians(tpdata[choose.genes,y.id],useNames = FALSE)
  t <- sign(a.m-b.m)*(abs(a.m-b.m))/(sqrt(a.v+b.v)+0.01)
  return(t)
}
calc.pc <-function(x,n.col,last.t,cor.method){
  y = list()
  y[["p"]]   <- apply(x,1,function(x){a <- cor.test(abs(x)[1:n.col],1:n.col,method=cor.method,alternative="g");return(a[["p.value"]])})
  y[["cor"]] <- apply(x,1,function(x){a <- cor.test(abs(x)[1:n.col],1:n.col,method=cor.method,alternative="g");return(a[["estimate"]])})
  y[["c"]] <- x[,last.t]*-log10(y[["p"]])
  y = as.data.frame(y)
  return(y)
}
neat.merge = function(x,y,...){
  merged.df = merge(x,y,by="row.names",...)
  rownames(merged.df) = merged.df[["Row.names"]]
  merged.df = merged.df[,-1]
  return(merged.df)
}
