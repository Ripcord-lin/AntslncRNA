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
