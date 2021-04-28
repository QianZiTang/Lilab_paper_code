library(MatchIt)
setwd('/Users/hsl/Desktop/forJinplot/addHiC/finalNEEcodes/')
data<-read.table(file='formatch.nochicken.nomouse.correct.filtered.xls',sep='\t',header=T,row.names=1)

m.out <- matchit(treat ~ mymeanexprs, data = data, method = "nearest",caliper=0.001,ratio=1)
## m.out <- matchit(treat ~ mymeanexprs, data = data, method = "nearest")
m.data <- match.data(m.out)

write.table(m.data,file="formatch_results.nochicken.nomouse.correct.filtered.xls",sep="\t",quote=F,row.names=T,col.names=T)
