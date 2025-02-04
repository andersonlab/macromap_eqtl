# read file 
#temp = list.files(pattern="*.txt")
#myfiles = lapply(temp, read.table)
#data_var_dec<-do.call(cbind.data.frame,myfiles)

#pdf()
#par(mar = c(7,12,2,2) + 0.1)
#boxplot(t(data_var_dec)[,order(apply(data_var_dec,1,median),decreasing=F)],las=2,horizontal = T)

library("RColorBrewer")
a<-readRDS("LFLMM.results.RDS")
source("Barplot.LFLMM.R")

pdf("var_dec.pdf")
par(mar = c(12,5,2,2) + 0.1)
Barplot(a,col=brewer.pal(8,"Set1")[2])
dev.off()





