library(data.table)
# enable parallel processing for computationally intensive operations.
library(doMC)
registerDoMC(cores = 5)
library("fst")

args = commandArgs(trailingOnly=TRUE) # 

df<-fread(paste("zcat",args[1]),nThread=5)

write.fst(df,paste0(dirname(args[1]),"/",paste0(paste(unlist(strsplit(basename(args[1]),"[.]"))[c(1,2)],collapse="."),".fst")),100)




