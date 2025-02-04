# add std error in significant QTLtools output
library(data.table)

args <- commandArgs(trailingOnly = TRUE)
opt_input = args[1]; # input all
opt_input_perm<-args[2] # input from standards QTL mapping to calculate degrees of freedom
opt_output = args[3]; # output file

DF<-read.table(opt_input_perm)
DF<-DF[!is.na(DF$V12),]
DF<-unique(DF$V12)

input_qtltools<-fread(input=paste("zcat",opt_input))


names(input_qtltools)<-c("Phenotype_ID","Chromosome_phe","TSS","TSS_end","Strand","Total_no_variants_cis","Distance_tss_variant","Variant_in_cis","Chromosome_var","Pos_variant","Pos_variant_end","Nominal_pvalue","Beta_regression","Best_hit")

# Calculate std error from pvalues and betas by calculating t-values and
#coef(y.x)[2,1]/(sign(coef(y.x)[2,1])*sqrt(qf(coef(y.x)[2,4],1,98,lower=F)))
#abs(beta/sqrt(qf(pvalue,1,degrees_freedom,lower=F)))


input_qtltools$std.err<-abs(input_qtltools$Beta_regression/sqrt(qf(input_qtltools$Nominal_pvalue,1,DF,lower=F)))
input_qtltools$std.err[(is.na(input_qtltools$std.err))] <-0
fout1=paste(opt_output, "stderr.rds", sep=".")
fout2=paste(opt_output, "stderr.txt", sep=".")
fout3=paste(opt_output, "stderr_5col.txt", sep=".")

#saveRDS(input_qtltools, fout1)

input_qtltools$Variant_in_cis<-paste(input_qtltools$Variant_in_cis,input_qtltools$Chromosome_var,input_qtltools$Pos_variant,sep="_")
input_qtltools$gene_snp<-paste(input_qtltools$Phenotype_ID,input_qtltools$Variant_in_cis,sep="_")

input_qtltools<-input_qtltools[!(duplicated(input_qtltools$gene_snp)),]


write.table(input_qtltools,fout2,sep="\t",quote=FALSE,col.names=TRUE,row.names=FALSE)
write.table(input_qtltools[,c(1,8,13,15,12)],fout3,sep="\t",quote=FALSE,col.names=TRUE,row.names=FALSE)

system (paste("gzip -f",fout2))
system (paste("gzip -f",fout3))
