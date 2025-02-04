# add std error in significant QTLtools output

args <- commandArgs(trailingOnly = TRUE)
opt_input = args[1]; # QTLtools all chunk file
opt_output = args[2];
input_type = args[3]; # This can be either 19 for raw output of all tests or 21 for the significant ones outputed from Olivier's R script

input_qtltools<-read.table(opt_input,header=F)

ifelse(input_type=="19",names(input_qtltools)<-c("Phenotype_ID","Chromosome_phe","TSS","TSS_end","Strand","Total_no_variants_cis","Distance_tss_variant","Best_variant_in_cis","Chromosome_var",
                         "Pos_variant","Pos_variant_end","DF","Dummy","Beta_dist_1","Beta_dist_2_number_of_ind_tests","Nominal_pvalue","Beta_regression","Empirical_pvalue",
                         "Corrected_pvalue"),
                         names(input_qtltools)<-c("Phenotype_ID","Chromosome_phe","TSS","TSS_end","Strand","Total_no_variants_cis","Distance_tss_variant","Best_variant_in_cis","Chromosome_var",
                                                  "Pos_variant","Pos_variant_end","DF","Dummy","Beta_dist_1","Beta_dist_2_number_of_ind_tests","Nominal_pvalue","Beta_regression","Empirical_pvalue",
                                                  "Corrected_pvalue","qvalue","nthreshold"))

# Calculate std error from pvalues and betas by calculating t-values and
#coef(y.x)[2,1]/(sign(coef(y.x)[2,1])*sqrt(qf(coef(y.x)[2,4],1,98,lower=F)))
#abs(beta/sqrt(qf(pvalue,1,degrees_freedom,lower=F)))

input_qtltools$std.err<-abs(input_qtltools$Beta_regression/sqrt(qf(input_qtltools$Nominal_pvalue,1,input_qtltools$DF,lower=F)))

fout1=paste(opt_output, "stderr.txt", sep=".")
write.table(input_qtltools, fout1, quote=FALSE, row.names=FALSE, col.names=TRUE,sep="\t")
