# Run this script from mac not from cluster 
library(locuscomparer)
library(stringr)
library(dplyr)

rm(list=ls())
## load data. The RDS is from figure 4 code 
df_ds_final_07_05<-readRDS("~/Desktop/FARM/myscratch//MacroMap/Data/MacroMap_manuscript_figures/data/df_ds_final_075macromap_05GTEx.RDS")
names(df_ds_final_07_05)<-lapply(df_ds_final_07_05, function (x) unique(x$gwas_trait))
setwd("/Users/np12/Documents/PROJECTS/MacroMap/coloc_results_df_ds_final_07_05")
source("~/Desktop/locus_compare_functions.R")

# in_fn2 must always be the QTL 
locuscompare2<-function (in_fn1, in_fn2, marker_col1 = "rsid", pval_col1 = "pval", 
                         title1 = "eQTL", marker_col2 = "rsid", pval_col2 = "pval", 
                         title2 = "GWAS", snp = NULL, population = "EUR", combine = TRUE, 
                         legend = TRUE, legend_position = c("bottomright", "topright", 
                                                            "topleft"), lz_ylab_linebreak = FALSE, genome = c("hg19","hg38"))
{
  d1 = read_metal(in_fn1, marker_col1, pval_col1)
  d2 = read_metal(in_fn2, marker_col2, pval_col2)
  merged = merge(d1, d2, by = "rsid", suffixes = c("1", "2"), 
                 all = FALSE)
  genome = match.arg(genome)
  merged = get_position(merged, genome)
  chr = unique(merged$chr)
  if (length(chr) != 1) 
    stop("There must be one and only one chromosome.")
  lead_snp<-data.frame(merge(d1[d1$rsid==snp,],d2[d2$rsid==snp,],by = "rsid", suffixes = c("1", "2"),all=FALSE),chr=substr(in_fn2$V9,4,6)[which(in_fn2$rsid==snp)],pos=in_fn2$V10[which(in_fn2$rsid==snp)])
  merged<-rbind(merged,lead_snp)
  merged<-merged[!duplicated(merged$rsid),]
  snp = get_lead_snp(merged, snp)
  ld = retrieve_LD(chr, snp, population)
  p = make_combined_plot(merged, title1, title2, ld, chr, snp, 
                         combine, legend, legend_position, lz_ylab_linebreak)
  return(p)
}

# debug 
#df_ds_final_07_05<-df_ds_final_07_05[c(1)]
#x<-df_ds_final_07_05$IBD

# you can plot all of them with the code below and check the CTSA and CAD plots 

#############################################################################################################################################
lapply(df_ds_final_07_05, function (x) {
  if  (length(x$study) !=0) {
    x1<-x %>% filter(study=="MacroMap") %>% filter (!condition_name %in% c("Ctrl_6","Ctrl_24"))  %>% filter (PP.H4.abf >0.5)
    x2<-x %>% filter(study=="MacroMap") %>% filter (condition_name %in% c("Ctrl_6","Ctrl_24"))
    x3<-rbind(x1,x2)
    trait<-unique(x$gwas_trait)
    genes<-unique(x$HGNC)
    system(paste("mkdir",trait))
    sapply(genes, function (x) dir.create(path=paste0(trait,"/",x),recursive = TRUE))
    
    y<-genes
    #i<-1
    sapply(genes, function (y) { 
      x4<-x3[x3$HGNC==y,]
      
      for (i in 1:dim(x4)[1]) {
        z<-x4[i,]
        qtl<-read.table(text=system(paste("/Users/np12/homebrew/bin/tabix",str_replace(z$qtl_path,"~","~/Desktop/FARM"),paste0(z$Chromosome_phe,":",z$TSS-100,"-",z$TSS+100)),intern = TRUE))
        gwas<-read.table(text=system(paste("/Users/np12/homebrew/bin/tabix",str_replace(z$gwas_path,"~","~/Desktop/FARM"),paste0(substr(z$Chromosome_phe,4,6),":",z$TSS-2000000,"-",z$TSS + 2000000)),intern = TRUE))
        gwas_header<-read.table(text=system(paste("/usr/bin/gzcat",str_replace(z$gwas_path,"~","~/Desktop/FARM"),paste0("|head -1")),intern = TRUE),header = TRUE)
        colnames(gwas)<-colnames(gwas_header)
        qtl$rsid<-as.character(str_split_fixed(qtl$V8, "_",n=2)[,1])
        if(unlist(str_split(z$gwas_path,"/"))[9]=="harmonized") {
          gwas<-gwas[paste(gwas$hm_chrom,gwas$hm_pos,sep="_") %in% paste(substr(qtl$V9,4,6),qtl$V10,sep="_"),]
          qtl<-qtl[paste(substr(qtl$V9,4,6),qtl$V10,sep="_") %in% paste(gwas$hm_chrom,gwas$hm_pos,sep="_"),]
          
          qtl<-qtl[!duplicated(qtl$V10),]
          gwas<-gwas[!duplicated(gwas$hm_pos),]
          
          gwas$rsid<-as.character(qtl$rsid)
          gwas$pval<-as.numeric(gwas$p_value)
          
          qtl$pval<-as.numeric(qtl$V12)
        } else {
          gwas<-gwas[paste(gwas$Chr,gwas$Pos,sep="_") %in% paste(substr(qtl$V9,4,6),qtl$V10,sep="_"),]
          qtl<-qtl[paste(substr(qtl$V9,4,6),qtl$V10,sep="_") %in% paste(gwas$Chr,gwas$Pos,sep="_"),]
          
          qtl<-qtl[!duplicated(qtl$V10),]
          gwas<-gwas[!duplicated(gwas$Pos),]
          
          gwas$rsid<-as.character(qtl$rsid)
          gwas$pval<-as.numeric(gwas$pval)
          
          qtl$pval<-as.numeric(qtl$V12)
        }
        
        print (paste0(trait,"_",z$condition_name,"_",y,".pdf"))
        p<-locuscompare2(in_fn1 = gwas, in_fn2 = qtl,title1 = paste(trait,"GWAS"), 
                         title2 = paste(z$condition_name, "eQTL"),snp=gwas$rsid[which.min(gwas$pval)],genome="hg38",combine = T)
        
        pdf(file=paste(getwd(),trait,y,(paste0(z$condition_name,"_",y,".pdf")),sep="/"))
        print(p)
        dev.off()
      }
    })
  }
})

# R version 3.6.2 (2019-12-12)
# Platform: x86_64-apple-darwin15.6.0 (64-bit)
# Running under: macOS  10.16
# 
# Matrix products: default
# LAPACK: /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRlapack.dylib
# 
# Random number generation:
#   RNG:     Mersenne-Twister 
# Normal:  Inversion 
# Sample:  Rounding 
# 
# locale:
#   [1] en_GB.UTF-8/en_GB.UTF-8/en_GB.UTF-8/C/en_GB.UTF-8/en_GB.UTF-8
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] dplyr_1.0.7         stringr_1.4.0       locuscomparer_1.0.0
# 
# loaded via a namespace (and not attached):
#   [1] rstudioapi_0.11  magrittr_1.5     cowplot_1.1.0    tidyselect_1.1.0 munsell_0.5.0    colorspace_1.4-1 R6_2.4.1         rlang_0.4.10     fansi_0.4.1      tools_3.6.2      grid_3.6.2       gtable_0.3.0     utf8_1.1.4      
# [14] DBI_1.1.0        ellipsis_0.3.2   assertthat_0.2.1 tibble_3.0.5     lifecycle_1.0.1  crayon_1.3.4     purrr_0.3.4      ggplot2_3.3.3    vctrs_0.3.8      glue_1.4.1       stringi_1.4.6    compiler_3.6.2   pillar_1.6.4    
# [27] generics_0.0.2   scales_1.1.0     pkgconfig_2.0.3 

