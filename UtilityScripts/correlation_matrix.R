library(RColorBrewer)
library(colorspace)
library(gplots)
library(tidyverse)
library("preprocessCore")


path_to_read_files<-"/nfs/users/nfs_n/np12/myscratch/MacroMap/Data/MacroMap_manuscript_figures/data/"
path_to_plot<-"/nfs/users/nfs_n/np12/myscratch/MacroMap/Data/MacroMap_manuscript_figures/plots/"


expression<-readRDS(paste0(path_to_read_files,"RDS_TPMs_filtd_macromap_fds.rds"))
covariates<-readRDS(paste0(path_to_read_files,"covariates_macromap_fds.rds"))

tpms_list <- split(data.frame(t(expression[-c(1:6)])), f=covariates$Stimulus_Hours)
tpms_list<-lapply(tpms_list, function(x) t(x) )
m<-lapply(tpms_list, function(x) rowMeans(x))
mm<-rbind.data.frame(m)
res<-cor(mm)


#write.table(res,file="correlation_matrix.txt",row.names = T,col.names=T,quote = F)

col_heat<-diverging_hcl(100, palette = "Blue-red")

pdf(paste0(path_to_plot,"correlation_matrix_gene_expression.pdf"))
heatmap.2(res,col=col_heat, trace="none")
dev.off()
