library("tidyverse")
library("RColorBrewer")
library("preprocessCore")
library("umap")
library("ggrepel")
library("ggplot2")
library("DESeq2")

source("/nfs/users/nfs_n/np12/bin/Rscripts/ggplot_themes.R")

# Load expression data and covariates 
covariates<-readRDS("/nfs/users/nfs_n/np12/myscratch/MacroMap/Analysis/QC/RDS_expression_covariates/Macromap_final_expression/protein_coding_lincRNA/covariates_macromap_fds.rds")
dds_expression<-readRDS("~/myscratch/MacroMap/Analysis/QC/RDS_expression_covariates/Macromap_final_expression/protein_coding_lincRNA/RDS_TPMs_filtd_macromap_fds.rds")

conditions<-c("CIL_24","CIL_6","Ctrl_24","Ctrl_6","IFNB_24","IFNB_6","IFNG_24","IFNG_6","IL4_24","IL4_6","LIL10_24","LIL10_6","P3C_24","P3C_6","R848_24",
              "R848_6","sLPS_24","sLPS_6","PIC_24","PIC_6","MBP_24","MBP_6","Prec_D0","Prec_D2")

conditions.order<-conditions[order(conditions)]


getPalette = colorRampPalette(brewer.pal(8, "Dark2"))(96) #  
getPalette<-getPalette[c(1,4,7,8,17,20,25,28,33,36,41,44,49,50,57,60,64,66,72,77,84,88,93,96)]

# log10 scale expression 
dds_expression.log<-log10(dds_expression[,-c(1:6)]+1)
# Sanity check  
plot(match(colnames(dds_expression.log),covariates$SampleID_RunID))
table(is.na(match(colnames(dds_expression.log),covariates$SampleID_RunID)))
table(match(colnames(dds_expression.log),covariates$SampleID_RunID) -1:length(covariates$SampleID_RunID))

##################################################################################################################################################################################################
# Correct data for all known covariates 

#This is to avoid NA in the model since I keep this sample from another run. Impute Differentiation_time_No_Days
covariates$Differentiation_time_No_Days[covariates$Line=="eipl"][which(is.na(covariates$Differentiation_time_No_Days[covariates$Line=="eipl"]))]<-30 


#Quantile normalizing log TPMs 
phen <- normalize.quantiles(dds_expression.log %>% as.matrix(),copy=TRUE) %>% as.data.frame()
rownames(phen) <- rownames(dds_expression.log)
colnames(phen) <- colnames(dds_expression.log)
print('Quantile normalized data..')

#Rank-based inverse normal transformation
rbint <- function(a,k,s){
  v <- (rank(a,tie="first")-k) /( s+1-(2*k))
  return(v)
}
s <- nrow(phen);
k <- 3/8

phen <- lapply( phen,  rbint, k=k,s=s) %>% lapply(qnorm) %>% as.data.frame()
rownames(phen) <- rownames(dds_expression.log)
colnames(phen) <- colnames(dds_expression.log)
print('Rank-based inverse normal transformed (rb-int) data..')


# regress out covariates 
res.expression <-t(residuals(lm (t(phen) ~ as.factor(as.character(covariates$RunID)) + 
                                   as.factor(as.character(covariates$Donor)) + as.factor(as.character(covariates$Library_prep)) + as.factor (covariates$Sex) + as.numeric(covariates$Purity_result_per) 
                                 + as.factor(covariates$Differentiation_media) + as.numeric(covariates$Purity_result_per) + as.numeric(covariates$Estimated_cell_diameter) 
                                 + as.factor(covariates$Differentiation_time_No_Days))))


# change umap std variables 
custom.config = umap.defaults
custom.config$n_neighbors <-100 # 15 defaults 
custom.config$n_epochs<-400 # large separates samples better # 200 defaults 

# umap with all conditions# submit this as job and save model 
umap_1<-umap(t(res.expression),custom.config)
umap_model = as.data.frame(umap_1$layout) 
umap_model$Stimulus_Hours<-covariates$Stimulus_Hours
umap_model$Lib_type<-covariates$Library_prep
umap_model$Donor<-covariates$Donor

saveRDS(umap_model,file="/nfs/users/nfs_n/np12/myscratch/MacroMap/Data/MacroMap_manuscript_figures/data/umap_model.RDS")

# you need to reverse coor *-1 to match previous version plots 
umap_model<-readRDS("/nfs/users/nfs_n/np12/myscratch/MacroMap/Data/MacroMap_manuscript_figures/data/umap_model.RDS")
umap_model$V1 <-umap_model$V1 *-1
umap_model$V2 <-umap_model$V2 *-1


umap_model$Stimulus<-sapply(strsplit(as.character(umap_model$Stimulus_Hours),"_"),"[[",1)


##########################################################################################################################################################################################v
# uncorrected data 

# change umap std variables 
custom.config = umap.defaults
custom.config$n_neighbors <-100
custom.config$n_epochs<-400 # large separates samples better 


umap_2<-umap(t(phen),custom.config)
umap_model2 = as.data.frame(umap_2$layout) 
umap_model2$Stimulus_Hours<-covariates$Stimulus_Hours
umap_model2$Lib_type<-covariates$Library_prep
umap_model2$Donor<-covariates$Donor


##########################################################################################################################################################################################
umap_plot_all<-function (obj_umap,toplot) {
  i=1
  j=2
  
  plot<- ggplot(data = obj_umap, aes_string(x = paste("V",i,sep=""), y = paste("V",j,sep=""), color = toplot)) + theme_bw()+ 
    geom_point(size = 1,shape=1)  + xlab(paste0("uMAP ",i)) + ylab(paste0("uMAP ",j)) + coord_fixed() + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank() ,
          axis.text=element_text(size=16),axis.title=element_text(size=20,face="bold"),legend.text=element_text(size=14),legend.title=element_text(size=15))
  return (plot)
}

##########################################################################################################################################################################################
# plot uncorrected data 

p<-umap_plot_all(umap_model2,"Lib_type")
pdf("/nfs/users/nfs_n/np12/myscratch/MacroMap/Data/MacroMap_manuscript_figures/plots/Suppl_Fig1e_umap_all_conditions_lib_type_uncorrected_data.pdf",useDingbats=FALSE)
print(p + theme(aspect.ratio=1))
dev.off() 
##
##########################################################################################################################################################################################
p<-umap_plot_all(umap_model,"Stimulus_Hours") +
  geom_label_repel(data= umap_model %>% slice(1:24),aes(label=Stimulus_Hours),show.legend = FALSE,box.padding   = 0.5,point.padding = 0.5,size=3) 

p<-p  + theme_Publication(20) 
p<- p + guides(color = guide_legend(override.aes = list(size = 4,shape=20),ncol=1,byrow=FALSE,title.position = "top")) + scale_color_manual(values=getPalette) + theme(legend.position="right")


pdf("/nfs/users/nfs_n/np12/myscratch/MacroMap/Data/MacroMap_manuscript_figures/plots/Fig1_support_umap_all_conditions_Stimulus_Hours.pdf",useDingbats=FALSE)
print(p)
dev.off()

getPalette = colorRampPalette(brewer.pal(8, "Dark2"))(96) #  
getPalette<-getPalette[c(1,4,7,8,17,20,25,28,33,36,41,44,49,50,57,60,64,66,72,77,84,88,93,96)]
get_odd<-getPalette[1:24%%2 == 1]

p1<-umap_plot_all(umap_model,"Stimulus") +  geom_label_repel(data= umap_model %>% slice(1:24),aes(label=Stimulus_Hours),show.legend = FALSE,box.padding   = 0.5,point.padding = 0.5,size=3) 
p1<-p1  + theme_Publication(20)
p1<- p1 + guides(color = guide_legend(override.aes = list(size = 4,shape=20),ncol=1,byrow=FALSE,title.position = "top")) + 
  scale_color_manual(values=get_odd) + theme(legend.position="right")

pdf("/nfs/users/nfs_n/np12/myscratch/MacroMap/Data/MacroMap_manuscript_figures/plots/Fig1a.umap_all_conditions_Stimulus.pdf",useDingbats=FALSE)
print(p1)
dev.off()

##Time point 

umap_model$Timepoint<-sapply(strsplit(as.character(umap_model$Stimulus_Hours),"_"),"[[",2)
umap_model$Timepoint[umap_model$Timepoint=="D0"]<-"Prec_D0"
umap_model$Timepoint[umap_model$Timepoint=="D2"]<-"Prec_D2"
umap_model$Timepoint<-factor(umap_model$Timepoint,levels=c("6","24","Prec_D0","Prec_D2"))
umap_model<-umap_model[rev(order(umap_model$Timepoint)),]

p1<-umap_plot_all(umap_model,"Timepoint") #+  stat_ellipse(aes(color=Time))
p1<-p1  + theme_Publication(20)
p1<- p1 + guides(color = guide_legend(override.aes = list(size = 4,shape=20),ncol=1,byrow=FALSE,title.position = "top")) + 
  scale_color_manual(values=brewer.pal(8, "Set1")) + theme(legend.position="right")


pdf("/nfs/users/nfs_n/np12/myscratch/MacroMap/Data/MacroMap_manuscript_figures/plots/Fig1b.umap_all_conditions_Stimulus_Hours_time.pdf",useDingbats=FALSE)
print(p1)
dev.off()


# abstract test 
fisher.test(matrix(c(424,1531,1954,6261),nrow=2,ncol=2))$p.val

# 
# colorRampPalette(brewer.pal(8, "Dark2"))
# umap_plot_stim_alpha<-function (obj_umap,toplot,sti) {
#   i=1
#   j=2
#   
#   plot<- ggplot(data = obj_umap, aes_string(x = paste("V",i,sep=""), y = paste("V",j,sep=""), color = toplot)) + theme_bw()+
#     geom_point(size = 1,alpha=0.05) + geom_point(size=2.5,data=umap_model[umap_model$Stimulus_Hours %in% sti,], aes_string(x = paste("V",i,sep=""), y = paste("V",j,sep=""))) + xlab(paste0("uMAP_",i)) + ylab(paste0("uMAP_",j)) + coord_fixed() +
#     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank() ,
#           axis.text=element_text(size=16),axis.title=element_text(size=20,face="bold"),legend.text=element_text(size=14),legend.title=element_text(size=15),
#           plot.title=element_text(color=getPalette[which(conditions.order %in% sti)],hjust = 0.5, size=24, face="bold")) + ggtitle(sti)  
#   return (plot)
# }
# # The function can plot more stimuli than 1 
# 
# for ( i in 1:length(conditions.order)) {
#   
#   p<-umap_plot_stim_alpha(umap_model,"Stimulus_Hours",conditions.order[i])
#   pdf(paste0("/nfs/users/nfs_n/np12/myscratch/MacroMap/Analysis/UMAP/Macromap_fds/plots.log10/umap_all_conditions_Stimulus_Hours_",conditions.order[i],".pdf"))
#   print(p + theme(aspect.ratio=1) + scale_color_manual(values=getPalette))
#   dev.off()
# } 
##################################################################################################################################################################################################
