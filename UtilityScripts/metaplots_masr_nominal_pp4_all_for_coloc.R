#Rscript metaplot_by_summary_stats.R eQTL.txt 
# Run it with R4 

suppressMessages(library(rmeta))
suppressMessages(library(RColorBrewer))
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(ggrepel))
suppressMessages(library("cowplot"))
suppressMessages(library(mashr))

################################################################################################################################################################################

path<-"/nfs/users/nfs_n/np12/myscratch/MacroMap/Analysis/eQTLs/Macromap_fds/analysis/mashR/analysis/analysisV2/8_flashR_best_eQTLs_random200K_commonbaseline/"

m.ed_common_baseline<-readRDS(paste0(path,"m.ed.RDS"))
Vhat<-readRDS(paste0(path,"Vhat.RDS"))


path<-"/nfs/users/nfs_n/np12/myscratch/MacroMap/Analysis/eQTLs/Macromap_fds/analysis/mashR/files/8_flashR_best_eQTLs_random200K/"
m.ed_all<-readRDS(paste0(path,"m.ed.RDS"))
Vhat_all <-readRDS("/nfs/users/nfs_n/np12/myscratch/MacroMap/Analysis/eQTLs/Macromap_fds/analysis/mashR/files/7_flashR_all_eQTLs_random200K/Vhat.RDS")


####################################################################################################################################################################################
source("/nfs/users/nfs_n/np12/bin/Rscripts/gencode_annotation.R")

args = commandArgs(trailingOnly=TRUE)

eQTL<-args[1]
disease<-args[2]

#eQTL<-"ENSG00000121594.11_rs9877891_chr3_119542019"
#disease<-"SSC"

gene<-sapply(strsplit(eQTL,split="_"),"[[",1)
snp<-sapply(strsplit(eQTL,split="_"),"[[",2)

gene_snp<-paste(gene,snp,sep="_")

getPalette = colorRampPalette(brewer.pal(8, "Dark2"))(96) #  
getPalette<-getPalette[c(1,4,7,8,17,20,25,28,33,36,41,44,49,50,57,60,64,66,72,77,84,88,93,96)]

conditions<-c("CIL_24","CIL_6","Ctrl_24","Ctrl_6","IFNB_24","IFNB_6","IFNG_24","IFNG_6","IL4_24","IL4_6","LIL10_24","LIL10_6","MBP_24","MBP_6","P3C_24","P3C_6","PIC_24",
              "PIC_6","Prec_D0","Prec_D2","R848_24","R848_6","sLPS_24","sLPS_6")

getPalette_no_Ctrl = colorRampPalette(brewer.pal(8, "Dark2"))(96) #  
getPalette_no_Ctrl<-getPalette_no_Ctrl[c(1,4,17,20,25,28,33,36,41,44,49,50,57,60,64,66,72,77,84,88,93,96)]

data<-read.table(text=system(paste('bash /nfs/users/nfs_n/np12/myscratch/MacroMap/Analysis/eQTLs/Macromap_fds/analysis/metaplots_mashR/extract_beta_std.sh',eQTL),intern = TRUE))


betasPlot<-data$V13
stdPlot<-data$V15
meta.ord<-order(betasPlot,decreasing=TRUE)

d<-data.frame(betasPlot,stdPlot,conditions)


####################################################################################################################################################################################
# Fit the 2 mash R models common baseline and all conditions 
bhat<-matrix(d$betasPlot[!d$conditions=="Ctrl_6"],nrow=1)
bhat<-t(matrix(rep(t(bhat),2),ncol=2))

shat<-matrix(d$stdPlot[!d$conditions=="Ctrl_6"],nrow=1)
shat<-t(matrix(rep(t(shat),2),ncol=2))

data_to_test_all_eQTLs = mash_set_data(bhat,shat) # load data that are significant (best snp per gene for every condition) 
data_to_test_all_eQTLs = mash_update_data(data_to_test_all_eQTLs,V=Vhat,ref=3)
m.ed.t_common_baseline<-mash(data_to_test_all_eQTLs,g=m.ed_common_baseline$fitted_g,fixg=TRUE) # fit the data to test 


bhat<-matrix(d$betasPlot,nrow=1)
bhat<-t(matrix(rep(t(bhat),2),ncol=2))

shat<-matrix(d$stdPlot,nrow=1)
shat<-t(matrix(rep(t(shat),2),ncol=2))

data_to_test_all_eQTLs_v2 = mash_set_data(bhat,shat) # load data that are significant (best snp per gene for every condition) 
data_to_test_all_eQTLs_v2 = mash_update_data(data_to_test_all_eQTLs_v2,V=Vhat_all)
m.ed.t_all<-mash(data_to_test_all_eQTLs_v2,g=m.ed_all$fitted_g,fixg=TRUE) # fit the data to test 

####################################################################################################################################################################################
# plot the common baseline model 

meta.ord<-order(m.ed.t_common_baseline$result$PosteriorMean[1,],decreasing=TRUE)
betasPlot<-m.ed.t_common_baseline$result$PosteriorMean[1,]
names(betasPlot)<-conditions[-c(3,4)]

stdPlot<-m.ed.t_common_baseline$result$PosteriorSD[1,]
names(stdPlot)<-conditions[-c(3,4)]

betasPlot<-betasPlot[meta.ord]
stdPlot<-stdPlot[meta.ord]

label<-names(betasPlot)

pdf(paste0(ref(gene),"_",gene_snp,".metaplot_mashR_baseline.pdf"))
metaplot(as.numeric(betasPlot),as.numeric(stdPlot),labels = label,xlab="PosteriorBetas",ylab="Conditions",
         col=meta.colors(box = as.character(getPalette_no_Ctrl[meta.ord])),main= paste(ref(gene),"\n",gene_snp))
dev.off()


betasPlot1  <- as.numeric(betasPlot)
sd <- as.numeric(stdPlot)
col<-as.character(getPalette_no_Ctrl[meta.ord])
conf.level<-0.95
ci.value <- -qnorm((1 - conf.level)/2)
lower.sd <- betasPlot1 - ci.value * sd
upper.sd <- betasPlot1 + ci.value * sd

df <- data.frame(label, betasPlot1, sd,col,lower.sd,upper.sd)



# reverses the factor level ordering for labels after coord_flip()
df$label <- factor(df$label, levels=(df$label))
xl<-extendrange(c(floor(range(df$lower.sd)[1]),ceiling(range(df$upper.sd)[2])),f=0.05)

fp_baseline <- ggplot(data=df, aes(x=label, y=betasPlot1, ymin=lower.sd, ymax=upper.sd,color=label)) +
  geom_pointrange() + scale_color_manual(values=col)  + ylim(xl) + 
  geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("Conditions") + ylab("PosteriorBetas") + theme_bw(base_size = 16) + theme(legend.position = "none") +
  theme(plot.title = element_text(hjust = 0.5,size = 16, face = "bold"))   + ggtitle(paste(ref(gene),"\n",gene_snp)) 

pdf(paste0(ref(gene),"_",gene_snp,".metaplot_mashR_ggplot_baseline_no_lfsr_thres.pdf"))
print(fp_baseline)
invisible(dev.off())


####################################################################################################################################################################################
# plot full mashR model 

meta.ord<-order(m.ed.t_all$result$PosteriorMean[1,],decreasing=TRUE)
betasPlot<-m.ed.t_all$result$PosteriorMean[1,]
names(betasPlot)<-conditions

stdPlot<-m.ed.t_all$result$PosteriorSD[1,]
names(stdPlot)<-conditions

betasPlot<-betasPlot[meta.ord]
stdPlot<-stdPlot[meta.ord]

label<-names(betasPlot)

pdf(paste0(ref(gene),"_",gene_snp,".metaplot_mashR_full_model.pdf"))
metaplot(as.numeric(betasPlot),as.numeric(stdPlot),labels = label,xlab="PosteriorBetas",ylab="Conditions",
         col=meta.colors(box = as.character(getPalette_no_Ctrl[meta.ord])),main= paste(ref(gene),"\n",gene_snp))
dev.off()


betasPlot1  <- as.numeric(betasPlot)
sd <- as.numeric(stdPlot)
col<-as.character(getPalette[meta.ord])
conf.level<-0.95
ci.value <- -qnorm((1 - conf.level)/2)
lower.sd <- betasPlot1 - ci.value * sd
upper.sd <- betasPlot1 + ci.value * sd

df <- data.frame(label, betasPlot1, sd,col,lower.sd,upper.sd)



# reverses the factor level ordering for labels after coord_flip()
df$label <- factor(df$label, levels=(df$label))
xl<-extendrange(c(floor(range(df$lower.sd)[1]),ceiling(range(df$upper.sd)[2])),f=0.05)

fp_full <- ggplot(data=df, aes(x=label, y=betasPlot1, ymin=lower.sd, ymax=upper.sd,color=label)) +
  geom_pointrange() + scale_color_manual(values=col)  + ylim(xl) + 
  geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("Conditions") + ylab("PosteriorBetas") + theme_bw(base_size = 16) + theme(legend.position = "none") +
  theme(plot.title = element_text(hjust = 0.5,size = 16, face = "bold"))   + ggtitle(paste(ref(gene),"\n",gene_snp)) 

pdf(paste0(ref(gene),"_",gene_snp,".metaplot_mashR_ggplot_full_model_no_lfsr_thres.pdf"))
print(fp_full)
invisible(dev.off())

##################################################################################################################################################
betasPlot<-data$V13
stdPlot<-data$V15
meta.ord<-order(betasPlot,decreasing=TRUE)

d<-data.frame(betasPlot,stdPlot,conditions)
d<-d[meta.ord,]

pdf(paste0(ref(gene),"_",eQTL,".metaplot_condition_by_condition_nominal.pdf"))
metaplot(as.numeric(d$betasPlot),as.numeric(d$stdPlot),labels = d$conditions,
         col=meta.colors(box = as.character(getPalette[meta.ord])),xlab="Nominal beta",ylab="Conditions",main= paste(ref(gene),"\n",gene_snp))
invisible(dev.off())

gene_snp<-sapply(strsplit(eQTL,"_"), function(x) paste(x[1:2], collapse = '_'))
cat(as.character(gene_snp))

label <- d$conditions
betasPlot1  <- as.numeric(d$betasPlot)
sd <- as.numeric(d$stdPlot)
col<-as.character(getPalette[meta.ord])
conf.level<-0.95
ci.value <- -qnorm((1 - conf.level)/2)
lower.sd <- betasPlot1 - ci.value * sd
upper.sd <- betasPlot1 + ci.value * sd

df <- data.frame(label, betasPlot1, sd,col,lower.sd,upper.sd)

# reverses the factor level ordering for labels after coord_flip()
df$label <- factor(df$label, levels=(df$label))
xl<-extendrange(c(floor(range(df$lower.sd)[1]),ceiling(range(df$upper.sd)[2])),f=0.05)


fp_nominal <- ggplot(data=df, aes(x=label, y=betasPlot1, ymin=lower.sd, ymax=upper.sd,color=label)) +
  geom_pointrange() + scale_color_manual(values=col)  + ylim(xl) + 
  geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("Conditions") + ylab("Nominal Beta") + theme_bw(base_size = 16) + theme(legend.position = "none") +
  theme(plot.title = element_text(hjust = 0.5,size = 16, face = "bold"))  + ggtitle(paste(ref(gene),"\n",gene_snp)) 

pdf(paste0(ref(gene),"_",eQTL,".metaplot_condition_by_condition_nominal_ggplot.pdf"))
print(fp_nominal)
invisible(dev.off())


################################################################################################################################################################################
#df_ds_final_09_05<-readRDS("/lustre/scratch119/realdata/mdt3/teams/gaffney/np12/MacroMap/Analysis/Coloc_analysis/Macromap_fds/coloc_new_gwas/analysis/df_ds_final_09_05.RDS")
df_ds_final_09_05<-readRDS("/nfs/users/nfs_n/np12/myscratch/MacroMap/Analysis/Coloc_analysis/Macromap_fds/coloc_v2/Analysis/df_ds_final_075_05.RDS")

lapply(df_ds_final_09_05, function (x)unique(data.frame(x)$gwas_trait))
names(df_ds_final_09_05)<-lapply(df_ds_final_09_05, function (x)unique(data.frame(x)$gwas_trait))

coloc_data<-as.data.frame(df_ds_final_09_05[names(df_ds_final_09_05)==disease],col.names="")
PP4<-coloc_data[coloc_data$phenotype_id==substr(gene,1,15),] %>% filter(study=="MacroMap") %>% distinct(PP.H4.abf) %>% pull()

df.betas_condition_by_condition<-data.frame(PP4=as.numeric(PP4[meta.ord]),Nominal_beta=betasPlot1,SD=sd,condition=d$conditions,lower.sd=lower.sd,upper.sd=upper.sd)
xl<-extendrange(c(floor(range(df.betas_condition_by_condition$lower.sd)[1]),ceiling(range(df.betas_condition_by_condition$upper.sd)[2])),f=0.1)
yl<-c(0,1)

p3<-ggplot(df.betas_condition_by_condition, aes(x=Nominal_beta, y=PP4,color =condition)) + geom_point(size=2.5) + geom_errorbarh(aes(xmin=lower.sd, xmax=upper.sd),height=0.01)+ 
  xlim(xl) + ylim(yl) + scale_color_manual(values=getPalette) +geom_label_repel(aes(label=condition),show.legend = FALSE,box.padding   = 0.6 ,point.padding = 0.6,max.overlaps=24) + 
  theme_bw(base_size = 16) + theme(legend.position = "none") + geom_vline(xintercept=0, linetype="dashed", color = "red",size=1.5) + 
  geom_hline(yintercept=0.75, linetype="dashed", color = "blue",size=1.5)

p4<-ggplot(df.betas_condition_by_condition, aes(x=Nominal_beta, y=PP4,color =condition)) + geom_point(size=2.5) + 
  geom_label_repel(aes(label=condition),show.legend = FALSE,box.padding   = 0.6 ,point.padding = 0.6,max.overlaps=24,direction= "both",nudge_x=ifelse(max(df.betas_condition_by_condition$Nominal_beta) >0,-1.5,1.5)) + 
  xlim( xl) + ylim(yl)+ scale_color_manual(values=getPalette) + theme_bw(base_size = 16) + theme(legend.position = "none") +
  geom_vline(xintercept=0, linetype="dashed", color = "red",size=1.5) +
  geom_hline(yintercept=0.75, linetype="dashed", color = "blue",size=1.5)



pdf(paste0(ref(gene),"_",eQTL,"_",disease,"_PP4vsBetas_condition_by_condition_summary_stats.pdf"))
print(p3)
dev.off()

pdf(paste0(ref(gene),"_",eQTL,"_",disease,"_PP4vsBetas_condition_by_condiion_summary_stats_no_SD.pdf"))
print(p4)
dev.off()

final_plot<-plot_grid(fp_nominal, p4, labels = "AUTO",rel_widths = c(1, 1))
pdf(paste0(ref(gene),"_",eQTL,"_",disease,"_metaplot_pp4vsBetas.pdf"),width = 11)
print(final_plot)
dev.off()


final_plot<-plot_grid(fp_baseline, p4, labels = "AUTO",rel_widths = c(1, 1))
pdf(paste0(ref(gene),"_",eQTL,"_",disease,"_metaplot_pp4vsPosteriorBetas_baseline.pdf"),width = 11)
print(final_plot)
dev.off()


final_plot<-plot_grid(fp_full, p4, labels = "AUTO",rel_widths = c(1, 1))
pdf(paste0(ref(gene),"_",eQTL,"_",disease,"_metaplot_pp4vsPosteriorBetas_all_data.pdf"),width = 11)
print(final_plot)
dev.off()


final_plot<-plot_grid(fp_full,fp_baseline, p4, labels = "AUTO",rel_widths = c(1, 1),nrow = 1)
pdf(paste0(ref(gene),"_",eQTL,"_",disease,"_metaplot_pp4vsPosteriorBetas_all_data_baseline.pdf"),width = 15)
print(final_plot)
dev.off()


