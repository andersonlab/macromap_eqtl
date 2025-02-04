library(stringr)
library(data.table)
library(qvalue)
library("colorspace")
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
##########################################################################################################################################################################################
source("/nfs/users/nfs_n/np12/myscratch/MacroMap/Data/MacroMap_manuscript_figures/scripts/ggplot_themes.R")

path_to_read_files<-"/nfs/users/nfs_n/np12/myscratch/MacroMap/Data/MacroMap_manuscript_figures/data/"
path_to_plot<-"/nfs/users/nfs_n/np12/myscratch/MacroMap/Data/MacroMap_manuscript_figures/plots/"

response_QTLs<-readRDS(file=paste0(path_to_read_files,"response_per_conditions.RDS"))

power_per_condition<-read.table(paste0(path_to_read_files,"power.conditions.txt"),h=F)

names(response_QTLs)<-power_per_condition$V1[-c(3,4)]

##########################################################################################################################################################################################
# prepare command to submit as jobs to extract nominal data from GTEx. Need to run only once 
files<-list.files(path = "~/myscratch/MacroMap/Data/GTEx_v8/", pattern = "*.tsv.gz", all.files = FALSE,full.names = TRUE)
is.odd <- function(x) x %% 2 != 0
files<-files[is.odd(1:98)]
tissues<-str_split_fixed(str_split_fixed(files,"//",2)[,2],"[.]",2)[,1]

list_of_list_commands<-lapply(response_QTLs, function (x) {
  
  test<-names(x)
  all_d<-str_split_fixed(test,"_",5)
  
  eQTLs<-data.frame(rbind(all_d[,c(1,2,3)][!all_d[,5]=="",],all_d[,c(1,3,4)][all_d[,5]=="",]))
  colnames(eQTLs)<-c("gene","chr","pos")
  eQTLs$chr<-substr(eQTLs$chr,4,6)
  commands<-NULL
  for (i in 1:length(files)) {
    commands[[i]]<-paste0("tabix ",files[i]," ",eQTLs$chr,":",eQTLs$pos,"-",eQTLs$pos,"|grep ",substr(eQTLs$gene,1,15),">>nominal_",tissues[i],".txt")
  }
  names(commands)<-tissues
  comm<-commands
  return(comm)
})

# # better to save the files and run the commands from bash and then load results to this once
# for(i in 1:length(list_of_list_commands)){
#   x<-list_of_list_commands[i]
#   condition<-str_split_fixed(str_split_fixed(names(x),"_sig",2)[,1],"n_",2)[,2]
#   system(paste0("mkdir ",condition))
#   for (j in 1:length(x[[1]])) {
#     tissue<-x[[1]][j]
#     write.table(tissue,file=paste0(condition,"/commands_",condition,"_",names(tissue),".txt"),col.names=F,row.names=F,quote = F)
#   }
# }

## DO the above once and them submit them as jobs .


## After jobs are done do : 

nonimal_per_condition<-NULL
conditions<-NULL

for(i in 1:length(list_of_list_commands)){
  x<-list_of_list_commands[i]
  condition<-names(x)
  conditions[i]<-condition

files_nom<-list.files(path = paste0(path_to_read_files,"GTEx_replication/",condition), pattern = "nominal*", all.files = FALSE,full.names = TRUE)
files_nom_gtex<-lapply(files_nom, read.table)
names(files_nom_gtex)<-names(x[[1]])
nonimal_per_condition[[i]]<-files_nom_gtex
}

names(nonimal_per_condition)<-conditions

pi1<-lapply(nonimal_per_condition, function (x) unlist(lapply(x, function (z) 1-qvalue(z[3])$pi0)))
p<-data.frame(pi1)

p1<-p[order(apply(p,1,median)),]
p2<-p1[,order(apply(p1,2,median))]

col_heat_red<-diverging_hcl(50, palette = "Blue-Red 3")
#col_heat<-diverging_hcl(100, palette = "Broc")
col_heat<-diverging_hcl(100, palette = "Vik")
#col_heat<-diverging_hcl(100, palette = "Purple-Brown" )

########################################################################################################################################################################################################
# plot response vs non response replication Supplementary figures 6a.b

response_QTLs_genes<-lapply(response_QTLs, function (x) { 
  substr(names(x[which(x)]),1,15)
})

no_response_QTLs_genes<-lapply(response_QTLs, function (x) { 
  substr(names(x[which(!x)]),1,15)
})

nonimal_per_condition_response<-NULL

for (i in 1:length(nonimal_per_condition)) {
  x<-nonimal_per_condition[[i]]
  y<-lapply(x, function (z) {
    z[z[[4]] %in% response_QTLs_genes[[i]],]
  }
         )
  nonimal_per_condition_response[[i]]<-y
}
names(nonimal_per_condition_response)<-conditions

nonimal_per_condition_no_response<-NULL

for (i in 1:length(nonimal_per_condition)) {
  x<-nonimal_per_condition[[i]]
  y<-lapply(x, function (z) {
    z[z[[4]] %in% no_response_QTLs_genes[[i]],]
  }
  )
  nonimal_per_condition_no_response[[i]]<-y
}
names(nonimal_per_condition_no_response)<-conditions


pi1_response<-lapply(nonimal_per_condition_response, function (x) unlist(lapply(x, function (z) 1-tryCatch(qvalue(z[3])$pi0, error=function(e) 0))))

pi1_response<-data.frame(pi1_response)
#pi1_response<-pi1_response[,-c(3,4)]

pi1_response<-pi1_response[order(apply(pi1_response,1,mean),decreasing = T),]
pi1_response<-pi1_response[,order(apply(pi1_response,2,mean),decreasing = T)]


pdf(paste0(path_to_plot,"Suppl.Fig_6b_pi1_response_heatmap.pdf"),width = 10,height = 10,useDingbats = FALSE)
pheatmap(pi1_response,cluster_cols = FALSE,cluster_rows = FALSE,col_heat_red,fontsize= 12, angle_col ="315",breaks= seq(0,1,0.02))
dev.off()


pi1_no_response<-lapply(nonimal_per_condition_no_response, function (x) unlist(lapply(x, function (z) 1-tryCatch(qvalue(z[3])$pi0, error=function(e) 0))))

pi1_no_response<-data.frame(pi1_no_response)

pi1_no_response<-pi1_no_response[match(rownames(pi1_response),rownames(pi1_no_response)),]
pi1_no_response<-pi1_no_response[, match(colnames(pi1_response),colnames(pi1_no_response))]

pdf(paste0(path_to_plot,"Suppl.Fig_6a_pi1_no_response_heatmap.pdf"),width = 10,height = 10,useDingbats = FALSE)
pheatmap(pi1_no_response,cluster_cols = FALSE,cluster_rows = FALSE,col_heat_red,fontsize= 12, angle_col ="315",breaks= seq(0,1,0.02))
dev.off()


########################################################################################################################################################################################################
pi1_no_response_box<-pi1_no_response
pi1_no_response_box$Tissue<-rownames(pi1_no_response)
pi1_no_response_box<-reshape2::melt(pi1_no_response_box,id.vars="Tissue",variable.name='Condition',value.name="Pi1")
pi1_no_response_box$eQTL<-"Non response eQTLs"

pi1_response_box<-pi1_response
pi1_response_box$Tissue<-rownames(pi1_response)
pi1_response_box<-reshape2::melt(pi1_response_box,id.vars="Tissue",variable.name='Condition',value.name="Pi1")
pi1_response_box$eQTL<-"Response eQTLs"

pi_boxplots_data<-rbind(pi1_no_response_box,pi1_response_box)
p<-ggplot(pi_boxplots_data, aes(y=Pi1,x=Condition,fill=eQTL)) + geom_boxplot()

p2<-ggplot(pi_boxplots_data, aes(y=Pi1,x=eQTL,fill=eQTL)) + geom_boxplot()


pdf(paste0(path_to_plot,"Suppl.Fig6c_boxplots_shared_vs_response.pdf"),useDingbats = FALSE)
p + scale_fill_brewer(palette="Accent") + theme_Publication(18)  + theme(axis.text.x=element_text(angle = 90, hjust = 0))
dev.off()

pdf(paste0(path_to_plot,"Suppl.Fig6d_boxplots_shared_vs_response_all.pdf"),useDingbats = FALSE)
p2+ scale_fill_brewer(palette="Accent") + theme_Publication(18)  + theme(axis.text.x=element_text(angle = 0))
dev.off()

w.test<-wilcox.test(pi1_no_response_box$Pi1,pi1_response_box$Pi1)
w.test$p.value

