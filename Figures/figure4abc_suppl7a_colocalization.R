# Run this script with R4 
library(data.table)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(tidyverse)
library(cowplot)
library(gplots)

source("/nfs/users/nfs_n/np12/bin/Rscripts/ggplot_themes.R")
source("/nfs/users/nfs_n/np12/myscratch/MacroMap/Data/MacroMap_manuscript_figures/scripts/gg.gap2.R")
source("/nfs/users/nfs_n/np12/myscratch/MacroMap/Data/MacroMap_manuscript_figures/scripts/gencode_annotation.R")
################################################################################################################################################################################################################################
# load coloc data and traits 

path_to_read_files<-"/nfs/users/nfs_n/np12/myscratch/MacroMap/Data/MacroMap_manuscript_figures/data/"
path_to_plot<-"/nfs/users/nfs_n/np12/myscratch/MacroMap/Data/MacroMap_manuscript_figures/plots/"

conditions<-c("CIL_24","CIL_6","Ctrl_24","Ctrl_6","IFNB_24","IFNB_6","IFNG_24","IFNG_6","IL4_24","IL4_6","LIL10_24","LIL10_6","MBP_24","MBP_6","P3C_24","P3C_6","PIC_24",
              "PIC_6","Prec_D0","Prec_D2","R848_24","R848_6","sLPS_24","sLPS_6")

# load coloc data for macromap and GTEx 
macromap<-fread(paste0(path_to_read_files,"colocalization/MacroMap.coloc.gz")) 
gtex<-fread(paste0(path_to_read_files,"colocalization/GTEx.coloc.gz"))

colocs<-rbind(macromap,gtex)

# load trait data 
traits<-read.table(paste0(path_to_read_files,"gwas_files_all.txt"),h=F,sep="\t")
traits<-traits[order(as.character(traits$V1)),]

# Correct spelling mistake 
traits$V7<-as.character(traits$V7)
traits$V7[which(traits$V7=="Autoimmune or Inflamatory disease")]<-"Autoimmune or Inflammatory disease"
traits$V7<-as.factor(traits$V7)


# add study trait tss and filter on gwas p-val and pp4
colocs$trait<-traits$V7[match(colocs$gwas_trait,traits$V1)]
colocs$phenotype_id<-substr(colocs$phenotype_id,1,15)
colocs<-colocs[order(colocs$phenotype_id),]
colocs$study<-ifelse(colocs$condition_name %in% conditions,"MacroMap","GTEx")

tss<-read.table(paste0(path_to_read_files,"gencode.v27.annotation_chr_pos_strand_geneid_gene_name.gtf"),h=F)
colnames(tss)<-c("Chromosome_phe","START_GENE","END_GENE","STRAND","Phenotype_ID","HGNC")

tss$TSS<-ifelse(tss$STRAND=="+",tss$START_GENE-1,tss$END_GENE-1)

colocs<-cbind(colocs,tss[-c(5)][match(colocs$phenotype_id,substr(tss$Phenotype_ID,1,15)),])


## remove them VITDL COFF  POS   UPD They have less than 10 significant regions  
remove_traits<-c("VITDL","COFF","POS","UPD")

colocs<-colocs[!colocs$gwas_trait %in% remove_traits,]
tr<-colocs  %>% distinct(gwas_trait) %>% pull() # 83 unique studies 

# 1955 genes 
length(colocs %>%  filter(study=="MacroMap" & PP.H4.abf >=0.75 & gwas_pval <=5e-08) %>% distinct(phenotype_id) %>% pull())
################################################################################################################################################################
# keep this section in case is needed 

#c<-colocs %>%  filter(study=="MacroMap" &  PP.H4.abf >0.75 & gwas_pval <=5e-08) 
#coloc_genes_in_ctrl<-c %>% filter(condition_name %in% c("Ctrl_6","Ctrl_24")) %>% distinct(phenotype_id) %>% pull()
#coloc_genes_not_in_ctrl<-c %>% filter(!phenotype_id %in% coloc_genes_in_ctrl) %>%  distinct(phenotype_id) %>% pull()  
#length(coloc_genes_not_in_ctrl) # 1324 

#d<-colocs %>%  filter(study=="MacroMap" & gwas_pval <=5e-08) 

#ctrl<-d %>% filter (phenotype_id %in% coloc_genes_not_in_ctrl) %>% filter (condition_name %in% c("Ctrl_6","Ctrl_24"))
#ctrl<-ctrl[!is.na(ctrl$PP.H4.abf),]

#length(unique(ctrl$phenotype_id[ctrl$PP.H4.abf >0.5 ])) # 303 
# 1324-303/1955 --> 52.22506 % 

#Relative to the naive condition, stimulation often increased the strength of evidence of colocalization with disease. 
#For example, ~68% of colocalized eGenes (n=1324) have PP4 >0.75 in the stimulated conditions but not in the naive cells(PP4<0.75). 
#This rate decreased to 52.2% depending on the choice of posterior probability cutoff for determining colocalization in naive cells (PP4 >0.75 in stimulated and PP4 <0.5 in naive). 
#The increased number of eGenes with colocalization evidence is likely to result from both the inclusion of more disease relevant stimulated conditions and from additional power due to biological replication.

################################################################################################################################################################
## Figure4A 
c<-colocs %>%  filter(study=="MacroMap" &  PP.H4.abf >0.75 & gwas_pval <=5e-08)  %>% group_by(condition_name) %>%  summarise(cluster = list(unique(phenotype_id)))

list_signif<-c$cluster
names(list_signif)<-c$condition_name
ctrl_main<- unique(c(c$cluster[[3]],c$cluster[[4]]))
ctrl_test<-ctrl_main

condi<-NULL
gene_le<-NULL

for (i in 1:22) {
  cond<-list_signif[names(list_signif) %in%  names(which.max(lapply(list_signif, function (x) length(which(!x %in% ctrl_test)))))]
  condi[[i]]<-names(cond)
  ctrl_test<-unique(c(ctrl_test,cond[[1]]))
  gene_le[[i]]<-length(ctrl_test)
}

gene_le<-append(append(gene_le,length(c$cluster[[3]]),after=F),length(ctrl_main),after=T)
condi<-append(append(condi,"Ctrl_6",after=F),"Ctrl_24",after=F)
ord<-match(condi,conditions)

getPalette = colorRampPalette(brewer.pal(8, "Dark2"))(96) #
getPalette<-getPalette[c(1,4,7,8,17,20,25,28,33,36,41,44,49,50,57,60,64,66,72,77,84,88,93,96)]

df<-data.frame(colocs=gene_le,Condition=condi,color=getPalette[ord])
df$Condition<- factor(df$Condition, levels = df$Condition)


pdf(paste0(path_to_plot,"Fig4.colocs_added_per_condition.pdf"),useDingbats=FALSE)
ggplot(df,aes(x=Condition, y = colocs,color =Condition,group=1))+ geom_point(size=4.5) + ylab("Cumulative number of\ncolocalized eGenes") + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 2000)) + geom_line() + scale_color_manual(values=getPalette[ord]) + 
  theme_Publication(base_size=20) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  geom_segment(x = 2,y = -10,xend = 2,yend = 632,col = "red",size=1)  + geom_segment(x = 0,y = 632,xend = 2,yend = 632,col = "red",size=1) + guides(color = FALSE) 
dev.off()

################################################################################################################################################################
# load response qtl data, create list of genes and list of egene-esnp pairs 

response_QTLs<-readRDS(file=paste0(path_to_read_files,"response_per_conditions.RDS"))
conditions.order<-conditions[order(conditions)]
names(response_QTLs)<-conditions.order[-c(3,4)]


response_QTLs_genes<-lapply(response_QTLs, function (x) { 
  substr(names(x[which(x)]),1,15)
})

# create list of gene-snp-condition for both reQTLs and shared QTLs 
list_of_list_response_QTLs_genes_snps<-lapply(response_QTLs, function (x) {
  x<-x[which(x)]
  test<-names(x)
  all_d<-str_split_fixed(test,"_",5)
  all_d[,1]<-substr(all_d[,1],1,15)
  egene_esnp<-c(apply(data.frame(all_d[,c(1,2)][all_d[,5]=="",]),1, paste,collapse="_"),apply(data.frame(all_d[,c(1,2,3)][!all_d[,5]=="",]),1, paste,collapse="_"))
  return(egene_esnp)
})

list_of_list_response_QTLs_genes_snps<-setNames(unlist(list_of_list_response_QTLs_genes_snps, use.names=F),rep(names(list_of_list_response_QTLs_genes_snps), lengths(list_of_list_response_QTLs_genes_snps)))
list_response_QTLs_genes_snps<-unique(paste0(as.character(list_of_list_response_QTLs_genes_snps),".",names(list_of_list_response_QTLs_genes_snps)))

list_of_list_shared_QTLs_genes_snps<-lapply(response_QTLs, function (x) {
  x<-x[which(!x)]
  test<-names(x)
  all_d<-str_split_fixed(test,"_",5)
  all_d[,1]<-substr(all_d[,1],1,15)
  egene_esnp<-c(apply(data.frame(all_d[,c(1,2)][all_d[,5]=="",]),1, paste,collapse="_"),apply(data.frame(all_d[,c(1,2,3)][!all_d[,5]=="",]),1, paste,collapse="_"))
  return(egene_esnp)
})

list_of_list_shared_QTLs_genes_snps<-setNames(unlist(list_of_list_shared_QTLs_genes_snps, use.names=F),rep(names(list_of_list_shared_QTLs_genes_snps), lengths(list_of_list_shared_QTLs_genes_snps)))
list_shared_QTLs_genes_snps<-unique(paste0(as.character(list_of_list_shared_QTLs_genes_snps),".",names(list_of_list_shared_QTLs_genes_snps)))

gene_condition_response<-paste0(substr(list_response_QTLs_genes_snps,1,15),".",str_split_fixed(list_response_QTLs_genes_snps,"[.]",2)[,2])
gene_condition_shared<-paste0(substr(list_shared_QTLs_genes_snps,1,15),".",str_split_fixed(list_shared_QTLs_genes_snps,"[.]",2)[,2])


# 617 to 424 check code below 
# Number of reQTL genes that coloc in macromap and condition in which we see the coloc is the same as the one in which we observe the reQTL
macro_response<-colocs %>% filter(study=="MacroMap" &  PP.H4.abf >0.75 & gwas_pval <=5e-08) %>%
  mutate("gene_cond"=paste0(phenotype_id,".",condition_name)) %>% filter(gene_cond %in% gene_condition_response) %>% distinct(phenotype_id)
# 424 reQTL genes 
#######################################################################################################################################################################################################
# Compare abs(betas) of these 424 reQTL genes 

response_eqtl_cond_coloc<-colocs %>% filter(study=="MacroMap" &  PP.H4.abf >0.75 & gwas_pval <=5e-08) %>% 
  mutate("gene_snp_cond"=paste0(paste0(phenotype_id,"_",qtl_lead),".",condition_name))%>% 
  filter(gene_snp_cond %in% list_response_QTLs_genes_snps) 
response_eqtl_cond_coloc_uniq<-unique(response_eqtl_cond_coloc$gene_snp_cond)

sign_eQTLs_betas<-readRDS(paste0(path_to_read_files,"sign_eQTLs_betas.RDS"))
sign_eQTLs_betas<-data.frame(sign_eQTLs_betas)
rownames(sign_eQTLs_betas)<-read.table(paste0(path_to_read_files,"rownames_sign_eQTLs_betas.txt"),h=F)$V1
egene_esnp_response<-sapply(strsplit(rownames(sign_eQTLs_betas),"_"), function (x) paste0(str_split_fixed(x[1],"\\.",2)[1],"_",x[2]))

list_response_eqtl_cond_coloc_uniq<-NULL

for (i in 1:length(response_eqtl_cond_coloc_uniq)) {
  egene_esnp_response_coloc<-sapply(strsplit(response_eqtl_cond_coloc_uniq[i],"[.]"),"[[",1)
  egene_esnp_response_coloc_conditions<-sapply(strsplit(response_eqtl_cond_coloc_uniq[i],"[.]"),"[[",2)
  
  df<-sign_eQTLs_betas[which(egene_esnp_response_coloc==egene_esnp_response),]
  col<-which(egene_esnp_response_coloc_conditions == conditions)
  list_response_eqtl_cond_coloc_uniq[i]<-ifelse(abs(df[,col]) > abs(df[,3]), egene_esnp_response_coloc,NA)
  
}

length(unique(substr(list_response_eqtl_cond_coloc_uniq[!is.na(list_response_eqtl_cond_coloc_uniq)],1,15)))
# 379 egenes  100*(379/424) 89.38679

#####################################################################################################################################################################
## Function to calculate colocs per trait based on the regions 
# another idea to plot pp4 
prep_data_response<-function (x)
{
  trait_regions<-read.table(paste0(path_to_read_files,"colocalization/gwas_regions/",x,".final_regions.bed"))
  names_traits<-trait_regions %>% unite(region, c(V1,V2,V3), sep = "_", remove = FALSE) %>% distinct(region) %>% pull() 
  all<-dim(trait_regions)[1]
  
  macro_response<-colocs %>% filter(gwas_trait==x) %>% filter(study=="MacroMap" &  PP.H4.abf >0.75 & gwas_pval <=5e-08) %>%
    mutate("gene_snp_cond"=paste0(paste0(phenotype_id,"_",qtl_lead),".",condition_name)) %>% filter(gene_snp_cond %in% list_response_QTLs_genes_snps) 
  
   
  macro_not_response<-colocs %>% filter(gwas_trait==x)  %>% filter(study=="MacroMap" &  PP.H4.abf >0.75 & gwas_pval <=5e-08) %>%
    mutate("gene_snp_cond"=paste0(paste0(phenotype_id,"_",qtl_lead),".",condition_name)) %>% filter(gene_snp_cond %in% list_shared_QTLs_genes_snps) 
  
  
  t.response<-macro_response[which(macro_response$phenotype_id %in% macro_not_response$phenotype_id)]
  t.shared<-macro_not_response[which(macro_not_response$phenotype_id %in% macro_response$phenotype_id)]
  
  
  
  if (dim(t.response)[1] !=0) {
    for (i in 1:length(unique(t.response$phenotype_id))) {
      gene<-unique(t.response$phenotype_id)[i]
      t.response.pp4<- t.response %>% filter (phenotype_id==gene) %>% pull(PP.H4.abf) %>% max()
      t.shared.pp4<- t.shared %>% filter (phenotype_id==gene) %>% pull(PP.H4.abf) %>% max()
      
      if( t.response.pp4 >t.shared.pp4  ) {
        macro_not_response<-macro_not_response[!macro_not_response$phenotype_id==gene]
      } else { 
        macro_response<-macro_response[!macro_response$phenotype_id==gene]
      }
    }
  }
  
  t.response_gwas_lead<-macro_response[which(macro_response$gwas_lead %in% macro_not_response$gwas_lead)]
  t.shared_gwas_lead<-macro_not_response[which(macro_not_response$gwas_lead %in% macro_response$gwas_lead)]
  
  if (dim(t.response_gwas_lead)[1] !=0) {
    for (i in 1:length(unique(t.response_gwas_lead$gwas_lead))) {
      snp<-unique(t.response_gwas_lead$gwas_lead)[i]
      t.response.pp4<- t.response_gwas_lead %>% filter (gwas_lead==snp) %>% pull(PP.H4.abf) %>% max()
      t.shared.pp4<- t.shared_gwas_lead %>% filter (gwas_lead==snp) %>% pull(PP.H4.abf) %>% max()
      
      if( t.response.pp4 >t.shared.pp4  ) {
        macro_not_response<-macro_not_response[!macro_not_response$gwas_lead==snp]
      } else { 
        macro_response<-macro_response[!macro_response$gwas_lead==snp]
      }
    }
  }
  
  macro_response_pp4<-macro_response %>% pull(PP.H4.abf) %>% max()
  macro_not_response_pp4<-macro_not_response %>% pull(PP.H4.abf) %>% max()
  
  
  
  return_snps<-function(x,study) {
    Chr<-x[1]
    pos1<-x[2]
    pos2<-x[3]
    res<- study %>% filter(substr(Chromosome_phe,4,6)==Chr & gwas_lead_pos>=pos1 & gwas_lead_pos<=pos2)  %>% distinct(gwas_lead) %>% pull() 
    return(res)
  }
  
  macro_response_snps<-apply(trait_regions,1,return_snps,macro_response)
  macro_not_response_snps<-apply(trait_regions,1,return_snps,macro_not_response)
  
  
  
  f_names<-function (x) {
    tryCatch(
      expr =  {
        names(x)<-names_traits
        return(x)
      },
      error= function (e) {
        return(NULL)
      }
      
    )
  }
  macro_response_snps<-f_names(macro_response_snps)
  macro_not_response_snps<-f_names(macro_not_response_snps)
  
  # if the regions apper in both macro_response_snps_regions and  macro_not_response_snps_regions 
  macro_response_snps_regions<- switch(is.null(macro_response_snps) +1, unique(reshape2::melt(macro_response_snps)$L1),NULL)
  macro_not_response_snps_regions<- switch(is.null(macro_not_response_snps) +1,unique(reshape2::melt(macro_not_response_snps)$L1),NULL)
  
  # if the regions apper in both macro_response_snps_regions and  macro_not_response_snps_regions remove them from macro_response_snps_regions
  macro_response_snps_regions<-macro_response_snps_regions[!macro_response_snps_regions %in% macro_not_response_snps_regions]
  
  not_coloc<-(all- length(unique(c(macro_response_snps_regions,macro_not_response_snps_regions))) )/all
  only_macro_response<-   length(macro_response_snps_regions)/all
  only_macro_not_response<-   length(macro_not_response_snps_regions)/all
  
  
  not_coloc_raw_n<-(all- length(unique(c(macro_response_snps_regions,macro_not_response_snps_regions))) )
  only_macro_response_raw_n<-   length(macro_response_snps_regions)
  only_macro_not_response_raw_n<-   length(macro_not_response_snps_regions)
  
  
  tmp.mat<-matrix(c(only_macro_response,only_macro_not_response,not_coloc))
  tmp.mat_raw<-matrix(c(only_macro_response_raw_n,only_macro_not_response_raw_n,not_coloc_raw_n))
  Trait<-c(rep(x,3))
  Study<-rep(c("Response_eQTLs","Non-response_eQTLs","Not colocalized"))
  coloc.stats<-tmp.mat*100
  coloc.stats_raw<-tmp.mat_raw
  #macro_response_pp4
  
  data<-data.frame(Trait,Study,Coloc=coloc.stats,Coloc_raw=coloc.stats_raw,mean_pp4_res=macro_response_pp4,mean_pp4_not_res=macro_not_response_pp4,total_loci=all)
  data$Study<-factor(data$Study,levels=c("Response_eQTLs","Non-response_eQTLs","Not colocalized"))
  return (data)
}


data_coloc_response<-NULL
for (j in 1:length(tr)){
  print(j)
  data_coloc_response[[j]]<-prep_data_response(tr[j])
}

data_coloc_response.bc<-data_coloc_response
data_coloc_response <- do.call("rbind",data_coloc_response)
data_coloc_response$Type<-traits$V7[match(data_coloc_response$Trait,traits$V1)]

data_coloc_response <- arrange(data_coloc_response, Study, desc(Coloc))

data_coloc_response$Trait<-factor(data_coloc_response$Trait,levels=unique(as.character(data_coloc_response$Trait)))

f<-data_coloc_response %>% group_by(Study) %>% summarize(M=sum(Coloc_raw),M1=sum(total_loci-Coloc_raw)) %>% filter(!Study=="Not colocalized")

mean_shared<-data_coloc_response %>% group_by(Type,Study) %>% summarize(M=mean(Coloc,na.rm=T)) %>% filter(Study=="Non-response_eQTLs")
mean_response<-data_coloc_response %>% group_by(Type,Study) %>% summarize(M=mean(Coloc,na.rm=T)) %>% filter(Study=="Response_eQTLs")

p<-ggplot(data_coloc_response) + geom_bar(aes(fill=Type, y=Coloc, x=Trait,alpha=Study),stat = "identity",position = position_stack(reverse = TRUE))+ 
  facet_grid(~Type, scales = "free",space="free_x") +  scale_alpha_discrete(range = c(1,0.2)) + scale_fill_brewer(palette = "Dark2",direction=-1) +
  geom_hline(data=mean_shared,aes(yintercept=M), size=1, color="grey70",linetype = "solid") + geom_hline(data=mean_response,aes(yintercept=M), size=1, color="grey35",linetype = "solid") + 
  theme(axis.text.y  = element_text(angle=90), axis.text.x = element_text(angle = 90,vjust =0.45),strip.text.x = element_blank(),legend.position="right", legend.direction="horizontal",
        legend.box = "vertical",panel.grid = element_blank(),panel.background = element_rect(fill = "white"),text = element_text(size = 10),legend.title = element_text(angle = 90),
        legend.text = element_text(angle = 90),axis.title.x = element_text(angle = 90,vjust = 0.5),legend.key.size= unit(0.5, "cm"),legend.spacing = unit(0.5, "cm"),
  ) +
  guides(fill = guide_legend(title.position = "left",label.position = "top",title.vjust = 1, label.hjust = 1)
         ,alpha = guide_legend(title.position = "left",label.position = "top"))


p1<-ggplot(data_coloc_response) + geom_bar(aes(fill=Type, y=Coloc_raw, x=Trait,alpha=Study),stat = "identity",position = position_stack(reverse = TRUE))+ 
  facet_grid(~Type, scales = "free",space="free_x") +  scale_alpha_discrete(range = c(1,0.2)) + scale_fill_brewer(palette = "Dark2",direction=-1) +
  geom_hline(data=mean_shared,aes(yintercept=M), size=1, color="grey70",linetype = "solid") + geom_hline(data=mean_response,aes(yintercept=M), size=1, color="grey35",linetype = "solid") + 
  theme(axis.text.y  = element_text(angle=90), axis.text.x = element_text(angle = 90,vjust =0.45),strip.text.x = element_blank(),legend.position="right", legend.direction="horizontal",
        legend.box = "vertical",panel.grid = element_blank(),panel.background = element_rect(fill = "white"),text = element_text(size = 10),legend.title = element_text(angle = 90),
        legend.text = element_text(angle = 90),axis.title.x = element_text(angle = 90,vjust = 0.5),legend.key.size= unit(0.5, "cm"),legend.spacing = unit(0.5, "cm"),
  ) +
  guides(fill = guide_legend(title.position = "left",label.position = "top",title.vjust = 1, label.hjust = 1)
         ,alpha = guide_legend(title.position = "left",label.position = "top"))

pdf(paste0(path_to_plot,"Fig4c_support_barplot_shared_response.pdf"),width =12,useDingbats=FALSE)
print(p)
dev.off()

pdf(paste0(path_to_plot,"Fig4c_support_barplot_shared_response_raw.pdf"),width =12,useDingbats=FALSE)
print(p1)
dev.off()


p1<-gg.gap::gg.gap(plot=p,segments=c(30,80),ylim=c(0,100),rel_heights=c(0.9,0,0.15))
p2<-gg.gap2(plot=p,segments=c(50,90),ylim=c(0,100),rel_heights=c(0.9,0,0.15),vjust=1)

pdf(paste0(path_to_plot,"Fig4c_barplot_shared_response_gg.gap.pdf"),width = 12,useDingbats=FALSE)
print(p2)
dev.off()
#####################################################################################################################################################################
toPlot<-data_coloc_response %>% filter (Study!="Not colocalized") %>% group_by(Trait)  %>% summarise(Coloc = sum(Coloc_raw))
toPlot$Type<-data_coloc_response$Type[1:83]
toPlot$Total<-data_coloc_response %>% filter (Study!="Not colocalized") %>% group_by(Trait) %>% distinct(total_loci) %>% pull(total_loci)

sum(toPlot$Coloc)
#2173 

p<-ggplot(toPlot,aes(x=Coloc,y=Total,col=Type)) + geom_point(size = 3)  + scale_color_brewer(palette = "Dark2",direction=-1) + 
  labs(x="Number of GWAS loci\n colocalized in MacroMap",y="Total GWAS loci") + geom_smooth(method=lm,aes(group=1)) + guides(col=guide_legend(nrow=3,byrow=F))

pdf(paste0(path_to_plot,"coloc_support_macromap_colocs_vs_total.pdf"))
p + theme_Publication(16) + theme(legend.position = "bottom",legend.box="vertical",plot.margin = margin(1,3,1,1, "cm"),legend.key.size= unit(0.3, "cm"))
dev.off()

data_coloc_box<-data_coloc_response %>% filter(Study!="Not colocalized") %>% group_by(Trait,Type) %>% summarize(Coloc=sum(Coloc)) 

p<-ggplot(data_coloc_box,aes(x=Type,y=Coloc,fill=Type)) + geom_boxplot() + coord_cartesian(ylim = c(0, 100)) + scale_fill_brewer(palette = "Dark2",direction=-1) + theme(axis.text.x = element_blank()) +
  scale_x_discrete(name="GWAS studies") + scale_y_continuous(name="Coloc %")

pdf(paste0(path_to_plot,"coloc_support_boxplot_all_eQTLs.pdf"))
p  + theme_Publication(16) + theme(axis.text.x = element_blank(),legend.position = "bottom",legend.box="vertical",plot.margin = margin(1,3,1,1, "cm"),legend.key.size= unit(0.3, "cm")
) + guides(fill=guide_legend(nrow=3,byrow=F))
dev.off()


data_coloc_box %>% group_by(Type) %>% summarise(M=mean(Coloc))
# 1 "Autoimmune or Inflammatory disease"  27.0
# 2 "Blood related diseases and traits "  32.2
# 3 "Cancer"                              25.9
# 4 "Heart related diseases and traits"   22.5
# 5 "Neuro related diseases"              21.2
# 6 "Other"                               23.4

mean( data_coloc_box %>% group_by(Type) %>% summarise(M=mean(Coloc)) %>% pull(M))
# 25.37102
range( data_coloc_box %>% group_by(Type) %>% summarise(M=mean(Coloc)) %>% pull(M))
#21.18837 32.15741


data_coloc_box %>% group_by(Type) %>% summarise(M=range(Coloc))
# 1 "Autoimmune or Inflammatory disease"  9.76
# 2 "Autoimmune or Inflammatory disease" 47.8 
# 3 "Blood related diseases and traits " 23.8 
# 4 "Blood related diseases and traits " 40   
# 5 "Cancer"                             19.3 
# 6 "Cancer"                             36.4 
# 7 "Heart related diseases and traits"   8   
# 8 "Heart related diseases and traits"  31.6 
# 9 "Neuro related diseases"              9.09
# 10 "Neuro related diseases"             41.4 
# 11 "Other"                              14.4 
# 12 "Other"                              40   


data_coloc_box<-data_coloc_response %>% filter(Study!="Not colocalized")

p<-ggplot(data_coloc_box,aes(x=Type,y=Coloc,fill=Type)) + geom_boxplot() + facet_wrap(~Study, scale="free") + 
  scale_fill_brewer(palette = "Dark2",direction=-1)  +scale_x_discrete(name="GWAS studies") + scale_y_continuous(name="Coloc %",limits = c(0, 100), breaks = seq(0, 100, by = 10))

pdf(paste0(path_to_plot,"Fig4b_boxplot_shared_response_eQTLs.pdf"),useDingbats=FALSE)
p + theme_Publication(16) + theme(axis.text.x = element_blank(),legend.position = "bottom",legend.box="vertical",plot.margin = margin(1,3,1,1, "cm"),
                                  legend.key.size= unit(0.3, "cm"),) + guides(fill=guide_legend(nrow=3,byrow=F))
dev.off()

################################################################################################################################################################################################
# calculate reQTL in macromap not in GTEx 

df<-colocs
order_gtex<-df %>% filter(study=="GTEx") %>% distinct(condition_name) %>% arrange(condition_name) %>% pull()
order_macro<-df %>% filter(study=="MacroMap") %>% distinct(condition_name) %>% arrange(condition_name) %>% pull()
df$condition_name <- factor(df$condition_name , levels=c(order_macro,order_gtex))


df_ds_final_075_05<-NULL
df_ds_final_075_05_only_conditions_pp4<-NULL
for (i in  1:length(unique(colocs$gwas_trait))) {
  print(i)
  ds<-unique(colocs$gwas_trait)[i]
  ds_genes<-colocs %>% filter(gwas_trait ==ds)
  macro_genes_ds<-ds_genes  %>% filter(study=="MacroMap" & PP.H4.abf >=0.75 & gwas_pval <=5e-08) %>% distinct(phenotype_id) %>% pull() 
  gtex_genes_ds<-ds_genes %>% filter(study=="GTEx") %>% filter(PP.H4.abf > 0.5 & gwas_pval <=5e-08) %>% distinct(phenotype_id) %>% pull()
  genes_ds_only_macromap<-macro_genes_ds[!macro_genes_ds %in% gtex_genes_ds]
  df_ds.genes<-colocs %>% filter(phenotype_id %in% genes_ds_only_macromap & gwas_trait ==ds & gwas_pval <=5e-08)
  df_ds.genes_2<-colocs %>% filter(phenotype_id %in% genes_ds_only_macromap & gwas_trait ==ds & gwas_pval <=5e-08 & PP.H4.abf >=0.75)

  df_ds_final_075_05[[i]]<-df_ds.genes
  df_ds_final_075_05_only_conditions_pp4[[i]]<-df_ds.genes_2
}

saveRDS(df_ds_final_075_05,paste0(path_to_read_files,"df_ds_final_075macromap_05GTEx.RDS"))


#1955
length(df %>%  filter(study=="MacroMap" & PP.H4.abf >=0.75 & gwas_pval <=5e-08) %>% distinct(phenotype_id) %>% pull() )
#998
length(unique(unlist(lapply(df_ds_final_075_05,function (x) unique(x$phenotype_id)))))

# 290 
#all_macro_respone_no_gtex<-unique(unlist(lapply(df_ds_final_075_05,function (x) unique(x$phenotype_id))))[unique(unlist(lapply(df_ds_final_075_05,function (x) unique(x$phenotype_id)))) %in% unique(as.character(unlist(response_QTLs_genes))) ]
#write.table(all_macro_respone_no_gtex,"all_macro_respone_no_gtex.txt")

all_macro<-df %>%  filter(study=="MacroMap" & PP.H4.abf >=0.75 & gwas_pval <=5e-08) %>% distinct(phenotype_id) %>% pull() 
####

#164  
all_macromap_condition_pp4_unique<-unique(unlist(lapply(df_ds_final_075_05_only_conditions_pp4, function (x) x %>% mutate("gene_cond"=paste0(phenotype_id,".",condition_name)) %>% distinct(gene_cond) %>% pull())))
length(unique(substr(all_macromap_condition_pp4_unique[all_macromap_condition_pp4_unique %in% gene_condition_response],1,15)))
all_macromap_condition_pp4_unique_no_gtex<-unique(substr(all_macromap_condition_pp4_unique[all_macromap_condition_pp4_unique %in% gene_condition_response],1,15))
#write.table(unique(substr(all_macromap_condition_pp4_unique[all_macromap_condition_pp4_unique %in% gene_condition_response],1,15)),"all_macro_respone_pp4_no_gtex.txt")

###########################################################################################################################################################################################################################
## Supplementary figure 7a
#df_ds_final_075_05<-readRDS(paste0(path_to_read_files,"df_ds_final_075macromap_05GTEx.RDS"))
df_to_plot<-do.call(data.frame,rbindlist(df_ds_final_075_05))
order_gtex<-df_to_plot %>% filter(study=="GTEx") %>% distinct(condition_name) %>% arrange(condition_name) %>% pull() %>% as.character()
order_macro<-df_to_plot %>% filter(study=="MacroMap") %>% distinct(condition_name) %>% arrange(condition_name) %>% pull() %>% as.character()
df_to_plot$condition_name <- factor(df_to_plot$condition_name , levels=c(order_macro,order_gtex))


plot_df<-function(df,name,path_to_plot,flag=0) { 
  traits=unique(df$condition_name)
  colocs1=data.frame(df[match(c(outer(sort(unique(df$condition_name)),unique(df$phenotype_id),paste)),paste(df$condition_name,df$phenotype_id)),])
  m=t(matrix(colocs1$PP.H4.abf ,length(traits)))
  colnames(m)=sort(unique(df$condition_name))
  #rownames(m)<-ref2(unique(df$phenotype_id))
  
  imax=function(x){seq(length(x))[x==max(x)][1]}
  Cor=function(x,y){x[is.na(x)]=0;cor(x,y)}
  m=m[order(apply(Cor(t(m[,1:24]),diag(24)),1,imax)),]
  
  zlim=c(-1,1)        
  X<-m
  X[X>zlim[2]]=zlim[2]
  X[X<zlim[1]]=zlim[1]
  N=ncol(X)
  M=nrow(X)
  r<-colorRampPalette(brewer.pal(8, "PRGn"))(100)
  
  pdf(paste0(path_to_plot,name,".pdf"),width  = 16,useDingbats=FALSE)
  heatmap.2(X, trace = "none", dendrogram = "none", Rowv = FALSE, Colv = FALSE, col = r[51:100], breaks = seq(0, 1, length.out = 51),
            key = TRUE,key.title ="PP4", keysize = 1 ,cexCol = 1.2, margins = c(20, 3), labRow = FALSE, labCol = colnames(m))
  
  dev.off()
}

plot_df(df_to_plot,"Supplemenatry_7a_pp4_per_tissue_image_all_macro",path_to_plot)



## plot CD80 Supplementary figure 10 This will create multiple plots 
## check the code in ~/myscratch/MacroMap/Analysis/Coloc_analysis/Macromap_fds/coloc_v2/Analysis/scripts/metaplots_masr_nominal_pp4_all_for_coloc.R you need to change the paths 

cd80<-"ENSG00000121594.11_rs9877891_chr3_119542019"

# create the metaplots 
plot_dis<-paste0("mkdir -p ",paste0(path_to_plot,"plots_075_05_response_metaplots/AST/CD80"))
system(plot_dis)
path_gene<-paste0(path_to_plot,"plots_075_05_response_metaplots/AST/CD80")
setwd(path_gene)
command_plot<-paste0("/software/R-4.1.0/bin/Rscript ",paste0("~/myscratch/MacroMap/Analysis/Coloc_analysis/Macromap_fds/coloc_v2/Analysis/scripts/metaplots_masr_nominal_pp4_all_for_coloc.R ",cd80," ","AST"))
system(command_plot)
setwd("../../")












