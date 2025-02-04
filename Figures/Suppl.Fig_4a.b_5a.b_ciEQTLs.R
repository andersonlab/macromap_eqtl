# Run it with R session 4 
library(dplyr)
library(tidyverse)
library(ggsci) # palettes 
library(ggpubr)
library(gridExtra)
library(grid)
library(mashr)
library(RColorBrewer)

##########################################################################################################################################################################################
source("/nfs/users/nfs_n/np12/myscratch/MacroMap/Data/MacroMap_manuscript_figures/scripts/ggplot_themes.R")
path_to_read_files<-"/nfs/users/nfs_n/np12/myscratch/MacroMap/Data/MacroMap_manuscript_figures/data/"
path_to_plot<-"/nfs/users/nfs_n/np12/myscratch/MacroMap/Data/MacroMap_manuscript_figures/plots/"

# load the data 
ci_eQTLs<-read.table(paste0(path_to_read_files,"all_sentinel_ci_eqtls.txt"),header=TRUE)
ci_eQTLs_cond<-ci_eQTLs %>% group_by(condition,rank,phe_id) %>% filter (rank >0) %>% 
  distinct(pos=paste0(phe_chr,":",phe_from,"-",phe_to),rs=var_id,eqtl=paste0(phe_id,"_",var_id,"_",phe_chr,"_",var_from))
##########################################################################################################################################################################################
# 1714 ci_eQTLs in stimulated 
ci_eQTLs_cond %>% ungroup %>% filter(!condition %in% c("Ctrl_6","Ctrl_24")) %>% distinct(phe_id)

ci_mm_eqtls <- ci_eQTLs %>%
  group_by(condition, phe_id, rank) %>%
  mutate(sen_rsid = paste("sen_", var_id[bwd_best_hit == "1"])) %>%
  ungroup(rank)  %>%
  mutate(cie_per_gn = factor(sum(bwd_best_hit == "1"))) %>%
  ungroup()%>% 
  unite(cond_gene, c(condition, phe_id), remove = F) %>%
  unite(cond_gene_rsid, c(condition, phe_id, var_id), remove = F) %>%
  mutate(across(rank, ~factor(., levels = c("0", "1", "2", "3"), labels = c("Primary eQTL", "Secondary eQTL", "Tertiary eQTL", "Quaternary eQTL"))))%>%
  mutate(pc = case_when(condition %in% c("IFNB_6", "IFNG_24", "IFNG_6", "LIL10_24", "Prec_D0", "Prec_D2") ~ "35",
                        condition %in% c("CIL_24", "CIL_6", "Ctrl_24", "Ctrl_6", "IFNB_24", "IL4_24", "MBP_24", "P3C_24", "PIC_24", "R848_24", "sLPS_24") ~ "40",
                        condition %in% c("IL4_6", "LIL10_6", "MBP_6", "P3C_6", "PIC_6", "R848_6", "sLPS_6") ~ "50"))


toplot<-ci_mm_eqtls %>% 
  ungroup() %>% 
  arrange(desc(cie_per_gn)) %>% 
  distinct(phe_id, `.keep_all` = T) %>% 
  add_count(cie_per_gn)


## Supplementary figures 4a 
p<-ggplot(data=toplot,aes(y =n, x = cie_per_gn)) + 
  geom_bar(position="dodge", stat = "identity",fill = pal_jama()(1)) + 
  geom_text(aes(label=n), size = 5, vjust=-0.5, hjust = 0.5) +
  labs(x = "Number of independent cis-eQTLs per gene", y = "Number of eGenes") + theme_classic(20) + 
  annotate("text", x = 3, y = 2600, label = "18% of eGenes have >1 CI cis-eQTL", size = 5.5, fontface = "bold") 

pdf(paste0(path_to_plot,"Suppl.Fig4a_CI_eQTLs_barplot.pdf"),useDingbats = FALSE)
print(p)
dev.off()
 

ci_mm_eqtls %>% 
  mutate(abs_tss = abs(dist_phe_var))%>% 
  compare_means(abs_tss ~ rank,  data = .)

my_comparisons <- list( c("Primary eQTL", "Secondary eQTL"), c("Primary eQTL", "Tertiary eQTL"), c("Primary eQTL", "Quaternary eQTL") )

p1<-ci_mm_eqtls %>% 
  mutate(abs_tss = abs(dist_phe_var)) %>%
  ggboxplot(., x = "rank", y = "abs_tss",
            fill = "rank", palette = "jco")+
  stat_compare_means(comparisons = my_comparisons) +
  labs(title = "Distance between CI cis-eQTL to TSS by rank", x = "CI cis-eQTLs rank", y = "Absolute distance to TSS") + 
  theme(text = element_text(size = 16),legend.text = element_text(size = 10),legend.position = "top") +
  guides(fill=guide_legend(nrow=2)) + scale_x_discrete(labels = c("Primary\neQTL","Secondary\neQTL","Tertiary\neQTL","Quaternary\neQTL"))
## Supplementary figures 4b 
pdf(paste0(path_to_plot,"Suppl.Fig4b_CI_eQTLs_tss_distance.pdf"),useDingbats = FALSE,width =7 ,height = 7)
print(p1)
dev.off()

pdf(paste0(path_to_plot,"Suppl.Fig4ab_CI_eQTLs_barplot_CI_eQTLs_tss_distance.pdf"),useDingbats = FALSE,width =14 ,height = 7)
grid.arrange(p,p1,ncol=2, widths = c(7, 7))
dev.off()

#########################################################################################################################################################
# Code for suppl Fig 5A 5B 


#########################################################################################################################################################
# Run the below once and save the object to extract info (betas std err) from the nominal pass summary stats to run mashR
# qtl<-NULL
# q<-NULL
# q1<-NULL
# 
# 
# for (i in 1:length(unique(ci_eQTLs$condition))) { 
#   cond<-unique(as.character(ci_eQTLs$condition))[i]
#   ci_eQTLs_cond_per_cond<-ci_eQTLs_cond %>% filter(condition==cond)
#   for (j in 1:dim(ci_eQTLs_cond_per_cond)[1]){
#     qtl<-read.table(text=system(paste0("tabix ~/myscratch/MacroMap/Analysis/eQTLs/Macromap_fds/analysis/eQTLs_per_TSS/1MB/nominal/",
#                                        condition,"/",condition,".nominal.geneTSS.sorted.txt.gz ",ci_eQTLs_cond_per_cond$pos[j]),intern = TRUE))
#     q1<-qtl[which(qtl$V16 %in% ci_eQTLs_cond_per_cond[j,]$eqtl),]
#     q[[cond]]<-rbind(q[[cond]],q1)
#   }
# }
# 
# saveRDS(paste0(path_to_read_files,"cieqtls_summary_stats.RDS"))

qtls<-readRDS(paste0(path_to_read_files,"cieqtls_summary_stats.RDS"))
m.ed_common_baseline<-readRDS(paste0(path_to_read_files,"m.ed.RDS"))
Vhat<-readRDS(paste0(path_to_read_files,"Vhat.RDS"))


secondary_eQTLs<-as.character(unlist(lapply(qtls,function (x) as.character(x$V16))))

# extract betas and std errors across all conditions 

get_beta<- function (eQTL) {
  
  data<-read.table(text=system(paste('bash ~/myscratch/MacroMap/Analysis/eQTLs/Macromap_fds/analysis/conditional_analysis_jesse/extract_beta_std.sh',eQTL),intern = TRUE))
  beta<-data$V13
  return(beta)  
  
}
get_std<- function (eQTL) {
  
  data<-read.table(text=system(paste('bash ~/myscratch/MacroMap/Analysis/eQTLs/Macromap_fds/analysis/conditional_analysis_jesse/extract_beta_std.sh',eQTL),intern = TRUE))
  std<-data$V15
  return(std)  
}

sign_betas<-NULL
sign_std<-NULL


# Do this once and save data 
######################################################################################################################################################################################################
#secondary_eQTLs<-gsub(";","';'",secondary_eQTLs)

# for (i in 1:length(secondary_eQTLs)) {
#   sign_betas[[i]]<-get_beta(secondary_eQTLs[i])
#   sign_std[[i]]<-get_std(secondary_eQTLs[i])
#     
# }

#saveRDS(sign_betas,"sign_betas.RDS")
#saveRDS(sign_std,"sign_std.RDS")
############################################################################################################################################################################################################
sign_betas<-readRDS(paste0(path_to_read_files,"cieqtls_sign_betas.RDS"))
sign_std<-readRDS(paste0(path_to_read_files,"cieqtls_sign_std.RDS"))

sign_betas<-do.call(rbind,sign_betas)
sign_std<-do.call(rbind,sign_std)

rownames(sign_betas)<-secondary_eQTLs
rownames(sign_std)<-secondary_eQTLs

colnames(sign_betas)<-names(qtls)
colnames(sign_std)<-names(qtls)

## There are eQTLs found as condition specific in multiple conditions, just remove the duplicates  no need to calculate mash stats 

sign_betas<-sign_betas[!duplicated(sign_betas),]
sign_std<-sign_std[!duplicated(sign_std),]


bhat<-as.matrix(sign_betas[,!colnames(sign_betas)=="Ctrl_6"])
shat<-as.matrix(sign_std[,!colnames(sign_std)=="Ctrl_6"])


data_to_test_all_eQTLs = mash_set_data(bhat,shat) # load data that are significant (best snp per gene for every condition) 
data_to_test_all_eQTLs = mash_update_data(data_to_test_all_eQTLs,V=Vhat,ref=3)
m.ed.t_common_baseline<-mash(data_to_test_all_eQTLs,g=m.ed_common_baseline$fitted_g,fixg=TRUE) # fit the data to test 

saveRDS(m.ed.t_common_baseline,file=paste0(path_to_read_files,"m.ed.t_baseline_secondary_eQTLs.RDS"))
########################################################################################################################################################################################################
# Analysis  and plots 
ci_eQTLs_cond_per_cond<-qtls
list_condition_by_condition_eQTLs_mashR<-lapply(ci_eQTLs_cond_per_cond,function (x) m.ed.t_common_baseline$result$lfsr[rownames(m.ed.t_common_baseline$result$lfsr) %in% x$V16,] )
list_condition_by_condition_eQTLs_mashR<-list_condition_by_condition_eQTLs_mashR[-c(3,4)]

response_QTLs<-NULL
for (i in 1:22) {
  response_QTLs[[i]]<-list_condition_by_condition_eQTLs_mashR[[i]][,i] <=0.05
} 
#saveRDS(response_QTLs,"response_inde_QTLs.RDS")

response_QTLs_plot<-sapply(response_QTLs, function (x) table(!x)) # this is just for ploting reasons 
colnames(response_QTLs_plot)<-names(ci_eQTLs_cond_per_cond)[-c(3,4)]
#response_QTLs_plot<-response_QTLs_plot[,order((apply(response_QTLs_plot,2,sum)),decreasing = TRUE)]

response_QTLs_plot<-response_QTLs_plot[,order(response_QTLs_plot[1,]/response_QTLs_plot[2,])]


pdf(paste0(path_to_plot,"Suppl.Fig5a_barplots_response_eQTLs_condition_by_condition.pdf"),useDingbats = FALSE)
par(mar = c(5,7, 3,2))
mp<-barplot(response_QTLs_plot,las=1,horiz=TRUE,col= brewer.pal(4, "Accent")[c(2,1)],xaxt="n",
            xlim=c(0,400),cex.names = 1.5,cex.main=1.7,main="Number of Conditionally independent\nresponse eQTLs")
axis(1, at=seq(0,400,100), labels=seq(0, 400, 100),las=1,cex.axis=1.5)
t<-round(apply(response_QTLs_plot,2, function(x) x[1]/(x[1]+x[2])*100),1)
text(response_QTLs_plot[1,],mp - 0.1, labels = paste0(t,"%"," ( ",response_QTLs_plot[1,]," ) "), pos = 4,cex=1.2,col="black",font=2)
dev.off()


unique_ci_egenes<-unique(substr(unique(as.character(unlist(lapply(list_condition_by_condition_eQTLs_mashR,rownames)))),1,15))
unique_ci_egenes_response<-unique(substr(unique(unlist(lapply(response_QTLs, function(x) names(x[x])))),1,15))

length(unique_ci_egenes_response)/length(unique_ci_egenes)
#######################################################################################################################################

response_QTLs_ind<-response_QTLs
response_QTLs<-readRDS(file=paste0(path_to_read_files,"response_per_conditions.RDS"))

conditions<-c("CIL_24","CIL_6","Ctrl_24","Ctrl_6","IFNB_24","IFNB_6","IFNG_24","IFNG_6","IL4_24","IL4_6","LIL10_24","LIL10_6","MBP_24","MBP_6","P3C_24","P3C_6","PIC_24",
              "PIC_6","Prec_D0","Prec_D2","R848_24","R848_6","sLPS_24","sLPS_6")

names(response_QTLs)<-conditions[-c(3,4)]

response_QTLs_genes_primary<-lapply(response_QTLs, function (x) substr(names(which (x)),1,15) )
response_QTLs_genes_secondary<-lapply(response_QTLs_ind, function (x) substr(names(which (x)),1,15) )

no_response_QTLs_genes_primary<-lapply(response_QTLs, function (x) substr(names(which (!x)),1,15) )
no_response_QTLs_genes_secondary<-lapply(response_QTLs_ind, function (x) substr(names(which (!x)),1,15) )


t1.1<-table(unique(unlist(response_QTLs_genes_secondary)) %in% unique(unlist(response_QTLs_genes_primary)))
t1.2<-table(unique(unlist(no_response_QTLs_genes_secondary)) %in% unique(unlist(response_QTLs_genes_primary)))

# I do have douple counts of CI eGenes being re or not re 
table(unique(unlist(no_response_QTLs_genes_secondary)) %in% unique(unlist(response_QTLs_genes_secondary)))

sec_reQTLs_primary_res<-t1.1[2]
sec_reQTLs_primary_nores<-t1.1[1]

no_sec_reQTLs_primary_res<-t1.2[2]
no_sec_reQTLs_primary_nores<-t1.2[1] - 106 # I have to remove the 106 genes that are response in some condition and not re in others from the Primary is not reQTL 

## Add this to text plus figure 
test<-fisher.test(matrix(c(sec_reQTLs_primary_res,no_sec_reQTLs_primary_res,sec_reQTLs_primary_nores,no_sec_reQTLs_primary_nores),nrow=2))
#Fisher's Exact Test for Count Data
##data:  table(df)
#p-value < 2.2e-16
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  7.693928 23.047551
#sample estimates:
#odds ratio 
#  12.92196 

# Create a data frame
df <- data.frame(Group = c("Primary eQTL is reQTL", "Primary eQTL is not reQTL", "Primary eQTL is reQTL", "Primary eQTL is not reQTL"),
                 Outcome = c("CI-eQTLs that are reQTLs", "CI-eQTLs that are reQTLs", "CI-eQTLs that are not reQTLs", "CI-eQTLs that are not reQTLs"),
                 Count = c(142, 17, 614, 951),
                 stringsAsFactors = FALSE)

# Calculate percentages
df$Percent <- df$Count / sum(df$Count) * 100

# Calculate Fisher's exact test
mat <- matrix(c(df[1, 3], df[2, 3], df[3, 3], df[4, 3]), nrow = 2, byrow = TRUE,
              dimnames = list(c("Primary eQTL is reQTL", "Primary eQTL is not reQTL"), 
                              c("CI-eQTLs that are reQTLs", "CI-eQTLs that are not reQTLs")))
res <- fisher.test(mat)

# Extract odds ratio, confidence interval, and p-value
odds_ratio <- res$estimate
conf_int <- res$conf.int
p_val <- res$p.value

# Calculate fold change
fe <- log2(odds_ratio)
#fc <- (142/756) / (17/968)  check this!! 

# Create the plot
p<-ggplot(df, aes(x = Outcome, y = Percent, fill = Group)) + 
  geom_col(position = "dodge") + 
  ylab("Percentage") + 
  scale_fill_manual(values = c("#8da0cb", "#fc8d62"), 
                    labels = c("Primary eQTL is not reQTL", "Primary eQTL is reQTL")) + 
  labs(title = paste("Comparison of CI-eQTLs that are reQTLs and\n CI-eQTLs that are not reQTLs"), 
       subtitle = paste("Fisher's exact test: p-value = ", round(p_val, 36)," OR = ", round(odds_ratio, 2), ", \n95% CI [", 
                        round(conf_int[1], 2), ", ", round(conf_int[2], 2), "], FE = ", round(fe, 2), sep = "")) + 
  theme_minimal() +
  theme(plot.title = element_text(size = 18, face = "bold"),
        plot.subtitle = element_text(size = 14),
        plot.caption = element_text(size = 10),
        legend.title = element_blank(),
        legend.position = "top",
        legend.text = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  guides(fill = guide_legend(reverse = TRUE)) + 
  geom_text(aes(label = paste0(round(Percent),"% (", Count, ")")), position = position_dodge(width = 0.9), vjust = -0.5, size = 4)

pdf("percentages_shared_response_primary_secondary_v2.pdf",useDingbats = FALSE)
pdf(paste0(path_to_plot,"Suppl.Fig5b_percentages_shared_response_primary_secondary_v2.pdf"),useDingbats = FALSE)

print(p)
dev.off()



