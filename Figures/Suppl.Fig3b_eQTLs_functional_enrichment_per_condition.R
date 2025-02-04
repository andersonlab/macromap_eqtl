library(ggplot2)
library(dplyr)
library(data.table)
library(RColorBrewer)
library(ggrepel)
library(grid)

path_to_read_files<-"/nfs/users/nfs_n/np12/myscratch/MacroMap/Data/MacroMap_manuscript_figures/data/"
path_to_plot<-"/nfs/users/nfs_n/np12/myscratch/MacroMap/Data/MacroMap_manuscript_figures/plots/"


filenames <- list.files(paste0(path_to_read_files,"/eqtl_functional_enrichment/"),
                        pattern="variant_level.bin", full.names=TRUE,recursive=TRUE)



read_data<-function (input) { 
  x=readBin(file(input,"rb"),"double",10000)
  x=matrix(x,(sqrt(length(x)*4+1)-1)/2)
  
  res=cbind(x[,1],sqrt(diag(solve(x[,-1]))))
  res[1,]=res[1,]/5
  
  res=cbind(res[,1],res[,1]-1.96*res[,2],res[,1]+1.96*res[,2],pchisq((res[,1]/res[,2])^2,1,lower=F))
  return(data.frame(res))
}


ldf <- lapply(filenames, read_data)
names(ldf) <- as.character(sapply(filenames,function (x) sapply(strsplit(x,"/"),"[",14)))

annot_name = c("TSS Prox (100kb)","CTCF_binding_site","enhancer","open_chromatin","promoter_flanking_region","promoter","TF_binding_site",
               "3_prime_UTR","5_prime_UTR","frameshift","intron","missense","NC_transcript","NC_transcript_exon","splice_acceptor","splice_donor"
               ,"splice_region","stop_gained","synonymous")


ldf<-lapply(ldf,function (x) { 
  rownames(x)<-annot_name
  colnames(x)=c("log_OR","CI_lower","CI_upper","P_value")
  return(x) }) 



data<-rbindlist(ldf, use.names=TRUE, fill=TRUE, idcol=TRUE)
data$Annotation<-rep(rownames(ldf[[1]]),24)
colnames(data)[1]<-"condition"

# Numbers for the paper  page 5
mean(p.adjust(data$P_value[data$Annotation=="promoter"],method="bonferroni"))
mean(p.adjust(data$P_value[data$Annotation=="TSS Prox (100kb)"],method="bonferroni"))
mean(data$log_OR[data$Annotation=="promoter"])


getPalette = colorRampPalette(brewer.pal(8, "Dark2"))(96) #  
getPalette<-getPalette[c(1,4,7,8,17,20,25,28,33,36,41,44,49,50,57,60,64,66,72,77,84,88,93,96)]

fp2 <-function (x) { ggplot(data=x, aes(x=Annotation, y=log_OR, ymin=CI_lower, ymax=CI_upper,col=condition)) +
    geom_pointrange(position=position_dodge(width=c(0.7)), size=0.2) +
    geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
    coord_flip() +  # flip coordinates (puts labels on y axis)
    xlab("Annotation") + ylab("log OR") + 
    theme_bw(base_size = 13) + scale_color_manual(values=getPalette,guide="none") + facet_grid( Annotation~ log_OR) + facet_wrap(~ condition,nrow=6,ncol=4)
  
}

 
grid_color<-function(gg) {
  g <- ggplot_gtable(ggplot_build(gg))
  strip_both <- which(grepl('strip-t', g$layout$name))
  fills <- getPalette[c(21:24,17:20,13:16,9:12,5:8,1:4)]
  k <- 1
  for (i in strip_both) {
    j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
    g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
    k <- k+1
  }
  return(g)
}

p.c<-grid_color(fp2(data  %>% filter(!Annotation %in% c("splice_acceptor","splice_donor","TF_binding_site","open_chromatin","NC_transcript_exon","intron"))))

pdf(paste0(path_to_plot,"Suppl.Fig3b_eQTL_functional_enrichment_per_condition.pdf"),width = 12,height = 10)
grid.draw(p.c)
dev.off()



