# locuscompare2 function in order to create plots with same ylim 

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
  p = make_combined_plot2(merged, title1, title2, ld, chr, snp, 
                         combine, legend, legend_position, lz_ylab_linebreak)
  return(p)
}


add_label<-function (merged, snp) 
{
  merged$label = ifelse(merged$rsid %in% snp, merged$rsid, "")
  return(merged)
}


make_locuszoom2<-function (metal, title, chr, color, shape, size,max_y, ylab_linebreak = FALSE) 
{
  p = ggplot(metal, aes(x = pos, logp)) + geom_point(aes(fill = rsid,size = rsid, shape = rsid), alpha = 0.8) + 
    geom_point(data = metal[metal$label !=  "", ], aes(x = pos, logp, fill = rsid, size = rsid, shape = rsid)) + 
    scale_fill_manual(values = color, guide = "none") + scale_shape_manual(values = shape, guide = "none") + scale_size_manual(values = size, guide = "none") + 
    coord_cartesian(ylim = c(0, max_y))+
    scale_x_continuous(labels = function(x) {
      sprintf("%.1f", x/1e+06)
    }) + ggrepel::geom_text_repel(aes(label = label)) + xlab(paste0("chr", 
                                                                    chr, " (Mb)")) + ylab(bquote(.(title) ~ -log[10] * "(P)")) + 
    theme_classic() + theme(plot.margin = unit(c(0.5, 1, 
                                                 0.5, 0.5), "lines"))
  if (ylab_linebreak == TRUE) {
    p = p + ylab(bquote(atop(.(title), -log[10] * "(P)")))
  }
  return(p)
}


make_combined_plot2<- function (merged, title1, title2, ld, chr, snp = NULL, combine = TRUE, 
          legend = TRUE, legend_position = c("bottomright", "topright", 
                                             "topleft"), lz_ylab_linebreak = FALSE) 
{
  snp = get_lead_snp(merged, snp)
  color = assign_color(merged$rsid, snp, ld)
  shape = ifelse(merged$rsid == snp, 23, 21)
  names(shape) = merged$rsid
  size = ifelse(merged$rsid == snp, 3, 2)
  names(size) = merged$rsid
  merged = add_label(merged, snp)
  
  p1 = make_scatterplot(merged, title1, title2, color, shape, 
                        size, legend, legend_position)
  metal1 = merged[, c("rsid", "logp1", "chr", "pos", "label")]
  colnames(metal1)[which(colnames(metal1) == "logp1")] = "logp"
  metal1<-metal1[!is.infinite(metal1$logp),]
  
  metal2 = merged[, c("rsid", "logp2", "chr", "pos", "label")]
  colnames(metal2)[which(colnames(metal2) == "logp2")] = "logp"
  metal2<-metal2[!is.infinite(metal2$logp),]
  
  max_y<-max(c(max(metal1$logp),max(metal2$logp)))
  
  
  p2 = make_locuszoom2(metal1, title1, chr, color, shape, size,max_y, 
                      lz_ylab_linebreak)
  p3 = make_locuszoom2(metal2, title2, chr, color, shape, size,max_y,
                      lz_ylab_linebreak)
  if (combine) {
    p2 = p2 + theme(axis.text.x = element_blank(), axis.title.x = element_blank())
    p4 = cowplot::plot_grid(p2, p3, align = "v", nrow = 2, 
                            rel_heights = c(0.8, 1))
    p5 = cowplot::plot_grid(p1, p4)
    return(p5)
  }
  else {
    return(list(locuscompare = p1, locuszoom1 = p2, locuszoom2 = p3))
  }
}

