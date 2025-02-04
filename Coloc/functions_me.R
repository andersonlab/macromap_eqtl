importVariantInformation <- function(file){
  info_col_names = c("chr","pos","snp_id","major","minor","MAF")
  into_col_types = "cicccn"
  snp_info = readr::read_delim(file, delim = "\t", col_types = into_col_types, col_names = info_col_names)
  snp_info = dplyr::mutate(snp_info, indel_length = pmax(nchar(major), nchar(minor))) %>%
    dplyr::mutate(is_indel = ifelse(indel_length > 1, TRUE, FALSE))
  return(snp_info)
}

constructQtlListForColoc <- function(qtl_root, sample_size_list){
  conditions = names(sample_size_list)
  min_pvalues = list()
  qtl_summary_list = list()

  #Iterate over conditions and fill lists
  for(condition in conditions){
    min_pvalue_path = file.path(qtl_root, "permuted", paste0(condition, ".permuted.txt.gz"))
    summary_path = file.path(qtl_root, "nominal", paste0(condition,"/",condition, ".nominal.sorted.txt.gz"))

    min_pvalues[[condition]] = importQTLtoolsTable(min_pvalue_path) %>% dplyr::select(Phenotype_ID, Best_variant_in_cis, qvalue)
    qtl_summary_list[[condition]] = summary_path
  }
  return(list(min_pvalues = min_pvalues, qtl_summary_list = qtl_summary_list, sample_sizes = sample_size_list))
}

importQTLtoolsTable<-function(file_path){
  col_types = "cciiciicciiidddddddd"
  table = readr::read_delim(file_path, col_names = TRUE, delim = "\t", col_types = col_types) %>%
    dplyr::filter(!is.na(Corrected_pvalue)) %>%
    dplyr::mutate(p_bonferroni = Nominal_pvalue*Beta_dist_2_number_of_ind_tests) %>%
    dplyr::mutate(p_bonferroni = pmin(p_bonferroni,1)) %>%
    dplyr::mutate(p_fdr = p.adjust(Corrected_pvalue, method = "BH")) %>%
    dplyr::mutate(qvalue = qvalue::qvalue(Corrected_pvalue)$qvalues) %>%
    dplyr::arrange(qvalue)
  return(table)
}

importGWASSummary <-function(summary_path){
  gwas_col_names = c("snp_id", "chr", "pos", "effect_allele", "MAF",
                     "p_nominal", "beta", "OR", "log_OR", "se", "z_score", "trait", "PMID", "used_file")
  gwas_col_types = c("ccicdddddddccc")
  gwas_pvals = readr::read_tsv(summary_path,
                               col_names = gwas_col_names, col_types = gwas_col_types, skip = 1)
  return(gwas_pvals)
}

prefilterColocCandidates <-function(qtl_min_pvalues, gwas_prefix, gwas_variant_info,
                                     fdr_thresh = 0.1, overlap_dist = 1e5, gwas_thresh = 1e-5){


  # debug
  #qtl_min_pvalues=phenotype_values$min_pvalues
  #gwas_variant_info = qtl_var_info
  #fdr_thresh = 1
  #overlap_dist = 1e7
  #gwas_thresh = 1e-5

  #Make sure that the qtl_df has neccessary columns
  assertthat::assert_that(assertthat::has_name(qtl_min_pvalues[[1]], "Phenotype_ID"))
  assertthat::assert_that(assertthat::has_name(qtl_min_pvalues[[1]], "Best_variant_in_cis"))
  assertthat::assert_that(assertthat::has_name(qtl_min_pvalues[[1]], "qvalue"))
  assertthat::assert_that(ncol(qtl_min_pvalues[[1]]) == 3)


  #Import top GWAS p-values
  gwas_pvals = importGWASSummary(paste0(gwas_prefix,".top_hits.txt.gz")) %>%
    dplyr::filter(p_nominal < gwas_thresh) %>%
    dplyr::transmute(chr = paste0('chr',chr), gwas_pos = pos)

  #Filter lead variants
  qtl_hits = purrr::map(qtl_min_pvalues, ~dplyr::filter(., qvalue < fdr_thresh))
  lead_variants = purrr::map_df(qtl_hits, identity) %>% unique()
  selected_variants = dplyr::filter(gwas_variant_info, snp_id %in% lead_variants$Best_variant_in_cis) %>%
    dplyr::select(chr, pos, snp_id) %>%
    dplyr::transmute(chr=chr,pos=pos,Best_variant_in_cis = snp_id)

  #Add GRCh37 coordinates (or #38 depends on used file gwas_variant_info )
  qtl_pos = purrr::map(qtl_hits, ~dplyr::left_join(., selected_variants, by = "Best_variant_in_cis") %>%
                         dplyr::filter(!is.na(pos)))

  #Identify genes that have associated variants nearby (ignoring LD)
  qtl_df_list = purrr::map(qtl_pos, ~dplyr::left_join(., gwas_pvals, by = "chr") %>%
                             dplyr::mutate(distance = abs(gwas_pos - pos)) %>%
                             dplyr::filter(distance < overlap_dist) %>%
                             dplyr::select(Phenotype_ID, Best_variant_in_cis) %>% unique())

}

qtltoolsTabixFetchPhenotypes<-function(phenotype_ranges, tabix_file){

  #Assertions about input
  assertthat::assert_that(class(phenotype_ranges) == "GRanges")
  assertthat::assert_that(assertthat::has_name(GenomicRanges::elementMetadata(phenotype_ranges), "phenotype_id"))

  #Set column names for rasqual
  qtltools_columns = c("phenotype_id","pheno_chr","pheno_start", "pheno_end",
                      "strand","n_snps", "distance", "snp_id", "snp_chr",
                      "snp_start", "snp_end", "p_nominal","beta", "is_lead","std.err","pheno_snp_id")
  qtltools_coltypes = "cciiciicciiddidc"

  result = list()
  for (i in seq_along(phenotype_ranges)){
    selected_phenotype_id = phenotype_ranges[i]$phenotype_id
    print(i)
    tabix_table = scanTabixDataFrame(tabix_file, phenotype_ranges[i],
                                     col_names = qtltools_columns, col_types = qtltools_coltypes)[[1]] %>%
      dplyr::filter(phenotype_id == selected_phenotype_id)

    #Add additional columns
    result[[selected_phenotype_id]] = tabix_table
  }
  return(result)
}

constructVariantRanges_gwas<-function(variant_df, variant_information, cis_dist){

  #Make key assertions
  assertthat::assert_that(assertthat::has_name(variant_df, "snp_id"))
  assertthat::assert_that(assertthat::has_name(variant_information, "snp_id"))
  assertthat::assert_that(assertthat::has_name(variant_information, "chr"))
  assertthat::assert_that(assertthat::has_name(variant_information, "pos"))

  #Filter variant information to contain only required snps
  var_info = dplyr::filter(variant_information, snp_id %in% variant_df$snp_id) %>%
    dplyr::select(snp_id, chr, pos, MAF)

  #Add variant info to variant df
  var_df = dplyr::left_join(variant_df, var_info, by = "snp_id") %>%
  dplyr::transmute(phenotype_id=phenotype_id,snp_id=snp_id,chr=substr(chr,4,10),pos=pos,MAF=MAF)

  #Make a ranges object
  var_ranges = var_df %>%
    dplyr::rename(seqnames = chr) %>%
    dplyr::mutate(start = pos - cis_dist, end = pos + cis_dist, strand = "*") %>%
    dataFrameToGRanges()

  return(var_ranges)
}

summaryReplaceSnpId<-function(summary_df, variant_information){

  #Make key assertions
  assertthat::assert_that(assertthat::has_name(summary_df, "snp_id"))
  assertthat::assert_that(assertthat::has_name(summary_df, "pos"))
  assertthat::assert_that(assertthat::has_name(summary_df, "chr"))

  #Filter variant information to contain only required snps
  var_info = dplyr::filter(variant_information, pos %in% summary_df$pos) %>%
    dplyr::select(snp_id, chr, pos, MAF) %>%   dplyr::transmute(snp_id,chr=gsub("chr","",chr),pos,MAF)

  #Remove MAF if it is present
  if(assertthat::has_name(summary_df, "MAF")){
    summary_df = dplyr::select(summary_df, -MAF)
  }

  #Add new coordinates
  new_coords = dplyr::select(summary_df, -snp_id) %>%
    dplyr::left_join(var_info, by = c("chr","pos")) %>%
    dplyr::filter(!is.na(snp_id)) %>%
    dplyr::arrange(pos)

  return(new_coords)
}

colocQtlGWAS_nikos<-function(qtl, gwas, N_qtl,N){

  #Check that QTL df has all correct names
  assertthat::assert_that(assertthat::has_name(qtl, "snp_id"))
  assertthat::assert_that(assertthat::has_name(qtl, "beta"))
  assertthat::assert_that(assertthat::has_name(qtl, "MAF"))
  assertthat::assert_that(assertthat::has_name(qtl, "p_nominal"))

  #Check that GWAS df has all correct names
  assertthat::assert_that(assertthat::has_name(gwas, "beta"))
  assertthat::assert_that(assertthat::has_name(gwas, "se"))
  assertthat::assert_that(assertthat::has_name(gwas, "snp_id"))
  assertthat::assert_that(assertthat::has_name(gwas, "log_OR"))
  assertthat::assert_that(assertthat::has_name(gwas, "MAF"))

  #Count NAs for log_OR and beta
  log_OR_NA_count = length(which(is.na(gwas$log_OR)))
  beta_NA_count = length(which(is.na(gwas$beta)))

  #Remove GWAS SNPs with NA std error
 if (is.null(N))
     {
 	 gwas = dplyr::filter(gwas, !is.na(se))
     }

  #If beta is not specified then use log_OR
  # For some GWAS I need to specify N in order coloc to run with just the pvalues
  if (is.null(N)) {
    if(beta_NA_count <= log_OR_NA_count){
    coloc_res = coloc::coloc.abf(dataset1 = list(pvalues = qtl$p_nominal,
                                                 N = N_qtl,
                                                 MAF = qtl$MAF,
                                                 type = "quant",
                                                 beta = qtl$beta,
                                                 snp = qtl$snp_id,
                                                 s=0.5),
                                 dataset2 = list(beta = gwas$beta,
                                                 varbeta = gwas$se^2,
                                                 type = "cc",
                                                 snp = gwas$snp_id,
                                                 MAF = gwas$MAF,
                                                s=0.5))
    } else {
      coloc_res = coloc::coloc.abf(dataset1 = list(pvalues = qtl$p_nominal,
                                                 N = N_qtl,
                                                 MAF = qtl$MAF,
                                                 type = "quant",
                                                 beta = qtl$beta,
                                                 snp = qtl$snp_id,
                                                 s=0.5),
                                 dataset2 = list(beta = gwas$log_OR,
                                                 varbeta = gwas$se^2,
                                                 type = "cc",
                                                 snp = gwas$snp_id,
                                                 MAF = gwas$MAF,
                                                 s=0.5))
    }
    } else {
      coloc_res = coloc::coloc.abf(dataset1 = list(pvalues = qtl$p_nominal,
                                                   N = N_qtl,
                                                   MAF = qtl$MAF,
                                                   type = "quant",
                                                   beta = qtl$beta,
                                                   snp = qtl$snp_id,
                                                   s=0.5),
                                   dataset2 = list(pvalues = gwas$p_nominal,
                                                   type = "cc",
                                                   snp = gwas$snp_id,
                                                   MAF = gwas$MAF,
                                                   s=0.5,N=N))
    }

  return(coloc_res)
  }


colocMolecularQTLs_nikos <- function(qtl_df, qtl_summary_path, gwas_summary_path,
                               gwas_variant_info, qtl_variant_info,
                               N_qtl = 84, cis_dist = 1e5,N = NULL,QTLTools = TRUE){

  #debug
  #qtl_df<-qtl_pairs[1,]
  #tl_summary_path<-phenotype_values$qtl_summary_list$IFNG_6
  #gwas_summary_path<-paste0(gwas_prefix, ".sorted.txt.gz")
  #gwas_variant_info<-qtl_var_info
  #qtl_variant_info<-qtl_var_info
  #N_qtl =phenotype_values$sample_sizes$IFNG_6
  #cis_dist = cis_window
  #QTLTools = TRUE
  #N=NULL
#print (N)
  #Assertions
  assertthat::assert_that(assertthat::has_name(qtl_df, "phenotype_id"))
  assertthat::assert_that(assertthat::has_name(qtl_df, "snp_id"))
  assertthat::assert_that(nrow(qtl_df) == 1)

  assertthat::assert_that(is.numeric(cis_dist))
  assertthat::assert_that(is.numeric(N_qtl))

  #Print for debugging
  message(qtl_df$phenotype_id)
  #print(qtl_summary_path)

  result = tryCatch({
    #Make GRanges object to fetch data
   qtl_ranges = constructVariantRanges(qtl_df, qtl_variant_info, cis_dist = cis_dist)
   gwas_ranges = constructVariantRanges_gwas(qtl_df, gwas_variant_info, cis_dist = cis_dist)

    #Fetch QTL summary stats
    if(QTLTools){
      qtl_summaries = qtltoolsTabixFetchPhenotypes(qtl_ranges, qtl_summary_path[[1]])[[1]] %>%
        dplyr::transmute(snp_id, chr = snp_chr, pos = snp_start, p_nominal, beta) %>%
        dplyr::mutate(snp_id=gsub("_.*","",snp_id))
    } else{
      qtl_summaries = fastqtlTabixFetchGenes(qtl_ranges, qtl_summary_path)[[1]]
    }

    #Fetch GWAS summary stats
    gwas_summaries = tabixFetchGWASSummary(gwas_ranges, gwas_summary_path)[[1]]

    #Substitute coordinate for the eqtl summary stats and add MAF
    qtl = summaryReplaceCoordinates(qtl_summaries, gwas_variant_info)

    #Substitute snp_id for the GWAS summary stats and add MAF
    gwas = summaryReplaceSnpId(gwas_summaries, gwas_variant_info)

    #Extract minimal p-values for both traits
    qtl_min = dplyr::arrange(qtl, p_nominal) %>% dplyr::filter(row_number() == 1)
    gwas_min = dplyr::arrange(gwas, p_nominal) %>% dplyr::filter(row_number() == 1)


    #Perform coloc analysis
    coloc_res = colocQtlGWAS_nikos(qtl, gwas, N_qtl = N_qtl,N=N)
    coloc_summary = dplyr::as_tibble(t(data.frame(coloc_res$summary))) %>%
      dplyr::mutate(qtl_pval = qtl_min$p_nominal, gwas_pval = gwas_min$p_nominal,
                    qtl_lead = qtl_min$snp_id, gwas_lead = gwas_min$snp_id,chr=qtl_min$chr,
                    gwas_lead_pos = gwas_min$pos,qtl_lead_pos=qtl_min$pos) #Add minimal pvalues

    #Summary list
    data_list = list(qtl = qtl, gwas = gwas)

    result = list(summary = coloc_summary, data = data_list)
  }, error = function(err) {
    print(paste("ERROR:",err))
    result = list(summary = NULL, data = NULL)
  }
  )
  return(result)
}

constructVariantRanges <-function(variant_df, variant_information, cis_dist){

  #Make key assertions
  assertthat::assert_that(assertthat::has_name(variant_df, "snp_id"))
  assertthat::assert_that(assertthat::has_name(variant_information, "snp_id"))
  assertthat::assert_that(assertthat::has_name(variant_information, "chr"))
  assertthat::assert_that(assertthat::has_name(variant_information, "pos"))

  #Filter variant information to contain only required snps
  var_info = dplyr::filter(variant_information, snp_id %in% variant_df$snp_id) %>%
    dplyr::select(snp_id, chr, pos, MAF)

  #Add variant info to variant df
  var_df = dplyr::left_join(variant_df, var_info, by = "snp_id")

  #Make a ranges object
  var_ranges = var_df %>%
    dplyr::rename(seqnames = chr) %>%
    dplyr::mutate(start = pos - cis_dist, end = pos + cis_dist, strand = "*") %>%
    dataFrameToGRanges()

  return(var_ranges)
}

colocMolecularQTLsByRow_nikos <- function(qtl_df,qtl_summary_path,gwas_summary_path,gwas_variant_info,qtl_variant_info,N_qtl,cis_dist,N){
  result = purrrlyr::by_row(qtl_df, ~colocMolecularQTLs_nikos(.,qtl_summary_path,gwas_summary_path,gwas_variant_info,qtl_variant_info,N_qtl,cis_dist,N)$summary, .collate = "rows")

}


dataFrameToGRanges<-function(df){
  #Convert a data.frame into a GRanges object

  gr = GenomicRanges::GRanges(seqnames = df$seqnames,
               ranges = IRanges::IRanges(start = df$start, end = df$end),
               strand = df$strand)

  #Add metadata
  meta = dplyr::select(df, -start, -end, -strand, -seqnames)
  GenomicRanges::elementMetadata(gr) = meta

  return(gr)
}

tabixFetchGWASSummary <- function(granges, summary_path){
  gwas_col_names = c("snp_id", "chr", "pos", "effect_allele", "MAF",
                     "p_nominal", "beta", "OR", "log_OR", "se", "z_score", "trait", "PMID", "used_file")
  gwas_col_types = c("ccicdddddddccc")
  gwas_pvalues = scanTabixDataFrame(summary_path, granges, col_names = gwas_col_names, col_types = gwas_col_types)
  return(gwas_pvalues)
}

scanTabixDataFrame <-function(tabix_file, param, ...){
  tabix_list = Rsamtools::scanTabix(tabix_file, param = param)
  df_list = lapply(tabix_list, function(x,...){
    if(length(x) > 0){
      if(length(x) == 1){
        #Hack to make sure that it also works for data frames with only one row
        #Adds an empty row and then removes it
        result = paste(paste(x, collapse = "\n"),"\n",sep = "")
        result = readr::read_delim(result, delim = "\t", ...)[1,]
      }else{
        result = paste(x, collapse = "\n")
        result = readr::read_delim(result, delim = "\t", ...)
      }
    } else{
      #Return NULL if the nothing is returned from tabix file
      result = NULL
    }
    return(result)
  }, ...)
  return(df_list)
}

summaryReplaceCoordinates <-function(summary_df, variant_information){

  #Make key assertions
  assertthat::assert_that(assertthat::has_name(summary_df, "snp_id"))
  assertthat::assert_that(assertthat::has_name(summary_df, "pos"))
  assertthat::assert_that(assertthat::has_name(summary_df, "chr"))

  #Filter variant information to contain only required snps
  var_info = dplyr::filter(variant_information, snp_id %in% summary_df$snp_id) %>%
    dplyr::select(snp_id, chr, pos, MAF)

  #Remove MAF if it is present
  if(assertthat::has_name(summary_df, "MAF")){
    summary_df = dplyr::select(summary_df, -MAF)
  }

  #Add new coordinates
  new_coords = dplyr::select(summary_df, -chr, -pos) %>%
    dplyr::left_join(var_info, by = "snp_id") %>%
    dplyr::filter(!is.na(pos)) %>%
    dplyr::arrange(pos)

  return(new_coords)
}
