suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("tidyr"))
suppressPackageStartupMessages(library("purrr"))
suppressPackageStartupMessages(library("coloc"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("devtools"))
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("purrrlyr"))
suppressPackageStartupMessages(library("Rsamtools"))

#detach("package:seqUtils", unload=TRUE)

#source("/nfs/users/nfs_n/np12/myscratch/MacroMap/Analysis/Coloc_analysis/Macromap_fds/coloc_new_gwas/IBD_Laura/MacroMap/scripts/functions_me.R")
#load_all("/nfs/users/nfs_n/np12/scratch_team170/ka8/projects/seqUtils")

#Parse command-line options
option_list <- list(
  make_option(c("-p", "--phenotype"), type="character", default=NULL,
              help="Type of QTLs used for coloc.", metavar = "type"),
  make_option(c("-w", "--window"), type="character", default=NULL,
              help="Size of the cis window.", metavar = "type"),
  make_option(c("--gwas"), type="character", default=NULL,
              help="Name of the GWAS trait", metavar = "type"),
  make_option(c("-d", "--dir"), type="character", default=NULL,
              help="Path to GWAS summary stats directory.", metavar = "type"),
  make_option(c("--qtl"), type="character", default=NULL,
              help="Path to the QTL directory.", metavar = "type"),
  make_option(c("-o", "--outdir"), type="character", default=NULL,
              help="Path to the output directory.", metavar = "type"),
  make_option(c("-s", "--samplesizes"), type="character", default=NULL,
              help="Path to the tab-separated text file with condition names and sample sizes.", metavar = "type"),
  make_option(c( "--gwasvarinfo"), type="character", default=NULL,
              help="Variant infromation file for the GWAS dataset.", metavar = "type"),
  make_option(c("--qtlvarinfo"), type="character", default=NULL,
              help="Variant information file for the QTL dataset.", metavar = "type"),
  make_option(c("--gwaslist"), type="character", default=NULL,
              help="Path to the list of GWAS studies.", metavar = "type"),
  make_option(c("--chunk"), type="integer", default=NULL,
            help="chunk to run", metavar = "type"),
  make_option(c("--numberGWAS"), type="integer", default=NULL,
              help="numberGWAS", metavar = "type"),
  make_option(c("--function_path_source"), type="character", default=NULL,
            help="Path to the functions_me file", metavar = "type")
)
opt <- parse_args(OptionParser(option_list=option_list))

#Debugging
# opt = list(gwas = "IBD", w = "2e6", p = "featureCounts", d = "/nfs/users/nfs_n/np12/myscratch/MacroMap/Analysis/Coloc_analysis/Macromap_fds/coloc_new_gwas/IBD_Laura/GWAS_data/",
#            o = "/nfs/users/nfs_n/np12/myscratch/MacroMap/Analysis/Coloc_analysis/Macromap_fds/coloc_new_gwas/IBD_Laura/MacroMap/coloc_out/",
#            qtl = "/nfs/users/nfs_n/np12/myscratch/MacroMap/Analysis/eQTLs/Macromap_fds/analysis/eQTLs_per_TSS/1MB/",
#            s = "/nfs/users/nfs_n/np12/myscratch/MacroMap/Analysis/eQTLs/Macromap_fds/analysis/eQTLs_per_TSS/sample_size_per_condition/IFNG_6",
#            gwasvarinfo = "/nfs/users/nfs_n/np12/myscratch/MacroMap/Data/Genotypes/vcf_hipsci_hg38/DS_Filt_M/chr_position_ID_major_minor_MAF_v37v2.txt.gz",
#            qtlvarinfo = "/nfs/users/nfs_n/np12/myscratch/MacroMap/Data/Genotypes/vcf_hipsci_hg38/DS_Filt_M/chr_position_ID_major_minor_MAF.txt.gz",
#            gwaslist = "/nfs/users/nfs_n/np12/myscratch/MacroMap/Analysis/Coloc_analysis/Macromap_fds/coloc_new_gwas/IBD_Laura/GWAS_data/IBD_gwas_files.txt",
#            chunk=6)

#Extract parameters for CMD options
gwas_id = opt$gwas
cis_window = as.numeric(opt$w)
phenotype = opt$p
gwas_dir = opt$d
qtl_dir = opt$qtl
outdir = opt$o
sample_size_path = opt$s
gwas_var_path = opt$gwasvarinfo
qtl_var_path = opt$qtlvarinfo
gwas_list = opt$gwaslist
chunk=opt$chunk
numberGWAS=opt$numberGWAS
function_path=opt$function_path_source

source(function_path)


#Import variant information
#gwas_var_info = importVariantInformation(gwas_var_path) # This set of variants is with v37 of the genome in case gwas summary stats are on 37
qtl_var_info = importVariantInformation(qtl_var_path) # v38
print("Variant information imported.")

#Import list of GWAS studies
gwas_stats_labeled = readr::read_tsv(gwas_list, col_names = c("trait","file_name","type"), col_type = "ccc")

#Import sample sizes
sample_sizes = readr::read_tsv(sample_size_path, col_names = c("condition_name", "sample_size"), col_types = "ci")
sample_sizes_list = as.list(sample_sizes$sample_size)
names(sample_sizes_list) = sample_sizes$condition_name


#Construct a new QTL list
phenotype_values = constructQtlListForColoc(qtl_dir, sample_sizes_list)

#Spcecify the location of the GWAS summary stats file
gwas_file_name = dplyr::filter(gwas_stats_labeled, trait == gwas_id)$file_name
gwas_prefix = file.path(gwas_dir, gwas_file_name)

#Prefilter coloc candidates
qtl_df_list = prefilterColocCandidates(phenotype_values$min_pvalues, gwas_prefix,
                                       gwas_variant_info = qtl_var_info, fdr_thresh = 1,
                                       overlap_dist = cis_window, gwas_thresh = 1e-5)
qtl_pairs = purrr::map_df(qtl_df_list, identity) %>% unique()
names(qtl_pairs)<-c("phenotype_id","snp_id")

# chunking
d<-1:length(qtl_pairs$phenotype_id)
d.split<-split(d, sort(d%%20))
qtl_pairs<-qtl_pairs[d.split[[chunk]],]

print(paste("Pre-filtering completed, number of genes to test",length(d.split[[chunk]]),"chunk",chunk))


#Test for coloc (if specified numberGWAS coloc is running with pvalues )

coloc_res_list = purrr::map2(phenotype_values$qtl_summary_list, phenotype_values$sample_sizes,
                             ~colocMolecularQTLsByRow_nikos(qtl_pairs, qtl_summary_path = .x,
                                                            gwas_summary_path = paste0(gwas_prefix, ".sorted.txt.gz"),
                                                            gwas_variant_info = qtl_var_info,
                                                            qtl_variant_info = qtl_var_info,
                                                            N_qtl = .y, cis_dist = cis_window,N=numberGWAS))

print("Coloc completed.")

#Export results
coloc_hits = purrr::map_df(coloc_res_list, identity, .id = "condition_name") %>% dplyr::arrange(gwas_lead) %>% dplyr::mutate(gwas_trait=gwas_id)
coloc_output = file.path(outdir, paste(gwas_id, phenotype, opt$w,sample_sizes$condition_name,"chunk",chunk,"txt", sep = "."))
write.table(coloc_hits, coloc_output, sep = "\t", quote = FALSE, row.names = FALSE)
