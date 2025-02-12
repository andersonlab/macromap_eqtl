#!/bin/bash

# tss window hash
declare -A WINDOW=( ["1MB"]="1000000" )

# declare conditions array for PCS

declare -A arr=(["CIL_24"]="40"	["CIL_6"]="40"	["Ctrl_24"]="40"		["Ctrl_6"]="40"	["IFNB_24"]="40"	["IFNB_6"]="35"	["IFNG_24"]="35"	["IFNG_6"]="35"	["IL4_24"]="40"	["IL4_6"]="50"	["LIL10_24"]="35"	["LIL10_6"]="50"	["MBP_24"]="40"	["MBP_6"]="50"	["P3C_24"]="40"	["P3C_6"]="50"	["PIC_24"]="40"	["PIC_6"]="50"	["Prec_D0"]="35"	["Prec_D2"]="35"	["R848_24"]="40"	["R848_6"]="50"	["sLPS_24"]="40" ["sLPS_6"]="50")


for TSS in "${!WINDOW[@]}"; do

  for cond in "${!arr[@]}";do
  #for I in A;do

    VCF=/lustre/scratch117/cellgen/team170/np12/MacroMap/Data/Genotypes/vcf_hipsci_hg38/DS_Filt_M/hipsci.wec.gtarray.HumanCoreExome.imputed_phased.INFO_0.4_filtered.20180102.genotypes.all_chr.hg38.sorted_unique_annotate_MAF5_filterMissing_added_DS.bcf.gz

    GENE=/lustre/scratch117/cellgen/team170/np12/MacroMap/Analysis/eQTLs/Macromap_fds/$cond/expression/exp_TPMs_fild_$cond\_QTLtools_input.txt.gz
    #DEL=/lustre/scratch117/cellgen/team170/np12/MacroMap/Analysis/eQTLs/Macromap_fds/exclude_phenotypes.txt
    # PC 10
    OUT_PC_10=/nfs/users/nfs_n/np12/myscratch/MacroMap/Analysis/eQTLs/Macromap_fds/$cond/analysis/$TSS/GENE_NOMINAL/
    LOG_PC_10=/nfs/users/nfs_n/np12/myscratch/MacroMap/Analysis/eQTLs/Macromap_fds/$cond/analysis/$TSS/GENE_NOMINAL/log
    COV=/lustre/scratch117/cellgen/team170/np12/MacroMap/Analysis/eQTLs/Macromap_fds/$cond/covariates/$cond\_${arr[$cond]}\.txt.gz


    echo "bsub -J \"QTL_$cond.PC_${arr[$cond]}[1-50]\"  -e $LOG_PC_10/chunk_$cond\_%I.err -o $LOG_PC_10/chunk_$cond\_%I.out -R \"select[mem>1500] rusage[mem=1500]\" -M 1500  \"/nfs/users/nfs_n/np12/src/QTLtools/qtltools_compiled/QTLtools_1.1_Ubuntu12.04_x86_64 cis --vcf $VCF --bed $GENE --out $OUT_PC_10/chunk_$cond.\\\$LSB_JOBINDEX.output.txt --log $LOG_PC_10/chunk_$cond.\\\$LSB_JOBINDEX.log --nominal 1  --window ${WINDOW[$TSS]} --normal --cov $COV  --chunk \\\$LSB_JOBINDEX 50\""

  done
done
