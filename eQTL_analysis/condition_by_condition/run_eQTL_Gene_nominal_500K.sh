#!/bin/bash

# tss window hash
declare -A WINDOW=( ["500K"]="500000" )

# declare conditions hash per PC
declare -A arr=( ["CIL_24"]="20" ["CIL_6"]="20" ["Ctrl_24"]="10" ["Ctrl_6"]="15" ["IFNB_24"]="20" ["IFNB_6"]="15" ["IFNG_24"]="15" ["IFNG_6"]="20" ["IL4_24"]="15" ["IL4_6"]="20" ["LIL10_24"]="20" ["LIL10_6"]="15" ["P3C_24"]="15" ["P3C_6"]="15" ["R848_24"]="15" ["R848_6"]="15" ["sLPS_24"]="20" ["sLPS_6"]="15" ["PIC_24"]="8" ["PIC_6"]="9" ["MBP_24"]="10" ["MBP_6"]="6" ["Prec_D0"]="1" ["Prec_D2"]="9")


for TSS in "${!WINDOW[@]}"; do

  for cond in "${!arr[@]}";do
  #for I in A;do

    VCF=/lustre/scratch117/cellgen/team170/np12/MacroMap/Data/Genotypes/vcf_hipsci_hg38/DS_Filt_M/hipsci.wec.gtarray.HumanCoreExome.imputed_phased.INFO_0.4_filtered.20180102.genotypes.all_chr.hg38.sorted_unique_annotate_MAF5_filterMissing_added_DS.bcf.gz

    GENE=/lustre/scratch117/cellgen/team170/np12/MacroMap/Analysis/eQTLs/Sep_all_2018/$cond/expression/exp_TPMs_fild_$cond\_QTLtools_input.txt.gz
    DEL=/lustre/scratch117/cellgen/team170/np12/MacroMap/Analysis/eQTLs/Sep_all_2018/exclude_phenotypes.txt
    # PC 10
    OUT_PC_10=/nfs/users/nfs_n/np12/myscratch/MacroMap/Analysis/eQTLs/Sep_all_2018/$cond/analysis/$TSS/GENE_NOMINAL/
    LOG_PC_10=/nfs/users/nfs_n/np12/myscratch/MacroMap/Analysis/eQTLs/Sep_all_2018/$cond/analysis/$TSS/GENE_NOMINAL/log
    COV=/lustre/scratch117/cellgen/team170/np12/MacroMap/Analysis/eQTLs/Sep_all_2018/$cond/covariates/$cond\_${arr[$cond]}\.txt.gz


    echo "bsub -J \"QTL_$cond.PC_${arr[$cond]}[1-50]\"  -e $LOG_PC_10/chunk_$cond\_%I.err -o $LOG_PC_10/chunk_$cond\_%I.out -R \"select[mem>500] rusage[mem=500]\" -M 500  \"/nfs/users/nfs_n/np12/src/QTLtools/qtltools_compiled/QTLtools_1.1_Ubuntu12.04_x86_64 cis --vcf $VCF --bed $GENE --out $OUT_PC_10/chunk_$cond.\\\$LSB_JOBINDEX.output.txt --log $LOG_PC_10/chunk_$cond.\\\$LSB_JOBINDEX.log --nominal 1  --window ${WINDOW[$TSS]} --normal --cov $COV --exclude-phenotypes $DEL --chunk \\\$LSB_JOBINDEX 50\""

  done
done
