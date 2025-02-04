#!/bin/bash

declare -a arr=("CIL_24" "CIL_6" "Ctrl_24" "Ctrl_6" "IFNB_24" "IFNB_6" "IFNG_24" "IFNG_6" "IL4_24" "IL4_6" "LIL10_24" "LIL10_6" "P3C_24" "P3C_6"  "R848_24" "R848_6" "sLPS_24" "sLPS_6"  "PIC_24" "PIC_6" "MBP_24" "MBP_6" "Prec_D0" "Prec_D2")


# tss window hash
declare -A WINDOW=( ["1MB"]="1000000" ["750K"]="750000" ["500K"]="500000" ["250K"]="250000" ["50K"]="50000")

for cond in "${arr[@]}"
do

for TSS in "${!WINDOW[@]}"; do
  VCF=/lustre/scratch117/cellgen/team170/np12/MacroMap/Data/Genotypes/vcf_hipsci_hg38/DS_Filt_M/hipsci.wec.gtarray.HumanCoreExome.imputed_phased.INFO_0.4_filtered.20180102.genotypes.all_chr.hg38.sorted_unique_annotate_MAF5_filterMissing_added_DS.bcf.gz
  GENE=/lustre/scratch117/cellgen/team170/np12/MacroMap/Analysis/eQTLs/Macromap_fds/$cond/expression/exp_TPMs_fild_$cond\_QTLtools_input.txt.gz
  # NO correction
  OUT_NO_PC=/nfs/users/nfs_n/np12/myscratch/MacroMap/Analysis/eQTLs/Macromap_fds/$cond/analysis/$TSS/NO_PC
  LOG_NO_PC=/nfs/users/nfs_n/np12/myscratch/MacroMap/Analysis/eQTLs/Macromap_fds/$cond/analysis/$TSS/NO_PC/log
  DEL=/lustre/scratch117/cellgen/team170/np12/MacroMap/Analysis/eQTLs/Macromap_fds/exlude_genes_50K.txt

  if [ "$TSS" == "50K" ]; then
    echo "bsub -J \"QTL_$TSS.no_pc.$cond[1-50]\"  -e $LOG_NO_PC/chunk_$TSS\_%I.err -o $LOG_NO_PC/chunk_$TSS\_%I.out -R \"select[mem>2000] rusage[mem=2000]\" -M 2000  \"QTLtools cis --vcf $VCF --bed $GENE --out $OUT_NO_PC/chunk_$TSS.\\\$LSB_JOBINDEX.output.txt --log $LOG_NO_PC/chunk_$TSS.\\\$LSB_JOBINDEX.log --permute 1000 --exclude-phenotypes $DEL --window ${WINDOW[$TSS]} --normal --chunk \\\$LSB_JOBINDEX 50\""

    for PC in 1 2 3 4 5 6 7 8 9 10 15 20 25 30 35 40 50 60 70 80 ; do
      #for PC in 0;do
        COV=/lustre/scratch117/cellgen/team170/np12/MacroMap/Analysis/eQTLs/Macromap_fds/$cond/covariates/$cond\_$PC\.txt.gz
        OUT=/nfs/users/nfs_n/np12/myscratch/MacroMap/Analysis/eQTLs/Macromap_fds/$cond/analysis/$TSS/PC$PC
        LOG=/nfs/users/nfs_n/np12/myscratch/MacroMap/Analysis/eQTLs/Macromap_fds/$cond/analysis/$TSS/PC$PC/log


       echo "bsub -J \"QTL_$TSS.$cond.$PC[1-50]\"  -e $LOG/chunk_$TSS\_%I.err -o $LOG/chunk_$TSS\_%I.out -R \"select[mem>2000] rusage[mem=2000]\" -M 2000  \"QTLtools cis --vcf $VCF --bed $GENE --out $OUT/chunk_$TSS.\\\$LSB_JOBINDEX.output.txt --log $LOG/chunk_$TSS.\\\$LSB_JOBINDEX.log --permute 1000  --cov $COV --exclude-phenotypes $DEL --window ${WINDOW[$TSS]} --normal --chunk \\\$LSB_JOBINDEX 50\""

  done
  else
    echo "bsub -J \"QTL_$TSS.no_pc.$cond[1-50]\"  -e $LOG_NO_PC/chunk_$TSS\_%I.err -o $LOG_NO_PC/chunk_$TSS\_%I.out -R \"select[mem>2000] rusage[mem=2000]\" -M 2000  \"QTLtools cis --vcf $VCF --bed $GENE --out $OUT_NO_PC/chunk_$TSS.\\\$LSB_JOBINDEX.output.txt --log $LOG_NO_PC/chunk_$TSS.\\\$LSB_JOBINDEX.log --permute 1000 --window ${WINDOW[$TSS]} --normal --chunk \\\$LSB_JOBINDEX 50\""

    for PC in 1 2 3 4 5 6 7 8 9 10 15 20 25 30 35 40 50 60 70 80 ; do
    #for PC in 0;do
      COV=/lustre/scratch117/cellgen/team170/np12/MacroMap/Analysis/eQTLs/Macromap_fds/$cond/covariates/$cond\_$PC\.txt.gz
      OUT=/nfs/users/nfs_n/np12/myscratch/MacroMap/Analysis/eQTLs/Macromap_fds/$cond/analysis/$TSS/PC$PC
      LOG=/nfs/users/nfs_n/np12/myscratch/MacroMap/Analysis/eQTLs/Macromap_fds/$cond/analysis/$TSS/PC$PC/log


     echo "bsub -J \"QTL_$TSS.$cond.$PC[1-50]\"  -e $LOG/chunk_$TSS\_%I.err -o $LOG/chunk_$TSS\_%I.out -R \"select[mem>2000] rusage[mem=2000]\" -M 2000  \"QTLtools cis --vcf $VCF --bed $GENE --out $OUT/chunk_$TSS.\\\$LSB_JOBINDEX.output.txt --log $LOG/chunk_$TSS.\\\$LSB_JOBINDEX.log --permute 1000  --cov $COV --window ${WINDOW[$TSS]} --normal --chunk \\\$LSB_JOBINDEX 50\""

   done
fi
 done
done
