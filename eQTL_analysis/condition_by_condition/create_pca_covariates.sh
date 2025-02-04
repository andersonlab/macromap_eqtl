#!/bin/bash

POP_PCS=/lustre/scratch117/cellgen/team170/np12/MacroMap/Data/Genotypes/vcf_hipsci_hg38/DS_Filt_M/stratification/eigenstrat/filtered_final/eigenstrat_hipsci_filtered_3PCs.txt

declare -a arr=("CIL_24" "CIL_6" "Ctrl_24" "Ctrl_6" "IFNB_24" "IFNB_6" "IFNG_24" "IFNG_6" "IL4_24" "IL4_6" "LIL10_24" "LIL10_6" "P3C_24" "P3C_6"  "R848_24" "R848_6" "sLPS_24" "sLPS_6"  "PIC_24" "PIC_6" "MBP_24" "MBP_6" "Prec_D0" "Prec_D2")

for cond in "${arr[@]}"
do
EXP_PCS_PATH=/lustre/scratch117/cellgen/team170/np12/MacroMap/Analysis/eQTLs/Macromap_fds/$cond/covariates
Rscript  /nfs/users/nfs_n/np12/myscratch/MacroMap/Analysis/eQTLs/Macromap_fds/scripts/create_strat_covariates.R $POP_PCS $EXP_PCS_PATH/$cond.pca $EXP_PCS_PATH/$cond.pop_3_PCs $EXP_PCS_PATH/$cond.pca.fltd

for pc in 1 2 3 4 5 6 7 8 9 10 15 20 25 30 35 40 50 60 70 80 ; do
  cat <(head -1 $EXP_PCS_PATH/$cond.pca.fltd) <(cat $EXP_PCS_PATH/$cond.pop_3_PCs|sed '1d'|head -n 3) <(cat $EXP_PCS_PATH/$cond.pca.fltd | sed '1d' | head -n $pc) | gzip -c > $EXP_PCS_PATH/$cond\_$pc\.txt.gz
done
done
