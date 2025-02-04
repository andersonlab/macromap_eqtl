#!/bin/bash

# tss window hash
declare -A WINDOW=( ["500K"]="500000" )

# declare conditions array for PCS
declare -A arr=( ["CIL_24"]="20" ["CIL_6"]="15" ["Ctrl_24"]="15" ["Ctrl_6"]="15" ["IFNB_24"]="25" ["IFNB_6"]="15" ["IFNG_24"]="20" ["IFNG_6"]="15" ["IL4_24"]="15" ["IL4_6"]="15" ["LIL10_24"]="15" ["LIL10_6"]="15" ["P3C_24"]="15" ["P3C_6"]="20" ["R848_24"]="15" ["R848_6"]="15" ["sLPS_24"]="20" ["sLPS_6"]="15" ["PIC_24"]="8" ["PIC_6"]="10" ["MBP_24"]="9" ["MBP_6"]="6" ["Prec_D0"]="6" ["Prec_D2"]="9" )

# old in order to rename the merged files 
declare -A arrB=( ["CIL_24"]="20" ["CIL_6"]="20" ["Ctrl_24"]="10" ["Ctrl_6"]="15" ["IFNB_24"]="20" ["IFNB_6"]="15" ["IFNG_24"]="15" ["IFNG_6"]="20" ["IL4_24"]="15" ["IL4_6"]="20" ["LIL10_24"]="20" ["LIL10_6"]="15" ["P3C_24"]="15" ["P3C_6"]="15" ["R848_24"]="15" ["R848_6"]="15" ["sLPS_24"]="20" ["sLPS_6"]="15" ["PIC_24"]="8" ["PIC_6"]="9" ["MBP_24"]="10" ["MBP_6"]="6" ["Prec_D0"]="1" ["Prec_D2"]="9")

for TSS in "${!WINDOW[@]}"; do

  for cond in "${!arr[@]}";do

  log=/nfs/users/nfs_n/np12/myscratch/MacroMap/Analysis/eQTLs/Sep_all_2018/$cond/analysis/$TSS/GENE_NOMINAL/log

#  echo "bsub -e $log/merge.nominal.err -o $log/merge.nominal.out 'cat /nfs/users/nfs_n/np12/myscratch/MacroMap/Analysis/eQTLs/Sep_all_2018/$cond/analysis/$TSS/GENE_NOMINAL/chunk*.txt |gzip -c > /nfs/users/nfs_n/np12/myscratch/MacroMap/Analysis/eQTLs/Sep_all_2018/$cond/analysis/$TSS/GENE_NOMINAL/$cond\_$TSS\_PC${arr[$cond]}_all.txt.gz;rm /nfs/users/nfs_n/np12/myscratch/MacroMap/Analysis/eQTLs/Sep_all_2018/$cond/analysis/$TSS/GENE_NOMINAL/chunk*.txt'"
echo " mv /nfs/users/nfs_n/np12/myscratch/MacroMap/Analysis/eQTLs/Sep_all_2018/$cond/analysis/$TSS/GENE_NOMINAL/$cond\_$TSS\_PC${arr[$cond]}_all.txt.gz /nfs/users/nfs_n/np12/myscratch/MacroMap/Analysis/eQTLs/Sep_all_2018/$cond/analysis/$TSS/GENE_NOMINAL/$cond\_$TSS\_PC${arrB[$cond]}_all.txt.gz"

  done
done
