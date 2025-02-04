declare -a arr=("CIL_24" "CIL_6" "Ctrl_24" "Ctrl_6" "IFNB_24" "IFNB_6" "IFNG_24" "IFNG_6" "IL4_24" "IL4_6" "LIL10_24" "LIL10_6" "P3C_24" "P3C_6"  "R848_24" "R848_6" "sLPS_24" "sLPS_6"  "PIC_24" "PIC_6" "MBP_24" "MBP_6" "Prec_D0" "Prec_D2")

for cond in "${arr[@]}"
do
  find $PWD/$cond/analysis/*  -path $PWD/$cond/analysis/results -prune -o -name "*GENE.txt.gz"  -exec ln -s {} /nfs/users/nfs_n/np12/myscratch/MacroMap/Analysis/eQTLs/Macromap_fds/$cond/analysis/results/ \;
  find $PWD/$cond/analysis/*  -path $PWD/$cond/analysis/results -prune -o -name "*significant.txt"  -exec ln -s {} /nfs/users/nfs_n/np12/myscratch/MacroMap/Analysis/eQTLs/Macromap_fds/$cond/analysis/results/ \;
done
