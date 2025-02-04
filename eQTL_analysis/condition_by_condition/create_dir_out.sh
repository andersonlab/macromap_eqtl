## declare an array variable

declare -a arr=("CIL_24" "CIL_6" "Ctrl_24" "Ctrl_6" "IFNB_24" "IFNB_6" "IFNG_24" "IFNG_6" "IL4_24" "IL4_6" "LIL10_24" "LIL10_6" "P3C_24" "P3C_6"  "R848_24" "R848_6" "sLPS_24" "sLPS_6"  "PIC_24" "PIC_6" "MBP_24" "MBP_6" "Prec_D0" "Prec_D2")


cd /nfs/users/nfs_n/np12/myscratch/MacroMap/Analysis/eQTLs/dir_stracture/

## loop through the above array
for cond in "${arr[@]}"
do
   find . -type d -exec mkdir -p -- /nfs/users/nfs_n/np12/myscratch/MacroMap/Analysis/eQTLs/Macromap_fds/$cond/{} \;

done
