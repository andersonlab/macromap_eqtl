
declare -a arr=("CIL_24" "CIL_6" "Ctrl_24" "Ctrl_6" "IFNB_24" "IFNB_6" "IFNG_24" "IFNG_6" "IL4_24" "IL4_6" "LIL10_24" "LIL10_6" "P3C_24" "P3C_6"  "R848_24" "R848_6" "sLPS_24" "sLPS_6"  "PIC_24" "PIC_6" "MBP_24" "MBP_6" "Prec_D0" "Prec_D2")

# tss window hash
declare -A WINDOW=( ["1MB"]="1000000" ["750K"]="750000" ["500K"]="500000" ["250K"]="250000" ["50K"]="50000")

for cond in "${arr[@]}"
do

for TSS in "${!WINDOW[@]}"; do

    echo find /nfs/users/nfs_n/np12/myscratch/MacroMap/Analysis/eQTLs/Macromap_fds/$cond/analysis/$TSS/NO_PC/log/  -name \"*out\" -exec grep -H \"Exited with\" {} \\\;

    for PC in 1 2 3 4 5 6 7 8 9 10 15 20 25 30 35 40 50 60 70 80 ; do

      #echo  find /nfs/users/nfs_n/np12/myscratch/MacroMap/Analysis/eQTLs/Sep2018/$cond/analysis/$TSS/PC$PC/log/ -name \"*out\" -exec grep -H \"Exited\" {} \\\;
      echo  find /nfs/users/nfs_n/np12/myscratch/MacroMap/Analysis/eQTLs/Macromap_fds/$cond/analysis/$TSS/PC$PC/log/ -name \"*out\" -exec grep -H \"Exited with\" {} \\\;

      done
    done
  done
