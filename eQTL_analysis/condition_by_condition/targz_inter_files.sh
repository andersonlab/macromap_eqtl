declare -a arr=("CIL_24" "CIL_6" "Ctrl_24" "Ctrl_6" "IFNB_24" "IFNB_6" "IFNG_24" "IFNG_6" "IL4_24" "IL4_6" "LIL10_24" "LIL10_6" "P3C_24" "P3C_6"  "R848_24" "R848_6" "sLPS_24" "sLPS_6"  "PIC_24" "PIC_6" "MBP_24" "MBP_6" "Prec_D0" "Prec_D2")

# tss window hash

for cond in "${arr[@]}"
  do
    find $cond/*  -type f -name "*out"  -printf '%h\n' |sort -u |while read line; do echo "tar -zcf $line/out.tgz $line/*out --remove-files"; done
    find $cond/*  -type f -name "*err"  -printf '%h\n' |sort -u |while read line; do echo "tar -zcf $line/err.tgz $line/*err --remove-files"; done
    find $cond/*  -type f -name "*log"  -printf '%h\n' |sort -u |while read line; do echo "tar -zcf $line/log.tgz $line/*log --remove-files"; done
done
