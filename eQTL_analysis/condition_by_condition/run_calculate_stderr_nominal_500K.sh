declare -A WINDOW=( ["500K"]="500000" )

# declare conditions array
declare -A arr=( ["CIL_24"]="20" ["CIL_6"]="20" ["Ctrl_24"]="10" ["Ctrl_6"]="15" ["IFNB_24"]="20" ["IFNB_6"]="15" ["IFNG_24"]="15" ["IFNG_6"]="20" ["IL4_24"]="15" ["IL4_6"]="20" ["LIL10_24"]="20" ["LIL10_6"]="15" ["P3C_24"]="15" ["P3C_6"]="15" ["R848_24"]="15" ["R848_6"]="15" ["sLPS_24"]="20" ["sLPS_6"]="15" ["PIC_24"]="8" ["PIC_6"]="9" ["MBP_24"]="10" ["MBP_6"]="6" ["Prec_D0"]="1" ["Prec_D2"]="9")

for TSS in "${!WINDOW[@]}"; do

  for cond in "${!arr[@]}";do
  #for I in A;do

    log=/nfs/users/nfs_n/np12/myscratch/MacroMap/Analysis/eQTLs/Sep_all_2018/$cond/analysis/$TSS/GENE_NOMINAL/log
    input=/nfs/users/nfs_n/np12/myscratch/MacroMap/Analysis/eQTLs/Sep_all_2018/$cond/analysis/$TSS/GENE_NOMINAL/$cond\_$TSS\_PC${arr[$cond]}_all.txt.gz
    input_QTLs=/nfs/users/nfs_n/np12/myscratch/MacroMap/Analysis/eQTLs/Sep_all_2018/$cond/analysis/$TSS/PC${arr[$cond]}/$TSS\_${arr[$cond]}\.GENE.txt.gz

    output=/nfs/users/nfs_n/np12/myscratch/MacroMap/Analysis/eQTLs/Sep_all_2018/$cond/analysis/$TSS/GENE_NOMINAL/$cond\_$TSS\_PC${arr[$cond]}_all

    echo "bsub -e $log/calculate_std.err -o $log/calculate_std.out  -R \"select[mem>40000] rusage[mem=40000]\" -M 40000  \"Rscript /nfs/users/nfs_n/np12/myscratch/MacroMap/Analysis/eQTLs/Sep_all_2018/scripts/calculate_stderr_nominal.R $input $input_QTLs $output\""


    done
    done
