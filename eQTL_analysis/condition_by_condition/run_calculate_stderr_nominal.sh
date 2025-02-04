# tss window hash
declare -A WINDOW=( ["1MB"]="1000000" )

# declare conditions array for PCS

declare -A arr=(["CIL_24"]="40"	["CIL_6"]="40"	["Ctrl_24"]="40"		["Ctrl_6"]="40"	["IFNB_24"]="40"	["IFNB_6"]="35"	["IFNG_24"]="35"	["IFNG_6"]="35"	["IL4_24"]="40"	["IL4_6"]="50"	["LIL10_24"]="35"	["LIL10_6"]="50"	["MBP_24"]="40"	["MBP_6"]="50"	["P3C_24"]="40"	["P3C_6"]="50"	["PIC_24"]="40"	["PIC_6"]="50"	["Prec_D0"]="35"	["Prec_D2"]="35"	["R848_24"]="40"	["R848_6"]="50"	["sLPS_24"]="40" ["sLPS_6"]="50")



for TSS in "${!WINDOW[@]}"; do

  for cond in "${!arr[@]}";do
  #for I in A;do

    log=/nfs/users/nfs_n/np12/myscratch/MacroMap/Analysis/eQTLs/Macromap_fds/$cond/analysis/$TSS/GENE_NOMINAL/log
    input=/nfs/users/nfs_n/np12/myscratch/MacroMap/Analysis/eQTLs/Macromap_fds/$cond/analysis/$TSS/GENE_NOMINAL/$cond\_$TSS\_PC${arr[$cond]}_all.txt.gz
    input_QTLs=/nfs/users/nfs_n/np12/myscratch/MacroMap/Analysis/eQTLs/Macromap_fds/$cond/analysis/$TSS/PC${arr[$cond]}/$TSS\_${arr[$cond]}\.GENE.txt.gz

    output=/nfs/users/nfs_n/np12/myscratch/MacroMap/Analysis/eQTLs/Macromap_fds/$cond/analysis/$TSS/GENE_NOMINAL/$cond\_$TSS\_PC${arr[$cond]}_all

    echo "bsub -e $log/calculate_std.err -o $log/calculate_std.out  -R \"select[mem>50000] rusage[mem=50000]\" -M 50000  \"Rscript /nfs/users/nfs_n/np12/myscratch/MacroMap/Analysis/eQTLs/Macromap_fds/scripts/calculate_stderr_nominal.R $input $input_QTLs $output\""


    done
    done
