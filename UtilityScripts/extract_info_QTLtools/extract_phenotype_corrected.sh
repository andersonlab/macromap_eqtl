while read line; do
   gene=$(echo "$line"| cut -f1 -d "_")
   snp=$(echo "$line"| cut -f2 -d "_")
   region=$(echo "$line"| cut -f3,4 -d "_"|tr "_" ":");
   cond=$(echo "$line"| cut -f5,6 -d "_");

echo $gene > significant_geneID
echo $snp > significant_snpID

declare -A WINDOW=( ["1MB"]="1000000" )

# declare conditions array for PCS

declare -A arr=(["CIL_24"]="40" ["CIL_6"]="40" ["Ctrl_24"]="40" ["Ctrl_6"]="40" ["IFNB_24"]="40" ["IFNB_6"]="35" ["IFNG_24"]="35" ["IFNG_6"]="35" ["IL4_24"]="40" ["IL4_6"]="50" ["LIL10_24"]="35" ["LIL10_6"]="50" ["MBP_24"]="40" ["MBP_6"]="50" ["P3C_24"]="40" ["P3C_6"]="50" ["PIC_24"]="40" ["PIC_6"]="50" ["Prec_D0"]="35" ["Prec_D2"]="35" ["R848_24"]="40" ["R848_6"]="50" ["sLPS_24"]="40" ["sLPS_6"]="50")


for TSS in "${!WINDOW[@]}"; do

  for cond in "${!arr[@]}";do
  #for I in A;do

    VCF=~/myscratch/MacroMap/Data/Genotypes/vcf_hipsci_hg38/DS_Filt_M/hipsci.wec.gtarray.HumanCoreExome.imputed_phased.INFO_0.4_filtered.20180102.genotypes.all_chr.hg38.sorted_unique_annotate_MAF5_filterMissing_added_DS.bcf.gz

    GENE=~/myscratch/MacroMap/Analysis/eQTLs/Macromap_fds/$cond/expression/exp_TPMs_fild_$cond\_QTLtools_input.txt.gz


    # PC 10
    OUT=./
    LOG=./
    COV=~/myscratch/MacroMap/Analysis/eQTLs/Macromap_fds/$cond/covariates/$cond\_${arr[$cond]}\.txt.gz

    echo "QTLtools correct --bed $GENE --out $OUT/exp_corrected_$cond.out --cov $COV --include-phenotypes significant_geneID --normal"
    echo "QTLtools extract --bed $GENE --vcf $VCF --out $OUT/snp_$cond.out --include-sites significant_snpID --region $region --include-phenotypes significant_geneID --cov $COV"

  done
done
done
