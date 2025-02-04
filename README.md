# **Gene expression QTL mapping in stimulated iPSC-derived macrophages provides insights into common complex diseases**

## **Overview**
This repository contains data, links to data, scripts, and analysis pipelines used for this study.


## **Abstract**
Many disease-associated variants are thought to be regulatory but are not present in existing eQTL catalogues. We hypothesize that these variants may regulate gene expression in **stimulated immune cells**. Using human iPSC-derived macrophages, we mapped eQTLs across **24 stimulation conditions** and found that:
- 76% of eQTLs detected in at least one stimulated condition were also present in naive cells.
- Response eQTLs (reQTLs) varied across conditions, ranging from **3.7% to 28.4%**.
- reQTLs specific to a single condition were rare (**1.11%**), but significantly overrepresented among disease-colocalizing eQTLs.
- We identified **21.7% additional disease effector genes** through reQTL colocalization, with **38.6% of these not found in the GTEx**.

Our findings underscore the importance of **context-specific regulatory variation** in understanding **common disease risk**.

---

## **Repository Structure**
```
.
├── Coloc                         # Colocalization analysis scripts and configs
├── Data                          # Processed/helper data for figures 
├── Differential_expression_analysis # DESeq2 scripts for differential expression
├── eQTL_analysis                 # eQTL mapping condition by condition and mashR
├── Figures                       # R scripts for figure generation
├── UtilityScripts                # Helper scripts for analysis and visualization
├── Variance_component_analysis   # Variance decomposition scripts of gene expression
```

---

## **Data Availability**
- **RNA-seq Data**  
  The dataset includes **4,698 RNA-seq libraries** from **209 iPSC lines**, spanning **24 stimulation conditions**.  
  - [Zenodo Repository](https://zenodo.org/records/11563707) (All `.txt.gz,tar` files are available here)

- **Data Visualization Portal**  
  Explore the data interactively through our visualization portal:  
  - [MacroMap Portal](https://www.macromapqtl.org.uk/)

### **Expression Data**
- MacroMap_raw_expression.txt.gz  
- MacroMap_TPMs_expression.txt.gz  
- MacroMap_TPMs_filtered_expression.txt.gz  

### **Colocalization Results**
- MacroMap_GTEx_coloc_results.txt.gz  

### **eQTL Summary Statistics**
- MacroMap_eQTL_nominal_summary_stats.tar  
- MacroMap_eQTL_permuted_summary_stats.tar  

## **Contact**
For questions  please contact **[Nikos Panousis]** at **[nikos.panousis@gmail.com]**.

---