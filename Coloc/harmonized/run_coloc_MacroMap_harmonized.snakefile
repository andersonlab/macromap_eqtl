localrules: run_all
rule run_all:
	input:
		expand("{coloc_output}/{gwas}.{phenotype}.{coloc_window}.{conditions}.chunk.{chunks}.txt", coloc_output = config["coloc_output"], gwas = config["gwas_traits"], phenotype = config["coloc_phenotypes"], coloc_window = config["coloc_window"],conditions =config["conditions"], chunks=config["chunks"],complete_run= config["complete_run"])
	output:
		expand("{path}/{file_name}.DONE",path=config["general_path"],file_name=config["complete_run"])
	resources:
		mem = 100
	threads: 1
	shell:
		"echo 'Done!' > {output}"


#Run coloc accross inflammatory traits
rule run_coloc:
	input:
		"/nfs/users/nfs_n/np12/myscratch/MacroMap/Analysis/eQTLs/Macromap_fds/analysis/eQTLs_per_TSS/1MB/nominal"
	output:
		"{coloc_output}/{gwas}.{phenotype}.{coloc_window}.{conditions}.chunk.{chunks}.txt"
	params:
		outdir = "{coloc_output}",
		qtl_dir = "/nfs/users/nfs_n/np12/myscratch/MacroMap/Analysis/eQTLs/Macromap_fds/analysis/eQTLs_per_TSS/1MB/",
		sample_sizes_dir = "/nfs/users/nfs_n/np12/myscratch/MacroMap/Analysis/eQTLs/Macromap_fds/analysis/eQTLs_per_TSS/sample_size_per_condition/",
		phenotype = "{phenotype}",
		gwas = "{gwas}",
		coloc_window = "{coloc_window}",
		conditions = "{conditions}",
		chunks="{chunks}"

	resources:

	shell:
		"Rscript {config[general_path]}/scripts/harmonized/GWAS_run_coloc.R "
		"--phenotype {wildcards.phenotype} --window {wildcards.coloc_window} "
		"--gwas {wildcards.gwas} --dir {config[gwas_dir]} --outdir {params.outdir} "
		"--qtl {params.qtl_dir} --samplesizes {params.sample_sizes_dir}{params.conditions} "
		"--gwasvarinfo {config[gwas_var_info]} --qtlvarinfo {config[qtl_var_info]} "
		"--gwaslist {config[gwas_list]} --chunk {params.chunks} --function_path_source {config[general_path]}/scripts/harmonized/functions_me.R "
		"--gwastype {config[gwastype]}"
