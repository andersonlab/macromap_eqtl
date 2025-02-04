snakefile="$(basename -- $1)";
configfile="$(basename -- $2)";
clusterfile="$(basename -- $3)";
dry=$4;

if [ $dry = "dry" ]; then
snakemake --cluster 'bsub -e {cluster.error} -o {cluster.output} -M {cluster.memory} -R {cluster.resources} -n {cluster.cores} -J {cluster.name} -q {cluster.queue} -m {cluster.select_queue}' --cluster-config $PWD/scripts/$clusterfile --jobs 4000 --keep-going --snakefile $PWD/scripts/$snakefile --configfile $PWD/scripts/$configfile  -np

else
snakemake --cluster 'bsub -e {cluster.error} -o {cluster.output} -M {cluster.memory} -R {cluster.resources} -n {cluster.cores} -J {cluster.name} -q {cluster.queue} -m {cluster.select_queue}' --cluster-config /$PWD/scripts/$clusterfile --jobs 4000 --keep-going --snakefile $PWD/scripts/$snakefile --configfile $PWD/scripts/$configfile
fi
