# copy files and keep dir stracture
find /nfs/users/nfs_n/np12/myscratch/MacroMap/Analysis/eQTLs_functional_enrichment/MacroMap_fds/fgwas/prepare_files/condition_by_condition/Analysis/Out/  -name "variant_level.bin" -exec sh -c 'dir=$(dirname "{}");mkdir -p "./${dir##*/}" ;cp "{}" "./${dir##*/}/$(basename "{}")"' \;

