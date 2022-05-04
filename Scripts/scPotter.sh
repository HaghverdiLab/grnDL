#!/bin/sh
path=/fast/AG_Haghverdi/Carla_Moelbert/Benchmarking_GRN_construction_DataDriven/grnDL
path_scripts="$path/Scripts"
output="$path/Data/GRNs_full"



for tag in $(ls -d $output/PBMC*| xargs -n1  basename) ; do

   # qsub -N $tag -l m_mem_free=16G -l h_rt=1:00:00  -pe smp 4 \
   # -o $output/$tag/out_scpotter -e $output/$tag/error_scpotter \
   sh  $path_scripts/scPotter_functionCall.sh $tag  $output $output $path_scripts/
done