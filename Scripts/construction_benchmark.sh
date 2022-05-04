#!/bin/sh
path=/fast/AG_Haghverdi/Carla_Moelbert/Benchmarking_GRN_construction_DataDriven/grnDL
path_scripts="$path/Scripts"

output="$path/Data/GRNs_subsets"
folder_test="$path/Data/PBMC_mono"
folder_train="$path/Data/subsets"


for score_threshold in 0 100; do # 0 
qsub -N GRN -l m_mem_free=16G -l h_rt=1:30:00 \
  $path_scripts/functioncall_construction_visual.sh $path/Data/PBMC_Mosaic $folder_test $path/Data/GRNs_full L1 Spearman  $score_threshold 5000
 
qsub -N GRN -l m_mem_free=16G -l h_rt=1:30:00 \
 $path_scripts/functioncall_construction_visual.sh $folder_test $folder_test $path/Data/GRNs_full L1 Spearman  $score_threshold 5000

#for folder in $folder_train/Mo*; do
#qsub -N GRN -l m_mem_free=16G -l h_rt=1:30:00 \
# $path_scripts/functioncall_construction_visual.sh $folder $folder_test $output L1 Spearman  $score_threshold 5000
#done 
done

