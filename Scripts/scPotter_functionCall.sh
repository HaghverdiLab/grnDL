#!/bin/sh

source  /home/cmoelbe/bin/anaconda3/bin/activate scPotter

start=`date +%s`

python3 $4/run_save.py \
                -i $2/ -tag $1 \
                --oo $3/$1/   \
                --classifier Chebnet \
                --n_hidden_GNN 32 \
                --n_hidden_FC 16 \
                --dropout 0.1 \
                --infer_graph False \
                --seed 0 \
                --K 2 \
                --epochs 30


end=`date +%s`

echo scPotter $1 `expr $end - $start` >> /fast/AG_Haghverdi/Carla_Moelbert/Incooperating_GRNs_Into_Celltype_Prediction/Results/runtime.txt



