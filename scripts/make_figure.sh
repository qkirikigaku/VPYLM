#!bin/bash
cancer_type=$1
scp -r idenken:project/VPYLM/result/${cancer_type}_* ./result/
python Drawing/draw_LL_ACC.py ${cancer_type}
python Drawing/find_best_index.py ${cancer_type}
python Drawing/draw_ngram.py ${cancer_type}
python Drawing/draw_likelihood.py ${cancer_type} Up
python Drawing/draw_likelihood.py ${cancer_type} Down
