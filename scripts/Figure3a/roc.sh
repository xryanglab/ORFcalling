#!/bin/env zsh
for i in ribocode rpbp riboORF orfrate
do
    grep -v None ${i}_roc.txt |awk '{print $1"_"$2"\t"$3"\t"$4}' >${i}_roc_noNA.txt
    Rscript roc.R ${i}_roc_noNA.txt ${i} > all_auc.txt
done
