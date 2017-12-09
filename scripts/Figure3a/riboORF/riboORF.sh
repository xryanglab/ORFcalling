#!/bin/env zsh
samtools view -h simulation_Aligned.sortedByCoord.out.bam  > HEK_293_Gao.sam
#step1 correct reads
perl /programs/ribORF.0.1/offsetCorrect.pl HEK_293_Gao.sam offset.corretion.parameters.txt corrected.mapping.sam

#step2 generate ORF genePred file
python generate_CDS_ORF_genePred.py > selectedORF.genePred
sed -i '1d' selectedORF.genePred
mkdir result

#step3 result
perl /programs/ribORF.0.1/ribORF.pl corrected.mapping.sam selectedORF.genePred result 10 10