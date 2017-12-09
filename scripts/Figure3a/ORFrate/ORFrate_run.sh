#!/bin/env zsh
gtfToGenePred -ignoreGroupsWithoutExons simulation.gtf stdout | genePredToBed stdin gencode.v19.annotation.bed
/programs/ORF-RATER-master/prune_transcripts.py --inbed gencode.v19.annotation.bed --summarytable tid_removal_summary.txt -p 16 -v hg19_genome.fa $1 > 1prune.log
/programs/ORF-RATER-master/make_tfams.py -v > 2make_tfams.log
/programs/ORF-RATER-master/find_orfs_and_types.py hg19_genome.fa --codons ATG -p 16 -v > 3find_ORF.log
/programs/ORF-RATER-master/psite_trimmed.py Ribo.bam --minrdlen 26 --maxrdlen 29 --subdir Ribo --tallyfile tallies.txt --cdsbed gencode.v19.annotation.bed -p 16 -v > 4psite.log
/programs/ORF-RATER-master/regress_orfs.py Ribo.bam --mincdsreads 10 --subdir Ribo -p 16 -v > 5regress.log
/programs/ORF-RATER-master/rate_regression_output.py Ribo -p 16 --CSV rate_regression.csv -v > 6rate.log
/programs/ORF-RATER-master/make_orf_bed.py --minrating 0 > 7make_orf.log
/programs/ORF-RATER-master/quantify_orfs.py Ribo.bam --subdir Ribo -p 16 -v --force --CSV quant.csv --minrating 0
