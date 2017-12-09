
#/data1/Programs/kent/bin/x86_64/gtfToGenePred -ignoreGroupsWithoutExons /workdata2/Projects/xiao/Share_RPFmethod/HEK_293_Gao_etal/simulation_CDS/modified.gtf stdout | /data1/Programs/kent/bin/x86_64/genePredToBed stdin gencode.v19.annotation.bed
/workdata/Share_RPFmethod/programs/ORF-RATER-master/prune_transcripts.py --inbed gencode.v19.annotation.bed --summarytable tid_removal_summary.txt -p 16 -v /workdata2/Projects/xiao/Share_RPFmethod/sharedFiles/annotation/gencodev19/hg19_genome.fa $1 > 1prune.log
/workdata/Share_RPFmethod/programs/ORF-RATER-master/make_tfams.py -v > 2make_tfams.log
/workdata/Share_RPFmethod/programs/ORF-RATER-master/find_orfs_and_types.py /workdata2/Projects/xiao/Share_RPFmethod/sharedFiles/annotation/gencodev19/hg19_genome.fa --codons ATG -p 16 -v > 3find_ORF.log
/workdata/Share_RPFmethod/programs/ORF-RATER-master/psite_trimmed.py $1 --minrdlen 26 --maxrdlen 29 --subdir Ribo --tallyfile tallies.txt --cdsbed gencode.v19.annotation.bed -p 16 -v > 4psite.log
/workdata/Share_RPFmethod/programs/ORF-RATER-master/regress_orfs.py $1 --mincdsreads 10 --subdir Ribo -p 16 -v > 5regress.log
/workdata/Share_RPFmethod/programs/ORF-RATER-master/rate_regression_output.py Ribo -p 16 --CSV rate_regression.csv -v > 6rate.log
/workdata/Share_RPFmethod/programs/ORF-RATER-master/make_orf_bed.py --minrating 0 > 7make_orf.log
/workdata/Share_RPFmethod/programs/ORF-RATER-master/quantify_orfs.py Ribo.bam --subdir Ribo -p 16 -v --force --CSV quant.csv --minrating 0