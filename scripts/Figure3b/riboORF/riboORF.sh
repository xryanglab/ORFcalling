#samtools view -h /workdata2/Projects/xiao/Share_RPFmethod/HEK_293_Gao_etal/remap/star_uniq/HEK_293_GaoAligned.sortedByCoord.out.bam  > HEK_293_Gao.sam
#step1 correct reads
#perl /workdata2/Projects/xiao/Share_RPFmethod/programs/ribORF.0.1/offsetCorrect.pl HEK_293_Gao.sam offset.corretion.parameters.txt corrected.mapping.sam
#step2 check 3-nt periodicity
#perl /workdata2/Projects/xiao/Share_RPFmethod/programs/ribORF.0.1/readDist.pl corrected.mapping.sam /workdata2/Projects/xiao/Share_RPFmethod/programs/ribORF.0.1/example.data/hg19.coding.gene.txt check_periodicity 1 30 50

#step3 generate ORF genePred file
python generate_CDS_ORF_genePred.py > selectedORF.genePred
sed -i '1d' selectedORF.genePred
mkdir result

#step4 result
perl /workdata2/Projects/xiao/Share_RPFmethod/programs/ribORF.0.1/ribORF.pl corrected.mapping.sam selectedORF.genePred result 10 10