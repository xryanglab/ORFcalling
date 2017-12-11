#!/usr/bin/env python
# -*- coding:UTF-8 -*-
__author__ = 'Zhengtao Xiao'

ribocode_path = "/lib/python2.7/site-packages/RiboCode/"
"""
sampling number of protein coding genes from counts file (counts number larger than specific number)
"""
PCG_SAMPLE_NUM = 1000
MIN_COUNTS = 5 # for selecting the genes 
accept_length = set([26,27,28,29]) # P-sites used for ORF calling
gtfFile = "gencode.v19.annotation.gtf"
bamFile1 = "rpf_Aligned.sortedByCoord.out.bam"
bamFile2 = "mRNA_Aligned.sortedByCoord.out.bam"
star_counts_file = 'rpf_ReadsPerGene.out.tab'

import random
from pybedtools import BedTool
import sys
sys.path.append(ribocode_path)
from prepare_transcripts import *

#step1 read start counts
gene_list = []
with open(star_counts_file) as fin:
	for line in fin:
		if line.startswith("N_"): continue
		tmp = line.strip().split("\t")
		if int(tmp[2]) > MIN_COUNTS:
			gene_list.append(tmp[0])

gene_list = set(gene_list)
coding_genes = set()
with open(gtfFile) as fin:
	for line in fin:
		if line.startswith("#"): continue
		field_dict = parsing_line(line)
		if field_dict["feature"] == "gene" and field_dict["attr"]["gene_id"] in gene_list:
			if field_dict["attr"]["gene_type"] == "protein_coding":
				coding_genes.add(field_dict["attr"]["gene_id"])


select_genes = random.sample(coding_genes,PCG_SAMPLE_NUM)
select_genes = set(select_genes)

#step2 read the gtf FIle:
with open("selected.gtf","w") as fout, open(gtfFile) as fin:
	for line in fin:
		if line.startswith("#"): continue
		field_dict = parsing_line(line)
		if field_dict["feature"] == "gene" and field_dict["attr"]["gene_id"] in select_genes:
			fout.write(line)

#step3 intersect with bamFile
select_bed = BedTool("selected.gtf")
bam1 = BedTool(bamFile1)
bam2 = BedTool(bamFile2)

#results_intersect1 = bam1.intersect(select_bed,s=True)
results_substract1 = bam1.intersect(select_bed,v=True,s=True)

results_intersect2 = bam2.intersect(select_bed,s=True)
#results_substract2 = bam1.intersect(select_bed,v=True,s=True)

#step4 output fastq
with open("simulated_rpf.fq","w") as fout:

	#write to RPF
	for i in results_substract1:
		fout.write("@%s\n%s\n+\n%s\n" % (i.fields[0],i.fields[9],i.fields[10]))
	for i in results_intersect2:
		seq = i.fields[9]
		qual = i.fields[10]
		if len(seq) not in accept_length:
			if len(seq) < min(accept_length):
				continue
			elif min(accept_length) <= len(seq) <= max(accept_length):
				pass
			else:
				trim = len(seq) - random.sample(accept_length,1)[0]
				trim_direction = random.choice([1,-1])
				if trim_direction > 0:
					seq = seq[trim:]
					qual = qual[trim:]
				else:
					seq = seq[:-trim]
					qual = qual[:-trim]
		fout.write("@%s\n%s\n+\n%s\n" % (i.fields[0],seq,qual))

