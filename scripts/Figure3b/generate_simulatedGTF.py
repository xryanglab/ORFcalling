#!/usr/bin/env python
# -*- coding:UTF-8 -*-
__author__ = 'Zhengtao Xiao'
#change the CDS to nonCDS

ribocode_path = "/data1/Programs/anaconda2/lib/python2.7/site-packages/RiboCode/"
gtfFile = "/workdata/DataBase/GRCh38R24/gencode.v24.primary_assembly.annotation.gtf"

import sys
sys.path.append(ribocode_path)
from prepare_transcripts import *
select_genes = set()
with open("filter_genes.gtf") as fin:
	for line in fin:
		if line.startswith("#"): continue
		field_dict = parsing_line(line)
		select_genes.add(field_dict["attr"]["gene_id"])
#
with open("modified.gtf","w") as fout,open(gtfFile) as fin:
	for line in fin:
		if line.startswith("#"): continue
		field_dict = parsing_line(line)
		if field_dict["attr"]["gene_id"] in select_genes:
			if field_dict["feature"] in ["CDS","start_codon","stop_codon"]:
				continue
			else:
				fout.write(line.replace("protein_coding","lincRNA"))
		else:
			fout.write(line)

