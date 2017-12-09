#!/usr/bin/env python
# -*- coding:UTF-8 -*-
__author__ = 'Zhengtao Xiao'

#Note: change the following path:
RiboCode_library_path="/Library/Python/2.7/site-packages/RiboCode/"
selected_gtf = "selected.gtf"
ribocode_annot_dir = "ribocode_annot/"
import sys
sys.path.append(RiboCode_library_path)
from prepare_transcripts import *
select_genes = []
with open(selected_gtf) as fin:
	for line in fin:
		if line.startswith("#"): continue
		field_dict = parsing_line(line)
		select_genes.append(field_dict["attr"]["gene_id"])

gene_dict,transcript_dict = load_transcripts_pickle(ribocode_annot_dir + "transcripts.pickle")
for g in select_genes:
	transcripts = gene_dict[g].transcripts
	for t in transcripts:
		tobj = transcript_dict[t]
		strand = tobj.genomic_iv.strand
		if not tobj.genomic_cds:
			continue
		else:
			line = []
			line.append(g)
			line.append(t)
			line.append(tobj.chrom)
			line.append(strand)
			if strand == "+":
				line.append(tobj.genomic_cds[0].start)
				line.append(tobj.genomic_cds[-1].end)
				line.append(tobj.genomic_cds[0].start)
				line.append(tobj.genomic_cds[-1].end)
			else:
				line.append(tobj.genomic_cds[-1].start)
				line.append(tobj.genomic_cds[0].end)
				line.append(tobj.genomic_cds[-1].start)
				line.append(tobj.genomic_cds[0].end)
			line.append(len(tobj.genomic_cds))
			if strand == "+":
				line.append(",".join(map(str,[i.start for i in tobj.genomic_cds]))+ ",")
				line.append(",".join(map(str,[i.end for i in tobj.genomic_cds]))+ ",")
			else:
				line.append(",".join(map(str,[i.start for i in tobj.genomic_cds[::-1]]))+ ",")
				line.append(",".join(map(str,[i.end for i in tobj.genomic_cds[::-1]]))+ ",")
			print "\t".join(map(str,line))
