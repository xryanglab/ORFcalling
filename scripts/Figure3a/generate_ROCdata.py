#!/usr/bin/env python
# -*- coding:UTF-8 -*-
__author__ = 'Zhengtao Xiao'
"""
generate the roc data using the results of each method
usage: python generate_ROCdata.py
"""
###Note: please change the path of following file
RiboCode_library_path="/Library/Python/2.7/site-packages/RiboCode/"
ribocode_annot_dir = "ribocode_annot"
selected_gtf = "selected.gtf"

#ribocode
ribocode_result_true = "RiboCode_results/results_true/RiboCode.txt"
ribocode_result_simulation = "RiboCode_results/results_simulation/RiboCode.txt"
#rpbp
rpbp_result_true = "RpBp_results/results_true/hek293.alignments-only-unique.length-26-27-28-29.offset-12-12-12-12.frac-0.2.rw-0.bayes-factors.bed.gz"
rpbp_result_simulation = "RpBp_results/results_simulation/hek293.alignments-only-unique.length-26-27-28-29.offset-12-12-12-12.frac-0.2.rw-0.bayes-factors.bed.gz"
#riboORF
riboORF_result_true = "riboORF_results/results_true/pred.pvalue.parameters.txt"
riboORF_result_simulation ="riboORF_results/results_simulation/pred.pvalue.parameters.txt" 
#ORFrate
orfrate_result_true = "ORFrate_results/results_true/quant.csv"
orfrate_result_simulation = "ORFrate_results/results_simulation/quant.csv"

#check the files:
import os,sys
for i in [selected_gtf,ribocode_result_true,ribocode_result_simulation,rpbp_result_true,rpbp_result_simulation,riboORF_result_true,riboORF_result_simulation,orfrate_result_true,orfrate_result_simulation]:
    if not os.path.exists(i):
        print "%s not found!" % i
        sys.exit()

import csv
import gzip
sys.path.append(RiboCode_library_path)
from prepare_transcripts import *
from collections import defaultdict

select_genes = set()
with open(selected_gtf) as fin:
    for line in fin:
        if line.startswith("#"): continue
        field_dict = parsing_line(line)
        select_genes.add(field_dict["attr"]["gene_id"])
if ribocode_annot_dir[-1]!="/": ribocode_annot_dir+="/"
gene_dict,transcript_dict = load_transcripts_pickle(ribocode_annot_dir + "transcripts.pickle")

##ribocode:
ribocode_dict = defaultdict(lambda:  ["None","None"])
with open(ribocode_result_true) as fin:
    csvreader = csv.DictReader(fin,delimiter="\t")
    for row in csvreader:
        gene_id = row["gene_id"]
        transcript_id = row["transcript_id"]
        pv = row["pval_combined"]
        if gene_id in select_genes:
            if not transcript_dict[transcript_id].cds: continue
            if row["ORF_type"] == "annotated":
                ribocode_dict[gene_id + "\t" + transcript_id][0] = pv
with open(ribocode_result_simulation) as fin:
    csvreader = csv.DictReader(fin,delimiter="\t")
    for row in csvreader:
        gene_id = row["gene_id"]
        transcript_id = row["transcript_id"]
        pv = row["pval_combined"]
        if gene_id in select_genes:
            if not transcript_dict[transcript_id].cds: continue
            if transcript_dict[transcript_id].cds.end  + 3 == int(row["ORF_tstop"]):
                ribocode_dict[gene_id + "\t" + transcript_id][1] = pv

#rpbp:
rpbp_dict = defaultdict(lambda:  ["None","None"])
with gzip.open(rpbp_result_true) as fin:
    csvreader = csv.DictReader(fin,delimiter="\t")
    for row in csvreader:
        tid = row["#id"].split("_")[0]
        gid = transcript_dict[tid].gene_id
        score = row["#bayes_factor_mean"]
        if score in ["inf","-inf"]: continue
        orf_start = row["#start"]
        orf_stop = row["#end"]
        if gid in select_genes:
            if not transcript_dict[tid].cds: continue
            if transcript_dict[tid].genomic_iv.strand == "+":
                cds_stop = transcript_dict[tid].genomic_cds[-1].end
                if cds_stop == int(orf_stop):
                    rpbp_dict[gid+"\t"+tid][0] = score
            else:
                cds_stop = transcript_dict[tid].genomic_cds[-1].start
                if cds_stop == int(orf_start):
                    rpbp_dict[gid+"\t"+tid][0] = score

with gzip.open(rpbp_result_simulation) as fin:
    csvreader = csv.DictReader(fin,delimiter="\t")
    for row in csvreader:
        tid = row["#id"].split("_")[0]
        gid = transcript_dict[tid].gene_id
        score = row["#bayes_factor_mean"]
        if score in ["inf","-inf"]: continue
        orf_start = row["#start"]
        orf_stop = row["#end"]
        if gid in select_genes:
            if not transcript_dict[tid].cds: continue
            if transcript_dict[tid].genomic_iv.strand == "+":
                cds_stop = transcript_dict[tid].genomic_cds[-1].end
                if cds_stop == int(orf_stop):
                    rpbp_dict[gid+"\t"+tid][1] = score
            else:
                cds_stop = transcript_dict[tid].genomic_cds[-1].start
                if cds_stop == int(orf_start):
                    rpbp_dict[gid+"\t"+tid][1] = score

#riboORF
riboORF_dict = defaultdict(lambda:  ["None","None"])
with open(riboORF_result_true) as fin:
    csvreader = csv.DictReader(fin,delimiter="\t")
    for row in csvreader:
        gid = row["geneID"]
        tid = row["transcriptID"]
        pvalue = row["pvalue"]
        if gid in select_genes:
            if not transcript_dict[tid].cds:
                continue
            else:
                riboORF_dict[gid + "\t" + tid][0] = pvalue
with open(riboORF_result_simulation) as fin:
    csvreader = csv.DictReader(fin,delimiter="\t")
    for row in csvreader:
        gid = row["geneID"]
        tid = row["transcriptID"]
        pvalue = row["pvalue"]
        if gid in select_genes:
            if not transcript_dict[tid].cds:
                continue
            else:
                riboORF_dict[gid + "\t" + tid][1] = pvalue

#orfrate
orfrate_dict = defaultdict(lambda:  ["None","None"])
with open(orfrate_result_true) as fin:
    csvreader = csv.DictReader(fin)
    for row in csvreader:
        tid = row["tid"]
        gid = transcript_dict[tid].gene_id
        score = row["orfrating"]
        orf_tstop = row["tstop"]

        if gid in select_genes:
            if not transcript_dict[tid].cds: continue
            if transcript_dict[tid].cds.end+3 == int(orf_tstop):
                orfrate_dict[gid+"\t"+tid][0] = score
with open(orfrate_result_simulation) as fin:
    csvreader = csv.DictReader(fin)
    for row in csvreader:
        tid = row["tid"]
        gid = transcript_dict[tid].gene_id
        score = row["orfrating"]
        orf_tstop = row["tstop"]

        if gid in select_genes:
            if not transcript_dict[tid].cds: continue
            if transcript_dict[tid].cds.end+3 == int(orf_tstop):
                orfrate_dict[gid+"\t"+tid][1] = score


#generate the ROCdata:
"""
using the intersection gene and transcript of ribocode,rpbp and riboORF
"""
intersect_ids = set(ribocode_dict.keys()).intersection(set(rpbp_dict.keys())).intersection(set(riboORF_dict.keys()))

all = defaultdict(lambda:["None\tNone"]*4)
with open("ribocode_roc.txt","w") as ribocode_f, open("rpbp_roc.txt","w") as rpbp_f, open("orfrate_roc.txt","w") as orfrate_f, open("riboORF_roc.txt","w") as riboORF_f:
    for k in intersect_ids:
        all[k][0] = "\t".join(ribocode_dict[k])
        all[k][1] = "\t".join(rpbp_dict[k])
        all[k][2] = "\t".join(riboORF_dict[k])
        if k in orfrate_dict:
            all[k][3] = "\t".join(orfrate_dict[k])
        ribocode_f.write("%s\t%s\t1\n%s\t%s\t0\n" % (k,ribocode_dict[k][0],k,ribocode_dict[k][1]))
        rpbp_f.write("%s\t%s\t1\n%s\t%s\t0\n" % (k,rpbp_dict[k][0],k,rpbp_dict[k][1]))
        riboORF_f.write("%s\t%s\t1\n%s\t%s\t0\n" % (k,riboORF_dict[k][0],k,riboORF_dict[k][1]))
        if k in orfrate_dict:
            orfrate_f.write("%s\t%s\t1\n%s\t%s\t0\n" % (k,orfrate_dict[k][0],k,orfrate_dict[k][1]))
