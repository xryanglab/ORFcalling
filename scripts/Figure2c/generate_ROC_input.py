#!/usr/bin/env python
# -*- coding:UTF-8 -*-
__author__ = 'Zhengtao Xiao'

import pandas as pd
import numpy as np

# read result
def read_resultccds(filename="data/ribotaper_data/results_ccds"):
	df = pd.read_table(filename,index_col=0,na_values='NA')
	return df
def read_ribocode(filename="ribocode_result.txt"):
	df = pd.read_table(filename,index_col=0,na_values='NAN')
	return df

def roc_output(df,outname="ROC_input.txt"):
	with open(outname,"w") as fout:
		fout.write("\t".join(["name","ribotaper","chisq","ORFscore","mymethod","truth"]) + "\n")
		for j in df.index:
			if np.nan not in df.ix[j][["pval_multit_3nt_ribo","chisq_ribo","ORF_score_ribo","my_ribo"]]:
				fout.write(j + "_ribo" + "\t" + "\t".join(map(str,df.ix[j][["pval_multit_3nt_ribo","chisq_ribo",		                                                            "ORF_score_ribo","my_ribo"]])) + "\t1\t\n")
			if np.nan not in df.ix[j][["pval_multit_3nt_rna","chisq_rna","ORF_score_rna","my_rna"]]:
				fout.write(j + "_rna" + "\t" + "\t".join(map(str,df.ix[j][["pval_multit_3nt_rna","chisq_rna",
				                                                                  "ORF_score_rna","my_rna"]])) + "\t0\t\n")


#main
#1. read result
ribotaper_df = read_resultccds("results_ccds")
ribocode = read_ribocode("ribocode_result.txt")

##Drop NA values
ribotaper_df_na = ribotaper_df[(ribotaper_df["pval_multit_3nt_ribo"] != np.nan) | (ribotaper_df["pval_multit_3nt_rna"] != np.nan)]
ribocode_df_na = ribocode[(ribocode["my_ribo"] != np.nan) | (ribocode['my_rna'] != np.nan)]

#2. merger
merger_df = ribotaper_df_na.merge(ribocode_df_na,how="inner",left_index=True,right_index=True)

#3. filter the psites_sum < 5:
merger_df_filter = merger_df[merger_df.P_sites_sum > 5]

#4. ROC results
roc_output(merger_df_filter)



