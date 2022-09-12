#!/usr/bin/env python3

import pandas as pd
import glob
import os
import argparse
import numpy as np



parser = argparse.ArgumentParser(description='Merge_vcf')

parser.add_argument('-v1', '--vcf1', dest='vcf_path1', type=str,
                        help='Sample name',nargs='?')
parser.add_argument('-v2', '--vcf2', dest='vcf_path2', type=str,
                        help='Sample name',nargs='?', const=1, default=1)
parser.add_argument('-v3', '--vcf3', dest='vcf_path3', type=str,
                        help='Sample name',nargs='?', const=1, default=1)

args = parser.parse_args()
files =args.vcf_path1, args.vcf_path2, args.vcf_path3


samp_files = list(files)


sample_name = []
for i in samp_files:
    i = "_".join(i.split("/")[-1].split("_")[0:3])
    if i not in sample_name:
        sample_name.append(i)
Sample_name = sample_name[0]


# joining files with concat and read_csv
Merge_df = pd.concat(map(pd.read_csv, files), ignore_index=True)
Merge_df = Merge_df.reset_index(drop=True)

df = Merge_df.replace(r'^\s*$', np.nan, regex=True)

####################

Merge_df['VarCal'] = Merge_df.groupby(["#CHROM","POS",'AA_change'])['Source'].transform(
                                              lambda x: ', '.join(x))
# drop duplicate data
Merge_df = Merge_df.drop_duplicates()
################
# add confidence and Avg VAF column
Merge_df_2 = Merge_df.groupby(["#CHROM","POS",'AA_change']).agg({'AA_change': 'count', 'VAF':'mean','DP' : 'mean'})
Merge_df_2.columns = ['Confidence', 'AVG_VAF', 'AVG_COV']
Merge_df_2 = Merge_df_2.reset_index()
#drop the
Merge_df = Merge_df.drop(['Source','QUAL','Unnamed: 0','AD','DP','AF','DP4','VAF','Annotation_Impact','Gene_ID'], axis=1)
Merge_df = Merge_df.drop_duplicates()

###################

result_1 = pd.merge(Merge_df_2, Merge_df,how='left',  on=["#CHROM","POS",'AA_change'])
result_1 = result_1.rename({'#CHROM' : 'CHROM'}, axis=1)


result_1.to_csv(Sample_name+"_merge"+'.csv')
