#!/usr/bin/env python3

import pandas as pd           ## Import Pandas library for processing dataframe as pd
import numpy as np            ## Import Numy for processing matrix as np
import sys                    ##Import sys library for maximizing the csv limit for large csv file which can be Geneious file.
import argparse


parser = argparse.ArgumentParser(description='Reads_coverage')

parser.add_argument('-S', '--Trimmed_stats', dest='reads_per_sample', type=str,
                        help='Trimmed_stats file')

parser.add_argument('-N', '--cov', dest='sam_depth', type=str,
                        help='coverage file')

args = parser.parse_args()
Cov_file = args.sam_depth
Trim_Stats = args.reads_per_sample

stats_file = pd.read_csv(Trim_Stats, delimiter = '\t',header=None, usecols=[0,1], nrows=2)

Stats = stats_file.T
new_header = Stats.iloc[0] #grab the first row for the header
Stats = Stats[1:] #take the data less the header row
Stats.columns = new_header #set the header row as the df header
#Stats = Stats.reset_index(drop=True)

new = Stats["#File"].str.split("_", n = 1, expand = True)
Stats['Sample_name'] = new[0]
Stats = Stats.rename(columns={'#Total':'Total_Raw_Reads'})

#############
cov_file = pd.read_csv(Cov_file, delimiter = '\t',usecols=[0,3])
# data transpose
cov = cov_file.T

#make first row as column names in pandas
new_header = cov.iloc[0] #grab the first row for the header
cov = cov[1:] #take the data less the header row
cov.columns = new_header #set the header row as the df header


cols = cov.columns

cov['Total_aligned_reads'] = cov[cols].sum(axis=1)
cov = cov.rename(columns={c: c+'_Reads_Aligned' for c in cov.columns if c in cols})


Sample_name = Cov_file.split("/")[-1].split("_")[0]
Sample_out = ("_".join(Cov_file.split("/")[-1].split("_")[0:-1]))
cov['Sample_name'] = Sample_name

cov["Result"] = np.where(cov['Total_aligned_reads'] > 100, 'Pass', 'Fail')
#cov["Result"] = np.where(cov['Sample_name'].str.contains('NTC').all() and cov['Total_aligned_reads'] > 3, 'Fail', 'pass')
# reorder columns
merge_df = pd.merge(Stats,cov, on='Sample_name')

merge_df = merge_df.reindex(sorted(merge_df.columns), axis=1)


merge_df.to_csv(Sample_out+'_readcoverage.csv')
