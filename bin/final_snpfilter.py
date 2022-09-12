#!/usr/bin/env python3

import pandas as pd
import numpy as np
import argparse
import sys

parser = argparse.ArgumentParser(description='snpfilter')

parser.add_argument('-C', '--WT_coverage', dest='WT_coverage', type=str,
                        help='Wildtype snps coverage file' )
parser.add_argument('-A', '--VCF_merge', dest='vcf_merge', type=str,
                        help='final merge vcf file')
parser.add_argument('-V', '--VOI_file', dest='VOI', type=str,
                        help='Variants of interest file')
args = parser.parse_args()
WT_cov = args.WT_coverage
vcf_merge = args.vcf_merge
VOI_file = args.VOI

Sample_out = ("_".join(vcf_merge.split("/")[-1].split("_")[0:3]))

# reading two csv files
# if coverage .csv file is empty, skip the code and exit
try:
    data1 = pd.read_csv(WT_cov)
except pd.errors.EmptyDataError:
    sys.exit()

# read the merged vcf file
data2 = pd.read_csv(vcf_merge)

# Replace DHPS_437Corrected to DHPS and rename chrom to CHROM
data1 = data1.replace({'DHPS_437Corrected': 'DHPS'})
data2 = data2.replace({'DHPS_437Corrected': 'DHPS'})
data1 = data1.rename({'Chrom': 'CHROM'}, axis=1)

# convert avg variant frequancy to percentage
data2["AVG_VAF"] = (100 * data2["AVG_VAF"]).round(1)
data2["VOI"] = data2["AA_change"]

# drop the rows where avg_vcf is less than 1%
data2 = data2[data2.AVG_VAF >= 1]
data2["AVG_VAF"] = data2["AVG_VAF"].astype(str) + '%'

# merge two datafrme (snp output and WT coverage out)
output1 = pd.concat([data2, data1])

# Creat a new column for mutation and wildtype
output1['Type'] = np.where(output1['AVG_VAF'].notnull(), "Mutation" , "WildType")

# covert coverge into integer
output1 = output1.astype({"AVG_COV":'int'})


##Drop duplicates meaning if the values are already in variants then drop it from the wildtypes
output1=output1.drop_duplicates(subset =["Sample_name", "VOI"] )

# Drop the  index column Unnamed and reset the index
output1.drop(output1.filter(regex="Unname"),axis=1, inplace=True)
output1.reset_index(drop=True, inplace=True)

#Drop two unwanted columns and replace empty cells with nan value
#output1.drop(['AA.pos', 'triplet'], axis=1, inplace=True)
output1.replace(r'^\s*$', np.nan, regex=True)


#Filter the snps based on VOI files creat a list
df_voi=pd.read_csv(VOI_file)
df_voi["VOI"]=df_voi["RefAA"]+df_voi["AAPos"].astype(str)+df_voi["AltAA"]
df_voi = df_voi.rename({'Gene':'CHROM'},axis=1)
voi_list=df_voi["VOI"].tolist()

#create a column snp_report if snps is present in voi file list, call it reportable snp otherwise novel
output1["SNP_Report"] = np.where(output1["VOI"].isin(voi_list), "Reportable" , "Novel")

#Reorder the columns
output1 = output1[['Sample_name', 'CHROM', 'POS', 'AA_change', 'AVG_VAF', 'AVG_COV', 'REF', 'ALT','VARTYPE', 'Annotation', 'VarCal'
                     ,'Confidence', 'VOI', 'Type', 'SNP_Report'  ]]
# final output in csv
output1.to_csv(Sample_out+'_final_snp.csv')

# separte reportable and novel mutations

Reportable = output1[output1['SNP_Report'] == 'Reportable']
Novel = output1[output1['SNP_Report'] == 'Novel']

Reportable.to_csv(Sample_out+'_Reportable.csv', index=False)
Novel.to_csv(Sample_out+'_Novel.csv', index=False)
