#!/usr/bin/env python3

import pandas as pd
import numpy as np
import argparse
import sys
min_cov = 5
min_VAF = 5

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

Sample_out = ("_".join(vcf_merge.split("/")[-1].split("_")[0:-1]))

# reading two csv files
# if coverage .csv file is empty, skip the code and exit
try:
    data1 = pd.read_csv(WT_cov)
except pd.errors.EmptyDataError:
    sys.exit()

# read the merged vcf file

try:
    data2 = pd.read_csv(vcf_merge)
except pd.errors.EmptyDataError:
    sys.exit()

# rename chrom to CHROM

data1 = data1.rename({'Chrom': 'CHROM'}, axis=1)

# convert avg variant frequancy to percentage
data2["AVG_VAF"] = (100 * data2["AVG_VAF"]).round(1)
data2["VOI"] = data2["AA_change"]

# drop the rows where avg_vcf is less than 1%
data2 = data2[data2.AVG_VAF > 5]

# merge two datafrme (snp output and WT coverage out)
output1 = pd.concat([data2, data1])

output1 = output1[~output1["Annotation"].str.contains("stop", na=False)]
output1 = output1[~output1["Annotation"].str.contains("splice", na=False)]


#output1 = output1.astype({"AVG_COV":'int'})
output1['AVG_VAF'] = output1['AVG_VAF'].fillna(0)



output1[['AVG_COV']] = output1[['AVG_COV']].astype(float).round(0)


#output1['AVG_VAF'] = output1['AVG_VAF'].str.replace(r'%', '')
output1['AVG_VAF'] = output1['AVG_VAF'].astype(float)



# Creat a new column for mutation and wildtype
#output1['Type'] = np.where(output1['AVG_VAF'].notnull(), "Mutation" , "WildType")
def snp_type(output1):

    if (output1["AVG_COV"] <= 2):
        return "No coverage"
    elif (output1["AVG_COV"] > 2) and (output1['AVG_VAF'] >= 95):
        return 'Mutant'
    elif (output1["AVG_COV"] > 2) and (1 < output1['AVG_VAF'] < 95):
        return 'Mixed'
    elif (output1["AVG_COV"] > 2) and (output1['AVG_VAF'] == 0.0):
        return "Wildtype"
    else:
        return None

output1["Type"] = output1.apply(snp_type, axis = 1)

output1['AVG_VAF'] = output1['AVG_VAF'].astype(str) + '%'
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

df_voi["AA_change"]=df_voi["RefAA"]+df_voi["AAPos"].astype(str)
df_voi = df_voi.rename({'Chr':'CHROM'},axis=1)
voi_list=df_voi["VOI"].tolist()


#create a column snp_report if snps is present in voi file list, call it reportable snp otherwise novel

def snp_report(output1):

    if (output1["AVG_COV"] != 0) and (output1["VOI"] in(voi_list)):
        return 'Reportable SNP'
    elif (output1["AVG_COV"] != 0) and (output1["VOI"] not in (voi_list)):
        return 'Novel SNP'

output1["SNP_Status"] = output1.apply(snp_report, axis = 1)


merge_df = pd.merge(output1, df_voi, on=['CHROM','VOI'], how="outer")

merge_df['SNP_REPORT'] = np.where(merge_df['snp_Report'].notnull(), merge_df['snp_Report'], merge_df['SNP_Status'])

merge_df = merge_df.rename({'AA_change_x': 'AA_change'}, axis=1)
merge_df['CHROM'] = merge_df['CHROM'].replace({'NC_009906.1': 'Pvcrt', 'NC_009910.1':'Pvdhfr', 'NC_009915.1': 'Pvmdr1', 'NC_009919.1': 'Pvdhps'})
#Reorder the columns
merge_df = merge_df[['Sample_name', 'CHROM', 'POS', 'AA_change', 'AVG_VAF', 'AVG_COV', 'REF', 'ALT','VARTYPE', 'Annotation', 'VarCal'
                    ,'Confidence', 'VOI', 'Type', 'SNP_REPORT']]

# final output in csv
merge_df.to_csv(Sample_out+'_final_snp.csv', index=False)
