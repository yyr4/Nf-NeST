#!/usr/bin/env python3

import pandas as pd
import glob
import os
import argparse
import numpy as np
import sys



parser = argparse.ArgumentParser(description='Merge_vcf')

parser.add_argument('-v1', '--vcf1', dest='vcf_path1', type=str,
                        help='Sample name',nargs='?')
parser.add_argument('-v2', '--vcf2', dest='vcf_path2', type=str,
                        help='Sample name',nargs='?', const=1, default=1)
parser.add_argument('-v3', '--vcf3', dest='vcf_path3', type=str,
                        help='Sample name',nargs='?', const=1, default=1)
parser.add_argument('-v4', '--vcf4', dest='vcf_path4', type=str,
                        help='Sample name',nargs='?', const=1, default=1)

args = parser.parse_args()
files =args.vcf_path1, args.vcf_path2, args.vcf_path3, args.vcf_path4

samp_files = list(files)
sample_name = []
for i in samp_files:
    i = "_".join(i.split("/")[-1].split("_")[0:-1])
    if i not in sample_name:
        sample_name.append(i)
#print(sample_name)
Sample_name = sample_name[0]

# joining files with concat and read_csv
Merge_df = pd.concat(map(pd.read_csv, files), ignore_index=True)
Merge_df = Merge_df.reset_index(drop=True)

df = Merge_df.replace(r'^\s*$', np.nan, regex=True)
#df.to_csv('file_merge.csv')
#print(df.columns)

Merge_df['Source'].str.strip()
####################

Merge_df['VarCal'] = Merge_df.groupby(["#CHROM","POS",'AA_change'])['Source'].transform(
                                              lambda x: ','.join(x))
# drop duplicate data
Merge_df = Merge_df.drop_duplicates()


################
# add confidence and Avg VAF column
Merge_df_2 = Merge_df.groupby(["#CHROM","POS",'AA_change']).agg({'AA_change': 'count', 'VAF':'mean','DP' : 'mean'})
Merge_df_2.columns = ['Confidence', 'AVG_VAF', 'AVG_COV']
Merge_df_2 = Merge_df_2.reset_index()
#drop the
Merge_df_3 = Merge_df.drop(['Source','QUAL','Unnamed: 0','AD','DP','AF','DP4','VAF','Annotation_Impact'], axis=1)
Merge_df_3 = Merge_df_3.drop_duplicates()



result_1 = pd.merge(Merge_df_2, Merge_df_3,how='left',  on=["#CHROM","POS",'AA_change'])
result_1 = result_1.rename({'#CHROM' : 'CHROM'}, axis=1)


################################
# drop duplicates while caling mnps and snps
#add new confidence

result_1.AVG_VAF = result_1.AVG_VAF.round(1)

df1 = result_1[result_1.duplicated(subset=['Sample_name','CHROM','POS','AVG_VAF'],keep=False)]

df_snps = df1.loc[df1["VARTYPE"] == 'SNP']


res = pd.merge(result_1,df_snps, indicator=True, how='outer').query('_merge=="left_only"').drop('_merge', axis=1)

res['codon'] = res.AA_change.str.extract('(\d+)').astype(str)

df2 = res[res.duplicated(subset=['Sample_name','CHROM','codon'], keep=False)]


## add new confiodence and varcal
df3 = df2.groupby(['Sample_name','codon'])['VarCal'].agg(','.join).reset_index()

df3['VarCal'] = df3['VarCal'].astype(str).str.split(',').apply(set).str.join(',')

df3['Confidence_2'] =  df3["VarCal"].apply(lambda x: len(x.split(',')))
df_merge = pd.merge(df2, df3, on=['Sample_name','codon'])

df_merge = df_merge.drop(columns=['VarCal_x','Confidence'])
df_merge = df_merge.rename(columns = {'VarCal_y': 'VarCal', 'Confidence_2':'Confidence' })


#df2_snps = df2.loc[df2["VARTYPE"].str.contains('SNP') ]
df2_snps = df2.loc[df2["VARTYPE"] =='SNP' ]

res_1 = pd.merge(res, df2_snps, indicator=True, how='outer').query('_merge=="left_only"').drop('_merge', axis=1)



res_1.set_index(['HGVS.c','POS', 'ALT'],inplace=True)
res_1.update(df_merge.set_index(['HGVS.c','POS','ALT']))
res_1.reset_index()

res_1.to_csv(Sample_name+"_merge"+'.csv')


###########Introns file###########


Introns = Merge_df[Merge_df['Annotation'] == "intron_variant"]

if Introns.empty:
        file = (Sample_name+"_introns"+'.csv')        # Create an empty file
        fb = open(file, "a")
        fb.write(",#CHROM,POS,HGVS.c,Confidence,AVG_VAF,AVG_COV,REF,ALT,Sample_name,VARTYPE,Annotation,VarCal")     # Add header
        fb.close()
        #warnings.filterwarnings("ignore")
        sys.exit()


Introns['VarCal'] = Introns.groupby(["#CHROM","POS",'HGVS.c'])['Source'].transform(lambda x: ', '.join(x))

Introns =Introns.drop_duplicates()
Introns_1 = Introns.groupby(["#CHROM","POS",'HGVS.c']).agg({'HGVS.c': 'count', 'VAF':'mean','DP' : 'mean'})
Introns_1.columns = ['Confidence', 'AVG_VAF', 'AVG_COV']
Introns_1 = Introns_1.reset_index()
Introns = Introns.drop(['Source','QUAL','Unnamed: 0','AD','DP','AF','DP4','VAF','Annotation_Impact','AA_change','HGVS.p'], axis=1)
Introns = Introns.drop_duplicates()

result_2 = pd.merge(Introns_1, Introns ,how='left',  on=["#CHROM","POS",'HGVS.c'])
result_2 = result_2[result_2.Confidence > 1]
#result_2 = result_2.reset_index()
result_2.to_csv(Sample_name+"_introns"+'.csv')
