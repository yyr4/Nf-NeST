#!/usr/bin/env python3

 # import the dependencies

import pandas as pd
import os
import glob
import numpy as np
import re
import argparse
import sys
import warnings

'''
This script process each vcf files Separately from each sample, extract the metadata  and annotation using the pandas dataframe.
The output will be csv file with
'''


#
parser = argparse.ArgumentParser(description='filename')
parser.add_argument('-n', dest='filename', type=str, help="name of vcf file")

args = parser.parse_args()
filename=args.filename


header = ''
informations = ''


##  If vcf file is not empty, split the file based on header and varinats information by tab Separated
##  If vcf file is empty, create an empty file and add header.

with open(filename, "r") as f:
            lines = f.readlines()

            if len(lines) > 0:
                chrom_index = [i for i, line in enumerate(lines) if line.strip().startswith("#CHROM")]
                data = lines[chrom_index[0]:]
                header = (data[0].strip().split("\t"))
                informations = [d.strip().split("\t") for d in data[1:]]

            else:
                file = ("_".join(filename.split("/")[-1].split("_vartype.vcf")[0:-1]))+'.csv'
                fb = open(file, "a")
                fb.write(",#CHROM,POS,REF,ALT,QUAL,Sample_name,Source,AD,DP,AF,DP4,VARTYPE,Annotation,Annotation_Impact,AA_change,VAF")
                fb.close()
                sys.exit()

## get source info from vcf file
for l in lines:
    if l.strip().startswith("##source"):
        #print(l)
        Source = l.strip().split("=")[1].split()[0]
        break
    else:
        Source = "Samtools"

## get sample name from filename
Sample_name = filename.split("/")[-1].split("_")[0]
Sample_out = ("_".join(filename.split("/")[-1].split("_vartype.vcf")[0:-1]))


vcf = pd.DataFrame(informations, columns=header)                            # create a dataframe called vcf with information and header
vcf.columns = [*vcf.columns[:-1], 'Genotype']                               # rename last column name to Genotype

## filter dataframe based on vartype

df_col = vcf.loc[:,['#CHROM', 'POS', 'REF', 'ALT', 'QUAL', 'INFO', 'FORMAT', 'Genotype']]    # Only select the columns of Interest exclude column ID and Filter
vartype = ["VARTYPE=SNP"]
df_rows = df_col[df_col['INFO'].str.contains('|'.join(vartype))]            # filter datafrma based on vartype = snp, MNP
df_rows = df_rows.reset_index(drop=True)                                    # reindex dataframe
df_rows['Sample_name'] = filename.split("/")[-1].split("_")[0]  # add column name sample name
df_rows['Source'] = Source            # add column name as Source

## Check for Empty File
## Check if vcf file does not contains vartype snp/MNP, then create an empty file with only header

if df_rows.empty:
        file = ("_".join(filename.split("/")[-1].split("_vartype.vcf")[0:-1]))+'.csv'        # Create an empty file
        fb = open(file, "a")
        fb.write(",#CHROM,POS,REF,ALT,QUAL,Sample_name,Source,AD,DP,AF,DP4,VARTYPE,Annotation,Annotation_Impact,AA_change,VAF")     # Add header
        fb.close()
        warnings.filterwarnings("ignore")
        sys.exit()

## process genotype data to get AD (Number of observation for each allele)

Genotype_1 =  df_rows['Genotype'].str.split(':', expand=True)
Genotype_1.columns = df_rows['FORMAT'].str.split(':', expand=True).iloc[0]

## Merge two dataframe
df_rows  = pd.merge(df_rows, Genotype_1, left_index=True, right_index=True)
df_rows["AD"] =np.nan if 'AD' not in df_rows.columns else df_rows['AD']

## filter the  INFO column to get the  annotation data
##'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length |

dict_info = {} # Create an empthy dictionary for 'INFO' data
dict_ANN = {}  # Create an empthy dictionary for Snpeff annotation  data
dict_info_merge = []

# Separate InFO by ';' and save in dictionary dict_info
for index, row in df_rows.iterrows():
    info_array = row["INFO"]
    for i in info_array.split(";"):
        if "=" in i:
            key = i.split("=")[0]
            value = i.split("=")[1]
            dict_info[key]=value

    # create a dictionary for annotation info
    my_Ann = dict_info["ANN"]
    keys = [ 'Allele' , 'Annotation' , 'Annotation_Impact',  'Gene_Name', 'Gene_ID' , 'Feature_Type' , 'Feature_ID' , 'Transcript_BioType' , 'Rank' , 'HGVS.c' , 'HGVS.p' , 'cDNA.pos/cDNA.length' , 'CDS.pos/CDS.length' , 'AA.pos/AA.length']
    values = my_Ann.split("|")[0:14]

    for i in range(len(keys)):
        dict_ANN[keys[i]] = values[i]

    # copy the content of dict_info to new dict so that we can remove annotation and merge the two dictonary
    New_info = dict_info
    New_info.pop("ANN")

    dict_merge = {**New_info, **dict_ANN}
    dict_info_merge.append(dict_merge)

# Create a new dataframe from Previous two merge list
df1 = pd.DataFrame(dict_info_merge)
# check if column exits or not: put NAN if column not exits
df1["POS"] = df_rows['POS']
df1["AF"] =np.nan if 'AF' not in df1.columns else df1['AF']
df1["DP4"] =np.nan if 'DP4' not in df1.columns else df1['DP4']


df1['Gene_Name'] = df1['Gene_Name'].replace(['DHFR-TS','DHPS','CYTB'],['DHFR','DHPS_437Corrected','mitochondrial_genome_-_CYTB_CDS'])
df1 = df1.rename({'Gene_Name': '#CHROM'}, axis=1)
df1['AA_change'] = [i[2:] for i in df1['HGVS.p']]


ANN_info = df1.loc[:,['DP',"AF",'DP4', 'VARTYPE','Annotation', 'Annotation_Impact', '#CHROM',
                      'POS','AA_change'] ]

VCF_info = df_rows.loc[:,['#CHROM', 'POS', 'REF', 'ALT', 'QUAL','Sample_name', 'Source','AD']]

result = pd.merge(VCF_info, ANN_info, on=["#CHROM","POS"])


VAF = []
try:
    for i in result['AD'].str.split(","):
        alt =  int(i[1])
        total = int(i[0])+ int(i[1])
        alfreq = round(alt/float(total),2)
        VAF.append(alfreq)

except:
    for i in result['DP4'].str.split(","):
        alt = int(i[2])+ int(i[3])
        total = int(i[0])+ int(i[1])+int(i[2])+ int(i[3])
        alfreq = round(alt/float(total),2)
        VAF.append(alfreq)

result["VAF"] = VAF
result = result.replace(r'^\s*$', np.nan, regex=True)

## Change the value of Pfcrt 75 codon to correcly annotate

# If #chrom = pfcrt, POS = 495 and AA_change = N75D; Then change the REF = AAT, ALT = GAT and AA_change= N75E

df_1 = result.loc[(result['#CHROM']=="PfCRT")&(result['POS']==495)|(result['AA_change']=='N75D')]

#df_1 = result[result['AA_change'].str.contains('N75D')]
df_1.replace({'AA_change':{'N75D' :'N75E'}, 'REF': {'A':'AAT'}, 'ALT': {'G':'GAT'}},inplace= True)


# If  #chrom = pfcrt, POS = 497 and AA_change = N75K : Then remove this row
result = result[(result.POS != 497) & (result.AA_change != 'N75K')  ]
result = result[(result.POS != 495) & (result.AA_change != 'N75D')  ]
result = result.reset_index(drop=True)


# COncat two df to get final dataframe and sort the POS value based on chrom and POS
dat1 = pd.concat([result, df_1], axis=0)
dat1.sort_values(by=['#CHROM','POS'],inplace=True)
dat1.to_csv(Sample_out+'.csv')
