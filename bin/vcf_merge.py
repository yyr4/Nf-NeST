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


#
parser = argparse.ArgumentParser(description='filename')
parser.add_argument('-n', dest='filename', type=str, help="name of vcf file")

args = parser.parse_args()
filename=args.filename


## get sample name from filename
Sample_name = filename.split("/")[-1].split("_")[0]
Sample_out = ("_".join(filename.split("/")[-1].split("_vartype.vcf")[0:-1]))


header = ''
informations = ''


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
              fb.write(",#CHROM,POS,REF,ALT,QUAL,Sample_name,Source,AD,DP,AF,DP4,VARTYPE,Annotation,Annotation_Impact,Rank,Gene_ID,HGVS.c,HGVS.p,CDS.length,AA.length,AA_change,VAF")
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

vcf = pd.DataFrame(informations, columns=header)
vcf = vcf.set_axis([*vcf.columns[:-1], 'Genotype'], axis=1, inplace=False)  # rename last column name to Genotype

## filter dataframe based on vartype.
df_col = vcf.loc[:,['#CHROM', 'POS', 'REF', 'ALT', 'QUAL', 'INFO', 'FORMAT', 'Genotype']]

#vartype = ["VARTYPE=SNP", "VARTYPE=DEL", 'VARTYPE=INS', 'VARTYPE=MNP']
vartype = ["VARTYPE=SNP", 'VARTYPE=MNP']
df_rows = df_col[df_col['INFO'].str.contains('|'.join(vartype))]   ## filter datafrma based on vartype = snp
df_rows = df_rows.reset_index(drop=True)
df_rows['Sample_name'] = Sample_name  # add column name sample name and source
df_rows['Source'] = Source

# check if vcf file does not contains vartype snp/MNP, then creat a empty file with header

if df_rows.empty:
        file = ("_".join(filename.split("/")[-1].split("_vartype.vcf")[0:-1]))+'.csv'
        fb = open(file, "a")
        fb.write(",#CHROM,POS,REF,ALT,QUAL,Sample_name,Source,AD,DP,AF,DP4,VARTYPE,Annotation,Annotation_Impact,Rank,Gene_ID,HGVS.c,HGVS.p,CDS.length,AA.pos,AA.length,AA_change,VAF")
        fb.close()
        warnings.filterwarnings("ignore")
        sys.exit()

## process genotype data to get AD
Genotype_1 =  df_rows['Genotype'].str.split(':', expand=True)
Genotype_1.columns = df_rows['FORMAT'].str.split(':', expand=True).iloc[0]

df_rows  = pd.merge(df_rows, Genotype_1, left_index=True, right_index=True)
df_rows["AD"] =np.nan if 'AD' not in df_rows.columns else df_rows['AD']

## filert INFO colm based on annotation
##'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length |

dict_info = {}
dict_ANN = {}
dict_info_merge = []


for index, row in df_rows.iterrows():
    #print (index, row)

    info_array = row["INFO"]
    for i in info_array.split(";"):
        if "=" in i:
            key = i.split("=")[0]
            value = i.split("=")[1]
            dict_info[key]=value
    my_Ann = dict_info["ANN"]

# creat a dcitonary for annotation info
    keys = [ 'Allele' , 'Annotation' , 'Annotation_Impact',  'Gene_Name', 'Gene_ID' , 'Feature_Type' , 'Feature_ID' , 'Transcript_BioType' , 'Rank' , 'HGVS.c' , 'HGVS.p' , 'cDNA.pos/cDNA.length' , 'CDS.pos/CDS.length' , 'AA.pos/AA.length']
    values = my_Ann.split("|")[0:14]

    for i in range(len(keys)):
        dict_ANN[keys[i]] = values[i]

## grap nucleoptide position from ANN
    val = my_Ann.split("|")

    r = re.compile("^c.[0-9]*[A-Z]{1}>[A-Z]{1}")
    newlist = list(filter(r.match,val))

    POS_ANN = [pos[2:-3] for pos in newlist]
    POS_ANN = (", ".join(POS_ANN))
    dict_ANN["POS_ANN"] = POS_ANN


# copy the content of dict_info to new dict so that we can remove annotation and merge the two dictonary
    New_info = dict_info
    New_info.pop("ANN")

    dict_merge = {**New_info, **dict_ANN}
    dict_info_merge.append(dict_merge)



df1 = pd.DataFrame(dict_info_merge)

# check if column exits or not: put NAN if column not exits
df1["POS"] = df_rows['POS']
df1["AF"] =np.nan if 'AF' not in df1.columns else df1['AF']
df1["DP4"] =np.nan if 'DP4' not in df1.columns else df1['DP4']


df1['Gene_Name'] = df1['Gene_Name'].replace(['DHFR-TS','DHPS','CYTB'],['DHFR','DHPS_437Corrected','mitochondrial_genome_-_CYTB_CDS'])

df1 = df1.rename({'Gene_Name': '#CHROM'}, axis=1)

df1['AA_change'] = [i[2:] for i in df1['HGVS.p']]

ANN_info = df1.loc[:,['DP',"AF",'DP4', 'VARTYPE','Annotation', 'Annotation_Impact', '#CHROM', 'Gene_ID',
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
print(result)

result.to_csv(Sample_out+'.csv')
