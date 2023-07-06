#!/usr/bin/env python3

import pandas as pd           ## Import Pandas library for processing dataframe as pd
import numpy as np            ## Import Numy for processing matrix as np
import sys                    ##Import sys library for maximizing the csv limit for large csv file which can be Geneious file.
import csv                    ##Import csv module to input the csv file to dataframe
import argparse


parser = argparse.ArgumentParser(description='merged_csv')
parser.add_argument('-f', dest='merged_csv', type=str, help="snpfilter output")


args = parser.parse_args()
filename=args.merged_csv



csv.field_size_limit(sys.maxsize)                                                                                                         ##Maximize the csv file input size
DF_1=pd.read_csv(filename)

def TreatmentDay(row): ##Set up a function for assigning TreatmentDay based on the values in the Sample_name column
    if row['Sample_name'][6:8]=="00":
          return '0'
    elif row['Sample_name'][6:8]!="00":
          return row['Sample_name'][6:8]

def POOLED(row):                                            ##Set up a function for POOLED based on the values in the Sample_name column
    if row['Sample_name'].replace(" ","")[12:13]=="p":      ##If position 12 has p in it
        return 'POOLED'                             ##then it is considered POOLED sample
    elif row['Sample_name'].replace(" ","")[12:13]!="p":    ##If position 12 has no p in it
        return 'individual'                         ##then it is considered individual sample

def year(row):  ##Set up a function for Year  based on the values in the Sample_name column
    return row['Sample_name'][0:2]

DF_1[['TREATMENT_DAY', "POOLED","YEAR"]] = DF_1.apply([TreatmentDay,POOLED,year], axis=1)

DF_1[['AVG_COV']] = DF_1[['AVG_COV']].astype(float)

DF_1['AVG_VAF'] = DF_1['AVG_VAF'].str.replace(r'%', '')
DF_1['AVG_VAF'] = DF_1['AVG_VAF'].astype('float')
DF_1 = DF_1.drop(DF_1[DF_1.AVG_VAF < 5 ].index)


# Remove rows which are stop gain/stop lost Annotation
DF_1 = DF_1[~DF_1["Annotation"].str.contains("stop", na=False)]
DF_1 = DF_1[~DF_1["Annotation"].str.contains("splice", na=False)]

# Assign coverage 0 as No coverage or No Amplification
DF_1.loc[(DF_1['AVG_COV']== 0,"Type")]  = 'No Coverage'

# fill NA as 0
DF_1['AVG_VAF'] = DF_1['AVG_VAF'].fillna(0)

# Assign mutatnt aixed and WT based on Varinat frequancy
DF_1.loc[DF_1['AVG_VAF'] < 95, 'Type'] = 'Mixed'
DF_1.loc[DF_1['AVG_VAF'] >= 95, 'Type'] = 'Mutant'
DF_1.loc[DF_1['AVG_VAF'] == 0, 'Type'] = 'Wildtype'


DF_1['AVG_VAF'] = DF_1['AVG_VAF'].astype(str) + '%'


DF_1.to_csv('All_final_snp.csv',index=False)

###### Repportable SNP
Reportable = DF_1.loc[DF_1['SNP_Report'] == "Reportable SNP"]

Reportable.drop(Reportable[Reportable.AVG_COV < 2 ].index)
Reportable.to_csv("Reportable_snps.csv", index=False)

######

Novel = DF_1.loc[DF_1['SNP_Report'] == "Novel SNP"]

Novel = Novel.drop(Novel[Novel.AVG_COV < 5 ].index)
Novel = Novel.drop(Novel[Novel.Confidence < 2 ].index)

Novel.to_csv("Novel_snps.csv", sep=',', index=False)
