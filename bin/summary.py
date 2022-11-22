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
df=pd.read_csv(filename)

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

df[['TREATMENT_DAY', "POOLED","YEAR"]] = df.apply([TreatmentDay,POOLED,year], axis=1)


df.to_csv('All_final_snp.csv',index=False)
