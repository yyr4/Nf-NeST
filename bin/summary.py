#!/usr/bin/env python3

import pandas as pd           ## Import Pandas library for processing dataframe as pd
import numpy as np            ## Import Numy for processing matrix as np
import sys                    ##Import sys library for maximizing the csv limit for large csv file which can be Geneious file.
import csv                    ##Import csv module to input the csv file to dataframe
import argparse
min_cov = 5
min_VAF = 10

parser = argparse.ArgumentParser(description='merged_csv')
parser.add_argument('-f', dest='merged_csv', type=str, help="snpfilter output")


args = parser.parse_args()
filename=args.merged_csv



csv.field_size_limit(sys.maxsize)                                                                                                         ##Maximize the csv file input size
DF_1=pd.read_csv(filename)


DF_1.to_csv('All_final_snp.csv',index=False)

###### Repportable SNP
Reportable = DF_1.loc[DF_1['SNP_REPORT'] != "Novel SNP"]

Reportable = Reportable.drop(Reportable[(Reportable['Type'] != "No coverage") & (Reportable['AVG_COV'] < min_cov) ].index) 

Reportable.to_csv("Reportable_snps.csv", index=False)


######

Novel = DF_1.loc[DF_1['SNP_REPORT'] == "Novel SNP"]

Novel = Novel.drop(Novel[Novel.AVG_COV < min_cov ].index)
Novel = Novel.drop(Novel[Novel.Confidence < 3 ].index)
#Novel = Novel.drop(Novel[Novel.AVG_VAF < min_VAF ].index)

Novel.to_csv("Novel_snps.csv", sep=',', index=False)
