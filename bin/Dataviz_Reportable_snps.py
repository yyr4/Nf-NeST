#!/usr/bin/env python3
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import sys
import argparse


parser = argparse.ArgumentParser(description='filename')
parser.add_argument('-r', dest='filename', type=str, help="name of snp summary file")

args = parser.parse_args()
filename=args.filename


try:
    DF = pd.read_csv(filename)

except pd.errors.EmptyDataError:
    DF = pd.DataFrame()
    print("Reportable SNPs are empty")
    sys.exit()

##########  Reportable_Per_SNP_depth
df1 = DF[['CHROM', 'VOI', 'AVG_COV']]

df1["variants"] = df1[['CHROM',"VOI"]].apply(":".join, axis=1)
df1 = df1.sort_values('variants')

#sns.set_style("whitegrid")


fig, ax = plt.subplots()
fig.set_size_inches(25, 25)
sns.boxplot(y = 'variants', x = 'AVG_COV', data = df1, hue="CHROM", dodge=False)
sns.despine()

fig.savefig('Reportable_Per_SNP_depth.pdf')

####################### Bar plot for Reportable_snps

df = DF.groupby(['CHROM','VOI','Type','SNP_Report']).size().reset_index(name='counts')

df_pv= df.pivot_table(values='counts', index=['CHROM','VOI',"SNP_Report"], columns='Type', aggfunc='first')
df_pv = df_pv.fillna(0).reset_index()

column_names = ['Mixed','Mutant','Wildtype']
df_pv['Total']= df_pv[column_names].sum(axis=1)
df_pv["Snps"] = df_pv["CHROM"] + ":" + df_pv["VOI"] + ":N=" + df_pv["Total"].astype(str)

df_pv = df_pv.rename(columns={'Mixed': 'Minor', 'Mutant':'Major'})
SNPvals=df_pv[["Snps",'Minor','Major','Wildtype','Total','SNP_Report']]

# from raw value to ratios

#Setup for loading
Totes = SNPvals.groupby('Snps')['Total'].sum().reset_index()
Minor = SNPvals.groupby('Snps')['Minor'].sum().reset_index()
Major = SNPvals.groupby('Snps')['Major'].sum().reset_index()
WT = SNPvals.groupby('Snps')['Wildtype'].sum().reset_index()

#Math and definition of SNPratio
Minor['SNPratio'] = [i / j for i,j in zip(Minor['Minor'], Totes['Total'])]
Major['SNPratio'] = [i / j for i,j in zip(Major['Major'], Totes['Total'])]
WT['SNPratio'] = [i / j  for i,j in zip(WT['Wildtype'], Totes['Total'])]

AllTogether = pd.concat([Minor.Snps, Minor.SNPratio, Major.SNPratio, WT.SNPratio], axis=1)
#AllTogether.columns = ['SNPs','Minor: AF < 50%','Major: AF >= 50%','WildType: AF=0%']
AllTogether.columns = ['Snps','Minor: AF < 95%', 'Major: AF >= 95%', 'WildType: AF=0%']
#AllTogether.to_csv("Tab_Table_snps.csv", index=False)

df_table_SNP=AllTogether.sort_values(by=['Snps'])

df_table_SNP["index"]=df_table_SNP.Snps.str.split(":").str[1].str[1:-1]
df_table_SNP["index"]=df_table_SNP["index"].astype(int)
df_table_SNP["index2"]=df_table_SNP.Snps.str.split(":").str[0]



plot = df_table_SNP.sort_values(by = ['index2', 'index'],ascending=False)[['Snps',  'Minor: AF < 95%','Major: AF >= 95%',  'WildType: AF=0%']].plot(x='Snps', kind='barh', stacked=True, title='Drug Resistance Reportable SNPs', figsize=(20,20), color={"Minor: AF < 95%": "#F3ABA8", "Major: AF >= 95%": "#98DAA7","WildType: AF=0%": "#5975A4"})
#plot.legend(bbox_to_anchor=(0.97, 0.1))

plot.legend(ncol = 2, loc = 'lower right')
#sns.despine(left = True, bottom = True)
plot.set(ylabel="SNPs")
plot.set(xlabel="SNP ratio")
plot.legend(loc=(1.04,0))
plt.savefig('SNPs-Reportable.pdf')
