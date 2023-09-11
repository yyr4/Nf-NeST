#!/usr/bin/env python3
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import sys
import argparse


parser = argparse.ArgumentParser(description='filename')
parser.add_argument('-n', dest='filename', type=str, help="name of snp summary file")

args = parser.parse_args()
filename=args.filename

try:
    Novel = pd.read_csv(filename)
except pd.errors.EmptyDataError:
    print('Empty Novel snps csv file!')
    sys.exit

df = Novel.groupby(['CHROM','VOI','Type','Annotation']).size().reset_index(name='counts')

df_pv= df.pivot_table(values='counts', index=['CHROM','VOI','Annotation'], columns='Type', aggfunc='first')
df_pv = df_pv.fillna(0).reset_index()

column_names = ['Mixed','Mutant']
df_pv['Total']= df_pv[column_names].sum(axis=1)
df_pv["Snps"] = df_pv["CHROM"] + ":" + df_pv["VOI"] + ":N=" + df_pv["Total"].astype(str)

df_pv = df_pv.rename(columns={'Mixed': 'Minor', 'Mutant':'Major'})
SNPvals=df_pv[["Snps",'Minor','Major','Total','Annotation']]

# Separate synonymous and missense:
SNPs_NS  = SNPvals[SNPvals['Annotation']  == "missense_variant"]


######################## Novel missense SNPS graph

#Setup for loading
Totes = SNPs_NS.groupby('Snps')['Total'].sum().reset_index()
Minor = SNPs_NS.groupby('Snps')['Minor'].sum().reset_index()
Major = SNPs_NS.groupby('Snps')['Major'].sum().reset_index()

#Math and definition of SNPratio
Minor['SNPratio'] = [i / j for i,j in zip(Minor['Minor'], Totes['Total'])]
Major['SNPratio'] = [i / j for i,j in zip(Major['Major'], Totes['Total'])]

AllTogether = pd.concat([Minor.Snps, Minor.SNPratio, Major.SNPratio], axis=1)

AllTogether.columns = ['Snps','Minor: AF < 95%', 'Major: AF >= 95%']

df_table_SNP=AllTogether.sort_values(by=['Snps'])

df_table_SNP["index"]=df_table_SNP.Snps.str.split(":").str[1].str[1:-1]
df_table_SNP["index"] = df_table_SNP["index"].str.extract('(\d+)').astype(int)
df_table_SNP["index2"]=df_table_SNP.Snps.str.split(":").str[0]

plot = df_table_SNP.sort_values(by = ['index2', 'index'],ascending=False)[['Snps',  'Minor: AF < 95%','Major: AF >= 95%']].plot(x='Snps', kind='barh', stacked=True, title='Novel missense Mutations', figsize=(20,20), color={"Minor: AF < 95%": "#F3ABA8", "Major: AF >= 95%": "#98DAA7"})

plot.legend(ncol = 2, loc = 'lower right')
plot.set(ylabel="SNPs")
plot.set(xlabel="SNP ratio")
plot.legend(loc=(1,0))
plt.savefig('SNPs-Novel-missense.pdf')




######################## Novel Synonymous SNPS graph

SNPs_S  = SNPvals[SNPvals['Annotation']  == "synonymous_variant"]

#Setup for loading
Totes = SNPs_S.groupby('Snps')['Total'].sum().reset_index()
Minor = SNPs_S.groupby('Snps')['Minor'].sum().reset_index()
Major = SNPs_S.groupby('Snps')['Major'].sum().reset_index()

#Math and definition of SNPratio
Minor['SNPratio'] = [i / j for i,j in zip(Minor['Minor'], Totes['Total'])]
Major['SNPratio'] = [i / j for i,j in zip(Major['Major'], Totes['Total'])]

AllTogether = pd.concat([Minor.Snps, Minor.SNPratio, Major.SNPratio], axis=1)

AllTogether.columns = ['Snps','Minor: AF < 95%', 'Major: AF >= 95%']

df_table_SNP=AllTogether.sort_values(by=['Snps'])

df_table_SNP["index"]=df_table_SNP.Snps.str.split(":").str[1].str[1:-1]
df_table_SNP["index"] = df_table_SNP["index"].str.extract('(\d+)').astype(int)
df_table_SNP["index2"]=df_table_SNP.Snps.str.split(":").str[0]

plot = df_table_SNP.sort_values(by = ['index2', 'index'],ascending=False)[['Snps', 'Minor: AF < 95%','Major: AF >= 95%']].plot(x='Snps', kind='barh', stacked=True, title='Novel synonymous Mutations', figsize=(20,20), color={"Minor: AF < 95%": "#F3ABA8", "Major: AF >= 95%": "#98DAA7"})
#plot.legend(bbox_to_anchor=(0.97, 0.1))

plot.legend(ncol = 2, loc = 'lower right')
#sns.despine(left = True, bottom = True)
plot.set(ylabel="SNPs")
plot.set(xlabel="SNP ratio")
plot.legend(loc=(1,0))
plt.savefig('SNPs-Novel-synonymous.pdf')
