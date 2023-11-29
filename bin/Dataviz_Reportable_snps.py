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

#######################  DMS EPI Report:

df_epi = DF[['Sample_name','CHROM','VOI','Type']]
def func(row):
    if row['Type'] == 'Mutant':
        return '1'
    elif row['Type'] =='Mixed':
        return '1'
    else:
        return '0'

df_epi['Total Mutation'] = df_epi.apply(func, axis=1).astype(int)

df_epi['codon'] = df_epi['VOI'].str[1:-1].astype(int)

df_epi['Gene_Variants'] =  df_epi['CHROM']  +":" +df_epi['VOI']

df2_epi =df_epi.groupby(['Sample_name','CHROM'])['Total Mutation'].sum().reset_index()
#print(df_epi)




pivot_table_1 = df_epi.pivot_table(index=['Gene_Variants','CHROM','codon'],columns=['Sample_name'], values='Type',sort=False ,aggfunc=lambda x: ' '.join(str(v) for v in x))
pivot_table_1 = pivot_table_1.sort_values(['CHROM','codon'])

pivot_table_1 = pivot_table_1.T



pivot_table_1 = pivot_table_1.droplevel([1,2], axis=1)
#print(pivot_table_1)



pivot_table_2 = df2_epi.pivot_table(index=['Sample_name'], columns=['CHROM'], values='Total Mutation').reset_index()
pivot_table_2 = pivot_table_2.add_suffix(': #Drug_resistant_mutations')
pivot_table_2 = pivot_table_2.rename(columns={"Sample_name: #Drug_resistant_mutations": "Sample_name"})

#print(pivot_table_2)


Merge_pivot = (pd.merge(pivot_table_2, pivot_table_1, how="outer", on=["Sample_name"])
             .set_index("Sample_name")
             .reset_index()

)

Merge_pivot.insert(loc=1, column='CaseID', value='NA' )
Merge_pivot.insert(loc=2, column='Travel_History', value='NA')
Merge_pivot = Merge_pivot.replace(['Wildtype', 'Mutant','Mixed','No coverage' ], ['WT', 'MT', 'MIX', 'NA'])

Merge_pivot.to_csv("DMS_EPI_report.csv", sep=',')


##########  Reportable_Per_SNP_depth
df1 = DF[['CHROM', 'VOI', 'AVG_COV']]

df1["variants"] = df1[['CHROM',"VOI"]].apply(":".join, axis=1)
df1 = df1.sort_values('variants')
print(df1)
#sns.set_style("whitegrid")
df1["index"]=df1.variants.str.split(":").str[1].str[1:-1]
df1["index"]=df1["index"].astype(int)
df1["index2"]=df1.variants.str.split(":").str[0]

df1 = df1.sort_values(by = ['index2', 'index'],ascending=True)

fig, ax = plt.subplots()
fig.set_size_inches(25, 25)
sns.boxplot(y = 'variants', x = 'AVG_COV', data = df1, hue="CHROM", dodge=False)
sns.despine()

fig.savefig('Reportable_Per_SNP_depth.pdf')

####################### Bar plot for Reportable_snps

df = DF.groupby(['CHROM','VOI','Type']).size().reset_index(name='counts')

df_pv= df.pivot_table(values='counts', index=['CHROM','VOI'], columns='Type', aggfunc='first')
df_pv = df_pv.fillna(0).reset_index()

column_names = ['Mixed','Mutant','Wildtype']
df_pv['Total']= df_pv[column_names].sum(axis=1)
df_pv["Snps"] = df_pv["CHROM"] + ":" + df_pv["VOI"] + ":N=" + df_pv["Total"].astype(str)

df_pv = df_pv.rename(columns={'Mixed': 'Minor', 'Mutant':'Major'})
SNPvals=df_pv[["Snps",'Minor','Major','Wildtype','Total']]

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
plot.legend(loc=(1,0))
plt.savefig('SNPs-Reportable.pdf')
