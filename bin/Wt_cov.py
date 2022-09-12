#!/usr/bin/env python3

import pandas as pd
import os
import glob
import numpy as np
import argparse
import sys
from pyfaidx import Fasta



parser = argparse.ArgumentParser(description='find_COV')

parser.add_argument('-R', '--referance', dest='Reference_fasta', type=str,
                        help='Ref file')
parser.add_argument('-B', '--Bed_file', dest='Bed_file', type=str,
                        help='Bed_file')
parser.add_argument('-N', '--cov', dest='sam_depth', type=str,
                        help='coverage file')
parser.add_argument('-V', '--VOI_file', dest='VOI', type=str,
                        help='Variants of interest file')

args = parser.parse_args()
Cov_file = args.sam_depth
ref = args.Reference_fasta
Bed_file = args.Bed_file
VOI_file = args.VOI

Sample_name = Cov_file.split("/")[-1].split("_")[0]
Sample_out = ("_".join(Cov_file.split("/")[-1].split("_")[0:3]))




#### parse the bed file
Bed = pd.read_csv(Bed_file, sep='\t', names=('chr', 'start', 'stop','Name','Phase', 'strand'), header=None)
df2 = Bed.loc[Bed['Name'].str.contains('CDS')]
#print(CDS)

## Here we are first align the each base of fasta file to per base coverage position and create dataframe

new = []
with Fasta(ref) as fasta:     ## fasta file with each genes
    with open(Cov_file, 'r') as nucleotides:   ## samtools depth output
            nucleotides.seek(0) # Ensure you're at the start of the file..
            first_char = nucleotides.read(1) # Get the first character
            if not first_char:
                sys.exit()
            else:
                nucleotides.seek(0) # The first character wasn't empty. Return to the start of the file.
                for line in nucleotides:
                    chrom, pos, cov = line.rstrip().split()
                    nuc = fasta[chrom][int(pos)-1].seq     ## pyfaidx module used to get fasta sequances from fasta file
                    i = ("{chrom}\t{pos}\t{nuc}\t{cov}".format(**locals()))
                    new.append(i)

Data = [i.split("\t") for i in new]

df = pd.DataFrame(Data, columns=['chr', 'pos','seq','cov'])
df[['pos','cov']] = df[['pos','cov']].astype(int)
#print(df)

df1 = df.loc[df['chr'].isin(['K13','DHFR','PfMDR1','mitochondrial_genome_-_CYTB_CDS'])]  ## find the rows with genes witout introns
df1 = df1.reset_index(drop=True)
#print(df1)
df1[['pos','cov']] = df1[['pos','cov']].astype(str)
df_ind1 = df1.groupby(df1.index // 3).agg(','.join)
###### DHPS
df3 = df.loc[df['chr'].isin(['DHPS_437Corrected'])]

df4 = df2.loc[df2['chr'].isin(['DHPS_437Corrected'])]

only_CDS = []
for i,j in df3.iterrows():
    for a,b in df4.iterrows():
        if (b['start'] < j["pos"]) & (j['pos'] < b['stop']+1):
            only_CDS.append(j['pos'])

dhps = df3.loc[df3['pos'].isin(only_CDS)]
dhps = dhps.reset_index(drop=True)
dhps[['pos','cov']] = dhps[['pos','cov']].astype(str)
df_ind2 = dhps.groupby(dhps.index // 3).agg(','.join)
#########PFCRT
df5 = df.loc[df['chr'].isin(['PfCRT'])]
df6 = df2.loc[df2['chr'].isin(['PfCRT'])]

only_CDS_CRT = []
for i,j in df5.iterrows():
    for a,b in df6.iterrows():
        if (b['start'] < j["pos"]) & (j['pos'] < b['stop']+1):
            only_CDS_CRT.append(j['pos'])

CRT = df5.loc[df5['pos'].isin(only_CDS_CRT)]
CRT = CRT.reset_index(drop=True)
CRT[['pos','cov']] = CRT[['pos','cov']].astype(str)
df_ind3 = CRT.groupby(CRT.index // 3).agg(','.join)
########### COncatenate all the dataframe

result = pd.concat([df_ind1,df_ind2, df_ind3], axis = 0)
df_ind = result.reset_index(drop=True)

df_ind['Chrom'] = df_ind['chr'].str.split(',', expand=True)[0]  # get gene name
df_ind["Codon_position"] = 1 + df_ind.groupby("Chrom").cumcount()  # get each codon position (3 nuc = 1 codon)
df_ind['AVG_COV'] = df_ind['cov'].apply(lambda x: sum(map(int, x.split(','))))   # sum of  cov of three consecutive bases
df_ind['AVG_COV'] = df_ind['AVG_COV']/3       # divide sum by 3
df_ind['triplet'] = df_ind['seq'].str.replace(',','')   # get three codon bases to convert into amino acid
df_ind.drop(['chr', 'seq','pos', 'cov'], axis=1, inplace=True)  # drop the two colums



AA_table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }

df_ind['AA'] = df_ind['triplet'].map(AA_table)
df_ind['Sample_name'] = Sample_name

df_ind['AA_change'] = df_ind['AA'] + df_ind['Codon_position'].astype(str)
column_names = ['Sample_name',"Chrom", "Codon_position", "AVG_COV","triplet", 'AA', 'AA_change']
df_ind = df_ind.reindex(columns=column_names)
df_ind = df_ind.replace({'DHPS_437Corrected': 'DHPS'})
df_ind = df_ind.rename({'Chrom': 'CHROM'}, axis=1)


df_voi=pd.read_csv(VOI_file)
df_voi["VOI"]=df_voi["RefAA"]+df_voi["AAPos"].astype(str)+df_voi["AltAA"]
df_voi["AA_change"]=df_voi["RefAA"]+df_voi["AAPos"].astype(str)
df_voi = df_voi.rename({'Gene':'CHROM'},axis=1)

merge_df = pd.merge(df_ind, df_voi, on=['CHROM','AA_change'])
merge_df.drop(['Chr', 'RefAA',  'AAPos', 'AltAA', 'Codon_position','AA' ], axis=1, inplace=True)


name = merge_df.to_csv(Sample_out+'_coverage.csv')
