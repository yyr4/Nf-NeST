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
parser.add_argument('-G', '--Gff file', dest='GFF_file', type=str,
                        help='GFF_file')
parser.add_argument('-N', '--cov', dest='sam_depth', type=str,
                        help='coverage file')
parser.add_argument('-V', '--VOI_file', dest='VOI', type=str,
                        help='Variants of interest file')

args = parser.parse_args()
Cov_file = args.sam_depth
ref = args.Reference_fasta
Gff_file = args.GFF_file
VOI_file = args.VOI

Sample_name = Cov_file.split("/")[-1].split("_")[0]
Sample_out = ("_".join(Cov_file.split("/")[-1].split("_")[0:-1]))


#### parse the bed file
#### parse the bed file

Bed = pd.read_csv(Gff_file, sep='\t', names=('chr','Source', 'type','start', 'stop','score','strand','Phase','attributes'), header=None)
Bed = Bed.dropna()

df2 = Bed.loc[Bed['type'].str.contains('CDS', na=False, case=False)]## Here we are first align the each base of fasta file to per base coverage position and create dataframe
df2[['start', 'stop']] = df2[['start', 'stop']].astype(int)
#df2['length'] = df2['stop'] - df2['start']



new = []
with Fasta(ref) as fasta:     ## fasta file with each genes
    with open(Cov_file, 'r') as nucleotides:   ## samtools depth output
            nucleotides.seek(0) # Ensure you're at the start of the file..
            first_char = nucleotides.read(1) # Get the first character
            if not first_char:
                raise Exception(Sample_out, "is empty!")

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



all =[]
for i,j in df.iterrows():
    for a,b in df2.iterrows():
        only_CDS_CRT= {}
        if (j['chr'] == b['chr']) & ((b['start'] <= j["pos"]) & (j['pos'] < b['stop']+1)):
            key = j['chr']
            value = j['pos']
            only_CDS_CRT[key]=value
            all.append(only_CDS_CRT)

all_keys =[i for s in [d.keys() for d in all] for i in s]
all_val = [i for s in [d.values() for d in all] for i in s]
df_new = pd.DataFrame()
df_new["chr"] = all_keys
df_new["pos"] = all_val


only_CDS = pd.merge(df_new, df, on=['chr','pos'])

only_CDS = only_CDS.reset_index(drop=True)
only_CDS[['pos','cov']] = only_CDS[['pos','cov']].astype(str)
df_ind3 = only_CDS.groupby(only_CDS.index // 3).agg(','.join)



df_ind = df_ind3.reset_index(drop=True)


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

df_ind = df_ind.rename({'Chrom': 'CHROM'}, axis=1)
df_ind.to_csv("1.csv")

df_voi=pd.read_csv(VOI_file)
df_voi["VOI"]=df_voi["RefAA"]+df_voi["AAPos"].astype(str)+df_voi["AltAA"]
df_voi["AA_change"]=df_voi["RefAA"]+df_voi["AAPos"].astype(str)
df_voi = df_voi.rename({'Chr':'CHROM'},axis=1)

merge_df = pd.merge(df_ind, df_voi, on=['CHROM','AA_change'])


merge_df = merge_df.drop(['Gene', 'RefAA',  'AAPos', 'AltAA', 'Codon_position','AA','snp_Report' ], axis=1)


name = merge_df.to_csv(Sample_out+'_coverage.csv')
