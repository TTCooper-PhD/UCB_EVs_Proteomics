#%%
import csv
import json
import math
import matplotlib as plt 
from matplotlib_venn import venn2,venn3
import os
import pandas as pd
import seaborn as sns
import statsmodels,statistics



protein=pd.read_csv('FragPipeAnalysis/combined_protein.tsv',sep='\t')
peptide=pd.read_csv('FragPipeAnalyses/combined_peptide.tsv',sep='\t')
ion=pd.read_csv('combined_ion.tsv',sep='\t')

#collect sample names--> list of sample names

samples=[x.split()[0] for x in protein.columns if "MaxLFQ" in x]
print(samples)

#Create Dictionary of all Data (DaD)
pro_root=['Protein', 'Protein ID', 'Entry Name', 'Gene', 'Protein Length',
    'Organism', 'Protein Existence', 'Description', 'Protein Probability',
    'Top Peptide Probability', 'Combined Total Peptides','Combined Spectral Count',
    'Combined Unique Spectral Count','Combined Total Spectral Count']
DAD={}

for sample in samples:
    temp=[x for x in protein.columns if sample in x]
    filterunit=temp+pro_root
    collection=protein[filterunit]
    DAD[sample]=collection

#create dictionaries that contain number of proteins identified and the IDs (pro_counts, pro_ids)
### Create Dictionary to Store Protiens ID'd in Each Fraction (Pro_Frac)
pro_frac={}
for sample in samples:
    frac=sample.split("_")[1]
    pro_frac[frac]=[]
pro_counts={}
pro_ids={}
id_grab="Protein ID"
for sample in samples:
    temp=DAD[sample]
    frac=sample.split("_")[1]
    f=f"{sample} Spectral Count"
    root_df=temp.loc[temp[f] > 0]
    pro_id=root_df[id_grab]
    count=0
    lst_id=[]
    for pro in pro_id:
        count+=1
        lst_id.append(pro)
    pro_counts[sample]=count
    pro_id[sample]=lst_id
    tmp_list=pro_frac[frac]+lst_id
    tmp_set=sorted(set(tmp_list))
    pro_frac[frac]=tmp_set

for i,j in pro_frac.items():
    print(i)
    print(len(j))

print(pro_counts)

fractions=['A','B','C','D','E','F',"Neat","Raw"]
pro_grouped_counts={}
for fraction in fractions:
    xi=[]
    print("____")
    for i,j in pro_counts.items():
        x=i.split("_")[1]
        if x == fraction:
            print("ya")
            xi.append(j)
    pro_grouped_counts[fraction]=xi
pro_gcount_df=pd.DataFrame(pro_grouped_counts)
print(pro_gcount_df)

pg_T=pro_gcount_df.T
pg_T.rename(
    columns={0: "UCB58", 1: "UCB78", 2:"UCB79"},
    inplace=True,
)

# sns.boxplot(data=pro_gcount_df)

print(pro_frac)
#Checjedf

samp1="C"
samp2="Raw"
a=set(pro_frac[samp1])
b=set(pro_frac[samp2])
a_lean=len(a.difference(b))
b_lean=len(b.difference(a))
mid=len(a.intersection(b))



# venn2(subsets = (a_lean,b_lean,mid), set_labels = (samp1, samp2))

v100=pd.read_excel('Vesicle_100.xlsx')
vpedia_100=set(v100["Vesicle_100"])
vp_overlap=a.intersection(vpedia_100)
# %%
