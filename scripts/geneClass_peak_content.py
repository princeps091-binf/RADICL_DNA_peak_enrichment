#%%
import pandas as pd
import numpy as np
import bioframe as bf
from os import listdir
import re
import hvplot.pandas
#%%
sample = "IPSC_replicate1"
peak_read_list_file = f"/home/vipink/Documents/FANTOM6/RADICL_DNA_peak_enrichment/data/processed/{sample}_peak_read_list.pkl"
clean_RADICL_by_RNA_folder = "/home/vipink/Documents/FANTOM6/RADICL_DNA_peak_enrichment/workflow/data/processed/chr/"
transcript_annotation_file = "/home/vipink/Documents/FANTOM6/data/annotation/FANTOM_CAT/FANTOM_CAT.lv3_robust.info_table.gene.tsv"
#%%
peak_read_list = pd.read_pickle(peak_read_list_file)

transcript_annotation_df = pd.read_csv(transcript_annotation_file,delimiter="\t")
loc_column = transcript_annotation_df.loc[:,"loc"].str.split(":|,")
gene_chromo = loc_column.apply(lambda x:x[0])
gene_start = loc_column.apply(lambda x:int(x[1].split('-')[0]))
gene_end = loc_column.apply(lambda x:int(x[1].split('-')[1]))
gene_strand = loc_column.apply(lambda x:x[2])
gene_df = pd.DataFrame({
    'chrom' : gene_chromo,
    'start' : gene_start,
    'end' : gene_end,
    'strand' : gene_strand,
    'ID' : transcript_annotation_df.geneID,
    'geneClass' : transcript_annotation_df.geneClass,
})
chr_set = np.unique(np.array([re.split("_|\.",f)[6] for f in listdir(clean_RADICL_by_RNA_folder)])) 
#%%
dfs = []
for chromo in chr_set:
    print(chromo)
    #%%
    chr_read_df = (pd.read_csv(f"{clean_RADICL_by_RNA_folder}{sample}_clean_RADICL_by_RNA_{chromo}.txt",
                             sep="\t",
                             header=None,
                             usecols=[0,1,2,3,5])
                             .rename(columns={
                                        0:'chrom',
                                        1:'start',
                                        2:'end',
                                        3:'RNA_ID',
                                        5:'strand'
                             }))
    peak_read_df = chr_read_df.query('RNA_ID in @peak_read_list')
    #%%
    all_count = bf.count_overlaps(gene_df.query("chrom == @chromo"),chr_read_df,on=['strand'])
    peak_count = bf.count_overlaps(gene_df.query("chrom == @chromo"),peak_read_df,on=['strand']).rename(columns={
        'count':'peak_count'
    })
    #%%
    chr_count_df = all_count.merge(peak_count,how='left')
    #%%
    dfs.append(chr_count_df)
#%%
genome_count = pd.concat(dfs)

(genome_count
 .assign(log_count = lambda df_:np.log10(df_.loc[:,'peak_count']))
 .assign(log_count = lambda df_:df_.log_count.mask(np.isinf,-5))
 .hvplot.kde('log_count',by='geneClass'))
# %%
(genome_count
 .query("count > 0")
 .assign(inpeak = lambda df_:df_.peak_count > 0)
 .assign(log_count = lambda df_:np.log10(df_.loc[:,'count']))
 .assign(log_count = lambda df_:df_.log_count.mask(np.isinf,-5))
 .hvplot.kde('log_count',by='inpeak'))

# %%
(genome_count
 .assign(log_peak_count = lambda df_:np.log10(df_.loc[:,'peak_count']))
 .assign(log_peak_count = lambda df_:df_.log_peak_count.mask(np.isinf,-1))
 .assign(log_all_count = lambda df_:np.log10(df_.loc[:,'count']))
 .assign(log_all_count = lambda df_:df_.log_all_count.mask(np.isinf,-1))

 .hvplot.scatter(x='log_all_count',y='log_peak_count',alpha=0.1,size=0.1))

# %%
