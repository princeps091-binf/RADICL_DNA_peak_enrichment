#%%
import pandas as pd
import numpy as np
import bioframe as bf
import pickle

#%%
RADICL_folder = "/home/vipink/Documents/FANTOM6/alternative_filter_blacklist_pipeline/workflow/data/processed/chr/"
read_folder = "/home/vipink/Documents/FANTOM6/alternative_filter_blacklist_pipeline/workflow/data/processed/DNA/chr/"
peak_folder = "/home/vipink/Documents/FANTOM6/alternative_filter_blacklist_pipeline/workflow/data/results/DNA/"
sample = "Neuron_replicate1"
out_file = f"/home/vipink/Documents/FANTOM6/RADICL_DNA_peak_enrichment/data/processed/{sample}_peak_read_list.pkl"
#%%
peak_df = (pd.read_csv(f"{peak_folder}{sample}/{sample}_all_peak.bed",sep="\t",header=None,usecols=[0,1,2],dtype={
    0:str,
    1:int,
    2:int
})
.rename(columns={
    0:'chrom',
    1:'start',
    2:'end'
}))
#%%

def extract_peak_read_ID(read_folder,sample,chr_set,peak_df):
    dfs = []
    for chromo in chr_set:
        print(chromo)
        read_df = (pd.read_csv(f"{read_folder}{sample}_clean_DNA_{chromo}.bed",
                            sep="\t",header=None,
                            usecols=[0,1,2,3],
                            dtype={
                                0:str,
                                1:int,
                                2:int,
                                3:str
                                    })
                    .rename(columns={
                        0:'chrom',
                        1:'start',
                        2:'end',
                        3:"ID"
                    }))
        chr_peak_read_ID = bf.overlap(read_df,peak_df.query("chrom == @chromo"),how='inner',return_input=True).loc[:,'ID'].to_list()
        dfs+=chr_peak_read_ID    
    return (dfs)
#%%
chr_set = peak_df.chrom.drop_duplicates().to_list()
#%%
peak_read_ID = extract_peak_read_ID(read_folder,sample,chr_set,peak_df)
#%%
with open(out_file, 'wb') as fp:
    pickle.dump(peak_read_ID,fp)