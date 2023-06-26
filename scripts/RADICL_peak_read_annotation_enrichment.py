#%%
import pandas as pd
import bioframe as bf
import scipy.stats as stats
import numpy as np
from scipy.stats import chi2_contingency
import hvplot.pandas

#%%
read_folder = "/home/vipink/Documents/FANTOM6/RADICL_hotspot/workflow/data/processed/DNA/chr/"
transcript_annotation_file = "/home/vipink/Documents/FANTOM6/data/annotation/FANTOM_CAT/FANTOM_CAT.lv3_robust.info_table.gene.tsv"
cre_file = "/home/vipink/Documents/FANTOM6/data/annotation/GRCh38-cCREs.bed"
enh_file = "/home/vipink/Documents/FANTOM6/data/annotation/GRCh38-ELS.bed"
peak_folder = "/home/vipink/Documents/FANTOM6/alternative_filter_blacklist_pipeline/workflow/data/results/DNA/"
sample = "IPSC_replicate1"
#%%
# load the annotation tables
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
    'geneType' : transcript_annotation_df.geneType,
    'geneCategory' : transcript_annotation_df.geneCategory
})
#%%
gene_df.groupby('geneClass').agg(nr=('start','count'))
#%%
enhancer_df = (pd.read_csv(enh_file,sep="\t",header=None,dtype = {
    0:str,
    1:int,
    2:int,
    3:str,
    4:str,
    5:str
})
.rename(columns={
    0:'chrom',
    1:'start',
    2:'end',
    3:'ID1',
    4:'ID2',
    5:'label'
}))
#%%
cre_df = (pd.read_csv(cre_file,sep="\t",header=None,dtype = {
    0:str,
    1:int,
    2:int,
    3:str,
    4:str,
    5:str
})
.rename(columns={
    0:'chrom',
    1:'start',
    2:'end',
    3:'ID1',
    4:'ID2',
    5:'label'
}))

# %%
# load the peak table
# %%
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
chr_set = peak_df.chrom.drop_duplicates().to_list()
#%%
def produce_count_df(read_df,peak_df,enhancer_df,chromo):
    peak_read = bf.count_overlaps(read_df,peak_df.query("chrom == @chromo")).query('count > 0')
    chr_peak_read_inter = bf.count_overlaps(peak_read,enhancer_df.query('chrom == @chromo')).query('count > 0').shape[0]
    chr_read_inter = bf.count_overlaps(read_df,enhancer_df.query('chrom == @chromo')).query('count > 0').shape[0]

    chr_count_tbl = pd.DataFrame({
                "chrom" : chromo,
                "set" : ["peak",'all','peak','all'],
                "io" : ["in","in","out","out"],
                "count" : np.array([chr_peak_read_inter,chr_read_inter,peak_read.shape[0]-chr_peak_read_inter,read_df.shape[0]-chr_read_inter])
            })
    return(chr_count_tbl)

def produce_anno_read_count(read_folder,sample,chr_set,peak_df,anno_df):
    dfs = []

    for chromo in chr_set:
        print(chromo)
        read_df = (pd.read_csv(f"{read_folder}{sample}_clean_DNA_{chromo}.bed",
                            sep="\t",header=None,
                            usecols=[0,1,2],
                            dtype={
                                0:str,
                                1:int,
                                2:int
                                    })
                    .rename(columns={
                        0:'chrom',
                        1:'start',
                        2:'end'
                    }))
        chr_count_tbl = produce_count_df(read_df,peak_df,anno_df,chromo)
        dfs.append(chr_count_tbl)
    count_tbl = pd.concat(dfs)
    return (count_tbl)

#%%
count_tbl = produce_anno_read_count(read_folder,sample,chr_set,peak_df,cre_df)
#%%
genome_summary = count_tbl.groupby(['set','io']).agg(nread=('count','sum')).reset_index()
plot = (genome_summary
 .merge(genome_summary.groupby('set').agg(read_sum= ('nread','sum')).reset_index())
 .assign(nratio= lambda df_:df_.nread/df_.read_sum)
 .loc[:,['set','io','nratio']]
 .pivot(index='set', columns='io',values='nratio')
 .hvplot.bar(y=["in","out"],cmap='Set1',stacked=True))
hvplot.save(plot,"test.html")
plot
# %%
peak_in = int(genome_summary.query("set == 'peak' & io == 'in'").nread.to_numpy())
peak_out = int(genome_summary.query("set == 'peak' & io == 'out'").nread.to_numpy())
all_in = int(genome_summary.query("set == 'all' & io == 'in'").nread.to_numpy())
all_out = int(genome_summary.query("set == 'all' & io == 'out'").nread.to_numpy())

obs = np.array([[peak_in, all_in], [peak_out, all_out]])
res = chi2_contingency(obs)
print((peak_in/(peak_in+peak_out))/(all_in/(all_in+all_out)))
print(res[1])
print(peak_in/(all_in+all_out))
# %%
