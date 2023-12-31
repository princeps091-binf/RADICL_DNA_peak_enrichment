#%%
import pandas as pd
import numpy as np
import bioframe as bf
from os import listdir
import re
import hvplot.pandas
from scipy.stats import chi2_contingency

#%%
sample = "IPSC_replicate1"
peak_read_list_file = f"/home/vipink/Documents/FANTOM6/RADICL_DNA_peak_enrichment/data/processed/{sample}_peak_read_list.pkl"
clean_RADICL_by_RNA_folder = "/home/vipink/Documents/FANTOM6/RADICL_DNA_peak_enrichment/workflow/data/processed/chr/"
transcript_annotation_file = "/home/vipink/Documents/FANTOM6/data/annotation/FANTOM_CAT/FANTOM_CAT.lv3_robust.info_table.gene.tsv"
# %%
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
    'geneType' : transcript_annotation_df.geneType,
    'geneCategory' : transcript_annotation_df.geneCategory
})

#%%
def get_gene_class_read_count(gene_df,gene_Class,read_df):
    tmp_gene_df = bf.merge(gene_df.groupby('geneClass').get_group(gene_Class),on=['strand'])
    read_in_count = bf.overlap(read_df,tmp_gene_df,on=['strand'],how='inner',return_index=True,return_input=False).shape[0]
    read_out = read_df.shape[0] - read_in_count
    return(pd.DataFrame({
        'io' : ['in','out'],
         'count' : [read_in_count,read_out],
         'anno' : gene_Class
    }))
    #%%

# %%

chr_set = np.unique(np.array([re.split("_|\.",f)[6] for f in listdir(clean_RADICL_by_RNA_folder)])) 
gene_Classes = gene_df.geneClass.drop_duplicates().to_list()
# %%
dfs = []
for chromo in chr_set:
    print(chromo)
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
    chr_gene_classes = gene_df.query("chrom == @chromo").geneClass.drop_duplicates().to_list()
    peak_count_df = (pd.concat([get_gene_class_read_count(gene_df.query("chrom == @chromo"),gene_class,peak_read_df) for gene_class in chr_gene_classes])
     .assign(read_set="peak",chrom=chromo))
    chr_count_df = (pd.concat([get_gene_class_read_count(gene_df.query("chrom == @chromo"),gene_class,chr_read_df) for gene_class in chr_gene_classes])
     .assign(read_set="all",chrom=chromo))
    chr_summary_df = pd.concat([peak_count_df,chr_count_df,]).reset_index()
    dfs.append(chr_summary_df)

# %%
genome_count_tbl = (pd.concat(dfs)
 .groupby(['io','anno','read_set'])
 .agg(nread=('count','sum'))
 .reset_index())
# %%
genome_geneClasses = genome_count_tbl.anno.drop_duplicates().to_list()
# %%
tmpClass = genome_geneClasses[0]
# %%
def compute_OR(tmpClass,genome_count_tbl):
    peak_in = int(genome_count_tbl.query("anno == @tmpClass").query("read_set == 'peak' & io == 'in'").nread.to_numpy())
    peak_out = int(genome_count_tbl.query("anno == @tmpClass").query("read_set == 'peak' & io == 'out'").nread.to_numpy())
    all_in = int(genome_count_tbl.query("anno == @tmpClass").query("read_set == 'all' & io == 'in'").nread.to_numpy())
    all_out = int(genome_count_tbl.query("anno == @tmpClass").query("read_set == 'all' & io == 'out'").nread.to_numpy())
    return((peak_in/(peak_in+peak_out))/(all_in/(all_in+all_out)))
#%%
pd.DataFrame({
    "GeneClass":genome_geneClasses,
    "OR":[compute_OR(tmpClass,genome_count_tbl) for tmpClass in genome_geneClasses]})
#%%
obs = np.array([[peak_in, all_in], [peak_out, all_out]])
res = chi2_contingency(obs)
print((peak_in/(peak_in+peak_out))/(all_in/(all_in+all_out)))
print(res[1])
print(peak_in/(all_in+all_out))

# %%
