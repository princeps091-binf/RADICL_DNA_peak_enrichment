#%%
import pandas as pd
import bioframe as bf
import hvplot.pandas
import numpy as np
#%%
gff_file = "/home/vipink/Documents/FANTOM6/data/annotation/hg38_gencode39_cat_ontcage_talon.gff3"
#%%
gff_df = pd.read_csv(gff_file,sep="\t",comment="#",header=None)

# %%
gene_tbl = (gff_df.iloc[:,[0,2,3,4,6,8]]
 .assign(meta_data=lambda df_:df_.iloc[:,5].str.extract(r'(gene_id=\w+\.?\w+)'))
 .assign(gene_meta=lambda df_:df_.iloc[:,5].str.extract(r'(gene_type=\w+)'))
 .assign(gene_meta_name=lambda df_:df_.iloc[:,5].str.extract(r'(gene_name=\w+)'))
 .assign(gene_id=lambda df_:df_.meta_data.str.split('=').str.get(1),
         gene_type=lambda df_:df_.gene_meta.str.split('=').str.get(1),
         gene_name=lambda df_:df_.gene_meta_name.str.split('=').str.get(1))
 .rename(columns={
     0:'chrom',
     2:'kind',
     3:'start',
     4:'end',
     6:'strand'
 })
 .groupby('gene_id')
 .agg(chrom=('chrom','first'),
      start=('start','min'),
      end=('end','max'),
      strand=('strand','first'),
      gene_type=('gene_type','first'),
      gene_name=('gene_name','first'))
 .reset_index()
)

# %%
(gene_tbl
 .assign(width=lambda df_:np.log10(df_.end-df_.start))
 .hvplot
 .kde('width',by='gene_type',subplots=True,shared_axes=False).cols(2))
# %%
gene_tbl[gene_tbl.gene_type.str.contains("lnc")]
# %%
