#!/usr/bin/env python
# coding: utf-8

# # Mid-gestation fetal cortex dataset: Cluster characterization
# 
# _**Single cell transcriptomics dataset from paper published by Trevino et al. (Cell 2021) characterizing human fetal cortex at mid-gestation**._
# 
# __Upstream Steps__
# 
# * Assemble adata
# * QC filter on cells
# * Expression filter on genes
# * Normalization and log10 transformation by Scanpy functions
# * Feature selection (HVG) by Scanpy functions
# * Dimensionality reduction
# * Batch correction by Harmony
# * Cluster identification
# 
# __This notebook__
# 
# * Cluster characterization
# 

# ----

# # 1. Environment Set Up
# 
# ## 1.1 Library upload

# In[1]:


import os
import sys

import numpy as np
import pandas as pd
import igraph as ig
import scanpy as sc
import scanpy.external as sce
from scipy.sparse import csr_matrix, isspmatrix
import gseapy as gp

#Plotting
import matplotlib.pyplot as plt
import seaborn as sns

#ultils
#import ipynbname
from datetime import datetime


# In[2]:


# Custom functions
sys.path.append('HelperFunctions')
import Day1Helper as fn


# In[3]:


sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=80, fontsize=12)


# ## 1.2 Start Computation time

# In[4]:


print(datetime.now())


# ----

# # 2. Read input files  

# In[5]:


path = '../DataDir/InputData/'
Id = 'Day1_2_TrevinoFiltNormAdata.h5ad'

input_file = path + Id

adata = sc.read(input_file)


# In[6]:


print('Loaded Normalizes AnnData object: number of cells', adata.n_obs)
print('Loaded Normalizes AnnData object: number of genes', adata.n_vars)

print('Available metadata for each cell: ', adata.obs.columns)


# # 3. Cluster characteristics

# ## 3.1 Check with original annotation

# In[7]:


sc.settings.set_figure_params(fontsize=9, figsize=[8, 6])

sc.pl.embedding(adata, basis="X_umap_harmony", color=['Leiden_Sel','cell_label', 'cell_class'],
                legend_loc='on data')


# In[8]:


pd.crosstab(adata.obs['cell_class'], adata.obs['Leiden_Sel'], normalize = 'columns', margins=True)


# > `_Try yourself_:` 
# > `what if you want to see total number of cells instead of %?`

# In[9]:


# write your code here


# > `_Try yourself_:` 
# > `what if you want to compare with the more granular annotation?`

# In[10]:


# write your code here


# #### Comparison outcome
# From the comparison above, we can hypothesize a first large-grain definition of our clusters.

# * __Excitatory Neurons__: C2, C3 (IPC), C4, C9
# * __Inhibitory Neurons__: C6, C8, C11
# * __Progenitors__: C1, C7
# * __Glia__: C0 (progenitors, OPC), C10 (microglia)
# * __Other__: C5

# ## 3.2 Cluster dendrogram

# In[11]:


sc.settings.set_figure_params(fontsize=10, figsize=[8, 6])

sc.tl.dendrogram(adata, groupby='Leiden_Sel', n_pcs=5, cor_method='spearman', linkage_method='average')
sc.pl.dendrogram(adata, groupby='Leiden_Sel', orientation='left')


# > `_Try yourself_:` 
# > `How to compare to a dendrogram created on the annotation from the original paper?`

# In[12]:


# write your code here


# ----

# # 4. Population markers
# 
# * Supervised approach: plot on the UMAP your own set of population markers
# * The definition of the markers strongly depends on biological knowledge and liteature survey, and it is context-dependent
# * As an example, we define below a dictionary of markers usually checked in cortex datasets (not comprehensive!!)
# * We use the __customUmap__ function defined in the helper file: check there to go back to the basic scanpy function, for the fine-tuning of the parameters or to create your own custom function
# 

# In[13]:


marker_dictionary = {

    'Proliferating_Progenitors': ['MKI67', 'CDC20', 'HMGB2', 'CCNB1', 'CCNB2', 'ASPM'], 
    'Radial_Glia': ['SOX2', 'PAX6', 'NES', 'VIM', 'HES1', 'GLI3'],
    'Intermediate_Progenitors': ['EOMES', 'ELAVL4', 'NHLH1', 'KCNQ3', 'INSM1', 'HES6'], 
    'oRG': ['FAM107A', 'HOPX',  'PTPRZ1', 'TNC', 'ITGB5'], 
    
    'Neurons': ['GAP43', 'DCX', 'STMN2', 'MAP2', 'SYT1', 'MEF2C'], 
    
    'Excitatory_Progenitors': ['EMX1', 'NEUROD1', 'NEUROD2', 'NEUROD6', 'NEUROG2', 'NEUROG1'],
    'Excitatory_Neurons': ['SLC17A6', 'SLC17A7', 'GRIN2A', 'GRIN2B', 'SLA'], 
    
    'Inhibitory_Early': ['NKX2-1', 'DLX5', 'DLX6', 'DLX1', 'DLX2', 'DLX6-AS1'],
    'Inhibitory_Neurons': ['GAD1', 'GAD2', 'SLC32A1', 'CALB1', 'CALB2', 'NPY', 'SST', 'PVALB', 'VIP'],
    
    'Astrocytes': ['GFAP', 'SLC1A3',  'S100B', 'AQP4', 'ALDH1L1', 'TNC'], 
    'Microglia': ['PTPRC', 'AIF1', 'CCL3', 'ITGAM', 'CX3CR1', 'CD74'], 
    
    'Endothelial_Pericytes': ['CLDN5', 'PECAM1', 'ABCG2', 'FOXC2', 'PDGFRB', 'IFITM2'], 
    }


# In[14]:


sc.settings.set_figure_params(dpi=80, fontsize=12)

for population in marker_dictionary: 
    print(f"{population.upper()}:")
    fn.customUmap(adata, marker_dictionary[population], size=8)
    print("\n\n\n")


# > `_Try yourself_:` 
# > `define below a small dictionary of genes you are interested in and plot them on the UMAP`

# In[15]:


# write your code here


# ------

# # 5. Cluster Top Markers

# ## 5.1 Identify Cluster top-markers
# 
# For each cluster, top-marker genes are identified by comparing the expression in that cluster versus all the others. Scanpy [tl.rank_genes_groups](https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.rank_genes_groups.html) function is used to this purpose.
# * groupby: cluster labels
# * method: statistical method for differential expression analysis. Here we specify Wilcoxon test (non-parametric). 

# In[16]:


adata.uns['log1p']['base'] = None
sc.tl.rank_genes_groups(adata, groupby='Leiden_Sel', method='wilcoxon', key_added='wilcox', 
                       use_raw=False, pts=True)


# In[17]:


GroupMarkers = pd.DataFrame(adata.uns['wilcox']['names']).head(101)
GroupMarkers.columns = 'Cl_' + GroupMarkers.columns

GroupMarkers.head(11)


# ## 5.2 Visualize marker genes
# 
# ### DotPlot

# In[18]:


sc.pl.rank_genes_groups_dotplot(adata, n_genes=4, key='wilcox', standard_scale="var")  
# standard scale: normalize each gene to range from 0 to 1


# ## 5.3 Filter top marker genes
# 
# * To improve the specificity of the markers, we can filter top-markers genes, e.g. considering ther % of expression in the target population and in the other cells.
# * We can then plot again the results with dotplot or matrixplot

# In[19]:


sc.tl.filter_rank_genes_groups(
    adata,
    min_in_group_fraction=0.3,
    max_out_group_fraction=0.2,
    key="wilcox", key_added="wilcox_filt",
)


# In[20]:


sc.pl.rank_genes_groups_dotplot(adata, n_genes=4, key='wilcox', 
                               standard_scale="var") 


# In[21]:


sc.pl.rank_genes_groups_dotplot(
    adata, key='wilcox', n_genes=4,
    values_to_plot="logfoldchanges",
    cmap='bwr', vmin=-4, vmax=4, min_logfoldchange=2,
    colorbar_title='log fold change',
)


# In[22]:


sc.pl.rank_genes_groups_matrixplot(adata, n_genes=5, key='wilcox_filt', standard_scale="var")


# ## 5.4 Marker overlap
# 
# * To gain further information on our clusters, we can overlap the top-genes of each cluster with external gene sets (e.g. reliable sets of population marker genes; genes identified by other studies, etc.).
# * As an example, here we overlap the marker dictionary defined above

# In[23]:


marker_overlap = sc.tl.marker_gene_overlap(adata, marker_dictionary, key='wilcox', 
                         top_n_markers=100, method='jaccard')

marker_overlap.columns = marker_overlap.columns.astype('category')


# In[24]:


sns.set_theme()
f, ax = plt.subplots(figsize=(10, 7))
sns.heatmap(marker_overlap, linewidths=.5, ax=ax, cmap = 'PuBu')


# ----

# # 6. Functional analysis by gseapy

# ## 6.1 Approaches
# 
# ### Overview
# Two main approaches can be employed
# 
# * __over-representation analysis__: test for enrichment of functionally relevant gene sets in the pool of the gene of interest (e.g. genes selected as cluster markers, DEGs) compared to the gene universe (all expressed/tested genes).  
# * __gsea-like approaches__: rank all tested genes according to statistical metrics (e.g. PValue, Fold-Change, etc.) and test the distribution of the gene sets along the list (e.g. more present in the top-part of the list vs evenly distributed).  
# 
# Each of these methods has advantages and drawbacks. In general, take into consideration that functional enrichment analyses are just a way to try to describe from a functional point of view your list of genes; they rely on vocabularies that can be variably complete and tuned to your biological system.
# Considering this, it is wiser to use the results as complementary to evidence coming from other approaches, rather than blindly and acritically rely on them.
# 
# ### Gseapy
# 
# * [Read the docs](https://gseapy.readthedocs.io/en/latest/index.html)
# * [Bioinformatics Paper](https://academic.oup.com/bioinformatics/article/39/1/btac757/6847088)
# * [GitHub repo](https://github.com/zqfang/GSEApy)
# 
# ### Other tools
# 
# * [GProfiler](https://pypi.org/project/gprofiler-official/)
# * [Decoupler](https://decoupler-py.readthedocs.io/en/latest/)
# * [TopGO](https://bioconductor.org/packages/devel/bioc/vignettes/topGO/inst/doc/topGO.pdf)

# In[25]:


# to check available gene sets
gp.get_library_name()[:5]


# ## 6.2 Pre-ranked GSEA
# 
# * the custom function __customGseapy__ (see helper file) first generates the list of marker genes ranked according to scanpy-calculated score for the cluster of interest
# * GSEA is then applied to the pre-ranked list for the gene sets of interests
# * as an example, here we test two genes sets (GO Biological Process; Cell Markers collection) on a subset of clusters.

# In[26]:


# C1: progenitors - C8: inhibitory neurons - C9: excitatory neurons - C10: microglia 

Cl = ['1', '8' , '9', '10']

for i in Cl: 
    print("\n\n {}".format(i)) 
    display(fn.customGseapy(adata, cluster=i, rank='wilcox', 
                            sets = ['GO_Biological_Process_2023'],
                            fdr_th=0.005, nes_th=1.75, show=5))


# # 7. Examine a gene list of interest

# ## 7.1 Define gene signatures
# 

# In[27]:


semaforins = fn.selectMarkers(adata, ["SEMA3A","SEMA3B","SEMA3C","SEMA3E","SEMA3F","SEMA4A","SEMA4B",
              "SEMA4C","SEMA4D","SEMA4F","SEMA4G","SEMA5A","SEMA5B","SEMA6A","SEMA6B","SEMA6C","SEMA6D","SEMA7A"])

#Kegg Steroid biosynthesis genes
steroid = fn.selectMarkers(adata, ["CEL","CYP27B1","CYP51A1","DHCR24","DHCR7","EBP","FDFT1","HSD17B7","LIPA","LSS","MSMO1",
                                   "NSDHL","SC5D","SOAT1","SOAT2","SQLE","TM7SF2"])

TFs = fn.selectMarkers(adata, ['DLX5', 'DLX1', 'SOX9', 'MEIS2', 
      'NEUROG2', 'NEUROD1', 'LHX5', 'LHX2', 
      'LHX9', 'EN1', 'EMX2', 'EMX1'])


# ## 7.2 Visualize levels across clusters

# #### Plot 'semaforins' across clusters by Stacked Violins 

# In[28]:


sc.settings.set_figure_params(figsize=[10, 8])
sc.pl.StackedViolin(adata, semaforins, groupby='Leiden_Sel').show()


# #### Plot 'TFs' across clusters and original annotation by dotplots anda matrixplot

# In[29]:


sc.pl.dotplot(adata, TFs, groupby='Leiden_Sel', swap_axes=True)


# In[30]:


sc.pl.dotplot(adata, TFs, groupby='cell_label', standard_scale="var", cmap='Blues', swap_axes=True)


# In[31]:


sc.pl.matrixplot(adata, TFs, groupby='Leiden_Sel', swap_axes=True, standard_scale="var")


# > `_Try yourself_:` 
# > `plot genes involved in steroid synthesis stratified by clusters or cell labels`

# In[32]:


# write your code here


# ## 7.3 Calculate and visualize gene scores
# 
# For a gene signature that you expect to be expressed coherently in the same cell population, you can employ the [tl.score.genes](https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.score_genes.html) to estimate an average expression score, and then plot it on the UMAP.  

# In[33]:


sc.tl.score_genes(adata, steroid, score_name='steroid_score')


# In[34]:


sc.pl.embedding(adata, basis="X_umap_harmony", 
                color='steroid_score')


# # 8. Save

# ## 8.1 Timestamp finished computations 

# In[35]:


print(datetime.now())


# ## 8.2 Save python and html versions

# In[36]:


nb_fname = '2_Clusters'
nb_fname


# In[37]:


get_ipython().run_cell_magic('bash', '-s "$nb_fname"', 'jupyter nbconvert "$1".ipynb --to="python"\njupyter nbconvert "$1".ipynb --to="html"\n')


# In[ ]:




