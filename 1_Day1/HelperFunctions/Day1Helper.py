### VERSION 0.2

### Modules

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import gseapy as gp
from scipy.sparse import csr_matrix, isspmatrix

#----------------------------------

### 1. densityQCs

def densityQCs(adataObj, hue=None, 
               min_counts=None, max_counts=None, 
               min_genes=None, max_genes=None,
               pct_counts_mito=None, pct_counts_ribo=None):   
    
    '''
    Returns density plots for the following QCs: 
        n_genes_by_counts (log10), total_counts (log10), pct_counts_mito (linear scale), pct_counts_ribo (linear scale).
    
    Parameters:
    - adataObj: the adata containing the QC metrics above specified. PS: at the moment, if the QC are not present, this raises an error. This point could be improved. 
    - hue: xxx. Default is None. 
    - min_counts
    - ** kwargs: I used it to be flexible in passing anything we want to filter on, arguments should be either:
        - valid arguments for sc.pp-filter_cells()
        - column names from adata.obs (and pass as corresponding value the label(s, list or tuple) that you want to filter)
        '''
    
    #Plot them in line so they take up less space
    fig, ax = plt.subplots(1, 4, figsize=(20,5))
    fig.tight_layout(pad=2)   #space between plots
    if hue != None:
        hue_s = adata.obs[hue].astype('string')
    else:
        hue_s = None

    ### Genes ---------------
    d1 = sns.kdeplot(np.log10(adataObj.obs['n_genes_by_counts']), fill=True, color='cornflowerblue', hue=hue_s, ax=ax[0])
    min_x, max_x = d1.get_xlim() 

    #Threshold lines and fill
    if min_genes != None:
        d1.axvline(np.log10(min_genes), 0, 1, c='red')  #set manually for chosen threshold
        d1.axvspan(min_x, np.log10(min_genes), alpha=0.2, color='red')
    if max_genes != None:
        d1.axvline(np.log10(max_genes), c='red')
        d1.axvspan(np.log10(max_genes), max_x, alpha=0.2, color='red')

    ### UMI ---------------
    d2 = sns.kdeplot(np.log10(adataObj.obs['total_counts']), fill=True, color='forestgreen', hue=hue_s, ax=ax[1])
    min_x, max_x = d2.get_xlim() 
        
    if min_counts != None:
        d2.axvline(np.log10(min_counts), 0, 1, c='red')  #set manually for chosen threshold
        d2.axvspan(min_x, np.log10(min_counts), alpha=0.2, color='red')
    if max_counts != None:
        d2.axvline(np.log10(max_counts), c='red')
        d2.axvspan(np.log10(max_counts), max_x, alpha=0.2, color='red')

    ### Mito % ---------------
    d3 = sns.kdeplot(adataObj.obs['pct_counts_mito'], fill=True, color='coral', hue=hue_s, ax=ax[2])
    min_x, max_x = d3.get_xlim() 

    #Threshold lines and fill
    if pct_counts_mito != None:
        d3.axvline(pct_counts_mito, 0, 1, c='red')  #set manually for chosen threshold
        d3.axvspan(pct_counts_mito, max_x, alpha=0.2, color='red')


    ### Ribo % ---------------
    d4 = sns.kdeplot(adataObj.obs['pct_counts_ribo'], fill=True, color='orchid', hue=hue_s, ax=ax[3])
    min_x, max_x = d4.get_xlim() 
    #ax[3].legend(loc='center left', bbox_to_anchor=(1.0, 1.0)) #upper right
    
    #Threshold lines and fill
    if pct_counts_ribo != None:
        d4.axvline(pct_counts_ribo, 0, 1, c='red')  #set manually for chosen threshold
        d4.axvspan(pct_counts_ribo, max_x, alpha=0.2, color='red')
    
    #Remove additional legends at need
    if hue != None:
        ax[0].get_legend().remove()
        ax[1].get_legend().remove()
        ax[2].get_legend().remove()
        
    # Remove all borders
    sns.despine(bottom = False, left = True)

#----------------------------------
    
### 2. filterCellBarplot
   
def filterCellBarplot(adataObj, **kwargs):
    
    '''
    Returns a barplot depicting cells present at the start of each spcified filtering step and how many cells get filtered.
    
    Parameters:
    - adataObj: the adata object to filter. The object original is NOT filtered (a copy is made), only the plot is returned.
    - ** kwargs: I used it to be flexible in passing anything we want to filter on, arguments should be either:
        - valid arguments for sc.pp-filter_cells()
        - column names from adata.obs (and pass as corresponding value the label(s, list or tuple) that you want to filter)
        '''
    
    sc.settings.verbosity = 1  #Don't show filtering messages
    
    #INIZIALIZE
    adata_filt = adataObj.copy()
    #Store the number of cells at each filtering step
    n_cells = [adata_filt.n_obs]
    #Compute how many cells are removed at each step
    removed = []
    #Store names for plot labels
    names = []   
    
    #FILTERING
    #I used kwargs and no named arguments to run everything in a for loop, maybe not too elegant
    #for each key-value preform filtering and append cell number and labels
    
    for key, value in kwargs.items(): 
        #Check exact arguments for sc.pp.filter_cells
        if key in ['min_counts', 'max_counts', 'min_genes', 'max_genes']:  
            sc.pp.filter_cells(adata_filt, **{key: value})
            n_cells.append(adata_filt.n_obs)
            names.append(key)    
        
        else:
            assert key in adataObj.obs.columns, 'Please specify valid adata.obs column names as arguments'
        
            #MORE FLEXIBLE: check for any argument containing mito/ribo
            #elif key in ['pct_counts_mito','pct_counts_ribo']:
            if 'pct_counts_mito' in key or 'pct_counts_ribo' in key:  
                 adata_filt = adata_filt[adata_filt.obs[key] < value, :]
                 n_cells.append(adata_filt.n_obs)
                 names.append(key)   
                
            else:
                #When we filter on more than one label in the same .obs column
                if type(value) is list or type(value) is tuple:
                    for el in value:
                        assert el in adataObj.obs[key].unique(), f"Please specify valid adata.obs['{key}'] values as parameter"
                        
                        adata_filt = adata_filt[adata_filt.obs[key] != el, :]
                        n_cells.append(adata_filt.n_obs)
                        names.append(el)   
                else:
                    assert value in adataObj.obs[key].unique(), f"Please specify valid adata.obs['{key}'] value as parameter"
                    
                    adata_filt = adata_filt[adata_filt.obs[key] != value, :]
                    n_cells.append(adata_filt.n_obs)
                    names.append(value) 
        
        
    #Update lists for plot
    names.append('end')
    removed = n_cells[1:] + [adata_filt.n_obs]
    #Use numpy arrays to subtract vectors
    rem = np.array(removed) - np.array(n_cells) #this order to have negative values
    
    #PLOT
    plt.figure(figsize=(20,8))
    plot1 = sns.barplot(x=names, y=n_cells, color='turquoise',  label = "starting cells")
    plot2 = sns.barplot(x=names, y=rem, color='crimson',  hatch='/',  label = "filtered cells").set(xlabel='Filtering Step', ylabel='Number of cells')
    plt.axhline(0, 0, 1, c='black') 
    plt.legend(frameon = False)
    
    # Annotate the bars
    for p in plot1.patches[:len(names)]:
        plot1.annotate(format(p.get_height(), '.0f'), 
                       (p.get_x() + p.get_width() / 2., p.get_height()), 
                       ha = 'center', va = 'center', 
                       xytext = (0, 10), 
                       textcoords = 'offset points') 
        
    for p in plot1.patches[len(names):]:
        plot1.annotate(format(p.get_height(), '.0f'), 
                       (p.get_x() + p.get_width() / 2., p.get_height()), 
                       ha = 'center', va = 'center', 
                       xytext = (0, -10), 
                       textcoords = 'offset points') 
    
    sc.settings.verbosity = 3
    plt.show()    

#----------------------------------

### 3. selectMarkers

def selectMarkers(adataObj, mList):
    
    """  
    From a list of gene names select only the genes that are present in adata.var
    """
    
    #Select markers present in adata
    p = adataObj.var_names[adataObj.var_names.isin(mList) == True]
    #Keep the same order as input list
    p = [x for x in mList if x in p]   
    
    #Select missing genes
    ab = set(mList).difference(set(adataObj.var_names))
    
    #Print message 
    if len(ab) == len(mList):
        print('\nAll markers are missing')
    else:
        print('\nThe following marker genes are missing: ', ab)
        
    return(p)

#----------------------------------
    
### 4. customUmap

def customUmap(adata, genes, size=10, ncols=3):
    genes = selectMarkers(adata, genes) 
    sc.pl.embedding(adata, basis='X_umap_harmony', color=genes, size=size, frameon=False,
               sort_order=False, vmin=0,  vmax='p99', layer='lognormcounts', ncols=ncols)


#----------------------------------
    
### 5. customGseapy

def customGseapy(adata, cluster, rank, 
                 sets = ['GO_Biological_Process_2023', 'CellMarker_Augmented_2021'],
                 fdr_th=0.005, nes_th=1.75, show=20):
    
    """  
    GSEA analysis with gseapy for cluster marker genes ranked by score.
    Custom wrapper function.
    """
    
    # create ranked list
    ranked_list = sc.get.rank_genes_groups_df(adata, cluster, key=rank).drop(columns=['logfoldchanges','pvals','pvals_adj',
                                     'pct_nz_group','pct_nz_reference']).set_index('names') 
    
    # perform gsea analysis 
    pre_res = gp.prerank(rnk=ranked_list, gene_sets=sets, 
                     min_size=20, max_size=250, ascending=False,
                     permutation_num=200, outdir=None, # don't write to disk
                     seed=6, verbose=False)
    
    # extract and organize results
    res = pd.DataFrame(pre_res.results).T.sort_values('nes', ascending=False).drop(columns=['name', 'matched_genes', 'hits', 'RES', 'es'])
    res = res[res.fdr < fdr_th] 
    res = res[res.nes > nes_th]                
    return res.head(show)

