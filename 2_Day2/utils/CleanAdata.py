import scanpy as sc


def CleanAdata(adata, obstokeep=[], vartokeep=[], obsmtokeep=[]):
	_adata = adata.copy()
	for i in [i for i in list(_adata.uns.keys()) if "_colors" not in i]:
		del _adata.uns[i]
	
	for i in [i for i in _adata.obs.columns.tolist() if i not in obstokeep]:
		del _adata.obs[i]
	
	for i in [i for i in _adata.var.columns.tolist() if i not in vartokeep]:
		del _adata.var[i]

	for i in [i for i in _adata.obsm if i not in obsmtokeep]:
		del _adata.obsm[i]
 
	if "leiden" in _adata.obs.columns:
		del _adata.obs["leiden"]
	if "louvain" in _adata.obs.columns:
		del _adata.obs["louvain"]
	if "highly_variable" in _adata.var.columns:
		del _adata.var["highly_variable"]
	del _adata.varm
	del _adata.obsp 
	#del _adata.layers
	#del _adata.obsm
	return _adata
