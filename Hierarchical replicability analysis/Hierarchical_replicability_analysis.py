import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import pymn
import os
seed=10
os.chdir("/.../aki_to_ckd/")
sc.logging.print_versions()
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
sns.set(style='white', font_scale=1.25)
plt.rc("axes.spines", top=False, right=False)
plt.rc('xtick', bottom=True)
plt.rc('ytick', left=True)




#========================================================================================================
#=========== Protocol 1: assessment of cell type replicability with unsupervised MetaNeighbor ===========
#========================================================================================================

#=========== #Step 2: Hierarchical cell type replicability analysis ===========

#load data
adata = sc.read("/.../merge_13k_shared_genes_x_641k_cells.h5ad")
adata
#AnnData object with n_obs × n_vars = 641212 × 13012
#    obs: 'orig.ident', 'nCount_RNA', 'nFeature_RNA', 'disease_l1', 'lab'
#    var: 'features'

print(adata.obs['disease_l1'].unique().tolist())
#['Control', 'AKI', 'Recovery', 'CKD']
print(adata.obs['lab'].unique().tolist())
#['balzer', 'doke', 'kirita', 'abedini', 'hinze', 'lake']

adata.obs['cell.type'] = adata.obs['disease_l1'].astype(str)
adata.obs['study_id'] = adata.obs['lab'].astype(str)
nan_count = pd.isna(adata.obs.disease_l1).sum()
nan_count #0

#get variable genes across datasets
pymn.variableGenes(adata, study_col='study_id')

#get hierarchical replicability AUROC scores
pymn.MetaNeighborUS(adata,
                    study_col='study_id',
                    ct_col='cell.type',
                    fast_version=True)

plt.figure()
dummy = ['#C00000', '#00B050', '#F00000', '#C00000', '#00B050', '#FFFF00', '#C00000', '#00B050', '#F00000', '#00B050', '#F00000', '#C00000', '#00B050', '#F00000', '#C00000', '#00B050']
name1 = list(adata.uns["MetaNeighborUS"])
col_colors1 = pd.Series(dummy, index=name1, name="lab")
pymn.plotMetaNeighborUS(adata, figsize=(14, 14.5), cmap='coolwarm', fontsize=22, 
                        col_colors=col_colors1, 
                        dendrogram_ratio=(0, .05),
                        cbar_pos=None) #tuple of (left, bottom, width, height), optional
plt.savefig('plotMetaNeighborUS3.svg')

adata.uns['MetaNeighborUS']
df = pd.DataFrame(adata.uns['MetaNeighborUS'])
df.to_csv("MN.csv")

pymn.topHits(adata, threshold=0.8)
adata.uns['MetaNeighborUS_topHits']
df = pd.DataFrame(adata.uns['MetaNeighborUS_topHits'])
df.to_csv("MN_topHits.csv")

#save outfile
adata.write('merge_13k_shared_genes_x_641k_cells_afterpymn.h5ad')


