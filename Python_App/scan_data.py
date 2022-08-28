import pandas , numpy as np , scanpy as sc , statistics as stat
import igraph as ig , louvain ,leidenalg
import logging
from traceback import format_exc

rodata="/home/izem/Desktop/STAGE/raw_feature_bc_matrix_mouse/"

try:
    odoto=sc.read_10x_mtx(rodata,var_names="gene_symbols",cache=True)


except FileNotFoundError:
    print(format_exc().split('\n')[-2]) #i want only the last message , the stack is not  important for biologist

# Multisample Analysis
#### In this analysis we're trying to see how to handle multi sample with scanpy

 # this is an earlier version of the dataset from the pbmc3k tutorial
adata_ref = sc.datasets.pbmc3k_processed()  # this is an earlier version of the dataset from the pbmc3k tutorial
adata_ref
adata = sc.datasets.pbmc68k_reduced()
print(adata_ref)

pasway="/home/izem/Desktop/STAGE/Izem_Single-cell/Patient_IMI/IMI_data_1/outs/filtered_feature_bc_matrix"

bassway="/home/izem/Desktop/STAGE/Izem_Single-cell/muscle_2_ctrl/outs/filtered_feature_bc_matrix/"

adata_isch=sc.read_10x_mtx(pasway,var_names='gene_symbols',cache=True)
adata_ctrl=sc.read_10x_mtx(bassway,var_names='gene_symbols',cache=True)

print(adata_isch)
print(adata_ctrl)

adata=adata_ctrl.concatenate(adata_isch)
adata

adata.var #rowdata equivalent
adata.var['mt']=adata.var_names.str.startswith('MT-')
adata.var[adata.var.mt==True]

sc.pp.filter_cells(adata,min_genes=150)
#the output n_genes is created here
# comment voir had les cellules qui expriment moins de 200
sc.pp.filter_genes(adata,min_cells=3)
#comment voir had les génes qui sont exprimés par trés peu de cellules
adata

sc.pp.calculate_qc_metrics(adata,qc_vars=['mt'],percent_top=None,log1p=False,inplace=True)
#ngenes_by_count is unclear and is different from n_genes + no answer from scanpy regarding this issue
adata=adata[adata.obs.pct_counts_mt< 10]
adata

sc.pp.normalize_total(adata,target_sum=1e4)# normalize every cell with 10 000 UMIs
sc.pp.log1p(adata) #change it to logcounts
adata

sc.pp.highly_variable_genes(adata)
sc.pp.highly_variable_genes(adata,min_mean=0.0125,max_mean=3,min_disp=0.5)

adata.var[adata.var.highly_variable]

adata.raw=adata #save raw adata
adata = adata [:, adata.var.highly_variable]
adata

sc.pp.regress_out(adata,['total_counts','pct_counts_mt'])

sc.pp.scale(adata,max_value=10)

sc.tl.pca(adata, svd_solver='arpack')

sc.pp.neighbors(adata, n_neighbors=20, n_pcs=40)

sc.tl.louvain(adata,flavor='vtraag')
#sc.tl.umap(adata,init_pos='louvain')
#sc.tl.paga(adata) #i need to run tl.leiden or tl.louvain before running tl.paga
#sc.pl.paga(adata, plot=False)  # remove `plot=False` if you want to see the coarse-grained graph
#sc.tl.umap(adata, init_pos='paga')
sc.tl.umap(adata)
sc.tl.leiden(adata)

sc.pl.umap(adata, color=['batch','louvain'],use_raw=True)

sc.tl.rank_genes_groups(adata, 'louvain', method='t-test')
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)