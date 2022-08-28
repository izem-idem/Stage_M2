import gseapy
import matplotlib.pyplot
import pandas as pd , numpy as np , scanpy as sc , statistics as stat
#import mnnpy
from traceback import format_exc
import scanorama
import scrublet as scr
import matplotlib as plt
import tkinter as tk
from tkinter import *
from tkinter import ttk, simpledialog
from tkinter import messagebox
import scrublet as scr
import anndata
from PIL import Image, ImageTk
import io
from contextlib import redirect_stdout
import os,time
import seaborn as sb
import math
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from datetime import datetime
import matplotlib.pyplot as plt
import loompy

'''percent_mt=10
min_genes=150
min_cells=3
adata=sc.read_10x_mtx("/home/izem/Desktop/STAGE/Izem_Single-cell/muscle_1_ctrl/outs/filtered_feature_bc_matrix/",var_names="gene_symbols",cache=True)
adata.obs['sample']="izem"

bdata=sc.read_10x_mtx("/home/izem/Desktop/STAGE/Izem_Single-cell/muscle_2_ctrl/outs/filtered_feature_bc_matrix/",var_names="gene_symbols",cache=True)
bdata.obs['sample']="hichem"

adata=adata.concatenate(bdata)


adata.var['mt'] = adata.var_names.str.startswith('MT-')
adata.var[adata.var.mt == True]
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
# ngenes_by_count is unclear and is different from n_genes + no answer from scanpy regarding this issue

adata = adata[adata.obs.pct_counts_mt < percent_mt]
sc.pp.filter_cells(adata, min_genes=min_genes)
sc.pp.filter_genes(adata, min_cells=min_cells)

adata.raw = adata  # save raw adata
sc.pp.normalize_total(adata, target_sum=1e4)  # normalize every cell with 10 000 UMIs
#print(adata.X[1:4,1015:1025].todense())
sc.pp.log1p(adata)  # change it to logcounts
#print("-----------------------------------------------------------------------------")
#print(adata.X[1:4,1015:1025].todense())
sc.pp.scale(adata)

####imagebox.grid_forget() #pour effacer l'image mais ça s'affiche pas

#################LET'S PASS TO HGVs ###########################################################
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5, batch_key="sample")
adata = adata[:, adata.var.highly_variable]
#print(adata.var.highly_variable)
var_genes = adata.var.highly_variable  # get the HVGs
plt.show(block=False)
sc.pl.highly_variable_genes(adata) #je vais tester les 2 plots
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')
sc.pl.pca_loadings(adata, components=[1,2,3])
print("ca va ?")
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=40)
sc.tl.louvain(adata, resolution=1.0, key_added="louvain_1.0")
sc.tl.louvain(adata, resolution=0.5, key_added="louvain_0.5")
print(adata)
sc.tl.umap(adata)
sc.tl.tsne(adata)
print(adata)'''
#IsADirectoryError
#toto=sc.read_h5ad("/home/izem/PycharmProjects/Explore/end_adata.h5ad")

#print(toto.obsm['X_umap'].shape)
#sc.pl.umap(toto,color='leiden_1.0')


'''adata=sc.read_10x_mtx("/home/izem/PycharmProjects/contexte/raw_feature_bc_matrix_2")
print(adata)
adata.write_h5ad("human_ctrl.h5ad")'''


if not os.path.exists("/home/izem"):
    print("ok")




'''adata=sc.read_loom("SCG71.loom")


percent_mt=10
min_genes=150
min_cells=3

adata=sc.read_10x_mtx("/home/izem/Desktop/STAGE/Izem_Single-cell/muscle_1_ctrl/outs/filtered_feature_bc_matrix/",var_names="gene_symbols",cache=True)
adata.obs['sample']="izem"

bdata=sc.read_10x_mtx("/home/izem/Desktop/STAGE/Izem_Single-cell/muscle_2_ctrl/outs/filtered_feature_bc_matrix/",var_names="gene_symbols",cache=True)
bdata.obs['sample']="hichem"

adata=adata.concatenate(bdata)

adata.var['mt'] = adata.var_names.str.startswith('MT-')
adata.var[adata.var.mt == True]
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
# ngenes_by_count is unclear and is different from n_genes + no answer from scanpy regarding this issue
adata = adata[adata.obs.pct_counts_mt < percent_mt]
sc.pp.filter_cells(adata, min_genes=min_genes)
sc.pp.filter_genes(adata, min_cells=min_cells)


scrub = scr.Scrublet(adata.X)
adata.obs['doublet_scores'], adata.obs['predicted_doublets'] = scrub.scrub_doublets()

adata.raw = adata  # save raw adata
sc.pp.normalize_total(adata, target_sum=1e4)  # normalize every cell with 10 000 UMIs
#print(adata.X[1:4,1015:1025].todense())
sc.pp.log1p(adata)  # change it to logcounts
#print("-----------------------------------------------------------------------------")
#print(adata.X[1:4,1015:1025].todense())
sc.pp.scale(adata)

####imagebox.grid_forget() #pour effacer l'image mais ça s'affiche pas

#################LET'S PASS TO HGVs ###########################################################
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5, batch_key="sample")
#adata = adata[:, adata.var.highly_variable]
#print(adata.var.highly_variable)
var_genes = adata.var.highly_variable  # get the HVGs
sc.pl.highly_variable_genes(adata) #je vais tester les 2 plots
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')
sc.pl.pca_loadings(adata, components=[1,2,3])
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=40)
sc.tl.louvain(adata, resolution=1.0, key_added="louvain_1.0")
sc.tl.louvain(adata, resolution=0.5, key_added="louvain_0.5")
sc.tl.umap(adata)
sc.tl.tsne(adata)
#sc.pl.umap(adata,color=['louvain_1.0','louvain_0.5'])
#plt.show()

sc.tl.rank_genes_groups(adata, "louvain_1.0",use_raw=True)
print(type(sc.pl.rank_genes_groups_matrixplot(adata, n_genes=5, groupby="louvain_1.0",show=False)))
plt.savefig("ubuntu.png")
'''


'''adata=sc.read_h5ad("before_marker.h5ad")

print(adata.obsm)'''


