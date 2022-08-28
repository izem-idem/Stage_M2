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



'''
adata=sc.read_10x_mtx("/home/izem/Desktop/Data_2/MTX_CTRL_TCDD/outs/raw_feature_bc_matrix/",var_names="gene_symbols",cache=True)
#adata=sc.read_10x_mtx("/home/izem/Desktop/Data_2/MTX_CTRL_TCDD/outs/filtered_feature_bc_matrix/",var_names="gene_symbols",cache=True)
adata.var['mt'] = adata.var_names.str.startswith('MT-')
adata.var[adata.var.mt == True]
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

bdata=sc.read_10x_mtx("/home/izem/Desktop/Data_2/MTX_CTRL_Nonane/outs/raw_feature_bc_matrix/",var_names="gene_symbols",cache=True)
#adata=sc.read_10x_mtx("/home/izem/Desktop/Data_2/MTX_CTRL_TCDD/outs/filtered_feature_bc_matrix/",var_names="gene_symbols",cache=True)
bdata.var['mt'] = bdata.var_names.str.startswith('MT-')
bdata.var[bdata.var.mt == True]
sc.pp.calculate_qc_metrics(bdata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

cdata=sc.read_10x_mtx("/home/izem/Desktop/Data_2/MTX_cKO_TCDD/outs/raw_feature_bc_matrix/",var_names="gene_symbols",cache=True)
#adata=sc.read_10x_mtx("/home/izem/Desktop/Data_2/MTX_CTRL_TCDD/outs/filtered_feature_bc_matrix/",var_names="gene_symbols",cache=True)
cdata.var['mt'] = cdata.var_names.str.startswith('MT-')
cdata.var[cdata.var.mt == True]
sc.pp.calculate_qc_metrics(cdata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

ddata=sc.read_10x_mtx("/home/izem/Desktop/Data_2/MTX_cKO_Nonane/outs/raw_feature_bc_matrix/",var_names="gene_symbols",cache=True)
#adata=sc.read_10x_mtx("/home/izem/Desktop/Data_2/MTX_CTRL_TCDD/outs/filtered_feature_bc_matrix/",var_names="gene_symbols",cache=True)
ddata.var['mt'] = ddata.var_names.str.startswith('MT-')
ddata.var[ddata.var.mt == True]
sc.pp.calculate_qc_metrics(ddata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)






sc.pp.filter_cells(adata, min_genes=220)
#sc.pp.filter_cells(adata, max_genes=3000)
sc.pp.filter_cells(adata, min_counts=200)
sc.pp.filter_genes(adata, min_cells=3)
print("nbre de cellules : ",adata.n_obs,"nbre de génes:" ,adata.n_vars)
sc.pp.filter_cells(adata, min_genes=220)
#sc.pp.filter_cells(adata, max_genes=3000)
sc.pp.filter_cells(bdata, min_counts=200)
sc.pp.filter_genes(bdata, min_cells=3)
print("nbre de cellules : ",bdata.n_obs,"nbre de génes:" ,bdata.n_vars)
sc.pp.filter_cells(cdata, min_genes=220)
#sc.pp.filter_cells(adata, max_genes=3000)
sc.pp.filter_cells(cdata, min_counts=200)
sc.pp.filter_genes(cdata, min_cells=3)
print("nbre de cellules : ",cdata.n_obs,"nbre de génes:" ,cdata.n_vars)
sc.pp.filter_cells(ddata, min_genes=220)
#sc.pp.filter_cells(adata, max_genes=3000)
sc.pp.filter_cells(ddata, min_counts=200)
sc.pp.filter_genes(ddata, min_cells=3)
print("nbre de cellules : ",ddata.n_obs,"nbre de génes:" ,ddata.n_vars)'''

tcdd=anndata.read_h5ad("/home/izem/ctrl_tcdd.h5ad")
nonane=anndata.read_h5ad("/home/izem/ctrl_nonane.h5ad")
tcdd.obs['sample']="ctrl_tcdd"
nonane.obs['sample']="ctrl_nonane"

sc.pp.normalize_total(tcdd,target_sum=1e4)# normalize every cell with 10 000 UMIs
sc.pp.log1p(tcdd) #change it to logcounts
print(tcdd)
#save raw matrix
tcdd.raw=tcdd
sc.pp.normalize_total(nonane,target_sum=1e4)# normalize every cell with 10 000 UMIs
sc.pp.log1p(nonane) #change it to logcounts
print(nonane)
#save raw matrix
nonane.raw=nonane

adata=tcdd.concatenate(nonane)

adata

adata_combat = sc.AnnData(X=adata.raw.X, var=adata.var, obs = adata.obs)

# first store the raw data
adata_combat.raw = adata_combat

# run combat
sc.pp.combat(adata_combat, key='sample')

sc.pp.highly_variable_genes(adata_combat)

sc.tl.pca(adata_combat, svd_solver='arpack')

sc.pp.neighbors(adata_combat, n_neighbors=20, n_pcs=40)

sc.tl.louvain(adata_combat,resolution=1.0,key_added="louvain")
#sc.tl.umap(adata,init_pos='louvain')
#sc.tl.paga(adata) #i need to run tl.leiden or tl.louvain before running tl.paga
#sc.pl.paga(adata, plot=False)  # remove `plot=False` if you want to see the coarse-grained graph
#sc.tl.umap(adata, init_pos='paga')
sc.tl.leiden(adata_combat,resolution=1.0,key_added="leiden")
sc.tl.umap(adata_combat)

sc.pl.umap(adata_combat, color=['louvain'])