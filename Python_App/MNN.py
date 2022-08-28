import os
import matplotlib.pyplot
import pandas as pd , numpy as np , scanpy as sc , statistics as stat
import mnnpy
from traceback import format_exc
import scanorama
import scrublet as scr
import matplotlib as plt
import tkinter as tk
from tkinter import *
from tkinter import  ttk
from tkinter import messagebox
import scrublet as scr
import anndata
from PIL import Image, ImageTk
import io
from contextlib import redirect_stdout
import os,time
import seaborn as sb
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.pyplot as plt
import sys





adata=anndata.read_h5ad(os.getcwd()+'/mnn.h5ad')

adata.obs['sample']=adata.obs['sample'].astype('category') #make it categorical
batches = adata.obs['sample'].cat.categories.tolist()
alldata = {}
for batch in batches:
    alldata[batch] = adata[adata.obs['sample'] == batch,]

listing=list(alldata.keys())

cdata = sc.external.pp.mnn_correct(xxx,batch_key = 'sample', save_raw = True) ## x will be replaced by the list chriki svd_dim=50

cor_data=cdata[0]

cor_data.write_h5ad(os.getcwd()+"/mnn.h5ad")



