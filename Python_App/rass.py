import tkinter as tk
from tkinter import *
from tkinter import ttk;
import matplotlib.pyplot as plt
from tkinter import messagebox
import scrublet as scr
import anndata
from PIL import Image
import io
from contextlib import redirect_stdout
import os
import time,scanpy as sc, anndata as ann
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import numpy as np
import random
import pandas as pd
import os
import subprocess
from subprocess import run
from matplotlib.figure import Figure
# OO backend (Tkinter) tkagg() function:
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import scrublet as scr
import Fonctions_tab1 as Ftab1
from tkinter.font import BOLD, Font
from itertools import chain
import gseapy





'''adata=anndata.read_h5ad('/home/izem/PycharmProjects/Explore/end_adata.h5ad')

adata.uns['log1p']['base']=None
#pd.DataFrame(adata.uns['rank_genes_groups']).to_csv("/home/izem/before.csv")
#csv_builder(adata,"all","leiden_1.0","before")


pio=["2","1"] #IndexError

x="3" #ValueError

sc.tl.rank_genes_groups(adata, 'leiden_1.0', groups=pio, reference= x, method='t-test')
csv_builder(adata,'all',0,'leiden_1.0',1,pio)'''

#print(adata.uns['rank_genes_groups']['scores'][0])

import tkinter as tk
from tkinter.filedialog import askopenfilename


# defining open_file_chooser function

glist=['IGKV4-1', 'CD55', 'IGKC', 'PPFIBP1', 'ABHD4', 'PCSK6', 'PGD', 'ARHGDIB', 'ITGB2', 'CARD6']


table=pd.read_csv("/home/izem/PycharmProjects/Explore/test.csv")
table.sort_values(by=table.columns[1], inplace=True, ascending=False)
table


tabla=pd.read_csv("/home/izem/Documents/test2.csv",header=None)

enr = gseapy.enrichr(gene_list=glist,
                 gene_sets='KEGG_2016',
                 background=tabla,
                 organism='Human', # don't forget to set organism to the one you desired! e.g. Yeast
                 description='test_name',
                 outdir='test/enrichr_kegg',
                 no_plot=False,
                 cutoff=0.5 # test dataset, use lower value from range(0,1)
                )

fig, ax = plt.subplots(1)
fig.tight_layout(w_pad=5, h_pad=5)
plt.subplots_adjust(left=0.05, bottom=0.1)

gseapy.dotplot(enr.res2d, title='KEGG_2013',cmap='viridis_r',ofname="resut_dot.png",top_term=15)

