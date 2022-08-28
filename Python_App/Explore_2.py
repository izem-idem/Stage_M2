import tkinter as tk
from tkinter import *
from tkinter import ttk, Tk
from tkinter import messagebox
import scrublet as scr
import anndata
from PIL import Image, ImageTk
import io
from contextlib import redirect_stdout
import os,time,scanpy as sc, anndata as ann
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import numpy as np
import random
import pandas as pd
import os
from matplotlib.figure import Figure
# OO backend (Tkinter) tkagg() function:
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import scrublet as scr
from tkinter.font import BOLD, Font
import matplotlib.pyplot as plt
import Fonctions_tab1 as Ftab1
import seaborn as sb
from tkinter.filedialog import askopenfilename
from tkinter.scrolledtext import ScrolledText
import gseapy
from gseapy.plot import barplot, dotplot


'''
percent_mt=10
min_genes=150
min_cells=3
adata=sc.read_10x_mtx("/home/izem/Desktop/STAGE/Izem_Single-cell/muscle_1_ctrl/outs/filtered_feature_bc_matrix/",var_names="gene_symbols",cache=True)
adata.obs['sample']="izem"


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
#sc.pl.highly_variable_genes(adata,show=False,save=".png") #je vais tester les 2 plots
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')
#sc.pl.pca_loadings(adata, components=[1,2,3,4,5,6,7,8]).
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=40)
sc.tl.louvain(adata, resolution=1.0, key_added="louvain")

sc.tl.rank_genes_groups(adata, 'louvain', method='t-test')
#sc.pl.pca_variance_ratio(adata, log=True,save=".png",show=False)
#sc.pl.rank_genes_groups(adata, n_genes=20, sharey=False,save=".png",show=False)

'''
'''def resize_func():
    image = Image.open("/home/izem/figures/violin.png")
    w = 400
    h = 400
    resize_img = image.resize((w, h))
    img = ImageTk.PhotoImage(resize_img)
    disp_img=Label(root)
    disp_img.config(image=img)
    disp_img.image = img
    disp_img.place(x=5,y=130)

'''






'''ws = Tk()
ws.title('PythonGuides')
ws.geometry('+100+50')
ws.config(bg='#4a7a8c')
root=Toplevel(ws)
root.geometry('850x500+10+500')
taga=Toplevel(ws)
taga.geometry('850x500+10+450')


frame = Frame(ws)
frame.pack()

Label(
    frame,
    text='Width'
    ).pack(side=LEFT)
width = Entry(frame, width=10)
width.insert(END, 300)
width.pack(side=LEFT)

Label(
    frame,
    text='Height'
    ).pack(side=LEFT)

height = Entry(frame, width=10)
height.insert(END, 350)
height.pack(side=LEFT)

resize_btn = Button(
    frame,
    text='Resize',
    command=lambda :resize_func()
)
resize_btn.pack(side=LEFT)

disp_img = Label()
disp_img.pack(pady=20)


ws.mainloop()

'''

'''txt = "--------------------------"
report = "{:s}\n The Marker gene detection plot is available : \n {:s} \n Number of clusters detected: {:d} \n Top 100 for each marker File: {:s} \n Comparaison method : t-test \n Correction method : Benjamini-Hocheberg \n p-value Cut-Off : 0.05 \n Clustering method : Louvain (Graph-based) \n{:s}".format(txt,"pathway",12,"pathway deuxiéme",txt)


root = Tk()

figure = Figure(figsize=(8, 5), dpi=100)
plot = figure.add_subplot(2, 1, 1)

x = [ 0.1, 0.2, 0.3, 0.4 ]
y = [ -0.1, -0.2, -0.3, -0.4 ]
plot.plot(x, y, color="red", marker="o",  linestyle="--")

canvas = FigureCanvasTkAgg(figure, root)
canvas.get_tk_widget().grid(row=0, column=0)

Label(root,text=report).pack()

root.mainloop() #to display the plots
'''

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
#sc.pl.highly_variable_genes(adata,show=False,save=".png") #je vais tester les 2 plots
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')
#sc.pl.pca_loadings(adata, components=[1,2,3,4,5,6,7,8]).
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=40)
sc.tl.louvain(adata, resolution=1.0, key_added="louvain")

tmp=pd.crosstab(adata.obs['sample'],adata.obs['louvain']) #normalize='index' for the percentages.
#tmp.columns.categories
dico={}
for i in range(0,len(tmp)):
    dico[tmp.index.categories[i]]={}
    for j in range(0,len(tmp.columns.categories)):
        dico[tmp.index.categories[i]][j]=tmp.iloc[i,j]

goose = ""

gos=""
for i in tmp.columns.categories:
    gos=gos+"\nCluster {0} = {1} cells ".format(i,tmp.iloc[:,int(i)].sum())

gas=""
for i in dico.keys():
    gas=gas+"\n\nReport for sample {0}\n".format(i)
    for j in dico[i]:
        gas=gas+"\nCluster {0} = {1} cells ".format(j,dico[i][j])

giz="\t Report for the Clustering:\n\nClustering Method: {0}\n\nCluster number: {1}\n\nTotal number of cell by cluster:{2}\n\nReport for each Sample:{3}".format(tmp.columns.name,len(tmp.columns),gos,gas)
'''


def txt_summary(togo,organism):
    gene_set_names = gseapy.get_library_name(organism=organism)
    txt=""
    for x in gene_set_names:
        txt=txt+x+"\n"

    window = Toplevel(togo)
    window.title("Report Clustering")
    window.geometry("")
    text_widget = tk.Text(window, height=30, width=70)
    # Create a scrollbar
    scroll_bar = tk.Scrollbar(window)
    # Pack the scroll bar
    # Place it to the right side, using tk.RIGHT
    scroll_bar.pack(fill=tk.BOTH, side=tk.RIGHT, expand=True)
    # Pack it into our tkinter application
    # Place the text widget to the left side
    text_widget.pack(side=tk.LEFT)

    # Insert text into the text widget
    text_widget.insert(tk.END, txt)
    text_widget["state"] = DISABLED




root = tk.Tk()

target=tk.StringVar()

background=tk.StringVar()
background.set("Optional")
taga = Toplevel(root)
taga.title("Enrichr Analysis")
taga.geometry('+%d+%d' % (800, 750))

gene7=tk.StringVar()
rapport=tk.StringVar()
var_organ=tk.StringVar()
pivalue = tk.DoubleVar()
pivalue.set(0.05)


#function that "replaces" text with canvas
import tkinter as tk
import klembord

root = tk.Tk()
text = tk.Text(root)
text.pack(fill='both', expand=True)

text.tag_configure('italic', font='TkDefaultFont 9 italic')
text.tag_configure('bold', font='TkDefaultFont 9 bold')

TAG_TO_HTML = {
    ('tagon', 'italic'): '<i>',
    ('tagon', 'bold'): '<b>',
    ('tagoff', 'italic'): '</i>',
    ('tagoff', 'bold'): '</b>',
}

def copy_rich_text(event):
    try:
        txt = text.get('sel.first', 'sel.last')
    except tk.TclError:
        # no selection
        return "break"
    content = text.dump('sel.first', 'sel.last', tag=True, text=True)
    html_text = []
    for key, value, index in content:
        if key == "text":
            html_text.append(value)
        else:
            html_text.append(TAG_TO_HTML.get((key, value), ''))
    klembord.set_with_rich_text(txt, ''.join(html_text))
    return "break"  # prevent class binding to be triggered

text.bind('<Control-c>', copy_rich_text)

text.insert("1.0", "Author et al. (2012). The title of the article. ")
text.insert("end", "Journal Name", "italic")
text.insert("end", ", ")
text.insert("end", "2", "bold")
text.insert("end", "(599), 1–5.")

text["state"] = DISABLED

root.mainloop()