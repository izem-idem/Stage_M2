import tkinter as tk
from tkinter import *
from tkinter import  ttk
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
import Main

#in this file i have all the work before i split up into many python file to lighten my scripts


okey_image="/home/izem/Pictures/done.png"
wrong_image="/home/izem/Pictures/wrong.png"


os.chdir("/home/izem") #in order to find the file cell_cycle in my directory


# <editor-fold desc="FUNCTIONS">

dico_path={} #to stock multiple sample

adata="will contain anndata"
dico_adata={} #will contain the many samples before i will concatenate
list_concat=[] #i will this list to concatenate the datas !!

global_varpath = "will contain one sample adata"
global_varname= "will contain the name of the unique sample data"


def check_path(tab):
    if radio_string.get() == "":
        messagebox.showerror("Input Error", "Error: Please precise if you have one sample or multi-sample")

    if samplename_count_tab1.get() == "" or fullpath_count_tab1.get() == "":
        messagebox.showerror("Input Error", "Error: Please You have to fill the two forms before going on ! ")

    # I have to precise to the user to not use blank space in sample name
    if radio_string.get() == "un" and not samplename_count_tab1.get() == "" and not fullpath_count_tab1.get() == "":
        global global_varpath, global_varname, adata
        var_path = fullpath_count_tab1.get().strip(" ")
        var_name = samplename_count_tab1.get()  # no need to do .strip because it's free from space.
        try:
            adata = sc.read_10x_mtx(var_path, var_names="gene_symbols", cache=True)
            messagebox.showinfo("Information Count Matrix", "Your sample was successfully loaded !")
        except FileNotFoundError:
            messagebox.showerror("File Error",
                                 "Error: The path you entered is not a directory\nOR\nThe path you entered does not contain the required files")
    if radio_string.get() == "deux" and not samplename_count_tab1.get() == "" and not fullpath_count_tab1.get() == "":
        global dico_path, list_concat
        dico_path = {}
        var_path = fullpath_count_tab1.get().strip(" ")
        var_name = samplename_count_tab1.get()  # no need to do .strip because it's free from space.
        try:
            sc.read_10x_mtx(var_path, var_names="gene_symbols", cache=True)
            if var_name not in dico_path.keys() and var_path not in dico_path.values():  # to avoid having the same sample name or the same sample (ie avoiding same fullpath)
                dico_path[var_name] = var_path
                messagebox.showinfo("Information Count Matrix",
                                    "Path entered correctly !!\nPlease add another path to add another sample\n\nYou have currently " + str(
                                        len(dico_path.keys())) + " sample(s) that are registered")
            elif var_name in dico_path.keys() or var_path in dico_path.values():  # print error of redundant values
                messagebox.showerror("Input Error",
                                     "Error: the full path or the sample name you entered are already used !!")
        except FileNotFoundError:
            messagebox.showerror("File Error",
                                 "Error: The path you entered is not a directory\nOR\nThe path you entered does not contain the required files")


def show_image(imagefile, tab, col, ro):
    global image, imagebox
    image = ImageTk.PhotoImage(file=imagefile)
    imagebox = ttk.Label(tab)
    imagebox.config(image=image)
    imagebox.image = image
    imagebox.grid(column=col, row=ro)


def log_window(title, var,x=100,y=50):
    if radio_report.get() == "non":
        print("report hided !")
    else:
        global pop
        #x=random.randint(1,root.winfo_screenwidth())
        #y=random.randint(1,root.winfo_screenheight())
        pop = Toplevel(root)
        pop.title(title)
        pop.geometry('+%d+%d' % (x, y))
        # pop.config(background="green")
        Label(pop, text=var).pack()


def plot_window(title, plot,single=True,x=450,y=50):
    if radio_report.get() == "non" :
        print("report hided !")
    else:
        global wind
        #x=random.randint(1,root.winfo_screenwidth())
        #y=random.randint(1,root.winfo_screenheight())
        wind = Toplevel(root)
        wind.title(title)
        wind.geometry('+%d+%d' % (x, y)) #i want to determine where my window pops up
        canvas = FigureCanvasTkAgg(plot, wind)
        if single ==True : #if i have only one plot so i will have to fix the size of the window
            canvas.get_tk_widget().config(width=700, height=400)
            canvas.get_tk_widget().pack()
        elif single == False or single == None: #if i have multiple plots or none (meaning i have doublet histogram) it automatically scales
            canvas.get_tk_widget().pack()

def png_window(title,name,w=300,h=300,x=100,y=50): #for the scanpy plot function who don't return figures !!
    if radio_report.get() == "non" :
        print("report hided !")
    else:
        global taga
        #x = random.randint(1, root.winfo_screenwidth()) #to popup randomly
        #y = random.randint(1, root.winfo_screenheight()) #to popup randomly
        taga=Toplevel(root)
        taga.title(title)
        taga.geometry('+%d+%d' % (x, y))
        image = Image.open(os.getcwd()+"/figures/"+name)
        resize_img = image.resize((w, h))
        img = ImageTk.PhotoImage(resize_img)
        disp_img = Label(taga)
        disp_img.config(image=img)
        disp_img.image = img
        disp_img.pack(pady=20)


def Analysis_Gate():
    if radio_string.get() == "un":  # One Sample
        print(global_varpath)
        print(global_varname)
        adata = sc.read_10x_mtx(global_varpath, var_names="gene_symbols", cache=True)
        adata.obs['sample'] = global_varname
        Analysis_Room(adata, percent_mt=5, min_genes=150,min_cells=3)  # je lance la vraie analyse (pour l'instant je code l'analyse rapide mais apr??s faudra passer des param??tres pour l'analyse d??taill??)

    else:  # Multi Sample
        for key in dico_path.keys():
            dico_adata[key] = sc.read_10x_mtx(dico_path[key], var_names="gene_symbols", cache=True)
            dico_adata[key].obs['sample'] = key  # name the samples
            list_concat.append(dico_adata[key])

        adata = list_concat[0].concatenate[list_concat[1:]]  # avec le 1er adata je concatene tous les autres adatas

        Analysis_Room(adata, percent_mt=5, min_genes=150, min_cells=3)  # je lance la vraie analyse , avec un adata unifi?? !!


def Analysis_Room(adata, percent_mt=5, min_genes=150, min_cells=3):
    ###show_image(okey_image,tab1,3,2)
    adata.var  # rowdata equivalent
    extract_qc_report(adata,percent_mt=percent_mt,min_cells=min_cells,min_genes=min_genes)
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    adata.var[adata.var.mt == True]
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    # ngenes_by_count is unclear and is different from n_genes + no answer from scanpy regarding this issue
    adata = adata[adata.obs.pct_counts_mt < percent_mt]
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)
    adata.raw = adata  # save raw adata
    sc.pp.normalize_total(adata, target_sum=1e4)  # normalize every cell with 10 000 UMIs
    sc.pp.log1p(adata)  # change it to logcounts
    sc.pp.scale(adata)

    cell_cycle_genes = [x.strip() for x in open('cell_cycle_genes.txt')]  # I have to have this file in my dir
    # print(len(cell_cycle_genes))
    # Split into 2 lists
    s_genes = cell_cycle_genes[:43]
    g2m_genes = cell_cycle_genes[43:]
    cell_cycle_genes = [x for x in cell_cycle_genes if x in adata.var_names]
    # print(len(cell_cycle_genes))

    sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)

    f = io.StringIO()
    with redirect_stdout(f):
        scrub = scr.Scrublet(adata.raw.X)
        adata.obs['doublet_scores'], adata.obs['predicted_doublets'] = scrub.scrub_doublets()  # c'est cette ligne qui produit l'output
        s = f.getvalue()
    log_window("Doublet Detection report", s)
    blot=scrub.plot_histogram()[0]
    plot_window("Doublet Detection plot", blot)
    sum(adata.obs['predicted_doublets'])
    # add in column with singlet/doublet instead of True/False
    adata.obs['doublet_info'] = adata.obs["predicted_doublets"].astype(str)
    adata
    adata_clone = adata  # I clone it in order to use it for Combat batch correction later
    ####imagebox.grid_forget() #pour effacer l'image mais ??a s'affiche pas

    #################LET'S PASS TO HGVs ###########################################################

    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5, batch_key="sample")
    adata = adata[:, adata.var.highly_variable]
    #print(adata.var.highly_variable)
    var_genes = adata.var.highly_variable  # get the HVGs
    sc.pl.highly_variable_genes(adata,save=".png",show=False) #je vais tester les 2 plots
    png_window("Report HGVs","filter_genes_dispersion.png",w=850,h=500,x=10,y=450)
    sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, svd_solver='arpack')
    #sc.pl.pca_loadings(adata, components=[1,2,3,4,5,6,7,8])
    sc.pl.pca_variance_ratio(adata, log=True,show=False,save=".png") #faut que je teste show = True si ??a me fait sortir les autres.
    png_window("Report Pca","pca_variance_ratio.png",w=400,h=400,x=900,y=450) #the figure is saved so i can call the function to display it !
    ######################### PLOT HEAT MAP #############################################################
    '''genes = adata.var['gene_ids']

    for pc in [1, 2, 3, 4]:
        g = adata.varm['PCs'][:, pc - 1]
        o = np.argsort(g)
        sel = np.concatenate((o[:10], o[-10:])).tolist()
        emb = adata.obsm['X_pca'][:, pc - 1]
        # order by position on that pc
        tempdata = adata[np.argsort(emb),]
        sc.pl.heatmap(tempdata, var_names=genes[sel].index.tolist(), groupby='predicted_doublets', swap_axes=True,
                      use_raw=False)'''
    ######################### PLOT HEAT MAP #############################################################
    #ploto=sc.pl.umap(adata, color=['sample','sample'], title="UMAP",return_fig=True)
    #plot_window("Umap Embedding",ploto,single=False)

    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=40)

    if radio_string.get() == "un": #if we have one cluster let's pass to the clustering
        sc.tl.louvain(adata, resolution=1.0, key_added="louvain")
        sc.tl.umap(adata)
        sc.tl.tsne(adata)
        sc.tl.rank_genes_groups(adata, 'louvain', method='t-test')
        sc.pl.rank_genes_groups(adata, n_genes=20, sharey=False,show=False,save=".png") # je vais voir comment on g??re ce plot !
        extract_marker_report(adata,"Top_100_genes_by_clusters.csv",os.getcwd()+"/figures/rank_genes_groups_louvain.png",x=1350,y=450)
        if radio_report.get() == "non":
            print("report hided !")
        else:
            im = Image.open(os.getcwd() + "/figures/rank_genes_groups_louvain.png") #to display the png !
            im.show()
    else: #so it's multisample !
        # create a new object with lognormalized counts
        adata_combat = sc.AnnData(X=adata_clone.X, var=adata_clone.var, obs=adata_clone.obs)
        # first store the raw data
        adata_combat.raw = adata_combat
        # Run Combat
        sc.pp.combat(adata_combat, key='sample')
        sc.pp.highly_variable_genes(adata_combat)
        sc.pl.highly_variable_genes(adata_combat,save=".png",show=False)
        png_window("Report HGVs (after Integration)", "filter_genes_dispersion.png", w=850, h=500,x=10,y=500)
        sc.pp.pca(adata_combat, n_comps=30, use_highly_variable=True, svd_solver='arpack')
        sc.pp.neighbors(adata_combat, n_pcs=30)
        sc.tl.umap(adata_combat)
        sc.tl.tsne(adata_combat, n_pcs=30)
        sc.tl.louvain(adata_combat, key_added="louvain",resolution=1.0)
        #tmp = pd.crosstab(adata_combat.obs['leiden_1.0'], adata_combat.obs['sample'], normalize='index')
        #tmp.plot.bar(stacked=True).legend(loc='upper right')
        sc.tl.rank_genes_groups(adata_combat, 'louvain', method='t-test')
        sc.pl.rank_genes_groups(adata_combat, n_genes=20, sharey=False,save=".png",show=False)
        extract_marker_report(adata_combat,"Top_100_genes_by_clusters.csv",os.getcwd()+"/figures/rank_genes_groups_louvain.png",x=1350,y=450)
        if radio_report.get() == "non":
            print("report hided !")
        else:
            im = Image.open(os.getcwd() + "/figures/rank_genes_groups_louvain.png") #to display the png !
            im.show()

def extract_marker_report(bdata,filename,png_name,x,y):
    nb=len(bdata.obs.louvain.unique()) #nb of clusters
    pd.DataFrame(bdata.uns['rank_genes_groups']['names']).head(100).to_csv(filename,index=False) #get the first 100 top genes
    txt = "--------------------------"
    report = "{:s}\n The Marker gene detection plot is available : \n {:s} \n Number of clusters detected: {:d} \n Top 100 for each marker File: {:s} \n Comparaison method : t-test \n Correction method : Benjamini-Hocheberg \n p-value Cut-Off : 0.05 \n Clustering method : Louvain (Graph-based) \n{:s}".format(txt,png_name,nb,filename,txt)
    log_window("Marker detection report",report,x,y)


def extract_qc_report(adata, percent_mt, min_cells, min_genes):
    bdata = sc.AnnData(X=adata.X, obs=adata.obs, var=adata.var)
    begin_cells = len(bdata.obs);begin_genes = len(bdata.var)
    cdata = sc.AnnData(X=adata.X, obs=adata.obs,var=adata.var)  # i could've worked with adata rather than creating a new var but i'm afraid that it modify the real adata variable
    sc.pp.filter_cells(bdata, min_genes=min_genes);diff_cells = begin_cells - len(bdata.obs);sc.pp.filter_genes(bdata, min_cells=3);
    diff_genes = begin_genes - len(bdata.var)
    cdata.var  # rowdata equivalent
    cdata.var['mt'] = cdata.var_names.str.startswith('MT-')
    cdata.var[cdata.var.mt == True]
    sc.pp.calculate_qc_metrics(cdata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    cdata = cdata[cdata.obs.pct_counts_mt <= 10];diff_cells_mt = begin_cells - len(cdata.obs)
    sc.pp.filter_cells(cdata, min_genes=min_genes);sc.pp.filter_genes(cdata, min_cells=3);total_cells = len(cdata.obs);total_genes = len(cdata.var)
    txt = "--------------------------"
    report = "{:s}\n Number of cells that express less than {:d} genes: {:d} \n Number of genes that are expressed by less than {:d} cells : {:d} \n Number of cells that express more than {:d}% of mitonchondrial genes : {:d} \n Total number of cells : {:d} \n Total number of genes : {:d} \n {:s}".format(txt, min_genes, diff_cells, min_cells, diff_genes, percent_mt, diff_cells_mt, total_cells, total_genes, txt)
    log_window("QC report",report)

#</editor-fold>

root = tk.Tk()
root.title("Tab Widget")

root.geometry("600x180")

tabControl = ttk.Notebook(root)

tab1 = ttk.Frame(tabControl)
tab2 = ttk.Frame(tabControl)

tabControl.add(tab1, text='Analyse Rapide')
tabControl.add(tab2, text='Analyse Appronfondie')
tabControl.pack(expand=1, fill="both")




# <editor-fold desc="TAB ONE">


Label(tab1, text="COUNT MATRIX").grid(column=0,row=0,sticky="NW")

fullpath_count_tab1 = tk.StringVar()
samplename_count_tab1=tk.StringVar()


fullpath_entry_tab1=Entry(tab1, textvariable=fullpath_count_tab1).grid(column=1, row=0, sticky="NW")
name_entry_tab1=Entry(tab1, textvariable=samplename_count_tab1).grid(column=1, row=1, sticky="NW")

fullpath_count_tab1.set("enter fullpath")
samplename_count_tab1.set("enter sample name")

ok_count_button_tab1 = ttk.Button(tab1, text="OKEY", command=lambda: check_path(tab1)).place(x=280, y=13) # j'ai utiliser .place mais je vais peut-??tre le regretter


radio_string=tk.StringVar()
Label(tab1,text="Analysis Mode :").grid(column=0,row=2,sticky="W")
radio_unisample=ttk.Radiobutton(tab1,text="Uni-sample",value="un",variable=radio_string).grid(column=2,row=2,sticky="W",ipadx=0)
radio_multisample=ttk.Radiobutton(tab1,text="Multi-sample",value="deux",variable=radio_string).grid(column=1,row=2,sticky="W")
Label(tab1,text="Show Reports :").grid(column=0,row=3,sticky="W")
radio_report=tk.StringVar()

radio_yes=ttk.Radiobutton(tab1,text="YES",value="oui",variable=radio_report).grid(column=2,row=3,sticky="W")
radio_no=ttk.Radiobutton(tab1,text="NO",value="non",variable=radio_report).grid(column=1,row=3,sticky="W")


tab1.columnconfigure(5,weight=1)
Lunch_tab1=ttk.Button(tab1, text="Launch Analysis",command=lambda : Analysis_Gate()).grid(column=5,row=10,sticky="SE")

# </editor-fold>




# <editor-fold desc="TAB TWO">
'''Label(tab2, text="COUNT MATRIX").grid(sticky="NW")

fullpath_count_tab2 = tk.StringVar()
Count_entry_tab2=Entry(tab2,textvariable=fullpath_count_tab2).grid(column=1,row=0)
ok_count_button_tab2 = ttk.Button(tab2, text="OKEY", command=check_path()).grid(column=2,row=0)


qc_tab2=ttk.Label(tab2,text="QUALITY CONTROL & FILTERING").grid(sticky="NW",padx=0,pady=30,column=0,row=0)

norm_tab2=ttk.Label(tab2,text="NORMALISATION").grid(sticky="NW",padx=0,pady=30,column=0,row=1)

feature_tab2=ttk.Label(tab2,text="FEATURE SELECTION").grid(sticky="NW",padx=0,pady=30)

dimreduc_tab2=ttk.Label(tab2,text="DIMENSIONALITY REDUCTION").grid(sticky="NW",padx=0,pady=30)

clustering_tab2=ttk.Label(tab2,text="CLUSTERING").grid(sticky="NW",padx=0,pady=30)

markergene_tab2=ttk.Label(tab2,text="MARKER GENE DETECTION").grid(sticky="NW",padx=0,pady=30)'''

# </editor-fold>




root.mainloop()

