import tkinter as tk
from itertools import chain
from tkinter.scrolledtext import ScrolledText
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
import matplotlib.pyplot as plt
from tkinter.font import BOLD, Font

place=os.getcwd()
okey_image="done.png"
wrong_image="wrong.png"


placa="/home/izem/Pictures/"

# <editor-fold desc="FUNCTIONS">

dico_path={} #to stock multiple sample

adata="will contain anndata"
dico_adata={} #will contain the many samples before i will concatenate
list_concat=[] #i will this list to concatenate the datas !!

var_path = "will contain one sample adata"
var_name= "will contain the name of the unique sample data"


okey_check=False

def check_path(radio_string,samplename_count_tab1,fullpath_count_tab1,radio_report,formate,output):
    if radio_string.get() == "" or radio_report.get() == "" or formate.get()=="":
        messagebox.showerror("Input Error", "Error: Please choose one option for the 3 propositions")

    if output.get()=="":
        messagebox.showerror("Input Error", "Error: Please provide a dir_name to store the output")
    if samplename_count_tab1.get() == "" or fullpath_count_tab1.get() == "":
        messagebox.showerror("Input Error", "Error: Please You have to fill the two forms before going on ! ")

    # I have to precise to the user to not use blank space in sample name
    if radio_string.get() == "un" and not samplename_count_tab1.get() == "" and not fullpath_count_tab1.get() == "" and not radio_report.get() == ""  and output.get()!="":
        os.system(f"mkdir {output.get().strip()}") #create output dir
        global var_path, var_name, adata,okey_check
        var_path = fullpath_count_tab1.get().strip(" ")
        var_name = samplename_count_tab1.get()  # no need to do .strip because it's free from space.
        okey_check=True
        if formate.get()=="10x_mtx":
            try:
                adata = sc.read_10x_mtx(var_path, var_names="gene_symbols", cache=True)
                messagebox.showinfo("Information Count Matrix", "Your sample was successfully loaded !")
            except FileNotFoundError:
                messagebox.showerror("File Error", "Error: The path you entered is not a directory\nOR\nThe path you entered does not contain the required files")
        elif formate.get()=="mtx":
            try:
                adata = sc.read_mtx(var_path)
                messagebox.showinfo("Information Count Matrix", "Your sample was successfully loaded !")
            except:
                messagebox.showerror("File Error", "Error: the file you provided does not fit")
        elif formate.get()=="h5ad":
            try:
                adata = sc.read_h5ad(var_path)
                messagebox.showinfo("Information Count Matrix", "Your sample was successfully loaded !")
            except:
                messagebox.showerror("File Error", "Error: the file you provided does not fit!")
        elif formate.get()=="loom":
            try:
                adata = sc.read_loom(var_path)
                messagebox.showinfo("Information Count Matrix", "Your sample was successfully loaded !")
            except:
                messagebox.showerror("File Error", "Error: the file you provided does not fit!")


    if radio_string.get() == "deux" and not samplename_count_tab1.get() == "" and not fullpath_count_tab1.get() == "" and output.get()!="":
        os.system(f"mkdir {output.get().strip()}")  # create output dir
        global dico_path, list_concat
        okey_check=True
        var_path = fullpath_count_tab1.get().strip(" ")
        var_name = samplename_count_tab1.get()  # no need to do .strip because it's free from space.
        if formate.get()=="10x_mtx":
            try:
                sc.read_10x_mtx(var_path, var_names="gene_symbols", cache=True)
                if var_name not in dico_path.keys() and var_path not in dico_path.values():  # to avoid having the same sample name or the same sample (ie avoiding same fullpath)
                    dico_path[var_name] = var_path
                    messagebox.showinfo("Information Count Matrix",
                                        "Path entered correctly !!\nPlease add another path to add another sample\n\nYou have currently " + str(
                                            len(dico_path.keys())) + " sample(s) that are registered")
                elif var_name in dico_path.keys() or var_path in dico_path.values():  # print error of redundant values
                    messagebox.showerror("Input Error","Error: the full path or the sample name you entered are already used !!")
            except :
                messagebox.showerror("File Error","Error: The file you provided does not fit")
        elif formate.get()=="mtx":
            try:
                sc.read_mtx(var_path)
                if var_name not in dico_path.keys() and var_path not in dico_path.values():  # to avoid having the same sample name or the same sample (ie avoiding same fullpath)
                    dico_path[var_name] = var_path
                    messagebox.showinfo("Information Count Matrix","Path entered correctly !!\nPlease add another path to add another sample\n\nYou have currently " + str(
                                            len(dico_path.keys())) + " sample(s) that are registered")
                elif var_name in dico_path.keys() or var_path in dico_path.values():  # print error of redundant values
                    messagebox.showerror("Input Error","Error: the full path or the sample name you entered are already used !!")
            except :
                messagebox.showerror("File Error","Error: The file you entered does not fit")
        elif formate.get()=="h5ad":
            try:
                sc.read_h5ad(var_path)
                if var_name not in dico_path.keys() and var_path not in dico_path.values():  # to avoid having the same sample name or the same sample (ie avoiding same fullpath)
                    dico_path[var_name] = var_path
                    messagebox.showinfo("Information Count Matrix","Path entered correctly !!\nPlease add another path to add another sample\n\nYou have currently " + str(
                                            len(dico_path.keys())) + " sample(s) that are registered")
                elif var_name in dico_path.keys() or var_path in dico_path.values():  # print error of redundant values
                    messagebox.showerror("Input Error","Error: the full path or the sample name you entered are already used !!")
            except :
                messagebox.showerror("File Error","Error: The file you entered does not fit")
        elif formate.get()=="loom":
            try:
                sc.read_loom(var_path)
                if var_name not in dico_path.keys() and var_path not in dico_path.values():  # to avoid having the same sample name or the same sample (ie avoiding same fullpath)
                    dico_path[var_name] = var_path
                    messagebox.showinfo("Information Count Matrix","Path entered correctly !!\nPlease add another path to add another sample\n\nYou have currently " + str(
                                            len(dico_path.keys())) + " sample(s) that are registered")
                elif var_name in dico_path.keys() or var_path in dico_path.values():  # print error of redundant values
                    messagebox.showerror("Input Error","Error: the full path or the sample name you entered are already used !!")
            except :
                messagebox.showerror("File Error","Error: The file you entered does not fit")



def show_image(imagefile, tab, col, ro):
    global image, imagebox
    image = ImageTk.PhotoImage(file=imagefile)
    imagebox = ttk.Label(tab)
    imagebox.config(image=image)
    imagebox.image = image
    imagebox.grid(column=col, row=ro)


def log_window(root,title, radio_report,var,x=100,y=50):
    if radio_report.get() == "non":
        print("report hided !")
    else:
        #global pop #no need to make it globally my friend
        #x=random.randint(1,root.winfo_screenwidth())
        #y=random.randint(1,root.winfo_screenheight())
        pop = Toplevel(root)
        pop.title(title)
        pop.geometry('+%d+%d' % (x, y))
        # pop.config(background="green")
        Label(pop, text=var).pack()


def plot_window(root,title, radio_report,plot,single=True,x=450,y=50):
    if str(type(radio_report))=="<class 'tkinter.StringVar'>":
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
            elif single == False or single == None: #if i have multiple plots or none (meaning i have doublet histogram or other things) it automatically scales
                canvas.get_tk_widget().pack()
    else:
        #x=random.randint(1,root.winfo_screenwidth())
        #y=random.randint(1,root.winfo_screenheight())
        wind = Toplevel(root)
        wind.title(title)
        wind.geometry('+%d+%d' % (x, y)) #i want to determine where my window pops up
        canvas = FigureCanvasTkAgg(plot, wind)
        if single ==True : #if i have only one plot so i will have to fix the size of the window
            canvas.get_tk_widget().config(width=700, height=400)
            canvas.get_tk_widget().pack()
        elif single == False or single == None: #if i have multiple plots or none (meaning i have doublet histogram or other things) it automatically scales
            canvas.get_tk_widget().pack()


def png_window(root,radio_report,title,name,w=300,h=300,x=100,y=50): #for the scanpy plot function who don't return figures !!
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


def Analysis_Gate(root,radio_string,radio_report,output):
    if okey_check==False:
        messagebox.showerror("Error","Please press Okey before continuing")

    if radio_string.get() == "un" and okey_check:  # One Sample
        print(var_path)
        print(var_name)

        adata = sc.read_10x_mtx(var_path, var_names="gene_symbols", cache=True)
        adata.obs['sample'] = var_name
        Analysis_Room(adata,root,radio_string,radio_report,percent_mt=5, min_genes=150,min_cells=3,outpute=output)  # je lance la vraie analyse (pour l'instant je code l'analyse rapide mais aprés faudra passer des paramétres pour l'analyse détaillé)

    elif radio_string.get()=="deux" and okey_check:  # Multi Sample
        for key in dico_path.keys():
            dico_adata[key] = sc.read_10x_mtx(dico_path[key], var_names="gene_symbols", cache=True)
            dico_adata[key].obs['sample'] = key  # name the samples
            list_concat.append(dico_adata[key])

        adata = list_concat[0].concatenate(list_concat[1:])  # avec le 1er adata je concatene tous les autres adatas

        Analysis_Room(adata,root,radio_string,radio_report,percent_mt=5, min_genes=150, min_cells=3,outpute=output)  # je lance la vraie analyse , avec un adata unifié !!


def Analysis_Room(adata,root,radio_string,radio_report,outpute,percent_mt=5, min_genes=150, min_cells=3):
    ###show_image(okey_image,tab1,3,2)
    adata.var  # rowdata equivalent
    extract_qc_report(adata,root,radio_report,adresse=f"{place}/{outpute.get().strip()}",percent_mt=percent_mt,min_cells=min_cells,min_genes=min_genes)
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
    log_window(root,"Doublet Detection report", radio_report,s,x=450,y=50)
    blot=scrub.plot_histogram()[0]
    plt.savefig(f"{place}/{outpute.get().strip()}/Doublet_Detection.png")
    plot_window(root,"Doublet Detection plot",radio_report,blot)
    sum(adata.obs['predicted_doublets'])
    # add in column with singlet/doublet instead of True/False
    adata.obs['doublet_info'] = adata.obs["predicted_doublets"].astype(str)
    adata
    adata_clone = adata  # I clone it in order to use it for Combat batch correction later
    ####imagebox.grid_forget() #pour effacer l'image mais ça s'affiche pas

    #################LET'S PASS TO HGVs ###########################################################

    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5, batch_key="sample")
    var_genes = adata.var.highly_variable  # get the HVGs
    sc.pl.highly_variable_genes(adata,save=".png",show=False) #je vais tester les 2 plots
    png_window(root,radio_report,"Report HGVs","filter_genes_dispersion.png",w=850,h=500,x=10,y=450)
    os.system(f"cp {os.getcwd()}/figures/filter_genes_dispersion.png {place}/{outpute.get().strip()}/") #mv the png to the output dir !!
    sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, svd_solver='arpack')
    #sc.pl.pca_loadings(adata, components=[1,2,3,4,5,6,7,8])
    sc.pl.pca_variance_ratio(adata, log=True,show=False,save=".png") #faut que je teste show = True si ça me fait sortir les autres.
    os.system(f"cp {os.getcwd()}/figures/pca_variance_ratio.png {place}/{outpute.get().strip()}/") #mv the png to the output dir !!
    png_window(root,radio_report,"Report Pca","pca_variance_ratio.png",w=400,h=400,x=900,y=450) #the figure is saved so i can call the function to display it !

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
        sc.pl.rank_genes_groups(adata, n_genes=20, sharey=False,show=False,save=".png") # je vais voir comment on gére ce plot !
        os.system(f"cp {os.getcwd()}/figures/rank_genes_groups_louvain.png {place}/{outpute.get().strip()}/") #mv the png to the output dir !!
        extract_marker_report(adata,root,radio_report,f"{place}/{outpute.get().strip()}/Top_100_genes_by_clusters.csv",f"{place}/{outpute.get().strip()}/rank_genes_groups_louvain.png",x=1350,y=450,adresse=f"{place}/{outpute.get().strip()}")
        if radio_report.get() == "non":
            print("report hided !")
        else:
            im = Image.open(os.getcwd() + "/figures/rank_genes_groups_louvain.png") #to display the png !
            im.show()
        plt.close("all")
        if radio_report.get() !="non":
            sc.pl.umap(adata, color='louvain')
            plt.show(block=False)
        adata.write_h5ad(f"{place}/{outpute.get().strip()}/end_adata.h5ad")
        os.system(f"mv {place}/{outpute.get().strip()} $HOME/Desktop") #move the output in a visible place
        messagebox.showinfo("Info","Quick Analysis Done !")

    else: #so it's multisample !
        # create a new object with lognormalized counts
        adata_combat = sc.AnnData(X=adata_clone.X, var=adata_clone.var, obs=adata_clone.obs)
        # first store the raw data
        adata_combat.raw = adata_combat
        # Run Combat
        sc.pp.combat(adata_combat, key='sample')
        sc.pp.highly_variable_genes(adata_combat)
        sc.pl.highly_variable_genes(adata_combat,save=".png",show=False)
        os.system(f"cp {os.getcwd()}/figures/filter_genes_dispersion.png {place}/{outpute.get().strip()}/")
        png_window(root,radio_report,"Report HGVs (after Integration)", "filter_genes_dispersion.png", w=850, h=500,x=10,y=500)
        sc.pp.pca(adata_combat, n_comps=30, use_highly_variable=True, svd_solver='arpack')
        sc.pp.neighbors(adata_combat, n_pcs=30)
        sc.tl.umap(adata_combat)
        sc.tl.tsne(adata_combat, n_pcs=30)
        sc.tl.louvain(adata_combat, key_added="louvain",resolution=1.0)
        #tmp = pd.crosstab(adata_combat.obs['leiden_1.0'], adata_combat.obs['sample'], normalize='index')
        #tmp.plot.bar(stacked=True).legend(loc='upper right')
        sc.tl.rank_genes_groups(adata_combat, 'louvain', method='t-test')
        sc.pl.rank_genes_groups(adata_combat, n_genes=20, sharey=False,save=".png",show=False)
        os.system(f"cp {os.getcwd()}/figures/rank_genes_groups_louvain.png {place}/{outpute.get().strip()}")
        extract_marker_report(adata_combat,root,radio_report,f"{place}/{outpute.get().strip()}/Top_100_genes_by_clusters.csv",f"{place}/{outpute.get().strip()}/rank_genes_groups_louvain.png",x=1350,y=450,adresse=f"{place}/{outpute.get().strip()}")
        if radio_report.get() == "non":
            print("report hided !")
        else:
            im = Image.open(os.getcwd() + "/figures/rank_genes_groups_louvain.png") #to display the png !
            im.show()
        plt.close("all")
        if radio_report.get() !="non":
            sc.pl.umap(adata_combat,color='louvain')
            plt.show(block=False)
        adata_combat.write_h5ad(f"{place}/{outpute.get().strip()}/end_adata.h5ad")
        os.system(f"mv {place}/{outpute.get().strip()} $HOME/Desktop")
        messagebox.showinfo("Info","Quick Analysis Done !")

def extract_marker_report(bdata,root,radio_report,filename,png_name,x,y,adresse):
    nb = len(bdata.obs.louvain.unique())  # nb of clusters
    txt = "--------------------------"
    report = "{:s}\n The Marker gene detection plot is available : \n {:s} \n Number of clusters detected: {:d} \n Top 100 for each marker File: {:s} \n Comparaison method : t-test \n Correction method : Benjamini-Hocheberg \n p-value Cut-Off : 0.05 \n Clustering method : Louvain (Graph-based) \n{:s}".format(
        txt, png_name, nb, filename, txt)

    with open(f'{adresse}/marker_report', 'w') as f:
        f.write(report)

    log_window(root, "Marker detection report", radio_report, report, x, y)

    lista = bdata.uns['rank_genes_groups']['scores'].tolist()
    listb = bdata.uns['rank_genes_groups']['pvals'].tolist()
    listc = bdata.uns['rank_genes_groups']['pvals_adj'].tolist()
    listd = bdata.uns['rank_genes_groups']['logfoldchanges'].tolist()
    liste = bdata.uns['rank_genes_groups']['names'].tolist()

    listo = []

    for w, x, y, z, p in zip(lista, listb, listc, listd, liste):
        listing = []
        for v in range(0, len(x)):
            listing.append(w[v])
            listing.append(x[v])
            listing.append(y[v])
            listing.append(z[v])
            listing.append(p[v])
            listing.append(" ")  # to put spaces between each cluster
        listo.append(listing)

    leicester = []
    for i in range(0, len(bdata.obs['louvain'].cat.categories)):
        lister = ['scores', 'p-value', 'adjusted_p-value', 'LogFC', 'marker_cluster_', ' ']
        lister[4] = lister[4] + str(i)
        leicester.append(lister)

    pd.DataFrame(listo, columns=list(chain(*leicester))).head(100).to_csv(filename, index=False)  # get the first 100 top genes



def extract_qc_report(adata, root,radio_report,percent_mt, min_cells, min_genes,adresse):
    bdata = sc.AnnData(X=adata.X, obs=adata.obs, var=adata.var)
    begin_cells = len(bdata.obs);begin_genes = len(bdata.var)
    cdata = sc.AnnData(X=adata.X, obs=adata.obs,var=adata.var)  # i could've worked with adata rather than creating a new var but i'm afraid that it modify the real adata variable
    sc.pp.filter_cells(bdata, min_genes=min_genes);diff_cells = begin_cells - len(bdata.obs);sc.pp.filter_genes(bdata, min_cells=min_cells);diff_genes = begin_genes - len(bdata.var)
    cdata.var['mt'] = cdata.var_names.str.startswith('MT-')
    cdata.var[cdata.var.mt == True]
    sc.pp.calculate_qc_metrics(cdata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    cdata = cdata[cdata.obs.pct_counts_mt < percent_mt];diff_cells_mt = begin_cells - len(cdata.obs)
    sc.pp.filter_cells(cdata, min_genes=min_genes);sc.pp.filter_genes(cdata, min_cells=min_cells);total_cells = len(cdata.obs);total_genes = len(cdata.var)
    txt = "--------------------------"
    report = "{0}\n Number of cells that express less than {1} genes: {2} \n Number of genes that are expressed by less than {3} cells : {4} \n Number of cells that express more than {5}% of mitonchondrial genes : {6} \n Total number of cells : {7} \n Total number of genes : {8} \n {9}"\
        .format(txt, min_genes, diff_cells, min_cells, diff_genes, percent_mt, diff_cells_mt, total_cells, total_genes, txt)
    with open(f'{adresse}/qc_report', 'w') as f:
        f.write(report)
    log_window(root,"QC report",radio_report,report)

#</editor-fold>
################################################################################Function tab2 (i can also use the functions that are above, for tab2)

def progress_log(root,listo,progress,num,heyt,step): #at which step of progress list do we have to stop
    window = Toplevel(root)
    window.title("Warning")
    window.geometry('420x%d+%d+%d' % (heyt,100, 250))
    #print(listo,progress)
    Label(window, text="Progression of your Analysis", font=Font(root, size=15, weight=BOLD)).pack()
    for i in range(0, len(progress[0:num])):
        Label(window, text=progress[i], font=Font(root, size=11, weight=BOLD)).pack()
        if listo[i] == False:
            Label(window, text="Not Done", font=Font(root, size=11, weight=BOLD), foreground="red").pack()
        else:
            Label(window, text='Done or Skipped', font=Font(root, size=11, weight=BOLD), foreground="green").pack()

    Label(window, text="Please complete the undone steps before \n{0}".format(step),
          font=Font(root, size=13, weight=BOLD)).pack()
    sola = Label(window)
    solb = Label(window)
    sola.image = ImageTk.PhotoImage(file=f"{placa}warn.png")
    solb.image = ImageTk.PhotoImage(file=f"{placa}warn.png")
    sola['image'] = sola.image
    solb['image'] = solb.image
    sola.place(x=0, y=0)
    solb.place(x=390, y=0)

def clust_report(tmp):
    dico = {}
    for i in range(0, len(tmp)):
        dico[tmp.iloc[i].name] = {}
        for j in range(0, len(tmp.columns.categories)):
            dico[tmp.index.categories[i]][j] = tmp.iloc[i, j]

    gos = ""
    for i in tmp.columns.categories:
        gos = gos + "\nCluster {0} = {1} cells ".format(i, tmp.iloc[:, int(i)].sum())

    gas = ""
    for i in dico.keys():
        gas = gas + "\n\nReport for sample {0} :\n".format(i)
        for j in dico[i]:
            gas = gas + "\nCluster {0} = {1} cells ".format(j, dico[i][j])

    return (gos,gas)