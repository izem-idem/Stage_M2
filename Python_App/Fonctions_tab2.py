import tkinter as tk
from tkinter import *
from tkinter import  ttk
from tkinter import messagebox
import scrublet as scr
import anndata
from PIL import Image, ImageTk
import scanorama
import io
from contextlib import redirect_stdout
import os,time,scanpy as sc, anndata as ann
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import numpy as np
import random
import pandas as pd
import Fonctions_tab1 as Ftab1
import matplotlib.pyplot as plt
import seaborn as sb
import math
from tkinter.font import BOLD, Font
from tkinter import ttk, simpledialog
from itertools import chain

placa="/home/izem/Pictures/"

place=os.getcwd()
adata_englobant="will englobe all the changes that will be done to adata !"
okey_check=False #if false then okey button is not clicked and i can't allow the user to go and click on the other buttons so for the 7 buttons I have to check if okey button is cliqued before !!

global_dico_adata={} #will contain the many samples before i will concatenate
global_list_concat=[] #i will this list to concatenate the datas !!

global_dico_path={}

qc_done=False #it's my flag for qc being done
norm_done=False
int_done=False
feat_done=False
dim_clust_done=False
mark_done=False


global_format=[]

flag_single=False
#### at each process i save my adata before the next process
txt="---------------"
mnn_flag=False
progress=['QC Filtering','Normalization','Integration','Feature Selection','Dim-Reduc Clustering','Marker Detection']
#listo=[qc_done,norm_done,int_done,feat_done,dim_done,clust_done,mark_done]

def check_path(radio_string, samplename_count_tab, fullpath_count_tab,formate,output):
    if radio_string.get() == "" or formate.get()=="":
        messagebox.showerror("Input Error", "Error: Please choose one option for the 2 propositions")

    if output.get()=="":
        messagebox.showerror("Input Error","Error:\nPlease provide a dir_name to store the output")

    if samplename_count_tab.get() == "" or fullpath_count_tab.get() == "":
        messagebox.showerror("Input Error", "Error: Please You have to fill the two forms before going on ! ")

    # I have to precise to the user to not use blank space in sample name
    if radio_string.get() == "un" and not samplename_count_tab.get() == "" and not fullpath_count_tab.get() == "" and output.get().strip()!="":
        global global_varpath, global_varname, adata
        global_varpath = fullpath_count_tab.get().strip(" ")
        global_varname = samplename_count_tab.get()  # no need to do .strip because it's free from space.

        if formate.get()=="10x_mtx":
            try:
                global okey_check
                adata = sc.read_10x_mtx(global_varpath, var_names="gene_symbols", cache=True)
                messagebox.showinfo("Information Count Matrix", "Your sample was successfully loaded !")
                okey_check=True
            except :
                messagebox.showerror("File Error", "Error: the file you provided does not fit")
        elif formate.get()=="mtx":
            try:
                adata = sc.read_mtx(global_varpath)
                messagebox.showinfo("Information Count Matrix", "Your sample was successfully loaded !")
                okey_check = True
            except:
                messagebox.showerror("File Error", "Error: the file you provided does not fit")
        elif formate.get()=="h5ad":
            try:
                adata = sc.read_h5ad(global_varpath)
                messagebox.showinfo("Information Count Matrix", "Your sample was successfully loaded !")
                okey_check = True
            except:
                messagebox.showerror("File Error", "Error: the file you provided does not fit!")
        elif formate.get()=="loom":
            try:
                adata = sc.read_loom(global_varpath)
                messagebox.showinfo("Information Count Matrix", "Your sample was successfully loaded !")
                okey_check = True
            except:
                messagebox.showerror("File Error", "Error: the file you provided does not fit!")


    if radio_string.get() == "deux" and not samplename_count_tab.get() == "" and not fullpath_count_tab.get() == "" and output.get()!="":
        global global_dico_path, global_list_concat,global_format
        global_varpath = fullpath_count_tab.get().strip(" ")
        global_varname = samplename_count_tab.get()  # no need to do .strip because it's free from space.
        flag=0
        print(global_varpath)
        if formate.get() == "10x_mtx":
            try:
                sc.read_10x_mtx(global_varpath, var_names="gene_symbols", cache=True)
                global_format.append(formate.get())
                flag=1
            except:
                messagebox.showerror("File Error", "Error: the file you provided does not fit!")
        elif formate.get() == "mtx":
            try:
                adata = sc.read_mtx(global_varpath)
                global_format.append(formate.get())
                flag = 1
            except:
                messagebox.showerror("File Error", "Error: the file you provided does not fit")
        elif formate.get() == "h5ad":
            try:
                adata = sc.read_h5ad(global_varpath)
                global_format.append(formate.get())
                flag = 1
            except:
                messagebox.showerror("File Error", "Error: the file you provided does not fit!")

        elif formate.get() == "loom":
            try:
                adata = sc.read_loom(global_varpath)
                global_format.append(formate.get())
                flag = 1
            except:
                messagebox.showerror("File Error", "Error: the file you provided does not fit!")

        if global_varname not in global_dico_path.keys() and global_varpath not in global_dico_path.values() and flag==1:  # to avoid having the same sample name or the same sample (ie avoiding same fullpath)
            global_dico_path[global_varname] = global_varpath
            messagebox.showinfo("Information Count Matrix","Path entered correctly !!\nPlease add another path to add another sample\n\nYou have currently " + str(len(global_dico_path.keys())) + " sample(s) that are registered")
            okey_check = True
            global qc_done,norm_done,int_done,feat_done,dim_done,clust_done,mark_done
            qc_done=False; norm_done=False; int_done=False; feat_done=False; dim_done=False; clust_done=False; mark_done=False# if the use clik on okey , we restart all the analysis

        elif global_varname in global_dico_path.keys() or global_varpath in global_dico_path.values() and flag==1:  # print error of redundant values
            messagebox.showerror("Input Error","Error: the full path or the sample name you entered are already used !!")

# <editor-fold desc="QC_STEP">
def QC_FILTER(root,radio_string,format,output):
    if okey_check==False:
        messagebox.showerror("Input Error","Please check your data by clicking on 'OKEY' before going on !")
    else:
        wind = Toplevel(root)
        wind.title("QC FILTERING")
        wind.geometry('+%d+%d' % (100, 50))

        percent_mt = tk.DoubleVar()
        min_cells = tk.IntVar()
        min_genes = tk.IntVar()
        yes_no_report = tk.StringVar()
        yes_no_doublet = tk.StringVar()


        Entry(wind, textvariable=percent_mt, width=7).grid(column=1, row=0, sticky="W")
        Label(wind, text="Maximum % of MT transcripts :").grid(column=0, row=0, sticky="W")
        Entry(wind, textvariable=min_cells, width=7).grid(column=1, row=1, sticky="W")
        Label(wind, text="Min number of genes expressed in a cell :").grid(column=0, row=1, sticky="W")
        Entry(wind, textvariable=min_genes, width=7).grid(column=1, row=2, sticky="W")
        Label(wind, text="Minimum number of cells expressing a gene :").grid(column=0, row=2, sticky="W")
        percent_mt.set(10);min_cells.set(3); min_genes.set(150)
        Label(wind, text="Doublet Detection :").grid(column=0, row=3, sticky="W")
        Radiobutton(wind, text="NO", variable=yes_no_doublet, value="No").place(x=130, y=69)
        Radiobutton(wind, text="YES", variable=yes_no_doublet, value="Yes").place(x=190, y=69)
        Label(wind, text="Show reports :").grid(column=0, row=4, sticky="W")
        Radiobutton(wind, text="NO", variable=yes_no_report, value="non").place(x=130, y=90)
        Radiobutton(wind, text="YES", variable=yes_no_report, value="oui").place(x=190, y=90)
        adata=Pre_Process(radio_string,format)
        before_qc=adata # i save my original adata before the qc process
        print("Begin QC_filter , I will work with this adata:",before_qc)
        os.system(f"mkdir {output.get().strip()}") #the dir, create it
        before_qc.write_h5ad(f"{place}/{output.get().strip()}/before_qc.h5ad")
        Button(wind, text="Visualize Plots",command=lambda : Visualize_plots(before_qc,root,"QC Plots",output)).grid(column=0, row=5, sticky="W", pady=5)
        Button(wind, text="Launch QC",command=lambda :QC_Analysis(percent_mt.get(),min_cells.get(),min_genes.get(),yes_no_doublet,yes_no_report,before_qc,root,output)).grid(column=1,row=5,sticky="E")
        Button(wind, text="Skip QC",command=lambda :Skip_QC(before_qc,output)).place(x=160,y=118)



def Pre_Process(radio_string,formate):
    if radio_string.get() == "un":  # One Sample
        if formate.get()=="10x_mtx":
            try:
                adata_englobant = sc.read_10x_mtx(global_varpath, var_names="gene_symbols", cache=True)
            except :
                messagebox.showerror("File Error", "Error: the file you provided does not fit")
        elif formate.get()=="mtx":
            try:
                adata_englobant = sc.read_mtx(global_varpath)
            except:
                messagebox.showerror("File Error", "Error: the file you provided does not fit")
        elif formate.get()=="h5ad":
            try:
                adata_englobant = sc.read_h5ad(global_varpath)
            except:
                messagebox.showerror("File Error", "Error: the file you provided does not fit!")
        elif formate.get()=="loom":
            try:
                adata_englobant = sc.read_loom(global_varpath)
            except:
                messagebox.showerror("File Error", "Error: the file you provided does not fit!")

        adata_englobant.obs['sample'] = global_varname
        return adata_englobant

    else:  # Multi Sample
        print(global_dico_path,global_list_concat)
        print(global_format)
        tint=0
        for key , format in zip(global_dico_path.keys(),global_format):
            if format == "10x_mtx":
                try:
                    print("le if de 10x_mtx")
                    global_dico_adata[key] = sc.read_10x_mtx(global_dico_path[key], var_names="gene_symbols", cache=True)
                    global_dico_adata[key].obs['sample'] = key  # name the samples
                    global_list_concat.append(global_dico_adata[key])
                    print(global_dico_adata[key])
                except:
                    messagebox.showerror("File Error", "Error: the file you provided does not fit!")
            elif format == "mtx":
                try:
                    global_dico_adata[key] = sc.read_mtx(global_dico_path[key])
                    global_dico_adata[key].obs['sample'] = key  # name the samples
                    global_list_concat.append(global_dico_adata[key])
                except:
                    messagebox.showerror("File Error", "Error: the file you provided does not fit!")
            elif format == "h5ad":
                try:
                    global_dico_adata[key] = sc.read_h5ad(global_dico_path[key])
                    global_dico_adata[key].obs['sample'] = key  # name the samples
                    global_list_concat.append(global_dico_adata[key])

                except:
                    messagebox.showerror("File Error", "Error: the file you provided does not fit!")
            elif format == "loom":
                try:
                    print("le if du loom")
                    global_dico_adata[key] = sc.read_loom(global_dico_path[key])
                    global_dico_adata[key].obs['sample'] = key  # name the samples
                    global_list_concat.append(global_dico_adata[key])
                    print(global_dico_adata[key])
                except:
                    messagebox.showerror("File Error", "Error: the file you provided does not fit!")

        adata_englobant = global_list_concat[0].concatenate(global_list_concat[1:])  # avec le 1er adata je concatene
        return adata_englobant


def Visualize_plots(before_qc_adata,root,title,output):
    clone_adata=None
    if qc_done==False: #there is no previous qc so i will use the "raw" adata
        before_qc_adata.var['mt'] = before_qc_adata.var_names.str.startswith('MT-')
        before_qc_adata.var[before_qc_adata.var.mt == True]
        sc.pp.calculate_qc_metrics(before_qc_adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
        clone_adata=before_qc_adata
    else: #i want to display the qcied adata so i use the before_norm
        clone_adata=anndata.read_h5ad("{0}/{1}/before_norm.h5ad".format(place,output.get().strip()))
    fig, ax = plt.subplots(2, 2, figsize=(10, 5))
    fig.tight_layout(w_pad=3, h_pad=3)
    plt.subplots_adjust(bottom=0.1, left=0.1)
    ax1 = ax[0,0];ax2 = ax[0,1];ax3 = ax[1,0]
    ax1.set_xlabel("Number of genes expressed")
    ax1.set_ylabel("Number of cells")
    ax2.set_xlim([0,30])
    ax2.set_xlabel("% of UMIs mapped to mitochondrial transcripts")
    ax2.set_ylabel("Number of cells")
    ax3.set_xlabel("Number of cells")
    ax3.set_ylabel("Number of genes")
    sb.histplot(clone_adata.obs['n_genes_by_counts'], kde=False, bins=100, ax=ax[0, 0])
    sb.histplot(clone_adata.obs['pct_counts_mt'], kde=False, bins=100, ax=ax[0, 1])
    sb.histplot(clone_adata.var['n_cells_by_counts'], kde=False, bins=100, ax=ax[1, 0])
    sb.histplot(clone_adata.obs['total_counts'], kde=False, bins=100, ax=ax[1, 1])
    Ftab1.plot_window(root, title, "oui" , plt.gcf(), single=None)

def QC_Analysis(percent_mt,min_cells,min_genes,doublet,report,adata,root,outpute):
    if doublet.get() == "" or report.get() == "":
        messagebox.showerror("Input Error","Either Show reports or Doublet Detection is missing.\n Please complete the 2 options before carrying on")
    else:
        global qc_done
        global before_norm #i use it at the end of the code
        qc_done=True #the qc_analysis is done , i want to save this information in a variable
        Ftab1.extract_qc_report(adata,root,report,percent_mt,min_cells,min_genes,adresse=f"{place}/{outpute.get().strip()}")
        adata.var['mt'] = adata.var_names.str.startswith('MT-')
        adata.var[adata.var.mt == True]
        sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
        # ngenes_by_count is unclear and is different from n_genes + no answer from scanpy regarding this issue
        adata = adata[adata.obs.pct_counts_mt < percent_mt]
        sc.pp.filter_cells(adata, min_genes=min_genes)
        sc.pp.filter_genes(adata, min_cells=min_cells)
        ##############some information about cell cycle genes
        cell_cycle_genes = [x.strip() for x in open('cell_cycle_genes.txt')]  # I have to have this file in my dir
        # print(len(cell_cycle_genes))
        # Split into 2 lists
        s_genes = cell_cycle_genes[:43]
        g2m_genes = cell_cycle_genes[43:]
        cell_cycle_genes = [x for x in cell_cycle_genes if x in adata.var_names]
        # print(len(cell_cycle_genes))
        s_genes = [x for x in s_genes if x in adata.var_names] #ensure to get only the genes that i have in my dataset
        g2m_genes = [x for x in g2m_genes if x in adata.var_names]
        if len(s_genes) and len(g2m_genes) != 0:
            sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)
        else:
            print("WARNING: all the cell_cycle are abscent from your dataset")
        if doublet.get() =="Yes":
            try:
                f = io.StringIO()
                with redirect_stdout(f):
                    scrub = scr.Scrublet(adata.X)
                    adata.obs['doublet_scores'], adata.obs['predicted_doublets'] = scrub.scrub_doublets()  # c'est cette ligne qui produit l'output
                    s = f.getvalue()
                Ftab1.log_window(root, "Doublet Detection report", report, s, x=450, y=50)
                blot = scrub.plot_histogram()[0]
                plt.savefig(f"{place}/{outpute.get().strip()}/Doublet_Detection.png")
                Ftab1.plot_window(root, "Doublet Detection plot", report, blot, single=None)
                sum(adata.obs['predicted_doublets'])
                # add in column with singlet/doublet instead of True/False
                adata.obs['doublet_info'] = adata.obs["predicted_doublets"].astype(str)
            except ZeroDivisionError: #TO AVOID THE ZERODIV ERROR
                messagebox.showwarning("Warning","Couldn't do the doublet detection\nbecause some element are divided by zero")

            #before_norm = sc.AnnData(X=adata.X,var=adata.var,obs=adata.obs)
            adata.write_h5ad("{0}/{1}/before_norm.h5ad".format(place,outpute.get().strip())) #the real before_norm is now craved in stone
            messagebox.showinfo("QC Info", "The QC Filtering process is finished !")
        else: #don't do doublet detection
            adata.write_h5ad("{0}/{1}/before_norm.h5ad".format(place,outpute.get().strip()))  # the real before_norm is now craved in stone
            messagebox.showinfo("QC Info", "The QC Filtering process is finished !")


def Skip_QC(adata,output):
    print("QC skipped,\nI print the before_Norm object:", adata)
    global qc_done #i have to redeclare the qc_done as being a global variable
    qc_done=True #i consider the qc done (so that the user can continue the analysis)
    #global before_norm
    #before_norm=adata #this adata is the before_qc so i pass it to befor_norm
    adata.write_h5ad("{0}/{1}/before_norm.h5ad".format(place,output.get().strip()))  # the real before_norm is now craved in stone
    messagebox.showwarning("QC Info","You have skipped the\nQC step")
    print("QC_done",qc_done)

# </editor-fold>


# <editor-fold desc="NORMALIZATION STEP">
def Normalizer(root,radio_string,output):
    if okey_check==False:
        messagebox.showerror("Input Error","Please check your data by clicking on 'OKEY' before going on !")
    elif qc_done == False:
        Ftab1.progress_log(root, [qc_done], progress, num=1,heyt=200,step="Normalization")
    else:
        qbel_norm=anndata.read_h5ad("{0}/{1}/before_norm.h5ad".format(place,output.get().strip()))
        print("Begin normalization: \n I will work with this adata : \n",qbel_norm)  #check if it's the child of qc_analysis
        normal_report = tk.StringVar()
        oui_non_report = tk.StringVar()
        taga = Toplevel(root)
        taga.title("NORMALIZATION")
        taga.geometry('+%d+%d' % (100, 50))

        Label(taga, text="Types of Normalization:", font=Font(root, size=11, weight=BOLD)).grid(column=1, row=0,sticky="W")
        Radiobutton(taga, text="CPM Normalization", value="cpm", variable=normal_report).grid(column=0, row=1)
        Radiobutton(taga, text="Normalization by Deconvolution", value="deconv", variable=normal_report).grid(column=1, row=1)
        Label(taga, text="Show reports :").grid(column=0, row=5, sticky="W")
        Radiobutton(taga, text="NO", variable=oui_non_report, value="non").place(x=90, y=50)
        Radiobutton(taga, text="YES", variable=oui_non_report, value="Yes").place(x=150, y=50)
        Button(taga, text="Normalize",command=lambda : Normal_Analysis(root,normal_report,oui_non_report,qbel_norm,radio_string,output)).grid(column=2, row=5, sticky="E")
        Button(taga, text="Skip Normalization",command=lambda : Skip_Norm(root,qbel_norm,radio_string,output)).grid(column=1, row=5,sticky="E",padx=20)

def Normal_Analysis(root,normal_report,oui_non_report,adata,radio_string,outpute): #must provide the before_norm
    global norm_done,int_done,flag_single #so that i can use it later
    if normal_report.get() == "" or oui_non_report.get() == "":
        messagebox.showerror("Input Error","Please choose a Normalization option and a Report option to carry on")
    elif normal_report.get() == "cpm":
        sc.pp.normalize_total(adata, target_sum=1e6)  # normalize every cell with 10^6 UMIs
        sc.pp.log1p(adata)  # to stabilize the variance let's transform it with log2
        lista = [adata.X[i, :].sum() for i in range(0, adata.X.shape[0])] #calculate the sum (ie col sum of matrix) of each cell (ie each col)
        fig, ax = plt.subplots(1)
        fig.tight_layout(w_pad=5, h_pad=5)
        plt.subplots_adjust(left=0.1, bottom=0.1)
        ax.set_xlabel('Library Size')
        ax.set_title("Distribution of Total UMIs by cell")
        ax.set_ylabel('Number of cells')
        sb.histplot(lista)
        Ftab1.plot_window(root,"Cpm Normalization",oui_non_report,plt.gcf(),single=None)
        plt.savefig(f"{place}/{outpute.get().strip()}/Cpm_normalization.png") #save as png
        #global before_integr
        adata.raw=adata #by raw i mean saving the logarithmized X in order to perform marker detection
        #before_integr=adata #save the last adata in a global variable
        print(adata.uns['log1p'])
        adata.write_h5ad("{0}/{1}/before_integr.h5ad".format(place,outpute.get().strip()))
        norm_done = True  # the flag is up , normalization has been done
        if radio_string.get() == "un":
            int_done = True
            flag_single=True
        messagebox.showinfo("Normalization Info", "The Normalization process is finished !")

    elif normal_report.get()=="deconv":

        adata.write_h5ad(f"{os.getcwd()}/test.h5ad") #save the current adata in h5ad
        filou = "{0}/deconv.R".format(os.getcwd()) #this path should also be automatic because all the machines do not have /home/izem/ dir !
        #file = "/home/izem/test.h5ad"
        os.system("Rscript --vanilla {:s} {:s} {:s}".format(filou,os.getcwd(),os.getcwd()+"/test.h5ad"))
        adata = sc.read_h5ad(os.getcwd()+"/complete.h5ad") #i get the adata which is normalized by deconv in R !
        #let's plot it

        fig, ax = plt.subplots(1, 3, figsize=(20, 5))
        fig.tight_layout(w_pad=5, h_pad=5)
        plt.subplots_adjust(left=0.05, bottom=0.1)
        #sc.pl.scatter(adata, 'sizeFactor', 'total_counts',color="sample",ax=ax[0],show=False,legend_loc='none')
        ax[0].scatter(adata.obs['libSizeFactor'], adata.obs['sizeFactor'], color="g", s=10)
        ax[0].set_xlabel("Library Size Factor")
        ax[0].set_title("Deconvolved Size Factor vs Library Size Factor")
        ax[0].set_ylabel("Deconvoluted Size Factor")
        ax[1].scatter(adata.obs['sizeFactor'],adata.obs['n_genes'], color="skyblue",s=10)
        ax[1].set_xlabel("Deconvoluted Size Factor")
        ax[1].set_ylabel("Number of genes")
        sb.histplot(adata.obs['sizeFactor'], bins=50, kde=False, ax=ax[2],color='salmon')
        ax[2].set_ylabel("Number of cells")
        ax[2].set_title("Distribution of Deconvolved Size Factor")
        Ftab1.plot_window(root, "Normalization by Deconvolution", oui_non_report, plt.gcf(), single=None)
        plt.savefig(f"{place}/{outpute.get().strip()}/Deconv_normalization.png") #save as png
        adata.raw = adata  # save the raw adata(log data)
        #before_integr = adata #save the last adata in a global variable to be ready for the next step (integration)
        adata.uns['log1p']={}
        adata.uns['log1p']['base']=None
        adata.write_h5ad("{0}/{1}/before_integr.h5ad".format(place,outpute.get().strip()))

        norm_done = True  # the flag is up , normalization has been done
        if radio_string.get() == "un":
            int_done = True
            flag_single=True
        messagebox.showinfo("Normalization Info", "The Normalization process is finished !")

def Skip_Norm(adata,radio_string,outpute):
    print("Norm skipped,\nI print the before_integer object:", adata)
    global norm_done,int_done,flag_single #i have to redeclare the qc_done as being a global variable
    norm_done=True #i consider the qc done (so that the user can continue the analysis)
    if radio_string.get() =="un" :
        int_done=True #also validate this step so that we can help each other
        flag_single=True
    global before_integr
    adata.raw=adata #save the raw for other usage
    #before_integr=adata #this adata is the before_qc so i pass it to befor_norm
    adata.write_h5ad("{0}/{1}/before_integr.h5ad".format(place,outpute.get().strip()))
    messagebox.showwarning("Normalization Info","You have skipped the\nNormalization step")

# </editor-fold>


# <editor-fold desc="INTEGRATION">

def Integrater(root,radio_string,output):
    if okey_check==False:
        messagebox.showerror("Input Error","Please check your data by clicking on 'OKEY' before going on !")
    else:
        if radio_string.get()== "un":
            messagebox.showwarning("Warning Process","You have only 1 sample \nYou don't need to Integrate.")
            global int_done
            int_done=True #validate this step for the unique sample
        if radio_string.get()=="deux":
            if qc_done==False  or norm_done==False :
                Ftab1.progress_log(root, [qc_done,norm_done], progress, num=2,heyt=200,step="Integration")
            elif qc_done and norm_done:
                qbel_integr = sc.read_h5ad(f"{place}/{output.get().strip()}/before_integr.h5ad")
                print("Begin Integration , I will work with this adata (before_interger) : ", qbel_integr)
                taga = Toplevel(root)
                taga.title("INTEGRATION")
                taga.geometry('+%d+%d' % (100, 300))

                integre_report = tk.StringVar()

                Label(taga, text="Integration Mode :", font=Font(taga, size=11, weight=BOLD)).grid(column=1, row=0,sticky="W")
                Radiobutton(taga, text="MNN", value="mnn", variable=integre_report).grid(column=0, row=1, sticky="W")
                Radiobutton(taga, text="Scanorama", value="scan", variable=integre_report).grid(column=1, row=1)
                Radiobutton(taga, text="Combat", value="combat", variable=integre_report).grid(column=2, row=1)
                Button(taga, text="Skip Integration",command=lambda :Skip_integration(qbel_integr,output)).grid(column=0, row=5, sticky="W")
                Button(taga, text="Integrate",command=lambda :Integration_Analysis(integre_report,qbel_integr,output)).grid(column=2, row=5, sticky="E")

def Skip_integration(adata,outpute):
    global int_done
    print("integration skipped,\nI print the before_feature object:",adata)
    messagebox.showwarning("Integration Info","You have skipped the\nIntegration step :( ")
    int_done=True
    global before_feature
    #before_feature=adata #here adata is actually before_integer , we pass directly to before feature.
    adata.write_h5ad("{0}/{1}/before_feature.h5ad".format(place,outpute.get().strip()))

def Integration_Analysis(integre_report,adata,outpute):
    bdata=sc.AnnData(X=adata.X,obs=adata.obs,var=adata.var,uns=adata.uns) #each function call , i want to refresh bdata and cdata
    cdata=sc.AnnData(X=adata.X,obs=adata.obs,var=adata.var,uns=adata.uns)
    global int_done
    int_done=True # that's if for the integration I raise the flag.
    if integre_report.get() == "mnn":
        global mnn_flag
        mnn_flag=True
        adata.write_h5ad(os.getcwd() + '/mnn.h5ad') #write the current adata which is befor_integer
        adata.obs['sample'] = adata.obs['sample'].astype('category')  # make it categorical
        batches = adata.obs['sample'].cat.categories.tolist()
        # print(batches)

        alldata = {}
        for batch in batches:
            alldata[batch] = adata[adata.obs['sample'] == batch,]
        listing = list(alldata.keys())
        lista = []
        for x in range(0, len(batches)): lista.append("alldata[listing[{:d}]],".format(x))
        phrase = "".join(lista)
        phrase = phrase[:len(phrase) - 1]  # on enléve la petite virgule

        #print(phrase)
        os.system("sed -i 's/xxx/{:s}/g' {:s}".format(phrase, os.getcwd() + "/MNN.py")) #will put the phrase in the function !!
        command = 'python3 '+os.getcwd() +'/MNN.py '+phrase  #problem of dir also for this file
        r = os.popen(command)  # Execute command
        info = r.readlines()  # read command output
        listo = [info[i].strip('\r\n') for i in range(0, len(info))]  # pour chaque element de lista , enlever les \n et \t
        if listo.__contains__("Done."):  # if Mnn is done correctly then !!
            ddata=anndata.read_h5ad(os.getcwd()+"/mnn.h5ad") # I recuperate my adata with the correction value on it.
            global before_feature
            #bdata.raw = adata.raw # I conserve the raw matrix (in adata.raw.X) (uncorrected) to perform marker detection
            #before_feature=bdata
            ddata.write_h5ad("{0}/{1}/before_feature.h5ad".format(place,outpute.get().strip()))
            messagebox.showinfo("Integration Info","The MNN correction\n was done successfully !!")
            #print("Anndata after MNN Integration:",before_feature)

            djoum = list(phrase)
            for i in range(0,len(djoum)):  # before each [ i will put a backslash to avoid the meta characters which is used in regex
                if djoum[i] == '[' or djoum[i] == ']':
                    djoum[i] = "\\" + djoum[i]

            phrase = "".join(djoum)

            os.system("sed -i 's/{:s}/xxx/g' {:s}".format(phrase, os.getcwd() + "/MNN.py")) #i want xxxx to be back in first arg of mnn_correct
        else:
            messagebox.showerror("Integration Problem","Something went wrong during the MNN integration")
            #maybe print some of log message
    elif integre_report.get() == "scan":
        bdata.obs['sample'] = bdata.obs['sample'].astype('category')  # make it categorical
        batches = bdata.obs['sample'].cat.categories.tolist()
        alldata = {}
        for batch in batches:
            alldata[batch] = bdata[bdata.obs['sample'] == batch,]
        # convert to list of AnnData objects
        adatas = list(alldata.values())
        # run scanorama.integrate
        scanorama.integrate_scanpy(adatas, dimred=50)
        #scanorame put it's corrected matrix in the .obsm of the adata (for each adata)
        scanorama_int = [ad.obsm['X_scanorama'] for ad in adatas] #put all the matrix for each sample in one list
        # make into one matrix.
        all_s = np.concatenate(scanorama_int)
        #print(all_s.shape)
        # add to the AnnData object
        bdata.obsm["Scanorama"] = all_s
        messagebox.showinfo("Integration Info", "The Scanorama correction\n was done successfully !!")
        #before_feature=bdata
        bdata.write_h5ad("{0}/{1}/before_feature.h5ad".format(place,outpute.get().strip()))
        #print("Anndata after Scanorama Integration:",before_feature)
    elif integre_report.get() == "combat":
        cdata.raw=cdata
        # Run Combat
        sc.pp.combat(cdata, key='sample')
        #before_feature=cdata
        cdata.write_h5ad("{0}/{1}/before_feature.h5ad".format(place,outpute.get().strip()))
        messagebox.showinfo("Integration Info", "The Combat correction\n was done successfully !!")
        #print("Anndata after Combat Integration:",before_feature)

# </editor-fold>


# <editor-fold desc="FEATURE SELECTION STEP">
def Selecter(root, radio_string,output):
    if okey_check==False:
        messagebox.showerror("Input Error","Please check your data by clicking on 'OKEY' before going on !")
    else:
        global qc_done
        if (qc_done==False  or norm_done==False or int_done==False) and radio_string.get()=="deux":
            Ftab1.progress_log(root,[qc_done,norm_done,int_done],progress,num=3,heyt=210,step="Feature Selection")
        elif (qc_done==False  or norm_done==False) and radio_string.get()=="un":
            Ftab1.progress_log(root,[qc_done,norm_done], progress, num=2, heyt=210,step="Feature Selection")
        else:
            if flag_single: #if this flag is true don't try to find before_feature.h5ad because before_integer is here
                qbel_integr = anndata.read_h5ad("{0}/{1}/before_integr.h5ad".format(place,output.get().strip()))
                qbel_integr.uns['log1p']['base'] = None
                print("Begin Feature selection , I will work with this adata ( before_integer) : ", qbel_integr)
            else:
                try:
                    qbel_feature=anndata.read_h5ad("{0}/{1}/before_feature.h5ad".format(place,output.get().strip()))
                    if 'log1p' in qbel_feature.uns.keys():
                        qbel_feature.uns['log1p']['base']=None
                    else: #if the h5ad comes from MNN it doesn't contain the log1p
                        qbel_feature.uns['log1p']={}
                        qbel_feature.uns['log1p']['base']=None
                    print("Begin Feature selection , I will work with this adata (before_feature) : ", qbel_feature)
                except FileNotFoundError: #ie we have only one sample so no before_feature adata has been generated. if try works it means that we have a multi-sample adata.
                    messagebox.showerror("error","Something went wrong \nthe file is not found")


            wind = Toplevel(root)
            wind.title("Feature Selection")
            wind.geometry('+%d+%d' % (100, 250))

            feature_report = tk.StringVar()
            ih_lela_report = tk.StringVar()

            seurat_a = tk.DoubleVar()
            seurat_b = tk.DoubleVar()
            seurat_c = tk.DoubleVar()

            seurat_a.set(0.0125)
            seurat_b.set(3)
            seurat_c.set(0.5)

            v3_a = tk.IntVar()

            pearson_a = tk.IntVar()
            pearson_b = tk.IntVar()
            pearson_c = tk.IntVar()

            pearson_a.set(100)
            pearson_b.set(3000)
            pearson_c.set(1000)

            Radiobutton(wind, text="Seurat Method", value="seurat", variable=feature_report).grid(column=0, row=0, columnspan=2, pady=5)
            Label(wind, text="Minimum Mean :").grid(column=0, row=1)
            Label(wind, text="Maximum Mean :").grid(column=0, row=2)
            Label(wind, text="Minimum Dispersions :").grid(column=0, row=3)
            Entry(wind, textvariable=seurat_a, width=9).grid(column=1, row=1)
            Entry(wind, textvariable=seurat_b, width=9).grid(column=1, row=2)
            Entry(wind, textvariable=seurat_c, width=9).grid(column=1, row=3)

            Radiobutton(wind, text="SeuratV3 Method", value="v3", variable=feature_report).grid(column=2, row=0)
            Label(wind, text="Top Number of genes \nto choose:").grid(column=2, row=1, rowspan=2)
            Entry(wind, textvariable=v3_a, width=9).grid(column=2, row=3)

            Radiobutton(wind, text="Pearson Residual Method", value="pearson", variable=feature_report).grid(column=3, row=0, columnspan=2)
            Label(wind, text="Theta value :").grid(column=3, row=1)
            Label(wind, text="Top Number of genes:").grid(column=3, row=2)
            Label(wind, text="Chunksize:").grid(column=3, row=3)
            Entry(wind, textvariable=pearson_a, width=9).grid(column=4, row=1)
            Entry(wind, textvariable=pearson_b, width=9).grid(column=4, row=2)
            Entry(wind, textvariable=pearson_c, width=9).grid(column=4, row=3)

            Label(wind, text="Show reports :").grid(column=0, row=4, sticky="W", pady=10)
            Radiobutton(wind, text="NO", variable=ih_lela_report, value="non").place(x=100, y=112)
            Radiobutton(wind, text="YES", variable=ih_lela_report, value="Yes").place(x=160, y=112)


            if radio_string.get()=="un": #one sample, i'll have to pass before_integer which the most recent adata.
                Button(wind, text="Perform Selection",command=lambda : Feature_Select_Analysis(feature_report,ih_lela_report,qbel_integr,seurat_a,seurat_b,seurat_c,pearson_a,pearson_b,pearson_c,v3_a,output)).grid(column=4, row=6, sticky="E")
                Button(wind, text="Skip Feature Selection", command=lambda: Skip_Feature_Selection(qbel_integr,output)).grid(column=0,row=6,sticky="W")
            else: #multi-sample i'll have to pass before_feature which is the most recent adata.
                Button(wind, text="Perform Selection",command=lambda: Feature_Select_Analysis(feature_report,ih_lela_report,qbel_feature,seurat_a,seurat_b,seurat_c,pearson_a,pearson_b,pearson_c,v3_a,output)).grid(column=4, row=6, sticky="E")
                Button(wind, text="Skip Feature Selection", command=lambda: Skip_Feature_Selection(qbel_feature,output)).grid(column=0, row=6, sticky="W")

def Feature_Select_Analysis(feat_rapport,ih_lela,adata,sirat_a,sirat_b,sirat_c,pirson_a,pirson_b,pirson_c,vtrois,outpute):
    if feat_rapport.get() =="":
        messagebox.showerror("Feature Info","You haven't choose any options\nPlease choose one of the 3 methods listed here")
    elif feat_rapport.get()=="pearson":
        #the pearson residuals uses the raw counts
        try:
            sc.experimental.pp.highly_variable_genes(adata, theta=pirson_a.get(), chunksize=pirson_c.get(), batch_key='sample', n_top_genes=pirson_b.get(),flavor='pearson_residuals')
        except ValueError:
            messagebox.showerror("Error Analysis","You didn't normalize your data\nso it's impossible to compute HVGs")
        plt.close("all") #close the figures before i'm thinking of running this line one time , because it will prevent two plt plot to be displayed.
        fig, ax = plt.subplots(1, 2, figsize=(20, 5),num="Pearson Residual Method")
        fig.tight_layout(w_pad=5, h_pad=5)
        plt.subplots_adjust(left=0.05, bottom=0.1)

        ax[0].scatter(adata.var['means'], adata.var['residual_variances'], c="lightgreen", s=5)
        ax[0].set_xlabel("mean expressions of genes")
        ax[0].set_ylabel("Pearson residual variance per gene")
        ax[1].scatter(adata.var['means'], [math.sqrt(x) for x in adata.var['variances'].values], c="salmon", s=5)
        ax[1].set_xlabel("mean expressions of genes")
        ax[1].set_ylabel("Standard deviation per gene")


        if ih_lela.get() =="non":
            print("report saved !")
            plt.savefig(f"{place}/{outpute.get().strip()}/filter_gene_dispersion_pearson.png")
        else:
            plt.show() # i don't block the execution of the code

        #adata = adata[:, adata.var.highly_variable] #if the use wants to do pca on all the genes they have to be available
        var_genes = adata.var.highly_variable
        global before_dimred_clust,feat_done
        feat_done=True
        #before_dimred_clust=adata #save the current adata to carry on the analysis.
        adata.write_h5ad("{0}/{1}/before_dimred_clust.h5ad".format(place,outpute.get().strip()))
        messagebox.showinfo("Progress Info","Feature Selection\nis done !!")

    elif feat_rapport.get()=="seurat":
        try:
            sc.pp.highly_variable_genes(adata, min_mean=sirat_a.get(), max_mean=sirat_b.get(), min_disp=sirat_c.get(),batch_key='sample')
        except ValueError:
            messagebox.showerror("Error Analysis","You didn't normalize your data\nso it's impossible to compute HVGs")

        plt.close("all")
        if ih_lela.get() == "non":
            print("report saved !")
            sc.pl.highly_variable_genes(adata,show=False,save=".png")
            os.system(f"cp {place}/figures/filter_genes_dispersion.png {place}/{outpute.get().strip()}") #i copy it into the output dir !
        else:
            sc.pl.highly_variable_genes(adata)
        #adata = adata[:, adata.var.highly_variable]
        var_genes = adata.var.highly_variable
        #i should see only this plot because i remove the others.
        #before_dimred_clust = adata
        adata.write_h5ad("{0}/{1}/before_dimred_clust.h5ad".format(place,outpute.get().strip()))
        feat_done=True
        messagebox.showinfo("Progress Info", "Feature Selection\nis done !!")

    elif feat_rapport.get()=="v3":
        try:
            sc.pp.highly_variable_genes(adata,batch_key='sample',n_top_genes=vtrois.get())
        except ValueError:
            messagebox.showerror("Error Analysis","You didn't normalize your data\nso it's impossible to compute HVGs")
        plt.close("all")

        if ih_lela == "non":
            sc.pl.highly_variable_genes(adata,show=False,save=".png")
            os.system(f"cp {place}/figures/filter_genes_dispersion.png {place}/{outpute.get().strip()}") #i copy it into the output dir !
        else:
            sc.pl.highly_variable_genes(adata)

        #adata = adata[:, adata.var.highly_variable]
        var_genes = adata.var.highly_variable
        #before_dimred_clust = adata
        adata.write_h5ad("{0}/{1}/before_dimred_clust.h5ad".format(place,outpute.get().strip()))
        feat_done=True
        messagebox.showinfo("Progress Info", "Feature Selection\nis done !!")


def Skip_Feature_Selection(adata,outpute):
    global before_dimred_clust
    #before_dimred_clust=adata
    adata.write_h5ad("{0}/{1}/before_dimred_clust.h5ad".format(place,outpute.get().strip()))
    global feat_done
    feat_done=True
    messagebox.showinfo("Info Analysis","You have skipped the\nfeature selection\nstep!!")

# </editor-fold>


# <editor-fold desc="DIM.REDUCION STEP">
def Reduc_clusterer(root,output):
    if okey_check==False:
        messagebox.showerror("Input Error","Please check your data by clicking on 'OKEY' before going on !")
    elif qc_done==False or norm_done==False or int_done==False or feat_done == False:
        Ftab1.progress_log(root,[qc_done,norm_done,int_done,feat_done],progress,num=4,heyt=250,step="Dim.Reduction/Clustering")
    else:
        wind = Toplevel(root)
        wind.title("DimReduc/Clustering")
        wind.geometry('+%d+%d' % (500, 250))

        pisiz = tk.StringVar()
        pisiz.set("50")

        hood = tk.StringVar()
        hood.set("15")

        koulech = tk.StringVar()
        clust_a = tk.StringVar()

        clust_b = tk.StringVar()
        clust_b.set("1.0")
        clust_c = tk.StringVar()
        clust_c.set("1.0")
        clust_d = tk.StringVar()
        clust_d.set("1.0/1.0")

        dimred = tk.StringVar()


        Label(wind, text="  Dimensionnality\nReduction", font=Font(root, size=11, weight=BOLD)).grid(column=0, row=0)
        Label(wind, text="PCA : ").grid(column=0, row=1)
        Radiobutton(wind, text="Only HVGs", variable=koulech, value="few").grid(column=0, row=2, sticky="W")
        Radiobutton(wind, text="All the genes", variable=koulech, value="all").grid(column=0, row=3, sticky="W")
        Label(wind, text="Compute \nneighbourhood", font=Font(root, size=11, weight=BOLD)).grid(column=0, row=4,
                                                                                                sticky="W")
        Label(wind, text="Number of PCs:").grid(column=0, row=5)
        Entry(wind, width=6, textvariable=pisiz).grid(column=0, row=6)
        Label(wind, text="Neighbourhood's size:").grid(column=0, row=7)
        Entry(wind, width=6, textvariable=hood).grid(column=0, row=8)
        Label(wind, text="Clustering\nMethod", font=Font(root, size=11, weight=BOLD)).grid(column=0, row=9)
        Radiobutton(wind, text="Leiden:", variable=clust_a, value="leiden").grid(column=0, row=10, sticky="w")
        Radiobutton(wind, text="Louvain:", variable=clust_a, value="louvain").grid(column=0, row=11, sticky="w")
        Radiobutton(wind, text="Both:", variable=clust_a, value="both").grid(column=0, row=12, sticky="w")
        Entry(wind, width=8, textvariable=clust_b).grid(column=0, row=10, sticky="E")
        Entry(wind, width=7, textvariable=clust_c).grid(column=0, row=11, sticky="E")
        Entry(wind, width=8, textvariable=clust_d).grid(column=0, row=12, sticky="E")
        Radiobutton(wind, text="both", variable=dimred, value="both", font=Font(root, size=7)).grid(column=0, row=13,
                                                                                                    sticky="W", pady=10)
        Radiobutton(wind, text="tsne", variable=dimred, value="tsne", font=Font(root, size=8)).grid(column=0, row=13,
                                                                                                    sticky="E")
        Radiobutton(wind, text="umap", variable=dimred, value="umap", font=Font(root, size=8)).grid(column=0, row=13)

        atoun = Image.open(f"{placa}/wecan.png")
        rat = ImageTk.PhotoImage(atoun)
        Button(wind, image=rat, fg="black", background="blue", width=50, height=400).place(x=170, y=0)

        toun = Image.open(f"{placa}/khanva.png")
        atr = ImageTk.PhotoImage(toun)
        Button(wind, image=atr, fg="black", background="white", width=50, height=400).place(x=220, y=0)

        oun = Image.open(f"{placa}/canvas.png")
        tra = ImageTk.PhotoImage(oun)
        Button(wind, image=tra, fg="black", background="red", width=50, height=400).place(x=275, y=0)

        kolor = tk.StringVar()
        kolor.set('geneA,sample')
        compo = tk.StringVar()
        compo.set("1,2/2,3")

        pcaavar = tk.BooleanVar()

        yumap = tk.StringVar()
        yumap.set('geneA,leiden_1.0')
        tisni = tk.StringVar()
        tisni.set('geneA,sample')
        stacked = tk.BooleanVar()


        pca_report = tk.BooleanVar()
        umap_report = tk.BooleanVar()
        tsne_report = tk.BooleanVar()

        Label(wind, text="Visualization", font=Font(root, size=11, weight=BOLD)).grid(column=3, row=0, padx=30,
                                                                                      sticky="w")
        Label(wind, text="Plot PCA", font=Font(root, size=9, weight=BOLD)).grid(column=3, row=1)
        Radiobutton(wind, text="yes", variable=pca_report, value=True).grid(column=3, row=2, sticky="w")
        Radiobutton(wind, text="No", variable=pca_report, value=False).grid(column=3, row=2, sticky="e")
        Label(wind, text="Color by :").grid(column=3, row=3, sticky="w")
        Entry(wind, width=12, textvariable=kolor).grid(column=3, row=3, sticky="e")
        Label(wind, text="N components:").grid(column=3, row=4, sticky="w")
        Entry(wind, width=8, textvariable=compo).grid(column=3, row=4, sticky="e")

        Label(wind, text="Plot Pca variance ratio", font=Font(root, size=9, weight=BOLD)).grid(column=3, row=5)
        Radiobutton(wind, text="yes", variable=pcaavar, value=True).grid(column=3, row=6, sticky="w")
        Radiobutton(wind, text="No", variable=pcaavar, value=False).grid(column=3, row=6, sticky="e")

        Label(wind, text="Plot UMAP", font=Font(root, size=9, weight=BOLD)).grid(column=3, row=7)
        Radiobutton(wind, text="yes", variable=umap_report, value=True).grid(column=3, row=8, sticky="w")
        Radiobutton(wind, text="No", variable=umap_report, value=False).grid(column=3, row=8, sticky="e")
        Label(wind, text="Color by :").grid(column=3, row=9, sticky="w")
        Entry(wind, width=12, textvariable=yumap).grid(column=3, row=9, sticky="e")
        Label(wind, text="Plot tSNE", font=Font(root, size=9, weight=BOLD)).grid(column=3, row=10, pady=5)
        Radiobutton(wind, text="yes", variable=tsne_report, value=True).grid(column=3, row=11, sticky="w")
        Radiobutton(wind, text="No", variable=tsne_report, value=False).grid(column=3, row=11, sticky="e")
        Label(wind, text="Color by :").grid(column=3, row=12, sticky="w")
        Entry(wind, width=12, textvariable=tisni).grid(column=3, row=12, columnspan=3, sticky="e")
        Label(wind, text="Display Cluster Report", font=Font(root, size=9, weight=BOLD)).grid(column=3, row=13)
        Radiobutton(wind, text="yes", variable=stacked, value=True).grid(column=3, row=14, sticky="w")
        Radiobutton(wind, text="No", variable=stacked, value=False).grid(column=3, row=14, sticky="e")

        qbel_dimred_clust=anndata.read_h5ad("{0}/{1}/before_dimred_clust.h5ad".format(place,output.get().strip()))
        print("Beging DimReduc_clusterization , i will work with this adata:",qbel_dimred_clust)
        #glob_adata=anndata.read_h5ad("{0}/global_adata.h5ad".format(os.getcwd()))
        Button(wind, text="Launch Clustering",command=lambda :DimReduc_Clustering_Analysis(
            qbel_dimred_clust,pisiz,hood,koulech,clust_a,clust_b,clust_c,clust_d,dimred,output)).grid(column=0, row=15, pady=5)
        Button(wind, text="Skip DimReduc/Clustering").grid(column=2, row=15, pady=5, sticky="E")
        Button(wind, text="Visualize Plots",command=lambda : Visualize_Clustering(root,anndata.read_h5ad("{0}/global_adata.h5ad".format(os.getcwd())),kolor,compo,pcaavar,yumap,tisni,stacked,pca_report,umap_report,tsne_report,clust_a,output)).grid(column=3, row=15, pady=5)


tmp1="make it global"
tmp2="make it global"
tmp="make it global"

def DimReduc_Clustering_Analysis(adata,pici,hodhod,kelch,clusta,clustb,clustc,clustd,dimblue,outpute):
    if pici.get()=="" or hodhod.get()=="" or kelch.get()=="" or clusta.get()=="" or dimblue.get() == "" or clustb.get()=="" or clustc.get()=="" or clustd.get()=="":
        messagebox.showwarning("Warning","You have to fill all the parameters\nin order to launch the analysis")
    else:
        flager=0 #if the flag is 0 then there is no error in the input format
        flagou=True
        try:
            sc.pp.pca(adata,use_highly_variable=kelch.get(),svd_solver='arpack',n_comps=50)
        except KeyError:
            messagebox.showerror("Error Info","Since you skipped the feature selection\nyou have no HVGs to compute\n")
            flagou=False

        if flagou:
            sc.pp.neighbors(adata, n_pcs=int(pici.get()), n_neighbors=int(hodhod.get()))
            global tmp
            tmp = {}
            if clusta.get() == "both":
                global tmp1, tmp2
                tmp1={}
                tmp2={}
                cluster=clustd.get().strip().split("/")
                for x in cluster[0].split(";"):
                    try:
                        sc.tl.leiden(adata,resolution=float(x) ,key_added="leiden_"+x)  # default resolution in 1.0
                        tmp1["leiden_"+x]=pd.crosstab(adata.obs['sample'], adata.obs['leiden_'+x]) #stock the cross tables
                    except ValueError:
                        messagebox.showerror("Error","Please correct the format of the input")
                        flager=1
                for y in cluster[1].split(";"):
                    try:
                        sc.tl.louvain(adata,resolution=float(y) ,key_added="louvain_"+y)
                        tmp2["louvain_" + y] = pd.crosstab(adata.obs['sample'], adata.obs['louvain_' +y])
                    except ValueError:
                        messagebox.showerror("Error", "Please correct the format of the input")
                        flager=1
            elif clusta.get() == "leiden":
                for x in clustb.get().split("/"):
                    try:
                        sc.tl.leiden(adata, resolution=float(x), key_added="leiden_" +x)
                        tmp["leiden_" +x]=pd.crosstab(adata.obs['sample'], adata.obs['leiden_'+x])
                    except ValueError:
                        messagebox.showerror("Error", "Please correct the format of the input")
                        flager=1
            elif clusta.get() == "louvain":
                for x in clustc.get().split("/"):
                   try:
                        sc.tl.louvain(adata,resolution=float(x),key_added="louvain_"+x)
                        tmp["louvain_" + x] = pd.crosstab(adata.obs['sample'], adata.obs['louvain_' + x])
                   except ValueError:
                        messagebox.showerror("Error", "Please correct the format of the input")
                        flager=1

            if flager==0:
                if dimblue.get() == "both":
                    sc.tl.umap(adata)
                    sc.tl.tsne(adata)
                elif dimblue.get() == "umap":
                    sc.tl.umap(adata)
                elif dimblue.get() == "tsne":
                    sc.tl.tsne(adata)
                #global adata_englobant
                #adata_englobant=adata #i'll pass this adata to the visualize plot although I could've used directly before_marker but , whatever
                adata.write_h5ad("{0}/global_adata.h5ad".format(os.getcwd()))
                global dim_clust_done,before_marker
                dim_clust_done=True
                #before_marker=adata
                adata.write_h5ad("{0}/{1}/before_marker.h5ad".format(place,outpute.get().strip()))
                messagebox.showinfo("Analysis Info","DimReduc/Clustering\ndone !!!")

def Visualize_Clustering(root,visu_adata,cooler,ocompo,pcastd,youmap,teasnea,stackover,pca_riport,umap_riport,tsne_riport,clusta_report,outpute):
    if cooler.get()=="" or ocompo.get()=="" or pcastd.get()=="" or youmap.get()=="" or teasnea.get() == "":
        messagebox.showwarning("Warning", "You have to fill all the parameters\nif you want to visualize")
        print("cooler {0} , ocompo {1} , pcavar {2} , youmap {3}, teasnea {4} , stackover {5} ".format(cooler.get(),ocompo.get(),pcastd.get(),youmap.get(),teasnea.get(),stackover.get()))
    if dim_clust_done == False:
        messagebox.showwarning("Warning", "Please run the DimReduc Clustering\nbefore visualizing")
    else:
        plt.close("all")
        cluter_sentence=[]#it will compile all the methods
        if pca_riport.get():
            try:
                sc.pl.pca(visu_adata,color=cooler.get().split(","),components=ocompo.get().split("/"))
            except KeyError:
                messagebox.showerror("Input Error","Either you miswrote one argument\nor\nthe key you are using doesn't exist !")
        if pcastd.get():
            sc.pl.pca_variance_ratio(visu_adata)
        if umap_riport.get():
            try:
                sc.pl.umap(visu_adata,color=youmap.get().split(",")) # i have to know which matrix i'm using.
            except KeyError:
                messagebox.showerror("Input Error","Either you miswrote one argument\nor\nthe key you are using doesn't exist !")
        if tsne_riport.get():
            try:
                sc.pl.tsne(visu_adata,color=teasnea.get().split(","))
            except KeyError:
                messagebox.showerror("Input Error", "Either you miswrote one argument\nor\nthe key you are using doesn't exist !")
        if stackover.get():
            #let's put some writen report.
            global giz
            giz=""
            if clusta_report.get() == "louvain":
                for x in tmp.keys():
                    giz = giz+"\t"+txt+"Clustering Report"+txt+"\n\n----Clustering Method: {0}\n\n----Cluster number: {1}\n\n----Total number of cell by cluster----{2}\n\n----Report for each Sample----{3}\n\n\n".format(
                    tmp[x].columns.name, len(tmp[x].columns),Ftab1.clust_report(tmp[x])[0],Ftab1.clust_report(tmp[x])[1])
                    cluter_sentence.append(x)
            elif clusta_report.get() == "leiden":
                for x in tmp.keys():
                    giz = giz+"\t"+txt+"Clustering Report"+txt+"\n\n----Clustering Method: {0}\n\n----Cluster number: {1}\n\n----Total number of cell by cluster----{2}\n\n----Report for each Sample----{3}\n\n\n".format(
                    tmp[x].columns.name, len(tmp[x].columns),Ftab1.clust_report(tmp[x])[0],Ftab1.clust_report(tmp[x])[1])
                    cluter_sentence.append(x)
            elif clusta_report.get() == "both":
                for x in tmp1.keys():
                    giz = giz + "\t"+txt+"Clustering Report"+txt+"\n\n----Clustering Method: {0}\n\n----Cluster number: {1}\n\n----Total number of cell by cluster----{2}\n\n----Report for each Sample----{3}\n\n\n".format(
                    tmp1[x].columns.name, len(tmp1[x].columns), Ftab1.clust_report(tmp1[x])[0], Ftab1.clust_report(tmp1[x])[1])
                    cluter_sentence.append(x)
                for x in tmp2.keys():
                    giz = giz + "\t"+txt+"Clustering Report"+txt+"\n\n----Clustering Method: {0}\n\n----Cluster number: {1}\n\n----Total number of cell by cluster----{2}\n\n----Report for each Sample----{3}\n\n\n".format(
                    tmp2[x].columns.name, len(tmp2[x].columns), Ftab1.clust_report(tmp2[x])[0], Ftab1.clust_report(tmp2[x])[1])
                    cluter_sentence.append(x)

            window=Toplevel(root)
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
            with open(f'{place}/{outpute.get().strip()}/clutering_report${"$".join(cluter_sentence)}', 'w') as f:
                f.write(giz)

            # Insert text into the text widget
            text_widget.insert(tk.END, giz)
            text_widget["state"] = DISABLED
# </editor-fold>


# <editor-fold desc="MARKER DETECTION STEP">
def Mark_Dectecter(root,output):
    if okey_check==False:
        messagebox.showerror("Input Error","Please check your data by clicking on 'OKEY' before going on !")
    elif qc_done==False or norm_done==False or int_done==False or feat_done == False or dim_clust_done==False:
        Ftab1.progress_log(root,[qc_done,norm_done,int_done,feat_done,dim_clust_done],progress,num=5,heyt=300,step="Marker Detection")
    else:
        wind = Toplevel(root)
        wind.title("Marker Detection")
        wind.geometry('+%d+%d' % (500, 250))

        churchill = tk.StringVar()
        correct = tk.StringVar()

        zefirst = tk.StringVar()

        zefirst.set("N")
        csv_report = tk.StringVar()
        cluster = tk.StringVar()

        bath=tk.StringVar()

        qbel_marker=anndata.read_h5ad("{0}/{1}/before_marker.h5ad".format(place,output.get().strip()))
        qbel_marker.uns['log1p']['base'] = None
        print("Begin Marker Detection , I will work with this adata (before marker):",qbel_marker)
        #print("I also want to see the raw adata of before marker:",qbel_marker.raw.to_adata().X[450:460,3000:3010])
        #print("I also want to see the adata before marker:",qbel_marker.X[450:460,3000:3010])

        Label(wind, text="Marker Gene Detection", font=Font(root, size=11, weight=BOLD)).grid(column=0, row=0, columnspan=3)
        Label(wind, text="Which Clustering ?", font=Font(root, size=11, weight=BOLD)).grid(column=0, row=1, sticky="w")
        Radiobutton(wind, text="Leiden", variable=cluster, value="leiden").grid(column=0, row=2)
        Radiobutton(wind, text="Louvain", variable=cluster, value="louvain").grid(column=1, row=2, sticky="w")
        Radiobutton(wind, text="Both", variable=cluster, value="both").grid(column=2, row=2, sticky="w")

        Label(wind, text="Statistical Test ", font=Font(root, size=11, weight=BOLD)).grid(column=0, row=3, sticky="w")
        Radiobutton(wind, text="Student t-test", variable=churchill, value="t-test").grid(column=0, row=4)
        Radiobutton(wind, text="Logistic Regression", variable=churchill, value="logreg").grid(column=1, row=4, sticky="w")
        Radiobutton(wind, text="Wilcoxon rank-sum test", variable=churchill, value="wilcoxon").grid(column=2, row=4, sticky="e")

        Label(wind, text="Correction Method", font=Font(root, size=11, weight=BOLD)).grid(column=0, row=5, sticky="W", pady=5)
        Radiobutton(wind, text="Benjamini-Hocheberg", variable=correct, value="benjamini-hochberg").grid(column=0, row=6, sticky="w")
        Radiobutton(wind, text="Bonferroni", variable=correct, value="bonferroni").grid(column=1, row=6, sticky="w")
        Label(wind, text="Report Section", font=Font(root, size=11, weight=BOLD)).grid(column=0, row=7, sticky="W", pady=10)
        Radiobutton(wind, text="Show plots", variable=bath, value="show").grid(column=0, row=8, sticky="w")
        Radiobutton(wind, text="Hide plots", variable=bath, value="hide").grid(column=1, row=8, sticky="w")
        Radiobutton(wind, text="Save all genes in csv", variable=csv_report, value="all").grid(column=0, row=9, pady=10)
        Radiobutton(wind, text="Save csv file, the first N genes", variable=csv_report, value="not_all").grid(column=1, row=9,pady=10,sticky="e")
        Radiobutton(wind, text="No csv", variable=csv_report, value="No").grid(column=2, row=9,pady=10,sticky="e")

        Entry(wind, textvariable=zefirst, width=13).grid(column=2, row=9, sticky="w")

        Button(wind, text="Detect Markers",command=lambda :Marker_Detection_Analysis(root,qbel_marker,churchill,correct,zefirst,csv_report,cluster,bath,output)).grid(column=0, row=10, pady=5)
        Button(wind, text="Skip Marker Detection",command=lambda : Skip_Marker_Detection(root,qbel_marker,output)).grid(column=1, row=10, pady=5, sticky="E")


def Marker_Detection_Analysis(root,mark_adata,church,korrect,zesecond,csvi,bluster,barth,outpute):
    flag = 0
    if church.get() == "" or korrect.get() == ""  or bluster.get() == "" or barth.get() == "" or csvi.get()=="":
        messagebox.showerror("Warning", "Please fill all the parameters\nbefore you proceed")
    elif (church.get() == "" or korrect.get() == ""  or bluster.get() == "" or barth.get() == "") or (csvi.get()=="not_all" and zesecond.get()==""):
        messagebox.showerror("Warning", "Please fill all the parameters\nbefore you proceed")
    else:
        global mark_done
        if bluster.get() == "louvain":
            lista = []
            plt.close("all")
            for x in mark_adata.obs.columns:
                if x.__contains__("louvain"):
                    lista.append(x)

            answer1 = simpledialog.askstring("Input",
                                             "You have one or multiple resolution for the louvain method:\n\n{0}\n\nIf you want to work with one resolution specify it in the input\n"
                                             "If you want to work with multiple resolutions specify it in the input with this format\n\nmethod1/method2\n\nIf you to work with them all,then type ALL".format(
                                                 lista), parent=root)

            if answer1 == "":
                messagebox.showerror("Error Input",
                                     "You haven't specify any resolution , you can't proceed to the analysis.")

            elif answer1.strip() == "ALL":  # the use want to marker detect on all of his clustering methods
                lista
                for x in lista:
                    if 'Scanorama' in mark_adata.obsm.keys():
                        sc.tl.rank_genes_groups(mark_adata, x.strip(), method=church.get(), corr_method=korrect.get())
                    else:
                        sc.tl.rank_genes_groups(mark_adata, x.strip(), method=church.get(), corr_method=korrect.get(),
                                                use_raw=True)

                    if csvi.get() != "No":
                        csv_builder(mark_adata, csvi, zesecond, x.strip(), flag=0, cluster=["nan", "nan"],outofput=f"{place}/{outpute.get().strip()}")
                    else:
                        print("No Csv generated")

                    if barth.get() == "show":
                        sc.pl.rank_genes_groups(mark_adata, n_genes=15, sharey=False, show=False)
                        sc.pl.rank_genes_groups_matrixplot(mark_adata, n_genes=5, groupby=x.strip(), title=x.strip(),
                                                           show=False)
                    else:
                        print("reports saved ! ")
                        sc.pl.rank_genes_groups(mark_adata, n_genes=15, sharey=False, show=False, save=".png")
                        os.system(
                            f"cp {place}/figures/rank_genes_groups_{x.strip()}.png {place}/{outpute.get().strip()}/")  # copy it in the output dir
                        plt.close("all")  # because the show is false so just in case , i clean the plt cache
                        sc.pl.rank_genes_groups_matrixplot(mark_adata, n_genes=5, groupby=x.strip(), title=x.strip(),
                                                           show=False)
                        plt.savefig(f"{place}/{outpute.get().strip()}/matrixplot_{x.strip()}.png")
                        plt.close("all")
                mark_done = True
                mark_adata.write_h5ad("{0}/end_adata.h5ad".format(f"{place}/{outpute.get().strip()}"))
                os.system(f"cp -r {place}/{outpute.get().strip()} $HOME/Desktop")
                messagebox.showinfo("Progress Analysis","Marker Detection Done !!!\ndata available in \n{0}/{1}/end_adata.h5ad".format("Desktop", outpute.get().strip()))

            elif len(answer1.strip().split("/")) > 1 and answer1.strip() != "ALL":
                for x in answer1.split("/"):
                    if x.strip() not in lista:
                        messagebox.showerror("Error Input", f"This resolution {x.strip()} do not exist in your data")
                        flag = 1
                        break

            if len(answer1.split(
                    "/")) > 1 and flag == 0:  # i have to put in if form (rather than elif) to force the program to evaluate it.
                for x in answer1.split("/"):
                    if 'Scanorama' in mark_adata.obsm.keys():
                        sc.tl.rank_genes_groups(mark_adata, x.strip(), method=church.get(), corr_method=korrect.get())
                    else:
                        sc.tl.rank_genes_groups(mark_adata, x.strip(), method=church.get(), corr_method=korrect.get(),
                                                use_raw=True)

                    if csvi.get() != "No":
                        csv_builder(mark_adata, csvi, zesecond, x.strip(), flag=0, cluster=["nan", "nan"],outofput=f"{place}/{outpute.get().strip()}")
                    else:
                        print("No Csv generated")

                    if barth.get() == "show":
                        sc.pl.rank_genes_groups(mark_adata, n_genes=15, sharey=False, show=False)
                        sc.pl.rank_genes_groups_matrixplot(mark_adata, n_genes=5, groupby=x.strip(), title=x.strip(),
                                                           show=False)
                    else:
                        print("reports saved ! ")
                        sc.pl.rank_genes_groups(mark_adata, n_genes=15, sharey=False, show=False, save=".png")
                        os.system(
                            f"cp {place}/figures/rank_genes_groups_{x.strip()}.png {place}/{outpute.get().strip()}/")  # copy it in the output dir
                        plt.close("all")  # because the show is false so just in case , i clean the plt cache
                        sc.pl.rank_genes_groups_matrixplot(mark_adata, n_genes=5, groupby=x.strip(), title=x.strip(),
                                                           show=False)
                        plt.savefig(f"{place}/{outpute.get().strip()}/matrixplot_{x.strip()}.png")
                        plt.close("all")
                mark_done = True
                mark_adata.write_h5ad("{0}/end_adata.h5ad".format(f"{place}/{outpute.get().strip()}"))
                os.system(f"cp -r {place}/{outpute.get().strip()} $HOME/Desktop")
                messagebox.showinfo("Progress Analysis","Marker Detection Done !!!\ndata available in \n{0}/{1}/end_adata.h5ad".format("Desktop", outpute.get().strip()))

            elif len(answer1.split("/")) == 1 and flag == 0 and answer1.strip() != "ALL":  # the use entered just one method
                if answer1.strip() not in lista:
                    messagebox.showerror("Error Input","This resolution {0} do not exist in your data".format(answer1.strip()))
                else:
                    if 'Scanorama' in mark_adata.obsm.keys():
                        sc.tl.rank_genes_groups(mark_adata, answer1.strip(), method=church.get(),
                                                corr_method=korrect.get())
                    else:
                        sc.tl.rank_genes_groups(mark_adata, answer1.strip(), method=church.get(),
                                                corr_method=korrect.get(), use_raw=True)

                    if csvi.get() != "No":
                        csv_builder(mark_adata, csvi, zesecond, answer1.strip(), flag=0, cluster=["nan", "nan"],outofput=f"{place}/{outpute.get().strip()}")
                    else:
                        print("No Csv generated")

                    if barth.get() == "show":
                        sc.pl.rank_genes_groups(mark_adata, n_genes=15, sharey=False, show=False)
                        sc.pl.rank_genes_groups_matrixplot(mark_adata, n_genes=5, groupby=answer1.strip(),title=answer1.strip(), show=False)
                    else:
                        print("reports saved ! ")
                        sc.pl.rank_genes_groups(mark_adata, n_genes=15, sharey=False, show=False, save=".png")
                        os.system(f"cp {place}/figures/rank_genes_groups_{answer1.strip()}.png {place}/{outpute.get().strip()}/")  # copy it in the output dir
                        plt.close("all")  # because the show is false so just in case , i clean the plt cache
                        sc.pl.rank_genes_groups_matrixplot(mark_adata, n_genes=5, groupby=answer1.strip(), title=answer1.strip(),show=False)
                        plt.savefig(f"{place}/{outpute.get().strip()}/matrixplot_{answer1.strip()}.png")
                        plt.close("all") # i use this here because at the end of the function I have a plt.show() so i clean all the plt cache so that plt.show() is useless
                    mark_done = True
                    mark_adata.write_h5ad("{0}/end_adata.h5ad".format(f"{place}/{outpute.get().strip()}"))
                    os.system(f"cp -r {place}/{outpute.get().strip()} $HOME/Desktop")
                    messagebox.showinfo("Progress Analysis","Marker Detection Done !!!\ndata available in \n{0}/{1}/end_adata.h5ad".format("Desktop", outpute.get().strip()))
            plt.show()  # to show the plot

        elif bluster.get() == "leiden":
            plt.close("all")
            lista = []
            for x in mark_adata.obs.columns:
                if x.__contains__("leiden"):
                    lista.append(x)

            answer1 = simpledialog.askstring("Input",
                                             "You have one or multiple resolution for the leiden method:\n\n{0}\n\nIf you want to work with one resolution specify it in the input\n"
                                             "If you want to work with multiple resolutions specify it in the input with this format\n\nmethod1/method2\n\nIf you to work with them all,then type ALL".format(
                                                 lista), parent=root)
            if answer1 == "":
                messagebox.showerror("Error Input",
                                     "You haven't specify any resolution , you can't proceed to the analysis.")

            elif answer1.strip() == "ALL":  # the use want to marker detect on all of his clustering methods
                lista
                for x in lista:
                    if 'Scanorama' in mark_adata.obsm.keys():
                        sc.tl.rank_genes_groups(mark_adata, x.strip(), method=church.get(), corr_method=korrect.get())
                    else:
                        sc.tl.rank_genes_groups(mark_adata, x.strip(), method=church.get(), corr_method=korrect.get(),
                                                use_raw=True)

                    if csvi.get() != "No":
                        csv_builder(mark_adata, csvi, zesecond, x.strip(), flag=0, cluster=["nan", "nan"],outofput=f"{place}/{outpute.get().strip()}")
                    else:
                        print("No Csv generated")

                    if barth.get() == "show":
                        sc.pl.rank_genes_groups(mark_adata, n_genes=15, sharey=False, show=False)
                        sc.pl.rank_genes_groups_matrixplot(mark_adata, n_genes=5, groupby=x.strip(), title=x.strip(),show=False)
                    else:
                        print("reports saved ! ")
                        sc.pl.rank_genes_groups(mark_adata, n_genes=15, sharey=False, show=False,save=".png")
                        os.system(f"cp {place}/figures/rank_genes_groups_{x.strip()}.png {place}/{outpute.get().strip()}/") #copy it in the output dir
                        plt.close("all") #because the show is false so just in case , i clean the plt cache
                        sc.pl.rank_genes_groups_matrixplot(mark_adata, n_genes=5, groupby=x.strip(), title=x.strip(),show=False)
                        plt.savefig(f"{place}/{outpute.get().strip()}/matrixplot_{x.strip()}.png")
                        plt.close("all")
                mark_done = True
                mark_adata.write_h5ad("{0}/end_adata.h5ad".format(f"{place}/{outpute.get().strip()}"))
                os.system(f"cp -r {place}/{outpute.get().strip()} $HOME/Desktop")
                messagebox.showinfo("Progress Analysis","Marker Detection Done !!!\ndata available in \n{0}/{1}/end_adata.h5ad".format("Desktop", outpute.get().strip()))

            elif len(answer1.strip().split("/")) > 1 and answer1.strip() != "ALL":
                for x in answer1.split("/"):
                    if x.strip() not in lista:
                        messagebox.showerror("Error Input",
                                             "This resolution {0} do not exist in your data".format(x.strip()))
                        flag = 1
                        break

            if len(answer1.split("/")) > 1 and flag == 0:
                for x in answer1.split("/"):
                    if 'Scanorama' in mark_adata.obsm.keys():
                        sc.tl.rank_genes_groups(mark_adata, x.strip(), method=church.get(), corr_method=korrect.get())
                    else:
                        sc.tl.rank_genes_groups(mark_adata, x.strip(), method=church.get(), corr_method=korrect.get(),use_raw=True)


                    if csvi.get() != "No":
                        csv_builder(mark_adata, csvi, zesecond, x.strip(), flag=0, cluster=["nan", "nan"],outofput=f"{place}/{outpute.get().strip()}")
                    else:
                        print("No Csv generated")

                    if barth.get() == "show":
                        sc.pl.rank_genes_groups(mark_adata, n_genes=15, sharey=False, show=False)
                        sc.pl.rank_genes_groups_matrixplot(mark_adata, n_genes=5, groupby=x.strip(), title=x.strip(),
                                                           show=False)
                    else:
                        print("reports saved ! ")
                        sc.pl.rank_genes_groups(mark_adata, n_genes=15, sharey=False, show=False, save=".png")
                        os.system(f"cp {place}/figures/rank_genes_groups_{x.strip()}.png {place}/{outpute.get().strip()}/")  # copy it in the output dir
                        plt.close("all")  # because the show is false so just in case , i clean the plt cache
                        sc.pl.rank_genes_groups_matrixplot(mark_adata, n_genes=5, groupby=x.strip(), title=x.strip(),show=False)
                        plt.savefig(f"{place}/{outpute.get().strip()}/matrixplot_{x.strip()}.png")
                        plt.close("all") # i use this here because at the end of the function I have a plt.show() so i clean all the plt cache so that plt.show() is useless
                mark_done = True
                mark_adata.write_h5ad("{0}/end_adata.h5ad".format(f"{place}/{outpute.get().strip()}"))
                os.system(f"cp -r {place}/{outpute.get().strip()} $HOME/Desktop")
                messagebox.showinfo("Progress Analysis","Marker Detection Done !!!\ndata available in \n{0}/{1}/end_adata.h5ad".format("Desktop", outpute.get().strip()))

            elif len(answer1.strip().split("/")) == 1 and flag == 0 and answer1.strip() != "ALL":
                if answer1.strip() not in lista:
                    messagebox.showerror("Error Input",
                                         "This resolution {0} do not exist in your data".format(answer1.strip()))
                else:
                    if 'Scanorama' in mark_adata.obsm.keys():  # no need to use raw because scanorama don't touch the X matrix .
                        sc.tl.rank_genes_groups(mark_adata, answer1.strip(), method=church.get(),
                                                corr_method=korrect.get())
                    else:
                        sc.tl.rank_genes_groups(mark_adata, answer1.strip(), method=church.get(),corr_method=korrect.get(), use_raw=True)

                    if csvi.get() != "No":
                        csv_builder(mark_adata, csvi, zesecond, answer1.strip(), flag=0, cluster=["nan", "nan"],outofput=f"{place}/{outpute.get().strip()}")
                    else:
                        print("No Csv generated")

                    if barth.get() == "show":
                        sc.pl.rank_genes_groups(mark_adata, n_genes=15, sharey=False, show=False)
                        sc.pl.rank_genes_groups_matrixplot(mark_adata, n_genes=5, groupby=answer1.strip(),
                                                           title=answer1.strip(), show=False)
                    else:
                        print("reports saved ! ")
                        sc.pl.rank_genes_groups(mark_adata, n_genes=15, sharey=False, show=False, save=".png")
                        os.system(
                            f"cp {place}/figures/rank_genes_groups_{answer1.strip()}.png {place}/{outpute.get().strip()}/")  # copy it in the output dir
                        plt.close("all")  # because the show is false so just in case , i clean the plt cache
                        sc.pl.rank_genes_groups_matrixplot(mark_adata, n_genes=5, groupby=answer1.strip(), title=x.strip(),
                                                           show=False)
                        plt.savefig(f"{place}/{outpute.get().strip()}/matrixplot_{answer1.strip()}.png")
                        plt.close("all")
                    mark_done = True
                    mark_adata.write_h5ad("{0}/end_adata.h5ad".format(f"{place}/{outpute.get().strip()}"))
                    os.system(f"cp -r {place}/{outpute.get().strip()} $HOME/Desktop")
                    messagebox.showinfo("Progress Analysis","Marker Detection Done !!!\ndata available in \n{0}/{1}/end_adata.h5ad".format("Desktop", outpute.get().strip()))

            plt.show()

        elif bluster.get() == "both":
            plt.close("all")
            lista = []
            for x in mark_adata.obs.columns:
                if x.__contains__("leiden") or x.__contains__("louvain"):
                    lista.append(x)

            answer1 = simpledialog.askstring("Input",
                                             "You have one or multiple resolution for the leiden/louvain method:\n\n{0}\n\nIf you want to work with one resolution specify it in the input\n"
                                             "If you want to work with multiple resolutions specify it in the input with this format\n\nmethod1/method2\n\nIf you to work with them all,then type ALL".format(
                                                 lista), parent=root)
            if answer1 == "":
                messagebox.showerror("Error Input",
                                     "You haven't specify any resolution , you can't proceed to the analysis.")

            elif answer1.strip() == "ALL":  # the use want to marker detect on all of his clustering methods
                lista
                for x in lista:
                    if 'Scanorama' in mark_adata.obsm.keys():
                        sc.tl.rank_genes_groups(mark_adata, x.strip(), method=church.get(), corr_method=korrect.get())
                    else:
                        sc.tl.rank_genes_groups(mark_adata, x.strip(), method=church.get(), corr_method=korrect.get(),
                                                use_raw=True)

                    if csvi.get() != "No":
                        csv_builder(mark_adata, csvi, zesecond, x.strip(), flag=0, cluster=["nan", "nan"],outofput=f"{place}/{outpute.get().strip()}")
                    else:
                        print("No Csv generated")

                    if barth.get() == "show":
                        sc.pl.rank_genes_groups(mark_adata, n_genes=15, sharey=False, show=False)
                        sc.pl.rank_genes_groups_matrixplot(mark_adata, n_genes=5, groupby=x.strip(), title=x.strip(),show=False)
                    else:
                        print("reports saved ! ")
                        sc.pl.rank_genes_groups(mark_adata, n_genes=15, sharey=False, show=False, save=".png")
                        os.system(f"cp {place}/figures/rank_genes_groups_{x.strip()}.png {place}/{outpute.get().strip()}/")  # copy it in the output dir
                        plt.close("all")  # because the show is false so just in case , i clean the plt cache
                        sc.pl.rank_genes_groups_matrixplot(mark_adata, n_genes=5, groupby=x.strip(), title=x.strip(),show=False)
                        plt.savefig(f"{place}/{outpute.get().strip()}/matrixplot_{x.strip()}.png")
                        plt.close("all")
                mark_done = True
                mark_adata.write_h5ad("{0}/end_adata.h5ad".format(f"{place}/{outpute.get().strip()}"))
                os.system(f"cp -r {place}/{outpute.get().strip()} $HOME/Desktop")
                messagebox.showinfo("Progress Analysis","Marker Detection Done !!!\ndata available in \n{0}/{1}/end_adata.h5ad".format("Desktop", outpute.get().strip()))

            elif len(answer1.split("/")) > 1 and answer1.strip() != "ALL":
                for x in answer1.split("/"):
                    if x.strip() not in lista:
                        messagebox.showerror("Error Input", f"This resolution {x.strip()} do not exist in your data")
                        flag = 1
                        break

            if len(answer1.split("/")) > 1 and flag == 0:
                for x in answer1.split("/"):
                    if 'Scanorama' in mark_adata.obsm.keys():  # no need to use raw in scanorama it will cause an error
                        sc.tl.rank_genes_groups(mark_adata, x.strip(), method=church.get(), corr_method=korrect.get())
                    else:
                        sc.tl.rank_genes_groups(mark_adata, x.strip(), method=church.get(), corr_method=korrect.get(),use_raw=True)

                    if csvi.get() != "No":
                        csv_builder(mark_adata, csvi, zesecond, x.strip(), flag=0, cluster=["nan", "nan"],outofput=f"{place}/{outpute.get().strip()}")
                    else:
                        print("No Csv generated")

                    if barth.get() == "show":
                        sc.pl.rank_genes_groups(mark_adata, n_genes=15, sharey=False, show=False)
                        sc.pl.rank_genes_groups_matrixplot(mark_adata, n_genes=5, groupby=x.strip(), title=x.strip(),show=False)
                    else:
                        print("reports saved ! ")
                        sc.pl.rank_genes_groups(mark_adata, n_genes=15, sharey=False, show=False, save=".png")
                        os.system(f"cp {place}/figures/rank_genes_groups_{x.strip()}.png {place}/{outpute.get().strip()}/")  # copy it in the output dir
                        plt.close("all")  # because the show is false so just in case , i clean the plt cache
                        sc.pl.rank_genes_groups_matrixplot(mark_adata, n_genes=5, groupby=x.strip(), title=x.strip(),show=False)
                        plt.savefig(f"{place}/{outpute.get().strip()}/matrixplot_{x.strip()}.png")
                        plt.close("all")

                mark_done = True
                mark_adata.write_h5ad("{0}/end_adata.h5ad".format(f"{place}/{outpute.get().strip()}"))
                os.system(f"cp -r {place}/{outpute.get().strip()} $HOME/Desktop")
                messagebox.showinfo("Progress Analysis","Marker Detection Done !!!\ndata available in \n{0}/{1}/end_adata.h5ad".format("Desktop", outpute.get().strip()))

            elif len(answer1.split("/")) == 1 and flag == 0 and answer1.strip() != "ALL":
                if answer1.strip() not in lista:
                    messagebox.showerror("Error Input",
                                         "This resolution {0} do not exist in your data".format(answer1.strip()))
                else:
                    if 'Scanorama' in mark_adata.obsm.keys():
                        sc.tl.rank_genes_groups(mark_adata, answer1.strip(), method=church.get(),corr_method=korrect.get())
                    else:
                        sc.tl.rank_genes_groups(mark_adata, answer1.strip(), method=church.get(),corr_method=korrect.get(), use_raw=True)

                    if csvi.get() != "No":
                        csv_builder(mark_adata, csvi, zesecond, answer1.strip(), flag=0, cluster=["nan", "nan"],outofput=f"{place}/{outpute.get().strip()}")
                    else:
                        print("No Csv generated")

                    if barth.get() == "show":
                        sc.pl.rank_genes_groups(mark_adata, n_genes=15, sharey=False, show=False)
                        sc.pl.rank_genes_groups_matrixplot(mark_adata, n_genes=5, groupby=answer1.strip(),title=answer1.strip(), show=False)
                    else:
                        print("reports saved ! ")
                        sc.pl.rank_genes_groups(mark_adata, n_genes=15, sharey=False, show=False, save=".png")
                        os.system(f"cp {place}/figures/rank_genes_groups_{answer1.strip()}.png {place}/{outpute.get().strip()}/")  # copy it in the output dir
                        plt.close("all")  # because the show is false so just in case , i clean the plt cache
                        sc.pl.rank_genes_groups_matrixplot(mark_adata, n_genes=5, groupby=answer1.strip(), title=answer1.strip(),show=False)
                        plt.savefig(f"{place}/{outpute.get().strip()}/matrixplot_{answer1.strip()}.png")
                        plt.close("all")
                    mark_done = True
                    mark_adata.write_h5ad("{0}/{1}/end_adata.h5ad".format(place,outpute.get().strip())) #had la ligne est différente dans docker.
                    os.system(f"cp -r {place}/{outpute.get().strip()} $HOME/Desktop")
                    messagebox.showinfo("Progress Analysis",
                                        "Marker Detection Done !!!\ndata available in \n{0}/{1}/end_adata.h5ad".format(
                                            "Desktop", outpute.get().strip()))
            plt.show()

def csv_builder(mark_adata,csvi,zethird,key,flag,cluster,outofput):
    lista = mark_adata.uns['rank_genes_groups']['scores'].tolist()
    listb = mark_adata.uns['rank_genes_groups']['pvals'].tolist()
    listc = mark_adata.uns['rank_genes_groups']['pvals_adj'].tolist()
    listd = mark_adata.uns['rank_genes_groups']['logfoldchanges'].tolist()
    liste = mark_adata.uns['rank_genes_groups']['names'].tolist()

    listo = []
    for w, x, y, z , p  in zip(lista, listb, listc, listd,liste):
        listing = []
        for v in range(0, len(x)):
            listing.append(w[v])
            listing.append(x[v])
            listing.append(y[v])
            listing.append(z[v])
            listing.append(p[v])
            listing.append(" ") #to put spaces between each cluster
        listo.append(listing)

    leicester = []
    if flag==1:
        for i in range(0, len(mark_adata.uns['rank_genes_groups']['scores'][0])):
            lister = ['scores', 'p-value', 'adjusted_p-value', 'LogFC', 'marker_cluster_',' ']
            lister[4] = lister[4] + cluster[i]
            leicester.append(lister)
    else:
        for i in range(0, len(mark_adata.obs[key].cat.categories)):
            lister = ['scores', 'p-value', 'adjusted_p-value','LogFC','marker_cluster_',' ']
            lister[4] = lister[4] + str(i)
            leicester.append(lister)

    if csvi.get() == "all":
        pd.DataFrame(listo, columns=list(chain(*leicester))).to_csv("{0}/all_genes_{1}.csv".format(outofput,key), index=False)
        messagebox.showinfo("Progress Analysis", "CSV CREATED!\n{0}/all_genes_{1}.csv".format(outofput,key))

    elif csvi.get() == "not_all":
        try:
            pd.DataFrame(listo, columns=list(chain(*leicester))).head(int(zethird.get().strip())).to_csv("{0}/{1}_genes_{2}.csv".format(outofput,zethird.get().strip(),key), index=False)
        except ValueError:
            messagebox.showerror("Error Format","Please provide numbers in this box")
        messagebox.showinfo("Progress Analysis", "CSV CREATED!\n{0}/{1}_genes_{2}.csv".format(outofput,zethird.get().strip(),key))

def Skip_Marker_Detection(root,adata,outpute):
    messagebox.showinfo("Info","You have skipped the\nmarker detection Step !")
    os.system(f"cp -r {place}/{outpute.get().strip()} $HOME/Desktop")
    #mark_done remains FALSE , to block the cluster compare analysis
    #WHICH ADATA ???

# </editor-fold>

# <editor-fold desc="CLUSTER COMPARAISON STEP">
def Comparer(root,output):
    if okey_check==False:
        messagebox.showerror("Input Error","Please check your data by clicking on 'OKEY' before going on !")
    elif qc_done==False or norm_done==False or int_done==False or feat_done == False or dim_clust_done==False or mark_done==False:
        Ftab1.progress_log(root,[qc_done,norm_done,int_done,feat_done,dim_clust_done,mark_done],progress,num=6,heyt=340,step="Cluster Comparaison")

    else:
        wind = Toplevel(root)
        wind.title("Cluster Analysis")
        wind.geometry('+%d+%d' % (1000, 250))

        masgid = tk.StringVar()
        sahih = tk.StringVar()

        group = tk.StringVar()

        tsv = tk.StringVar()

        ref = tk.StringVar()

        douche = tk.StringVar()

        awal=tk.StringVar()
        awal.set("N")

        group.set("5,2,3")



        ref.set("1")

        try:
            after_marker=anndata.read_h5ad(f"{place}/{output.get().strip()}/end_adata.h5ad".format()) #le end_adata rahou yezghoud !! faut le capter mais attention kayen mode pycharm launch et desktop launch
            after_marker.uns['log1p']['base'] = None
            drapeau=0
        except FileNotFoundError:
            messagebox.showerror("Error Input","Something wrong happened\n the data are lost !!\nyou can't do this analysis")
            drapeau=1

        if drapeau==0:
            Label(wind, text="Compare Clusters", font=Font(root, size=13, weight=BOLD)).grid(column=0, row=0, columnspan=3)
            Label(wind, text="Cluster(s)", font=Font(root, size=11, weight=BOLD)).grid(column=0, row=2, sticky="w", pady=10)
            Entry(wind, textvariable=group, width=10).grid(column=0, row=2, sticky="e")
            Label(wind, text="VS Reference", font=Font(root, size=11, weight=BOLD)).grid(column=1, row=2, sticky="w")
            Entry(wind, textvariable=ref, width=10).grid(column=1, row=2, sticky="e")

            Label(wind, text="Statistical Test ", font=Font(root, size=11, weight=BOLD)).grid(column=0, row=3, sticky="w")
            Radiobutton(wind, text="Student t-test", variable=masgid, value="t-test").grid(column=0, row=4)
            Radiobutton(wind, text="Logistic Regression", variable=masgid, value="logreg").grid(column=1, row=4, sticky="w")
            Radiobutton(wind, text="Wilcoxon rank-sum test", variable=masgid, value="wilcoxon").grid(column=2, row=4,sticky="e")

            Label(wind, text="Correction Method", font=Font(root, size=11, weight=BOLD)).grid(column=0, row=5, sticky="W", pady=5)

            Radiobutton(wind, text="Benjamini-Hocheberg", variable=sahih, value="benjamini-hochberg").grid(column=0, row=6,sticky="w")

            Radiobutton(wind, text="Bonferroni", variable=sahih, value="bonferroni").grid(column=1, row=6, sticky="w")

            Label(wind, text="Report Section", font=Font(root, size=11, weight=BOLD)).grid(column=0, row=7, sticky="W", pady=10)
            Radiobutton(wind, text="Show plots", variable=douche, value="show").grid(column=0, row=8, sticky="w")
            Radiobutton(wind, text="Hide plots", variable=douche, value="hide").grid(column=1, row=8, sticky="w")

            Radiobutton(wind, text="Save all genes in csv", variable=tsv, value="all").grid(column=0, row=9, pady=10)
            Radiobutton(wind, text="Save csv file, the first N genes :", variable=tsv, value="not_all").grid(column=1,row=9,pady=10,sticky="e")
            Entry(wind, textvariable=awal, width=13).grid(column=2, row=9, sticky="w")

            Radiobutton(wind, text="No csv ", variable=tsv, value="No").grid(column=2,row=9,sticky="e")

            Button(wind, text="Compare Clusters",command=lambda : Cluster_Compare_Analysis(root,after_marker,group,ref,masgid,sahih,awal,tsv,douche,output)).grid(column=0, row=10, pady=5)
            Button(wind, text="Skip Comparaison").grid(column=1, row=10, pady=5, sticky="E")

def Cluster_Compare_Analysis(root,clust_adata,fariq,reffe,cherch,haqiq,thani,tisvi,bakht,output):
    if cherch.get() == "" or haqiq.get() == "" or thani.get() == "" or bakht.get() == "" or tisvi.get() == "" or fariq.get()=="" or reffe.get()=="" or len(reffe.get().split(",")) >1 :
        messagebox.showerror("Warning", "Please fill (correctly) all the parameters\nbefore you proceed")
    elif (cherch.get() == "" or haqiq.get() == "" or fariq.get()=="" or bakht.get() == "" or fariq.get()=="" or reffe.get()=="" ) or (tisvi.get() == "not_all" and thani.get() == ""):
        messagebox.showerror("Warning", "Please fill all the parameters\nbefore you proceed")
    else:

        lista = []
        for x in clust_adata.obs.columns:
            if x.__contains__("leiden") or x.__contains__("louvain"):
                lista.append(x)
        idjaba= simpledialog.askstring("Input","You have one or multiple resolution for the leiden/louvain method:\n\n{0}\n\nFrom which method do you want to compare your clusters ?\n"
                                       "Please specify a method:".format(lista), parent=root)
        try:
            sc.tl.rank_genes_groups(clust_adata,groupby=idjaba.strip(),groups=fariq.get().strip().split(","),reference=reffe.get().strip(),corr_method=haqiq.get(),method=cherch.get())
            flager=0
        except IndexError:
            messagebox.showerror("Error Input","Please specify a correct cluster name\n one the cluster don't exist")
            flager=1
        except ValueError:
            messagebox.showerror("Reference Error","The reference should be\n from among the clusters")
            flager=1
        except : #other errors
            messagebox.showerror("Error Input","Something went wrong\n your input is not in the \n correct format")
            flager=1

        if flager==0: #we can continue the analysis
            if tisvi.get() != "No":
                csv_builder(clust_adata,tisvi,thani,key="{0}_vs_{1}".format(fariq.get().strip().split(","),reffe.get().strip()),flag=1,cluster=fariq.get().strip().split(","),outofput=f"{place}/{output.get().strip()}")
            else:
                print("No Csv generated")
            if bakht.get() == "show":
                plt.close("all")
                sc.pl.rank_genes_groups(clust_adata, groups=fariq.get().strip().split(","), n_genes=20)
                plt.show()
            else:
                print("report saved!")
                plt.close("all")
                sc.pl.rank_genes_groups(clust_adata, groups=fariq.get().strip().split(","), n_genes=20,show=False)
                plt.savefig(f"{place}/{output.get().strip()}/rank_genes_groups_[{'-'.join(fariq.get().strip().split(','))}]VS[{reffe.get().strip()}].png")

            os.system(f"mv {place}/{output.get().strip()}/* $HOME/Desktop/{output.get().strip()} && rmdir") # à la fin ge3 je déplace entiérement le dossier dans le bureau.

# </editor-fold>