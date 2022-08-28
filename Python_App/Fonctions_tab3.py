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
from tkinter.filedialog import askopenfilename
from tkinter.font import BOLD, Font
from tkinter import ttk, simpledialog
from itertools import chain
import gseapy
from gseapy.plot import barplot, dotplot

endroit = os.getcwd()



def txt_summary(togo, organism):
    gene_set_names = gseapy.get_library_name(organism=organism)
    txt = ""
    for x in gene_set_names:
        txt = txt + x + "\n"

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
    # text_widget["state"] = DISABLED


# <editor-fold desc="ORA Method">
def ORA(root):

    def KEGG():
        lala.grid_forget()
        BP.grid_forget()
        MF.grid_forget()
        CC.grid_forget()


    def GOO():
        lala.grid(column=0, row=6, sticky="w")
        BP.grid(column=0, row=6)
        MF.grid(column=0, row=6, sticky="e")
        CC.grid(column=1, row=6, sticky="w")

    target = tk.StringVar()

    background = tk.StringVar()
    background.set("Optional")  # the function enrichGO/KEGG provide a default background to do this analysis.
    taga = Toplevel(root)
    target.set("Gene Symbol List")
    taga.title("OR Analysis")
    taga.geometry('+%d+%d' % (800, 750))

    saheh = tk.StringVar()

    pivalue = tk.DoubleVar()

    pivalue.set(0.05)

    GO=tk.StringVar()
    anto = tk.StringVar()

    orga = tk.StringVar()

    anto.set("MF")

    Label(taga, text="Over-Representation Analysis", font=Font(root, size=13, weight=BOLD)).grid(column=0, row=0,
                                                                                                 columnspan=3)
    Radiobutton(taga, text="GO Mode", command=lambda: GOO(),variable=GO,value="GO").grid(column=0, row=1)
    Radiobutton(taga, text="KEGG Mode", command=lambda: KEGG(),variable=GO,value="KEGG").grid(column=1, row=1, sticky="w")
    Label(taga, text="Target Gene list separated by <ENTER>").grid(column=0, row=2)
    Entry(taga, textvariable=target, width=15).grid(column=1, row=2, sticky="E", padx=10)
    Label(taga, text="Background Gene list separated by <ENTER>").grid(column=0, row=3)
    Entry(taga, textvariable=background, width=15).grid(column=1, row=3, sticky="E", padx=10)
    Label(taga, text="P-value Cutoff").grid(column=0, row=4)
    Entry(taga, textvariable=pivalue, width=15).grid(column=1, row=4, sticky="E", padx=10)
    Label(taga, text="Correction Method :").grid(column=0, row=5, sticky="w")
    Radiobutton(taga, text="Benjamini-Hocheberg", value="BH", variable=saheh).grid(column=0, row=5, sticky="E")
    Radiobutton(taga, text="Bonferroni", value="bonferroni", variable=saheh).grid(column=1, row=5, sticky="w")
    lala = ttk.Label(taga, text="Ontology:")
    lala.grid(column=0, row=6, sticky="w")
    BP = ttk.Radiobutton(taga, text="BP", value="BP", variable=anto)
    BP.grid(column=0, row=6)
    MF = ttk.Radiobutton(taga, text="MF", value="MF", variable=anto)
    MF.grid(column=0, row=6, sticky="e")
    CC = ttk.Radiobutton(taga, text="CC", value="CC", variable=anto)
    CC.grid(column=1, row=6, sticky="w")
    Label(taga, text="Organism:").grid(column=0, row=7, sticky="w")
    Radiobutton(taga, text="Homo Sapiens", value="human", variable=orga).grid(column=0, row=7)
    Radiobutton(taga, text="Mus Musculus", value="mouse", variable=orga).grid(column=1, row=7, sticky="w")
    # j'ai enlevé le save report car s'il choisit hide report ça va rien garder/montrer donc ça ne sert à rien.

    Button(taga, text="Enrich!", command=lambda: ORA_analysis(root, target, background, saheh, pivalue, anto, orga,GO),
           font=Font(root, size=10, weight=BOLD)).grid(column=0, row=9, columnspan=3)


def ORA_analysis(root, cible, contexte, korrect, pval, onto, organizem,GO_flag):
    flager = 0
    if cible.get == "" or contexte.get() == "" or korrect.get() == "" or pval.get() == "" or organizem.get() == "" or GO_flag.get()=="":
        messagebox.showerror("Error Params", "Please Fill all the \nparameters before continuing")
    else:
        with open("{0}/target".format(os.getcwd()), mode="w") as file:
            file.write(cible.get())

        if contexte.get().strip() != "Optional":  # meaning that the user has entered a background list (we don't care which length)
            with open("{0}/background".format(os.getcwd()), mode="w") as file:
                file.write(contexte.get())
                flag_contexte = 1
        else:
            with open("{0}/background".format(os.getcwd()), mode="w") as file:
                file.write("nothing")  # i just create the file but it contains nothing
            flag_contexte = 0  # i will have to use the background provided by the methods enrich !

        os.system("rm Enrich_result_GO.csv") #remove the previous one in case the enrichment is empty
        os.system("rm Enrich_result_KEGG.csv")

        place = simpledialog.askstring("Output dir","The results are being generated\nplease provide a dir_name to store the output")  # we are in a docker don't ask for the location

        os.system(f"mkdir {place.strip()}") #create the directory

        os.system("Rscript --vanilla {0} {1} {2} {3} {4} {5} {6} {7} {8} {9}".format(
            os.getcwd() + "/GO_lab.R", os.getcwd(), GO_flag.get(), os.getcwd() + "/target", flag_contexte,
            os.getcwd() + "/background", organizem.get(), pval.get(), korrect.get(), onto.get()))


        print("Flag=", GO_flag.get())

        if GO_flag.get() == "GO":
            try:
                data = pd.read_csv("Enrich_result_GO.csv")
                data.to_csv(f"{place.strip()}/Enrich_result_GO.csv") #re write the file in the good output_dir
            except FileNotFoundError:
                messagebox.showerror("Error", "Something went wrong the csv is not found\nMaybe 0 enrichment terms were found")
                flager = 1
        else:
            try:
                data = pd.read_csv("Enrich_result_KEGG.csv")
                data.to_csv(f"{place.strip()}/Enrich_result_KEGG.csv") #re write the file in the good output_dir
            except FileNotFoundError:
                messagebox.showerror("Error", "Something went wrong the csv is not found\nMaybe 0 enrichment terms were found")
                flager = 1

        if flager == 0:  # no problem carry on the process
            data = data[['Description', 'GeneRatio', 'p.adjust', 'qvalue']]
            data.columns = ['Term', 'Overlap', 'Adjusted P-value', 'Combined Score']
            data['Combined Score'] = 1 / data['Combined Score']  # to have correct values

            gseapy.barplot(data, ofname=f"{place.strip()}/result_bar.png", top_term=15)
            gseapy.dotplot(data, ofname=f"{place.strip()}/result_dot.png", top_term=15)

            wind = Toplevel(root)
            wind.title("End report")
            if GO_flag.get() == "GO":
                txt = f"The results has been generated:\n{endroit}/{place.strip()}/result_dot.png\n{endroit}/{place.strip()}/result_bar.png\n{endroit}/{place.strip()}/Enrich_result_GO.csv"
            else:
                txt = f"The results has been generated:\n{endroit}/{place.strip()}/result_dot.png\n{endroit}/{place.strip()}/result_bar.png\n{endroit}/{place.strip()}/Enrich_result_KEGG.csv"

            Label(wind, text=txt, font=Font(weight=BOLD)).pack()
            messagebox.showinfo("ORA Info", "The Enrichment Analysis\nis done !!")


# </editor-fold>


# <editor-fold desc="GSEA Method">
def GSEA(root):
    def open_file_chooser():
        global filename, flag_filename
        filename = askopenfilename()
        flag_filename = messagebox.askyesno("Format Question", "Does your file contain a header ?")
        print(filename)

    target = tk.StringVar()
    logf = tk.StringVar()
    taga = Toplevel(root)
    taga.title("GSE Analysis")
    taga.geometry('+%d+%d' % (800, 750))

    gene7 = tk.StringVar()
    var_organ = tk.StringVar()

    porte = tk.StringVar()
    perma = tk.StringVar()
    perma.set("100")

    target.set("BTNL9\nRNASE1\nITGA1\nHLA-DRB1\nPTPRB\nEGFL7")
    logf.set("8.73524761199951\n7.2530255317688\n7.0928897857666\n6.97891712188721\n6.95244693756104\n6.95154476165772")

    Label(taga, text="Gene Set enrichment analysis\n(GSEA)", font=Font(root, size=13, weight=BOLD)).grid(column=0,
                                                                                                         row=0,
                                                                                                         columnspan=3)
    ttk.Radiobutton(taga, text="Enter gene list in text bar", value="entry", variable=porte).grid(column=0, row=1,sticky="w")
    Label(taga, text="Ordered Gene list separated by <ENTER>").grid(column=0, row=2)
    Entry(taga, textvariable=target, width=15).grid(column=1, row=2, sticky="E", padx=10)
    Label(taga, text=" LogFC for each Gene separated by <ENTER>").grid(column=0, row=3)
    Entry(taga, textvariable=logf, width=15).grid(column=1, row=3, sticky="E", padx=10)
    ttk.Radiobutton(taga, text="Enter gene list in a file", value="csv_txt", variable=porte).grid(column=0, row=4,
                                                                                                  sticky="w")
    Label(taga, text=" Enter ordered gene list with a CSV file").grid(column=0, row=5)
    Button(taga, text="Browse", command=lambda: open_file_chooser()).grid(column=1, row=5, sticky="W", padx=10)
    Label(taga, text="Which Organism :").grid(column=0, row=6, sticky="w")
    ttk.Radiobutton(taga, text="Human", value="Human", variable=var_organ).grid(column=0, row=6)
    ttk.Radiobutton(taga, text="Mouse", value="Mouse", variable=var_organ).grid(column=0, row=6, sticky="E")
    ttk.Radiobutton(taga, text="Yeast", value="Yeast", variable=var_organ).grid(column=1, row=6, sticky="w")
    ttk.Radiobutton(taga, text="Fly", value="Fly", variable=var_organ).grid(column=1, row=6, sticky="e")
    Label(taga, text="Which Gene Set").grid(column=0, row=7, sticky="w")
    gene7.set("GO_Biological_Process_2021")
    Entry(taga, textvariable=gene7, width=10).grid(column=0, row=7)
    Button(taga, text="List", command=lambda: txt_summary(taga, organism=var_organ.get())).grid(column=0, row=7,sticky="e")
    Label(taga, text="Permutation").grid(column=1, row=7, sticky="w")
    Entry(taga, textvariable=perma, width=7).grid(column=1, row=7, sticky="e")

    def gate():
        if porte.get() == "entry":
            GSEA_Analysis(root, target, logf, "None", porte, var_organ, gene7, "None", perma)
        elif porte.get() == "csv_txt":
            GSEA_Analysis(root, target, logf, filename, porte, var_organ, gene7, flag_filename, perma)
        elif porte.get() == "" or perma.get() == "":
            messagebox.showerror("Error Input", "Please fill all the parameters before continuing")

    Button(taga, text="Enrich!", font=Font(root, size=10, weight=BOLD), command=lambda: gate()).grid(column=0, row=10,columnspan=3)


def GSEA_Analysis(root, cible, foldchange, filou, flag, organizem, gene_set, header, permutation):
    flager = 0
    if gene_set.get() == "" or organizem.get() == "":
        messagebox.showerror("Error Input", "Please fill all the parameters before continuing")
    else:
        if flag.get() == "csv_txt":
            print(filou)
            if filou == "":
                messagebox.showerror("Input error", "You didn't provide any file!!")
            elif filou.endswith(".csv") == False:
                messagebox.showerror("Format Error", "The file you provided is not in format csv\nProvide a .csv file")

            else:
                if header:
                    table = pd.read_csv(filou)
                    table.sort_values(by=table.columns[1], inplace=True, ascending=False)
                else:
                    table = pd.read_csv(filou, header=None)
                    table.sort_values(by=[1], inplace=True, ascending=False)

                place = simpledialog.askstring("Output dir",
                                               "The results are being generated\nplease provide a dir_name to store the output")
                try:
                    gsea_res = gseapy.prerank(rnk=table, gene_sets=gene_set.get().strip(),
                                              permutation_num=int(permutation.get()),
                                              outdir="{0}/{1}".format(endroit, place.strip()))
                except ValueError:
                    messagebox.showerror("Error",
                                         "Please check your csv file\nYou should have only 2 columns:\ngene name in first column\nlogf in second column")
                    flager = 1
                except TypeError:
                    messagebox.showerror("Error","Please check your csv file\n verify the format(str and float) of the 2 columns")
                    flager = 1

                if flager == 0:
                    terms = gsea_res.res2d.index
                    gseapy.gseaplot(rank_metric=gsea_res.ranking, term=terms[0], **gsea_res.results[terms[0]])

                    wind = Toplevel(root)
                    wind.title("End report")
                    txt = "The results has been generated and can be found in:\n{0}/{1}".format(endroit, place.strip())
                    Label(wind, text=txt, font=Font(weight=BOLD)).pack()
                    messagebox.showinfo("GSEA Info", "The Enrichment Analysis\nis done !!")

        elif flag.get() == "entry":
            if cible.get() == "" or foldchange.get() == "":
                messagebox.showerror("Input error", "You have to fill the 2 first entries\nbefore continuing")
            elif not len(cible.get().strip().split("\n")) == len(foldchange.get().strip().split("\n")):
                messagebox.showerror("Input Error", "The 2 lists do not have the same length")
            else:
                try:
                    float_logf = [float(item) for item in foldchange.get().strip().split(
                        "\n")]  # the value are in type string so i convert them into float
                    zipped = list(zip(cible.get().strip().split("\n"), float_logf))
                    table = pd.DataFrame(zipped, columns=['genes', 'LogFC'])
                    table.sort_values(by=['LogFC'], inplace=True,ascending=False)  # if the user give a ordered list or no , i will order it and then process it
                    place = simpledialog.askstring("Output dir","The results are being generated\nplease provide a directory name to store the output")
                except ValueError:
                    messagebox.showerror("Error",
                                         "Please check your input\n the first entry contain only gene name \n the second entry , only floats")

            try:
                gsea_res = gseapy.prerank(rnk=table, gene_sets=gene_set.get().strip(),
                                          permutation_num=int(permutation.get()),
                                          outdir="{0}/{1}".format(endroit, place.strip()))
            except TypeError:
                messagebox.showerror("Error",
                                     "Please check your input\n the first entry contain only gene name \n the second entry , only floats")
                flager = 1

            if flager == 0:
                terms = gsea_res.res2d.index
                gseapy.gseaplot(rank_metric=gsea_res.ranking, term=terms[0], **gsea_res.results[terms[0]])
                wind = Toplevel(root)
                wind.title("End report")
                txt = "The results has been generated and can be found in:\n{0}/{1}".format(endroit, place.strip())
                Label(wind, text=txt, font=Font(weight=BOLD)).pack()
                messagebox.showinfo("GSEA Info", "The Enrichment Analysis\nis done !!")


# </editor-fold>


# <editor-fold desc="Enrichr Method">
def Enrichr(root):
    target = tk.StringVar()
    background = tk.StringVar()
    background.set("Optional")
    taga = Toplevel(root)
    taga.title("Enrichr Analysis")
    taga.geometry('+%d+%d' % (800, 750))

    gene7 = tk.StringVar()
    var_organ = tk.StringVar()
    pivalue = tk.StringVar()
    pivalue.set("0.05")

    Label(taga, text="Enrichr analysis", font=Font(root, size=13, weight=BOLD)).grid(column=0, row=0, columnspan=3)
    Label(taga, text="Gene list separated by <ENTER>").grid(column=0, row=2)
    Entry(taga, textvariable=target, width=15).grid(column=1, row=2, sticky="E", padx=10)
    Label(taga, text=" Background gene list separated by <ENTER>").grid(column=0, row=3)
    Entry(taga, textvariable=background, width=15).grid(column=1, row=3, sticky="E", padx=10)
    Label(taga, text="P-value Cutoff").grid(column=0, row=4)
    Entry(taga, textvariable=pivalue, width=15).grid(column=1, row=4, sticky="E", padx=10)
    Label(taga, text="Which Organism :").grid(column=0, row=6, sticky="w")
    ttk.Radiobutton(taga, text="Human", value="Human", variable=var_organ).grid(column=0, row=6)
    ttk.Radiobutton(taga, text="Mouse", value="Mouse", variable=var_organ).grid(column=0, row=6, sticky="E")
    ttk.Radiobutton(taga, text="Yeast", value="Yeast", variable=var_organ).grid(column=1, row=6, sticky="w")
    ttk.Radiobutton(taga, text="Fly", value="Fly", variable=var_organ).grid(column=1, row=6, sticky="e")
    Label(taga, text="Which Gene Set").grid(column=0, row=7, sticky="w", pady=10)
    Entry(taga, textvariable=gene7, width=10).grid(column=0, row=7)
    Button(taga, text="List", command=lambda: txt_summary(taga, organism=var_organ.get())).grid(column=0, row=7,
                                                                                                sticky="e")
    Button(taga, text="Enrich!", font=Font(root, size=10, weight=BOLD),
           command=lambda: Enrichr_Analysis(root, target, background, var_organ, gene7, pivalue)).grid(column=0, row=10,
                                                                                                       columnspan=3)


def Enrichr_Analysis(root, cible, contexte, organism, gene_set, pval):
    if cible.get == "" or contexte.get() == "" or pval.get() == "" or organism.get() == "":
        messagebox.showerror("Error Params", "Please Fill all the \nparameters before continuing")
    elif cible.get == "" and (contexte.get() == "" or contexte.get() == "Optional"):
        print("rak hna ?")
        messagebox.showerror("Error Params", "Please Fill all the \nparameters before continuing")
    else:
        list_target = cible.get().strip().split("\n")
        if contexte != "Optional":
            list_background = contexte.get().strip().split("\n")
            flag = 1
        else:
            flag = 0

        if flag == 1:
            place = simpledialog.askstring("Output dir",
                                           "The results are being generated\nplease provide a directory name to store the output")
            enr = gseapy.enrichr(gene_list=list_target,
                                 gene_sets=gene_set.get().strip(),
                                 background=list_background,
                                 organism=organism.get().strip(),
                                 # don't forget to set organism to the one you desired! e.g. Yeast
                                 outdir='{0}/{1}'.format(endroit, place),
                                 no_plot=False,
                                 cutoff=float(pval.get())  # test dataset, use lower value from range(0,1)
                                 )
            wind = Toplevel(root)
            wind.title("End report")
            txt = "The results has been generated and can be found in:\n{0}/{1}".format(endroit, place.strip())
            Label(wind, text=txt, font=Font(weight=BOLD)).pack()
            messagebox.showinfo("Enrichr Info", "The Enrichment Analysis\nis done !!")


        elif flag == 0:
            place = simpledialog.askstring("Output dir",
                                           "The results are being generated\nplease provide a directory name to store the output")
            enr = gseapy.enrichr(gene_list=list_target,
                                 gene_sets=gene_set.get().strip(),
                                 organism=organism.get().strip(),
                                 outdir='{0}/{1}'.format(endroit, place),
                                 no_plot=False,
                                 cutoff=float(pval.get())  # test dataset, use lower value from range(0,1)
                                 )
            wind = Toplevel(root)
            wind.title("End report")
            txt = "The results has been generated and can be found in:\n{0}/{1}".format(endroit, place.strip())
            Label(wind, text=txt, font=Font(weight=BOLD)).pack()
            messagebox.showinfo("Enrichr Info", "The Enrichment Analysis\nis done !!")
# </editor-fold>