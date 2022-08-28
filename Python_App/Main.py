import matplotlib
#matplotlib.use('Agg') #Anti-grain geometry engine. use for docker image
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
import Fonctions_tab1 as Ftab1
import Fonctions_tab2 as Ftab2
import Fonctions_tab3 as Ftab3
from tkinter.font import BOLD, Font
import matplotlib.pyplot as plt


# <editor-fold desc="Main Window">
splash_win= tk.Tk()

splash_win.geometry('+%d+%d' % (300,200))
txt="~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"

Label(splash_win, text= txt,font= ('Times New Roman', 30,"bold italic"),foreground="black").pack()
Label(splash_win, text="Welcome To",font= ('Times New Roman', 30,"bold italic"),foreground="black").pack()
Label(splash_win, text="Wahed_cell", font=('Times New Roman',40,"bold italic"),foreground="salmon").pack()
Label(splash_win, text= "A User-Friendly tool to Perform Single-cell RNA-seq Analysis",font= ('Times New Roman', 30,"bold italic"),foreground="black").pack()
Label(splash_win, text=txt,font= ('Times New Roman', 30,"bold italic"),foreground="black").pack()

placa="/home/izem/Pictures/"

image = Image.open(f"{placa}imrb.jpeg")
resize_img = image.resize((300, 100))
sola = Label(splash_win)
sola.image = ImageTk.PhotoImage(resize_img)
sola['image'] = sola.image
sola.place(x=0, y=40)

image_a = Image.open(f"{placa}upec.png")
resize_img_a = image_a.resize((200, 130))
solb = Label(splash_win)
solb.image = ImageTk.PhotoImage(resize_img_a)
solb['image'] = solb.image
solb.place(x=900, y=35)




splash_win.overrideredirect(True)

def mainWin():
   splash_win.destroy()
   root = tk.Tk()
   root.title("Wahed_Cell")

   root.geometry("690x165")

   tabControl = ttk.Notebook(root)
   tab1 = ttk.Frame(tabControl)
   tab2 = ttk.Frame(tabControl)
   tab3 = ttk.Frame(tabControl)

   tabControl.add(tab1, text='Quick Analysis')
   tabControl.add(tab2, text='Advanced Analysis')
   tabControl.add(tab3, text="Functional Enrichment")
   tabControl.pack(expand=1, fill="both")

   # <editor-fold desc="TAB ONE">


   Label(tab1, text="COUNT MATRIX").grid(column=0, row=0, sticky="NW")

   fullpath_count_tab1 = tk.StringVar()
   samplename_count_tab1 = tk.StringVar()

   Entry(tab1, textvariable=fullpath_count_tab1).grid(column=1, row=0, sticky="NW")
   Entry(tab1, textvariable=samplename_count_tab1).grid(column=1, row=1, sticky="NW")

   fullpath_count_tab1.set("enter fullpath")
   samplename_count_tab1.set("enter sample name")

   format_tab1 = tk.StringVar()
   format_tab1.set("10x_mtx")
   Label(tab1, text="File Format :").grid(column=0, row=2, sticky="W")
   Radiobutton(tab1, text="10x_mtx", value="10x_mtx", variable=format_tab1).grid(column=1, row=2, sticky="w")
   Radiobutton(tab1, text="mtx", value="mtx", variable=format_tab1).grid(column=1, row=2, sticky="e")
   Radiobutton(tab1, text="h5ad", value="h5ad", variable=format_tab1).grid(column=2, row=2, sticky="w")
   Radiobutton(tab1, text="loom", value="loom", variable=format_tab1).grid(column=3, row=2, sticky="w")

   radio_string = tk.StringVar()
   Label(tab1, text="Analysis Mode :").grid(column=0, row=3, sticky="W")
   ttk.Radiobutton(tab1, text="Uni-sample", value="un", variable=radio_string).grid(column=2, row=3,sticky="W")

   ttk.Radiobutton(tab1, text="Multi-sample", value="deux", variable=radio_string).grid(column=1,row=3,sticky="W")
   Label(tab1, text="Show Reports :").grid(column=0, row=4, sticky="W")
   radio_report = tk.StringVar()

   ttk.Button(tab1, text="OKEY",command=lambda: Ftab1.check_path(radio_string, samplename_count_tab1, fullpath_count_tab1, radio_report,format_tab1,output_tab1)).place(x=280, y=13)  # j'ai utiliser .place mais je vais peut-être le regretter

   ttk.Radiobutton(tab1, text="YES", value="oui", variable=radio_report).grid(column=2, row=4, sticky="W")
   ttk.Radiobutton(tab1, text="NO", value="non", variable=radio_report).grid(column=1, row=4, sticky="W")

   output_tab1=tk.StringVar()

   Label(tab1,text="Output Dir ").grid(column=3,row=4,sticky="w")
   Entry(tab1, textvariable=output_tab1,width=10).grid(column=4, row=4,sticky="w")

   ttk.Button(tab1, text="Launch Analysis", command=lambda: Ftab1.Analysis_Gate(root, radio_string, radio_report,output_tab1)).place(x=150, y=110)  # i should remove this button

   # </editor-fold>

   # <editor-fold desc="TAB TWO">
   Label(tab2, text="COUNT MATRIX").grid(column=0, row=0, sticky="W")

   fullpath_count_tab2 = tk.StringVar()
   samplename_count_tab2 = tk.StringVar()

   Entry(tab2, textvariable=fullpath_count_tab2).grid(column=1, row=0, sticky="NW")
   Entry(tab2, textvariable=samplename_count_tab2).grid(column=1, row=1, sticky="NW")

   fullpath_count_tab2.set("/home/izem/Desktop/STAGE/Izem_Single-cell/muscle_1_ctrl/outs/filtered_feature_bc_matrix/")
   samplename_count_tab2.set("izem")

   radio_string_tab2 = tk.StringVar()
   Label(tab2, text="Analysis Mode :").grid(column=0, row=2, sticky="W")
   ttk.Radiobutton(tab2, text="Uni-sample", value="un", variable=radio_string_tab2).grid(
       column=2, row=2, sticky="W", ipadx=0)
   ttk.Radiobutton(tab2, text="Multi-sample", value="deux", variable=radio_string_tab2).grid(
       column=1, row=2, sticky="W")

   format_tab2 = tk.StringVar()
   format_tab2.set("10x_mtx")

   ttk.Button(tab2, text="OKEY",command=lambda: Ftab2.check_path(radio_string_tab2, samplename_count_tab2,fullpath_count_tab2, format_tab2,output_tab2)).place(x=284,y=13)
   ttk.Button(tab2, text="QC Filtering", command=lambda: Ftab2.QC_FILTER(root, radio_string_tab2,format_tab2,output_tab2)).place(x=420, y=5)
   ttk.Button(tab2, text="Normalization", command=lambda: Ftab2.Normalizer(root, radio_string_tab2,output_tab2)).place(x=520, y=5)
   ttk.Button(tab2, text="Integration", command=lambda: Ftab2.Integrater(root, radio_string_tab2,output_tab2)).place(x=420, y=40)
   ttk.Button(tab2, text="Feature Select", command=lambda: Ftab2.Selecter(root, radio_string_tab2,output_tab2)).place(x=530, y=40)
   ttk.Button(tab2, text="Dim. reduc/Clustering", command=lambda: Ftab2.Reduc_clusterer(root,output_tab2)).place(x=465, y=75)
   ttk.Button(tab2, text="Marker Detection", command=lambda: Ftab2.Mark_Dectecter(root,output_tab2)).place(x=420, y=110)
   ttk.Button(tab2, text="Compare Cluster", command=lambda: Ftab2.Comparer(root,output_tab2)).place(x=550, y=110)


   Label(tab2, text="File Format :").grid(column=0, row=3, sticky="W")
   Radiobutton(tab2, text="10x_mtx", value="10x_mtx", variable=format_tab2).grid(column=1, row=3, sticky="w")
   Radiobutton(tab2, text="mtx", value="mtx", variable=format_tab2).grid(column=1, row=3, sticky="e")
   Radiobutton(tab2, text="h5ad", value="h5ad", variable=format_tab2).grid(column=1, row=4, sticky="E")
   Radiobutton(tab2, text="loom", value="loom", variable=format_tab2).grid(column=1, row=4, sticky="w")

   output_tab2 = tk.StringVar()

   Label(tab2, text="Output Directory ").grid(column=0, row=5, sticky="w")
   Entry(tab2, textvariable=output_tab2, width=10).grid(column=1, row=5, sticky="w")


   Button(tab2,fg="red", background="red",activebackground="red" ,width=1, height=200).place(x=378, y=0)

   # </editor-fold>

   # <editor-fold desc="TAB THREE">

   Label(tab3, text="Type of Enrichment Analysis", font=Font(root, size=13, weight=BOLD)).grid(column=0, row=0,columnspan=3)
   Button(tab3, text="Over-Representation Analysis\n(ORA)", command=lambda: Ftab3.ORA(root)).grid(column=0, row=1, sticky="W", pady=20)
   Button(tab3, text="Gene Set enrichment analysis\n(GSEA)", command=lambda: Ftab3.GSEA(root)).grid(column=1, row=1,sticky="E", pady=10, padx=10)
   Button(tab3, text="Enrichr\nMethod", command=lambda: Ftab3.Enrichr(root)).grid(column=2, row=1, pady=10, padx=10)

   # </editor-fold>



splash_win.after(3000, mainWin)

mainloop()
# </editor-fold>

