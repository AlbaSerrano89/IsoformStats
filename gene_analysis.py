#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 29 17:01:31 2018

@author: aserrano
"""
import argparse
import os
import sys
import pprint
import Functions
import pandas as pd

parser = argparse.ArgumentParser()

parser.add_argument("--data", required = True, help = "A csv or gzip file with the data you want to analyse.")
parser.add_argument("--genename", required = False, type = str, help = "The gene to analyse.")
parser.add_argument("--genepos", required = False, type = int, help = "The position of the gene to analyse.")
parser.add_argument("--minexp", required = False, type = float, default = 0.8, help = "A threshold used to put the least expressed isoform as 'other'.")
parser.add_argument("--minsamps", required = False, type = int, default = 10, help = "A threshold to determine the minimum of samples to take a gene into account.")
parser.add_argument("--plotfile", required = False, help = "The name of the plot file.")

parser.add_argument("-GN", "--genenameF", nargs = "*", help = "Returns the name of the gene.")
parser.add_argument("-GP", "--geneposF", nargs = "*", help = "Returns the position of the gene in the list.")
parser.add_argument("-GI", "--geneinfo", nargs = "*", help = "Returns all the information contained in 'data' file about the gene.")
parser.add_argument("-GS", "--genestats", nargs = "*", help = "Returns some statistics of each isoform of the gene, reversed ordered by the mean of the expression of each one.")
parser.add_argument("-GFP", "--genefilprop", nargs = "*", help = "Returns the same as GS, but filtering by the minimum expression you want to be reached.")
parser.add_argument("-GC", "--geneclass", nargs = "*", help = "Returns a list with different proportions for each isoform.")
parser.add_argument("-GBox", "--geneboxplot", nargs = "*", help = "Returns a box plot showing the first statistics of the isoforms of the gene you select.")
parser.add_argument("-GM", "--genematrix", nargs = "*", help = "Returns a plot showing all the isoforms expressed in all the samples of the gene you select, in a matrix way.")
parser.add_argument("-GBar", "--genebarplot", nargs = "*", help = "Returns a general view of the gene you select, using an stacked barplot.")
parser.add_argument("-GFB", "--genefiltbarplot", nargs = "*", help = "Returns a stacked bar plot showing just the most expressed isoforms of the gene you select.")

args = parser.parse_args()

data = args.data
genename = args.genename
genepos = args.genepos
minexp = args.minexp
minsamps = args.minsamps
plotfile = args.plotfile

genenameF = args.genenameF
geneposF = args.geneposF
geneinfo = args.geneinfo
genestats = args.genestats
genefilprop = args.genefilprop
geneclass = args.geneclass
geneboxplot = args.geneboxplot
genematrix = args.genematrix
genefiltbarplot = args.genefiltbarplot
genebarplot = args.genebarplot

# ---------------------- Starting parameters ---------------------- #
if genename == None and genepos == None:
    sys.exit('ERROR: You must specify the gene you are interested in.')
elif genename == None:
    gene = genepos
else:
    gene = genename

# ---------------------- Control ---------------------- #
if genenameF == None and geneposF == None and geneinfo == None and genestats == None and genefilprop == None and geneclass == None and geneboxplot == None and genematrix == None and genebarplot == None and genefiltbarplot == None:
    sys.exit("ERROR: You did not write any function.")

# ---------------------- Functions ---------------------- #
DATA = Functions.reading_data(data)

plotdir = 'plots'
pd.set_option('max_columns', 100)

if genenameF != None:
    a = Functions.gene_name(DATA, gene)
    
elif geneposF != None:
    a = Functions.gene_pos(DATA, gene)

elif geneinfo != None:
    a = Functions.gene_info(DATA, gene, False)
    
elif genestats != None:
    a = Functions.gene_statistics(DATA, gene, False)
    
elif genefilprop != None:
    a = Functions.gene_filtered_proportions(DATA, gene, False, minexp)
    
elif geneclass != None:
    a = Functions.gene_classification(DATA, gene, False, minexp, minsamps)
    
elif genematrix != None:
    if plotfile == None:
        genename = Functions.gene_name(DATA, gene)
        
        if not os.path.isdir(plotdir):
            os.mkdir(plotdir)
        
        plotfile = '{}/{}_Matrix.png'.format(plotdir, genename)
    
    Functions.gene_matrix(DATA, gene, plotfile)

elif genebarplot != None:
    if plotfile == None:
        genename = Functions.gene_name(DATA, gene)
        
        if not os.path.isdir(plotdir):
            os.mkdir(plotdir)
        
        plotfile = '{}/{}_Barplot.png'.format(plotdir, genename)
    
    Functions.gene_barplot(DATA, gene, plotfile)

elif geneboxplot != None:
    if plotfile == None:
        genename = Functions.gene_name(DATA, gene)
        
        if not os.path.isdir(plotdir):
            os.mkdir(plotdir)
        
        plotfile = '{}/{}_Boxplot.png'.format(plotdir, genename)
    
    Functions.gene_boxplot(DATA, gene, plotfile)

elif genefiltbarplot != None:
    if plotfile == None:
        genename = Functions.gene_name(DATA, gene)
        
        if not os.path.isdir(plotdir):
            os.mkdir(plotdir)
        
        plotfile = '{}/{}_{}_FilteredBarplot.png'.format(plotdir, genename, minexp)
    
    Functions.gene_filtered_barplot(DATA, gene, plotfile, minexp)

try:
    if type(a) == str:
        print(a)
    else:
        pprint.pprint(a)
        print('\n')
except NameError:
    pass
