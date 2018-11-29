#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 29 17:01:31 2018

@author: aserrano
"""

import argparse
import pprint
import Functions
import pandas as pd

parser = argparse.ArgumentParser()

parser.add_argument("--data", required = True, help = "A csv or gzip file with the data you want to analyse.")
parser.add_argument("--genename", required = False, help = "The gene to analyse.")
parser.add_argument("--genepos", required = False, help = "The position of the gene to analyse.")
parser.add_argument("--dftype", required = False, help = "Type of answer: dictionary or dataframe?")
parser.add_argument("--ordered", required = False, help = "Do you want the barplot ordered?")
parser.add_argument("--minexp", required = False, help = "A threshold used to put the least expressed isoform as 'other'.")
parser.add_argument("--minsamps", required = False, help = "A threshold to determine the minimum of samples to take a gene into account.")

parser.add_argument("-GN", "--genenameF", nargs = "*", help = "Returns the name of the gene.")
parser.add_argument("-GP", "--geneposF", nargs = "*", help = "Returns the name of the gene.")
parser.add_argument("-GI", "--geneinfo", nargs = "*", help = "Returns all the information contained in 'data' file about the gene.")
parser.add_argument("-GS", "--genestats", nargs = "*", help = "Returns the mean proportion of each isoform of the gene.")
parser.add_argument("-GFP", "--genefilprop", nargs = "*", help = "Returns a list of the isoforms and their respective total proportions, reversed ordered by these proportions.")
parser.add_argument("-GC", "--geneclass", nargs = "*", help = "Returns a list with different proportions for each isoform.")

parser.add_argument("-GBox", "--geneboxplot", nargs = "*", help = "Returns a box plot showing the first statistics of the isoforms of the gene you select.")
parser.add_argument("-GM", "--genematrix", nargs = "*", help = "Returns a plot showing all the isoforms expressed in all the samples of the gene you select, in a matrix way.")
parser.add_argument("-GBar", "--genebarplot", nargs = "*", help = "Returns a general view of the gene you select, using an stacked barplot.")
parser.add_argument("-GFB", "--genefiltbarplot", nargs = "*", help = "Returns a stacked bar plot showing just the most expressed isoforms of the gene you select.")

args = parser.parse_args()

data = args.data
genename = args.genename
genepos = args.genepos
dftype = args.dftype
ordered = args.ordered
minexp = args.minexp
minsamps = args.minsamps

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

if genename == None and genepos == None:
    print('ERROR: You must specify the gene you are interested in.')
elif genename == None:
    gene = int(genepos)
else:
    gene = genename

if dftype == None:
    dftype = 'T'

if ordered != None:
    if ordered == ('False' or 'F'):
        ordered = False
    elif ordered == ('True' or 'T'):
        ordered = True
else:
    ordered = True

if minexp != None:
    minexp = float(minexp)
else:
    minexp = 0.1

if minsamps != None:
    minsamps = int(minsamps)
else:
    minsamps = 10

DATA = Functions.reading_data(data)

pd.set_option('max_columns', 100)

if genenameF != None:
    a = Functions.gene_name(DATA, gene)
    
elif geneposF != None:
    a = Functions.gene_pos(DATA, gene)

elif geneinfo != None:
    a = Functions.gene_info(DATA, gene, dftype)
    
elif genestats != None:
    a = Functions.gene_statistics(DATA, gene, dftype)
    
elif genefilprop != None:
    a = Functions.gene_filtered_proportions(DATA, gene, dftype, minexp)
    
elif geneclass != None:
    a = Functions.gene_classification(DATA, gene, dftype, minexp, minsamps)
    
elif geneboxplot != None:
    Functions.gene_boxplot(DATA, gene, ordered)

elif genematrix != None:
    Functions.gene_matrix(DATA, gene)

elif genebarplot != None:
    Functions.gene_barplot(DATA, gene)
    
elif genefiltbarplot != None:
    Functions.gene_filtered_barplot(DATA, gene, minexp)

try:
    if type(a) == str:
        print(a)
    else:
        pprint.pprint(a)
except NameError:
    pass