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
parser.add_argument("--genepos", required = False, help = "Returns the name of the gene.")
parser.add_argument("--allsamps", required = False, help = "Do you want too see the info of all the samples, or just the expressed ones?.")
parser.add_argument("--ordered", required = False, help = "Do you want the barplot ordered?.")
parser.add_argument("--outfile", required = False, help = "The name of the bigsummary file.")
parser.add_argument("--minexp", required = False, help = "A threshold used to put the least expressed isoform as 'other'.")
parser.add_argument("--mintotexp", required = False, help = "A second threshold, in this case, to select the most significative isoforms.")
parser.add_argument("--minsamps", required = False, help = "A third threshold, in this case, to determine the minimum of samples to take a gene into account.")

parser.add_argument("-GI", "--geneinfo", nargs = "*", help = "Returns the name of the gene.")
parser.add_argument("-GN", "--genenameF", nargs = "*", help = "Returns the name of the gene.")
parser.add_argument("-GS", "--genestats", nargs = "*", help = "Returns the mean proportion of each isoform of the gene.")
parser.add_argument("-GFP", "--genefilprop", nargs = "*", help = "Returns a list of the isoforms and their respective total proportions, reversed ordered by these proportions.")
parser.add_argument("-GC", "--geneclass", nargs = "*", help = "Returns a list with different proportions for each isoform.")

#parser.add_argument("-M", "--matrix", nargs = "*", help = "Returns a plot showing all the isoforms expressed in all the samples of the gene you select, in a matrix way.")
#parser.add_argument("-Box", "--boxplot", nargs = "*", help = "Returns a box plot showing the first statistics of the isoforms of the gene you select.")
#parser.add_argument("-G", "--general", nargs = "*", help = "Returns a general view of the gene you select, using an stacked barplot.")
#parser.add_argument("-Bar", "--barplot", nargs = "*", help = "Returns a stacked bar plot showing just the most expressed isoforms of the gene you select.")

args = parser.parse_args()

data = args.data
genename = args.genename
genepos = args.genepos
allsamps = args.allsamps
ordered = args.ordered
outfile = args.outfile
minexp = args.minexp
mintotexp = args.mintotexp
minsamps = args.minsamps

geneinfo = args.geneinfo
genenameF = args.genenameF
genestats = args.genestats
genefilprop = args.genefilprop
geneclass = args.geneclass

#matrix = args.matrix
#boxplot = args.boxplot
#general = args.general
#barplot = args.barplot

if genename == None and genepos == None:
    print('ERROR: You must specify the gene you want to analyse.')
elif genename == None:
    gene = int(genepos)
else:
    gene = genename

if allsamps == None:
    allsamps = 'F'

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

if mintotexp != None:
    mintotexp = float(mintotexp)
else:
    mintotexp = 0.7

if minsamps != None:
    minsamps = int(minsamps)
else:
    minsamps = 10

DATA = Functions.reading_data(data)

pd.set_option('max_columns', 100)

if geneinfo != None:
    a = Functions.gene_info(DATA, gene, allsamps)
    
elif genenameF != None:
    a = Functions.gene_name(DATA, gene)
    
elif genestats != None:
    a = Functions.gene_statistics(DATA, gene, 'T')
    
elif genefilprop != None:
    a = Functions.gene_filtered_proportions(DATA, gene, 'T', minexp)
    
elif geneclass != None:
    a = Functions.gene_classification(DATA, gene, 'T', minexp, mintotexp, minsamps)
    
#elif matrix != None:
#    a = Functions.visual_matrix(DATA, gene)
#
#elif boxplot != None:
#    a = Functions.gene_boxplot(DATA, gene, ordered)
#
#elif general != None:
#    a = Functions.gene_general_view(DATA, gene)
#    
#elif barplot != None:
#    a = Functions.stacked_barplot(DATA, gene, minexp)

if type(a) == str:
    print(a)
else:
    pprint.pprint(a)