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
parser.add_argument("--thres1", required = False, help = "A threshold used to put the least expressed isoform as 'other'.")
parser.add_argument("--thres2", required = False, help = "A second threshold, in this case, to select the most significative isoforms.")
parser.add_argument("--thres3", required = False, help = "A third threshold, in this case, to determine the minimum of samples to take a gene into account.")

parser.add_argument("-GI", "--geneinfo", nargs = "*", help = "Returns the name of the gene.")
parser.add_argument("-GN", "--genenameF", nargs = "*", help = "Returns the name of the gene.")
parser.add_argument("-GP", "--geneprop", nargs = "*", help = "Returns the mean proportion of each isoform of the gene.")
parser.add_argument("-GFP", "--genefilprop", nargs = "*", help = "Returns a list of the isoforms and their respective total proportions, reversed ordered by these proportions.")
parser.add_argument("-GC", "--geneclass", nargs = "*", help = "Returns a list with different proportions for each isoform.")

parser.add_argument("-G", "--general", nargs = "*", help = "Returns a general view of the gene you select, using an stacked barplot.")
parser.add_argument("-B", "--stackedbarplot", nargs = "*", help = "Returns a stacked bar plot showing just the most expressed isoforms of the gene you select.")

args = parser.parse_args()

data = args.data
genename = args.genename
genepos = args.genepos
allsamps = args.allsamps
ordered = args.ordered
outfile = args.outfile
thres1 = args.thres1
thres2 = args.thres2
thres3 = args.thres3

geneinfo = args.geneinfo
genenameF = args.genenameF
geneprop = args.geneprop
genefilprop = args.genefilprop
geneclass = args.geneclass

general = args.general
stackedbarplot = args.stackedbarplot

if genename == None and genepos == None:
    print('ERROR: You must specify the gene you want to analyse.')
elif genename == None:
    gene = genepos
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

if thres1 != None:
    thres1 = float(thres1)
else:
    thres1 = 0.1

if thres2 != None:
    thres2 = float(thres2)
else:
    thres2 = 0.7

if thres3 != None:
    thres3 = int(thres3)
else:
    thres3 = 10

DATA = Functions.reading_data(data)

pd.set_option('max_columns', 100)

if geneinfo != None:
    try:
        gene = int(gene)
    except ValueError:
        gene = gene
    a = Functions.gene_info(DATA, gene, allsamps)
    
    if 'is not in your data' in a:
        print('ERROR: ' + a)
    else:
        pprint.pprint(a)

elif genenameF != None:
    try:
        gene = int(gene)
    except ValueError:
        gene = gene

    a = Functions.formatting_gene(DATA, gene)
    
    if 'is not in your data' in a:
        print('ERROR: ' + a)
    else:
        print(a)

elif geneprop != None:
    try:
        gene = int(gene)
    except ValueError:
        gene = gene

    a = Functions.gene_proportions(DATA, gene)
    
    if 'not expressed' in a or 'is not in your data' in a:
        print(a)
    else:
        pprint.pprint(a)

elif genefilprop != None:
    try:
        gene = int(gene)
    except ValueError:
        gene = gene

    a = Functions.gene_filtered_proportions(DATA, gene, thres1)
    
    if 'not expressed' in a or 'is not in your data' in a:
        print(a)
    else:
        pprint.pprint(a)

elif geneclass != None:
    try:
        gene = int(gene)
    except ValueError:
        gene = gene

    a = Functions.gene_classification(DATA, gene, thres1, thres2, thres3)
    
    if 'not expressed' in a or 'is not in your data' in a:
        print(a)
    else:
        pprint.pprint(a)
    
elif general != None:
    try:
        gene = int(gene)
    except ValueError:
        gene = gene

    a = Functions.gene_general_view(DATA, gene)
    
    if type(a) == str:
        print(a)
    
elif stackedbarplot != None:
    try:
        gene = int(gene)
    except ValueError:
        gene = gene

    a = Functions.stacked_barplot(DATA, gene, thres1)

    if type(a) == str:
        print(a)
