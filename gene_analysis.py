#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 29 17:01:31 2018

@author: aserrano
"""

import argparse
import pprint
import Functions

parser = argparse.ArgumentParser()

parser.add_argument("--data", required = True, help = "A csv or gzip file with the data you want to analyse.")
parser.add_argument("--genename", required = False, help = "The gene to analyse.")
parser.add_argument("--genepos", required = False, help = "Returns the name of the gene.")
parser.add_argument("--ncols", required = False, help = "Returns the name of the gene.")
parser.add_argument("--outfile", required = False, help = "The name of the bigsummary file.")
parser.add_argument("--thres1", required = False, help = "A threshold used to put the least expressed isoform as 'other'.")
parser.add_argument("--thres2", required = False, help = "A second threshold, in this case, to select the most significative isoforms.")
parser.add_argument("--thres3", required = False, help = "A third threshold, in this case, to determine the minimum of samples to take a gene into account.")

parser.add_argument("-GI", "--geneinfo", nargs = "*", help = "Returns the name of the gene.")
parser.add_argument("-GN", "--genenameF", nargs = "*", help = "Returns the name of the gene.")
parser.add_argument("-TP", "--totprop", nargs = "*", help = "Returns the mean proportion of each isoform of the gene.")
parser.add_argument("-SP", "--sumprop", nargs = "*", help = "Returns a list of the isoforms and their respective total proportions, reversed ordered by these proportions.")
parser.add_argument("-IS", "--isofstats", nargs = "*", help = "Returns a list with different proportions for each isoform.")
parser.add_argument("-TG", "--typegene", nargs = "*", help = "Returns a list with different proportions for each isoform.")

parser.add_argument("-G", "--general", nargs = "*", help = "Returns a general view of the gene you select, using a barplot.")
parser.add_argument("-P", "--pieplot", nargs = "*", help = "Returns a pie plot showing the isoforms of the gene you select in order of expression.")
parser.add_argument("-B", "--stackedbarplot", nargs = "*", help = "Returns a stacked bar plot showing just the most expressed isoforms of the gene you select.")

args = parser.parse_args()

data = args.data
genename = args.genename
genepos = args.genepos
ncols = args.ncols
outfile = args.outfile
thres1 = args.thres1
thres2 = args.thres2
thres3 = args.thres3

geneinfo = args.geneinfo
genenameF = args.genenameF
totprop = args.totprop
sumprop = args.sumprop
isofstats = args.isofstats
typegene = args.typegene

general = args.general
pieplot = args.pieplot
stackedbarplot = args.stackedbarplot

if genename == None and genepos == None:
    print('ERROR: You must specify the gene you want to analyse.')
elif genename == None:
    gene = genepos
else:
    gene = genename

if ncols != None:
    ncols = int(ncols)
else:
    ncols = 10

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

if geneinfo != None:
    # python gene_analysis.py '--data csv_gz_file' '--genename genename' '--genepos geneposition' -GI '--ncols numberofcolumns'
    try:
        gene = int(gene)
    except ValueError:
        gene = gene

    a = Functions.gene_info(DATA, gene, ncols)
    
    if 'is not in your data' in a:
        print('ERROR: ' + a)
    else:
        pprint.pprint(a)

elif genenameF != None:
    # python gene_analysis.py '--data csv_gz_file' '--genename genename' '--genepos geneposition' -GN 
    try:
        gene = int(gene)
    except ValueError:
        gene = gene

    a = Functions.formatting_gene(DATA, gene)
    
    if 'is not in your data' in a:
        print('ERROR: ' + a)
    else:
        print(a)

elif totprop != None:
    # python gene_analysis.py '--data csv_gz_file' '--genename genename' '--genepos geneposition' -TP
    try:
        gene = int(gene)
    except ValueError:
        gene = gene

    a = Functions.total_prop(DATA, gene)
    
    if 'not expressed' in a or 'is not in your data' in a:
        print(a)
    else:
        pprint.pprint(a)

elif sumprop != None:
    # python gene_analysis.py '--data csv_gz_file' '--genename genename' '--genepos geneposition' -SP '--thres1 thres1'
    try:
        gene = int(gene)
    except ValueError:
        gene = gene

    a = Functions.summarized_props(DATA, gene, thres1)
    
    if 'not expressed' in a or 'is not in your data' in a:
        print(a)
    else:
        pprint.pprint(a)

elif isofstats != None:
    # python gene_analysis.py '--data csv_gz_file' '--genename genename' '--genepos geneposition' -IS '--thres1 thres1' '--thres2 thres2'
    try:
        gene = int(gene)
    except ValueError:
        gene = gene

    a = Functions.isoform_stats(DATA, gene, thres1, thres2)
    
    if 'not expressed' in a or 'is not in your data' in a:
        print(a)
    else:
        pprint.pprint(a)

elif typegene != None:
    # python gene_analysis.py '--data csv_gz_file' '--genename genename' '--genepos geneposition' -TG '--thres1 thres1' '--thres2 thres2' '--thres3 thres3'
    try:
        gene = int(gene)
    except ValueError:
        gene = gene

    a = Functions.type_of_gene(DATA, gene, thres1, thres2, thres3)
    
    if 'not expressed' in a or 'is not in your data' in a:
        print(a)
    else:
        pprint.pprint(a)
    
elif general != None:
    # python gene_analysis.py '--data csv_gz_file' '--genename genename' '--genepos geneposition' -G
    try:
        gene = int(gene)
    except ValueError:
        gene = gene

    a = Functions.general_view(DATA, gene)
    
    if type(a) == str:
        print(a)
    
elif pieplot != None:
    # python gene_analysis.py '--data csv_gz_file' '--genename genename' '--genepos geneposition' -P '--thres1 thres1'
    try:
        gene = int(gene)
    except ValueError:
        gene = gene

    a = Functions.pie_plot(DATA, gene, thres1)

    if type(a) == str:
        print(a)
    
elif stackedbarplot != None:
    # python gene_analysis.py '--data csv_gz_file' '--genename genename' '--genepos geneposition' -B '--thres1 thres1'
    try:
        gene = int(gene)
    except ValueError:
        gene = gene

    a = Functions.stacked_barplot(DATA, gene, thres1)

    if type(a) == str:
        print(a)