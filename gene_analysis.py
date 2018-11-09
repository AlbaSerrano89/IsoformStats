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

parser.add_argument("data", help = "A csv or gzip file with the data you want to analyse.")
parser.add_argument("gene", help = "The gene to analyse.")
parser.add_argument("--outfile", dest = "outfile", required = False, help = "The name of the bigsummary file.")
parser.add_argument("--thres1", dest = "thres1", required = False, help = "A threshold used to put the least expressed isoform as 'other'.")
parser.add_argument("--thres2", dest = "thres2", required = False, help = "A second threshold, in this case, to select the most significative isoforms.")
parser.add_argument("--thres3", dest = "thres3", required = False, help = "A third threshold, in this case, to determine the minimum of samples to take a gene into account.")

parser.add_argument("-GN", "--genename", nargs = "*", dest = "genename", help = "Returns the name of the gene.")
parser.add_argument("-TP", "--totprop", nargs = "*", dest = "totprop", help = "Returns the mean proportion of each isoform of the gene.")
parser.add_argument("-SP", "--sumprop", nargs = "*", dest = "sumprop", help = "Returns a list of the isoforms and their respective total proportions, reversed ordered by these proportions.")
parser.add_argument("-IS", "--isofstats", nargs = "*", dest = "isofstats", help = "Returns a list with different proportions for each isoform.")
parser.add_argument("-TG", "--typegene", nargs = "*", dest = "typegene", help = "Returns a list with different proportions for each isoform.")

parser.add_argument("-G", "--general", nargs = "*", dest = "general", help = "Returns a general view of the gene you select, using a barplot.")
parser.add_argument("-P", "--pieplot", nargs = "*", dest = "pieplot", help = "Returns a pie plot showing the isoforms of the gene you select in order of expression.")
parser.add_argument("-B", "--stackedbarplot", nargs = "*", dest = "stackedbarplot", help = "Returns a stacked bar plot showing just the most expressed isoforms of the gene you select.")

args = parser.parse_args()

data = args.data
gene = args.gene
outfile = args.outfile
thres1 = args.thres1
thres2 = args.thres2
thres3 = args.thres3

genename = args.genename
totprop = args.totprop
sumprop = args.sumprop
isofstats = args.isofstats
typegene = args.typegene

general = args.general
pieplot = args.pieplot
stackedbarplot = args.stackedbarplot

DATA = Functions.reading_data(data)

if genename != None:
    # python gene_analysis.py csv_file genename/geneposition -GN 
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
    # python gene_analysis.py csv_file genename/geneposition -TP
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
    # python gene_analysis.py csv_file genename/geneposition -SP '--threshold thres'
    try:
        gene = int(gene)
    except ValueError:
        gene = gene

    if thres1 != None:
        thres1 = float(thres1)
    else:
        thres1 = 0.9

    a = Functions.summarized_props(DATA, gene, thres1)
    
    if 'not expressed' in a or 'is not in your data' in a:
        print(a)
    else:
        pprint.pprint(a)

elif isofstats != None:
    # python gene_analysis.py csv_file genename/geneposition -IS '--threshold1 thres1' '--threshold2 thres2'
    try:
        gene = int(gene)
    except ValueError:
        gene = gene

    if thres1 != None:
        thres1 = float(thres1)
    else:
        thres1 = 0.9

    if thres2 != None:
        thres2 = float(thres2)
    else:
        thres2 = 0.7

    a = Functions.isoform_stats(DATA, gene, thres1, thres2)
    
    if 'not expressed' in a or 'is not in your data' in a:
        print(a)
    else:
        pprint.pprint(a)

elif typegene != None:
    # python gene_analysis.py csv_file genename/geneposition -TG '--threshold1 thres1' '--threshold2 thres2' '--threshold3 thres3'
    try:
        gene = int(gene)
    except ValueError:
        gene = gene

    if thres1 != None:
        thres1 = float(thres1)
    else:
        thres1 = 0.9

    if thres2 != None:
        thres2 = float(thres2)
    else:
        thres2 = 0.7

    if thres3 != None:
        thres3 = int(thres3)
    else:
        thres3 = 10

    a = Functions.type_of_gene(DATA, gene, thres1, thres2, thres3)
    
    if 'not expressed' in a or 'is not in your data' in a:
        print(a)
    else:
        pprint.pprint(a)
    
elif general != None:
    # python gene_analysis.py csv_file genename/geneposition -G
    try:
        gene = int(gene)
    except ValueError:
        gene = gene

    a = Functions.general_view(DATA, gene)
    
    if type(a) == str:
        print(a)
    
elif pieplot != None:
    # python gene_analysis.py csv_file genename/geneposition -P '--threshold thres'
    try:
        gene = int(gene)
    except ValueError:
        gene = gene

    if thres1 != None:
        thres1 = float(thres1)
    else:
        thres1 = 0.9
        
    a = Functions.pie_plot(DATA, gene, thres1)

    if type(a) == str:
        print(a)
    
elif stackedbarplot != None:
    # python gene_analysis.py csv_file genename/geneposition -B '--threshold thres'
    try:
        gene = int(gene)
    except ValueError:
        gene = gene

    if thres1 != None:
        thres1 = float(thres1)
    else:
        thres1 = 0.9
    
    a = Functions.stacked_barplot(DATA, gene, thres1)

    if type(a) == str:
        print(a)