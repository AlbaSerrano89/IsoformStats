#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 29 17:01:31 2018

@author: aserrano
"""

import argparse
import Functions

#if __name__ == "__main__":
parser = argparse.ArgumentParser()

parser.add_argument("data", help = "A csv or gzip file with the data you want to analyse.")
parser.add_argument("--gene", dest = "gene", required = False, help = "The gene to analyse.")
parser.add_argument("--outfile", dest = "outfile", required = False, help = "The name of the bigsummary file.")
parser.add_argument("--threshold", dest = "thres", required = False, help = "A threshold used to put the least expressed isoform as 'other'.")
parser.add_argument("--threshold2", dest = "thres2", required = False, help = "A second threshold, in this case, to select the most significative isoforms.")
parser.add_argument("--tissuefile", dest = "tissuefile", required = False, help = "The name of the bigsummary file.")
parser.add_argument("--tissue", dest = "tissue", required = False, help = "The tissue of the data.")

parser.add_argument("-G", "--general", nargs = "*", dest = "general", help = "Returns a general view of the gene you select, using a barplot.")
parser.add_argument("-P", "--pieplot", nargs = "*", dest = "pieplot", help = "Returns a pie plot showing the isoforms of the gene you select in order of expression.")
parser.add_argument("-B", "--stackedbarplot", nargs = "*", dest = "stackedbarplot", help = "Returns a stacked bar plot showing just the most expressed isoforms of the gene you select.")
parser.add_argument("-BS", "--bigsummary", nargs = "*", dest = "bigsummary", help = "Returns a csv file with a conclusion for each gene.")
parser.add_argument("-SB", "--summarybarplot", nargs = "*", dest = "sumbar", help = "Returns a barplot with the statistics of the dataframe of 'bigsummary'.")
parser.add_argument("-E", "--expression", nargs = "*", dest = "expsum", help = "Returns a list with the statistics of the dataframe of 'bigsummary' and a pie plot for an easy visualization.")
parser.add_argument("-St", "--stats", nargs = "*", dest = "stats", help = "Returns a list with the statistics of the dataframe of 'bigsummary' and a pie plot for an easy visualization.")

args = parser.parse_args()

data = args.data
general = args.general
pieplot = args.pieplot
stackedbarplot = args.stackedbarplot
bigsummary = args.bigsummary
sum_bar = args.sumbar
expr_sum = args.expsum
stats = args.stats

gene = args.gene
outfile = args.outfile
tissuefile = args.tissuefile
thres = args.thres
thres2 = args.thres2
tissue = args.tissue

DATA = Functions.reading_data(data)

if general != None:
    # python analyse_rnaseq.py csv_file -G --gene genename/geneposition
    # G = general view of the expression of 1 gene using a barplot
    if gene != None:
        try:
            gene1 = int(gene)
        except ValueError:
            gene1 = gene
    
        a = Functions.general_view(DATA, gene1)
        
        if type(a) == str:
            print(a)
    else:
        parser.error("You have to specify the gene name or the gene position in your list of genes.")
    
elif pieplot != None:
    # python analyse_rnaseq.py csv_file -P --gene genename/geneposition '--threshold thres'
    # P = pie plot of the distribution of the expression of 1 gene
    if gene != None:
        try:
            gene1 = int(gene)
        except ValueError:
            gene1 = gene
    
        if thres != None:
            thres1 = float(thres)
        else:
            thres1 = 0.9
            
        a = Functions.pie_plot(DATA, gene1, thres1)
    
        if type(a) == str:
            print(a)
    else:
        parser.error("You have to specify the gene name or the gene position in your list of genes.")
    
elif stackedbarplot != None:
    # python analyse_rnaseq.py csv_file -B --gene genename/geneposition '--threshold thres'
    # B = sorted stacked barplot of the distribution of the expression of 1 gene
    if gene != None:
        try:
            gene1 = int(gene)
        except ValueError:
            gene1 = gene
    
        if thres != None:
            thres1 = float(thres)
        else:
            thres1 = 0.9
        
        a = Functions.stacked_barplot(DATA, gene1, thres1)
    
        if type(a) == str:
            print(a)
    else:
        parser.error("You have to specify the gene name or the gene position in your list of genes.")
    
elif bigsummary != None:
    # python analyse_rnaseq.py csv_file -BS '--threshold thres' '--threshold2 thres2' '--outfile bigsummary_filename'
    # BS = dataframe with a summary of the classification of all the genes: NotExpressed, Mono, Bi, Tri or Multi?
    if thres != None:
        thres1 = float(thres)
    else:
        thres1 = 0.9

    if thres2 != None:
        thres3 = float(thres2)
    else:
        thres3 = 0.7
    
    if outfile != None:
        filename = outfile
    else:
        filename = 'classif_genes_'
    
    Functions.big_summary(DATA, thres1, thres3, filename)

if sum_bar != None:
    # python analyse_rnaseq.py csv_file -SB '--tissuefile bigsummary_filename' '--tissue tissue_name'
    # SB = barplot with the classification of the genes
    if tissuefile != None:
        filename = tissuefile
    else:
        filename = 'classif_genes_90.0_70.0.csv'
    
    if tissue != None:
        tissue1 = tissue
    else:
        if data.endswith('gz'):
            tissue1 = data[:-3]
        else:
            tissue1 = data[:-4]
    
    Functions.stats_barplot(filename, tissue1)
    
elif expr_sum != None:
    # python analyse_rnaseq.py csv_file -E '--tissuefile bigsummary_filename' '--tissue tissue_name'
    # E = classification of the genes: expressed or not expressed? and a pieplot
    if tissuefile != None:
        filename = tissuefile
    else:
        filename = 'classif_genes_90.0_70.0.csv'
    
    if tissue != None:
        tissue1 = tissue
    else:
        if data.endswith('gz'):
            tissue1 = data[:-3]
        else:
            tissue1 = data[:-4]
    
    Functions.notexp_stats(filename, tissue1)
    
elif stats != None:
    # python analyse_rnaseq.py csv_file -St '--tissuefile bigsummary_filename' '--tissue tissue_name'
    # St = statistics of the classification of all the genes: Mono, Bi, Tri or Multi?
    if tissuefile != None:
        filename = tissuefile
    else:
        filename = 'classif_genes_90.0_70.0.csv'
    
    if tissue != None:
        tissue1 = tissue
    else:
        if data.endswith('gz'):
            tissue1 = data[:-3]
        else:
            tissue1 = data[:-4]
    
    Functions.statistics(filename, tissue1)