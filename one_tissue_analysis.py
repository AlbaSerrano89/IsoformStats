#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  6 12:54:19 2018

@author: aserrano
"""

import argparse
import pprint
import Functions

parser = argparse.ArgumentParser()

parser.add_argument("data", help = "A csv or gzip file with the data you want to analyse.")
parser.add_argument("--thres1", dest = "thres1", required = False, help = "A threshold used to put the least expressed isoform as 'other'.")
parser.add_argument("--thres2", dest = "thres2", required = False, help = "A second threshold, in this case, to select the most significative isoforms.")
parser.add_argument("--thres3", dest = "thres3", required = False, help = "A third threshold, in this case, to determine the minimum of samples to take a gene into account.")
parser.add_argument("--bstissuefile", dest = "bstissuefile", required = False, help = "The name of the bigsummary file.")
parser.add_argument("--tissuestatsfile", dest = "tissuestatsfile", required = False, help = "The name of the tissue statistics file.")
parser.add_argument("--dffile", dest = "dffile", required = False, help = "True or False: do you want a file with the dataframe?.")

parser.add_argument("-BS", "--bigsummary", nargs = "*", dest = "bigsummary", help = "Returns a csv file with a conclusion for each gene.")
parser.add_argument("-St", "--stats", nargs = "*", dest = "stats", help = "Returns a list with the statistics of the dataframe of 'bigsummary'.")
parser.add_argument("-SB", "--summarybarplot", nargs = "*", dest = "sumbar", help = "Returns a barplot with the statistics of the dataframe of 'bigsummary'.")
parser.add_argument("-NP", "--notexppie", nargs = "*", dest = "notexppie", help = "Returns a pie plot for an easy visualization.")
parser.add_argument("-EP", "--exprpie", nargs = "*", dest = "exprpie", help = "Returns a list with the statistics of the dataframe of 'bigsummary' and a pie plot for an easy visualization.")

args = parser.parse_args()
tissuestatsfile = args.tissuestatsfile
thres1 = args.thres1
thres2 = args.thres2
thres3 = args.thres3
bstissuefile = args.bstissuefile

data = args.data
bigsummary = args.bigsummary
exprpie = args.exprpie
sum_bar = args.sumbar
notexppie = args.notexppie
stats = args.stats
dffile = args.dffile

DATA = Functions.reading_data(data)

if bigsummary != None:
    # python analyse_rnaseq.py csv_file -BS '--threshold thres' '--threshold2 thres2' '--threshold3 thres3' '--bstissuefile bigsummary_filename'
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
    
    if bstissuefile != None:
        bstissuefile = bstissuefile
    else:
        bstissuefile = '_genes_'
    
    Functions.big_summary(DATA, thres1, thres2, thres3, bstissuefile)

if stats != None:
    # python analyse_rnaseq.py csv_file -St '--bstissuefile bigsummary_filename' '--tissuestatsfile tissue_statistics_filename' '--dffile ¿'T'/'F'?'
    
    if tissuestatsfile != None:
        tissuestatsfile = tissuestatsfile
    else:
        tissuestatsfile = '_statistics.csv'
    
    if dffile == None or dffile == 'T':
        dffile = 'T'
    elif dffile == 'F':
        dffile = 'F'
    else:
        dffile = dffile
    
    a = Functions.statistics(DATA, bstissuefile, tissuestatsfile, dffile)
    
    if type(a) == dict:
        pprint.pprint(a)

if sum_bar != None:
    # python analyse_rnaseq.py csv_file -SB '--bstissuefile bigsummary_filename'
    Functions.stats_barplot(DATA, bstissuefile)
    
elif notexppie != None:
    # python analyse_rnaseq.py csv_file -NP '--bstissuefile bigsummary_filename'
    Functions.notexp_pie(DATA, bstissuefile)
    
elif exprpie != None:
    # python analyse_rnaseq.py csv_file -NE '--bstissuefile bigsummary_filename'
    Functions.expr_pie(DATA, bstissuefile)