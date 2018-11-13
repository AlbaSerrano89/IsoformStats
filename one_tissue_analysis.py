#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  6 12:54:19 2018

@author: aserrano
"""

import argparse
import pandas
import pprint
import Functions

parser = argparse.ArgumentParser()

parser.add_argument("data", help = "A csv or gzip file with the data you want to analyse.")
parser.add_argument("--thres1", dest = "thres1", nargs = "+", required = False, help = "A threshold used to put the least expressed isoform as 'other'.")
parser.add_argument("--thres2", dest = "thres2", nargs = "+", required = False, help = "A second threshold, in this case, to select the most significative isoforms.")
parser.add_argument("--thres3", dest = "thres3", nargs = "+", required = False, help = "A third threshold, in this case, to determine the minimum of samples to take a gene into account.")
parser.add_argument("--out_bsdir", dest = "out_bsdir", required = False, help = "The path where the bigsummary file will be saved.")
parser.add_argument("--out_bsfile", dest = "out_bsfile", required = False, help = "The name of the bigsummary file.")
parser.add_argument("--dfbsfile", dest = "dfbsfile", required = False, help = "True or False: do you want a file with the dataframe?")
parser.add_argument("--in_bsfile", dest = "in_bsfile", required = False, help = "The name of the bigsummary file.")
parser.add_argument("--out_statsfile", dest = "out_statsfile", required = False, help = "The name of the tissue statistics file.")
parser.add_argument("--dfstatsfile", dest = "dfstatsfile", required = False, help = "True or False: do you want a file with the dataframe?")

parser.add_argument("-BS", "--bigsummary", nargs = "*", dest = "bigsummary", help = "Returns a csv file with a conclusion for each gene.")
parser.add_argument("-CS", "--compstats", nargs = "*", dest = "compstats", help = "Returns a list with the statistics of the dataframe of 'bigsummary'.")
parser.add_argument("-St", "--stats", nargs = "*", dest = "stats", help = "Returns a list with the statistics of the dataframe of 'bigsummary'.")
parser.add_argument("-SB", "--summarybarplot", nargs = "*", dest = "sumbar", help = "Returns a barplot with the statistics of the dataframe of 'bigsummary'.")
parser.add_argument("-EB", "--exprbarplot", nargs = "*", dest = "exprbar", help = "Returns a list with the statistics of the dataframe of 'bigsummary' and a pie plot for an easy visualization.")

args = parser.parse_args()

data = args.data
thres1 = args.thres1
thres2 = args.thres2
thres3 = args.thres3
out_bsdir = args.out_bsdir
out_bsfile = args.out_bsfile
dfbsfile = args.dfbsfile
in_bsfile = args.in_bsfile
out_statsfile = args.out_statsfile
dfstatsfile = args.dfstatsfile

bigsummary = args.bigsummary
compstats = args.compstats
stats = args.stats
sum_bar = args.sumbar
exprbar = args.exprbar

DATA = Functions.reading_data(data)

if thres1 == None:
    thres1 = 0.1

if thres2 == None:
    thres2 = 0.7

if thres3 == None:
    thres3 = 10

if bigsummary != None:
    # python analyse_rnaseq.py csv_file -BS '--threshold thres' '--threshold2 thres2' '--threshold3 thres3' '--out_bsfile bigsummary_filename'
    thres1 = float(thres1)
    thres2 = float(thres2)
    thres3 = int(thres3)
    
    if out_bsfile != None:
        out_bsfile = out_bsfile
    else:
        out_bsfile = '_genes_'
    
    if out_statsfile != None:
        out_statsfile = out_statsfile
    else:
        out_statsfile = '_statistics.csv'
    
    if dfbsfile == None or dfbsfile == 'T':
        if out_bsdir == None:
            parser.error('An output directory is required.')
        
        dfbsfile = 'T'
    elif dfbsfile == 'F':
        dfbsfile = 'F'
    else:
        dfbsfile = dfbsfile
    
    a = Functions.big_summary(DATA, out_bsdir, dfbsfile, out_bsfile, thres1, thres2, thres3)

    if type(a) == pandas.core.frame.DataFrame:
        pprint.pprint(a)

if compstats != None:
    # python analyse_rnaseq.py csv_file -CS '--threshold thres' '--threshold2 thres2' '--threshold3 thres3'
    Functions.compare_thres(DATA, thres1, thres2, thres3)

if stats != None:
    # python analyse_rnaseq.py csv_file -St '--in_bsfile bigsummary_filename' '--out_statsfile tissue_statistics_filename' '--dffile ¿'T'/'F'?'
    if dfstatsfile == None or dfstatsfile == 'T':
        if out_statsfile != None:
            out_statsfile = out_statsfile
        else:
            out_statsfile = '_statistics.csv'
        
        dfstatsfile = 'T'
    elif dfstatsfile == 'F':
        dfstatsfile = 'F'
    else:
        dfstatsfile = dfstatsfile

    a = Functions.statistics(DATA, in_bsfile, out_statsfile, dfstatsfile)
    
    if type(a) == pandas.core.frame.DataFrame:
        pprint.pprint(a)

if sum_bar != None:
    # python analyse_rnaseq.py csv_file -SB '--in_bsfile bigsummary_filename'
    Functions.stats_barplot(DATA, in_bsfile)
    
elif exprbar != None:
    # python analyse_rnaseq.py csv_file -EB '--in_bsfile bigsummary_filename'
    Functions.expr_barplot(DATA, in_bsfile)