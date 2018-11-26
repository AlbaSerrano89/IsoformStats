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
parser.add_argument("--minexp", dest = "minexp", required = False, help = "A threshold used to put the least expressed isoform as 'other'.")
parser.add_argument("--mintotexp", dest = "mintotexp", required = False, help = "A second threshold, in this case, to select the most significative isoforms.")
parser.add_argument("--minsamps", dest = "minsamps", required = False, help = "A third threshold, in this case, to determine the minimum of samples to take a gene into account.")
parser.add_argument("--out_bsdir", dest = "out_bsdir", required = False, help = "The path where the bigsummary file will be saved.")
parser.add_argument("--out_bsfile", dest = "out_bsfile", required = False, help = "The name of the bigsummary file.")
parser.add_argument("--dfbsfile", dest = "dfbsfile", required = False, help = "True or False: do you want a file with the dataframe?")
parser.add_argument("--in_bsfile", dest = "in_bsfile", required = False, help = "The name of the bigsummary file.")
parser.add_argument("--out_statsfile", dest = "out_statsfile", required = False, help = "The name of the tissue statistics file.")
parser.add_argument("--freq_type", dest = "freq_type", required = False, help = "The name of the bigsummary file.")
parser.add_argument("--dfstatsfile", dest = "dfstatsfile", required = False, help = "True or False: do you want a file with the dataframe?")

parser.add_argument("-BS", "--bigsummary", nargs = "*", dest = "bigsummary", help = "Returns a csv file with a conclusion for each gene.")
#parser.add_argument("-St", "--stats", nargs = "*", dest = "stats", help = "Returns a list with the statistics of the dataframe of 'bigsummary'.")
#parser.add_argument("-SB", "--summarybarplot", nargs = "*", dest = "sumbar", help = "Returns a barplot with the statistics of the dataframe of 'bigsummary'.")
#parser.add_argument("-EB", "--exprbarplot", nargs = "*", dest = "exprbar", help = "Returns a list with the statistics of the dataframe of 'bigsummary' and a pie plot for an easy visualization.")

args = parser.parse_args()

data = args.data
minexp = args.minexp
mintotexp = args.mintotexp
minsamps = args.minsamps
out_bsdir = args.out_bsdir
out_bsfile = args.out_bsfile
dfbsfile = args.dfbsfile
in_bsfile = args.in_bsfile
out_statsfile = args.out_statsfile
freq_type = args.freq_type
dfstatsfile = args.dfstatsfile

bigsummary = args.bigsummary
#stats = args.stats
#sum_bar = args.sumbar
#exprbar = args.exprbar

#DATA = Functions.reading_data(data)

if out_bsfile != None:
    out_bsfile = out_bsfile
else:
    out_bsfile = '_genes_'

if minexp == None:
    minexp = 0.1

if mintotexp == None:
    mintotexp = 0.7

if minsamps == None:
    minsamps = 10

minexp = float(minexp)
mintotexp = float(mintotexp)
minsamps = int(minsamps)

if out_statsfile != None:
    out_statsfile = out_statsfile
else:
    out_statsfile = '_statistics.csv'

if freq_type == None:
    freq_type = 'rel'

if dfstatsfile == None:
    dfstatsfile = 'F'

if dfstatsfile == 'T':
    if out_statsfile == None:
        out_statsfile = '_statistics.csv'

if bigsummary != None:
    if out_bsdir == None:
        parser.error('An output directory is required.')

    #a = Functions.big_summary(DATA, out_bsdir, out_bsfile, minexp, mintotexp, minsamps)
    a = Functions.big_summary(args.data, out_bsdir, out_bsfile, minexp, mintotexp, minsamps)

#if stats != None:
#    a = Functions.statistics(DATA, in_bsfile, freq_type, out_statsfile, dfstatsfile)
#
#if sum_bar != None:
#    Functions.stats_barplot(DATA, in_bsfile)
#    
#elif exprbar != None:
#    Functions.expr_barplot(DATA, in_bsfile)

try:
    if type(a) == pandas.core.frame.DataFrame:
        pprint.pprint(a)
except NameError:
    pass