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
parser.add_argument("--thres1", dest = "thres1", required = False, help = "A threshold used to put the least expressed isoform as 'other'.")
parser.add_argument("--thres2", dest = "thres2", required = False, help = "A second threshold, in this case, to select the most significative isoforms.")
parser.add_argument("--thres3", dest = "thres3", required = False, help = "A third threshold, in this case, to determine the minimum of samples to take a gene into account.")
parser.add_argument("--out_bsdir", dest = "out_bsdir", required = False, help = "The path where the bigsummary file will be saved.")
parser.add_argument("--out_bsfile", dest = "out_bsfile", required = False, help = "The name of the bigsummary file.")
parser.add_argument("--dfbsfile", dest = "dfbsfile", required = False, help = "True or False: do you want a file with the dataframe?")
parser.add_argument("--in_bsfile", dest = "in_bsfile", required = False, help = "The name of the bigsummary file.")
parser.add_argument("--out_statsfile", dest = "out_statsfile", required = False, help = "The name of the tissue statistics file.")
parser.add_argument("--freq_type", dest = "freq_type", required = False, help = "The name of the bigsummary file.")
parser.add_argument("--dfstatsfile", dest = "dfstatsfile", required = False, help = "True or False: do you want a file with the dataframe?")

parser.add_argument("-BS", "--bigsummary", nargs = "*", dest = "bigsummary", help = "Returns a csv file with a conclusion for each gene.")
parser.add_argument("-St", "--stats", nargs = "*", dest = "stats", help = "Returns a list with the statistics of the dataframe of 'bigsummary'.")

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
freq_type = args.freq_type
dfstatsfile = args.dfstatsfile

bigsummary = args.bigsummary
stats = args.stats

DATA = Functions.reading_data(data)

if out_bsfile != None:
    out_bsfile = out_bsfile
else:
    out_bsfile = '_genes_'

if thres1 == None:
    thres1 = 0.1

if thres2 == None:
    thres2 = 0.7

if thres3 == None:
    thres3 = 10

thres1 = float(thres1)
thres2 = float(thres2)
thres3 = int(thres3)

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
    if dfbsfile == None or dfbsfile == 'T':
        if out_bsdir == None:
            parser.error('An output directory is required.')
        
        dfbsfile = 'T'
    elif dfbsfile == 'F':
        dfbsfile = 'F'
    else:
        dfbsfile = dfbsfile

    a = Functions.big_summary(DATA, out_bsdir, dfbsfile, out_bsfile, thres1, thres2, thres3)

elif stats != None:
    a = Functions.statistics(DATA, in_bsfile, freq_type, out_statsfile, dfstatsfile)

try:
    if type(a) == pandas.core.frame.DataFrame:
        pprint.pprint(a)
except NameError:
    pass
