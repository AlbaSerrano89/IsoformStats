#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  5 16:14:47 2018

@author: aserrano
"""

import argparse
import Functions

parser = argparse.ArgumentParser()

parser.add_argument("--originaldir", dest = "originaldir", help = "A folder where all the tissue to analyse are located.")
parser.add_argument("--pref", dest = "prefix", help = "In case there were more than one type of data and a prefix diffenciates them.")
parser.add_argument("--bsdir", dest = "bigsummariesdir", help = "The directory where the big summaries will be saved.")
parser.add_argument("--statsfile", dest = "statsfile", help = "The name of the csv file where the statistics will be saved.")

parser.add_argument("--thres1", dest = "thres1", required = False, help = "A threshold used to put the least expressed isoform as 'other'.")
parser.add_argument("--thres2", dest = "thres2", required = False, help = "A second threshold, in this case, to select the most significative isoforms.")
parser.add_argument("--thres3", dest = "thres3", required = False, help = "A third threshold, in this case, to determine the minimum of samples to take a gene into account.")

parser.add_argument("-ATSum", "--alltissuesummaries", nargs = "*", dest = "alltissuesummaries", help = "Analyses all the RNA-seq files of a directory and returns a new directory with a big summary for each tissue.")
parser.add_argument("-ATSta", "--alltissuestats", nargs = "*", dest = "alltissuestats", help = "Returns the statistics of all tissues.")

args = parser.parse_args()

originaldir = args.originaldir
bigsummariesdir = args.bigsummariesdir
statsfile = args.statsfile

prefix = args.prefix
thres1 = args.thres1
thres2 = args.thres2
thres3 = args.thres3

alltissuesummaries = args.alltissuesummaries
alltissuestats = args.alltissuestats

if prefix == None:
    prefix = ''

if bigsummariesdir == None:
    bigsummariesdir = ''

if statsfile == None:
    statsfile = 'all_tissues_statistics.csv'

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

if alltissuesummaries != None:
    # python all_tissue_analysis.py -ATS --originaldir '--pref prefix' '--bsdir big_summaries_dir' '--thres1 thres1' '--thres2 thres2' '--thres3 thres3'
    Functions.all_tissues_analysis(originaldir, thres1, thres2, thres3, prefix, bigsummariesdir)
    
elif alltissuestats != None:
    # python all_tissue_analysis.py -TS --bsdir big_summaries_dir' '--statsfile statistics_filename'
    Functions.all_tissues_stats(bigsummariesdir, statsfile)