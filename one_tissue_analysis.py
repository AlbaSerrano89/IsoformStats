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

parser.add_argument("--data", help = "A csv or gzip file with the data you want to analyse.")
parser.add_argument("--minexp", required = False, help = "A threshold used to put the least expressed isoform as 'other'.")
parser.add_argument("--minsamps", required = False, help = "A threshold to determine the minimum of samples to take a gene into account.")
parser.add_argument("--out_tsdir", required = False, help = "The name of the directory where the tissue summary file will be saved.")
parser.add_argument("--out_tsfile", required = False, help = "The name of the tissue summary file.")
parser.add_argument("--out_thresdir", required = False, help = "The name of the directory where all the different tissue summary files will be saved.")
parser.add_argument("--seqexp", required = False, help = "The difference between the minexp thresholds you want to compare.")
parser.add_argument("--seqnumsamps", nargs = '+', required = False, help = "The maximum number of samples you want to compare.")
parser.add_argument("--in_thresdir", required = False, help = "The name of the directory where all the different tissue summary files has been saved by tissuediffthressum.")
parser.add_argument("--in_tsfile", required = False, help = "The name of the tissue summary file.")
parser.add_argument("--savefile", required = False, help = "T/F: do you want to save a file with the statistics?")
parser.add_argument("--out_statsfile", required = False, help = "If savefile = T, the name of the tissue statistics file.")
parser.add_argument("--in_statsfile", required = False, help = "The name of the tissue statistics file.")
parser.add_argument("--genetype", nargs = '+', required = False, help = "The concrete gene type you are interested in.")
parser.add_argument("--drop_tsfile", required = False, help = "T/F: do you want to remove the tissue summary file?")

parser.add_argument("-TSu", "--tissuesummary", nargs = "*", help = "Returns a csv file with a conclusion for each gene.")
parser.add_argument("-TSt", "--tissuestats", nargs = "*", help = "Returns a list with the statistics of the dataframe of 'tissuesummary'.")
parser.add_argument("-TDSu", "--tissuediffthressum", nargs = "*", help = "Creates a directory with different tissue summaries, changing the thresholds.")
parser.add_argument("-TDSt", "--tissuediffthresstats", nargs = "*", help = "Creates a directory with the analysis of all the different tissue summaries created by tissuediffthressum.")
parser.add_argument("-TAB", "--tissueallbarplot", nargs = "*", help = "Returns a barplot with the statistics of all the genes of the dataframe of 'tissuesummary'.")
parser.add_argument("-TEB", "--tissueexprbarplot", nargs = "*", help = "Returns a barplot with the statistics of the expressed genes of the dataframe of 'tissuesummary'.")

args = parser.parse_args()

data = args.data
minexp = args.minexp
minsamps = args.minsamps
out_tsdir = args.out_tsdir
out_tsfile = args.out_tsfile
out_thresdir = args.out_thresdir
seqexp = args.seqexp
seqnumsamps = args.seqnumsamps
in_thresdir = args.in_thresdir
in_tsfile = args.in_tsfile
savefile = args.savefile
out_statsfile = args.out_statsfile
in_statsfile = args.in_statsfile
genetype = args.genetype
drop_tsfile = args.drop_tsfile

tissuesummary = args.tissuesummary
tissuestats = args.tissuestats
tissuediffthressum = args.tissuediffthressum
tissuediffthresstats = args.tissuediffthresstats
tissueallbarplot = args.tissueallbarplot
tissueexprbarplot = args.tissueexprbarplot

if minexp == None:
    minexp = 0.1

if minsamps == None:
    minsamps = 10

if out_tsdir == None:
    out_tsdir = 'results'

if out_tsfile == None:
    out_tsfile = ''

if out_thresdir == None:
    out_thresdir = 'different_thresholds'

if seqexp == None:
    seqexp = 0.05

minexp = float(minexp)
minsamps = int(minsamps)
seqexp = float(seqexp)

if seqnumsamps != None:
    seqnumsamps = [int(x) for x in seqnumsamps]

if savefile == None or savefile == 'T':
    savefile = 'T'
    
    if out_statsfile == None:
        out_statsfile = ''

if genetype == None:
    genetype = ''

if drop_tsfile == None:
    drop_tsfile = 'T'

if tissuesummary != None:
    Functions.tissue_summary(data, out_tsdir, out_tsfile, minexp, minsamps)

elif tissuestats != None:
    a = Functions.tissue_statistics(in_tsfile, savefile, out_statsfile, genetype, drop_tsfile)

elif tissueallbarplot != None:
    Functions.tissue_all_barplot(in_statsfile)
    
elif tissueexprbarplot != None:
    Functions.tissue_exp_barplot(in_statsfile)

elif tissuediffthressum != None:
    Functions.tissue_difthres_summaries(data, seqnumsamps, out_thresdir, seqexp)

elif tissuediffthresstats != None:
    Functions.tissue_difthres_statistics(in_thresdir, genetype, drop_tsfile)

try:
    if type(a) == pandas.core.frame.DataFrame:
        pprint.pprint(a)
except NameError:
    pass