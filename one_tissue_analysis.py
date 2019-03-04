#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  6 12:54:19 2018

@author: aserrano
"""
import os
import sys
import argparse
import pandas as pd
import pprint
import Functions

parser = argparse.ArgumentParser()

parser.add_argument("--data", help = "A csv or gzip file with the data you want to analyse.")
parser.add_argument("--minexp", required = False, type = float, default = 0.8, help = "A threshold used to put the least expressed isoform as 'other'.")
parser.add_argument("--minsamps", required = False, type = int, default = 10, help = "A threshold to determine the minimum of samples to take a gene into account.")
parser.add_argument("--out_tsdir", required = False, type = str, default = os.getcwd(), help = "The name of the directory where the tissue summary file will be saved.")
parser.add_argument("--out_tsfile", required = False, type = str, default = '', help = "The name of the tissue summary file.")
parser.add_argument("--out_thresdir", required = False, type = str, default = '_DiffThres_Summaries', help = "The name of the directory where all the different tissue summary files will be saved.")
parser.add_argument("--seqexp", required = False, type = float, default = 0.05, help = "The difference between the minexp thresholds you want to compare.")
parser.add_argument("--seqnumsamps", nargs = '+', required = False, type = int, help = "A sequence of maximum number of samples you want to compare.")
parser.add_argument("--ncpus", required = False, type = int, default = 1, help = "The number of CPUs you will use.")
parser.add_argument("--num_cores", required = False, type = int, default = 1, help = "The number of cores you will use.")
parser.add_argument("--in_thresdir", required = False, help = "The name of the directory where all the different tissue summary files has been saved by tissuediffthressum.")
parser.add_argument("--in_tsfile", required = False, help = "The name of the tissue summary file.")
parser.add_argument("--savefile", required = False, action = 'store_true', help = "T/F: do you want to save a file with the statistics?")
parser.add_argument("--out_statsdir", required = False, type = str, default = os.getcwd(), help = "If savefile = T, the name of the directory where we want to save the tissue statistics file.")
parser.add_argument("--out_statsfile", required = False, type = str, default = '', help = "If savefile = T, the name of the tissue statistics file.")
parser.add_argument("--in_statsfile", required = False, help = "The name of the tissue statistics file.")
parser.add_argument("--genetype", nargs = '+', required = False, default = '', help = "The concrete gene type you are interested in.")
parser.add_argument("--drop_tsfile", required = False, action = 'store_true', help = "T/F: do you want to remove the tissue summary file?")
parser.add_argument("--in_thresfile", required = False, help = "The name of the file where the tissue different thresholds statistics file has been saved by tissuediffthresstats.")
parser.add_argument("--expressed", required = False, action = 'store_true', help = "T/F: do you want to see just the results of the expressed genes?")
parser.add_argument("--samplots", nargs = '+', required = False, default = '', help = "For how many samples do you want to see the histograms?")
parser.add_argument("--plotfile", required = False, help = "The name of the plot file.")

parser.add_argument("-TSu", "--tissuesummary", action = 'store_true', help = "Returns a csv file with a conclusion for each gene.")
parser.add_argument("-TSt", "--tissuestats", action = 'store_true', help = "Returns a list with the statistics of the dataframe of 'tissuesummary'.")
parser.add_argument("-TDSu", "--tissuediffthressum", action = 'store_true', help = "Creates a directory with different tissue summaries, changing the thresholds.")
parser.add_argument("-TDSt", "--tissuediffthresstats", action = 'store_true', help = "Creates a directory with the analysis of all the different tissue summaries created by tissuediffthressum.")
parser.add_argument("-TB", "--tissuebarplot", action = 'store_true', help = "Returns a barplot with the statistics of all the genes of the dataframe of 'tissuesummary'.")
parser.add_argument("-TDB", "--tissuediffthresplot", action = 'store_true', help = "Returns a set of histograms with a comparison of the different threshold results.")

args = parser.parse_args()

data = args.data
minexp = args.minexp
minsamps = args.minsamps
out_tsdir = args.out_tsdir
out_tsfile = args.out_tsfile
out_thresdir = args.out_thresdir
seqexp = args.seqexp
seqnumsamps = args.seqnumsamps
ncpus = args.ncpus
num_cores = args.num_cores
in_thresdir = args.in_thresdir
in_tsfile = args.in_tsfile
savefile = args.savefile
out_statsdir = args.out_statsdir
out_statsfile = args.out_statsfile
in_statsfile = args.in_statsfile
genetype = args.genetype
drop_tsfile = args.drop_tsfile
in_thresfile = args.in_thresfile
expressed = args.expressed
samplots = args.samplots
plotfile = args.plotfile

tissuesummary = args.tissuesummary
tissuestats = args.tissuestats
tissuediffthressum = args.tissuediffthressum
tissuediffthresstats = args.tissuediffthresstats
tissuebarplot = args.tissuebarplot
tissuediffthresplot = args.tissuediffthresplot

#print('seqnumsamps: {}.'.format(seqnumsamps))
# ---------------------- Functions ---------------------- #
pd.set_option('max_columns', 100)

if tissuesummary:
    Functions.tissue_summary(data, out_tsdir, out_tsfile, minexp, minsamps)

elif tissuestats:
    a = Functions.tissue_statistics(in_tsfile, savefile, out_statsdir, out_statsfile, genetype, drop_tsfile)

elif tissuediffthressum:
    Functions.tissue_difthres_summaries(data, seqnumsamps, out_thresdir, seqexp, ncpus, num_cores)

elif tissuediffthresstats:
    Functions.tissue_difthres_statistics(in_thresdir, out_statsdir, out_statsfile, genetype, drop_tsfile)

elif tissuebarplot:
    if in_statsfile != None:
        in_file = in_statsfile
    elif in_thresfile != None:
        in_file = in_thresfile
    
    Functions.tissue_barplot(in_file, plotfile, expressed, minexp, minsamps)

elif tissuediffthresplot:
    Functions.tissue_difthres_barplot(in_thresfile, plotfile, expressed, genetype, samplots)
else:
    sys.exit("ERROR: You did not write any function.")

try:
    if type(a) == pd.core.frame.DataFrame:
        pprint.pprint(a)
except NameError:
    pass