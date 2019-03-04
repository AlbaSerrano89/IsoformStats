#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  5 16:14:47 2018

@author: aserrano
"""
import sys
import argparse
import Functions

parser = argparse.ArgumentParser()

parser.add_argument("--initial_dir", required = False, help = "The name of the directory where all initial data is saved.")
parser.add_argument("--seqexp", required = False, type = float, default = 0.05, help = "The difference between the minexp thresholds you want to compare.")
parser.add_argument("--seqnumsamps", nargs = '+', required = False, type = int, help = "A sequence of maximum number of samples you want to compare.")
parser.add_argument("--ncpus", required = False, type = int, default = 1, help = "The number of CPUs you will use.")
parser.add_argument("--num_cores", required = False, type = int, default = 1, help = "The number of cores you will use.")
parser.add_argument("--out_thresdir", required = False, type = str, default = 'AllTissues_Summaries_DifThres/', help = "The name of the directory where all the different thresholds tissue summary files will be saved.")
parser.add_argument("--in_thresdir", required = False, type = str, default = 'AllTissues_Summaries_DifThres/', help = "The name of the directory where all the different tissue summary files has been saved.")
parser.add_argument("--out_statsdir", required = False, type = str, default = 'AllTissues_Statistics_DifThres/', help = "A folder where all the tissue statistics files are located.")
parser.add_argument("--in_statsdir", required = False, type = str, default = 'AllTissues_Statistics_DifThres/', help = "A folder where all the tissue statistics files are located.")
parser.add_argument("--genetype", nargs = '+', required = False, default = '', help = "The concrete gene type you are interested in.")
parser.add_argument("--drop_tsfiles", required = False, action = 'store_true', help = "T/F: do you want to remove the tissue summary files?")
parser.add_argument("--minexp", required = False, type = float, default = 0.8, help = "A threshold used to put the least expressed isoform as 'other'.")
parser.add_argument("--minsamps", required = False, type = int, default = 10, help = "A threshold to determine the minimum of samples to take a gene into account.")
parser.add_argument("--expressed", required = False, action = 'store_true', help = "T/F: do you want to see just the results of the expressed genes?")
parser.add_argument("--plotfile", required = False, type = str, default = 'AllTissuesBarplot', help = "The name of the plot file.")

parser.add_argument("-ATDTSum", "--alltissdiffthressum", nargs = "*", help = "Creates the directory and the files where all the statistics per tissue will be saved.")
parser.add_argument("-ATDTSta", "--alltissdiffthresstats", nargs = "*", help = "Creates the directory and the files where all the statistics per tissue will be saved.")
parser.add_argument("-ATB", "--alltissuesbarplot", nargs = "*", help = "Draws an horizontal barplot to show the differences between all the tissues.")

args = parser.parse_args()

initial_dir = args.initial_dir
seqnumsamps = args.seqnumsamps
out_thresdir = args.out_thresdir
seqexp = args.seqexp
ncpus = args.ncpus
num_cores = args.num_cores
in_thresdir = args.in_thresdir
out_statsdir = args.out_statsdir
in_statsdir = args.in_statsdir
genetype = args.genetype
drop_tsfiles = args.drop_tsfiles
minexp = args.minexp
minsamps = args.minsamps
expressed = args.expressed
plotfile = args.plotfile

alltissdiffthressum = args.alltissdiffthressum
alltissdiffthresstats = args.alltissdiffthresstats
alltissuesbarplot = args.alltissuesbarplot

if alltissdiffthressum != None:
    Functions.all_tissues_difthres_summaries(initial_dir, seqnumsamps, out_thresdir, seqexp, ncpus, num_cores)

elif alltissdiffthresstats != None:
    Functions.all_tissues_difthres_statistics(in_thresdir, out_statsdir, genetype, drop_tsfiles)

elif alltissuesbarplot != None:
    Functions.all_tissues_barplot(in_statsdir, minexp, minsamps, expressed, plotfile, genetype)
else:
    sys.exit("ERROR: You did not write any function.")