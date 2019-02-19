#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  5 16:14:47 2018

@author: aserrano
"""
import argparse
import Functions

parser = argparse.ArgumentParser()

parser.add_argument("--initial_dir", required = False, help = "The name of the directory where all initial data is saved.")
parser.add_argument("--seqnumsamps", nargs = '+', required = False, help = "A sequence of maximum number of samples you want to compare.")
parser.add_argument("--out_threstsdir", required = False, help = "The name of the directory where all the different tissue summary files has been saved.")
parser.add_argument("--seqexp", required = False, help = "The difference between the minexp thresholds you want to compare.")
parser.add_argument("--ncpus", required = False, help = "The number of CPUs you will use.")
parser.add_argument("--num_cores", required = False, help = "The number of cores you will use.")
parser.add_argument("--in_thresdir", required = False, help = "The name of the directory where all the different tissue summary files has been saved.")
parser.add_argument("--out_statsdir", required = False, help = "A folder where all the tissue statistics files are located.")
parser.add_argument("--in_statsdir", required = False, help = "A folder where all the tissue statistics files are located.")
parser.add_argument("--genetype", nargs = '+', required = False, help = "The concrete gene type you are interested in.")
parser.add_argument("--drop_tsfiles", required = False, help = "T/F: do you want to remove the tissue summary files?")
parser.add_argument("--minexp", required = False, help = "A threshold used to put the least expressed isoform as 'other'.")
parser.add_argument("--minsamps", required = False, help = "A threshold to determine the minimum of samples to take a gene into account.")
parser.add_argument("--plotfile", required = False, help = "The name of the plot file.")

parser.add_argument("-ATDTSum", "--alltissdiffthressum", nargs = "*", help = "Creates the directory and the files where all the statistics per tissue will be saved.")
parser.add_argument("-ATDTSta", "--alltissdiffthresstats", nargs = "*", help = "Creates the directory and the files where all the statistics per tissue will be saved.")
parser.add_argument("-ATB", "--alltissuesbarplot", nargs = "*", help = "Draws an horizontal barplot to show the differences between all the tissues.")
parser.add_argument("-ATBE", "--alltissuesbarplotexpressed", nargs = "*", help = "Draws an horizontal barplot to show the differences between the expressed genes of all the tissues.")

args = parser.parse_args()

initial_dir = args.initial_dir
seqnumsamps = args.seqnumsamps
out_threstsdir = args.out_threstsdir
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
plotfile = args.plotfile

alltissdiffthressum = args.alltissdiffthressum
alltissdiffthresstats = args.alltissdiffthresstats
alltissuesbarplot = args.alltissuesbarplot
alltissuesbarplotexpressed = args.alltissuesbarplotexpressed

if out_threstsdir == None:
    out_threstsdir = 'AllTissues_Summaries_DifThres/'

if seqexp == None:
    seqexp = 0.05

seqexp = float(seqexp)

if ncpus == None:
    ncpus = 1
else:
    int(ncpus)

if num_cores == None:
    num_cores = 1
else:
    int(num_cores)

if in_thresdir == None:
    in_thresdir = 'AllTissues_Summaries_DifThres/'

if out_statsdir == None:
    out_statsdir = 'AllTissues_Statistics_DifThres/'

if in_statsdir == None:
    in_statsdir = 'AllTissues_Statistics_DifThres/'

if genetype == None:
    genetype = ''

if drop_tsfiles == None:
    drop_tsfiles = 'T'
 
if minexp != None:
    minexp = float(minexp)

if minsamps != None:
    minsamps = int(minsamps)

if alltissdiffthressum != None:
    Functions.all_tissues_difthres_summaries(initial_dir, seqnumsamps, out_threstsdir, seqexp, ncpus, num_cores)

elif alltissdiffthresstats != None:
    Functions.all_tissues_difthres_statistics(in_thresdir, out_statsdir, genetype, drop_tsfiles)

elif alltissuesbarplot != None:
    if plotfile == None:
        plotfile = 'AllTissuesBarplot'
    
    Functions.all_tissues_barplot(in_statsdir, minexp, minsamps, plotfile, genetype)

elif alltissuesbarplotexpressed != None:
    if plotfile == None:
        plotfile = 'AllTissuesBarplotExpr'
    
    Functions.all_tissues_barplot_expr(in_statsdir, minexp, minsamps, plotfile, genetype)