#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  5 16:14:47 2018

@author: aserrano
"""
import argparse
import Functions

parser = argparse.ArgumentParser()

parser.add_argument("--in_thresdir", required = False, help = "The name of the directory where all the different tissue summary files has been saved.")
parser.add_argument("--statsdir", required = False, help = "A folder where all the tissue statistics files are located.")
parser.add_argument("--genetype", nargs = '+', required = False, help = "The concrete gene type you are interested in.")
parser.add_argument("--drop_tsfiles", required = False, help = "T/F: do you want to remove the tissue summary files?")
parser.add_argument("--minexp", required = False, help = "A threshold used to put the least expressed isoform as 'other'.")
parser.add_argument("--minsamps", required = False, help = "A threshold to determine the minimum of samples to take a gene into account.")
parser.add_argument("--plotfile", required = False, help = "The name of the plot file.")

parser.add_argument("-ATDT", "--alltissdiffthres", nargs = "*", help = "Creates the directory and the files where all the statistics per tissue will be saved.")
parser.add_argument("-ATB", "--alltissuesbarplot", nargs = "*", help = "Draws an horizontal barplot to show the differences between all the tissues.")
parser.add_argument("-ATBE", "--alltissuesbarplotexpressed", nargs = "*", help = "Draws an horizontal barplot to show the differences between the expressed genes of all the tissues.")

args = parser.parse_args()

in_thresdir = args.in_thresdir
statsdir = args.statsdir
genetype = args.genetype
drop_tsfiles = args.drop_tsfiles
minexp = args.minexp
minsamps = args.minsamps
plotfile = args.plotfile

alltissdiffthres = args.alltissdiffthres
alltissuesbarplot = args.alltissuesbarplot
alltissuesbarplotexpressed = args.alltissuesbarplotexpressed

if statsdir == None:
    statsdir = 'diffthres_stats_alltissues/'

if genetype == None:
    genetype = ''

if drop_tsfiles == None:
    drop_tsfiles = 'T'
 
if minexp != None:
    minexp = float(minexp)

if minsamps != None:
    minsamps = int(minsamps)

if alltissdiffthres != None:
    Functions.all_tissues_difthres_statistics(in_thresdir, statsdir, genetype, drop_tsfiles)

elif alltissuesbarplot != None:
    if plotfile == None:
        plotfile = 'AllTissuesBarplot.pdf'
       
    Functions.all_tissues_barplot(statsdir, plotfile, minexp, minsamps, genetype)

elif alltissuesbarplotexpressed != None:
    if plotfile == None:
        plotfile = 'AllTissuesBarplot_Expressed.pdf'
       
    Functions.all_tissues_barplot_expr(statsdir, plotfile, minexp, minsamps, genetype)
