#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 29 17:01:31 2018

@author: aserrano
"""

import sys
import argparse
import pickle
import Functions

#if __name__ == "__main__":
parser = argparse.ArgumentParser()
parser.add_argument("-U", "--upload", dest = "upload", nargs = 1, help = "Uploads the original file and creates a temporal binary file to make the loading of the data faster.")
parser.add_argument("-G", "--general", dest = "general", nargs = 1, help = "Returns a general view of the gene you select, using a barplot")
parser.add_argument("-P1", "--pieplot1", dest = "pieplot1", nargs = "+", help = "Returns a pie plot showing the isoforms of the gene you select in order of expression.\nThe threshold is used to put the least expressed isoform as 'other'.")
parser.add_argument("-B", "--barplot", dest = "barplot", nargs = "+", help = "Returns a bar plot showing just the most expressed isoforms of the gene you select.\nThe threshold is used to put the least expressed isoform as 'other'.")
parser.add_argument("-SB", "--stackedbarplot", dest = "stackedbarplot", nargs = "+", help = "Returns a stacked bar plot showing just the most expressed isoforms of the gene you select.\nThe threshold is used to put the least expressed isoform as 'other'.")
parser.add_argument("-BS", "--bigsummary", dest = "bigsummary", nargs = "+", help = "Returns a dataframe with a conclusion for each gene.\nThe threshold is used to put the least expressed isoform as 'other'.\nThe threshold2 is the value you consider is significative.\n")
parser.add_argument("-St", "--stats", dest = "stats", nargs = "+", help = "Returns a list with the statistics of the dataframe of 'bigsummary' and a pie plot for an easy visualization.\nThe threshold is used to put the least expressed isoform as 'other'.\nThe threshold2 is the value you consider is significative.\n")
args = parser.parse_args()

upload = args.upload
general = args.general
pieplot1 = args.pieplot1
barplot = args.barplot
stackedbarplot = args.stackedbarplot
bigsummary = args.bigsummary
stats = args.stats

# U = upload the initial data --> analyse_rnaseq.py U file
if sys.argv[1] == 'U' or sys.argv[1] == 'upload':
    DATA = Functions.reading_data(sys.argv[2])
    f = open('temp_file','wb')
    pickle.dump(DATA, f)
else:
    try:
        pickle_file = open('temp_file', 'rb')
        DATA = pickle.load(pickle_file)
        
        if sys.argv[1] == 'G' or sys.argv[1] == 'general':
            # G = general view of the expression of 1 gene using barplots --> analyse_rnaseq.py G gene
            Functions.general_view(DATA, sys.argv[2])
            
        elif sys.argv[1] == 'P1' or sys.argv[1] == 'pieplot1':
            # P1 = pie plot of the distribution of the expression of 1 gene --> analyse_rnaseq.py P1 gene 'thres'
            if len(sys.argv) == 4:
                Functions.pie_plot(DATA, sys.argv[2], float(sys.argv[3]))
            else:
                Functions.pie_plot(DATA, sys.argv[2])
        
        elif sys.argv[1] == 'B' or sys.argv[1] == 'barplot':
            # B = sorted barplot of the distribution of the expression of 1 gene --> analyse_rnaseq.py B gene 'thres'
            if len(sys.argv) == 4:
                Functions.sorted_barplot(DATA, sys.argv[2], float(sys.argv[3]))
            else:
                Functions.sorted_barplot(DATA, sys.argv[2])
        
        elif sys.argv[1] == 'SB' or sys.argv[1] == 'stackedbarplot':
            # SB = sorted stacked barplot of the distribution of the expression of 1 gene --> analyse_rnaseq.py SB gene 'thres'
            if len(sys.argv) == 4:
                Functions.sorted_stacked_barplot(DATA, sys.argv[2], float(sys.argv[3]))
            else:
                Functions.sorted_stacked_barplot(DATA, sys.argv[2])
        
        elif sys.argv[1] == 'BS' or sys.argv[1] == 'bigsummary':
            # BS = dataframe with a summary of the classification of all the genes: Mono, Bi, Tri or Multi? --> analyse_rnaseq.py BS 'thres' 'thres2' 
            Functions.big_summary(DATA, int(sys.argv[2]), int(sys.argv[3]))
        
        elif sys.argv[1] == 'St' or sys.argv[1] == 'statistics':
            # St = statistics of the classification of all the genes: Mono, Bi, Tri or Multi? --> analyse_rnaseq.py St 'thres' 'thres2' 'plot'
            df = Functions.big_summary(DATA, sys.argv[2], sys.argv[3]).T
            df = df.astype('category')
            Functions.statistics(df, sys.argv[2])
        
        else:
            print("That's all!")
        
    except FileNotFoundError:
        print('You must charge the data before using any other function.\nPlease, type "analyse_rnaseq.py U file".')