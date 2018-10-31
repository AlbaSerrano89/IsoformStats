#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 29 17:01:31 2018

@author: aserrano
"""

import argparse
import pickle
import Functions

#if __name__ == "__main__":
parser = argparse.ArgumentParser()
parser.add_argument("-U", "--upload", dest = "upload", nargs = 1, help = "Uploads the original file and creates a temporal binary file to make the loading of the data faster.")
parser.add_argument("-G", "--general", dest = "general", nargs = 1, help = "Returns a general view of the gene you select, using a barplot")
parser.add_argument("-P", "--pieplot", dest = "pieplot", nargs = "+", help = "Returns a pie plot showing the isoforms of the gene you select in order of expression.\nThe threshold is used to put the least expressed isoform as 'other'.")
parser.add_argument("-B", "--barplot", dest = "barplot", nargs = "+", help = "Returns a bar plot showing just the most expressed isoforms of the gene you select.\nThe threshold is used to put the least expressed isoform as 'other'.")
parser.add_argument("-SB", "--stackedbarplot", dest = "stackedbarplot", nargs = "+", help = "Returns a stacked bar plot showing just the most expressed isoforms of the gene you select.\nThe threshold is used to put the least expressed isoform as 'other'.")
parser.add_argument("-BS", "--bigsummary", dest = "bigsummary", nargs = "*", help = "Returns a csv file with a conclusion for each gene.\nThe threshold is used to put the least expressed isoform as 'other'.\nThe threshold2 is the value you consider is significative.\n")
parser.add_argument("-E", "--expression", dest = "expsum", nargs = "+", help = "Returns a list with the statistics of the dataframe of 'bigsummary' and a pie plot for an easy visualization.\nThe threshold is used to put the least expressed isoform as 'other'.\nThe threshold2 is the value you consider is significative.\n")
parser.add_argument("-St", "--stats", dest = "stats", nargs = "+", help = "Returns a list with the statistics of the dataframe of 'bigsummary' and a pie plot for an easy visualization.\nThe threshold is used to put the least expressed isoform as 'other'.\nThe threshold2 is the value you consider is significative.\n")
args = parser.parse_args()

upload = args.upload
general = args.general
pieplot = args.pieplot
barplot = args.barplot
stackedbarplot = args.stackedbarplot
bigsummary = args.bigsummary
expr_sum = args.expsum
stats = args.stats

# U = upload the initial data --> analyse_rnaseq.py U file
if upload != None:
    DATA = Functions.reading_data(upload[0])
    f = open('temp_file','wb')
    pickle.dump(DATA, f)
else:
    try:
        pickle_file = open('temp_file', 'rb')
        DATA = pickle.load(pickle_file)
        
        if general != None:
            # G = general view of the expression of 1 gene using barplots --> analyse_rnaseq.py G gene
            try:
                gene = int(general[0])
            except ValueError:
                gene = general[0]
            
            Functions.general_view(DATA, gene)
            
        elif pieplot != None:
            # P = pie plot of the distribution of the expression of 1 gene --> analyse_rnaseq.py P1 gene 'thres'
            try:
                gene = int(pieplot[0])
            except ValueError:
                gene = pieplot[0]
            
            try:
                thres = float(pieplot[1])
            except IndexError:
                thres = 0.9
            
            Functions.pie_plot(DATA, gene, thres)
        
        elif barplot != None:
            # B = sorted barplot of the distribution of the expression of 1 gene --> analyse_rnaseq.py B gene 'thres'
            try:
                gene = int(barplot[0])
            except ValueError:
                gene = barplot[0]
            
            try:
                thres = float(barplot[1])
            except IndexError:
                thres = 0.9
            
            Functions.sorted_barplot(DATA, gene, thres)

        elif stackedbarplot != None:
            # SB = sorted stacked barplot of the distribution of the expression of 1 gene --> analyse_rnaseq.py SB gene 'thres'
            try:
                gene = int(stackedbarplot[0])
            except ValueError:
                gene = stackedbarplot[0]
            
            try:
                thres = float(stackedbarplot[1])
            except IndexError:
                thres = 0.9
            
            Functions.sorted_stacked_barplot(DATA, gene, thres)

        elif bigsummary != None:
            # BS = dataframe with a summary of the classification of all the genes: Mono, Bi, Tri or Multi? --> analyse_rnaseq.py BS 'thres' 'thres2' 
            try:
                thres = float(bigsummary[0])
            except IndexError:
                thres = 0.9

            try:
                thres2 = float(bigsummary[1])
            except IndexError:
                thres2 = 0.7
            
            print('holi')
            Functions.big_summary(DATA, thres, thres2)

        elif expr_sum != None:
            # E = classification of the gene: expressed or not expressed? and a pieplot --> analyse_rnaseq.py E csv_file 'thres' 'thres2'
            try:
                thres = float(expr_sum[1])
            except IndexError:
                thres = 0.9

            try:
                thres2 = float(expr_sum[2])
            except IndexError:
                thres2 = 0.7
            
            Functions.notexp_stats(expr_sum[0])
        
        elif stats != None:
            # St = statistics of the classification of all the genes: Mono, Bi, Tri or Multi? --> analyse_rnaseq.py St csv_file 'thres' 'thres2'
            try:
                thres = float(stats[1])
            except IndexError:
                thres = 0.9

            try:
                thres2 = float(stats[2])
            except IndexError:
                thres2 = 0.7
            
            Functions.statistics(stats[0])
#        
    except FileNotFoundError:
        print('You must charge the data before using any other function.\nPlease, type "python analyse_rnaseq.py -U csv_file".')