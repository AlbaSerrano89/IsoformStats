#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  9 07:55:04 2018

@author: aserrano
"""
# I ------------------------ IMPORTING LIBRARIES ------------------------
import Functions
import pickle
# F ------------------------ IMPORTING LIBRARIES ------------------------

# I ------------------------ READING AND FORMATING THE FILE ------------------------
DATA = Functions.reading_data('SMTSD_Fallopian_Tube.csv')
# F ------------------------ READING AND FORMATING THE FILE ------------------------

f = open('temp_file','wb')
pickle.dump(DATA, f)


pickle_file = open('temp_file2', 'rb')
todo = pickle.load(pickle_file)
print(todo)
Functions.general_view(todo, 0)
# I ------------------------ INITIAL PLOT: DISTRIBUTION PER GENE ------------------------
Functions.general_view(DATA, 'ENSG00000000003.10')
Functions.general_view(DATA, 0)
Functions.general_view(DATA, 1)
Functions.general_view(DATA, 15)
Functions.general_view(DATA, 100)
# F ------------------------ INITIAL PLOT: DISTRIBUTION PER GENE ------------------------
Functions.total_prop(DATA, 0)
# I ------------------------ SORTING PROPORTIONS OF EACH ISOFORM PER GENE ------------------------
Functions.summarized_props(DATA, 0, 0.8)
Functions.summarized_props(DATA, 0)
Functions.summarized_props(DATA, 1)
Functions.summarized_props(DATA, 15)
Functions.summarized_props(DATA, 100)
# F ------------------------ SORTING PROPORTIONS OF EACH ISOFORM PER GENE ------------------------

# I ------------------------ PROPORTIONS PIE PLOT: DISTRIBUTION PER GENE ------------------------
Functions.pie_plot(DATA, 0)
Functions.pie_plot(DATA, 0, 0.8)
Functions.pie_plot(DATA, 1)
Functions.pie_plot(DATA, 15)
Functions.pie_plot(DATA, 100)
# F ------------------------ PROPORTIONS PIE PLOT: DISTRIBUTION PER GENE ------------------------

# I ------------------------ PROPORTIONS BARPLOT: DISTRIBUTION PER GENE ------------------------
Functions.sorted_barplot(DATA, 0)
Functions.sorted_barplot(DATA, 0, 0.8)
Functions.sorted_barplot(DATA, 1)
Functions.sorted_barplot(DATA, 15)
Functions.sorted_barplot(DATA, 100)
# F ------------------------ PROPORTIONS BARPLOT: DISTRIBUTION PER GENE ------------------------

# I ------------------------ PROPORTIONS BARPLOT: DISTRIBUTION PER GENE ------------------------
Functions.sorted_stacked_barplot(DATA, 0)
Functions.sorted_stacked_barplot(DATA, 0, 0.8)
Functions.sorted_stacked_barplot(DATA, 1)
Functions.sorted_stacked_barplot(DATA, 15)
Functions.sorted_stacked_barplot(DATA, 100)
# F ------------------------ PROPORTIONS BARPLOT: DISTRIBUTION PER GENE ------------------------

# I ------------------------ NUMBER OF ISOFORMS EXPRESSED PER GENE ------------------------
Functions.number_of_isoforms(DATA, 0)
Functions.number_of_isoforms(DATA, 0, thres2 = 0.9)
Functions.number_of_isoforms(DATA, 1)
Functions.number_of_isoforms(DATA, 15)
Functions.number_of_isoforms(DATA, 100)
# F ------------------------ NUMBER OF ISOFORMS EXPRESSED PER GENE ------------------------

# I ------------------------ CLASSIFYING EACH GENE ------------------------
Functions.type_of_gene(DATA, 0)
Functions.type_of_gene(DATA, 0, thres2 = 0.9)
Functions.type_of_gene(DATA, 1)
Functions.type_of_gene(DATA, 15)
Functions.type_of_gene(DATA, 100)
# F ------------------------ CLASSIFYING EACH GENE ------------------------

# I ------------------------ BIG SUMMARY ------------------------
Functions.big_summary(DATA, 4)

df = Functions.big_summary(DATA).T
df = df.astype('category')

# Statistics
Functions.statistics(df, plot = True)
# F ------------------------ BIG SUMMARY ------------------------

# I ------------------------ TEC ------------------------
inters = Functions.getting_interactions(DATA)

Functions.TEC_ij(inters, 0)
Functions.TEC_ij(inters, 1)
Functions.TEC_ij(inters, 10)
Functions.TEC_ij(inters, 100)
Functions.TEC_ij(inters, 2000)

Functions.TEClevels_ij(inters, 0)
Functions.TEClevels_ij(inters, 1)
Functions.TEClevels_ij(inters, 10)
Functions.TEClevels_ij(inters, 100)
Functions.TEClevels_ij(inters, 2000)
# F ------------------------ TEC ------------------------