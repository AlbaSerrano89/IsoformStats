#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 19 12:06:44 2018

@author: aserrano
"""
import os
import numpy
import matplotlib.pyplot as plt
import pandas as pd

def reading_data(csv_file):
    lastbar = csv_file.rfind('/')
    file = csv_file[(lastbar + 1):]
    
    ni = 0
    nf = len(file)
    
    if file.startswith('SMTS_'):
        ni = 5
    elif file.startswith('SMTSD_'):
        ni = 6
    
    if file.endswith('gz'):
        nf = -7
        data = pd.read_csv(csv_file, compression = 'gzip', header = 0, sep = '\t')
    elif file.endswith('csv'):
        nf = -4
        data = pd.read_csv(csv_file, header = 0, sep = '\t')
    
    tissue = file[ni:nf]
    genes = list(data['gene_name'].unique())
    
    return([genes, data, tissue])

def formatting_gene(all_data, gene):
    genes = all_data[0]
    
    if type(gene) == int:
        gene2 = genes[gene]
    elif type(gene) == str or type(gene) == numpy.str_:
        if gene in genes:
            gene2 = gene
        else:
            return('The gene ' + gene + ' is not in your data.')
    else:
        return('Incorrect type of element for the gene: ' + str(type(gene)))
    
    return(gene2)

def gene_info(all_data, gene, all_samps = 'F'):
    gene2 = formatting_gene(all_data, gene)
    genes = all_data[0]
        
    if gene2 in genes:
        df = all_data[1]
        gene_df = df.loc[df['gene_name'] == gene2].drop(['gene_name', 'gene_id'], axis = 1).set_index('transcript_id')
        
        if all_samps == 'F':
            columns = list(gene_df)[2:]
            
            for col in columns:
                if gene_df[col][0] == -1:
                    gene_df = gene_df.drop([col], axis = 1)
        
        return(gene_df)
    else:
        return(gene2)

def gene_statistics(all_data, gene):
    gene2 = formatting_gene(all_data, gene)
    genes = all_data[0]
        
    if gene2 in genes:
        geneinfo = gene_info(all_data, gene2)
        
        if len(list(geneinfo)) == 2:
            return('The gene ' + gene2 + ' is not expressed in any of the analysed samples.')
        else:
            gene_stats = geneinfo.T[2:].astype('float64').describe().T.sort_values(by = 'mean', ascending = False)
            return(gene_stats)
    else:
        return(gene2)

def gene_proportions(all_data, gene):
    gene2 = formatting_gene(all_data, gene)
    genes = all_data[0]
    
    if gene2 in genes:
        gene_stats = gene_statistics(all_data, gene2)
        
        if 'not expressed' in gene_stats:
            return(gene_stats)
        else:
            sorted_mean = gene_stats['mean']
            
            acum_props = []
            acum = 0
            for iso in list(sorted_mean):
                acum = acum + iso
                if acum < 1:
                    acum_props.append(acum)
                else:
                    acum_props.append(1.0)
            
            df = pd.DataFrame({'mean' : sorted_mean, 'cumulation' : acum_props})
            
            return(df)
    else:
        return(gene2)

def gene_filtered_proportions(all_data, gene, thres = 0.1):
    gene2 = formatting_gene(all_data, gene)
    genes = all_data[0]
        
    if gene2 in genes:
        gene_props = gene_proportions(all_data, gene)
        
        if 'not expressed' in gene_props:
            return(gene_props)
        else:
            means = gene_props['mean']
            cumul = gene_props['cumulation']
            
            high_props = []
            cumu_props = []
            for i in range(0, len(means)):
                if means[i] >= thres:
                    high_props.append(means[i])
                    cumu_props.append(cumul[i])
                else:
                    if i == 0:
                        return('There is not any isoform with a mean expression higher than ' + str(thres))
            
            nhigh = len(high_props)
            high_props.append(sum(means[nhigh:]))
            
            columns = list(gene_props[:nhigh].index)
            columns.append('other')
            
            if round(high_props[-1] + cumu_props[-1], 5) == 1:
                cumu_props.append(1.0)
            else:
                cumu_props.append(high_props[-1] + cumu_props[-1])
            
            df = pd.DataFrame({'mean' : high_props, 'cumulation' : cumu_props}, index = columns)
            
            return(df)
    else:
        return(gene2)

def gene_filtered_statistics(all_data, gene, thres = 0.1):
    gene2 = formatting_gene(all_data, gene)
    genes = all_data[0]
    
    if gene2 in genes:
        gene_props = gene_filtered_proportions(all_data, gene, thres)
        
        if 'not expressed' in gene_props or 'higher than' in gene_props:
            return(gene_props)
        else:
            gene_props = gene_props.drop(['other'])
            total_sum = gene_props['cumulation'][-1]
            
            new_means = list(gene_props['mean'])/total_sum
            new_cumus = list(gene_props['cumulation'])/total_sum
            
            gene_props['new_means'] = new_means
            gene_props['new_cumulation'] = new_cumus

            return(gene_props)
    else:
        return(gene2)

def gene_classification(all_data, gene, thres = 0.1, thres2 = 0.7, thres3 = 10):
    gene2 = formatting_gene(all_data, gene)
    genes = all_data[0]
    
    if gene2 in genes:
        isos_stats = gene_filtered_statistics(all_data, gene, thres)
        
        geneinfo = gene_info(all_data, gene2)
        colnames = list(gene_info(all_data, gene2))
        genetype = geneinfo[colnames[0]][0]
        
        samp_names = colnames[2:]
        nsamps = len(samp_names)
        
        Types = {}
        if 'not expressed' in isos_stats or 'higher than' in isos_stats or nsamps < thres3:
            isos = transtype = ini_props = acum_props = new_ini_props = new_acum_props = '-'
            
            if 'not expressed' in isos_stats:
                cog = 'NotExpressed'
            elif 'higher than' in isos_stats:
                cog = 'LowExpressedTranscripts'
            else:
                cog = 'FewSamples'
        else:
            cols = ['new_cumulation']
            nisos = len(isos_stats[isos_stats[cols] <= thres2][cols].dropna()) + 1
            
            if nisos == 1:
                cog = 'Monoform'
            elif nisos == 2:
                cog = 'Biform'
            elif nisos == 3:
                cog = 'Triform'
            else:
                cog = 'Multiform'
            
            new_isos = isos_stats[:nisos]
            
            isos = list(new_isos.index)
            ini_props = list(new_isos['mean'])
            acum_props = list(new_isos['cumulation'])
            new_ini_props = list(new_isos['new_means'])
            new_acum_props = list(new_isos['new_cumulation'])
            
            transtype = []
            for iso in isos:
                transtype.append(geneinfo[colnames[1]][iso])
                
        Types[gene2] = [genetype, cog, isos, transtype, ini_props, acum_props, new_ini_props, new_acum_props]
        return(Types)
    else:
        return(gene2)

def big_summary(all_data, out_bsdir, dfbsfile = 'T', out_bsfile = '_genes_', thres = 0.1, thres2 = 0.7, thres3 = 10):
    tissue = all_data[2]
    
    all_types = {}
    for i in range(0, len(all_data[0])):
        res = gene_classification(all_data, i, thres, thres2, thres3)
        all_types.update(res)
    
    df = pd.DataFrame(all_types, dtype = 'category')
    df = df.T
    
    df.columns = ['Gene', 'GeneType', 'GeneClassification', 'Transcripts', 'TranscriptTypes', 'InitialProportion', 'CumulativeProportion', 'NewProportion', 'NewCumulativeProportion']
    
    if dfbsfile == 'T':
        if not os.path.isdir(out_bsdir):
            os.mkdir(out_bsdir)
        
        if out_bsdir[-1] != '/':
            out_bsdir = out_bsdir + '/'
        
        if out_bsfile != '_genes_':
            bstissuefile = out_bsdir + out_bsfile
        else:
            bstissuefile = out_bsdir + tissue + out_bsfile + str(thres * 100) + '_' + str(thres2 * 100) + '_' + str(thres3) + '.csv'
        
        df.to_csv(bstissuefile)
    elif dfbsfile == 'F':
        return(df)
    else:
        print(dfbsfile + ' is not an acceptable option for --dfstatsfile.')

def statistics(all_data, in_bsfile, freq_type = 'rel', tissuestatsfile = '_statistics.csv', dfstatsfile = 'F'):
    df = pd.read_csv(in_bsfile, index_col = 0)
    
    df['GeneType'] = df.GeneType.astype('category')
    df['GeneClassification'] = df.GeneClassification.astype('category')
    df['TranscriptTypes'] = df.TranscriptTypes.astype('category')
    
    numbs = df.GeneClassification.value_counts()
    
    df2 = pd.DataFrame({'TotalCounts' : numbs})
    df3 = pd.crosstab(df.GeneClassification, df.GeneType)
    
    df4 = df2.join(df3)
    
    if freq_type == 'abs':
        final_df = df4
    elif freq_type == 'rel':
        final_df = df4 / sum(numbs)
    
    if dfstatsfile == 'T':
        if tissuestatsfile != '_statistics.csv':
            tissuestatsfile = tissuestatsfile
        else:
            tissue = all_data[2]
            tissuestatsfile = tissue + tissuestatsfile
        
        final_df.to_csv(tissuestatsfile)
        
    elif dfstatsfile == 'F':
        return(final_df)
    else:
        print(dfstatsfile + ' is not an acceptable option for --dfstatsfile.')

def all_tissues_analysis(original_dir, pref = '', new_dir = '', thres = 0.1, thres2 = 0.7, thres3 = 10):
    if pref != '' and pref[-1] != '_':
        pref = pref + '_'
    
    if new_dir == '':
        if original_dir[-1] == '/':
            bsumm_dir = original_dir + pref + 'big_summaries'
        else:
            bsumm_dir = original_dir + '/' + pref + 'big_summaries'
    else:
        bsumm_dir = new_dir
    
    if not os.path.isdir(bsumm_dir):
        os.mkdir(bsumm_dir)
    
    os.chdir(original_dir)
    
    nfiles = 0
    for file in os.listdir():
        if file.startswith(pref) and file.endswith('gz') and 'all' not in file:
            nfiles += 1
    
    print('There are ' + str(nfiles) + ' tissues to analyse.\n')
    
    i = 0
    for file in os.listdir():
        if 'all' not in file.lower():
            if file.endswith('gz'):
                nf = -7
            elif file.endswith('csv'):
                nf = -4
            
            if pref != '':
                if pref == 'SMTS_':
                    ni = 5
                elif pref == 'SMTSD_':
                    ni = 6
                
                if file.startswith(pref) and file.endswith('gz'):
                    tissue = file[ni:nf]
                    
                    if original_dir[-1] == '/':
                        tissuepath = original_dir + file
                    else:
                        tissuepath = original_dir + '/' + file
                    
                    print(tissuepath)
                    DATA = reading_data(tissuepath)
                    big_summary(DATA, bsumm_dir, dfbsfile = 'T', out_bsfile = 'BS_' + tissue + '.csv', thres = thres, thres2 = thres2, thres3 = thres3)
                    
                    i += 1
                    print(str(i) + '.- Big summary of the tissue ' + tissue + ': done.')
            else:
                if file.endswith('gz'):
                    tissue = file[:nf]
                    
                    if original_dir[-1] == '/':
                        tissuepath = original_dir + file
                    else:
                        tissuepath = original_dir + '/' + file
                    
                    print(tissuepath)
                    DATA = reading_data(tissuepath)
                    big_summary(DATA, bsumm_dir, dfbsfile = 'T', out_bsfile = 'BS_' + tissue + '.csv', thres = thres, thres2 = thres2, thres3 = thres3)
                    
                    i += 1
                    print(str(i) + '.- Big summary of the tissue ' + tissue + ': done.')















# ----------------------------------------------- PLOTS ----------------------------------------------- #
def gene_general_view(all_data, gene):
    gene2 = formatting_gene(all_data, gene)
    genes = all_data[0]
    
    if gene2 in genes:
        gene_props = gene_info(all_data, gene)
        
        if 'not expressed' in gene_props:
            return(gene_props)
        else:
            ax = gene_props.T[2:].plot.bar(stacked = True, fontsize = 8)
            ax.set_title('Distribution of samples of the gene ' + gene2)
            ax.set_xlabel("Samples")
            ax.set_ylabel("Proportion")
            plt.subplots_adjust(bottom = 0.22, top = 0.95)
            plt.show()
    else:
        return(gene2)

def stacked_barplot(all_data, gene, thres = 0.1):
    gene2 = formatting_gene(all_data, gene)
    genes = all_data[0]
    
    if gene2 in genes:
        gene_props = gene_filtered_proportions(all_data, gene, thres)
        
        if 'not expressed' in gene_props:
            return(gene_props)
        else:
            selected_isos = list(gene_props.index)
            selected_isos.pop()
            geneinfo = gene_info(all_data, gene).drop(['gene_type', 'transcript_type'], axis = 1)
            
            df = geneinfo.T[selected_isos]
            other_df = geneinfo.T.drop(selected_isos, axis = 1)
            df['other'] = other_df.sum(axis = 1)
            
            ax = df.plot.bar(stacked = True, fontsize = 8)
            ax.set_title('Distribution of samples of the gene ' + gene2)
            ax.set_xlabel("Samples")
            ax.set_ylabel("Proportion")
            plt.subplots_adjust(bottom = 0.22, top = 0.95)
            plt.show()
    else:
        return(gene2)
