#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 19 12:06:44 2018

@author: aserrano
"""
import os
import gzip
import statistics
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def reading_data(csv_file):
    Data = dict()
    
    lastbar = csv_file.rfind('/')
    file = csv_file[(lastbar + 1):]
    
    if file.startswith('SMTS_'):
        ni = 5
    elif file.startswith('SMTSD_'):
        ni = 6
    
    if file.endswith('gz'):
        nf = -7
    elif file.endswith('csv'):
        nf = -4
    
    tissue = file[ni:nf]
    
    with gzip.open(csv_file, 'rt') if csv_file.endswith('gz') else open(csv_file) as rd:
        nsamps = 0
        samp_names = ''
        
        for line in rd:
            ff = line.strip().split('\t')
            if line.startswith('transc'):
                columns = tuple(ff)
                samp_names = columns[5:]
                nsamps = len(samp_names)
                
            else:
                # here I parse all lines to generate the Data dictionary
                # keys = 'gene_names'
                # Data[gene] = list of isos
                # each iso  = [ENST, samp]
                gene_name = ff[1]
                gene_type = ff[2]
                nums  = list()
                samps = list()
                for i in range(0, nsamps):
                    num = float(ff[i + 5])
                    if num >= 0:
                        samps.append(samp_names[i])
                        nums.append(num)
                # here I update or generate the data[gene_id] value
                tmp = Data.get(gene_name, None)
                if tmp is None:
                    Data[gene_name] = [gene_type, tuple(samps), [ff[0], ff[3], tuple(nums)]]
                else:
                    tmp.append([ff[0], ff[3], tuple(nums)])
                    Data[gene_name] = tmp

    return([list(Data.keys()), Data, tissue, samp_names])

def gene_name(all_data, gene):
    genes = all_data[0]
    
    if type(gene) == int:
        gene2 = genes[gene]
    elif type(gene) == (str or np.str_):
        if gene in genes:
            gene2 = gene
        else:
            return('The gene ' + gene + ' is not in your data.')
    else:
        return('Incorrect type of element for the gene: ' + str(type(gene)))
    
    return(gene2)

def gene_info(all_data, gene):
    gene2 = gene_name(all_data, gene)
    genes = all_data[0]
        
    if gene2 in genes:
        geneinfo1 = all_data[1][gene2]
        
        nsamps = len(geneinfo1[1])
        for trans in geneinfo1[2:]:
            if len(trans[2]) == 0:
                trans[2] = tuple(np.repeat(0.0, nsamps))
        
        geneinfo = {}
        geneinfo[gene2] = geneinfo1
        return(geneinfo)
    else:
        return(gene2)

def gene_data(all_data, gene):
    gene2 = gene_name(all_data, gene)
    genes = all_data[0]
        
    if gene2 in genes:
        geneinfo = gene_info(all_data, gene2)[gene2]
        
        gene_type = geneinfo[0]
        samp_names = geneinfo[1]
        count = len(samp_names)
        
        genedata = {}
        genedata[gene2] = [gene_type, samp_names, count]
        return(genedata)
    else:
        return(gene2)

def gene_statistics(all_data, gene, pretty = 'F'):
    gene2 = gene_name(all_data, gene)
    genes = all_data[0]
    
    if gene2 in genes:
        geneinfo = gene_info(all_data, gene2)[gene2]
        genedata = gene_data(all_data, gene2)[gene2]
        
        if genedata[2] == 0:
            return('The gene ' + gene2 + ' is not expressed in any of the analysed samples.')
        else:
            def takeMean(elem):
                return elem[2]
            
            isoforms = []
            list_stats = []
            for trans in geneinfo[2:]:
                iso = trans[0]
                isoforms.append(iso)
                isotype = trans[1]
                mn = np.mean(trans[2])
                
                if len(trans[2]) == 1:
                    sd = 0.0
                else:
                    sd = statistics.stdev(trans[2])
                Q1 = np.percentile(trans[2], 25)
                med = np.median(trans[2])
                Q3 = np.percentile(trans[2], 75)
                minim = min(trans[2])
                maxim = max(trans[2])
                
                list_stats.append([iso, isotype, mn, sd, Q1, med, Q3, minim, maxim])
            
            list_stats = sorted(list_stats, key = takeMean, reverse = True)
            
            for i in range(0, len(list_stats)):
                if i == 0:
                    new = list_stats[i][2]
                else:
                    new = list_stats[i][2] + list_stats[i - 1][3]
                
                if round(new, 5) == 1:
                    new = 1.0
                
                list_stats[i].insert(3, new)
            
            gene_stats = {}
            gene_stats[gene2] = genedata + list_stats
            
            if pretty == 'T':
                print('\n\n\t' + gene2 + ' --> ' + genedata[0] + '\n\n')
                pretty_df = gene_stats[gene2][3:]
                df = pd.DataFrame(pretty_df, index = isoforms, columns = ['Transcript', 'TranscriptType', 'Mean', 'CumulativeMean', 'StandardDeviation', 'Q1', 'Median', 'Q3', 'Minimum', 'Maximum'])
                df2 = df.drop(['Transcript'], axis = 1)
                return(df2)
            else:
                return(gene_stats)
    else:
        return(gene2)

def gene_filtered_proportions(all_data, gene, pretty = 'F', minexp = 0.1):
    gene2 = gene_name(all_data, gene)
    genes = all_data[0]
        
    if gene2 in genes:
        gene_stats = gene_statistics(all_data, gene2, 'F')
        
        if 'not expressed' in gene_stats:
            return(gene_stats)
        else:
            genestats = gene_stats[gene2]
            genedata = gene_data(all_data, gene2)[gene2]
            
            means = [trans[2] for trans in genestats[3:]]
            means_array = np.asarray(means)
            high_props = list(means_array[means_array > minexp])
            ntrans = len(high_props)
            
            if ntrans == 0:
                return('There is not any isoform of the gene ' + gene2 + ' with a mean expression higher than ' + str(minexp))
            else:
                transcripts = [trans[0] for trans in genestats[3:][:ntrans]]
                trans_type = [trans[1] for trans in genestats[3:][:ntrans]]
                cumu_props = [trans[3] for trans in genestats[3:][:ntrans]]
                
                total_sum = np.float64(cumu_props[-1])
                
                new_means = list(means[:ntrans] / total_sum)
                new_cumus = list(cumu_props / total_sum)
                
                if ntrans != len(genestats[3:]):
                    transcripts.append('other')
                    trans_type.append('-')
                    high_props.append(sum(means[ntrans:]))
                
                    last_cum = high_props[-1] + cumu_props[-1]
                    if round(last_cum, 5) == 1:
                        cumu_props.append(1.0)
                    else:
                        print('The gene ' + gene2 + ' is just explained in a ' + str(round(last_cum + 100, 2)) + '%.')
                        cumu_props.append(last_cum)
                    
                    new_means.append('-')
                    new_cumus.append('-')
            
            results = [transcripts, trans_type, high_props, cumu_props, new_means, new_cumus]
            
            fin_res = []
            for i in range(0, len(results[0])):
                fin_res.append([lst[i] for lst in results])
                
            gene_filt_props = {}
            gene_filt_props[gene2] = genedata + fin_res
            
            if pretty == 'T':
                print('\n\n\t' + gene2 + ' --> ' + genedata[0] + '\n\n')
                pretty_df = gene_filt_props[gene2][3:]
                df = pd.DataFrame(pretty_df, index = results[0], columns = ['Transcript', 'TranscriptType', 'Mean', 'CumulativeMean', 'NewMean', 'NewCumulativeMean'])
                df2 = df.drop(['Transcript'], axis = 1)
                return(df2)
            else:
                return(gene_filt_props)
    else:
        return(gene2)

def gene_classification(all_data, gene, pretty = 'F', minexp = 0.1, mintotexp = 0.7, minsamps = 10):
    gene2 = gene_name(all_data, gene)
    genes = all_data[0]
    
    if gene2 in genes:
        gene_filt_props = gene_filtered_proportions(all_data, gene, 'F', minexp)
        genedata = gene_data(all_data, gene2)[gene2]
        
        if type(gene_filt_props) == str:
            transcripts = transcript_types = means = cumu_means = new_means = new_cumu_means = '-'
            
            if 'not expressed' in gene_filt_props:
                cog = 'NotExpressed'
            else:
                cog = 'LowExpressedTranscripts'
        else:
            if genedata[2] < minsamps:
                transcripts = transcript_types = means = cumu_means = new_means = new_cumu_means = '-'
                cog = 'FewSamples'
            else:
                gene_filt = gene_filt_props[gene2][3:]
                gene_filt.pop()
                
                if len(gene_filt) == 1:
                    cog = 'Monoform'
                elif len(gene_filt) == 2:
                    cog = 'Biform'
                elif len(gene_filt) == 3:
                    cog = 'Triform'
                else:
                    cog = 'Multiform'
                
                transcripts = [trans[0] for trans in gene_filt]
                transcript_types = [trans[1] for trans in gene_filt]
                means = [trans[2] for trans in gene_filt]
                cumu_means = [trans[3] for trans in gene_filt]
                new_means = [trans[4] for trans in gene_filt]
                new_cumu_means = [trans[5] for trans in gene_filt]
                
        Types = {}
        Types[gene2] = [genedata[0], cog, transcripts, transcript_types, means, cumu_means, new_means, new_cumu_means]
        
        if pretty == 'T':
            print('\n\n\t' + gene2 + ' --> ' + genedata[0] + '\n\n')
            df = pd.DataFrame(gene_filt, index = transcripts, columns = ['Transcript', 'TranscriptType', 'Mean', 'CumulativeMean', 'NewMean', 'NewCumulativeMean'])
            df2 = df.drop(['Transcript'], axis = 1)
            return(df2)
        else:
            return(Types)
    else:
        return(gene2)

def big_summary(all_data, out_bsdir, out_bsfile = '_genes_', minexp = 0.1, mintotexp = 0.7, minsamps = 10):
    tissue = all_data[2]
    
    all_types = {}
    for i in range(0, len(all_data[0])):
        res = gene_classification(all_data, i, 'F', minexp, mintotexp, minsamps)
        all_types.update(res)
    
    df = pd.DataFrame(all_types, dtype = 'category')
    df = df.T
    
    df.columns = ['GeneType', 'GeneClassification', 'Transcripts', 'TranscriptTypes', 'Mean', 'CumulativeMean', 'NewMean', 'NewCumulativeMean']
    
    if out_bsdir[-1] != '/':
        out_bsdir = out_bsdir + '/'
    
    if out_bsfile != '_genes_':
        bstissuefile = out_bsdir + out_bsfile
    else:
        bstissuefile = out_bsdir + tissue + out_bsfile + str(minexp * 100) + '_' + str(mintotexp * 100) + '_' + str(minsamps) + '.csv'
    
    if not os.path.isdir(out_bsdir):
        os.mkdir(out_bsdir)
    
    df.to_csv(bstissuefile)

def tissue_statistics(all_data, in_bsfile, freq_type = 'rel', tissuestatsfile = '_statistics.csv', dfstatsfile = 'F'):
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

def all_tissues_analysis(original_dir, pref = '', new_dir = '', minexp = 0.1, mintotexp = 0.7, minsamps = 10):
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
                    big_summary(DATA, bsumm_dir, dfbsfile = 'T', out_bsfile = 'BS_' + tissue + '.csv', minexp = minexp, mintotexp = mintotexp, minsamps = minsamps)
                    
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
                    big_summary(DATA, bsumm_dir, dfbsfile = 'T', out_bsfile = 'BS_' + tissue + '.csv', minexp = minexp, mintotexp = mintotexp, minsamps = minsamps)
                    
                    i += 1
                    print(str(i) + '.- Big summary of the tissue ' + tissue + ': done.')

def all_tissues_stats(bigsummaries_dir, allstats_file = 'all_tissues_statistics.csv'):
    types = ['NotExpressed', 'FewSamples', 'Monoform', 'Biform', 'Triform', 'Multiform']
    
    stats = {}
    for file in os.listdir(bigsummaries_dir):
        if file != 'Statistics':
            df = pd.read_csv(bigsummaries_dir + file, engine='python')
            tissue = file[3:-4]
            
            nnexp = df.loc[df['Type'] == 'NotExpressed']['Type'].count()
            nfew = df.loc[df['Type'] == 'FewSamples']['Type'].count()
            nmono = df.loc[df['Type'] == 'Monoform']['Type'].count()
            nbi = df.loc[df['Type'] == 'Biform']['Type'].count()
            ntri = df.loc[df['Type'] == 'Triform']['Type'].count()
            nmulti = df.loc[df['Type'] == 'Multiform']['Type'].count()
            numbs = [nnexp, nfew, nmono, nbi, ntri, nmulti]
            
            stats[tissue] = numbs
    
    df = pd.DataFrame(stats, index = types)
    df = df.T
    df.to_csv(allstats_file)















# ----------------------------------------------- PLOTS ----------------------------------------------- #
def visual_matrix(all_data, gene):
    gene2 = gene_name(all_data, gene)
    genes = all_data[0]
    
    if gene2 in genes:
        gene_props = gene_info(all_data, gene, 'T')
        
        samp_names = list(gene_props)[2:]
        nsamps = len(samp_names)
        indx = np.arange(nsamps)
        
        isoforms = list(gene_props.index)
        nisos = len(isoforms)
        
        indy = np.arange(nisos) + 0.5
        
        fig, ax = plt.subplots()
        width = 1
        
        acum_vect = [0 for i in range(0, nsamps)]
        sum_vect = [1 for i in range(0, nsamps)]
            
        if 'not expressed' in gene_props:
            for iso in isoforms:
                ax.bar(indx, sum_vect, width, bottom = acum_vect, color = 'white', edgecolor = 'gray')
                acum_vect = tuple(sum(x) for x in zip(acum_vect, sum_vect))
        else:
            gene_props2 = gene_props.T
            for iso in isoforms:
                colors = []
                for numb in gene_props2[iso][2:]:
                    if numb > 0:
                        colors.append('black')
                    else:
                        colors.append('white')
                    
                ax.bar(indx, sum_vect, width, bottom = acum_vect, color = colors, edgecolor = 'gray')
                acum_vect = tuple(sum(x) for x in zip(acum_vect, sum_vect))
        
        ax.set_xlabel('Samples')
        ax.set_ylabel('Isoforms')
        ax.set_title('Gene ' + gene2)
        ax.set_xticks(indx)
        ax.set_xticklabels(samp_names, rotation = 'vertical', fontsize = 8)
        ax.set_yticks(indy)
        ax.set_yticklabels(isoforms, fontsize = 8)
        plt.subplots_adjust(bottom = 0.22, top = 0.95, left = 0.2)
        plt.show()
    else:
        return(gene2)

def gene_boxplot(all_data, gene, ordered = True):
    gene2 = gene_name(all_data, gene)
    genes = all_data[0]
    
    if gene2 in genes:
        gene_props = gene_info(all_data, gene)
        
        if 'not expressed' in gene_props:
            return(gene_props)
        else:
            samp_names = list(gene_props)[2:]
            gene_props = gene_props[samp_names]
            
            if ordered:
                gene_props['median'] = gene_props.median(axis = 1)
                df = gene_props.sort_values(by = 'median', ascending = False).T
                df2 = df.drop(['median'])
                ax = df2.plot.box(rot = 90, fontsize = 8)
            else:
                ax = gene_props.plot.box(rot = 90, fontsize = 8)
            
            ind = list(np.arange(0, 1.1, 0.1))
            ax.set_yticks(ind)
            ax.set_title('Distribution of samples of the gene ' + gene2)
            axes = plt.gca()
            axes.set_ylim([-0.05, 1])
            ax.set_xlabel("Isoforms")
            plt.subplots_adjust(bottom = 0.22, top = 0.95)
            plt.show()
    else:
        return(gene2)

def gene_general_view(all_data, gene):
    gene2 = gene_name(all_data, gene)
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

def stacked_barplot(all_data, gene, minexp = 0.1):
    gene2 = gene_name(all_data, gene)
    genes = all_data[0]
    
    if gene2 in genes:
        gene_props = gene_filtered_proportions(all_data, gene, minexp)
        
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

def stats_barplot(all_data, csv_file):
    stats = statistics(all_data, csv_file, dffile = 'F')
    tissue = all_data[2]
    
    types = list(stats.index)
    props = list(stats.Typeprop)
    
    numbs = []
    for count in list(stats.Counts):
        numbs.append(int(count))
    
    labs = []
    for i in range(0, len(types)):
        labs.append(types[i] + ': ' + str(numbs[i]))
    
    ntypes = len(types)
    colors = ['C0', 'C9', 'C1', 'C2', 'C3', 'C5']
    
    fig, ax = plt.subplots()
    
    acum = props[0]
    ax.bar(1, props[0], 1, label = labs[0], color = colors[0])
    for i in range(1, ntypes):
        ax.bar(1, props[i], 1, bottom = acum, label = labs[i], color = colors[i])
        acum =  acum + props[i]
    
    ax.set_ylabel('Proportion')
    ax.set_title('Isoform distribution of ' + tissue)
    ax.set_xticks([0])
    ax.legend()
    plt.subplots_adjust(bottom = 0.1, top = 0.9)
    plt.show()

def expr_barplot(all_data, csv_file):
    stats = statistics(all_data, csv_file, dffile = 'F')
    tissue = all_data[2]
    
    types = list(stats.index)[2:]
    props = list(stats.Typeprop)[2:]
    
    numbs = []
    for count in list(stats.Counts)[2:]:
        numbs.append(int(count))
    
    total = sum(numbs)
    
    props = []
    for numb in numbs:
        props.append(numb / float(total))
    
    labs = []
    for i in range(0, len(types)):
        labs.append(types[i] + ': ' + str(numbs[i]))
    
    ntypes = len(types)
    colors = ['C1', 'C2', 'C3', 'C5']
    
    fig, ax = plt.subplots()
    
    acum = props[0]
    ax.bar(1, props[0], 1, label = labs[0], color = colors[0])
    for i in range(1, ntypes):
        ax.bar(1, props[i], 1, bottom = acum, label = labs[i], color = colors[i])
        acum =  acum + props[i]
    
    ax.set_ylabel('Proportion')
    ax.set_title('Isoform distribution of ' + tissue)
    ax.set_xticks([0])
    ax.legend()
    plt.subplots_adjust(bottom = 0.1, top = 0.9)
    plt.show()

def all_tissues_barplot(stats_filename = 'all_tissues_statistics.csv'):
    df1 = pd.read_csv(stats_filename, index_col = 0)
    df2 = df1.T
    
    types = list(df1.columns)
    tissues = list(df2.columns)
    
    total_col = {}
    for tissue in tissues:
        total_col[tissue] = sum(df2[tissue])
    
    props_dict = {}
    for tissue in tissues:
        props_tiss = []
        for numbs in list(df2[tissue]):
            props_tiss.append(numbs / total_col[tissue])
        
        props_dict[tissue] = props_tiss
    
    props_df1 = pd.DataFrame(props_dict, index = types)
    props_df2 = props_df1.T
    
    fig, ax = plt.subplots()
    ind = np.arange(len(tissues))
    width = 0.5
    colors = ['C0', 'C9', 'C1', 'C2', 'C3', 'C5']
    
    props_vect = []
    for typegene in types:
        props_vect.append(list(props_df2[typegene]))
    
    acum_vect = props_vect[0]

    ax.barh(ind, props_vect[0], width, label = types[0], color = colors[0])
    for i in range(1, len(types)):
        ax.barh(ind, props_vect[i], width, left = acum_vect, label = types[i], color = colors[i])
        acum_vect = tuple(sum(x) for x in zip(acum_vect, props_vect[i]))
    
    ax.set_xlabel('Proportion')
    ax.set_ylabel('Tissues')
    ax.set_title('Counts of the analysis')
    ax.set_yticks(ind)
    ax.set_yticklabels(tissues, fontsize = 8)
    ax.legend()
    plt.subplots_adjust(bottom = 0.2, top = 0.9)
    plt.show()















# ----------------------------------------------- TEC analysis ----------------------------------------------- #

def getting_interactions(all_data):
    data_interactions = {}
    for gene, dat in all_data[1].items():
        iso_gene = []
        for i in range(0, len(dat) - 1):
            for j in range(i + 1, len(dat)):
                iso_gene.append([(dat[i][0], dat[j][0]), [dat[i][1], dat[j][1]]])
    
        data_interactions[gene] = iso_gene
    
    return(data_interactions)

def TEC_ij(ints, gene):
    genes = list(ints.keys())
    if type(gene) == int:
        gene2 = genes[gene]
    elif type(gene) == str or type(gene) == np.str_:
        if gene in genes:
            gene2 = gene
        else:
            print('The gene ' + gene2 + ' is not in your data.')
    else:
        print('Incorrect type of element for the gene: ' + str(type(gene)))
    
    TEC = {}
    TEC_list = []
    for inter in ints[gene2]:
        isoform1 = list(inter[1][0])
        isoform2 = list(inter[1][1])
        
        d1 = d2 = t1 = t2 = 0
        for i in range(0, len(isoform1)):
            if isoform1[i] == 0 and isoform2[i] == 0:
                d1 += 0
                d2 += 0
            elif isoform1[i] == 0:
                d2 += 1
                t2 += 1
            elif isoform2[i] == 0:
                d1 += 1
                t1 += 1
            else:
                t1 += 1
                t2 += 1
        
        if t1 != 0:
            m1 = d1 / float(t1)
        else:
            m1 = 0
    
        if t2 != 0:
            m2 = d2 / float(t2)
        else:
            m2 = 0
    
        TECij = round((m1 + m2) / 2.0, 3)
        TEC_list.append([inter[0], [(d1, t1, d2, t2), TECij]])
    
    TEC[gene2] = TEC_list
    return(TEC)

def TEClevels_ij(ints, gene):
    genes = list(ints.keys())
    if type(gene) == int:
        gene2 = genes[gene]
    elif type(gene) == str or type(gene) == np.str_:
        if gene in genes:
            gene2 = gene
        else:
            print('The gene ' + gene2 + ' is not in your data.')
    else:
        print('Incorrect type of element for the gene: ' + str(type(gene)))
    
    TEC = {}
    TEC_list = []
    for inter in ints[gene2]:
        acum1 = 0
        acum2 = 0
        sum1 = 0
        sum2 = 0
        for i in range(0, len(inter[1][0])):
            sum1 = sum1 + inter[1][0][i]
            sum2 = sum2 + inter[1][1][i]
            
            if inter[1][0][i] > inter[1][1][i]:
                acum1 = acum1 + (inter[1][0][i] - inter[1][1][i])
            else:
                acum2 = acum2 + (inter[1][1][i] - inter[1][0][i])
        
        if sum1 == 0 and sum2 == 0:
            tec_ij = 0
        elif sum1 == 0:
            tec_ij = (acum2 / float(sum2)) / 2.0
        elif sum2 == 0:
            tec_ij = (acum1 / float(sum1)) / 2.0
        else:
            tec_ij = (acum1 / float(sum1) + acum2 / float(sum2)) / 2.0
        
        TECij = round(tec_ij, 3)
        TEC_list.append([inter[0], TECij])
    
    TEC[gene2] = TEC_list
    return(TEC)
