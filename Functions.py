#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 19 12:06:44 2018
@author: aserrano
"""
import os
import gzip
import glob
import re
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
#import cProfile

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
                genename = ff[1]
                genetype = ff[2]
                geneid = ff[4]
                nums  = list()
                samps = list()
                for i in range(0, nsamps):
                    num = float(ff[i + 5])
                    if num >= 0:
                        samps.append(samp_names[i])
                        nums.append(num)
                # here I update or generate the data[gene_id] value
                tmp = Data.get(geneid, None)
                if tmp is None:
                    Data[geneid] = [genename, genetype, tuple(samps), [ff[0], ff[3], tuple(nums)]]
                else:
                    tmp.append([ff[0], ff[3], tuple(nums)])
                    Data[geneid] = tmp

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

def gene_pos(all_data, gene):
    genes = all_data[0]
    
    if type(gene) == int:
        gene2 = gene
    elif type(gene) == (str or np.str_):
        if gene in genes:
            gene2 = genes.index(gene)
        else:
            return('The gene ' + gene + ' is not in your data.')
    else:
        return('Incorrect type of element for the gene: ' + str(type(gene)))
    
    return(gene2)

def gene_info(all_data, gene, pretty = 'F'):
    gene2 = gene_name(all_data, gene)
    genes = all_data[0]
        
    if gene2 in genes:
        geneinfo = {}
        geneinfo[gene2] = all_data[1][gene2]
        
        if pretty == 'T':
            print('\n' + gene2 + ' - ' + geneinfo[gene2][0] + ' ---> ' + geneinfo[gene2][1] + '\n')
            transinfo = geneinfo[gene2][3:]
            
            transcripts = [trans[0] for trans in transinfo]
            
            final_list = []
            for trans in transinfo:
                new_trans = [trans[0], trans[1]]
                
                for num in trans[2]:
                    new_trans.append(num)
                
                final_list.append(new_trans)
            
            pretty_df = {}
            pretty_df[gene2] = final_list
            
            cols = ['Transcript', 'TranscriptType']
            cols.extend(geneinfo[gene2][2])
            
            pretty_df = pretty_df[gene2]
            df = pd.DataFrame(pretty_df, index = transcripts, columns = cols)
            df2 = df.drop(['Transcript'], axis = 1)
            return(df2)
        else:
            return(geneinfo)
    else:
        return(gene2)

def gene_data(all_data, gene):
    gene2 = gene_name(all_data, gene)
    genes = all_data[0]
        
    if gene2 in genes:
        geneinfo = gene_info(all_data, gene2)[gene2]
        
        genename = geneinfo[0]
        genetype = geneinfo[1]
        sampnames = geneinfo[2]
        nsamps = len(sampnames)
        
        genedata = {}
        genedata[gene2] = [genename, genetype, sampnames, nsamps]
        return(genedata)
    else:
        return(gene2)

def gene_statistics(all_data, gene, pretty = 'F'):
    gene2 = gene_name(all_data, gene)
    genes = all_data[0]
    
    if gene2 in genes:
        geneinfo = gene_info(all_data, gene2)[gene2]
        genedata = gene_data(all_data, gene2)[gene2]
        
        if genedata[3] == 0:
            return('The gene ' + gene2 + ' is not expressed in any of the analysed samples.')
        else:
            def takeMean(elem):
                return elem[2]
            
            isoforms = []
            list_stats = []
            for trans in geneinfo[3:]:
                iso = trans[0]
                isoforms.append(iso)
                isotype = trans[1]
                mn = np.mean(trans[2])
                
                if len(trans[2]) == 1:
                    sd = 0.0
                else:
                    sd = np.std(np.array(trans[2]))
                try: 
                    minim, Q1, med, Q3, maxim = np.percentile(trans[2], [0, 25, 50, 75, 100])
                except:
                    #print(trans[2])
                    #raise
                    minim, Q1, med, Q3, maxim = [-1]*5
                
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
                print('\n' + gene2 + ' - ' + genedata[0] + ' ---> ' + genedata[1] + '\n')
                pretty_df = gene_stats[gene2][4:]
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
            
            means = [trans[2] for trans in genestats[4:]]
            means_array = np.asarray(means)
            high_props = list(means_array[means_array > minexp])
            ntrans = len(high_props)
            
            if ntrans == 0:
                return('There is not any isoform of the gene ' + gene2 + ' with a mean expression higher than ' + str(minexp))
            else:
                transcripts = [trans[0] for trans in genestats[4:][:ntrans]]
                trans_type = [trans[1] for trans in genestats[4:][:ntrans]]
                cumu_props = [trans[3] for trans in genestats[4:][:ntrans]]
                
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
                        print('The gene ' + gene2 + ' is just explained in a ' + str(round(last_cum * 100, 2)) + '%.')
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
                print('\n' + gene2 + ' - ' + genedata[0] + ' ---> ' + genedata[1] + '\n')
                pretty_df = gene_filt_props[gene2][4:]
                df = pd.DataFrame(pretty_df, index = results[0], columns = ['Transcript', 'TranscriptType', 'Mean', 'CumulativeMean', 'NewMean', 'NewCumulativeMean'])
                df2 = df.drop(['Transcript'], axis = 1)
                return(df2)
            else:
                return(gene_filt_props)
    else:
        return(gene2)

def gene_classification(all_data, gene, pretty = 'F', minexp = 0.1, minsamps = 10):
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
            if genedata[3] < minsamps:
                transcripts = transcript_types = means = cumu_means = new_means = new_cumu_means = '-'
                cog = 'FewSamples'
            else:
                gene_filt = gene_filt_props[gene2][4:]
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
        Types[gene2] = genedata + [cog, transcripts, transcript_types, means, cumu_means, new_means, new_cumu_means]
        
        if pretty == 'T':
            print('\n' + gene2 + ' - ' + genedata[0] + ', ' + genedata[1] + ' ---> ' + cog)
            
            if cog == 'Monoform' or cog == 'Biform' or cog == 'Triform' or cog == 'Multiform':
                print('\n')
                df = pd.DataFrame(gene_filt, index = transcripts, columns = ['Transcript', 'TranscriptType', 'Mean', 'CumulativeMean', 'NewMean', 'NewCumulativeMean'])
                df2 = df.drop(['Transcript'], axis = 1)
                return(df2)
            else:
                return('')
        else:
            return(Types)
    else:
        return(gene2)

def tissue_summary(csv_file, out_tsdir = 'results', out_tsfile = '', minexp = 0.1, minsamps = 10):
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
    
    if not os.path.isdir(out_tsdir):
        os.mkdir(out_tsdir)
    
    if out_tsfile == '':
        minexp_txt = str(round(minexp, 2))
        minsamps_txt = str(int(minsamps))
        
        if len(minexp_txt) < 4:
            minexp_list = list(minexp_txt)
            minexp_list.append('0')
            minexp_txt = ''.join(minexp_list)
        
        bstissuefile = os.path.abspath(out_tsdir) + '/' + os.path.basename(csv_file[ni:nf]) + '_' + minexp_txt + '_' + minsamps_txt + '.csv'
    else:
        bstissuefile = os.path.abspath(out_tsdir) + '/' + os.path.basename(out_tsfile)
    
    header = ';'.join(['GeneId', 'GeneName', 'GeneType', 'ExpressedSamples', 'NumberOfExpressedSamples', 'GeneClassification', 'Transcripts', 'TranscriptTypes', 'Mean', 'CumulativeMean', 'NewMean', 'NewCumulativeMean'])
#    bstissuefile=bstissuefile.split('/')[-1] 
    with gzip.open(csv_file, 'rt') if csv_file.endswith('gz') else open(csv_file) as rd, open(bstissuefile,'w') as wr:
        wr.write(header + '\n')
        
        nsamps = 0
        samp_names = ''
        
        for line in rd:
            ff = line.strip().split('\t')
            if line.startswith('transc'):
                columns = ff
                samp_names = columns[5:]
                nsamps = len(samp_names)
            else:
                genename = ff[1]
                genetype = ff[2]
                geneid = ff[4]
                samps = list()
                samps = [samp_names[i] for i in range(0, nsamps)  if float(ff[i + 5]) >= 0]
                nums  = [float(i) for i in ff[5:]  if float(i) >= 0] 
                
                tmp = Data.get(geneid, None)
                if tmp is None:
                    # process me
                    if len(Data.keys()) > 0:
                        res = gene_classification([list(Data.keys()), Data, tissue, samp_names] , 0, 'F', minexp, minsamps)
                        
                        for k, v in (res.items()):
                            outline = [k]
                            tmp_list = [str(x) for x in v]
                            outline.extend(tmp_list)
                            wr.write(';'.join(outline) + '\n')
                        
                        if len(Data.keys()) > 1:
                            print(Data.keys())
                            raise
                        Data = dict()
                    Data[geneid] = [genename, genetype, samps, [ff[0], ff[3], nums]]
                else:
                    tmp.append([ff[0], ff[3], nums])
                    Data[geneid] = tmp
        
        res = gene_classification([list(Data.keys()),Data, tissue, samp_names] , 0, 'F', minexp, minsamps)
        for k, v in (res.items()):
            outline = [k]
            tmp_list = [str(x) for x in v]
            outline.extend(tmp_list)
            wr.write(';'.join(outline) + '\n')
    
    return None

def tissue_statistics(in_tsfile, savefile = 'T', out_statsfile = '', genetype = '', drop_tsfile = 'T'):
    df = pd.read_csv(in_tsfile, index_col = 0, sep = ';')
    df2 = pd.crosstab(df.GeneClassification, df.GeneType)
    
    nrows = df2.shape[0]
    ncols = df2.shape[1]
    
    if nrows != 8:
        if 'Biform' not in df2.index:
            biform = pd.DataFrame([[0 for i in range(0, ncols)]], dtype = "int", columns = df2.columns)
            df2 = df2.append(biform)
            df2.rename(index = {0:'Biform'}, inplace = True)
            
        if 'FewSamples' not in df2.index:
            fewsamp = pd.DataFrame([[0 for i in range(0, ncols)]], dtype = "int", columns = df2.columns)
            df2 = df2.append(fewsamp)
            df2.rename(index = {0:'FewSamples'}, inplace = True)
        
        if 'LowExpressedTranscripts' not in df2.index:
            lowexp = pd.DataFrame([[0 for i in range(0, ncols)]], dtype = "int", columns = df2.columns)
            df2 = df2.append(lowexp)
            df2.rename(index = {0:'LowExpressedTranscripts'}, inplace = True)
        
        if 'Monoform' not in df2.index:
            monoform = pd.DataFrame([[0 for i in range(0, ncols)]], dtype = "int", columns = df2.columns)
            df2 = df2.append(monoform)
            df2.rename(index = {0:'Monoform'}, inplace = True)
        
        if 'Multiform' not in df2.index:
            multiform = pd.DataFrame([[0 for i in range(0, ncols)]], dtype = "int", columns = df2.columns)
            df2 = df2.append(multiform)
            df2.rename(index = {0:'Multiform'}, inplace = True)
        
        if 'NotExpressed' not in df2.index:
            notexp = pd.DataFrame([[0 for i in range(0, ncols)]], dtype = "int", columns = df2.columns)
            df2 = df2.append(notexp)
            df2.rename(index = {0:'NotExpressed'}, inplace = True)
        
        if 'Triform' not in df2.index:
            triform = pd.DataFrame([[0 for i in range(0, ncols)]], dtype = "int", columns = df2.columns)
            df2 = df2.append(triform)
            df2.rename(index = {0:'Triform'}, inplace = True)
    
    df2.sort_index(inplace = True)
    
    if genetype != '':
        df2 = df2[genetype]
    
    df2['Total'] = df2.sum(axis = 1)
    
    if drop_tsfile == 'T':
        os.remove(in_tsfile) 
    
    if savefile == 'T':
        if out_statsfile == '':
            out_statsfile = in_tsfile[:-4] + '_statistics.csv'
        
        df2.to_csv(out_statsfile)
    else:
        return(df2)

def tissue_difthres_summaries(csv_file, seqnumsamps, out_thresdir = 'comparing_thresholds', seqexp = 0.05,ncpus=1):
    minexp_seq = np.arange(0, 1, seqexp)
    from joblib import Parallel, delayed
    import multiprocessing
    num_cores = 20
    def dummy_func(i,j):
        tissue_summary(csv_file, out_tsdir = out_thresdir, out_tsfile = '', minexp = i, minsamps = j)
        
    import itertools    
    #os.system("taskset -p 0xff %d" % os.getpid())
    result =  Parallel(n_jobs=num_cores,prefer="processes")(delayed(dummy_func)(i,j) for i,j in list(itertools.product(minexp_seq,seqnumsamps)))

    
    return None

def tissue_difthres_statistics(in_thresdir, genetype = '', drop_tsfile = 'T'):
    if in_thresdir[-1] != '/':
        in_thresdir = in_thresdir + '/'
    
    for file in os.listdir(in_thresdir):
        infile = in_thresdir + file
        tissue_statistics(infile, 'T', '', genetype, drop_tsfile)

def tissue_difthres_stats_1file(in_thresdir):
    if in_thresdir[-1] != '/':
        in_thresdir = in_thresdir + '/'
    
    all_files = glob.glob(in_thresdir + '*statistics.csv')
    return(sorted(all_files))
#    for file in all_files:
#        print(file)










# ----------------------------------------------- PLOTS ----------------------------------------------- #
# --- GENE PLOTS --- #
def gene_boxplot(all_data, gene, ordered = True):
    gene2 = gene_name(all_data, gene)
    genes = all_data[0]
    
    if gene2 in genes:
        geneinfo = gene_info(all_data, gene)[gene2]
        
        nsamps = len(geneinfo[2])
        
        if nsamps == 0:
            return('The gene ' + gene2 + ' is not expressed in any of the analysed samples.')
        else:
            if ordered == True:
                def takeMedian(elem):
                    return np.median(elem[2])
                
                transcripts = sorted(geneinfo[3:], key = takeMedian, reverse = True)
            else:
                transcripts = geneinfo[3:]
                
            transnames = [trans[0] for trans in transcripts]
            numbers = [trans[2] for trans in transcripts]
            
            fig, ax = plt.subplots()
            ax.boxplot(numbers)
            ax.set_xlabel('Transcripts')
            ax.set_ylabel('Expression')
            ax.set_title('Gene ' + gene2)
            ax.set_xticklabels(transnames, rotation = 'vertical', fontsize = 8)
            ax.yaxis.grid(True, linestyle = '-', which = 'major', color = 'lightgrey', alpha = 0.5)
            axes = plt.gca()
            axes.set_ylim([-0.05, 1])
            plt.subplots_adjust(bottom = 0.22, top = 0.95)
            plt.show()
    else:
        return(gene2)

def gene_matrix(all_data, gene):
    gene2 = gene_name(all_data, gene)
    genes = all_data[0]
    
    if gene2 in genes:
        geneinfo = gene_info(all_data, gene)[gene2]
        
        sampnames = all_data[3]
        nsamps = len(sampnames)
        indx = np.arange(nsamps)
        
        transcripts = [trans[0] for trans in geneinfo[3:]]
        ntrans = len(transcripts)
        indy = np.arange(ntrans) + 0.5
        
        fig, ax = plt.subplots()
        width = 1
        
        acum_vect = [0 for i in range(0, nsamps)]
        sum_vect = [1 for i in range(0, nsamps)]
        
        if len(geneinfo[2]) == 0:
            for iso in transcripts:
                ax.bar(indx, sum_vect, width, bottom = acum_vect, color = 'white', edgecolor = 'gray')
                acum_vect = tuple(sum(x) for x in zip(acum_vect, sum_vect))
        else:
            for trans in geneinfo[3:]:
                i = 0
                colors = []
                for samp in sampnames:
                    if samp in geneinfo[2]:
                        if trans[2][i] > 0:
                            colors.append('blue')
                        else:
                            colors.append('white')
                        i += 1
                    else:
                        colors.append('white')
                        
                ax.bar(indx, sum_vect, width, bottom = acum_vect, color = colors, edgecolor = 'gray')
                acum_vect = tuple(sum(x) for x in zip(acum_vect, sum_vect))
                
        
        ax.set_xlabel('Samples')
        ax.set_ylabel('Transcripts')
        ax.set_title('Gene ' + gene2)
        ax.set_xticks(indx)
        ax.set_xticklabels(sampnames, rotation = 'vertical', fontsize = 8)
        ax.set_yticks(indy)
        ax.set_yticklabels(transcripts, fontsize = 8)
        plt.subplots_adjust(bottom = 0.22, top = 0.95, left = 0.2)
        plt.show()
    else:
        return(gene2)

def gene_barplot(all_data, gene):
    gene2 = gene_name(all_data, gene)
    genes = all_data[0]
    
    if gene2 in genes:
        geneinfo = gene_info(all_data, gene)[gene2]
        sampnames = geneinfo[2]
        nsamps = len(sampnames)
        
        if nsamps == 0:
            return('The gene ' + gene2 + ' is not expressed in any of the analysed samples.')
        else:
            indx = np.arange(nsamps)
            
            fig, ax = plt.subplots()
            width = 0.75
            
            acum_vect = [0 for i in range(0, nsamps)]
        
            for trans in geneinfo[3:]:
                ax.bar(indx, trans[2], width, bottom = acum_vect, label = trans[0])
                acum_vect = tuple(sum(x) for x in zip(acum_vect, trans[2]))
            
            ax.set_title('Gene ' + gene2)
            ax.set_xlabel("Samples")
            ax.set_ylabel("Proportion")
            ax.set_xticks(indx)
            ax.set_xticklabels(sampnames, rotation = 'vertical', fontsize = 8)
            ax.set_ylim([-0.05, 1])
            plt.subplots_adjust(bottom = 0.22, top = 0.95)
            plt.legend()
            plt.show()
    else:
        return(gene2)

def gene_filtered_barplot(all_data, gene, minexp = 0.1):
    gene2 = gene_name(all_data, gene)
    genes = all_data[0]
    
    if gene2 in genes:
        genefilt = gene_filtered_proportions(all_data, gene, 'F', minexp)
        
        if 'not expressed' in genefilt:
            return(genefilt)
        else:
            genefilt = genefilt[gene2]
            geneinfo = gene_info(all_data, gene)[gene2]
            sampnames = geneinfo[2]
            nsamps = len(sampnames)
            
            selected_trans = [trans[0] for trans in genefilt[4:-1]]
            
            selected_numbs = []
            for transcript in geneinfo[3:]:
                if transcript[0] in selected_trans:
                    selected_numbs.append(transcript)
            
            indx = np.arange(nsamps)
            
            fig, ax = plt.subplots()
            width = 0.75
            
            acum_vect = [0 for i in range(0, nsamps)]
        
            for trans in selected_numbs:
                ax.bar(indx, trans[2], width, bottom = acum_vect, label = trans[0])
                acum_vect = tuple(sum(x) for x in zip(acum_vect, trans[2]))
            
            sum_vect = [1 for i in range(0, nsamps)]
            other = tuple(round(num, 8) for num in np.subtract(sum_vect, acum_vect))
            
            ax.bar(indx, other, width, bottom = acum_vect, label = 'other')
                
            ax.set_xlabel('Samples')
            ax.set_ylabel('Expression')
            ax.set_title('Gene ' + gene2)
            ax.set_xticks(indx)
            ax.set_xticklabels(sampnames, rotation = 'vertical', fontsize = 8)
            axes = plt.gca()
            axes.set_ylim([-0.05, 1])
            plt.subplots_adjust(bottom = 0.22, top = 0.95)
            plt.legend()
            plt.show()
    else:
        return(gene2)

# --- TISSUE PLOTS --- #

def tissue_all_barplot(in_statsfile):
    stats = pd.read_csv(in_statsfile, index_col = 0)
    
    types = list(stats.index)
    counts = list(stats.Total)
    total = sum(counts)
    
    props = [count / total for count in counts]
    
    labs = []
    for i in range(0, len(types)):
        labs.append(types[i] + ': ' + str(counts[i]))
    
    ntypes = len(types)
    
    fig, ax = plt.subplots()
    
    acum = 0
    for i in range(0, ntypes):
        ax.bar(1, props[i], 1, bottom = acum, label = labs[i])
        acum =  acum + props[i]
    
    ax.set_xlabel(in_statsfile)
    ax.set_ylabel('Proportion')
    ax.set_title('All genes distribution')
    ax.set_xticks([0])
    ax.legend()
    ax.set_ylim([-0.05, 1])
    plt.subplots_adjust(bottom = 0.1, top = 0.9)
    plt.show()

def tissue_exp_barplot(in_statsfile):
    stats = pd.read_csv(in_statsfile, index_col = 0)
    df = stats.loc[['Biform', 'Monoform', 'Multiform', 'Triform']]
    
    types = list(df.index)
    counts = list(df.Total)
    total = sum(counts)
    
    props = [count / total for count in counts]
    
    labs = []
    for i in range(0, len(types)):
        labs.append(types[i] + ': ' + str(counts[i]))
    
    ntypes = len(types)
    colors = ['C0', 'C3', 'C4', 'C6']
    
    fig, ax = plt.subplots()
    
    acum = props[0]
    ax.bar(1, props[0], 1, label = labs[0], color = colors[0])
    for i in range(1, ntypes):
        ax.bar(1, props[i], 1, bottom = acum, label = labs[i], color = colors[i])
        acum =  acum + props[i]
    
    ax.set_xlabel(in_statsfile)
    ax.set_ylabel('Proportion')
    ax.set_title('Expressed gene distribution')
    ax.set_xticks([0])
    ax.legend()
    ax.set_ylim([-0.05, 1])
    plt.subplots_adjust(bottom = 0.1, top = 0.9)
    plt.show()

#def tissue_difthres_barplot(in_thresdir, nsamps):
#    if in_thresdir[-1] != '/':
#        in_thresdir = in_thresdir + '/'
#    
#    all_files = glob.glob(in_thresdir + '*statistics.csv')
#    
#    for file in all_files:
#        underscores = [x.start() for x in re.finditer('_', file)]
#        minsamps = file[underscores[-2] + 1 : underscores[-1]]
#        
#        stats = pd.read_csv(file, index_col = 0)
#        
#        types = list(stats.index)
#        counts = list(stats.Total)
#        total = sum(counts)
#        
#        props = [count / total for count in counts]
#        
#        labs = []
#        for i in range(0, len(types)):
#            labs.append(types[i] + ': ' + str(counts[i]))
#        
#        ntypes = len(types)
#        
#        ax = plt.subplot2grid((nsamps, 1), (0, 0))
#        
#        acum = 0
#        for i in range(0, ntypes):
#            ax.bar(1, props[i], 1, bottom = acum, label = labs[i])
#            acum =  acum + props[i]
#        
#    ax.set_xlabel('Different thresholds')
#    ax.set_ylabel('Proportion')
#    ax.set_title('All genes distribution')
#    ax.set_xticks([0])
#    ax.legend()
#    ax.set_ylim([-0.05, 1])
#    plt.subplots_adjust(bottom = 0.1, top = 0.9)
#    plt.show()

# --- ALL TISSUES PLOTS --- #

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
