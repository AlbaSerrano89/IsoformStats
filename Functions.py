#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 19 12:06:44 2018
@author: aserrano
"""
import os
import gzip
from joblib import Parallel, delayed
import itertools
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pylab import gcf

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
            return('ERROR: The gene ' + gene + ' is not in your data.')
    else:
        return('Incorrect type of element for the gene: ' + str(type(gene)))
    
    return(gene2)

def gene_pos(all_data, gene):
    genes = all_data[0]
    
    if type(gene) == int:
        gene2 = str(gene)
    elif type(gene) == (str or np.str_):
        if gene in genes:
            gene2 = str(genes.index(gene))
        else:
            return('ERROR: The gene ' + gene + ' is not in your data.')
    else:
        return('Incorrect type of element for the gene: ' + str(type(gene)))
    
    return(gene2)

def gene_info(all_data, gene, dicttype = True):
    gene2 = gene_name(all_data, gene)
    genes = all_data[0]
    
    if gene2 in genes:
        geneinfo = {}
        geneinfo[gene2] = all_data[1][gene2]
        
        if dicttype:
            return(geneinfo)
        else:
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

def gene_statistics(all_data, gene, dicttype = True):
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
            
            list_stats = []
            for trans in geneinfo[3:]:
                iso = trans[0]
                isotype = trans[1]
                mn = np.mean(trans[2])
                
                if len(trans[2]) == 1:
                    sd = 0.0
                else:
                    sd = np.std(np.array(trans[2]))
                try: 
                    minim, Q1, med, Q3, maxim = np.percentile(trans[2], [0, 25, 50, 75, 100])
                except:
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
            
            if dicttype:
                return(gene_stats)
            else:
                print('\n' + gene2 + ' - ' + genedata[0] + ' ---> ' + genedata[1] + '\n')
                pretty_df = gene_stats[gene2][4:]
                isoforms = [isoinfo[0] for isoinfo in pretty_df]
                df = pd.DataFrame(pretty_df, index = isoforms, columns = ['Transcript', 'TranscriptType', 'Mean', 'CumulativeMean', 'StandardDeviation', 'Q1', 'Median', 'Q3', 'Minimum', 'Maximum'])
                df2 = df.drop(['Transcript'], axis = 1)
                return(df2)
    else:
        return(gene2)

def gene_filtered_proportions(all_data, gene, dicttype = True, minexp = 0.8):
    gene2 = gene_name(all_data, gene)
    genes = all_data[0]
        
    if gene2 in genes:
        gene_stats = gene_statistics(all_data, gene2)
        
        if 'not expressed' in gene_stats:
            return(gene_stats)
        else:
            genestats = gene_stats[gene2]
            genedata = gene_data(all_data, gene2)[gene2]
            
            cum_means = [trans[3] for trans in genestats[4:]]
            ind_thres = np.searchsorted(cum_means, minexp, side = 'right')
            ntrans = ind_thres + 1
            
            transcripts = [trans[0] for trans in genestats[4:][:ntrans]]
            trans_type = [trans[1] for trans in genestats[4:][:ntrans]]
            high_props = [trans[2] for trans in genestats[4:][:ntrans]]
            cumu_props = [trans[3] for trans in genestats[4:][:ntrans]]
            
            total_sum = np.float64(cumu_props[-1])
            
            new_means = list(high_props / total_sum)
            new_cumus = list(cumu_props / total_sum)
            
            if ntrans != len(genestats[3:]):
                transcripts.append('other')
                trans_type.append('-')
                high_props.append(1 - sum(high_props))
                
                last_cum = cum_means[-1]
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
            
            if dicttype:
                return(gene_filt_props)
            else:
                print('\n' + gene2 + ' - ' + genedata[0] + ' ---> ' + genedata[1] + '\n')
                pretty_df = gene_filt_props[gene2][4:]
                df = pd.DataFrame(pretty_df, index = results[0], columns = ['Transcript', 'TranscriptType', 'Mean', 'CumulativeMean', 'NewMean', 'NewCumulativeMean'])
                df2 = df.drop(['Transcript'], axis = 1)
                return(df2)
    else:
        return(gene2)

def gene_classification(all_data, gene, dicttype = True, minexp = 0.8, minsamps = 10):
    gene2 = gene_name(all_data, gene)
    genes = all_data[0]
    
    if gene2 in genes:
        gene_filt_props = gene_filtered_proportions(all_data, gene, True, minexp)
        genedata = gene_data(all_data, gene2)[gene2]
        
        if type(gene_filt_props) == str:
            transcripts = transcript_types = means = cumu_means = new_means = new_cumu_means = '-'
            cog = 'NotExpressed'
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
        
        if dicttype:
            return(Types)
        else:
            print('\n' + gene2 + ' - ' + genedata[0] + ', ' + genedata[1] + ' ---> ' + cog)
            
            if cog == 'Monoform' or cog == 'Biform' or cog == 'Triform' or cog == 'Multiform':
                print('\n')
                df = pd.DataFrame(gene_filt, index = transcripts, columns = ['Transcript', 'TranscriptType', 'Mean', 'CumulativeMean', 'NewMean', 'NewCumulativeMean'])
                df2 = df.drop(['Transcript'], axis = 1)
                return(df2)
            else:
                return('')
    else:
        return(gene2)

def tissue_summary(csv_file, out_tsdir = os.getcwd(), out_tsfile = '', minexp = 0.8, minsamps = 10):
    Data = dict()    
    lastbar = csv_file.rfind('/')
    file = csv_file[(lastbar + 1):]    
    
    ni = 0
    if file.startswith('SMTS_'):
        ni = 5
    elif file.startswith('SMTSD_'):
        ni = 6
    
    nf = len(file)
    if file.endswith('gz'):
        nf = -7
    elif file.endswith('csv'):
        nf = -4
        
    tissue = file[ni:nf]
    
    if not os.path.isdir(out_tsdir):
        os.mkdir(out_tsdir)
    
    if out_tsdir[-1] != '/':
        out_tsdir = out_tsdir + '/'
        
    if out_tsfile == '':
        minexp_txt = str(round(minexp, 2))
        minexp_txt = minexp_txt.ljust(5, '0')
        
        bstissuefile = '{}{}_{}_{}.csv'.format(out_tsdir, file[:nf], minexp_txt, minsamps)
    else:
        bstissuefile = out_tsdir + out_tsfile
    
    header = ';'.join(['GeneId', 'GeneName', 'GeneType', 'MinimumExpression', 'MinimumSamples', 'ExpressedSamples', 'NumberOfExpressedSamples', 'GeneClassification', 'Transcripts', 'TranscriptTypes', 'Mean', 'CumulativeMean', 'NewMean', 'NewCumulativeMean'])
    
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
                        res = gene_classification([list(Data.keys()), Data, tissue, samp_names] , 0, True, minexp, minsamps)
                        
                        for k, v in (res.items()):
                            outline = [k]
                            tmp_list = [str(x) for x in v]
                            outline.extend(tmp_list)
                            wr.write(';'.join(outline[:3] + [str(round(minexp, 2)), str(round(minsamps))] + outline[3:]) + '\n')
                        
                        if len(Data.keys()) > 1:
                            print(Data.keys())
                            raise
                        Data = dict()
                    Data[geneid] = [genename, genetype, samps, [ff[0], ff[3], nums]]
                else:
                    tmp.append([ff[0], ff[3], nums])
                    Data[geneid] = tmp
        
        res = gene_classification([list(Data.keys()), Data, tissue, samp_names] , 0, True, minexp, minsamps)
        for k, v in (res.items()):
            outline = [k]
            tmp_list = [str(x) for x in v]
            outline.extend(tmp_list)
            wr.write(';'.join(outline[:3] + [str(round(minexp, 2)), str(round(minsamps))] + outline[3:]) + '\n')
    
    return None

def tissue_statistics(in_tsfile, savefile = True, out_statsdir = os.getcwd(), out_statsfile = '', genetype = '', drop_tsfile = True):
    df = pd.read_csv(in_tsfile, index_col = 0, sep = ';')
    
    minexp = df['MinimumExpression'][0]
    minsamp = df['MinimumSamples'][0]
    
    df2 = pd.crosstab(df.GeneClassification, df.GeneType)
    
    nrows = df2.shape[0]
    ncols = df2.shape[1]
    
    if nrows != 8:
        if 'Monoform' not in df2.index:
            monoform = pd.DataFrame([[0 for i in range(0, ncols)]], dtype = "int", columns = df2.columns)
            df2 = df2.append(monoform)
            df2.rename(index = {0:'Monoform'}, inplace = True)
        
        if 'Biform' not in df2.index:
            biform = pd.DataFrame([[0 for i in range(0, ncols)]], dtype = "int", columns = df2.columns)
            df2 = df2.append(biform)
            df2.rename(index = {0:'Biform'}, inplace = True)
        
        if 'Triform' not in df2.index:
            triform = pd.DataFrame([[0 for i in range(0, ncols)]], dtype = "int", columns = df2.columns)
            df2 = df2.append(triform)
            df2.rename(index = {0:'Triform'}, inplace = True)
        
        if 'Multiform' not in df2.index:
            multiform = pd.DataFrame([[0 for i in range(0, ncols)]], dtype = "int", columns = df2.columns)
            df2 = df2.append(multiform)
            df2.rename(index = {0:'Multiform'}, inplace = True)
        
        if 'NotExpressed' not in df2.index:
            notexp = pd.DataFrame([[0 for i in range(0, ncols)]], dtype = "int", columns = df2.columns)
            df2 = df2.append(notexp)
            df2.rename(index = {0:'NotExpressed'}, inplace = True)
        
        if 'FewSamples' not in df2.index:
            fewsamp = pd.DataFrame([[0 for i in range(0, ncols)]], dtype = "int", columns = df2.columns)
            df2 = df2.append(fewsamp)
            df2.rename(index = {0:'FewSamples'}, inplace = True)
    
    df2.sort_index(inplace = True)
    
    if genetype != '':
        df2 = df2[genetype]
    
    df2['Total'] = df2.sum(axis = 1)
    df2.insert(0, 'MinimumExpression', [minexp for i in range(0, 6)])
    df2.insert(1, 'MinimumSamples', [minsamp for i in range(0, 6)])
    
    df2 = df2.T
    df2 = df2[['Monoform', 'Biform', 'Triform', 'Multiform', 'NotExpressed', 'FewSamples']]
    df2 = df2.T
    
    if savefile:
        if drop_tsfile:
            os.remove(in_tsfile)
        
        if not os.path.isdir(out_statsdir):
            os.mkdir(out_statsdir)
    
        if out_statsdir[-1] != '/':
            out_statsdir = out_statsdir + '/'
        
        if out_statsfile == '':
            out_statsfile = os.path.basename(in_tsfile)[:-4] + '_statistics.csv'
        
        final_file = out_statsdir + out_statsfile
        
        df2.to_csv(final_file)
    else:
        return(df2)

def tissue_difthres_summaries(csv_file, seqnumsamps, out_thresdir = '_DiffThres_Summaries', seqexp = 0.05, ncpus = 1, num_cores = 1):
    minexp_seq = np.arange(0, 1, seqexp)
    
    maxlen = max([len(str(x)) for x in seqnumsamps])
    seqnumsamps_new = [int(str(x).rjust(maxlen, '0')) for x in seqnumsamps]
    
    if out_thresdir == '_DiffThres_Summaries':
        out_thresdir = os.path.basename(csv_file)[:-7] + out_thresdir
    
    if not os.path.isdir(out_thresdir):
        os.mkdir(out_thresdir)
    
    def dummy_func(i, j):
        tissue_summary(csv_file, out_tsdir = out_thresdir, out_tsfile = '', minexp = i, minsamps = j)
    
    Parallel(n_jobs = num_cores, prefer = "processes")(delayed(dummy_func)(i, j) for i, j in list(itertools.product(minexp_seq, seqnumsamps_new)))
    
    return None

def tissue_difthres_statistics(in_thresdir, out_statsdir = os.getcwd(), out_statsfile = '', genetype = '', drop_tsfile = True):
    if in_thresdir[-1] != '/':
        in_thresdir = in_thresdir + '/'
    
    if not os.path.isdir(out_statsdir):
        os.mkdir(out_statsdir)

    if out_statsdir[-1] != '/':
        out_statsdir = out_statsdir + '/'
    
    df = pd.DataFrame()
    for file in sorted(os.listdir(in_thresdir)):
        infile = in_thresdir + file
        df = df.append(tissue_statistics(infile, False, '', '', genetype, drop_tsfile))
        
        if drop_tsfile:
            os.remove(infile)
        
    tissue_list = os.path.basename(infile).split('_')[:-2]
    tissue = '_'.join(tissue_list)
    
    if out_statsfile == '':
        out_statsfile = tissue + '_statistics.csv'
    
    df.to_csv(out_statsdir + out_statsfile)

def all_tissues_difthres_summaries(initial_dir, seqnumsamps, out_threstsdir = 'AllTissues_Summaries_DifThres/', seqexp = 0.05, ncpus = 1, num_cores = 1):
    if initial_dir[-1] != '/':
        initial_dir = initial_dir + '/'
    
    for file in sorted(os.listdir(initial_dir)):
        print(file)
        tissue_difthres_summaries(initial_dir + file, seqnumsamps, out_threstsdir, seqexp, ncpus, num_cores)

def all_tissues_difthres_statistics(in_thresdir = 'AllTissues_Summaries_DifThres/', stats_dir = 'AllTissues_Statistics_DifThres/', genetype = '', drop_tsfile = True):
    if in_thresdir[-1] != '/':
        in_thresdir = in_thresdir + '/'
    
    if not os.path.isdir(stats_dir):
        os.mkdir(stats_dir)
    
    beginnings = []
    for file in sorted(os.listdir(in_thresdir)):
        tissue = '_'.join(file.split('_')[:-2])
        beginnings.append(tissue)
    
    beginnings = sorted(list(set(beginnings)))
    for tissue in beginnings:
        df = pd.DataFrame()
        for file in sorted(os.listdir(in_thresdir)):
            if file.startswith(tissue + '_0'):
                infile = in_thresdir + file
                df = df.append(tissue_statistics(infile, False, '', '', genetype, drop_tsfile))
        
        df.to_csv(stats_dir + tissue + '_statistics.csv')













# ----------------------------------------------- PLOTS ----------------------------------------------- #

# --- GENE PLOTS --- #
def gene_matrix(all_data, gene, plotfile):
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
        
        fig2, ax2 = plt.subplots()
        width = 1
        
        acum_vect = [0 for i in range(0, nsamps)]
        sum_vect = [1 for i in range(0, nsamps)]
        
        if len(geneinfo[2]) == 0:
            for iso in transcripts:
                ax2.bar(indx, sum_vect, width, bottom = acum_vect, color = 'white', edgecolor = 'gray')
                acum_vect = tuple(sum(x) for x in zip(acum_vect, sum_vect))
        else:
            for trans in geneinfo[3:]:
                i = 0
                colors = []
                for samp in sampnames:
                    if samp in geneinfo[2]:
                        if trans[2][i] > 0:
                            colors.append('black')
                        else:
                            colors.append('white')
                        i += 1
                    else:
                        colors.append('white')
                        
                ax2.bar(indx, sum_vect, width, bottom = acum_vect, color = colors, edgecolor = 'gray')
                acum_vect = tuple(sum(x) for x in zip(acum_vect, sum_vect))
                
        
        ax2.set_xlabel('Samples')
        ax2.set_ylabel('Transcripts')
        ax2.set_title('Gene ' + gene2)
        ax2.set_xticks(indx)
        ax2.set_xticklabels(sampnames, rotation = 'vertical', fontsize = 8)
        ax2.set_yticks(indy)
        ax2.set_yticklabels(transcripts, fontsize = 8)
        plt.subplots_adjust(bottom = 0.41, top = 0.95, left = 0.25)
        plt.savefig(plotfile)
    else:
        return(gene2)

def gene_barplot(all_data, gene, plotfile):
    gene2 = gene_name(all_data, gene)
    genes = all_data[0]
    
    if gene2 in genes:
        geneinfo = gene_info(all_data, gene)[gene2]
        sampnames = geneinfo[2]
        nsamps = len(sampnames)
        
        gene_stats = gene_statistics(all_data, gene2)[gene2]
        
        final_gene_info = []
        for gnsts in gene_stats[4:]:
            for i in range(len(geneinfo[3:])):
                if gnsts[0] == geneinfo[3:][i][0]:
                    final_gene_info.append(geneinfo[3 + i])
        
        if nsamps == 0:
            return('The gene ' + gene2 + ' is not expressed in any of the analysed samples.')
        else:
            indx = np.arange(nsamps)
            
            fig3, ax3 = plt.subplots()
            width = 0.75
            
            acum_vect = [0 for j in range(0, nsamps)]
            
            for trans in final_gene_info:
                ax3.bar(indx, trans[2], width, bottom = acum_vect, label = trans[0])
                acum_vect = tuple(sum(x) for x in zip(acum_vect, trans[2]))
            
            ax3.set_title('Gene ' + gene2)
            ax3.set_xlabel("Samples")
            ax3.set_ylabel("Proportion")
            ax3.set_xticks(indx)
            ax3.set_xticklabels(sampnames, rotation = 'vertical', fontsize = 8)
            ax3.set_ylim([-0.05, 1])
            plt.subplots_adjust(bottom = 0.41, top = 0.95)
            plt.legend()
            plt.savefig(plotfile)
    else:
        return(gene2)

def gene_boxplot(all_data, gene, plotfile):
    gene2 = gene_name(all_data, gene)
    genes = all_data[0]
    
    if gene2 in genes:
        geneinfo = gene_info(all_data, gene)[gene2]
        
        nsamps = len(geneinfo[2])
        
        if nsamps == 0:
            return('The gene ' + gene2 + ' is not expressed in any of the analysed samples.')
        else:
            def takeMedian(elem):
                return np.median(elem[2])
            
            transcripts = sorted(geneinfo[3:], key = takeMedian, reverse = True)
            
            transnames = [trans[0] for trans in transcripts]
            numbers = [trans[2] for trans in transcripts]
            
            fig1, ax1 = plt.subplots()
            ax1.boxplot(numbers)
            ax1.set_xlabel('Transcripts')
            ax1.set_ylabel('Expression')
            ax1.set_title('Gene ' + gene2)
            ax1.set_xticklabels(transnames, rotation = 'vertical', fontsize = 8)
            ax1.yaxis.grid(True, linestyle = '-', which = 'major', color = 'lightgrey', alpha = 0.5)
            ax1.set_ylim([-0.05, 1])
            plt.subplots_adjust(bottom = 0.32, top = 0.95)
            plt.savefig(plotfile)
    else:
        return(gene2)

def gene_filtered_barplot(all_data, gene, plotfile, minexp = 0.8):
    gene2 = gene_name(all_data, gene)
    genes = all_data[0]
    
    if gene2 in genes:
        genefilt = gene_filtered_proportions(all_data, gene, True, minexp)
        
        if 'not expressed' in genefilt:
            return(genefilt)
        else:
            genefilt = genefilt[gene2]
            geneinfo = gene_info(all_data, gene)[gene2]
            sampnames = geneinfo[2]
            nsamps = len(sampnames)
            
            selected_trans = [trans[0] for trans in genefilt[4:-1]]
            
            final_gene_info = []
            for trans in selected_trans:
                for i in range(len(geneinfo[3:])):
                    if trans == geneinfo[3:][i][0]:
                        final_gene_info.append(geneinfo[3 + i])
            
            indx = np.arange(nsamps)
            
            fig4, ax4 = plt.subplots()
            width = 0.75
            
            acum_vect = [0 for i in range(0, nsamps)]
            
            for trans in final_gene_info:
                ax4.bar(indx, trans[2], width, bottom = acum_vect, label = trans[0])
                acum_vect = tuple(sum(x) for x in zip(acum_vect, trans[2]))
            
            sum_vect = [1 for i in range(0, nsamps)]
            other = tuple(round(num, 8) for num in np.subtract(sum_vect, acum_vect))
            
            ax4.bar(indx, other, width, bottom = acum_vect, color = 'gray', label = 'other')
                
            ax4.set_xlabel('Samples')
            ax4.set_ylabel('Expression')
            ax4.set_title('Gene {}\nMinimum expression: {}'.format(gene2, minexp))
            ax4.set_xticks(indx)
            ax4.set_xticklabels(sampnames, rotation = 'vertical', fontsize = 8)
            ax4.set_ylim([-0.05, 1])
            plt.subplots_adjust(bottom = 0.41, top = 0.9)
            plt.legend()
            plt.savefig(plotfile)
    else:
        return(gene2)

# --- TISSUE PLOTS --- #

def tissue_barplot(in_statsfile, plotfile = None, expressed = False, minexp = 0.8, minsamps = 10):
    stats = pd.read_csv(in_statsfile, index_col = 0)
    
    filename = os.path.basename(in_statsfile)[:-15]
    
    if stats.shape[0] > 6:
        stats = stats.loc[stats.MinimumExpression == minexp]
        stats = stats.loc[stats.MinimumSamples == minsamps]
    else:
        minexp = stats.MinimumExpression[0]
        minsamps = int(stats.MinimumSamples[0])
        filename = '_'.join(filename.split('_')[:-2])
    
    if plotfile == None:
        if expressed:
            plotfile = '{}_{}_{}_ExpressedGenesBarplot.png'.format(filename, minexp, minsamps)
        else:
            plotfile = '{}_{}_{}_AllGenesBarplot.png'.format(filename, minexp, minsamps)
    
    stats = stats.drop(['MinimumExpression', 'MinimumSamples'], axis = 1)
    
    if expressed:
        stats = stats.loc[['Monoform', 'Biform', 'Triform', 'Multiform']]
    
    types = list(stats.index)
    counts = list(stats.Total)
    total = sum(counts)
    
    props = [count / total for count in counts]
    
    labs = []
    for i in range(0, len(types)):
        labs.append('{}: {}'.format(types[i], int(counts[i])))
    
    ntypes = len(types)
    
    fig5, ax5 = plt.subplots()
    
    acum = 0
    for i in range(0, ntypes):
        ax5.bar(1, props[i], 1, bottom = acum, label = labs[i])
        acum =  acum + props[i]
    
    ax5.set_xlabel('MinimumExpression: {}\nMinimumSamples: {}'.format(minexp, minsamps))
    ax5.set_ylabel('Proportion')
    
    if expressed:
        ax5.set_title('Expressed gene distribution')
    else:
        ax5.set_title('All genes distribution')
    
    ax5.set_xticks([0])
    ax5.legend()
    ax5.set_ylim([-0.05, 1])
    plt.subplots_adjust(bottom = 0.1, top = 0.9)
    plt.savefig(plotfile)

def tissue_difthres_barplot(in_thresfile, plotfile = None, expressed = False, genetype = '', samplots = ''):
    stats = pd.read_csv(in_thresfile, index_col = 0)
    
    if genetype != '':
        stats1 = stats[['MinimumSamples', 'MinimumExpression'] + genetype]
        stats2 = stats[genetype]
        stats3 = stats2.sum(axis = 1)
        stats = pd.concat([stats1, stats2], axis = 1, sort = False)
        stats['Total'] = stats3
    
    if samplots == '':
        samps = stats.MinimumSamples.unique()
    else:
        samps = [int(x) for x in samplots]
    
    if expressed:
        stats = stats.loc[['Monoform', 'Biform', 'Triform', 'Multiform']]
        colors = ['C0', 'C1', 'C2', 'C3']
        ind = np.arange(4)
    else:
        colors = ['C0', 'C1', 'C2', 'C3', 'C4', 'C5']
        ind = np.arange(6)
    
    width = 0.85
    
    if plotfile == None:
        filename = os.path.basename(in_thresfile)[:-14]
        
        if genetype != '':
            gtps = '-'.join(genetype) + '_'
        else:
            gtps = ''
        
        if samplots != '':
            smplts = '-'.join(samplots) + '_'
        else:
            smplts = ''
        
        filename = '{}{}{}'.format(filename, gtps, smplts)
        if expressed:
            plotfile = '{}ExpressedDiffThresHistograms.png'.format(filename)
        else:
            plotfile = '{}AllDiffThresHistograms.png'.format(filename)
    
    fig6 = plt.figure()
    
    l = 0
    for i in samps:
        stats3 = stats.loc[stats.MinimumSamples == i]
        expres = stats3.MinimumExpression.unique()
        
        lim = (len(expres) - 1) / 2
        
        ax6 = fig6.add_subplot(len(samps), 1, l + 1)
        
        m = 0
        for j in expres:
            stats4 = stats3.loc[stats3.MinimumExpression == j]
            
            types = list(stats4.index)
            counts = stats4.Total
            total = sum(counts)
            
            props = [count / total for count in counts]
            
            for k in ind:
                ax6.bar(ind[k] + (m - lim) * width/len(expres), props[k], width/len(expres), color = colors[k])
            
            m += 1
        
        ax6.set_ylabel('Min. samples: ' + str(i), rotation = 'horizontal', fontsize = 10, labelpad = 50)
        
        if i == samps[-1]:
            ax6.set_xticks(ind)
            ax6.set_xticklabels(types, rotation = 'vertical', fontsize = 8)
        else:
            ax6.set_xticks(ind)
            ticks = ['' for x in range(7)]
            ax6.set_xticklabels(ticks, rotation = 'vertical', fontsize = 8)
        
        l += 1
    
    fig = gcf()
    
    if genetype == '':
        fig.suptitle('Comparing the different thresholds', fontsize = 14)
        top_marg = 0.93
    else:
        gtps = '\n'.join(genetype)
        fig.suptitle('Comparing the different thresholds for:\n' + gtps, fontsize = 10)
        ngtps = len(genetype)
        top_marg = 0.93 - 0.03 * ngtps
    
    plt.subplots_adjust(bottom = 0.2, top = top_marg, right = 0.98, left = 0.27)
    plt.savefig(plotfile)

# --- ALL TISSUES PLOTS --- #

def all_tissues_barplot(stats_directory, minexp, minsamps, expressed = False, plotfile = 'AllTissuesBarplot', genetype = ''):
    if stats_directory[-1] != '/':
        stats_directory = stats_directory + '/'
    
    df = pd.DataFrame()
    for file in os.listdir(stats_directory):
        df1 = pd.read_csv(stats_directory + file, index_col = 0)
        
        if genetype != '':
            stats1 = df1[['MinimumSamples', 'MinimumExpression']]
            stats2 = df1[genetype]
            stats3 = stats2.sum(axis = 1)
            df1 = pd.concat([stats1, stats2], axis = 1, sort = False)
            df1['Total'] = stats3
            
            gtps = '_' + '-'.join(genetype)
        else:
            gtps = ''
        
        df1 = df1.loc[df1.MinimumExpression == minexp]
        df1 = df1.loc[df1.MinimumSamples == minsamps]
        df1 = df1.Total
        df = df.append(df1.T)
        df.rename(index = {'Total': file[:-15]}, inplace = True)
        
    if expressed:
        df = df[['Monoform', 'Biform', 'Triform', 'Multiform']]
        colors = ['C0', 'C1', 'C2', 'C3']
        filename = 'Expr'
    else:
        df = df[['Monoform', 'Biform', 'Triform', 'Multiform', 'NotExpressed', 'FewSamples']]
        colors = ['C0', 'C1', 'C2', 'C3', 'C4', 'C5']
        filename = ''
    
    df2 = df.T
    
    types = list(df)
    tissues = list(df.index)
    
    props_dict = {}
    counts_dict = {}
    for tissue in tissues:
        counts_dict[tissue] = [int(x) for x in list(df2[tissue])]
        total_col = sum(df2[tissue])
        
        props_tiss = []
        for numbs in list(df2[tissue]):
            if total_col != 0:
                props_tiss.append(numbs / total_col)
            else:
                props_tiss.append(0)
        
        props_dict[tissue] = props_tiss
    
    props_df1 = pd.DataFrame(props_dict, index = types)
    props_df2 = props_df1.T
    
    if plotfile == 'AllTissuesBarplot':
        plotfile = '{}{}_{}_{}{}.png'.format(plotfile, filename, minexp, minsamps, gtps)
    
    fig, ax8 = plt.subplots()
    ind = np.arange(len(tissues))
    width = 0.5
    
    props_vect = []
    for typegene in types:
        props_vect.append(list(props_df2[typegene]))
    
    acum_vect = tuple(0 for x in props_vect[0])
    
    for i in range(0, len(types)):
        ax8.barh(ind, props_vect[i], width, left = acum_vect, label = types[i], color = colors[i])
        acum_vect = tuple(sum(x) for x in zip(acum_vect, props_vect[i]))
        
        for j in range(len(acum_vect)):
            if props_vect[i][j] != 0:
                pos_text = acum_vect[j] - 0.005
                text = str(counts_dict[tissues[j]][i])
            else:
                pos_text = 0
                text = ''
            
            plt.text(pos_text, j, text, color = 'black', fontsize = 'xx-small', ha = 'right', va = 'center')
    
    ax8.set_xlabel('Minimum expression: ' + str(minexp) + '\nMinimum samples: ' + str(minsamps))
    ax8.set_ylabel('Tissues')
    
    if genetype == '':
        fig.suptitle('Counts of the analysis', fontsize = 14)
        top_marg = 0.93
    else:
        gtps = '\n'.join(genetype)
        fig.suptitle('Counts of the analysis for:\n' + gtps, fontsize = 10)
        ngtps = len(genetype)
        top_marg = 0.93 - 0.03 * ngtps
    
    ax8.set_yticks(ind)
    ax8.set_yticklabels(tissues, fontsize = 8)
    ax8.legend()
    ax8.set_xlim([-0.05, 1])
    plt.subplots_adjust(bottom = 0.18, top = top_marg, left = 0.25, right = 0.98)
    plt.savefig(plotfile)