#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 19 12:06:44 2018

@author: aserrano
"""
import os
import gzip
import numpy
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
                nsamps = len(ff[2:])
                samp_names = tuple(ff[2:])

            else:
                # here I parse all lines to generate the Data dictionary
                # keys = 'gene_names'
                # Data[gene] = list of isos
                # each iso  = [ENST, samp]
                gene_id = ff[1]
                nums  = list()
                samps = list()
                for i in range(0, nsamps):
                        num = float(ff[i + 2])
                        if num >= 0:
                            samps.append(samp_names[i])
                            nums.append(num)
                # here I update or generate the data[gene_id] value
                tmp  = Data.get(gene_id, None)
                if tmp is None:
                    Data[gene_id] = [[ff[0], tuple(samps), tuple(nums)]]
                else:
                    tmp.append([ff[0], tuple(samps), tuple(nums)])
                    Data[gene_id] = tmp

    return([Data.keys(), Data, tissue])

def formatting_gene(all_data, gene):
    genes = list(all_data[0])
    
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

def gene_info(all_data, gene, ncols = 10):
    gene2 = formatting_gene(all_data, gene)
    genes = list(all_data[0])
        
    if gene2 in genes:
        print('Data of the gene: ' + gene2)
        
        gene_info = {}
        for iso_info in all_data[1][gene2]:
            gene_info[iso_info[0]] = iso_info[2]
        
        df = pd.DataFrame(gene_info, index = (iso_info[1]))
        df = df.T
        pd.set_option('max_columns', ncols)
        return(df)
    else:
        return(gene2)

def total_prop(all_data, gene):
    gene2 = formatting_gene(all_data, gene)
    genes = list(all_data[0])
        
    if gene2 in genes:
        dat = all_data[1]
        samp_names = dat[gene2][0][1]
        nsamps = len(samp_names)
        
        if nsamps == 0:
            return('The gene ' + gene2 + ' is not expressed in any of the analysed samples.')
        else:
            list_props = []
            Proportions = {}
            for iso, samps, props in dat[gene2]:
                iso_prop = sum(props) / float(nsamps)
                list_props.append([iso, iso_prop])
            
            Proportions[gene2] = list_props
            
            return(Proportions)
    else:
        return(gene2)

def summarized_props(all_data, gene, thres = 0.9):
    thres = float(thres)
    
    gene2 = formatting_gene(all_data, gene)
    genes = list(all_data[0])
        
    if gene2 in genes:
        gene_props = total_prop(all_data, gene)
        
        if type(gene_props) == dict:
            def takeSecond(elem):
                return elem[1]
        
            for gene, iso_props in gene_props.items():
                iso_props.sort(reverse = True, key = takeSecond)
                
                acum = [iso_props[0][1]]
                for i in range(1, len(iso_props)):
                    acum.append(iso_props[i][1] + acum[i - 1])
        
            if acum[0] == 0:
                final_props = gene_props[gene]
                print('This gene is not expressed in these samples.')
            else:
                ind_thres = numpy.searchsorted(acum, thres, side = 'right')
                if ind_thres == len(acum):
                    perc = round(acum[ind_thres - 1] * 100, 2)
                    print('The expressed isoforms just explain the ' + str(perc) + '% of the gene expression of these samples.')
                else:
                    rem = round(1 - acum[ind_thres], 3)
                
                    if rem == iso_props[-1][1] and rem != 0:
                        final_props = iso_props
                    elif rem == 0:
                        final_props = [iso_props[0:(ind_thres + 1)]]
                    else:
                        final_props = [iso_props[0:(ind_thres + 1)], ['other', rem]]
        
            return(final_props)
        else:
            return(gene_props)
    else:
        return(gene2)

def isoform_stats(all_data, gene, thres = 0.9, thres2 = 0.7):
    thres = float(thres)
    thres2 = float(thres2)
    
    gene2 = formatting_gene(all_data, gene)
    genes = list(all_data[0])
    
    if gene2 in genes:
        gene_props = summarized_props(all_data, gene, thres)
        
        if 'not expressed' in gene_props:
            return(gene_props)
        else:
            if (len(gene_props) == 2 and gene_props[1][0] == 'other') or (len(gene_props) == 1):
                keeped_isos = gene_props[0]
            else:
                keeped_isos = gene_props
            
            isos = [keeped_isos[0][0]]
            values = [keeped_isos[0][1]]
            acum = [keeped_isos[0][1]]
            for i in range(1, len(keeped_isos)):
                isos.append(keeped_isos[i][0])
                values.append(keeped_isos[i][1])
                acum.append(acum[i - 1] + keeped_isos[i][1])
            
            total = acum[-1]
            new_prop = []
            for num in acum:
                new_num = num / float(total)
                new_prop.append(new_num)

            return([isos, values, acum, new_prop])
    else:
        return(gene2)

def type_of_gene(all_data, gene, thres = 0.9, thres2 = 0.7, thres3 = 2):
    thres = float(thres)
    thres2 = float(thres2)
    thres3 = int(thres3)
    
    gene2 = formatting_gene(all_data, gene)
    genes = list(all_data[0])
    
    if gene2 in genes:
        isos_stats = isoform_stats(all_data, gene, thres, thres2)
        
        Types = {}
        if 'not expressed' in isos_stats:
            tog = 'NotExpressed'
            isos = '-'
            stats = '-'
        else:
            dat = all_data[1]
            samp_names = dat[gene2][0][1]
            nsamps = len(samp_names)
        
            if nsamps < thres3:
                tog = 'FewSamples'
                isos = '-'
                stats = '-'
            else:
                i = 0
                while isos_stats[3][i] < thres2:
                    i += 1
                
                if i == 0:
                    tog = 'Monoform'
                elif i == 1:
                    tog = 'Biform'
                elif i == 2:
                    tog = 'Triform'
                else:
                    tog = 'Multiform'
            
                isos = isos_stats[0][:(i + 1)]
                stats = isos_stats[1][:(i + 1)]
        
        Types[gene2] = [tog, isos, stats]
        return(Types)
    else:
        return(gene2)

def big_summary(all_data, thres = 0.9, thres2 = 0.7, thres3 = 10, bstissuefile = '_genes_'):
    thres = float(thres)
    thres2 = float(thres2)
    thres3 = int(thres3)
    
    tissue = all_data[2]
    
    all_types = {}
    for i in range(0, len(all_data[0])):
        res = type_of_gene(all_data, i, thres, thres2, thres3)
        all_types.update(res)
    
    df = pd.DataFrame(all_types, index = ['Type', 'Isoforms', 'IsoformsProportion'], dtype = 'category')
    df = df.T
    
    if bstissuefile != '_genes_':
        bstissuefile = bstissuefile
    else:
        bstissuefile = tissue + bstissuefile + str(thres * 100) + '_' + str(thres2 * 100) + '_' + str(thres3) + '.csv'
    
    df.to_csv(bstissuefile)

def statistics(all_data, bstissuefile, tissuestatsfile = '_statistics.csv', dffile = 'T'):
    df = pd.read_csv(bstissuefile)
    
    nnexp = df.loc[df['Type'] == 'NotExpressed']['Type'].count()
    nfew = df.loc[df['Type'] == 'FewSamples']['Type'].count()
    nmono = df.loc[df['Type'] == 'Monoform']['Type'].count()
    nbi = df.loc[df['Type'] == 'Biform']['Type'].count()
    ntri = df.loc[df['Type'] == 'Triform']['Type'].count()
    nmulti = df.loc[df['Type'] == 'Multiform']['Type'].count()
    numbs = [nnexp, nfew, nmono, nbi, ntri, nmulti]
    
    ntotal = sum(numbs)
    
    pnexp = nnexp / float(ntotal)
    pfew = nfew / float(ntotal)
    pmono = nmono / float(ntotal)
    pbi = nbi / float(ntotal)
    ptri = ntri / float(ntotal)
    pmulti = nmulti / float(ntotal)
    props = [pnexp, pfew, pmono, pbi, ptri, pmulti]
    
    types = ['NotExpressed', 'FewSamples', 'Monoform', 'Biform', 'Triform', 'Multiform']
    
    stats = {}
    for i in range(0, len(types)):
        stats[types[i]] = [numbs[i], props[i]]
    
    if dffile == 'T':
        df = pd.DataFrame(stats, index = ['Counts', 'Typeform'])
        df = df.T
        
        if tissuestatsfile != '_statistics.csv':
            tissuestatsfile = tissuestatsfile
        else:
            tissue = all_data[2]
            tissuestatsfile = tissue + tissuestatsfile
        
        df.to_csv(tissuestatsfile)
    elif dffile == 'F':
        return(stats)
    else:
        print(dffile + 'is not an acceptable option for --dffile.')

def all_tissues_analysis(original_dir, pref = '', new_dir = '', thres = 0.9, thres2 = 0.7, thres3 = 10):
    thres = float(thres)
    thres2 = float(thres2)
    thres3 = int(thres3)
    
    if pref != '':
        pref = pref + '_'
    
    if new_dir == '':
        bsumm_dir = original_dir + pref + 'big_summaries'
    else:
        bsumm_dir = new_dir
    
    os.mkdir(bsumm_dir)
    os.chdir(original_dir)
    
    nfiles = 0
    for file in os.listdir():
        if file.startswith(pref) and file.endswith('gz'):
            nfiles += 1
    
    print('There are ' + str(nfiles) + ' tissues to analyse.\n')
    
    i = 0
    for file in os.listdir():
        if pref != '':
            print("pref != ''")
            if pref == 'SMTS_':
                ni = 5
            elif pref == 'SMTSD_':
                ni = 6
            
            if file.endswith('gz'):
                nf = -7
            elif file.endswith('csv'):
                nf = -4
            
            if file.startswith(pref) and file.endswith('gz'):
                tissuepath = original_dir + file
                tissue = file[ni:nf]
                
                DATA = reading_data(tissuepath)
                
                os.chdir(bsumm_dir)
                big_summary(DATA, thres, thres2, thres3, filename = 'BS_' + tissue + '.csv')
                
                os.chdir(original_dir)
                
                i += 1
                print(str(i) + '.- Big summary of the tissue ' + tissue + ': done.')
        else:
            tissuepath = original_dir + file
            tissue = file[:-3]
            
            DATA = reading_data(tissuepath)
            
            os.chdir(bsumm_dir)
            big_summary(DATA, thres, thres2, thres3, filename = 'BS_' + tissue + '.csv')
            
            os.chdir(original_dir)
        
            i += 1
            print(str(i[0]) + '.- Big summary of the tissue ' + tissue + ': done.')

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
def general_view(all_data, gene):
    gene2 = formatting_gene(all_data, gene)
    genes = list(all_data[0])
    
    if gene2 in genes:
        dat = all_data[1]
        samp_names = dat[gene2][0][1]
        nsamps = len(samp_names)
        
        if nsamps == 0:
            return('The gene ' + gene2 + ' is not expressed in any of the analysed samples.')
        else:
            n_isos = len(dat[gene2])
            lim = (n_isos - 1) / 2.0
            ind = numpy.arange(nsamps)
            
            fig, ax = plt.subplots()
            width = 0.75
            
            for i in range(0, n_isos):
                info_gene = dat[gene2][i]
                ax.bar(ind + (i - lim) * width / float(n_isos), info_gene[2], width / float(n_isos), label = info_gene[0])
        
            ax.set_ylabel('Proportion')
            ax.set_title('Distribution of samples of the gene ' + gene2)
            ax.set_xticks(ind)
            ax.set_xticklabels(samp_names, rotation = 'vertical', fontsize = 8)
            ax.legend()
            plt.subplots_adjust(bottom = 0.2, top = 0.9)
            plt.show()
    else:
        return(gene2)

def pie_plot(all_data, gene, thres = 0.9):
    thres = float(thres)
    
    gene2 = formatting_gene(all_data, gene)
    genes = list(all_data[0])
        
    if gene2 in genes:
        gene_props = summarized_props(all_data, gene, thres)
        
        if 'not expressed' in gene_props:
            return(gene_props)
        else:
            final_isos = []
            final_props = []
            if len(gene_props) == 2:
                isos = gene_props[0]
                other = gene_props[1]
            
                for iso in isos:
                    final_isos.append(iso[0])
                    final_props.append(iso[1])
            
                final_isos.append(other[0])
                final_props.append(other[1])
            else:
                for iso in gene_props[0]:
                    final_isos.append(iso[0])
                    final_props.append(iso[1])
                
            fig1, ax1 = plt.subplots()
            
            ax1.pie(final_props, labels = final_isos, autopct = '%.2f%%', shadow = True, startangle = 0)
            ax1.set_title('Distribution of samples of the gene ' + gene2)
            ax1.axis('equal')
            plt.show()
    else:
        return(gene2)

def stacked_barplot(all_data, gene, thres = 0.9):
    thres = float(thres)
    
    gene2 = formatting_gene(all_data, gene)
    genes = list(all_data[0])
    
    if gene2 in genes:
        gene_props = summarized_props(all_data, gene, thres)
        
        if 'not expressed' in gene_props:
            return(gene_props)
        else:
            dat = all_data[1]
            samp_names = dat[gene2][0][1]
            nsamps = len(samp_names)
            
            if (len(gene_props) == 2 and gene_props[1][0] == 'other') or (len(gene_props) == 1):
                keeped_isos = gene_props[0]
            else:
                keeped_isos = gene_props
            
            sorted_iso = []
            for iso in keeped_isos:
                for iso_data in dat[gene2]:
                    if iso[0] == iso_data[0]:
                        sorted_iso.append(iso_data)
            
            n_isos = len(sorted_iso)
            ind = numpy.arange(nsamps)
            
            fig, ax = plt.subplots()
            width = 0.5
            
            acum_vect = sorted_iso[0][2]
            ax.bar(ind, acum_vect, width, label = sorted_iso[0][0])
            for i in range(1, n_isos):
                ax.bar(ind, sorted_iso[i][2], width, bottom = acum_vect, label = sorted_iso[i][0])
                acum_vect = tuple(sum(x) for x in zip(acum_vect, sorted_iso[i][2]))
            
            ax.set_ylabel('Proportion')
            ax.set_title('Distribution of samples of the gene ' + gene2)
            ax.set_xticks(ind)
            ax.set_xticklabels(samp_names, rotation = 'vertical', fontsize = 8)
            ax.legend()
            plt.subplots_adjust(bottom = 0.2, top = 0.9)
            plt.show()
    else:
        return(gene2)

def stats_barplot(all_data, csv_file):
    stats = statistics(all_data, csv_file, dffile = 'F')
    tissue = all_data[2]
    
    types = list(stats.keys())
    
    numbs = []
    props = []
    for value in stats.values():
        numbs.append(int(value[0]))
        props.append(float(value[1]))
    
    labs = []
    for i in range(0, len(types)):
        labs.append(types[i] + ': ' + str(numbs[i]))
    
    ntypes = len(stats)
    acum = props[0]
    colors = ['C0', 'C9', 'C1', 'C2', 'C3', 'C5']
    
    fig, ax = plt.subplots()
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

def notexp_pie(all_data, csv_file):
    stats = statistics(all_data, csv_file, dffile = 'F')
    tissue = all_data[2]
    
    numbs = []
    props = []
    for value in stats.values():
        numbs.append(int(value[0]))
        props.append(float(value[1]))
    
    ini_types = list(stats.keys())
    types = ini_types[0:2]
    types.append('Expressed')
    
    new_numbs = [numbs[0], numbs[1], sum(numbs[2:])]
    labs = []
    for i in range(0, len(types)):
        labs.append(types[i] + ': ' + str(new_numbs[i]))
    
    new_props = [props[0], props[1], sum(props[2:])]
    
    fig1, ax1 = plt.subplots()
    colors = ['C0', 'C9', 'C6']
    ax1.pie(new_props, labels = labs, autopct = '%1.1f%%', shadow = True, startangle = 0, colors = colors)
    ax1.set_title('Isoform distribution of ' + tissue + '\nExpressed and not expressed isoforms')
    ax1.axis('equal')
    plt.show()

def expr_pie(all_data, csv_file):
    stats = statistics(all_data, csv_file, dffile = 'F')
    tissue = all_data[2]
    
    numbs = []
    props = []
    for value in stats.values():
        numbs.append(int(value[0]))
        props.append(float(value[1]))
    
    types = list(stats.keys())[2:]
    new_numbs = numbs[2:]
    labs = []
    for i in range(0, len(types)):
        labs.append(types[i] + ': ' + str(new_numbs[i]))
    
    new_total = sum(new_numbs)
    new_props = []
    for numb in new_numbs:
        new_prop = numb / float(new_total)
        new_props.append(new_prop)
    
    colors = ['C1', 'C2', 'C3', 'C5']
    
    fig1, ax1 = plt.subplots()
    ax1.pie(new_props, labels = labs, autopct = '%1.1f%%', shadow = True, startangle = 0, colors = colors)
    ax1.set_title('Isoform distribution of ' + tissue + '\nExpressed isoforms')
    ax1.axis('equal')
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
    ind = numpy.arange(len(tissues))
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
    elif type(gene) == str or type(gene) == numpy.str_:
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
    elif type(gene) == str or type(gene) == numpy.str_:
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
