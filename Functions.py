#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 19 12:06:44 2018

@author: aserrano
"""
import numpy
import matplotlib.pyplot as plt
import pandas as pd
import gzip

def reading_data(csv_file):
    Data = dict()
    with gzip.open(csv_file) if csv_file.endswith('gz') else open(csv_file) as rd:
        #data = rd.read() # pfile.read()
        #pfile.close()
        
        nsamps = 0
        samp_names = ''
        
        for line in rd:
            ff = line.strip().split('\t')
            if line.startswith('transc'):
                        nsamps = len(ff[2:])
                        samp_names = tuple(ff[2:])

            else:
                # here I parse all lines
                # to generate the Data dictionary
                # keys are 'gene_names'
                # Data[gene] = list of isos
                # each iso  = [ ENST, samp]
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

    return([Data.keys(), Data, nsamps, samp_names])

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

def total_prop(all_data, gene):
    gene2 = formatting_gene(all_data, gene)
    genes = list(all_data[0])
        
    if gene2 in genes:
        dat = all_data[1]
        samp_names = dat[gene2][0][1]
        nsamps = len(samp_names)
        
        if nsamps == 0:
            return(['The gene ' + gene2 + ' is not expressed in any of the analysed samples.'])
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
        elif type(gene_props) == list:
            return(gene_props[0])
    else:
        return(gene2)

def pie_plot(all_data, gene, thres = 0.9):
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
            
            ax1.pie(final_props, labels = final_isos, autopct = '%1.1f%%', shadow = True, startangle = 0)
            ax1.set_title('Distribution of samples of the gene ' + gene2)
            ax1.axis('equal')
            plt.show()
    else:
        return(gene2)

def stacked_barplot(all_data, gene, thres = 0.9):
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

def isoform_stats(all_data, gene, thres = 0.9, thres2 = 0.7):
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

def type_of_gene(all_data, gene, thres = 0.9, thres2 = 0.7):
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

def big_summary(all_data, thres = 0.9, thres2 = 0.7, filename = 'classif_genes_'):
    all_types = {}
    for i in range(0, len(all_data[0])):
        res = type_of_gene(all_data, i, thres, thres2)
        all_types.update(res)
    
    df = pd.DataFrame(all_types, index = ['Type', 'Isoforms', 'IsoformsProportion'], dtype = "category")
    df = df.T
    
    if filename != 'classif_genes_':
        filename1 = filename
    else:
        filename1 = filename + str(thres * 100) + '_' + str(thres2 * 100) + '.csv'
    
    df.to_csv(filename1)
    return(df)

def stats_barplot(csv_file = 'classif_genes_90.0_70.0.csv', tissue_name = 'Tissue'):
    df = pd.read_csv(csv_file)
    ntotal = df['Type'].count()
    
    pnexp = df.loc[df['Type'] == 'NotExpressed']['Type'].count() / float(ntotal)
    pmono = df.loc[df['Type'] == 'Monoform']['Type'].count() / float(ntotal)
    pbi = df.loc[df['Type'] == 'Biform']['Type'].count() / float(ntotal)
    ptri = df.loc[df['Type'] == 'Triform']['Type'].count() / float(ntotal)
    pmulti = df.loc[df['Type'] == 'Multiform']['Type'].count() / float(ntotal)
    
    stats = [['NotExpressed', 'Monoform', 'Biform', 'Triform', 'Multiform'], [pnexp, pmono, pbi, ptri, pmulti]]
    
    fig, ax = plt.subplots()
    ntypes = len(stats[1])    
    acum = stats[1][0]
    ax.bar(1, acum, 1, label = stats[0][0])
    for i in range(1, ntypes):
        ax.bar(1, stats[1][i], 1, bottom = acum, label = stats[0][i])
        acum =  acum + stats[1][i]
    
    ax.set_ylabel('Proportion')
    ax.set_title('Isoform distribution of ' + tissue_name)
    ax.set_xticks([0])
    ax.legend()
    plt.subplots_adjust(bottom = 0.1, top = 0.9)
    plt.show()
    
    return stats

def notexp_stats(csv_file = 'classif_genes_90.0_70.0.csv', tissue_name = 'Tissue'):
    df = pd.read_csv(csv_file)
    ntotal = df['Type'].count()
    
    expr = df.loc[df['Type'] != 'NotExpressed']['Type'].count() / float(ntotal)
    notexp = df.loc[df['Type'] == 'NotExpressed']['Type'].count() / float(ntotal)
    
    stats = [['Expressed', 'NotExpressed'], [expr, notexp]]
    
    fig1, ax1 = plt.subplots()
    ax1.pie(stats[1], labels = stats[0], autopct = '%1.1f%%', shadow = True, startangle = 0)
    ax1.set_title('Isoform distribution of ' + tissue_name + '\nExpressed and not expressed isoforms')
    ax1.axis('equal')
    plt.show()
    
    return stats

def statistics(csv_file = 'classif_genes_90.0_70.0.csv', tissue_name = 'Tissue'):
    df = pd.read_csv(csv_file)
    nmono = df.loc[df['Type'] == 'Monoform']['Type'].count()
    nbi = df.loc[df['Type'] == 'Biform']['Type'].count()
    ntri = df.loc[df['Type'] == 'Triform']['Type'].count()
    nmulti = df.loc[df['Type'] == 'Multiform']['Type'].count()
    ntotal = nmono + nbi + ntri + nmulti
    
    pmono = nmono / float(ntotal)
    pbi = nbi / float(ntotal)
    ptri = ntri / float(ntotal)
    pmulti = nmulti / float(ntotal)
    
    stats = [['Monoform', 'Biform', 'Triform', 'Multiform'], [pmono, pbi, ptri, pmulti]]
    
    fig1, ax1 = plt.subplots()
    ax1.pie(stats[1], labels = stats[0], autopct = '%1.1f%%', shadow = True, startangle = 0)
    ax1.set_title('Isoform distribution of ' + tissue_name + '\nExpressed isoforms')
    ax1.axis('equal')
    plt.show()
    
    return stats

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