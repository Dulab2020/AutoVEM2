#=================================================================================================================
# AutoVEM2 V1.0
# Author: Xi Binbin
# Email: 201766841276@mail.scut.edu.cn
#=================================================================================================================
import os
import sys
import time
import numpy as np
import pandas as pd 
import call
import matplotlib.pyplot as plt 
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.patches as mpatches

def snp_filter(file, directory, length, sites=None, fre=None):
    '''
    filter SNP sites

    :param file: snp_merged.tsv
    :param directory: output directory
    :param sites: snp sites that you are interested
    :param length: int, length of reference genome sequence
    :param fre: frequency of snp sites that more than fre will be obtained, default 0.05
    :returns snp_pos, snp_ref_alt(dict)
    '''
    prefix = directory[:]
    suffix = "snp_sites.tsv"
    snp_sites = os.path.join(directory, 'snp_sites.tsv')
    if os.path.exists(snp_sites):
        print("%s already exists."%snp_sites)
        i = 0
        while True:
            i = i + 1
            temp = os.path.join(prefix, r'#%s.%s#'%(suffix,i))
            if os.path.exists(temp):
                continue
            else :
                os.rename(snp_sites, temp)
                print("Back up %s to %s"%(snp_sites, temp))
                break


    df = pd.read_csv(file, sep='\t')
    ids = df['Id'].unique().tolist()
    n_genome = len(ids)
    #print(n_genome)


    counts = df['Position'].value_counts()
    frequency = counts/n_genome
    frequency = frequency.round(decimals=4)
    frequency = frequency.sort_index(ascending=True)
    snp_dict = dict(zip(frequency.index.tolist(), frequency.values.tolist()))
    for item in list(range((length+1))):
        if item not in snp_dict:
            snp_dict[item] = 0.0
    if 0 in snp_dict:
        snp_dict.pop(0)


    snp_pos = list()    # important sites
    if sites is None:
        if fre is None:
            fre = 0.05
        for key,item in snp_dict.items():
            if item>=fre:
                snp_pos.append(key)
    else:
        if fre is None:
            for item in sites:
                if item in snp_dict:
                    snp_pos.append(item)       
        else:
            for item in sites:
                if ((item in snp_dict) and (snp_dict[item]>=fre)):
                    snp_pos.append(item)
                
    if len(snp_pos)<=1:
        print('There are no or too little sites that meet your requirements(<2).')
        sys.exit()
    print("The following %s mutations meet you requirements:"%len(snp_pos), end=' ')
    for item in snp_pos:
        print(item, end=' ')
    print()
    print("For more information about these sites, please view the %s file"%snp_sites)
    time.sleep(5)


    snp_pos.sort()
    snp_ref_alt = dict()
    for item in snp_pos:
        df2 = df[df['Position']==item]
        if len(df2.index)==0:
            Ref = "N"
            Alt = "N"
        else:
            Ref = df2['Ref'].value_counts().index.tolist()[0]
            Alt = df2['Alt'].value_counts().index.tolist()[0]
        snp_ref_alt[item] = [Ref, Alt, snp_dict[item]]
    header = 'Position\tRef\tAlt\tFrequency\n'
    with open(snp_sites, 'a') as fhand:
        fhand.write(header)
        for key, item in snp_ref_alt.items():
            record = str(key)+'\t'+str(item[0])+'\t'+str(item[1])+'\t'+str(item[2])+'\n'
            fhand.write(record)


    n_snp = len(snp_ref_alt) 
    if n_snp>=500:
        print('There are too many sites meet your requirements (there are %s, should <500).'%n_snp)
        print('Please increase the frequency to filter out more sites or specify no more than 500 sites')
        sys.exit()


    return snp_pos, snp_ref_alt

def ref_haplotype(position, refsequence):
    '''
    reference haplotype sequence

    :param position：snp positions, int
    :param refsequence: reference genome sequence
    :returns referenceHaplotype
    '''
    with open(refsequence, 'r') as fhand:
        referenceGenomeSequence = str()
        for line in fhand.readlines():
            line = line.strip()
            line = line.replace(' ','')
            if len(line)==0:
                continue
            if line[0] == '>':
                continue
            else :
                referenceGenomeSequence = referenceGenomeSequence + line

    referenceHaplotype = str()
    for item in position:
        referenceHaplotype = referenceHaplotype + referenceGenomeSequence[item-1]
    return referenceHaplotype

def genome_haplotype(file, position, positionAlt, referenceHaplotype, directory):
    '''
    get haplotype sequence of genome

    :param file: snp_merged.tsv
    :param position: mutation positions
    :param positionAlt: position_ref_alt(dict)
    :param referenceHaplotype: reference haplotype sequence
    :param directory: output directory
    :returns filePath: data.tsv
    '''
    filePath = os.path.join(directory, 'data.tsv')
    prefix = directory[:]
    suffix = 'data.tsv'
    if os.path.exists(filePath):
        print("%s already exists."%filePath)
        i = 0
        while True:
            i = i + 1
            tmp = os.path.join(prefix,r"#%s.%s#"%(suffix, i))
            if os.path.exists(tmp):
                continue
            else:
                os.rename(filePath, tmp)
                print("Back up %s to %s"%(filePath, tmp))
                break

    df = pd.read_csv(file, sep='\t')
    ids = df['Id'].unique().tolist()

    Date = list()
    Country = list()
    Case_id = list()
    snp_positions = list()

    for item in ids:
        df1 = df[df['Id']==item]
        case = item
        date = df1['Date'].value_counts().index.tolist()
        if len(date)==0:
            date = "NA"
        else:
            date = date[0]
        country = df1['Country'].value_counts().index.tolist()
        if len(country)==0:
            country = "NA"
        else:
            country = country[0]
        Case_id.append(case)
        Date.append(date)
        Country.append(country)
        mutation_positions = df1['Position'].tolist()
        mutation_positions = list([int(x) for x in mutation_positions])
        snp_positions.append(mutation_positions)

    haplotypes = list()
    for item in snp_positions:
        sites = str()
        i = -1
        for site in position:
            i = i + 1
            site = int(site)
            if site in item:
                sites = sites + positionAlt[site][1]
            else:
                sites = sites + referenceHaplotype[i]
        haplotypes.append(sites)
    data = pd.DataFrame(data={'Id': Case_id,
                                   'Date': Date,
                                   'Country':Country,
                                   'Hap': haplotypes})
    data.to_csv(filePath, sep='\t', index=False)

    
    return filePath

def block_file(position, directory):
    '''
    block.txt

    :param position: snp_pos
    :param directory: output directory
    :returns blockFile：block.txt
    '''
    num = len(position)
    blockFile = os.path.join(directory, 'block.txt')
    prefix = directory[:]
    suffix = 'block.txt'
    if os.path.exists(blockFile):
        print("%s already exists."%blockFile)
        i = 0
        while True:
            i = i + 1
            tmp = os.path.join(prefix,r"#%s.%s#"%(suffix, i))
            if os.path.exists(tmp):
                continue
            else:
                os.rename(blockFile, tmp)
                print("Back up %s to %s"%(blockFile, tmp))
                break
    with open(blockFile, 'a') as fhand:
        for i in list(range(num)):
            fhand.write(str(i+1))
            fhand.write('\t')
    
    return blockFile

def map_file(position, directory):
    '''
    snp.info

    :param position: snp_ref_alt
    :returns mapFile: snp.info
    '''
    mapFile = os.path.join(directory, 'snp.info')
    prefix = directory[:]
    suffix = 'snp.info'
    if os.path.exists(mapFile):
        print("%s already exists."%mapFile)
        i = 0
        while True:
            i = i + 1
            tmp = os.path.join(prefix,r"#%s.%s#"%(suffix, i))
            if os.path.exists(tmp):
                continue
            else:
                os.rename(mapFile, tmp)
                print("Back up %s to %s"%(mapFile, tmp))
                break
    with open(mapFile, 'a') as fhand:
        for key, item in position.items():
            name = str(item[0])+str(key)+str(item[1])
            record = name + '\t' + str(key) + '\n'
            fhand.write(record)
    
    return mapFile

def ped_file(file, directory):
    '''
    生成ped格式文件

    :param file: data.tsv
    :param directory: output directory
    :returns pedFile: snp.ped
    '''
    pedFile = os.path.join(directory, 'snp.ped')
    prefix = directory[:]
    suffix = 'snp.ped'
    if os.path.exists(pedFile):
        print("%s already exists."%pedFile)
        i = 0
        while True:
            i = i + 1
            tmp = os.path.join(prefix,r"#%s.%s#"%(suffix, i))
            if os.path.exists(tmp):
                continue
            else:
                os.rename(pedFile, tmp)
                print("Back up %s to %s"%(pedFile, tmp))
                break
    with open(pedFile, 'a') as f:
        with open(file, 'r') as fhand:
            i = 0
            for line in fhand.readlines():
                i = i + 1
                if i == 1:
                    continue
                line = line.rstrip()
                line = line.split('\t')
                record = line[0]+'\t'+line[0]+'\t'+'0'+'\t'+'0'+'\t'+'0'+'\t'+'0'
                for item in line[-1]:
                    record = record + '\t' + item + '\t' + item
                record = record + '\n'
                f.write(record)
    
    return pedFile

def linkage_analysis(ped, mapf, block, directory):
    '''
    linkage analysis

    :param ped: snp.ped
    :param mapf: snp.info
    :param block: block.txt
    :param directory: output directory
    :returns haplotypesFile: plot.CUSTblocks
    '''

    temp = os.path.join(directory, 'plot')
    haplotypesFile = os.path.join(directory, 'plot.CUSTblocks')
    prefix = directory[:]
    suffix = 'plot.CUSTblocks'
    if os.path.exists(haplotypesFile):
        print("%s already exists."%haplotypesFile)
        i = 0
        while True:
            i = i + 1
            tmp = os.path.join(prefix,r"#%s.%s#"%(suffix, i))
            if os.path.exists(tmp):
                continue
            else:
                os.rename(haplotypesFile, tmp)
                print("Back up %s to %s"%(haplotypesFile, tmp))
                break

    os.system('java -jar %s -n -skipcheck -pedfile %s -info %s -blocks %s -png -out %s' % (call.haploview, ped, mapf, block, temp))

    return haplotypesFile

def haplotyper(file, dataFile, directory):
    '''
    hyplotype

    :param file: plot.CUSTblocks
    :param dataFile: data.tsv
    :param directory: output directory
    :returns dataPlot, haplotypes: data_plot.tsv
    '''
    dataPlot = os.path.join(directory, 'data_plot.tsv')
    prefix = directory[:]
    suffix = 'data_plot.tsv'
    if os.path.exists(dataPlot):
        print("%s already exists."%dataPlot)
        i = 0
        while True:
            i = i + 1
            tmp = os.path.join(prefix,r"#%s.%s#"%(suffix, i))
            if os.path.exists(tmp):
                continue
            else:
                os.rename(dataPlot, tmp)
                print("Back up %s to %s"%(dataPlot, tmp))
                break
    haplotypes = os.path.join(directory, 'haplotypes.tsv')
    prefix = directory[:]
    suffix = 'haplotypes.tsv'
    if os.path.exists(haplotypes):
        print("%s already exists."%haplotypes)
        i = 0
        while True:
            i = i + 1
            tmp = os.path.join(prefix,r"#%s.%s#"%(suffix, i))
            if os.path.exists(tmp):
                continue
            else:
                os.rename(haplotypes, tmp)
                print("Back up %s to %s"%(haplotypes, tmp))
                break
    nt_dict = {'1': 'A', '2': 'C', '3': 'G', '4': 'T'}
    
    hap_dict = dict()
    with open(file, 'r') as fhand:
        i = -1
        for line in fhand.readlines():
            sequence = ''
            i = i + 1
            if i == 0:
                continue 
            hap = 'H' + str(i)
            sequence_list = line.split()
            sequence_num = sequence_list[0]
            for num in sequence_num:
                sequence = sequence + nt_dict[num]
            hap_dict[sequence] = hap

    header = 'Id\tDate\tCountry\tHap\tName\n'
    with open(dataPlot, 'a') as f:
        f.write(header)
        with open(dataFile, 'r') as fhand:
            i = 0
            for line in fhand.readlines():
                i = i + 1
                if i == 1:
                    continue
                line = line.rstrip()
                line = line.split('\t')
                hap = line[-1]
                record = line[0]+'\t'+line[1]+'\t'+line[2]+'\t'+line[3]+'\t'
                if hap in hap_dict:
                    name = hap_dict[hap]
                else:
                    name = 'other'
                record = record + name + '\n'
                f.write(record)
    
    df = pd.read_csv(dataPlot, sep='\t')
    num_case = len(df.index)
    haps = df.Name.unique().tolist()

    haps_sorted = []
    for i in list(range(len(haps)-1)):
        h = 'H' + str(i+1)
        haps_sorted.append(h)
    haps_sorted.append("other")
    print(haps_sorted)

    header = 'Name\tSequence\tFrequency\n'
    with open(haplotypes, 'a') as fhand:
        fhand.write(header)
        for item in haps_sorted:
            df2 = df[df.Name==item]
            number = len(df2.index)
            frequency = number/num_case
            frequency = round(frequency, 4)
            seq = df2.Hap.unique().tolist()[0]
            if item == "other":
                seq = "NA"
            record = str(item)+'\t'+str(seq)+'\t'+str(frequency)+'\n'
            fhand.write(record)

    return dataPlot

def hap_plot(datafile, snpfile, snp_position, directory, length):
    '''
    :param datafile: data_plot.tsv
    :param snpfile: snp_merged.tsv
    :param snp_position: mutation position
    :param directory: output directory
    :param length: reference genome length
    '''
    df = pd.read_csv(datafile, sep='\t')
    haps_tmp = df.Name.unique().tolist()
    num = len(haps_tmp)
    haps = list()
    i = 0
    for x in haps_tmp:
        i = i + 1
        if i == num:
            break
        else :
            item = 'H'+ str(i)
        haps.append(item)
    
    haps_id = dict()
    for item in haps:
        df2 = df[df.Name==item]
        id_list = df2.Id.tolist()
        haps_id[item] = id_list

    df = pd.read_csv(snpfile, sep='\t')
    data = dict()
    for hap,ids in haps_id.items():
        df2 = df.loc[df.Id.isin(ids)]
        num = len(df2.Id.unique().tolist())
        tmp = df2.Position.value_counts()
        fre = tmp/num
        fre = fre.round(decimals=4)
        fre = fre[fre.values>0.0001]
        fre = fre.sort_index(ascending=True)
        pos = fre.index.tolist()
        frequency = fre.values.tolist()
        # pos_fre = dict(zip(pos, frequency))

        # pos_fre_dict = dict()
        # for item in list(range(length+1)):
        #     if item not in pos_fre:
        #         pos_fre_dict[item] = 0.0
        #     else:
        #         pos_fre_dict[item] = pos_fre[item]
        # pos = list()
        # frequency = list()
        # for key, value in pos_fre_dict.items():
        #     pos.append(key)
        #     frequency.append(value)
        if pos[0] == 0:
            del pos[0]
            del frequency[0]
        data[hap] = [pos,frequency]


    flag = len(data)
    positions = snp_position[:]
    lx = len(positions)
    pdf_file = os.path.join(directory, 'hap_mutations.pdf')
    prefix = directory[:]
    suffix = 'hap_mutations.pdf'
    if os.path.exists(pdf_file):
        print("%s already exists."%pdf_file)
        i = 0
        while True:
            i = i + 1
            tmp = os.path.join(prefix,r"#%s.%s#"%(suffix, i))
            if os.path.exists(tmp):
                continue
            else:
                os.rename(pdf_file, tmp)
                print("Back up %s to %s"%(pdf_file, tmp))
                break
    pdf = PdfPages(pdf_file)
    fig, axes = plt.subplots(flag,1,figsize=(10,5*flag))
    i = -1
    for hap,values in data.items():
        if flag==1:
            ax = axes
        else:
            i = i + 1
            ax = axes[i]
        x = values[0]
        height = values[1]
        colors = list()
        for pos in x:
            if pos in positions:
                colors.append('blueviolet')
            else:
                colors.append('tomato')
        ax.set_xlim((int(-0.05*length),int(1.05*length)))
        ax.set_ylim((-0.05,1.05))
        ax.set_ylabel(str(hap))
        ax.set_xlabel('Position')
        y_ticks = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
        ax.set_yticks(y_ticks)
        ax.set_yticklabels(y_ticks)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.set_facecolor('whitesmoke')
        if i==0:
            blueviolet_patch = mpatches.Patch(color='blueviolet', label='the %s specific sites'%lx)
            tomato_patch = mpatches.Patch(color='tomato', label='other mutation sites')
            ax.legend(handles=[blueviolet_patch, tomato_patch], frameon=False, bbox_to_anchor=(0,1), ncol=2, loc="lower left")
        ax.bar(x=x, height=height,color=colors,width=0.5)
    pdf.savefig()
    pdf.close()
    plt.close()

def module2(file, directory, refsequence, sites=None, frequency=None):
    '''
    mutation analysis

    :param file: snp_merged.tsv file produced by AutoVEM2 call
    :param directory: output directory
    :param refsequence: reference genome sequence
    :param sites: mutation sites that you are interested
    :param frequecy: mutations with mutation frequecy lower than the given frequency will be filtered out, default 0.05
    '''
    # print(sites)
    # print(directory)
    # print(refsequence)
    # print(sites)
    # print(frequency)
    if not os.path.exists(os.path.abspath(file)):
        print("Error: can't find the %s file"%(os.path.abspath(file)))
        sys.exit()
    if not os.path.exists(os.path.abspath(refsequence)):
        print("Error: can't find the %s file"%(os.path.abspath(file)))
        sys.exit()
    
    path = os.path.abspath(directory)
    if os.path.exists(path):
        print("%s already exists."%path)
        dirname, filename = os.path.split(path)
        i = 0
        while True:
            i = i + 1
            temp = os.path.join(dirname, r"#%s.%s#"%(filename,i))
            if os.path.exists(temp):
                continue
            else :
                os.rename(path, temp)
                print("Back up %s to %s"%(path, temp))
                break
    os.mkdir(path)


    with open(refsequence, 'r') as fhand:
        reference = ''
        for line in fhand.readlines():
            line = line.strip()
            line = line.replace(' ','')
            if len(line)==0:
                continue
            if line[0]==">":
                continue
            else:
                reference = reference + line
    length_genome = len(reference)
    # print(length_genome)
    if length_genome == 0:
        print('The %s file has no genome sequence'%refsequence)
        sys.exit()

    print('Obtaining specific sites...')
    snp_position, snp_ref_alt = snp_filter(file, directory, length_genome, sites=sites, fre=frequency)
    print('Done!')

    print('Obtaining reference haplotype sequence and data.tsv file...')
    ref_haplotype_sequence = ref_haplotype(snp_position, refsequence)
    dataFile = genome_haplotype(file, snp_position, snp_ref_alt, ref_haplotype_sequence, directory)
    print('Done!')

    print('Obtaining block.txt file...')
    blockFile = block_file(snp_position, directory)
    print('Done!')

    print('Obtaining map file...')
    mapFile = map_file(snp_ref_alt, directory)
    print('Done!')

    print('Obtaining snp.ped file...')
    pedFile = ped_file(dataFile, directory)
    print('Done!')

    print('Linkage analyzing...')
    haplotypesFile = linkage_analysis(pedFile, mapFile, blockFile, directory)
    os.remove(mapFile)
    os.remove(pedFile)
    os.remove(blockFile)
    print('Done!')
    
    print('Haplotyping......')
    data_plot = haplotyper(haplotypesFile, dataFile, directory)
    os.remove(haplotypesFile)
    os.remove(dataFile)
    print('Done!')

    print('Plotting hap_mutations...')
    hap_plot(data_plot, file, snp_position, directory, length_genome)
    print('Done!')

    return data_plot
