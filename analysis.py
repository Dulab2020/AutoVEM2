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

def snp_filter(file, directory, sites=None, fre=None):
    '''
    filter SNP sites

    :param file: the absolute path of the snp_merged.tsv
    :param directory: the absolute path of the output directory
    :param sites: snp sites that of interest
    :param fre: frequency of snp sites that more than fre will be obtained, default 0.05
    :returns snp_pos, snp_ref_alt(dict)
    '''
    ## storage the retained SNP sites
    ## check and create the file
    snp_sites = os.path.join(directory, 'snp_sites.tsv')
    if os.path.exists(snp_sites):
        os.system(f"rm -rf {snp_sites}")
    ## read the SNV files
    df = pd.read_csv(file, sep='\t')
    ## count the number of sequences
    ids = df['Id'].unique().tolist()
    n_genome = len(ids)
    ## calculate the mutation frequency at each position
    counts = df['Position'].value_counts()
    frequency = counts/n_genome
    frequency = frequency.round(decimals=4)
    frequency = frequency.sort_index()
    frequency = frequency[frequency.index!=0]
    snp_dict = dict(zip(frequency.index.tolist(), frequency.values.tolist()))
    ## get the filtered SNP positions
    snp_pos = list()
    ## sites with mutation frequency bigger than the given frequency or default frequency (0.05) will be retained.
    if sites is None:
        if fre is None:
            snpSites = frequency[frequency.values>=0.05]
            snp_pos = snpSites.index.tolist()
        else:
            snpSites = frequency[frequency.values>=fre]
            snp_pos = snpSites.index.tolist()
    ## if provided sites, these sites will be retained
    else:
        snp_pos = sites
    ## whether the number of sites bigger than 1 and less than 500
    if len(snp_pos)<=1:
        print('There are no or too little sites that meet your requirements.')
        sys.exit()
    elif len(snp_pos)>=500:
        print(f'There are too many sites meet your requirements (should be less than 500).')
        print('Please increase the frequency to filter out more sites or specify no more than 500 sites')
        sys.exit()
    else:
        snp_pos.sort()
        ## print the site retained
        out_line = ""
        for site in snp_pos:
            out_line = out_line + " " + str(site)
        print(f"The following sites will be retained: {out_line}")

        ## get the mutation information of the retained sites
        snp_ref_alt = dict()
        for snp in snp_pos:
            df2 = df[df['Position']==snp]
            Ref = df2['Ref'].value_counts().index.tolist()[0]
            Alt = df2['Alt'].value_counts().index.tolist()[0]
            snp_ref_alt[snp] = (Ref, Alt, snp_dict[snp])
        ## write the information of retained sites to the record file
        header = 'Position\tRef\tAlt\tFrequency\n'
        with open(snp_sites, 'a') as fhand:
            fhand.write(header)
            for pos, (Ref, Alt, Fre) in snp_ref_alt.items():
                record = str(pos) + '\t' + Ref + '\t' + Alt + '\t' + str(Fre) + '\n'
                fhand.write(record)
            
        return snp_pos, snp_ref_alt

def ref_haplotype(position, refsequence):
    '''
    reference haplotype sequence

    :param position: snp positions, int
    :param refsequence: reference genome sequence
    :returns referenceHaplotype
    '''
    with open(refsequence, 'r') as fhand:
        referenceGenomeSequence = str()
        for line in fhand.readlines():
            line = line.strip()
            line = line.replace(' ','')
            line = line.replace('\t','')
            if len(line)==0:
                pass
            elif line[0] == '>':
                pass
            else:
                referenceGenomeSequence = referenceGenomeSequence + line

    referenceHaplotype = str()
    for pos in position:
        referenceHaplotype = referenceHaplotype + referenceGenomeSequence[pos-1]
    return referenceHaplotype

def genome_haplotype(file, positions, referenceHaplotype, snp_alt, directory):
    '''
    get haplotype sequence of genome

    :param file: snp_merged.tsv
    :param positions: mutation positions
    :param referenceHaplotype: reference haplotype sequence
    :param snp_alt: the snp_ref_alt_fre dict, {pos:[ref, alt, fre]}
    :param directory: output directory
    :returns filePath: data.tsv
    '''
    ## the file stores the haplotype sequence of each sequence
    filePath = os.path.join(directory, 'data.tsv')
    if os.path.exists(filePath):
        os.system(f"rm -rf {filePath}")
    ## read the mutation file
    df = pd.read_csv(file, sep='\t')
    ## used to store the meta information of the sequences
    Date = []
    Country = []
    Case_id = []
    snp_positions = []
    Haplotypes = []
    ## group according to the Id columns
    group = df.groupby("Id")
    j = 0
    for idx, df_s in group:
        ## print processing information
        j = j + 1
        print(f"[INFO]: Obtaining the haplotype sequence of the {j}th sequence: {idx}")
        date = df_s['Date'].value_counts().index.tolist()[0]
        country = df_s['Country'].value_counts().index.tolist()[0]
        mutation_positions = set(df_s["Position"].tolist())
        hap_seq = ""
        for i, snp in enumerate(positions):
            if snp in mutation_positions:
                hap_seq = hap_seq + snp_alt[snp][1]
            else:
                hap_seq = hap_seq + referenceHaplotype[i]

        Case_id.append(idx)
        Date.append(date)
        Country.append(country)
        Haplotypes.append(hap_seq)

    data_df = pd.DataFrame(data={'Id': Case_id,
                                   'Date': Date,
                                   'Country':Country,
                                   'Hap': Haplotypes})
    data_df.to_csv(filePath, sep='\t', index=False)

    return filePath

def block_file(positions, directory):
    '''
    block.txt

    :param positions: snp_pos
    :param directory: output directory
    :returns blockFile: block.txt
    '''
    num = len(positions)
    blockFile = os.path.join(directory, 'block.txt')
    if os.path.exists(blockFile):
        os.system(f"rm -rf {blockFile}")
    ## write the output block
    with open(blockFile, 'a') as fhand:
        for i in list(range(num)):
            fhand.write(str(i+1))
            fhand.write('\t')
    
    return blockFile

def map_file(positions, directory):
    '''
    snp.info

    :param positions: snp_ref_alt
    :returns mapFile: snp.info
    '''
    mapFile = os.path.join(directory, 'snp.info')
    if os.path.exists(mapFile):
        os.system(f"rm -rf {mapFile}")
    ## write the maker to the snp.info file
    with open(mapFile, 'a') as fhand:
        for key, item in positions.items():
            name = str(item[0])+str(key)+str(item[1])
            record = name + '\t' + str(key) + '\n'
            fhand.write(record)
    
    return mapFile

def ped_file(file, directory):
    '''
    snp.ped

    :param file: data.tsv
    :param directory: output directory
    :returns pedFile: snp.ped
    '''
    ## the snp.ped file
    pedFile = os.path.join(directory, 'snp.ped')
    if os.path.exists(pedFile):
        os.system(f"rm -rf {pedFile}")

    df = pd.read_table(file)
    Ids = df.Id.tolist()
    HaplotypeSequence = df.Hap.tolist()
    with open(pedFile, 'a') as f:
        record = ""
        constantString = '\t0\t0\t0\t0\t'
        for i,idx in enumerate(Ids):
            hap = ""
            hap_seq = HaplotypeSequence[i]
            genotype = list()
            for base in hap_seq:
                genotype.append(base)
                genotype.append(base)
            hap = "\t".join(genotype)
            record = str(i) + "\t" + str(idx) + constantString + hap + "\n"
            f.write(record)
            record = ""

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
    if directory[-1] == "/":
        temp = directory + "plot"
    else:
        temp = directory + "/plot"
    haplotypesFile = os.path.join(directory, 'plot.CUSTblocks')
    if os.path.exists(haplotypesFile):
        os.system(f"rm -rf {haplotypesFile}")

    os.system(f'java -jar {call.haploview} -n -skipcheck -pedfile {ped} -info {mapf} -blocks {block} -png -out {temp}')

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
    if os.path.exists(dataPlot):
        os.system(f"rm -rf {dataPlot}")
    
    haplotypes = os.path.join(directory, 'haplotypes.tsv')
    if os.path.exists(haplotypes):
        os.system(f"rm -rf {haplotypes}")
        
    nt_dict = {'1': 'A', '2': 'C', '3': 'G', '4': 'T'}
    ## obtain the haplotype sequence
    hap_dict = dict()
    with open(file, 'r') as fhand:
        for i, line in enumerate(fhand.readlines()):
            sequence = ''
            if i == 0:
                pass
            else:
                hap = 'H' + str(i)
                sequence_num = line.split()[0]
                proportion = float(line.split()[1][1:-1])
                ## if proportion less than 0.01, it means the variant population is too small
                ## so this variant population will be filtered out
                if proportion < 0.01:
                    break
                else:
                    ## get the haplotype sequence and name
                    for num in sequence_num:
                        sequence = sequence + nt_dict[num]
                    hap_dict[sequence] = [hap, proportion]

    df = pd.read_table(dataFile)
    Name = list()
    for hapSeq in df.Hap.tolist():
        if hapSeq in hap_dict:
            Name.append(hap_dict[hapSeq][0])
        else:
            Name.append("other")
    df["Name"] = Name
    df.to_csv(dataPlot, sep="\t", index=None)
    
    headers = "Name\tSequence\tFrequency\n"
    with open(haplotypes, 'a') as fhand:
        fhand.write(headers)
        for key, value in hap_dict.items():
            tem = ''
            tem = value[0] + '\t' + key + '\t' + str(value[1]) + "\n"
            fhand.write(tem)

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
        df2 = df[np.isin(df.Id,ids)]
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

    :param file: the absolute path of the snp_merged.tsv file produced by AutoVEM2 call
    :param directory: the absolute path of the output directory
    :param refsequence: the absolute path of the reference genome sequence
    :param sites: mutation sites that of interest
    :param frequecy: mutations with mutation frequecy lower than the given frequency will be filtered out, default 0.05
    '''
    ## get the reference genome sequence and the length of the reference genome sequence
    reference = ''
    with open(refsequence, 'r') as fhand:
        for line in fhand.readlines():
            line = line.strip()
            line = line.replace(' ','')
            line = line.replace("\t", "")
            if len(line)==0:
                continue
            elif line[0]==">":
                continue
            else:
                reference = reference + line
    length_genome = len(reference)
    ## check whether the reference genome sequence exists
    if length_genome == 0:
        print(f'The {refsequence} file has no genome sequence')
        sys.exit()
    ## check whether the values of the --sites is valid
    if sites is None:
        pass
    else:
        sites.sort()
        minimum = sites[0]
        if minimum<=0 or sites[-1]>length_genome:
            print("Value of the --sites argument should be bigger than 0 and no bigger than the length of the reference genome.")
            sys.exit()
    ## obtaining the specific sites
    print('Obtaining specific sites...')
    snp_position, snp_ref_alt = snp_filter(file, directory, sites=sites, fre=frequency)
    print('Done!')
    ## get the reference haplotype sequence
    ## get the haplotype sequence of each genome seqeunce
    print('Obtaining the reference haplotype sequence ...')
    ref_haplotype_sequence = ref_haplotype(snp_position, refsequence)
    print("Done!")
    print("Obtaining the haplotype sequence of each genome sequence")
    dataFile = genome_haplotype(file, snp_position, ref_haplotype_sequence, snp_ref_alt, directory)
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
    print('Done!')
    
    print('Haplotyping......')
    data_plot  = haplotyper(haplotypesFile, dataFile, directory)
    print('Done!')
    try:
        os.remove(mapFile)
        os.remove(pedFile)
        os.remove(blockFile)
        os.remove(haplotypesFile)
        os.remove(dataFile)
        print("Having done haplotyping and linkage analysis!")
    except:
        print("Having done haplotyping and linkage analysis!")

    return data_plot
