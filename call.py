#=================================================================================================================
# AutoVEM2
# Author: Xi Binbin
# Email: 201766841276@mail.scut.edu.cn
#=================================================================================================================
import os
import sys
import shutil
import copy

softDirectory = os.path.dirname(os.path.abspath(sys.argv[0]))   # AutoVEM2 directory
ploidy = os.path.join(softDirectory, 'ploidy.txt')  # ploidy.txt file, used in call function
haploview = os.path.join(softDirectory, 'Haploview.jar')    # Haploview tool
logging = os.path.join(softDirectory, 'out.log')

def extract_sequence(directory, genomeDirectory):
    '''
    :param directory: input directory
    :param genomeDirectory: genome directory
    :return filesPath
    '''
    files = os.listdir(directory)
    files_list = list()
    for file in files:
        file = os.path.join(directory, file)
        if os.path.isfile(file):
            files_list.append(file)


    pointer = 1
    for file in files_list:
        with open(file, 'r') as fhand:
            id_path = str()
            for line in fhand.readlines():
                line = line.rstrip()
                if len(line)==0:
                    continue
                if line[0] == '>':
                    pointer = 0
                    tmp = line[:]
                    line = line.split('|')
                    if len(line) != 4:
                        pointer = 1
                        continue
                    id = line[1]
                    id = id.split()[0]
                    id = id + '.fa'
                    id_path = os.path.join(genomeDirectory, id)
                    if os.path.exists(id_path):
                        pointer = 1
                        continue
                    else :
                        with open(id_path, 'a') as f:
                            f.write(tmp+'\n')
                            pointer = 0
                            continue

                if pointer == 0:
                    with open(id_path, 'a') as f:
                        f.write(line+'\n')
                        continue
                else:
                    continue

    filesPath = list()
    filesList = os.listdir(genomeDirectory)
    for path in filesList:
        path = os.path.join(genomeDirectory, path)
        filesPath.append(path)

    return filesPath

def create_index(path, refsequence):
    '''
    create index
    :param path: index directory
    :param refsequence: fasta fomat reference sequence file
    '''
    os.system('bowtie2-build -f %s %s/index' % (refsequence, path))

def genome_quality_control(file, referenceLength=None, a=None, b=None, c=None):
    '''
    quality control for the first time

    :param file: a fasta format file of a genome sequence
    :param a: length
    :param b: number of unknown bases, number_n
    :param c: number of degenerate bases, number_db
    :returns flag
    '''
    if not os.path.getsize(file):
        return -1

    with open(file, 'r') as fhand:
        genomeSequence = str()
        for line in fhand.readlines():
            line = line.strip()
            line = line.replace(" ","")
            if len(line)==0:
                continue
            if line[0]=='>':
                continue
            genomeSequence = genomeSequence + line


    len_sequence = len(genomeSequence)
    number_n = 0
    number_db = 0
    letters = ['A', 'T', 'G', 'C', 'N']
    for item in genomeSequence:
        item = item.capitalize()
        if item == 'N':
            number_n = number_n + 1
        elif item not in letters:
            number_db = number_db + 1
        else :
            pass
    
    if a is None:
        pass
    else :
        if a > 1:
            if len_sequence < a:
                return -1
        else:
            if (len_sequence/referenceLength) < a:
                return -1

    if b is None:
        pass
    else:
        if b > 1:
            if number_n > b:
                return -1
        else:
            if (number_n/referenceLength) > b:
                return -1

    if c is None:
        pass
    else:
        if c > 1:
            if number_db > c:
                return -1
        else:
            if (number_db/referenceLength) > c:
                return -1  
  
    return 0

def split_sequence(file, path, filter="yes"):
    '''
    split genome sequence to reads(~80nt)

    :param file: genome sequence
    :param path: output directory
    :param filter: wether filter out sequences with unclear collection time, default yes
    :returns bassicMessage(dict), splitFile(file), tempAnalysisDirectory(temp directory):
    '''
    China = ['taiwan', 'guangzhou', 'fujian', 'sichuan','wuhan', 'hangzhou', 'jiangsu', 
             'shanghai','shandong', 'guangdong', 'foshan', 'nanchang', 'yingtan','hunan', 
             'shangrao', 'yichun', 'beijing','jingzhou','pingxiang', 'shenzhen', 'lishui', 
             'zhejiang', 'fuzhou', 'shaoxing', 'xinyu','yunnan', 'jiujiang', 'chongqing', 
             'henan', 'hefei', 'fuyang', 'changzhou', 'jiangxi', 'ganzhou', 'hong kong', 
             'macau', 'qingdao', 'liaoning', 'harbin', 'tianmen', 'jian']   
    deletes = ['lion', 'cat', 'env', 'canine', 'tiger', 'mink']


    with open(file, 'r') as fhand:
        head = str()
        sequence = str()
        for line in fhand.readlines():
            line = line.rstrip()
            if len(line)==0:
                continue
            if line[0]==">":
                header = line[:]
                continue
            line = line.replace(" ","")
            sequence = sequence + line


    length = len(sequence)
    # section = length//80
    subSequence = list()
    # subSequence.append(head)
    for i in list(range(0,length-30)):
        if (i+100) >= length:
            sub = sequence[i:]
        else:
            sub = sequence[i:(i+100)]
        sub = sub.upper()
        subSequence.append(sub)

    # header = subSequence[0]
    header = header.lstrip('>')
    tmp = copy.deepcopy(header)
    _, id, date, region = tmp.split('|')
    id = id.split()[0]

    tmp1 = copy.deepcopy(date)
    tmp1 = tmp1.lower()
    tmp2 = copy.deepcopy(region)
    tmp2 = tmp2.lower()
    if filter == "yes":
        if ((tmp1=='na') or (tmp2=='na')):
            return -1, -1, -1
    if tmp2 in deletes:
        return -1, -1, -1
    if tmp2 in China:
        region = 'China'
        

    region = region.title()
    basicMessage = dict()
    basicMessage['Id'] = id
    basicMessage['Country'] = region
    basicMessage['Date'] = date


    splitFastaFileName = str(basicMessage['Id']) + '_split.fa'
    tempAnalysisDirectory = os.path.join(path, str(id))
    if os.path.exists(tempAnalysisDirectory):
        shutil.rmtree(tempAnalysisDirectory)
    os.mkdir(tempAnalysisDirectory)
    splitedFile = os.path.join(tempAnalysisDirectory, splitFastaFileName)
    with open(splitedFile, 'a') as fhand:
        for i, subsequence in enumerate(subSequence):
            readName = ">read" + str(i) + "\n" + str(subsequence) + "\n"
            fhand.write(readName)


    return basicMessage, splitedFile, tempAnalysisDirectory

def align(file, index, directory):
    '''
    alignment

    :param file: XXX_split.fa
    :param index: index directory
    :param directory: temp directory
    :return temp.sam
    '''
    samFile = os.path.join(directory, 'temp.sam')
    os.system('bowtie2 -f -x %s/index -U %s -S %s' % (index, file, samFile))


    return samFile
    
def sort(file, directory):
    '''
    sort reads

    :param file: temp.sam
    :param directory: temp directory
    :return temp.bam
    '''
    bamFile = os.path.join(directory, 'temp.bam')
    os.system('samtools sort %s > %s' % (file, bamFile))


    return bamFile

def mpileup(file, directory, refsequence):
    '''
    bam format -> vcf format

    :param file: temp.bam
    :param directory: temp directory
    :returns vcfFile: temp.vcf
    '''
    vcfFile = os.path.join(directory, 'temp.vcf')
    os.system('bcftools mpileup %s --fasta-ref %s > %s' % (file, refsequence, vcfFile))


    return vcfFile

def call(vcfFile, directory, ploidy, a=None):
    '''
    call SNPs

    :param vcfFile: temp.vcf
    :param directory: temp directory
    :param ploidy: ploidy.txt
    :param path: output directory
    :param a: number of INDELs
    :returns flag, snpFile
    '''
    snp_indel_file_path = os.path.join(directory, 'snp_indel.vcf')
    snp_file_path = os.path.join(directory, 'snp.vcf')
    indel_file_path = os.path.join(directory, 'indel.vcf')


    # can add [--threads <int>] to use multithreading
    os.system('bcftools call --ploidy 1 -vm %s -o %s' % (vcfFile, snp_indel_file_path))
    if a is None:
        os.system('vcftools --vcf %s --recode --remove-indels --stdout > %s' % (snp_indel_file_path, snp_file_path))
        return 0, snp_file_path
    else:
        os.system('vcftools --vcf %s --recode --keep-only-indels --stdout > %s' % (snp_indel_file_path, indel_file_path))
        n_indels = 0
        with open(indel_file_path, 'r') as fhand:
            for line in fhand.readlines():
                line = line.rstrip()
                if len(line)==0:
                    continue
                if line[0]=='#':
                    continue
                else :
                    n_indels = n_indels + 1
        if n_indels > a:
            return -1, -1
        else :
            os.system('vcftools --vcf %s --recode --remove-indels --stdout > %s' % (snp_indel_file_path, snp_file_path))
            return 0, snp_file_path

def snp_mutation_information(file):
    '''
    obtain SNPs information

    :param file: snp.vcf
    :returns mutationInformation(dict): key=['Position', 'Ref', 'Alt']
    '''

    mutation_lines = list()
    with open(file, 'r') as fhand:
        for line in fhand.readlines():
            line = line.rstrip()
            if len(line)==0:
                pass
            elif line[0]=="#":
                pass
            else:
                mutation_lines.append(line)
    
    mutationInformation = list()
    if len(mutation_lines) == 0:
        record = {'Position':0, 'Ref':'NA', 'Alt':'NA'}
        mutationInformation.append(record)
    else:
        for item in mutation_lines:
            item = item.split()
            pos = int(item[1])
            ref = str(item[3])
            alt = str(item[4])
            snpMutationMessage = dict()
            snpMutationMessage['Position'] = pos
            snpMutationMessage['Ref'] = ref
            snpMutationMessage['Alt'] = alt
            mutationInformation.append(snpMutationMessage)

    return mutationInformation

def module1(inputDirectory, outputDirectory, reference, collection_time="yes", length=None, number_n=None, number_db=None, number_indels=None):
    '''
    Call SNVs

    :param inputDirectory: raw genome sequences directory
    :param reference: reference genome sequence
    :param collection_time: whether filter out sequence that is without clear collection time or country or not
    :param length: length of genome sequence
    :param number_n: number of unknown bases
    :param number_db: number of degenerate bases
    :param number_indels: number of INDELs
    :param outputDirectory: output directory
    :return snp_merged.tsv, sequences_information.tsv
    '''
    path = os.path.abspath(outputDirectory)
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

    snpFile = os.path.join(path, 'snp_merged.tsv')
    seq_info = os.path.join(path, 'sequences_information.tsv')
    indexDirectory = os.path.join(path, 'index')
    os.mkdir(indexDirectory)
    refsequence = shutil.copy(reference, indexDirectory)
    genomeDirectory = os.path.join(path, 'genome')
    os.mkdir(genomeDirectory)


    if not os.path.exists(inputDirectory):
        print("%s doesn't exist."%inputDirectory)
        sys.exit()
    
    print('Extracting genome sequences...')
    files_path = extract_sequence(inputDirectory, genomeDirectory)
    print('Done!')
    total = len(files_path)
    print('There are %s sequences.'%total)

    ref_str = str()
    with open(refsequence, 'r') as fhand:
        for line in fhand.readlines():
            line = line.strip()
            if len(line) == 0:
                continue
            elif line[0] == ">":
                continue
            else:
                line = line.replace(" ", "")
                ref_str = ref_str + line
    length_ref = len(ref_str)
    print("Length of virus reference genome sequence is %s" % length_ref)


    create_index(indexDirectory, refsequence)

    columns = '''Id\tDate\tCountry\tPosition\tRef\tAlt\n'''
    with open(snpFile, 'a') as fhand:
        fhand.write(columns)
    

    pass_n = 0
    for file in files_path:
        flag = genome_quality_control(file, referenceLength=length_ref,a=length, b=number_n, c=number_db)
        if flag == -1:
            continue
        else:
            basicmessage, splitfasta, tempdirectory = split_sequence(file, path, filter=collection_time)
            if basicmessage == -1:
                continue
            else :
                samFile = align(splitfasta, indexDirectory, tempdirectory)
                bamFile = sort(samFile, tempdirectory)
                vcfFile = mpileup(bamFile, tempdirectory, refsequence)
                filter, vcfSnpFile = call(vcfFile, tempdirectory, ploidy, a=number_indels)
                if filter == -1:
                    if os.path.exists(tempdirectory):
                        shutil.rmtree(str(tempdirectory))
                    continue
                else:
                    pass_n = pass_n + 1


                snpMutationInformation = snp_mutation_information(vcfSnpFile)
                snp_mutation = dict()
                snp_mutation['Id'] = basicmessage['Id']
                snp_mutation['Date'] = basicmessage['Date']
                snp_mutation['Country'] = basicmessage['Country']


                with open(snpFile, 'a', encoding='utf-8') as fhand:
                    for snp in snpMutationInformation:
                        if (("Position" in snp) and ("Ref" in snp) and ("Alt" in snp)):
                            snp_mutation['Position'] = snp['Position']
                            snp_mutation['Ref'] = snp['Ref']
                            snp_mutation['Alt'] = snp['Alt']
                        else:
                            continue
                        record = str(snp_mutation['Id']) + "\t" + str(snp_mutation['Date']) + "\t" + str(snp_mutation['Country']) + "\t" + str(snp_mutation['Position']) + "\t" + str(snp_mutation['Ref']) + "\t" + str(snp_mutation['Alt']) +'\n'
                        fhand.write(record)
        if os.path.exists(logging):
            os.remove(logging)
        if os.path.exists(tempdirectory):
            shutil.rmtree(tempdirectory)
    if os.path.exists(indexDirectory):
        shutil.rmtree(indexDirectory)
    if os.path.exists(genomeDirectory):
        shutil.rmtree(genomeDirectory)


    not_pass_n = total - pass_n
    with open(seq_info, 'a') as fhand:
        line1 = 'Total genome sequences: %s'%total + '\n'
        line2 = 'Pass quality control: %s'%pass_n + '\n'
        line3 = 'Not pass quality control: %s'%not_pass_n + '\n'
        fhand.write(line1)
        fhand.write(line2)
        fhand.write(line3)
    print('Have successfully obtained all SNV mutations.')

    return snpFile