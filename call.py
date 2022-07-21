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
    extract each genome sequence from the multi-sequences files

    :param directory: the input directory
    :param genomeDirectory: the absolute path of output genome files directory

    :return filesPath
    '''
    ## get the absolute path of the multi-sequences fasta files
    in_dir = os.path.abspath(directory)
    files = os.listdir(in_dir)
    files_list = list()
    for file in files:
        temp = os.path.join(in_dir, file)
        if os.path.isfile(temp):
            files_list.append(temp)
        else:
            pass
    ## pointer: indicates whether the sequence passed the quality control
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
                    idx = line[1]
                    idx = "_".join(idx.split())
                    idx = idx + '.fa'
                    id_path = os.path.join(genomeDirectory, idx)
                    if os.path.exists(id_path):
                        pointer = 1
                        continue
                    else:
                        with open(id_path, 'a') as f:
                            f.write(tmp+'\n')
                            pointer = 0
                            continue
                if pointer == 0:
                    line = line.replace(" ", "")
                    line = line.replace("\t", "")
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

def create_index(refsequence):
    '''
    create index

    :param refsequence: the fasta fomat reference sequence file
    '''
    ## bowtie2 index and .fai index
    os.system(f'bowtie2-build -f {refsequence} {refsequence}')
    os.system(f"samtools faidx {refsequence}")
    ## dict index
    prefix, suffix = os.path.split(refsequence)
    fileName = ".".join(suffix.split(".")[:-1]) + ".dict"
    outFile = os.path.join(prefix, fileName)
    os.system(f"samtools dict {refsequence} -o {outFile}")

def genome_quality_control(file, referenceLength=None, a=None, b=None, c=None):
    '''
    quality control for the first time

    :param file: a fasta format file of a genome sequence
    :param a: length
    :param b: number of unknown bases, number_n
    :param c: number of degenerate bases, number_db

    :returns flag
    '''
    ## test whether the genome sequence is readable
    try:
        os.path.getsize(file)
    except:
        return -1
    ## get the nucleotide sequence of the genome
    with open(file, 'r') as fhand:
        genomeSequence = str()
        for line in fhand.readlines():
            line = line.strip()
            if len(line)==0:
                continue
            if line[0]=='>':
                continue
            genomeSequence = genomeSequence + line
    len_sequence = len(genomeSequence)

    ## count the number of the unknown bases and the degenerate bases
    number_n = 0
    number_db = 0
    letters = {'A', 'T', 'G', 'C', 'N'}
    for item in genomeSequence:
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
    ## get the header line of the genome sequence
    ## get the nucleotide sequence of the genome sequence
    with open(file, 'r') as fhand:
        header = str()
        sequence = str()
        for line in fhand.readlines():
            line = line.rstrip()
            if len(line)==0:
                continue
            if line[0]==">":
                header = line[:]
            else:
                sequence = sequence + line

    ## length of the genome sequence
    length = len(sequence)
    subSequence = list()

    ## split the genome sequence into 100bp reads
    ## the splice window is 100bp
    ## the step is 1bp at each time
    for i in list(range(0,length-30)):
        if (i+100) >= length:
            sub = sequence[i:]
        else:
            sub = sequence[i:(i+100)]
        subSequence.append(sub)
    
    ## get the metadata of the genome sequence
    ## and test whether the region or collection date of the sequence is unknown
    header = header[1:]
    _, idx, date, region = header.split('|')
    if filter == "yes":
        if ((date.lower()=='na') or (region.lower()=='na')):
            return -1, -1, -1  
    ## get the basic information of the genome sequence
    basicMessage = dict()
    idx = "_".join(idx.split())
    basicMessage['Id'] = idx
    basicMessage['Country'] = region
    basicMessage['Date'] = date

    ## write the spliced reads to the fasta file
    splitFastaFileName = str(basicMessage['Id']) + '_split.fa'
    tempAnalysisDirectory = os.path.join(path, str(idx))
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
    :param index: the reference genome sequence
    :param directory: temp directory
    :return temp.sam
    '''
    samFile = os.path.join(directory, 'temp.sam')
    os.system(f'bowtie2 -f -x {index} -U {file} -S {samFile}')

    return samFile
    
def sort(file, directory):
    '''
    sort reads

    :param file: temp.sam
    :param directory: temp directory
    :return temp.bam
    '''
    bamFile = os.path.join(directory, 'temp.bam')
    os.system(f'samtools sort {file} > {bamFile}')
    os.system(f"samtools index {bamFile}")

    return bamFile

def mpileup(file, directory):
    '''
    add @RG to the bam file

    :param file: temp.bam
    :param directory: temp directory

    :returns vcfFile: temp.vcf
    '''
    
    addHeader = os.path.join(directory, "temp_addheader.bam")
    os.system(f"picard AddOrReplaceReadGroups -I {file} -O {addHeader} --RGID Sample --RGLB AMPLICON --RGPL ILLUMINA --RGPU unit1 --RGSM Sample")
    os.system(f"samtools index {addHeader}")

    return addHeader

def call(vcfFile, directory, refSeq, a=None):
    '''
    call SNPs

    :param vcfFile: temp_addheader.bam
    :param directory: temp directory
    :param path: output directory
    :param a: number of INDELs

    :returns flag, snpFile
    '''
    ## the calling variant files
    snp_indel_file_path = os.path.join(directory, 'snp_indel.vcf')
    snp_file_path = os.path.join(directory, 'snp.vcf')
    indel_file_path = os.path.join(directory, 'indel.vcf')

    ## call variants
    os.system(f'gatk HaplotypeCaller -I {vcfFile} -O {snp_indel_file_path} -R {refSeq} -ploidy 1')
    ## filter the raw variants according the number of INDELs
    if a is None:
        os.system(f'vcftools --vcf {snp_indel_file_path} --recode --remove-indels --stdout > {snp_file_path}')
        return 0, snp_file_path
    else:
        os.system(f'vcftools --vcf {snp_indel_file_path} --recode --keep-only-indels --stdout > {indel_file_path}')
        n_indels = 0
        with open(indel_file_path, 'r') as fhand:
            for line in fhand.readlines():
                line = line.rstrip()
                if len(line)==0:
                    pass
                elif line[0]=='#':
                    pass
                else:
                    n_indels = n_indels + 1
        if n_indels > a:
            return -1, -1
        else:
            os.system(f'vcftools --vcf {snp_indel_file_path} --recode --remove-indels --stdout > {snp_file_path}')
            return 0, snp_file_path

def snp_mutation_information(file):
    '''
    obtain SNPs information

    :param file: snp.vcf
    :returns mutationInformation(dict): key=['Position', 'Ref', 'Alt']
    '''
    ## obtain the all mutation records of a genome sequence
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
    ## 
    mutationInformation = list()
    ## if there isn't SNP mutations of a sequence, the mutation position will be designated as 0
    ## and the reference base and alternate base will be designated as NA
    if len(mutation_lines) == 0:
        record = {'Position':0, 'Ref':'NA', 'Alt':'NA'}
        mutationInformation.append(record)
    else:
        ## handle each mutation records
        for item in mutation_lines:
            snpMutationMessage = dict()
            item = item.split()
            pos = int(item[1])
            ref = str(item[3])
            alt = str(item[4])
            snpMutationMessage['Position'] = pos
            snpMutationMessage['Ref'] = ref
            snpMutationMessage['Alt'] = alt
            mutationInformation.append(snpMutationMessage)

    return mutationInformation

def module1(inputDirectory, outputDirectory, reference, collection_time="yes", length=None, number_n=None, number_db=None, number_indels=None):
    '''
    Call SNVs

    :param inputDirectory: raw genome sequences directory
    :param outputDirectory: output directory
    :param reference: reference genome sequence
    :param collection_time: whether filter out sequence that is without clear collection time or country or not
    :param length: length of genome sequence
    :param number_n: number of unknown bases
    :param number_db: number of degenerate bases
    :param number_indels: number of INDELs
    :return snp_merged.tsv, sequences_information.tsv
    '''
    ## if the input directory doesn't exists, exit
    if not os.path.exists(inputDirectory):
        print(f"{inputDirectory} doesn't exist.")
        sys.exit()

    ## path: the absolute path of the output directory
    ## decide whether the output directory has exists,
    ## if exists, backup the already exists directory,
    ## and make a new empty output directory.
    path = os.path.abspath(outputDirectory)
    if os.path.exists(path):
        print(f"{path} already exists.")
        dirname, filename = os.path.split(path)
        i = 0
        while True:
            i = i + 1
            temp = os.path.join(dirname, r"#%s.%s#"%(filename,i))
            if os.path.exists(temp):
                continue
            else :
                os.rename(path, temp)
                print(f"Back up {path} to {temp}")
                break
    os.mkdir(path)

    ## make files and directories
    snpFile = os.path.join(path, 'snp_merged.tsv')
    columns = '''Id\tDate\tCountry\tPosition\tRef\tAlt\n'''
    with open(snpFile, 'a') as fhand:
        fhand.write(columns)
    seq_info = os.path.join(path, 'sequences_information.tsv')

    ## make the index directory and copy the reference file to the index directory
    indexDirectory = os.path.join(path, 'index')
    os.mkdir(indexDirectory)
    ref_path = os.path.abspath(reference)
    refsequence = shutil.copy(ref_path, indexDirectory)

    ## make the genome directory in the output directory
    ## this directory will store the extracted genome sequences
    genomeDirectory = os.path.join(path, 'genome')
    os.mkdir(genomeDirectory)
    
    ## extract genome sequences and put each sequence into a single file
    ## the extracted genome sequences will be stored in the above genome directory
    ## files_path: the list of absolute file path of the undercalling files
    print('Extracting genome sequences...')
    files_path = extract_sequence(inputDirectory, genomeDirectory)
    print('Done!')
    total = len(files_path)
    print(f'There are {total} sequences.')

    ## get the reference sequence
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
                line = line.replace("\t", "")
                line = line.upper()
                ref_str = ref_str + line
    length_ref = len(ref_str)
    print(f"Length of virus reference genome sequence is {length_ref}")

    ## create the index of the reference genome sequence
    create_index(refsequence)
    ## call variants
    pass_n = 0
    for file in files_path:
        flag = genome_quality_control(file, referenceLength=length_ref,a=length, b=number_n, c=number_db)
        if flag == -1:
            continue
        else:
            basicmessage, splitfasta, tempdirectory = split_sequence(file, path, filter=collection_time)
            if basicmessage == -1:
                continue
            else:
                samFile = align(splitfasta, refsequence, tempdirectory)
                bamFile = sort(samFile, tempdirectory)
                vcfFile = mpileup(bamFile, tempdirectory)
                filter, vcfSnpFile = call(vcfFile, tempdirectory, refsequence, a=number_indels)
                if filter == -1:
                    if os.path.exists(tempdirectory):
                        shutil.rmtree(str(tempdirectory))
                    continue
                else:
                    pass_n = pass_n + 1
                ## write the SNP variant information to the mutation records file
                snpMutationInformation = snp_mutation_information(vcfSnpFile)
                Idx = basicmessage['Id']
                Date = basicmessage['Date']
                Country = basicmessage['Country']
                with open(snpFile, 'a', encoding='utf-8') as fhand:
                    for snp in snpMutationInformation:
                        Position = snp['Position']
                        Ref = snp['Ref']
                        Alt = snp['Alt']
                        record = str(Idx) + "\t" + str(Date) + "\t" + str(Country) + "\t" + str(Position) + "\t" + str(Ref) + "\t" + str(Alt) +'\n'
                        fhand.write(record)
        ## remove the temporary file and directory
        if os.path.exists(logging):
            os.remove(logging)
        if os.path.exists(tempdirectory):
            shutil.rmtree(tempdirectory)
    ## remove the index and genome directory
    if os.path.exists(indexDirectory):
        shutil.rmtree(indexDirectory)
    if os.path.exists(genomeDirectory):
        shutil.rmtree(genomeDirectory)
    ## write the quality control message to file
    not_pass_n = total - pass_n
    with open(seq_info, 'a') as fhand:
        line1 = f'Total genome sequences: {total}\n'
        line2 = f'Pass quality control: {pass_n}\n'
        line3 = f'Not pass quality control: {not_pass_n}\n'
        fhand.write(line1)
        fhand.write(line2)
        fhand.write(line3)
    print('Have successfully obtained all SNV mutations.')

    return snpFile