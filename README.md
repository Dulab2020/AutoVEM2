# README

### 0. Basic Message

Program: AutoVEM2 (works on Linux machine)

Description: Virus key mutations and epidemic trends analysis tool

Version: V2.1

Authors: Xibinbin

​				 Duhongli

Contacts: 201766841276@mail.scut.edu.cn 

​					hldu@scut.edu.cn

Year: 2021

License: Released under GNU General Public License

**Citation**

AutoVEM2 can be cited with the following paper:

Xi, B., Chen, Z., Li, S., Liu, W., Jiang, D., Bai, Y., ... & Du, H. (2021). AutoVEM2: A flexible automated tool to analyze candidate key mutations and epidemic trends for virus. *Computational and structural biotechnology journal*, *19*, 5029-5038.

### 1. Prerequisite

- Python(v3.8.6 or higher) with pandas, numpy and matplotlib packages (latest version)
- Java Runtime Environment (should meet the requirements of Haploview tool)

**Question**: How could I test whether the Java Runtime Environment meet the requirements of Haploview tool?

**Answer**: Run the following command and see if the help information could be outputted correctly:

```
cd ${PATH}/AutoVEM2
java -jar Haploview.jar -n -h
```

### 2. Dependencies

Bowtie2 (v2.4.2), SAMtools (v1.10), VCFtools (v0.1.16), Picard AddOrReplaceReadGroups (v2.27.4), and GATK4 HaplotypeCaller are required. Please make sure you have installed these tools and add them to the environment variables globally before using AutoVEM2. 

Haploview.jar (v4.2) has been provided. And please make sure the Java Runtime Environment meets the requirements of Haploview. For more details about Haploview, please visit: https://www.broadinstitute.org/haploview/tutorial

### 3. Installation

AutoVEM2 needn't to be installed and can be used directly. Once your have installed the dependencies, you can download the source code and run AutoVEM2 as follows:

```
git clone https://github.com/Dulab2020/AutoVEM2.git
cd ${PATH}/AutoVEM2
python run.py --help
python run.py call --help
python run.py analysis --help
python run.py plot --help
python run.py pipeline --help
```

### 4. Format of Genome Sequences

All genome sequences can be put into a fasta format file or several fasta format files. And put the file(s) into a directory. The directory will be as the input of AutoVEM2.

Format of virus genome sequence is as follows:

``` 
>virus_name|sequence_id|collection_date|country_or_region
Whole_genome_sequence
```

***Explanations in detail***

` virus_name`: Name of the virus, customizable. 

` sequence_id`: Identifier of the genome sequence, customizable. Make sure it is unique and do not have any spaces or tabs.

` collection_date`: Collection time of virus sequence. Format: ` year-month-day`, the following three examples are all valid ` 2020-12-13` , ` 2020-12` , ` 2020` . If collection time is not clear, use ` NA`.

` country_or_region`: Country or region where the sequence comes from. If there are several words, use spaces as separators. If not clear, use ` NA` .

***Examples***

```
>HPV|JN565301|2001-12-13|Croatia
ACTACAATAATTCATGTATAAAACTAAGGGCGTAACCGAAATCGGTTGA...ATGC
>HPV|JN565304|2001-12|Croatia
ACTACAATAATTCATGTATAAAACTAAGGGCGTAACCGAAATCGGTTGA...ATGC
>HPV|JN565305|2001|Croatia
ACTACAATAATTCATGTATAAAACTAAGGGCGTAACCGAAATCGGTTGA...ATGC
>HPV|JN565306|NA|Croatia
ACTACAATAATTCATGTATAAAACTAAGGGCGTAACCGAAATCGGTTGA...ATGC
>HPV|JN565307|NA|NA
ACTACAATAATTCATGTATAAAACTAAGGGCGTAACCGAAATCGGTTGA...ATGC
>HPV|JN565308|2001-12-13|NA
ACTACAATAATTCATGTATAAAACTAAGGGCGTAACCGAAATCGGTTGA...ATGC
>HPV|JN565309|2001-12-13|United States
ACTACAATAATTCATGTATAAAACTAAGGGCGTAACCGAAATCGGTTGA...ATGC
```

### 5. Usage

AutoVEM2 has three modules, including `call` module, `analysis` module, and `plot` module, which can be used modularly or as a whole. The `pipeline`  combines the three modules to perform a complete analysis function. But we recommend using AutoVEM2 modularly which can separate the whole analysis  into the following three steps:

- ` call`: **The first step**. This module will carry out quality control of the genomes and find all **SNV** mutations and stored them into the ` snp_merged.tsv` file which is the input file of the `analysis` module. And this module will also produce a ` sequences_information.tsv` file. It stores the summary result of quality control. 

  For more details, please use ` python run.py call -h` command line.

- ` analysis`: **The second step**. This module will obtain key mutation sites and find haplotype of every sequence. Haplotype information is stored in the ` data_plot.tsv` file. This file is the input file of the `plot` module.

  For more details, please use ` python run.py analysis -h` command line.

- ` plot`: **The third step.** Visualize the epidemic trends of the haplotypes.

  For more details, please use ` python run.py plot -h` command line.

### 6. Output Files

` snp_merged.tsv` 	Produced by the ` call` module, ` \t` delimited. It stores the information of all **SNV** mutations of all genome sequences. 

> **Id**: the unique identifier of a sequence 
>
> **Date**: collection date of the sequence 
>
> **Country**: country or region where the sequence was collected 
>
> **Position**: mutation position 
>
> **Ref**: reference base 
>
> **Alt**: mutation base

​           ![image-20210504091125418](https://github.com/Dulab2020/AutoVEM2/blob/main/images/image-20210504091125418.png)                               	

` sequences_information.tsv`	Produced by the ` call` module. It stores the summary information about quality control result.

![image-20210504091253330](https://github.com/Dulab2020/AutoVEM2/blob/main/images/image-20210504091253330.png)



` snp_sites.tsv`	Produced by the ` analysis` module , `\t` delimited. It stores the information of the specific mutation sites that meet your requirement.

> **Position:** mutation position
>
> **Ref:** reference base
>
> **Alt:** mutation base
>
> **Frequency:** mutation frequency

![image-20210504091534370](https://github.com/Dulab2020/AutoVEM2/blob/main/images/image-20210504091534370.png)



`plot.LD.PNG`	Produced by the ` analysis` module , linkage analysis results of the specific sites.

![plot.LD](https://github.com/Dulab2020/AutoVEM2/blob/main/images/plot.LD.PNG)



` haplotypes.tsv`	Produced by the `analysis` module, `\t` delimited. It stores the information about the haplotypes. Haplotypes whose frequency less than 1% are named as `other` .

> **Name**: the name of each haplotype.
>
> **Sequence**: the sequence of each haplotype.
>
> **Frequency**:  the frequency of each haplotype in the population.

![image-20210504091834774](https://github.com/Dulab2020/AutoVEM2/blob/main/images/image-20210504091834774.png)

​							

` data_plot.tsv`	Produced by the `analysis`  module, `\t` delimited. This file is the input file of the `plot` module. If you want to plot the epidemic trends of every haplotypes, when `call` SNVs, don't use 

`--region_date_filter no` command.

> **Id**: the identifier of the genome sequence
>
> **Date**: the date when the sequence was collected
>
> **Country**: the country or region where the sequence was collected
>
> **Hap**: the haplotype sequence of the genome sequence
>
> **Name**: the name of the haplotype

![image-20210504092628923](https://github.com/Dulab2020/AutoVEM2/blob/main/images/image-20210504092628923.png)

` hap_date.pdf`	Produced by the `plot` module . Show the epidemic trends of haplotypes.
![image-20210504092426244](https://github.com/Dulab2020/AutoVEM2/blob/main/images/image-20210504092426244.png)

## 7. Examples

An example has been provided, you can download it from the github repo to your local computer.

```shell
git clone https://github.com/Dulab2020/AutoVEM2.git
cd ../AutoVEM2
# call module
# This command can produced the snp_merged.tsv and sequences_information.tsv files
# All output by call command will be stored in the Example/call folder
python run.py call --input Example/genomes --ref Example/reference/ref_SARS-CoV-2.fa --length 29000 --number_n 15 --number_db 50 --number_indels 2 --output Example/call

# analysis module
# This command can produced the snp_sites.tsv, plot.LD.PNG, haplotypes.tsv, 
# hap_mutations.pdf and data_plot.tsv files
# All output by analysis command will be stored in the Example/analysis folder
python run.py analysis --input Example/call/snp_merged.tsv --ref Example/reference/ref_SARS-CoV-2.fa --frequency 0.2 --output Example/analysis

# plot module
# This command can produced the hap_date.pdf file
# All output by plot command will be stored in the Example/plot folder
python run.py plot --input Example/analysis/data_plot.tsv --days 7 --output Example/plot

# integrate the three into the pipeline module
# This command can produced all files mentioned above
# All output by pipeline command will be stored in the Example/pipeline folder
python run.py pipeline --input Example/genomes --ref Example/reference/ref_SARS-CoV-2.fa --length 29000 --number_n 15 --number_db 50 --number_indels 2 --frequency 0.2 --days 7 --output Example/pipeline
```
