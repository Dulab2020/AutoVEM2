#!/usr/bin/env python

#=================================================================================================================
# AutoVEM2
# Author: Xi Binbin
# Email: 201766841276@mail.scut.edu.cn
#=================================================================================================================
import os
import sys
import argparse
import call 
import analysis
import plot

def call_func(args):
    '''
    call SNVs
    '''
    inputDir = args.input
    outDir = args.output
    ref = args.ref
    length = args.length
    number_n = args.number_n
    number_db = args.number_db
    number_indels = args.number_indels
    collection_filter = args.region_date_filter
    if not os.path.exists(inputDir):
        print(r"Can not find directory %s."%inputDir)
        sys.exit()
    if not os.path.isdir(inputDir):
        print(r'Error: bad values of the "--input" argument, should be a directory rather than a file.')
        sys.exit()
    if not os.path.exists(ref):
        print(r"Can not find the %s file."%ref)
        sys.exit()
    if not os.path.isfile(ref):
        print(r"Error: %s is not a file."%ref)
        sys.exit()
    if length is None:
        pass
    else:
        if length<=0:
            print(r"Length should be bigger than 0.")
            sys.exit()
    if number_n is None:
        pass
    else:
        if number_n<=0:
            print(r'Value of the "--number_n" argument should be bigger than 0.')
            sys.exit()
    if number_db is None:
        pass
    else:
        if number_db<=0:
            print(r'Value of the "--number_db" argument should be bigger than 0.')
            sys.exit()
    if number_indels is None:
        pass
    else:
        if number_indels<=0:
            print(r'Value of the "--number_indels" argument should be bigger than 0.')
            sys.exit()
        
    call.module1(inputDir, outDir, ref, collection_time=collection_filter, length=length, number_n=number_n, number_db=number_db, number_indels=number_indels)

def analysis_func(args):
    '''
    analysis
    '''
    inputFile = args.input
    outDir = args.output
    ref = args.ref
    frequency = args.frequency
    sites = args.sites
    if sites is None:
        pass
    else:
        sites.sort()
    if not os.path.exists(inputFile):
        print("Can not find the %s file."%inputFile)
        sys.exit()
    if not os.path.isfile(inputFile):
        print("Error: %s is not a file."%inputFile)
        sys.exit()
    if not os.path.exists(ref):
            print("Can not find the %s file."%ref)
            sys.exit()
    if not os.path.isfile(ref):
        print("Error: %s is not a file."%ref)
        sys.exit()
    if frequency is None:
        pass
    else:
        if ((frequency<=0) or (frequency>1)):
            print('Value of the "--frequency" argument should not less than 0 and bigger than 1.')
            sys.exit()
    if sites is None:
        pass
    else:
        minimum = sites[0]
        if minimum<=0:
            print("There is(are) negative value(s) in your giver sites.")
            sys.exit()

    path = os.path.abspath(outDir)
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

    analysis.module2(inputFile, outDir, ref, sites=sites, frequency=frequency)

def plot_func(args):
    '''
    plot 
    '''
    dataPlot = args.input
    outDir = args.output
    days = args.days
    if not os.path.exists(dataPlot):
        print("Can not find file %s."%dataPlot)
        sys.exit()
    if not os.path.isfile(dataPlot):
        print("Error: %s is not a file."%dataPlot)
        sys.exit()
    if days<=0:
        print("Error: bad values of '--days' argument, should be a positive integer.")
        sys.exit()
    
    path = os.path.abspath(outDir)
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
    
    plot.module3(dataPlot, outDir, days)

def pipeline_func(args):
    '''
    pipeline analysis
    '''
    inputDir = args.input
    outDir = args.output
    ref = args.ref
    days = args.days
    length = args.length
    number_n = args.number_n
    number_db = args.number_db
    number_indels = args.number_indels
    frequency = args.frequency
    collection_filter = args.region_date_filter
    sites = args.sites
    if sites is None:
        pass
    else:
        sites.sort()
    if not os.path.exists(inputDir):
        print("Can not find directory %s."%inputDir)
        sys.exit()
    if not os.path.isdir(inputDir):
        print('Error: bad values of the "--input" argument, should be a directory rather than a file.')
        sys.exit()
    if not os.path.exists(ref):
        print("Can not find the %s file."%ref)
        sys.exit()
    if not os.path.isfile(ref):
        print("Error: %s is not a file."%ref)
        sys.exit()
    if length is None:
        pass
    else:
        if length<0:
            print("Length should be bigger than or equal to 0.")
            sys.exit()
    if number_n is None:
        pass
    else:
        if number_n<0:
            print('Value of the "--number_n" argument should be bigger than or equal to 0.')
            sys.exit()
    if number_db is None:
        pass
    else:
        if number_db<0:
            print('Value of the "--number_db" argument should be bigger than or equal to 0.')
            sys.exit()
    if number_indels is None:
        pass
    else:
        if number_indels<0:
            print('Value of the "--number_indels" argument should be bigger than or equal to 0.')
            sys.exit()
    if frequency is None:
        pass
    else:
        if ((frequency<=0) or (frequency>1)):
            print('Value of the "--frequency" argument should not less than 0 and bigger than 1.')
            sys.exit()
    if sites is None:
        pass
    else:
        minimum = sites[0]
        if minimum<=0:
            print("There is(are) negative value(s) in your given sites.")
            sys.exit()
    if days<=0:
        print("Error: bad values of '--days' argument, should be a positive integer.")
        sys.exit()
    snpFile = call.module1(inputDir, outDir, ref, collection_time=collection_filter, length=length, number_n=number_n, number_db=number_db, number_indels=number_indels)
    dataPlot = analysis.module2(snpFile, outDir, ref, sites=sites, frequency=frequency)
    plot.module3(dataPlot, outDir, days)

parser = argparse.ArgumentParser(prog='AutoVEM')
parser.add_argument('--version', '-v', action='version', version='AutoVEM2 V2.1')
subparsers = parser.add_subparsers()

call_parser = subparsers.add_parser('call', help='Quality control of genomes and call all SNVs')
call_parser.add_argument('--input', type=str, required=True, help='Directory stores fasta format virus genome sequeces file(s).')
call_parser.add_argument('--ref', type=str, required=True, help='Fasta format reference genome file.')
call_parser.add_argument('--length', type=float, default=None, help="Sequence whose length is less than the given length(int) or percent(0-1) of reference length will be filtered out.")
call_parser.add_argument('--number_n', type=float, default=None, help="Sequence with unknown base exceeding the given value(int) or percent(0-1) of reference length will be filtered out.")
call_parser.add_argument('--number_db', type=float, default=None, help="Sequence with degenerate base exceeding the given value(int) or percent(0-1) of reference length will be filtered out.")
call_parser.add_argument('--number_indels', type=int, default=None, help="Sequence with INDEL counts exceeding the given value will be filtered out.")
call_parser.add_argument('--region_date_filter', choices=['yes','no'], default="yes", help="if yes, sequence with unclear collection time or unclear region will be filtered out, default yes.")
call_parser.add_argument('--output', type=str, required=True, help='Output Directory.')
call_parser.set_defaults(func=call_func)

analysis_parser = subparsers.add_parser('analysis', help='Screen out candidate key mutations and acquire haplotye of each genomes.')
analysis_parser.add_argument('--input', type=str, required=True, help="snp_merged.tsv file produced by the call module.")
analysis_parser.add_argument('--ref', type=str, required=True, help='Fasta format reference genome file.')
analysis_parser.add_argument('--frequency', type=float, default=None, help='''Find SNVs with mutation frequency equal to or more than the given frequency (0,1]''')
analysis_parser.add_argument('--sites', type=int, nargs='+', default=None, help='''Space delimited list of integer numbers.''')
analysis_parser.add_argument('--output', type=str, required=True, help='Output directory')
analysis_parser.set_defaults(func=analysis_func)

plot_parser = subparsers.add_parser('plot', help='Visualize virus epidemic trends.')
plot_parser.add_argument('--input', type=str, required=True, help="data_plot.tsv file.")
plot_parser.add_argument('--days', type=int, required=True, help='Day interval.')
plot_parser.add_argument('--output', type=str, required=True, help='Output directory')
plot_parser.set_defaults(func=plot_func)

pipeline_parser = subparsers.add_parser('pipeline', help="Run the analysis completely")
pipeline_parser.add_argument('--input', type=str, required=True, help='Directory stores virus genome sequeces.')
pipeline_parser.add_argument('--ref', type=str, required=True, help='Fasta format reference genome file.')
pipeline_parser.add_argument('--frequency', type=float, default=None, help='''Find SNVs with mutation frequency equal to or more than the given frequency (0,1]''')
pipeline_parser.add_argument('--sites', type=int, nargs='+', default=None, help='''Space delimited list of integer numbers.''')
pipeline_parser.add_argument('--number_n', type=float, default=None, help="Sequence with unknown base exceeding the given value(int) or percent(0-1) of reference length will be filtered out.")
pipeline_parser.add_argument('--number_db', type=float, default=None, help="Sequence with degenerate base exceeding the given value(int) or percent(0-1) of reference length will be filtered out.")
pipeline_parser.add_argument('--length', type=float, default=None, help="Sequence whose length is less than the given length(int) or percent(0-1) of reference length will be filtered out.")
pipeline_parser.add_argument('--number_indels', type=int, default=None, help="Sequence with INDEL counts exceeding the given value will be filtered out.")
pipeline_parser.add_argument('--region_date_filter', choices=['yes','no'], default=None, help="if yes, sequence with unclear collection time or unclear region will be filtered out, default yes.")
pipeline_parser.add_argument('--days', type=int, required=True, help='Day interval.')
pipeline_parser.add_argument('--output', type=str, required=True, help='Output directory')
pipeline_parser.set_defaults(func=pipeline_func)

args = parser.parse_args()
args.func(args)
