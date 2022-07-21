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

def arguments():
    parser = argparse.ArgumentParser(prog='AutoVEM2', description="Variant calling of virus whole genome sequence")
    parser.add_argument('--version', '-v', action='version', version='AutoVEM2 v1.2')
    subparsers = parser.add_subparsers()
    ## call module
    call_parser = subparsers.add_parser('call', help='Quality control of genomes and call SNVs')
    call_parser.add_argument('--input', type=str, required=True, help='Directory that stores the fasta format virus genome sequeces file(s)')
    call_parser.add_argument('--ref', type=str, required=True, help='File of fasta format reference genome sequence')
    call_parser.add_argument('--length', type=float, default=None, help="Sequence whose length is less than the given length(int) or percent(0-1) of reference length will be filtered out")
    call_parser.add_argument('--number_n', type=float, default=None, help="Sequence with the number of unknown bases exceeding the given value(int) or the percent(0-1) of reference length will be filtered out")
    call_parser.add_argument('--number_db', type=float, default=None, help="Sequence with the number of degenerate bases exceeding the given value(int) or the percent(0-1) of reference length will be filtered out")
    call_parser.add_argument('--number_indels', type=int, default=None, help="Sequence with the number of INDELs exceeding the given value will be filtered out")
    call_parser.add_argument('--region_date_filter', choices=['yes','no'], default="yes", help="if yes, sequence with unclear collection time or unclear region will be filtered out, default yes")
    call_parser.add_argument('--output', type=str, required=True, help='Output Directory')
    call_parser.set_defaults(func=call_func)
    ## analysis module
    analysis_parser = subparsers.add_parser('analysis', help='Screen out candidate key mutations and acquire haplotye of each genomes')
    analysis_parser.add_argument('--input', type=str, required=True, help="the snp_merged.tsv file produced by the call module")
    analysis_parser.add_argument('--ref', type=str, required=True, help='Fasta format reference genome file')
    analysis_parser.add_argument('--frequency', type=float, default=None, help='''SNVs with mutation frequency less than the given frequency will be filtered out''')
    analysis_parser.add_argument('--sites', type=int, nargs='+', default=None, help='''Space delimited list of positions that of interest''')
    analysis_parser.add_argument('--output', type=str, required=True, help='Output directory')
    analysis_parser.set_defaults(func=analysis_func)
    ## plot module
    plot_parser = subparsers.add_parser('plot', help='Visualize virus epidemic trends')
    plot_parser.add_argument('--input', type=str, required=True, help="the data_plot.tsv file")
    plot_parser.add_argument('--days', type=int, default=7, help='time interval, default 7 days')
    plot_parser.add_argument('--output', type=str, required=True, help='Output directory')
    plot_parser.set_defaults(func=plot_func)
    ## pipeline
    pipeline_parser = subparsers.add_parser('pipeline', help="Run the analysis pipeline completely")
    pipeline_parser.add_argument('--input', type=str, required=True, help='Directory which stores virus genome sequeces')
    pipeline_parser.add_argument('--ref', type=str, required=True, help='Fasta format reference genome file')
    pipeline_parser.add_argument('--frequency', type=float, default=None, help='''SNVs with mutation frequency less than the given frequency will be filtered out''')
    pipeline_parser.add_argument('--sites', type=int, nargs='+', default=None, help='''Space delimited list of positions that of interest''')
    pipeline_parser.add_argument('--number_n', type=float, default=None, help="Sequence with the number of unknown bases exceeding the given value(int) or the percent(0-1) of reference length will be filtered out")
    pipeline_parser.add_argument('--number_db', type=float, default=None, help="Sequence with the number of degenerate bases exceeding the given value(int) or the percent(0-1) of reference length will be filtered out")
    pipeline_parser.add_argument('--length', type=float, default=None, help="Sequence whose length is less than the given length(int) or the percent(0-1) of reference length will be filtered out")
    pipeline_parser.add_argument('--number_indels', type=int, default=None, help="Sequence with the number of INDELs exceeding the given value will be filtered out")
    pipeline_parser.add_argument('--region_date_filter', choices=['yes','no'], default=None, help="if yes, sequence with unclear collection time or unclear region will be filtered out, default yes")
    pipeline_parser.add_argument('--days', type=int, default=7, help='time interval, default 7 days')
    pipeline_parser.add_argument('--output', type=str, required=True, help='Output directory')
    pipeline_parser.set_defaults(func=pipeline_func)

    args = parser.parse_args()

    return args

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
    ## get the absolute path of the input directory
    ## check whether the input directory exists
    inputDir = os.path.abspath(inputDir)
    if not os.path.exists(inputDir):
        print(f"Can not find directory {inputDir}.")
        sys.exit()
    if not os.path.isdir(inputDir):
        print('Error: bad values of the --input argument, should be a directory rather than a file.')
        sys.exit()
    ## get the absolute path of the output directory
    outDir = os.path.abspath(outDir)
    ## get the absolute path of the reference genome file
    ## check whether the reference genome file exists
    ref = os.path.abspath(ref)
    if not os.path.exists(ref):
        print(f"Can not find the {ref} file.")
        sys.exit()
    if not os.path.isfile(ref):
        print(f"Error: {ref} is not a file.")
        sys.exit()

    if length is None:
        pass
    else:
        if length<=0:
            print("Length should be bigger than 0.")
            sys.exit()
    if number_n is None:
        pass
    else:
        if number_n<=0:
            print('Value of the --number_n argument should be bigger than 0.')
            sys.exit()
    if number_db is None:
        pass
    else:
        if number_db<=0:
            print('Value of the --number_db argument should be bigger than 0.')
            sys.exit()
    if number_indels is None:
        pass
    else:
        if number_indels<=0:
            print('Value of the --number_indels argument should be bigger than 0.')
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
    ## get the absolute path of the input file
    ## check if the input file exists
    inputFile = os.path.abspath(inputFile)
    if not os.path.exists(inputFile):
        print(f"Can not find the {inputFile} file.")
        sys.exit()
    if not os.path.isfile(inputFile):
        print(f"Error: {inputFile} is not a file.")
        sys.exit()
    ## chech the output argument
    path = os.path.abspath(outDir) 
    ## get the absolute path of the reference genome sequence
    ## check whether the reference genome file exists
    ref = os.path.abspath(ref)
    if not os.path.exists(ref):
        print(f"Can not find the {ref} file.")
        sys.exit()
    if not os.path.isfile(ref):
        print(f"Error: {ref} is not a file.")
        sys.exit()
    ## check whether the value of the --frequency is valid
    if frequency is None:
        pass
    else:
        if ((frequency<=0) or (frequency>1)):
            print('Value of the --frequency argument should not less than 0 and bigger than 1.')
            sys.exit()
    ## the validation of the --sites argument is checked in the module2() function

    analysis.module2(inputFile, path, ref, sites=sites, frequency=frequency)

def plot_func(args):
    '''
    plot 
    '''
    dataPlot = args.input
    outDir = args.output
    days = args.days
    if days is None:
        days = 7
    if not os.path.exists(dataPlot):
        print("Can not find file %s."%dataPlot)
        sys.exit()
    if not os.path.isfile(dataPlot):
        print("Error: %s is not a file."%dataPlot)
        sys.exit()
    if days<=0:
        print("Error: bad values of '--days' argument, should be a positive integer.")
        sys.exit()
    
    outDir = os.path.abspath(outDir)

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
    if days is None:
        days = 7
    snpFile = call.module1(inputDir, outDir, ref, collection_time=collection_filter, length=length, number_n=number_n, number_db=number_db, number_indels=number_indels)
    dataPlot = analysis.module2(snpFile, outDir, ref, sites=sites, frequency=frequency)
    plot.module3(dataPlot, outDir, days)

if __name__ == "__main__":
    args = arguments()
    args.func(args)
