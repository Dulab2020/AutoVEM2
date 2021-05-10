#=================================================================================================================
# AutoVEM2
# Author: Xi Binbin
# Email: 201766841276@mail.scut.edu.cn
#=================================================================================================================

import os
import sys
import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt 
import datetime
from math import ceil
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.patches as mpatches
import matplotlib.colors as mcolors


def plot(file, directory, days):
    '''
    结果可视化

    :param file: data_plot.tsv
    :param directory: ouput directory
    :param days: interval days
    '''

    pdfFile = os.path.join(directory, 'hap_date.pdf')
    df = pd.read_csv(file, sep='\t')
    df.index = pd.to_datetime(df["Date"])
    df = df.sort_index(ascending=True)

    n_genomes = df.shape[0]
    haplotypes = df["Name"].unique().tolist()
    countries = df["Country"].unique().tolist()
    n_haplotypes=len(haplotypes)
    n_countries = len(countries)

    n=days
    if n_genomes <= 1:
        print("Too little genomes.")
        sys.exit()
    start_day = df.index[0]
    end_day = df.index[-1]
    delta = end_day - start_day
    m = ceil(delta.days/n)
    n_periods = list(range(m))
    date_point = list()
    date_point.append(start_day.strftime('%Y/%m/%d'))
    for c in n_periods:
        b = pd.to_datetime(start_day+datetime.timedelta(days = int((c+1)*n)))
        date_point.append(b.strftime('%Y/%m/%d'))

    def name(n):
            haps = list()
            for i in list(range(n-1)):
                s = 'H' + str(i+1)
                haps.append(s)
            haps.append('other')
            return haps

    def size(n):
        m=0
        if(n>=1000):
            m=1
        elif(n>=100):
            m=0.85
        elif(n>=10):
            m=0.7
        else:
            m=0.5
        return m

    sort_hap = name(n_haplotypes)
    colors_all=list(mcolors.CSS4_COLORS.keys())
    colors=colors_all[10:(n_haplotypes+10)]
    colors_dict = dict(zip(sort_hap, colors))


    x_size = 2.5  + 0.396 + m*0.396 + 0.396 + 0.1 + 0.396*7
    y_size = 1.2 + 0.396 + n_countries*0.396 + 0.396 + 0.396*2
    a = x_origin = 2.5/x_size
    b = y_origin = 1.2/y_size
    c = x_len = (m+2)*0.396/x_size
    d = y_len = (n_countries+2)*0.396/y_size
    lx = 0.396/x_size
    ly = 0.396/y_size

    # pdf = PdfPages(pdfFile)
    plt.figure(figsize=(x_size, y_size))
    axes_main = plt.axes([a,b,c,d])
    axes_main.set_xlim(0, m+2)
    axes_main.set_ylim(0, n_countries+2)
    xticks = list(range(1, m+2))
    yticks = list(np.arange(1.5, n_countries+1, 1))
    xticklabels = date_point[:]
    yticklabels = countries[::-1]
    axes_main.set_xticks(xticks)
    axes_main.set_yticks(yticks)
    axes_main.set_xticklabels(xticklabels)
    axes_main.xaxis.set_tick_params(rotation=30, labelsize=12)
    axes_main.set_yticklabels(yticklabels, size=12)


    for c in n_periods:
        a = pd.to_datetime(start_day+datetime.timedelta(days = int(c*n)))
        b = pd.to_datetime(start_day+datetime.timedelta(days = int((c+1)*n)))
        if b == end_day:
            df2 = df[(df.index >=a) & (df.index <= b)]
        else :
            df2 = df[(df.index >=a) & (df.index < b)]
        for r, j in enumerate(yticklabels):
            df3 = df2[df2["Country"] == j]
            haps_count = list()
            haps_color = list()
            for i, k in enumerate(sort_hap):
                df4 = df3[df3["Name"] == k]
                hap_num = len(df4.index)
                haps_count.append(hap_num)
                haps_color.append(colors_dict[k])
            num_all = sum(haps_count)
            if num_all == 0:
                continue
            p1 = x_origin + lx + lx*c
            p2 = y_origin + ly + ly*r
            p3 = lx
            p4 = ly
            ax = plt.axes([p1, p2, p3, p4])
            ax.axis("off")
            ax.set_xlim(0, 1)
            ax.set_ylim(0, 1)
            radius = size(num_all)
            ax.pie(haps_count, colors=haps_color, radius=radius, normalize=True)

    handle = list()      
    for hap in sort_hap:
        hap_patch = mpatches.Patch(color=colors_dict[hap], label=str(hap))
        handle.append(hap_patch)
    p1 = (x_size-7*0.396)/x_size
    p2 = 1.2/y_size
    p3 = 4*0.396/x_size
    p4 = (y_size-1.2-0.396*2)/y_size
    axes_main = plt.axes([p1,p2,p3,p4])
    axes_main.axis("off")
    axes_main.set_xticks([])
    axes_main.set_yticks([])
    axes_main.legend(handles=handle, loc="upper center", ncol=1, frameon=False, fontsize=12)

    p2 = (x_size-4*0.396)/x_size
    en = (y_size-0.396*3)/y_size
    ax4=plt.axes([p2, en, lx, ly])
    ax4.pie([1],labels=[""],colors=['white'],rotatelabels=True, wedgeprops={'linewidth': 1, 'edgecolor': "black"}, radius=1.0)
    ax4.annotate("n>=1000",xy=(p1, en), fontsize=12)
    ax5=plt.axes([p2, en-ly, lx, ly])
    ax5.pie([1],labels=[""],colors=['white'],rotatelabels=True, wedgeprops={'linewidth': 1, 'edgecolor': "black"}, radius=0.85)
    ax5.annotate("n>=100",xy=(p1, en-ly), fontsize=12)
    ax6=plt.axes([p2, en-2*ly, lx, ly])
    ax6.pie([1],labels=[""],colors=['white'],rotatelabels=True, wedgeprops={'linewidth': 1, 'edgecolor': "black"}, radius=0.7)
    ax6.annotate("n>=10",xy=(p1, en-2*ly), fontsize=12)
    ax7=plt.axes([p2, en-3*ly, lx, ly])
    ax7.pie([1],labels=[""],colors=['white'],rotatelabels=True, wedgeprops={'linewidth': 1, 'edgecolor': "black"}, radius=0.5)
    ax7.annotate("n>=1",xy=(p1, en-3*ly), fontsize=12)
    plt.savefig(pdfFile)
    plt.close()

    
def module3(dataPlot, outDirectory, days):
    file = dataPlot
    if not os.path.exists(os.path.abspath(file)):
        print("Error: can't find the %s file"%(os.path.abspath(file)))
        sys.exit()

    print('Plotting...')
    plot(dataPlot, outDirectory, days)
    print('Done!')
