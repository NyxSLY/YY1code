"""
make line plots of insulation score
wholegenome(), make line plots for each chromosome, input folder path
single_chr(), make line plot for the chromosome, can plot region if given coords
boxplot(), no statistics yet
"""

from mirnylib import h5dict, genome
from mirnylib import plotting
import sys
import numpy as np
import matplotlib.pyplot as plt
from os import listdir
import re
import os
import seaborn as sns
from scipy.stats import wilcoxon


def read_data(fi):
    data = []
    with open(fi, 'r') as f:
        for line in f:
            if not line.startswith('bin'):
                continue
            # NA to 0
            # data.append([float(x) if x != 'NA' else 0 for x in line.split('\t')[-4:-1]])

            # delete NA
            if 'NA' in line.split('\t')[-4:-1]:
                continue
            else:
                data.append([float(x) for x in line.split('\t')[-4:-1]])
    data = np.asarray(data)
    x = data.T[0]
    y = data.T[1]
    z = data.T[2]
    return x, y, z


def get_resolution(fi):
    with open(fi, 'r') as f:
        for line in f:
            if not line.startswith('bin'):
                continue
            return int(line.split('\t')[2]) - int(line.split('\t')[1])


def plot(fi1, fi2, chrom, threshold=0.2, **kwargs):
    resolution = get_resolution(fi1)
    x1, y1, z1 = read_data(fi1)
    x2, y2, z2 = read_data(fi2)
    x = x1
    plt.figure(figsize=(20, 4))
    plt.plot(x, y1, label='P5')
    plt.plot(x, y2, label='P15')
    plt.plot([0, x[-1]], [0, 0], color='k', linestyle='--')
    # plt.legend()

    if 'diff' in kwargs:
        diff = y1-y2
        plt.plot(x, diff, label='diff', color='k')

    plt.legend()
    axes = plt.gca()
    axes.set_ylim([-2, 2])

    # add boundaires indication
    if 'boundaries' in kwargs:
        data1 = []
        data2 = []
        fi1_boundaries = fi1+'.boundaries.bed'
        fi2_boundaries = fi2+'.boundaries.bed'

        with open(fi1_boundaries, 'r') as f:
            for line in f:
                if not line.startswith('chr'):
                    continue
                tokens = line.rstrip().split('\t')
                if float(tokens[-1]) > threshold:
                    data1.append([int(tokens[1]), int(tokens[2]), float(tokens[-1])])

        with open(fi2_boundaries, 'r') as f:
            for line in f:
                if not line.startswith('chr'):
                    continue
                tokens = line.rstrip().split('\t')
                if float(tokens[-1]) > threshold:
                    data2.append([int(tokens[1]), int(tokens[2]), float(tokens[-1])])

        data1 = np.asarray(data1)
        data2 = np.asarray(data2)
        midpoint1 = (data1.T[1] - data1.T[0])/resolution/2 + data1.T[0]/resolution
        midpoint2 = (data2.T[1] - data2.T[0])/resolution/2 + data2.T[0]/resolution
        values1 = data1.T[-1]
        values2 = data2.T[-1]
        y1 = [-1.8 for x in values1]
        y2 = [-1.9 for x in values2]
        plt.plot(midpoint1, values1, 'ob', label='P5 boundary')
        plt.plot(midpoint2, values2, 'og', label='P15 boundary')
        #for i, txt in enumerate(values1):
        #    plt.annotate(round(txt,2), (midpoint1[i], values1[i]))

        if 'labelBoundaries' in kwargs:
            x = []
            pct_change = []
            with open('/Volumes/Luyang1/hic/domain/20k/motif/data2_more_insulate.bed', 'r') as f:
                for line in f:
                    if re.search('^{}\t'.format(chrom), line):
                        tokens = line.rstrip().split('\t')
                        midpoint = (int(tokens[2]) - int(tokens[1]))/resolution + int(tokens[1])/resolution
                        x.append(midpoint)
                        pct_change.append(round(float(tokens[4]),2))
            y = [1 for a in x]

            for a,b,v in zip(x,y,pct_change):
                c = b+0.5
                plt.annotate(v, xy=(a, b), xytext=(a, c), arrowprops=dict(facecolor='blue', shrink=0.05, alpha=0.5))
            #plt.plot(x, y, 'oy')


    plt.legend()

    if 'region' in kwargs:
        region = np.asarray([float(x) for x in kwargs.pop('region')])/resolution
        axes.set_xlim(list(region))
        start = int(float(region[0])*resolution)
        end = int(float(region[1])*resolution)
        plt.show()
        plt.savefig(os.getcwd()+'/'+'chr{}_{}_{}_curve.png'.format(chrom, start, end), dpi=150)
    else:
        plt.show()
        plt.savefig(os.getcwd()+'/'+'chr{}_curve.png'.format(chrom), dpi=150)

    # plt.close()


def wholegenome(path):
    """
    path: the full path of the fold contain all .insulation file
    :return: figures of each chromosome
    """
    for fi in listdir(path):
        if fi.startswith('P5') and re.search('insulation$', fi):
            print(fi)
            fi1 = fi
            fi2 = re.sub('P5', 'P15', fi)
            chrom = re.search('chr(\d+)', fi1).group(1)
            plot(path+"/"+fi1, path+"/"+fi2, chrom, diff=1)


def single_chr(chrom, **kwargs):
    fi1 = '/Volumes/Luyang1/hic/domain/40k/P5_chr{}_ismatrix.is2000001.ids400001.insulation'.format(chrom)
    fi2 = '/Volumes/Luyang1/hic/domain/40k/P15_chr{}_ismatrix.is2000001.ids400001.insulation'.format(chrom)
    chrom = 'chr'+str(chrom)
    if 'region' in kwargs:
        region = kwargs.pop('region')
        #plot(fi1, fi2, chrom, region=region, diff=1, boundaries=1, labelBoundaries=1)
        plot(fi1, fi2, chrom, region=region, diff=1)
    # TODO: need to change to match kwargs with proper command

    else:
        #plot(fi1, fi2, chrom, diff=1, boundaries=1, labelBoundaries=1)
        plot(fi1, fi2, chrom, diff=1)


def boxplot(fi1, fi2, threshold=0.7, **kwargs):
    """
    make boxplot of the insulation strength for a chr or region
    fi1 and fi2 should be "XXXX.is2000001.ids400001.insulation.boundaries.bed"
    :return:
    """

    # read file
    data1 = []
    data2 = []
    with open(fi1, 'r') as f:
        for line in f:
            if not line.startswith('chr'):
                continue
            tokens = line.rstrip().split('\t')
            if float(tokens[-1]) > threshold:
                data1.append(float(tokens[-1]))

    with open(fi2, 'r') as f:
        for line in f:
            if not line.startswith('chr'):
                continue
            tokens = line.rstrip().split('\t')
            if float(tokens[-1]) > threshold:
                data2.append(float(tokens[-1]))
    data = [data1, data2]
    # boxplot
    #plt.boxplot(data,showfliers=False)
    #plt.plot([0, 3], [np.median(data1), np.median(data1)], color='k', linestyle='--')
    plot = sns.boxplot(data=data, showfliers=False)

    # calculate p-value
    # t, pvalue = wilcoxon(data1, data2)


    # statistical annotation
    """
    http://stackoverflow.com/questions/36578458/how-does-one-insert-statistical-annotations-stars-or-p-values-into-matplotlib
    """
    x1, x2 = 0, 1
    h = ((np.percentile(data1, 75) - np.percentile(data1, 25))*1.5+np.percentile(data1, 75))/8
    y, h, col = max([
        ((np.percentile(data1, 75) - np.percentile(data1, 25))*1.5+np.percentile(data1, 75))
        , ((np.percentile(data2, 75) - np.percentile(data2, 25))*1.5+np.percentile(data2, 75))
                    ]) + h, h, 'k'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*.5, y+h, 'P-value', ha='center', va='bottom', color=col)


    # histogram
    # plt.hist(np.log(data1), 20, alpha=0.4)
    # plt.hist(np.log(data2), 20, alpha=0.4)
    # plt.legend(['P5', 'P15'])
    plt.show()


def test():
    chrom = 19
    fi1 = '/Volumes/Luyang1/hic/domain/40k/P5_chr{}_ismatrix.is2000001.ids400001.insulation.boundaries.bed'.format(chrom)
    fi2 = '/Volumes/Luyang1/hic/domain/40k/P15_chr{}_ismatrix.is2000001.ids400001.insulation.boundaries.bed'.format(chrom)
    boxplot(fi1, fi2, threshold=0.2)

if __name__ == '__main__':
    chrom = 9
    fi1 = '/Volumes/Luyang1/hic/domain/40k/P5_chr{}_ismatrix.is2000001.ids400001.insulation.boundaries.bed'.format(chrom)
    fi2 = '/Volumes/Luyang1/hic/domain/40k/P15_chr{}_ismatrix.is2000001.ids400001.insulation.boundaries.bed'.format(chrom)
    #boxplot(fi1, fi2, threshold=0.2)
    #single_chr(14, region=(0*20000, 1000*20000), diff=1, boundaries=1)
    single_chr(1, region=(230000000,245000000))
    # wholegenome('/Volumes/Luyang1/hic/domain/20k')

