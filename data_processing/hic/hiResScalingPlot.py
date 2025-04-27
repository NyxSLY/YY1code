from hiclib.highResBinnedData import HiResHiC
from mirnylib import genome
import matplotlib.pyplot as plt
import sys
import re
import numpy as np
from mirnylib import numutils
import logging
import mirnylib
from mirnylib import plotting



"""
make sacling plot from high resolution heatmap
"""

logging.basicConfig(level='INFO')


def plotScaling(res, label="BLA", color=None, plotUnit=1000):
        "plots scaling of a heatmap,treating arms separately"

        data = res.data
        bins = numutils.logbins(
            2, res.genome.maxChrmArm // res.resolution, 1.17)

        chroms = []
        masks = []
        for i in range(res.chromosomeCount):
            logging.info('Start chr{}'.format(i+1))
            chr_data = data[(i, i)].getData()
            s = np.sum(chr_data, axis=0) > 0
            mask = s[:, None] * s[None, :]
            if i != 0:
                end = res.centromerePositions[i] - res.chromosomeEnds[i-1]
            else:
                end = res.centromerePositions[i]
            chroms.append(chr_data[:end, :end])
            masks.append(mask[:end, :end])

            beg = end
            chroms.append(chr_data[beg:, beg:])
            masks.append(mask[beg:, beg:])
        observed = []
        expected = []
        for i in range(len(bins) - 1):
            low = bins[i]
            high = bins[i + 1]
            obs = 0
            exp = 0
            for j in range(len(chroms)):
                if low > len(chroms[j]):
                    continue
                high2 = min(high, len(chroms[j]))
                for k in range(low, high2):
                    obs += np.sum(np.diag(chroms[j], k))
                    exp += np.sum(np.diag(masks[j], k))
            observed.append(obs)
            expected.append(exp)
        observed = np.array(observed, float)
        expected = np.array(expected, float)
        values = observed / expected
        bins = np.array(bins, float)
        bins2 = 0.5 * (bins[:-1] + bins[1:])
        norm = np.sum(values * (bins[1:] - bins[:-1]) * (
            res.resolution / float(plotUnit)))
        args = [res.resolution * bins2 / plotUnit, values / (1. * norm)]
        if color is not None:
            args.append(color)
        plt.plot(*args, label=label, linewidth=2)


def single_chr_plotScaling(res, chrom, resolution, label="BLA", color=None, plotUnit=1000):
        "plots scaling of a heatmap,treating arms separately"

        data = res.data
        i = chrom
        chrmArmLength = max(res.centromerePositions[i] - res.chromosomeStarts[i],
                            res.chromosomeEnds[i] - res.centromerePositions[i])*resolution

        bins = numutils.logbins(
            2, chrmArmLength // resolution, 1.17)

        chroms = []
        masks = []
        i = chrom
        logging.info('Start chr{}'.format(i+1))
        chr_data = data[(i, i)].getData()
        s = np.sum(chr_data, axis=0) > 0
        mask = s[:, None] * s[None, :]
        end = res.centromerePositions[i] - res.chromosomeStarts[i]
        chroms.append(chr_data[:end, :end])
        masks.append(mask[:end, :end])

        beg = end
        chroms.append(chr_data[beg:, beg:])
        masks.append(mask[beg:, beg:])
        observed = []
        expected = []
        for i in range(len(bins) - 1):
            low = bins[i]
            high = bins[i + 1]
            obs = 0
            exp = 0
            for j in range(len(chroms)):
                if low > len(chroms[j]):
                    continue
                high2 = min(high, len(chroms[j]))
                for k in range(low, high2):
                    obs += np.sum(np.diag(chroms[j], k))
                    exp += np.sum(np.diag(masks[j], k))
            observed.append(obs)
            expected.append(exp)
        observed = np.array(observed, float)
        expected = np.array(expected, float)
        print(observed)
        print(expected)
        values = observed / expected
        bins = np.array(bins, float)
        bins2 = 0.5 * (bins[:-1] + bins[1:])
        norm = np.sum(values * (bins[1:] - bins[:-1]) * (
            resolution / float(plotUnit)))
        args = [resolution * bins2 / plotUnit, values / (1. * norm)]
        if color is not None:
            args.append(color)
        plt.plot(*args, label=label, linewidth=2)


def single_chrom_plot(db, db15, resolution):
    """
    make scaling plots of each chromosome
    :return: saved figs
    """

    for i in range(db.chromosomeCount):
        figname = 'chr{} High res distance-dependent interaction frequency.png'.format(i+1)
        plt.figure(i)
        single_chr_plotScaling(db, i, resolution, label="P5", color="#A7A241", plotUnit=1000)
        single_chr_plotScaling(db15, i, resolution, label="P15", color="#344370", plotUnit=1000)
        ax = plt.gca()
        mirnylib.plotting.removeAxes()
        fs = 6
        plt.xlabel("Genomic distance (KB)", fontsize=6)
        plt.ylabel("Contact probability", fontsize=6)
        for xlabel_i in ax.get_xticklabels():
            xlabel_i.set_fontsize(fs)
        for xlabel_i in ax.get_yticklabels():
            xlabel_i.set_fontsize(fs)
        legend = plt.legend(loc=0, prop={"size": 6})
        legend.draw_frame(False)
        plt.xscale("log")
        plt.yscale("log")
        plt.savefig(figname, dpi=300)

        ax = plt.gca()
        mirnylib.plotting.removeAxes()
        fs = 6
        plt.xlabel("Genomic distance (KB)", fontsize=6)
        plt.ylabel("Contact probability", fontsize=6)
        for xlabel_i in ax.get_xticklabels():
            xlabel_i.set_fontsize(fs)
        for xlabel_i in ax.get_yticklabels():
            xlabel_i.set_fontsize(fs)
        legend = plt.legend(loc=0, prop={"size": 6})
        legend.draw_frame(False)
        plt.xscale("log")
        plt.yscale("log")
        plt.savefig(figname, dpi=300)


def whole_genome_plot():
    """
    make scaling plot with combinded data of all chromosomes
    :return: a saved fig
    """
    plotScaling(db, label="P5", color="#A7A241", plotUnit=1000)
    plotScaling(db15, label="P15", color="#344370",plotUnit=1000)

    ax = plt.gca()
    mirnylib.plotting.removeAxes()
    fs = 6
    plt.xlabel("Genomic distance (KB)", fontsize=6)
    plt.ylabel("Contact probability", fontsize=6)
    for xlabel_i in ax.get_xticklabels():
        xlabel_i.set_fontsize(fs)
    for xlabel_i in ax.get_yticklabels():
        xlabel_i.set_fontsize(fs)
    legend = plt.legend(loc=0, prop={"size": 6})
    legend.draw_frame(False)
    plt.xscale("log")
    plt.yscale("log")
    plt.savefig(figname, dpi=300)


figname = 'High res distance-dependent interaction frequency.png'
genome_db = genome.Genome('/Volumes/LACIE/Human_database/hg19/bowtie2_index/fasta/hg19',
                      gapFile='/Volumes/LACIE/Human_database/hg19/bowtie2_index/fasta/hg19/gap.txt',
                      readChrms=['#', 'X'])
genome_db.setEnzyme('MboI')
heatmap = sys.argv[1]
name = heatmap.split('/')[-1].split('.')[0]

# check suffix
if heatmap.split('.')[-1] != 'byChr':
    raise IOError("{} do not have byChr suffix, maybe its not high resolution by chromosome h5dict format"\
        .format(heatmap))

# Read resolution
resolution = int(re.search('-(\d+)k', name).group(1))*1000
logging.info('resolution = {}'.format(resolution))
db = HiResHiC(genome_db, resolution, name, mode='w')
logging.info('Start to load heatmap {}'.format(heatmap))
db.loadData(dictLike=heatmap)

heatmap15 = sys.argv[2]
name15 = heatmap15.split('/')[-1].split('.')[0]
db15 = HiResHiC(genome_db, resolution, name15, mode='w')
logging.info('Start to load heatmap {}'.format(heatmap))
db15.loadData(dictLike=heatmap15)

single_chrom_plot(db, db15, resolution)