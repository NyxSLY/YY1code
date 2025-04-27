from mirnylib import h5dict, genome
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import re
import os
import os.path
import argparse
import logging
from scipy.stats import poisson


"""
change log:
Mon Oct  8 12:50:23 CDT 2018
1) add one more option, to generate log(P15/P5,2) heatmap
   Also, if 2 isfile provided, label both of them in the ratio heatmap
2) add a filter to processInsulationScoreResult, filter out TADs with low signal. These TADs
   usually are centromere or regions with low mappability

Wed Oct 17 16:23:32 CDT 2018
1) add an quick method to filter TADs, which is based the gap.txt file of hg19, but it's not perfect.
   some regions not in the gaps do not have mapped reads. So this is a quick and dirty way
2) set the orignial method as slow method
3) add a reset option, if reset is invoke, it will make new final boundaires file. This should be used
   when quick method was performed and want to try with slow method again.

Mon Aug 26 12:57:45 CDT 2019
modifed to plot super enhancer around regions
    1) don't need TAD labeling, but need to add line for annotation
    2) plot multiple plots at once, read bed file
    3) when plot ratio, plot pvalue -- line 88 
"""


logging.basicConfig(format='%(message)s', level=logging.INFO)
path = os.getcwd()


def get_parameter():
    parser = argparse.ArgumentParser(description='make region plot from highRes byChr h5dict file')
    parser.add_argument('-hm', required=True, help='the heatmap byChr h5dict file')
    parser.add_argument('-bed', required=True, help='the 6col bed file of SE region')
    parser.add_argument('-width', default=10, type=int,help='number of bins plot up-stream and dn-stream. Default=10')
    parser.add_argument('-r', help='only plot the coords region if provided. Default: None', nargs='+', type=int)
    parser.add_argument('-isfile', help='label tads on the heatmap if provided. File should be generated'
                                        'from insulation score result "xxxxx.insulation.boundaries.bed"')
    parser.add_argument('-is_threshold', type=float, default=0.2, help='the threshold of '
                                                                      'insulation score. Default=0.2')
    parser.add_argument('-b', help='will plot a curve flat under the heatmap, if provided')
    parser.add_argument('-s', type=str, help='generate log(hm heatmap/s heatmap) ratio heatmap.'
                                              'for example, -hm p15.heatmap and -s p5.heatmap will generate '
                                              'log(p15/p5) ratio heatmap')
    parser.add_argument('-isfile1', help='the second isfile used to label tabs on heatmap. Should only invoke'
                                         'when -s is enabled')
    parser.add_argument('-reset', help='if invoked, will make new final boundaries file and replace the existing one')

    return parser.parse_args()

def region_plot(hm_file, chr, resolution, **kwargs):

    if 'region' in kwargs:
        region = kwargs.pop('region')
        region_start = int(region[0])
        region_end = int(region[1])
        start_bin = int(region_start/resolution)
        end_bin = int(region_end/resolution)
    else:
        start_bin = 0
        end_bin = h5dict.h5dict(hm_file)['{} {}'.format(chr, chr)].shape[0]
        region_start = 0
        region_end = end_bin*resolution

    if 'ratio' not in kwargs:
        if re.search('.byChr', hm_file):
            hm = h5dict.h5dict(hm_file, 'r')
            region_hm = hm['{} {}'.format(chr, chr)][start_bin:end_bin, start_bin:end_bin]
        else:  # not byChr, should be matrix directly
            # TODO: detect if h5dict format
            region_hm = np.loadtxt(hm_file)[start_bin:end_bin, start_bin:end_bin]
    else:
        if re.search('.byChr', hm_file):
            hm = h5dict.h5dict(hm_file, 'r')
            region_hm = hm['{} {}'.format(chr, chr)][start_bin:end_bin, start_bin:end_bin]
            hm1 = h5dict.h5dict(args.s, 'r')
            region_hm1 = hm1['{} {}'.format(chr, chr)][start_bin:end_bin, start_bin:end_bin]
            # region_hm = region_hm/region_hm1
            # calculate p value

            di = np.diag_indices(region_hm.shape[0])
            di_1 = (di[0][1:], di[1][:-1]) # line under diagnal
            di_2 = (di[0][:-1], di[1][1:]) # line upon diagnal
            region_hm[di] = 0
            region_hm1[di] = 0
            region_hm[di_1] = 0
            region_hm[di_2] = 0
            region_hm1[di_1] = 0
            region_hm1[di_2] = 0
            mask1 = region_hm>10
            mask2 = region_hm1>10
            maskk = np.any([mask1, mask2], axis=0)
            np.set_printoptions(suppress=True)          

            region_hm = region_hm.astype(int)
            region_hm1 = region_hm1.astype(int)
            np.savetxt('{}_matrix.csv'.format(sample), region_hm, delimiter=',')
            np.savetxt('{}_matrix1.csv'.format(sample), region_hm1, delimiter=',')
            pos_pvalue1 = poisson.sf(region_hm+1,region_hm1+1)
            pos_pvalue2 = poisson.sf(region_hm1+1,region_hm+1)
            pos_pvalue = pos_pvalue1
            mask = np.array(region_hm1)>np.array(region_hm)
            pos_pvalue[mask] = pos_pvalue2[mask]
            pos_pvalue[np.invert(maskk)]= 1
            region_hm = pos_pvalue

        else:  # not byChr, should be matrix directly
            # TODO: detect if h5dict format
            region_hm = np.loadtxt(hm_file)[start_bin:end_bin, start_bin:end_bin]
            region_hm1 = np.loadtxt(args.s)[start_bin:end_bin, start_bin:end_bin]
            region_hm = region_hm/region_hm1

    # only export matrix
    if 'only_matrix' in kwargs:
        return region_hm

    print(region_hm.shape)
    fig = plt.figure(figsize=(8,8))
    ax1 = fig.add_subplot(111)

    # set absolute coords
    # TODO: axis has problem
    #ax1.set_xticks([x for x in xrange(start_bin, end_bin, (end_bin-start_bin)/5)])
    #ax1.set_xticklabels([str(x*resolution/1000)+'kb' for x in xrange(start_bin, end_bin, (end_bin-start_bin)/5)])
    #ax1.set_yticks([x for x in xrange(start_bin, end_bin, (end_bin-start_bin)/5)])
    #ax1.set_yticklabels([str(x*resolution/1000)+'kb' for x in xrange(start_bin, end_bin, (end_bin-start_bin)/5)])
    ax1.spines['right'].set_visible(False)
    ax1.spines['top'].set_visible(False)

    if 'ratio' not in kwargs:
        plot_matrix(ax1, np.log(region_hm), cmap='afmhot_r')
    else:
        plot_matrix(ax1, np.log(region_hm), cmap='afmhot')

    # add TADs
    if 'tads' in kwargs:
        regions = []
        tad_bed = open(kwargs.pop('tads'), 'r')
        for line in tad_bed:
            if not line.startswith('chr'):
                continue
            region = line.rstrip().split('\t')[1:3]
            regions.append(tuple([int(x) for x in region]))
        bin_regions = [np.asarray(x)/resolution for x in regions if x[1] <= region_end]
        addTads(ax1, bin_regions, resolution, region_start)

    if 'tads1' in kwargs:
        regions = []
        tad_bed = open(kwargs.pop('tads1'), 'r')
        for line in tad_bed:
            if not line.startswith('chr'):
                continue
            region = line.rstrip().split('\t')[1:3]
            regions.append(tuple([int(x) for x in region]))
        bin_regions = [np.asarray(x)/resolution for x in regions if x[1] <= region_end]
        addTads(ax1, bin_regions, resolution, region_start, color='r')


    # add line -- label SE
    lableSE(fig, start_bin, end_bin)

    # add diagnoal
    addline(fig, [0, end_bin-start_bin], [0, end_bin-start_bin])

    if 'ratio' in kwargs:
        pass
        plt.savefig('{}_chr{}_{}_{}_{}_ratio.png'.format(sample, chr+1, region_start, region_end, str(resolution/1000)+'k'), dpi=1200)
    else:
        plt.savefig('{}_chr{}_{}_{}_{}.png'.format(sample, chr+1, region_start, region_end, str(resolution/1000)+'k'), dpi=1200)


def read_bedgraph(fil):
    data = {}
    with open(fil, 'r') as f:
        for line in f:
            tokens = line.rstrip().split('\t')
            try:
                data[tokens[0]].append(tokens[1:])
            except KeyError:
                data[tokens[0]] =[]
                data[tokens[0]].append(tokens[1:])
    return data


def addTads(fig, regions, res, region_start, color='black'):
    for region in regions:
        fig.add_patch(
            patches.Rectangle(
                tuple([region[0]-region_start/res, region[0]-region_start/res]),
                region[1]-region[0],
                region[1]-region[0],
                fill=False,
                alpha=0.8,
                linestyle='-',
                linewidth=0.3,
                color=color,

        )
    )


def addline(fig, pos1, pos2, color='black'):
    plt.plot(pos1, pos2, 'k--', lw=0.5)


def lableSE(fig, start_bin, end_bin):
    """
    SE should always be in the middle
    """
    width = end_bin-start_bin
    mid_bin = int((end_bin-start_bin)/2)

    addline(fig, [0, width], [mid_bin, mid_bin])
    addline(fig, [mid_bin, mid_bin], [0, width])


def processInsulationScoreResult_old(result_file, strength_threshold, res):
    """
    process the insulation score result 'xxxxx.insulation.boundaries.bed'
    the last column is insulation strength. Filter the data based on that.
    Default = 0.2
    """
    data = []
    with open(result_file, 'r') as f:
        for line in f:
            if not line.startswith('chr'):
                continue
            tokens = line.rstrip().split('\t')
            if float(tokens[4]) >= strength_threshold:
                bin_s = float(tokens[1])/res
                bin_e = float(tokens[2])/res
                bin_m = (bin_e - bin_s)/2 + bin_s
                data.append([tokens[0], bin_m, tokens[4]])

    # write filtered bed file in 3 column bed format (chr   start   end)
    new_data = []
    for i in range(len(data)-1):
        new_data.append([data[i][0], str(int(float(data[i][1])*res)), str(int(float(data[i+1][1])*res))])

    outfile = '_'.join(result_file.split('.')[:-2])+'.finalDomain.bed'
    os.chdir(os.getcwd())
    with open(outfile, 'w') as f:
        for i in new_data:
            f.write('\t'.join(i)+'\n')

    return outfile


def processInsulationScoreResult_slow(result_file, strength_threshold, res, hm_file, **kwargs):
    """
    process the insulation score result 'xxxxx.insulation.boundaries.bed'
    the last column is insulation strength. Filter the data based on that.
    Default = 0.2

    Changes:
    Tue Oct  9 18:16:07 CDT 2018
    Because I add a step to filter empty TAD, the whole process takes too long time. Now
    add a step to check if the insulation file is already exsit, if so, skip this function

    """
    outfile = '_'.join(result_file.split('.')[:-2])+'.finalDomain.bed'

    if args.reset is not None:
        logging.info('Will generate new boundaires files: {}'.format(outfile))
    else:
        if os.path.isfile(outfile):
            logging.info('File: {} already exit. Do not need to regenerate.'.format(outfile))
            return outfile

    hm = h5dict.h5dict(hm_file, 'r')
    data = []
    with open(result_file, 'r') as f:
        for line in f:
            if not line.startswith('chr'):
                continue
            tokens = line.rstrip().split('\t')
            if float(tokens[4]) >= strength_threshold:
                bin_s = float(tokens[1])/res
                bin_e = float(tokens[2])/res
                bin_m = (bin_e - bin_s)/2 + bin_s
                data.append([tokens[0], bin_m, tokens[4]])

    # write filtered bed file in 3 column bed format (chr   start   end)
    new_data = []
    test_tad_data = {}
    for i in range(len(data)-1):
        tad_bin_s = int(float(data[i][1]))
        tad_bin_e = int(float(data[i+1][1]))
        chrom = int(re.search('chr(\d+)', data[i][0]).group(1)) - 1
        fraction = test_tad(hm, tad_bin_s, tad_bin_e, chrom)
        id = str(chrom) + '_' + str(tad_bin_s) + '_' + str(tad_bin_e)
        test_tad_data[id] = fraction
        if fraction > 0.1:
            new_data.append([data[i][0], str(int(float(data[i][1])*res)), str(int(float(data[i+1][1])*res))])

    os.chdir(path)
    with open(outfile, 'w') as f:
        for i in new_data:
            f.write('\t'.join(i)+'\n')

    for i in test_tad_data:
        print(i+'\t'+str(test_tad_data[i]))
    return outfile


def processInsulationScoreResult(result_file, strength_threshold, res, gap):
    """
    process the insulation score result 'xxxxx.insulation.boundaries.bed'
    the last column is insulation strength. Filter the data based on that.
    Default = 0.2

    Use a simple method to filter the TADs. Import gap.txt
    """
    outfile = '_'.join(result_file.split('.')[:-2])+'.finalDomain.bed'
    data = []
    with open(result_file, 'r') as f:
        for line in f:
            if not line.startswith('chr'):
                continue
            tokens = line.rstrip().split('\t')
            if float(tokens[4]) >= strength_threshold:
                s = float(tokens[1])
                e = float(tokens[2])
                bin_m = (e - s)/2 + s
                data.append([tokens[0], bin_m, tokens[4]])

    # write filtered bed file in 3 column bed format (chr   start   end)
    new_data = []
    test_tad_data = {}
    boundaries_data = []
    for i in range(len(data)-1):
        s = float(data[i][1])
        e = float(data[i+1][1])
        query_region = [s, e]
        chrom = re.search('(chr\d+)', data[i][0]).group(1)

        overlaps = []
        for region in gap[chrom]:
            overlaps.append(getOverlap(query_region, region))
        if len([x for x in overlaps if x > 0.3]) == 0:
            new_data.append([data[i][0], str(int(float(data[i][1]))), str(int(float(data[i+1][1])))])

    os.chdir(path)
    with open(outfile, 'w') as f:
        for i in new_data:
            f.write('\t'.join(i)+'\n')

    for i in test_tad_data:
        print(i+'\t'+str(test_tad_data[i]))
    return outfile


def getOverlap(a, b):
    overlap = max(0, min(a[1], b[1]) - max(a[0], b[0]))
    return float(overlap) / float(a[1]-a[0])


def test_tad(hm, start_bin, end_bin, chrom):
    """
    test if the tad is empty
    :param hm:
    :param start_bin:
    :param end_bin:
    :param chrom:
    :return:
    """
    # data = {}
    region_hm = hm['{} {}'.format(chrom, chrom)][start_bin:end_bin, start_bin:end_bin]
    if np.size(region_hm) == 0:
        return False
    faction = float(np.count_nonzero(region_hm))/np.size(region_hm)
    # id = str(chrom) + '_' + str(start_bin) + '_' + str(end_bin)
    # data[id] = faction
    return faction


def plot_matrix(ax, matrix, **kwargs):
    """Plot a 2D array with a colorbar.

    Parameters
    ----------

    matrix : a 2d numpy array
        A 2d array to plot
    clip_min : float, optional
        The lower clipping value. If an element of a matrix is <clip_min, it is
        plotted as clip_min.
    clip_max : float, optional
        The upper clipping value.
    label : str, optional
        Colorbar label
    ticklabels1 : list, optional
        Custom tick labels for the first dimension of the matrix.
    ticklabels2 : list, optional
        Custom tick labels for the second dimension of the matrix.
    """
    clip_min = kwargs.pop('clip_min', -np.inf)
    clip_max = kwargs.pop('clip_max', np.inf)

    if 'ticklabels1' in kwargs:
        ax.yticks(list(range(matrix.shape[0])))
        ax.gca().set_yticklabels(kwargs.pop('ticklabels1'))

    if 'ticklabels2' in kwargs:
        ax.xticks(list(range(matrix.shape[1])))
        ax.gca().set_xticklabels(kwargs.pop('ticklabels2'))

    img = ax.imshow(
        np.clip(matrix, a_min=clip_min, a_max=clip_max),
        interpolation='nearest',
        **kwargs)
    if 'label' not in kwargs:
        cb = plt.colorbar(img, ax=ax)
    else:
        cb = plt.colorbar(img, ax=ax).set_label(kwargs['label'])

    # hide ticks of colorbar
    if 'hide_colorbar_ticks' in kwargs:
        cb.set_ticks([])


def debug():
    chrom = 1
    hm = '/Volumes/Luyang1/hic/map_result/hm/byChr/IC-P5_heatmap-40k.byChr'
    resolution = h5dict.h5dict(hm)['resolution']
    sample = re.search('(P\d+)', hm).group(1)
    isfile = '/Volumes/Luyang1/hic/domain/40k/P5_chr1_ismatrix.is2000001.ids400001.insulation.boundaries.bed'
    tads = processInsulationScoreResult(isfile, 0.2, resolution)
    region_plot(hm, chrom-1, resolution, region=[230000000,250000000], tads=tads,
                bedgraph='/Volumes/Luyang1/hic/PCA/200k/hasDiagonal/P5/P5-200k-hasDiagonal.PC1.bedGraph')


# make new
if __name__ == '__main__':
    args = get_parameter()
    resolution = h5dict.h5dict(args.hm)['resolution']
    sample = re.search('(P\d+)', args.hm).group(1)
    tads = None
    tads1 = None
    gap = {}
    with open("/Volumes/LACIE/Human_database/hg19/bowtie2_index/fasta/hg19/gap.bed", 'r') as f:
        for line in f:
            a = line.rstrip().split('\t')
            if a[0] == 'chrX':
                a[0] = 'chr23'
            if a[0] =='chrY':
                a[0] = 'chr24'
            if a[0] not in gap:
                gap[a[0]] = []
            gap[a[0]].append([int(a[1]), int(a[2])])

    with open(args.bed, 'r') as f:
        for line in f:
            a = line.rstrip().split('\t')
            if a[0] == 'chrX':
                a[0] = 'chr23'
            if a[0] =='chrY':
                a[0] = 'chr24'
            chrom = int(re.sub('chr','',a[0]))
            start = int(a[1])
            end = int(a[2])
            name = a[3]
            region = [start-args.width*resolution, end+args.width*resolution]
            print(chrom, region)
            chrom = chrom - 1

            if args.s is None:
                if args.isfile is not None:
                    # tads = processInsulationScoreResult(args.isfile, args.is_threshold, resolution, gap)
                    tads = processInsulationScoreResult_slow(args.isfile, args.is_threshold, resolution, args.hm)
                # if args.isfile1 is not None:
                    # tads1 = processInsulationScoreResult(args.isfile1, args.is_threshold, resolution, gap)
                #    tads1 = processInsulationScoreResult_slow(args.isfile1, args.is_threshold, resolution, args.s)
                if tads is None and tads1 is None:
                    region_plot(args.hm, chrom, resolution, region=region)
                if tads is None and tads1 is not None:
                    logging.warning('isfile1 is invoked without -s, exiting')
                    exit()
                    # region_plot(args.hm, chrom, resolution, region=args.r, tads1=tads1)
                if tads is not None and tads1 is None:
                    region_plot(args.hm, chrom, resolution, region=region, tads=tads)
                if tads is not None and tads1 is not None:
                    logging.warning('isfile1 is invoked without -s, exiting')
                    exit()
                    # region_plot(args.hm, chrom, resolution, region=args.r, tads=tads, tads1=tads1)

            else:
                if args.isfile is not None:
                    # tads = processInsulationScoreResult(args.isfile, args.is_threshold, resolution, gap)
                    tads = processInsulationScoreResult_slow(args.isfile, args.is_threshold, resolution, args.hm)
                if args.isfile1 is not None:
                    # tads1 = processInsulationScoreResult(args.isfile1, args.is_threshold, resolution, gap)
                    tads1 = processInsulationScoreResult_slow(args.isfile1, args.is_threshold, resolution, args.s)

                if tads is None and tads1 is None:
                    region_plot(args.hm, chrom, resolution, region=region, ratio=True)
                if tads is None and tads1 is not None:
                    region_plot(args.hm, chrom, resolution, region=region, tads1=tads1, ratio=True)
                if tads is not None and tads1 is None:
                    region_plot(args.hm, chrom, resolution, region=region, tads=tads, ratio=True)
                if tads is not None and tads1 is not None:
                    region_plot(args.hm, chrom, resolution, region=region, tads=tads, tads1=tads1, ratio=True)


