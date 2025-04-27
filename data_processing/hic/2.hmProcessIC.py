import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import re
from mirnylib import genome
from mirnylib import h5dict
from mirnylib import plotting
from hiclib.binnedData import binnedDataAnalysis, binnedData, experimentalBinnedData
import mirnylib.systemutils
import mirnylib.plotting
import scipy.stats
import scipy.ndimage
import logging
cr = scipy.stats.spearmanr


"""
usage: script.py <hm>
"""


logging.basicConfig(level=logging.DEBUG)


# genome_path = '/media/wd/Data/Luyang/hg19/new/fasta/hg19'
raw_hm = sys.argv[1]
path = '/'.join(raw_hm.split('/')[:-1]) + '/'
if path == '/':
    path = os.getcwd() + '/'
prefix = raw_hm.split('/')[-1].split('.')[0]
name = re.search(('P\d+'), raw_hm).group(0)

# parameters
intermedia = False
heatmap = False
remove_diagonal = False

genome_db = genome.Genome('/Volumes/LACIE/Human_database/hg19/bowtie2_index/fasta/hg19',
                          gapFile='/Volumes/LACIE/Human_database/hg19/bowtie2_index/fasta/hg19/gap.txt',
                          readChrms=['#', 'X', "Y"])
genome_db.setEnzyme('MboI')

# Read resolution from the dataset.
raw_heatmap = h5dict.h5dict(raw_hm, mode='r')
resolution = int(raw_heatmap['resolution'])

if resolution <= 100*1000:
    # do not generate whole genome heatmap. No enough memory
    heatmap = False

# Create a binnedData object, load the data.
BD = binnedData(resolution, genome_db, readChrms=["#", "X", "Y"])
BD.simpleLoad(raw_hm, name)

# Plot the raw heatmap.
if heatmap is True:
    plt.figure()
    plotting.plot_matrix(np.log(BD.dataDict[name]))
    plt.savefig(path+name+'raw_{}.png'.format(prefix))

# Remove 1% of regions with low coverage.
BD.removePoorRegions(cutoff=1)

# Plot the raw heatmap directly.
if intermedia is True and heatmap is True:
    plt.figure()
    plotting.plot_matrix(np.log(BD.dataDict[name]), clip_min=0)
    plt.savefig(path+name+'_PR_{}.png'.format(prefix))

# Remove the contacts between loci located within the same bin.
if remove_diagonal is True:
    BD.removeDiagonal()

# Plot the heatmap.
if intermedia is True and heatmap is True:
    plt.figure()
    plotting.plot_matrix(np.log(BD.dataDict[name]), clip_min=0)
    plt.savefig(path+name+'PR_RD_{}.png'.format(prefix))

# Remove bins with less than half of a bin sequenced.
BD.removeBySequencedCount(0.5)

# Plot the heatmap.
if intermedia is True and heatmap is True:
    plt.figure()
    plotting.plot_matrix(np.log(BD.dataDict[name]), clip_min=0)
    plt.savefig(path+name+'PR_RD_RC_{}.png'.format(prefix))

# Truncate top 0.05% of interchromosomal counts (possibly, PCR blowouts).
BD.truncTrans(high=0.0005)

# Plot the heatmap.
if intermedia is True and heatmap is True:
    plt.figure()
    plotting.plot_matrix(np.log(BD.dataDict[name]), clip_min=0)
    plt.savefig(path+name+'PR_RD_RC_RT{}.png'.format(prefix))

# Perform iterative correction.
BD.iterativeCorrectWithoutSS(tolerance=1e-20)

# Save the iteratively corrected heatmap.
outfile = path+'IC-{}.hm'.format(prefix)
print(outfile)
BD.export(name, outfile)

# Plot the corrected heatmap directly.
if heatmap is True:
    plt.figure()
    plotting.plot_matrix(np.log(BD.dataDict[name]))
    plt.savefig(path+'_IC_{}.png'.format(prefix))