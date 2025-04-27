"""
Perform filter and IC for high resolution (<=100 kb) data (.byChr)
All set tolerance = 1e-5, same with low reolution setting
only perform removePoorRegion and removeDiagonal, because those the only 2 provided, and the same with example
"""

from hiclib.highResBinnedData import HiResHiC
from mirnylib import genome
import sys
import os
import re


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

db = HiResHiC(genome_db, resolution, name, mode='w')
db.loadData(dictLike=heatmap)
db.removePoorRegions(1)
db.removeDiagonal()
db.iterativeCorrection(1e-5)
db.export('IC-'+name+'.byChr')
