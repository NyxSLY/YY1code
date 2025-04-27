"""
count the number of long range contact

UsageL

The script <int(count contact less than this range)>
"""
from hiclib import fragmentHiC
from mirnylib import h5dict, genome
import numpy as np
import matplotlib.pyplot as plt
import sys


length = int(sys.argv[1])
filename = 'fragment_dataset_P5.hdf5'
filename1 = 'fragment_dataset_P15.hdf5'
path = '/Volumes/Luyang1/hic/map_result/'
genome_db = genome.Genome('/Volumes/LACIE/Human_database/hg19/bowtie2_index/fasta/hg19',
                          gapFile='/Volumes/LACIE/Human_database/hg19/bowtie2_index/fasta/hg19/gap.txt',
                          readChrms=['#', 'X', "Y"])
genome_db.setEnzyme('MboI')

Dset5 = fragmentHiC.HiCdataset(filename=path+filename, mode='r', genome=genome_db)
Dset15 = fragmentHiC.HiCdataset(filename=path+filename1, mode='r', genome=genome_db)
data = {}
for i in range(0, 24):
	maskN5 = (abs(Dset5.cuts1 - Dset5.cuts2) < length)*(Dset5.chrms1 == Dset5.chrms2)*(Dset5.chrms1==i)
	maskN15 = (abs(Dset15.cuts1 - Dset15.cuts2) < length)*(Dset15.chrms1 == Dset15.chrms2)*(Dset15.chrms1==i)
	N = maskN5.sum()
	N15 = maskN15.sum()
	i = str(i+1)
	data['chr{}'.format(i)] = (N, N15)
	print "%s chr%s Cis contact < 100kb: %s" % (filename, i, N)
	print "%s chr%s Cis contact < 100kb: %s" % (filename, i, N15)

outfile = 'chr_wide_short_contact_%s.csv' %(length)
with open(outfile, 'w') as f:
	for i in data:
		line = '{},{},{}\n'.format(i, data[i][0], data[i][1])
		f.write(line)

