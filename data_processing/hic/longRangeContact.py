"""
count the number of long range contract
"""
from hiclib import fragmentHiC
from mirnylib import h5dict, genome
import numpy as np
import matplotlib.pyplot as plt


filename = 'fragment_dataset_P5.hdf5'
path = '/Volumes/Luyang1/hic/map_result/'
genome_db = genome.Genome('/Volumes/LACIE/Human_database/hg19/bowtie2_index/fasta/hg19',
                          gapFile='/Volumes/LACIE/Human_database/hg19/bowtie2_index/fasta/hg19/gap.txt',
                          readChrms=['#', 'X', "Y"])
genome_db.setEnzyme('MboI')

Dset = fragmentHiC.HiCdataset(filename=path+filename, mode='r', genome=genome_db)
maskN = (abs(Dset.cuts1 - Dset.cuts2) < 100000)*(Dset.chrms1 == Dset.chrms2)
N = maskN.sum()
mask = (Dset.chrms1 == Dset.chrms2)
T = mask.sum()

print "%s Cis contact < 100kb: %s" % (filename, N)
print "Total Cis: %s" % (T)


filename = 'fragment_dataset_P15.hdf5'
path = '/Volumes/Luyang1/hic/map_result/'
genome_db = genome.Genome('/Volumes/LACIE/Human_database/hg19/bowtie2_index/fasta/hg19',
                          gapFile='/Volumes/LACIE/Human_database/hg19/bowtie2_index/fasta/hg19/gap.txt',
                          readChrms=['#', 'X', "Y"])
genome_db.setEnzyme('MboI')

Dset = fragmentHiC.HiCdataset(filename=path+filename, mode='r', genome=genome_db)
maskN = (abs(Dset.cuts1 - Dset.cuts2) < 100000)*(Dset.chrms1 == Dset.chrms2)
N = maskN.sum()
mask = (Dset.chrms1 == Dset.chrms2)
T = mask.sum()

print "%s Cis contact < 100kb: %s" % (filename, N)
print "Total Cis: %s" % (T)