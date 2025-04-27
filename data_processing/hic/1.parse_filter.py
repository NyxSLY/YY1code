import os
import logging
from hiclib import mapping
from mirnylib import h5dict, genome
from hiclib import fragmentHiC


logging.basicConfig(level=logging.DEBUG)


wholeGenomeResolutionsKb = [1000, 500, 200, 100]
byChromosomeResolutionsKb = [40, 20]

outpath = '/Volumes/Luyang1/hic/tmp/senecence_Criscione2016/sen/'
sample = 'HiC-Sen-rep2'  # should be P5 or P15

genome_db = genome.Genome('/Volumes/LACIE/Human_database/hg19/bowtie2_index/fasta/hg19',
                          gapFile='/Volumes/LACIE/Human_database/hg19/bowtie2_index/fasta/hg19/gap.txt',
                          readChrms=['#', 'X'])
genome_db.setEnzyme('HindIII')


# B. Parse the mapped sequences into a Python data structure,
#    assign the ultra-sonic fragments to restriction fragments.

mapped_reads = h5dict.h5dict(outpath+'{}_mapped_reads.hdf5'.format(sample))
mapping.parse_sam(
    sam_basename1=outpath+'{}_1.bam'.format(sample),
    sam_basename2=outpath+'{}_2.bam'.format(sample),
    out_dict=mapped_reads,
    genome_db=genome_db,
    enzyme_name='HindIII')


# fragment level filter
fragments = fragmentHiC.HiCdataset(
    filename=outpath+'fragment_dataset_{}.hdf5'.format(sample),
    genome=genome_db,
    maximumMoleculeLength=600,
    mode='w')

if os.path.exists(outpath+'{}_mapped_reads.hdf5'.format(sample)):
    fragments.parseInputData(dictLike=outpath+'{}_mapped_reads.hdf5'.format(sample))
else:
    if sample.isupper():
        if os.path.exists(outpath+'{}_mapped_reads.hdf5'.format(sample.lower())):
            fragments.parseInputData(dictLike=outpath+'{}_mapped_reads.hdf5'.format(sample.lower()))
        else:
            raise IOError('\nNo parsed mapping result avaialbe!')
    else:
        if os.path.exists(outpath+'{}_mapped_reads.hdf5'.format(sample.upper())):
            fragments.parseInputData(dictLike=outpath+'{}_mapped_reads.hdf5'.format(sample.upper()))
        else:
            raise IOError('\nNo parsed mapping result avaialbe!')


# filter
fragments.filterRsiteStart(offset=5)
fragments.filterDuplicates()
fragments.filterLarge(100000, 10)
fragments.filterExtreme(cutH=0.001, cutL=0)

# print filter stat
fragments.writeFilteringStats()
fragments.printMetadata(saveTo=outpath+"{}_statistics.txt".format(sample))

fragments._sortData()

# save heatmap
outfile = outpath+'heatmap-{}k.hm'
for res in wholeGenomeResolutionsKb:
    fragments.saveHeatmap(outpath+'heatmap-{}k.hm'.format(res), res*1000)
for res in byChromosomeResolutionsKb:
    fragments.saveByChromosomeHeatmap(outpath+'heatmap-{}k.byChr'.format(res), res*1000)
