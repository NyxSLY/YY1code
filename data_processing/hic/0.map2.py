import os
import logging
from hiclib import mapping

logging.basicConfig(level=logging.DEBUG)

if not os.path.exists('/Volumes/Luyang1/hic/tmp/senecence_Criscione2016/temp/'):
    os.mkdir('/Volumes/Luyang1/hic/tmp/senecence_Criscione2016/temp/')

# A. Map the reads iteratively.
mapping.iterative_mapping(
    bowtie_path='/Users/Luyang/Downloads/bowtie2-2.2.4/bowtie2',
    bowtie_index_path='/Volumes/LACIE/Human_database/hg19/bowtie2_index/hg19',
    fastq_path='/Volumes/Luyang1/hic/tmp/senecence_Criscione2016/HiC-Sen-rep2_2.fastq',
    out_sam_path='/Volumes/Luyang1/hic/tmp/senecence_Criscione2016/HiC-Sen-rep2_2.bam',
    min_seq_len=25,
    len_step=5,
    seq_start=0,
    seq_end=100,
    nthreads=4,  # on intel corei7 CPUs 4 threads are as fast as
                 # 8, but leave some room for you other applications
    #max_reads_per_chunk = 10000000,  #optional, on low-memory machines
    temp_dir='/Volumes/Luyang1/hic/tmp/senecence_Criscione2016/temp/',  # optional, keep temporary files here
    bowtie_flags='--very-sensitive')
