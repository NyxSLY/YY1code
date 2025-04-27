# read sample name and fastq information from csv file, and run fastp and hisat2 based on that
# 2022-5-7
# add parallel ability

import pandas as pd
import subprocess as sp
import re
import sys
import logging
import os
import numpy as np
from contextlib import suppress
from queue import Queue
from threading import Thread
import time

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(message)s')
q = Queue(maxsize=0)
num_threads = 10

def job(i, q):
    while True:
        # function here
        start_time = time.time()
        cmds = q.get()
        for cmd in cmds:
            logging.info('worker{}:\n{}\n'.format(i, cmd))
            sp.call(cmd, shell=True)
        q.task_done()
        logging.info('{} FINISHED: time elapsed: {}'.format(cmd, time.time()-start_time))

def map_pipeline(fqfile1, fqfile2, prefix, genome_index):
    trim_output1 = '{}_1.trimmed.fq'.format(prefix)
    trim_output2 = '{}_2.trimmed.fq'.format(prefix)
    trim_cmd = f'fastp -i {fqfile1} -I {fqfile2} -o {trim_output1} -O {trim_output2} -j {prefix}.trim.json -h {prefix}.trim.html'
    map_cmd = f'bowtie2 -x {genome_index} -1 {fqfile1} -2 {fqfile2} | samtools view -bSh - > {prefix}.filtered.bam'
    index_cmd = f'samtools sort {prefix}.filtered.bam -o {prefix}.sorted.filtered.bam; samtools index {prefix}.sorted.filtered.bam'
    cmds = [map_cmd]
    return(cmds)

def map_pipeline_SE(fqfile1, prefix, genome_index):
    trim_output1 = '{}_1.trimmed.fq'.format(prefix)
    trim_cmd = f'fastp -i {fqfile1} -o {trim_output1} -j {prefix}.trim.json -h {prefix}.trim.html'
    map_cmd = f'hisat2 -p 8 -x {genome_index} -U {trim_output1} --summary-file {prefix}.mappingSummary.txt | samtools view -bSh - | samtools sort - -o {prefix}.sorted.bam'

    index_cmd = f'samtools index {prefix}.sorted.bam'
    cmds = [trim_cmd, map_cmd, index_cmd]
    return(cmds)

x = pd.read_csv(sys.argv[1])

fq1_list = x['read1']
fq2_list = x['read2']
prefix_list = x['prefix']
genome_list = x['genome']

queue_list = []
for f1, f2, y, z in zip(fq1_list, fq2_list, prefix_list, genome_list):
    queue_list.append([f1, f2, y, z])

for i in range(num_threads):
    worker = Thread(target=job, args=(i, q,))
    worker.setDaemon(True)
    worker.start()

for i in queue_list:
    if i[1] == '' or pd.isnull(i[1]):
        i[0] = '\"{}\"'.format(i[0])
        a = map_pipeline_SE(i[0], i[2], i[3])
    else:
        i[0] = '\"{}\"'.format(i[0])
        i[1] = '\"{}\"'.format(i[1])
        a = map_pipeline(i[0], i[1], i[2], i[3])
    q.put(a)
q.join()

