"""
following 3.highResIC.py, should combine to one script eventually.

Exact matrix file from h5dict format to txt format. Also add 3 columns (chr, start, end) at the start of each line
in order to agree with homer format and run TADs caller from Ren Lab
"""

from mirnylib.h5dict import h5dict
import numpy as np
import re


#TODO: option: export all combination, export all cis, export all trans
#TODO: option: export specific chr-chr matirx

heatmap_file = '/Volumes/Luyang1/hic/map_result/hm/P5_heatmap-10k.byChr'
heatmap = h5dict(heatmap_file, mode='r')
try:
    resolution = int(heatmap['resolution'])
except KeyError:
    resolution = int(re.search('-(\d+)k', heatmap_file.split("/")[-1]).group(1)) * 1000
chr21 = heatmap['17 17']  # chrs are 0-based
bins = chr21.shape[0]
chrom = 18
path = '/Volumes/Luyang1/hic/map_result/hm/'
c1 = np.full((bins, 1), chrom, dtype=int)
c2 = np.arange(bins).reshape(bins, 1) * resolution
c3 = c2 + resolution
row_head = np.c_[c1, c2, c3]
new_data = np.c_[c1, c2, c3, chr21]
np.savetxt(path+'chr18.matrix', new_data, delimiter='\t', fmt='%i\t%i\t%i'+'\t%g'*bins)
