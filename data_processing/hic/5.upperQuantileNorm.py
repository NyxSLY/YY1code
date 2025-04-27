import numpy as np
from mirnylib import h5dict
import sys
import re
import os


heatmap_file1 = sys.argv[1]
heatmap_file2 = sys.argv[2]
#heatmap_file = '/Volumes/Luyang1/hic/map_result/hm/byChr/IC-P15_heatmap-40k.byChr'


def calculate_f(heatmap_file):
    hm = h5dict.h5dict(heatmap_file, 'r')
    hm_list = list(hm)[:-1]
    all = []
    for i in hm_list:
        mask = hm[i] > 0
        filtered = hm[i][mask].reshape(1, -1)
        if i == '0 0':
            all = filtered
        else:
            all = np.concatenate((all, filtered), axis=1)
        print(i)
    f = np.percentile(all, 75)
    print(f)
    return f


hm = h5dict.h5dict(heatmap_file1, 'r')
if re.search('.hm$', heatmap_file1):
    name = '.'.join(os.path.basename(1).split('.')[:-1])
    shapee = hm['heatmap'].shape
    one_D_heatmap = hm['heatmap'].reshape(1, shapee[0]*shapee[1])
    mask = one_D_heatmap > 0
    non_0_heatmap = one_D_heatmap[mask]
    upper_quartile_factor = np.percentile(non_0_heatmap, 75)
    normalized_one_D_heatmap = one_D_heatmap / upper_quartile_factor
    normalized_heatmap = normalized_one_D_heatmap.reshape(shapee[0], shapee[1])
    new_hm = h5dict.h5dict(name+'.norm.hm', 'w')
    new_hm.update(hm)
    new_hm['heatmap'] = normalized_heatmap

elif re.search('.byChr$', heatmap_file1):
    name = '.'.join(os.path.basename(heatmap_file1).split('.')[:-1])
    f1 = calculate_f(heatmap_file1)
    f2 = calculate_f(heatmap_file2)

    hm = h5dict.h5dict(heatmap_file1, 'r')
    hm_list = list(hm)
    mydict = h5dict.h5dict(name+'.norm.byChr')
    for i in hm_list:
        data = hm[i]/f1 * float(f1+f2) * 0.5
        mydict[i] = data
    try:
        mydict["resolution"] = hm["resolution"]
    except KeyError:
        mydict["resolution"] = int(re.search(".*-(\d+)k", name).group(1))*1000

    name = '.'.join(os.path.basename(heatmap_file2).split('.')[:-1])
    hm = h5dict.h5dict(heatmap_file2, 'r')
    hm_list = list(hm)[:-1]
    mydict = h5dict.h5dict(name+'.norm.byChr')
    for i in hm_list:
        data = hm[i]/f2 * float(f1+f2) * 0.5
        mydict[i] = data
    try:
        mydict["resolution"] = hm["resolution"]
    except KeyError:
        mydict["resolution"] = int(re.search(".*-(\d+)k", name).group(1))*1000