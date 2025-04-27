"""
normalize the total read depth to 1M
"""


import numpy as np
from mirnylib import h5dict
import sys
import re
import os


#heatmap_file1 = sys.argv[1]
#heatmap_file2 = sys.argv[2]

heatmap_file1 = '/Volumes/Luyang1/hic/map_result/hm/byChr/P5_heatmap-40k.byChr'
heatmap_file2 = '/Volumes/Luyang1/hic/map_result/hm/byChr/P15_heatmap-40k.byChr'
#heatmap_file = '/Volumes/Luyang1/hic/map_result/hm/byChr/IC-P15_heatmap-40k.byChr'


def calculate_f(heatmap_file):
    hm = h5dict.h5dict(heatmap_file, 'r')
    hm_list = list(hm)
    try:
        hm_list.remove('resolution')
    except:
        pass
    all = []
    for i in hm_list:
        summ = hm[i].sum()
        if i == '0 0':
            all = summ
        else:
            all = all + summ
        print(i)
    f = all
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
    if f1 > f2:
        f = float(f2)/float(f1)

        hm = h5dict.h5dict(heatmap_file1, 'r')
        hm_list = list(hm)
        try:
            hm_list.remove('resolution')
        except:
            pass
        mydict = h5dict.h5dict(name+'.totalnorm.byChr')

        for i in hm_list:
            data = hm[i]
            mydict[i] = data
        try:
            mydict["resolution"] = hm["resolution"]
        except KeyError:
            mydict["resolution"] = int(re.search(".*-(\d+)k", name).group(1))*1000
    #--------------------------------------------------------------------------------------
        name = '.'.join(os.path.basename(heatmap_file2).split('.')[:-1])
        hm = h5dict.h5dict(heatmap_file2, 'r')
        hm_list = list(hm)
        try:
            hm_list.remove('resolution')
        except:
            pass

        mydict = h5dict.h5dict(name+'.totalnorm.byChr')
        for i in hm_list:
            data = hm[i]/f
            mydict[i] = data
        try:
            mydict["resolution"] = hm["resolution"]
        except KeyError:
            mydict["resolution"] = int(re.search(".*-(\d+)k", name).group(1))*1000

    else:
        f = float(f1)/float(f2)

        hm = h5dict.h5dict(heatmap_file1, 'r')
        hm_list = list(hm)
        try:
            hm_list.remove('resolution')
        except:
            pass
        mydict = h5dict.h5dict(name+'.totalnorm.byChr')

        for i in hm_list:
            data = hm[i]/f
            mydict[i] = data
        try:
            mydict["resolution"] = hm["resolution"]
        except KeyError:
            mydict["resolution"] = int(re.search(".*-(\d+)k", name).group(1))*1000
    #--------------------------------------------------------------------------------------
        name = '.'.join(os.path.basename(heatmap_file2).split('.')[:-1])
        hm = h5dict.h5dict(heatmap_file2, 'r')
        hm_list = list(hm)
        try:
            hm_list.remove('resolution')
        except:
            pass

        mydict = h5dict.h5dict(name+'.totalnorm.byChr')
        for i in hm_list:
            data = hm[i]
            mydict[i] = data
        try:
            mydict["resolution"] = hm["resolution"]
        except KeyError:
            mydict["resolution"] = int(re.search(".*-(\d+)k", name).group(1))*1000