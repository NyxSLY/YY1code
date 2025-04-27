chrom=$1
region=$2
pa=$3

python2 /Users/Luyang/scripts/hic/script/highResRegionPlot.py -chr $chrom -r $region -hm /Volumes/Luyang1/hic/map_result/hm/byChr/noDiag/IC-P${pa}_heatmap-40k.byChr -isfile /Volumes/Luyang1/hic/domain/40k/P${pa}_chr${chrom}_ismatrix.is2000001.ids400001.insulation.boundaries.bed
