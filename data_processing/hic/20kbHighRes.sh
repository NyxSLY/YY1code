# 20kb version with region
chr=$1
region=$2

python2 /Users/Luyang/scripts/hic/script/highResRegionPlot.py -hm /Volumes/Luyang1/hic/map_result/hm/byChr/noDiag/IC-P15_heatmap-20k.byChr -chr $chr -isfile /Volumes/Luyang1/hic/domain/20k/P15_chr${chr}_20000_ismatrix.is1000001.ids200001.insulation.boundaries.bed -r $region &

python2 /Users/Luyang/scripts/hic/script/highResRegionPlot.py -hm /Volumes/Luyang1/hic/map_result/hm/byChr/noDiag/IC-P5_heatmap-20k.byChr -chr $chr -isfile /Volumes/Luyang1/hic/domain/20k/P5_chr${chr}_20000_ismatrix.is1000001.ids200001.insulation.boundaries.bed -r $region &

python2 /Users/Luyang/scripts/hic/script/highResRegionPlot.py -hm /Volumes/Luyang1/hic/map_result/hm/byChr/noDiag/IC-P15_heatmap-20k.byChr -chr $chr -isfile /Volumes/Luyang1/hic/domain/20k/P15_chr${chr}_20000_ismatrix.is1000001.ids200001.insulation.boundaries.bed -r $region -s /Volumes/Luyang1/hic/map_result/hm/byChr/noDiag/IC-P5_heatmap-20k.byChr -isfile1 /Volumes/Luyang1/hic/domain/20k/P5_chr${chr}_20000_ismatrix.is1000001.ids200001.insulation.boundaries.bed