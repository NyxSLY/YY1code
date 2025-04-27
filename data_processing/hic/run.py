from plots import plots
from hiclib.binnedData import binnedDataAnalysis, binnedData
import matplotlib.pyplot as plt
from mirnylib import h5dict, genome
import mirnylib
from mirnylib import plotting


def singleResHeatmap():
    # for one heatmap
    heatmap = '/Volumes/Luyang1/hic/map_result/hm/P5_heatmap-200k.hm'
    db = plots(heatmap)
    db.singleChrHeatmap()


def wholeGenomePlot():
    heatmap = '/Volumes/Luyang1/hic/map_result/hm/remvoeDiagonal/IC-P5_heatmap-1000k.hm'
    db = plots(heatmap)
    db.wholeGenomeHeatmap()


def heatmapAllRes():
    # for all resolution, and young and old
    samples = ['P5', 'P15']
    resolutions = ['1000k', '500k', '200k', '100k']

    for sample in samples:
        for res in resolutions:
            heatmap = 'IC-{}_heatmap-{}.hm'.format(sample, res)
            db = plots(heatmap)
            db.singleChrHeatmap()


def scalingPlot():
    "Distance-dependent interaction frequency"
    genome_db = genome.Genome('/Volumes/LACIE/Human_database/hg19/bowtie2_index/fasta/hg19',
                          gapFile='/Volumes/LACIE/Human_database/hg19/bowtie2_index/fasta/hg19/gap.txt',
                          readChrms=['#', 'X', "Y"])
    genome_db.setEnzyme('MboI')
    hmP5 = '/Volumes/Luyang1/hic/map_result/hm/remvoeDiagonal/IC-P5_heatmap-200k.hm'
    hmP15 = '/Volumes/Luyang1/hic/map_result/hm/remvoeDiagonal/IC-P15_heatmap-200k.hm'
    rawhmP5 = '/Volumes/Luyang1/hic/map_result/hm/P5_heatmap-500k.hm'
    rawhmP15 = '/Volumes/Luyang1/hic/map_result/hm/P15_heatmap-500k.hm'

    heatmapP5 = h5dict.h5dict(hmP5, mode='r')
    resolutionP5 = int(heatmapP5['resolution'])

    heatmapP15 = h5dict.h5dict(hmP15, mode='r')
    resolutionP15 = int(heatmapP15['resolution'])

    figname = 'distance-dependent interaction frequency.png'
    plt.figure(figsize=(4, 4))
    p5 = binnedDataAnalysis(resolutionP5, genome_db, ["#", "X", "Y"])
    p15 = binnedDataAnalysis(resolutionP15, genome_db, ["#", "X", "Y"])
    # rawp5 = binnedDataAnalysis(resolutionP5, genome_db, ["#", "X", "Y"])
    # rawp15 = binnedDataAnalysis(resolutionP5, genome_db, ["#", "X", "Y"])
    p5.simpleLoad(hmP5, "P5")
    p15.simpleLoad(hmP15, "P15")
    # rawp5.simpleLoad(rawhmP5, 'raw P5')
    # rawp15.simpleLoad(rawhmP15, 'raw P15')
    """
    p5.removePoorRegions(cutoff=1)
    p5.removeDiagonal()
    p5.removeBySequencedCount(0.5)
    p5.truncTrans(high=0.0005)
    p5.iterativeCorrectWithoutSS()

    p15.removePoorRegions(cutoff=1)
    p15.removeDiagonal()
    p15.removeBySequencedCount(0.5)
    p15.truncTrans(high=0.0005)
    p15.iterativeCorrectWithoutSS()
    """
    p5.plotScaling("P5", label="P5", color="#A7A241", plotUnit=1000,)
    p15.plotScaling("P15", label="P15", color="#344370",plotUnit=1000)
    # rawp5.plotScaling("raw P5", label="P5", color="blue")
    # rawp15.plotScaling("raw P15", label="P15", color="red")

    ax = plt.gca()
    mirnylib.plotting.removeAxes()
    fs = 6
    plt.xlabel("Genomic distance (KB)", fontsize=6)
    plt.ylabel("Contact probability", fontsize=6)
    for xlabel_i in ax.get_xticklabels():
        xlabel_i.set_fontsize(fs)
    for xlabel_i in ax.get_yticklabels():
        xlabel_i.set_fontsize(fs)
    legend = plt.legend(loc=0, prop={"size": 6})
    legend.draw_frame(False)
    plt.xscale("log")
    plt.yscale("log")
    plt.savefig(figname, dpi=300)

if __name__ == '__main__':
    scalingPlot()
