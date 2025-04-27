rm(list = ls())
load('D:\\BaiduNetdiskWorkspace\\OneDrive\\work\\3dGenome\\Find_enhancer_target\\BASE.RData')
setwd('D:\\BaiduNetdiskWorkspace\\OneDrive\\work\\3dGenome\\manuscript\\writing\\Weiwei_comments_based_on_GSA\\YY1_redistribution')
library(rtracklayer)
library(GenomicRanges)
library(regioneR)

getwd()
yy1
pros

# TAD boundaries
tad_boundaries_1 = import('P5_wholeGenome.finalboundaries.bed','BED')
tad_boundaries_2 = import('P15_wholeGenome.finalboundaries.bed','BED')

tad_boundaries = GenomicRanges::union(tad_boundaries_2, tad_boundaries_1)
data.frame(tad_boundaries) -> tad_boundaries
tad_boundaries$seqnames = as.character(tad_boundaries$seqnames)
tad_boundaries$seqnames[tad_boundaries$seqnames=='chr23'] = 'chrX'
makeGRangesFromDataFrame(tad_boundaries) -> tad_boundaries


# specific YY1
yy1.p15.sp = yy1_p15[!overlapsAny(yy1_p15, yy1_p5)]
yy1.p5.sp = yy1_p5[!overlapsAny(yy1_p5, yy1_p15)]
yy1.common = yy1_p5[overlapsAny(yy1_p5, yy1_p15)]
export(yy1.common, 'YY1_common.bed','BED')

# stringent universe
univ = import('universe_strigent.bed', format='BED')

data.frame(g$mask) -> mask
mask[!grepl('_', mask$seqnames),] -> mask
makeGRangesFromDataFrame(mask) -> mask

seqlevels(mask, pruning.mode="coarse") <- data.frame(univ)$seqnames %>% unique() %>% as.character()


perm.yy1p15sp.tadb = permTest(A=tad_boundaries, B=yy1.p15.sp,
                          randomize.function=randomizeRegions,
                          genome=univ, mask=mask, evaluate.function=numOverlaps, count.once=T, ntimes=1000)

perm.yy1p5sp.tadb = permTest(A=tad_boundaries, B=yy1.p5.sp,
                         randomize.function=randomizeRegions,
                         genome=univ, mask=mask, evaluate.function=numOverlaps, count.once=T, ntimes=1000)


# SE
# SE also change with age

# In Y, test P5.sp.YY1 and p5 SE(P5 sp SE + common SE)
# In O, test P15.sp.YY1 and p15 SE(P15 sp SE + common SE)
# then compare with them
# In summary, they all compare with their own age


# In Y
perm.yy1p5sp.SE = permTest(A=se[!grepl('p15',se$name),], B=yy1.p5.sp,
                                  randomize.function=randomizeRegions,
                                  genome=univ, mask=mask, evaluate.function=numOverlaps, count.once=T, ntimes=1000)

# In O
perm.yy1p15sp.SE = permTest(A=se[!grepl('p5',se$name),], B=yy1.p15.sp,
                                 randomize.function=randomizeRegions,
                                 genome=univ, mask=mask, evaluate.function=numOverlaps, count.once=T, ntimes=1000)


# use all age SE, not just age-specific SE
perm.yy1p5sp.SE1 = permTest(A=se[!grepl('p15',se$name),], B=yy1_p5,
                           randomize.function=randomizeRegions,
                           genome=univ, mask=mask, evaluate.function=numOverlaps, count.once=T, ntimes=100)

# In O
perm.yy1p15sp.SE1 = permTest(A=se[!grepl('p5',se$name),], B=yy1_p15,
                            randomize.function=randomizeRegions,
                            genome=univ, mask=mask, evaluate.function=numOverlaps, count.once=T, ntimes=100)


# calculate the distribution - could make a table in the supplemental figure.

# yy1.p5.sp
sum(!overlapsAny(yy1.p5.sp, se)) # not on SE
sum(overlapsAny(yy1.p5.sp, se.p5)) # on P5 SE
sum(overlapsAny(yy1.p5.sp, se.p15)) # on P15 SE
sum(overlapsAny(yy1.p5.sp, se.common)) # on common SE


# yy1.p15.sp
sum(!overlapsAny(yy1.p15.sp, se)) # not on SE
sum(overlapsAny(yy1.p15.sp, se.p5)) # on P5 SE
sum(overlapsAny(yy1.p15.sp, se.p15)) # on P15 SE
sum(overlapsAny(yy1.p15.sp, se.common)) # on common SE



plot(perm.yy1p15sp.tadb)
plot(perm.yy1p5sp.tadb)







sum(overlapsAny(yy1.p5.sp, tad_boundaries))/length(yy1.p5.sp)
sum(overlapsAny(yy1.p15.sp, tad_boundaries))/length(yy1.p15.sp)

sum(overlapsAny(yy1.p5.sp, se))/length(yy1.p5.sp)
sum(overlapsAny(yy1.p15.sp, se))/length(yy1.p15.sp)

sum(overlapsAny(yy1.p5.sp, tad_boundaries))/length(yy1.p5.sp)
sum(overlapsAny(yy1.p15.sp, tad_boundaries))/length(yy1.p15.sp)

sum(overlapsAny(yy1_p5, se))/length(se)
sum(overlapsAny(yy1_p15, se))/length(se)







