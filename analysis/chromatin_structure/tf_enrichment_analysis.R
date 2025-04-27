setwd(dirname(rstudioapi::getSourceEditorContext()$path))
load('BASE.RData')



# make universe
union(p5_target_region, p15_target_region) -> all_target_region
all_target_region_extend=GenomicRanges::resize(all_target_region, 2000000)
all_target_region_extend=GenomicRanges::shift(all_target_region_extend, -1000000)
all_target_region_extend = intersect(all_target_region_extend, g$genome)
all_target_region_extend = reduce(all_target_region_extend)
bin = 40*1000
windows.s = slidingWindows(all_target_region_extend, width=bin, step=bin)
windows.s = unlist(windows.s)
windows.s
export(windows.s, 'universe_strigent.bed', 'BED')

univ = windows.s

# perm_ctcf = permTest(A=all_target_region, B=ctcf, 
#                randomize.function=randomizeRegions, 
#                genome=univ, mask=g$mask, evaluate.function=numOverlaps, count.once=T, ntimes=100)
# 
# perm_yy1 = permTest(A=all_target_region, B=yy1, 
#                randomize.function=randomizeRegions, 
#                genome=univ, mask=g$mask, evaluate.function=numOverlaps, count.once=T, ntimes=100)
# 
# plot(perm_ctcf)
# perm_ctcf
# plot(perm_yy1)
# perm_yy1


union(makeGRangesFromDataFrame(se_target_pro.p5),makeGRangesFromDataFrame(se_target_pro.p15)) -> all_target_pro


perm_ctcf = permTest(A=all_target_pro, B=ctcf, 
                     randomize.function=randomizeRegions, 
                     genome=univ, mask=g$mask, evaluate.function=numOverlaps, count.once=T, ntimes=100)

perm_yy1 = permTest(A=all_target_pro, B=yy1, 
                    randomize.function=randomizeRegions, 
                    genome=univ, mask=g$mask, evaluate.function=numOverlaps, count.once=T, ntimes=100)

plot(perm_ctcf)
perm_ctcf
plot(perm_yy1)
perm_yy1



output_ppt(plot(perm_ctcf), 'ctcf enrichment score.pptx')
output_ppt(plot(perm_yy1), 'YY1 enrichment score.pptx')


#-----------------------------------------------
# enhancer region - YY1 is more enriched than CTCF

se_region = import('super_enhancer.bed','BED')

perm_ctcf.se = permTest(A=se_region, B=ctcf, 
                        randomize.function=randomizeRegions, 
                        genome=univ, mask=g$mask, evaluate.function=numOverlaps, count.once=T, ntimes=1000)

perm_yy1.se = permTest(A=se_region, B=yy1, 
                       randomize.function=randomizeRegions, 
                       genome=univ, mask=g$mask, evaluate.function=numOverlaps, count.once=T, ntimes=1000)

plot(perm_ctcf.se)
perm_ctcf.se
plot(perm_yy1.se)
perm_yy1.se


output_ppt(plot(perm_ctcf.se), 'ctcf enrichment score with SE.pptx')
output_ppt(plot(perm_yy1.se), 'YY1 enrichment score with SE.pptx')





#-----------------------------------------------
# TAD boundaries region - YY1 is more enriched than CTCF
tad_boundaries_1 = import('P5_wholeGenome.finalboundaries.bed','BED')
tad_boundaries_2 = import('P15_wholeGenome.finalboundaries.bed','BED')

tad_boundaries = union(tad_boundaries_2, tad_boundaries_1)
data.frame(tad_boundaries) -> tad_boundaries
tad_boundaries$seqnames = as.character(tad_boundaries$seqnames)
tad_boundaries$seqnames[tad_boundaries$seqnames=='chr23'] = 'chrX'
makeGRangesFromDataFrame(tad_boundaries) -> tad_boundaries

perm_ctcf.tadb = permTest(A=tad_boundaries, B=ctcf, 
                          randomize.function=randomizeRegions, 
                          genome=univ, mask=g$mask, evaluate.function=numOverlaps, count.once=T, ntimes=100)

perm_yy1.tadb = permTest(A=tad_boundaries, B=yy1, 
                         randomize.function=randomizeRegions, 
                         genome=univ, mask=g$mask, evaluate.function=numOverlaps, count.once=T, ntimes=100)

plot(perm_ctcf.tadb)
perm_ctcf.tadb
plot(perm_yy1.tadb)
perm_yy1.tadb


output_ppt(plot(perm_ctcf.tadb), 'ctcf enrichment score with TAD boundaires.pptx')

output_ppt(plot(perm_yy1.tadb), 'YY1 enrichment score with TAD boundaires.pptx')


