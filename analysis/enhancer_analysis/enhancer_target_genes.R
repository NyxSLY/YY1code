---
title: "Find out target genes of super enhancer"
output: html_notebook
---


```{r, setup=T, message=F}
library(data.table)
library(stringr)
library(dplyr)
library(rtracklayer)
library(GenomicRanges)
library(regioneR)
library(BSgenome.Hsapiens.UCSC.hg19.masked)
g = getGenomeAndMask(BSgenome.Hsapiens.UCSC.hg19.masked)
library(ggplot2)
session::restore.session('BASE.RData')

library(commonAPI)
#source('/Users/Luyang/scripts/R/commonAPI.R')
#source('/Users/Luyang/scripts/R/annotation.R')

#source('D:\\BaiduNetdiskWorkspace\\OneDrive\\scripts\\R\\commonAPI.R')
source('D:\\BaiduNetdiskWorkspace\\OneDrive\\scripts\\R\\annotation.R')
session::restore.session('functions.RData')

#install.packages('ggbeeswarm')
library('ggbeeswarm')


library(officer)
library(rvg)

trace(grDevices::png, quote({
  if (missing(type) && missing(antialias)) {
    type <- "cairo-png"
    antialias <- "subpixel"
  }
}), print = FALSE)

getwd()
```


# read results
```{r, warning=F}
dir='total_normalized'

p5_target_region = GRanges()
p15_target_region = GRanges()

for(d in list.files(dir)){
    if(str_detect(d, "P5_")){
        files = list.files(file.path(dir, d))
        #print(files)
        file = file.path(dir, d, files[grepl('1e-3.bed', files)])
        target_region = import(file, format='BED')
        p5_target_region = c(p5_target_region, target_region)
    }
    else if(str_detect(d, "P15_")){
        files = list.files(file.path(dir, d))
        #print(files)
        file = file.path(dir, d, files[grepl('1e-3.bed', files)])
        target_region = import(file, format='BED')
        p15_target_region = c(p15_target_region, target_region)
    }
}


# filter target regions - p5 target gene only keep p5 se and common se

data.frame(p5_target_region) %>% dplyr::filter(grepl('common_se', name) | grepl('p5_se', name)) %>% makeGRangesFromDataFrame(keep.extra.columns = T) -> p5_target_region

data.frame(p15_target_region) %>% dplyr::filter(grepl('common_se', name) | grepl('p15_se', name)) %>% makeGRangesFromDataFrame(keep.extra.columns = T) -> p15_target_region

```


# distance between SE and target region
```{r}
data.frame(p5_target_region)$name -> p5.info
data.frame(p15_target_region)$name -> p15.info


c(p5.info, p15.info) %>% length()
c(p5.info, p15.info) -> info
info[!duplicated(info)] %>% length()

!duplicated(info) %>% sum()
```


# identify genes in the region
```{r}
library(EnsDb.Hsapiens.v75)
hg19=EnsDb.Hsapiens.v75
pros = promoters(hg19, upstream = 1000, downstream=1000)
pros = data.frame(pros)
pros$seqnames = as.character(pros$seqnames)
pros$seqnames = paste0('chr',pros$seqnames)
pros = makeGRangesFromDataFrame(pros, keep.extra.columns = T)

get_target_tss <- function(target_region){
    target_region.df = data.frame(target_region)
    target_region.df$id = sapply(str_split(target_region.df$name, ':'), '[', 1)
    # find out genes
    m = findOverlaps(target_region, pros)
    m = data.frame(m)
    out = tapply(m$subjectHits, m$queryHits, as.vector)
    
    genes_in_region = list()
    for(i in names(out)){
        target = data.frame(pros[out[[i]]])
        target$region_id = target_region.df$id[as.numeric(i)]
        target$target_region_chr = target_region.df$seqnames[as.numeric(i)]
        target$target_region_start = target_region.df$start[as.numeric(i)]
        target$target_region_end = target_region.df$end[as.numeric(i)]
        target$info = target_region.df$name[as.numeric(i)]
        genes_in_region[[i]] = target
    }
    do.call('rbind', genes_in_region) -> se_target
    return(se_target)
}

get_target_tss(p5_target_region) -> se_target_pro.p5
get_target_tss(p15_target_region) -> se_target_pro.p15

se_target_pro.p5[!duplicated(se_target_pro.p5$gene_id), ] -> se_target_gene.p5
se_target_pro.p15[!duplicated(se_target_pro.p15$gene_id), ] -> se_target_gene.p15

se_target_gene.p5 -> se_target_gene.p5.save
se_target_gene.p15 -> se_target_gene.p15.save
```


# how many genes in the target region per SE
```{r}
new = se_target_gene.p5

count_ = function(new){
  new$pair = paste(new$region_id, new$target_region_chr, new$target_region_start, new$target_region_end, sep='_')
  tapply(new$gene_id,new$pair, function(X) length(unique(X)))
}

count_both = function(new, new1){
  new = rbind(new, new1)
  new$pair = paste(new$region_id, new$target_region_chr, new$target_region_start, new$target_region_end, sep='_')
  tapply(new$gene_id,new$pair, function(X) length(unique(X)))
}

count_(se_target_gene.p5) -> gene_count_per_se.p5
count_(se_target_gene.p15) -> gene_count_per_se.p15

count_both(se_target_gene.p5, se_target_gene.p15) -> gene_count_per_se

barplot(table(gene_count_per_se.p5), xlim=c(1,21))
barplot(table(gene_count_per_se.p15), xlim=c(1,21))

barplot(table(gene_count_per_se))

output_ppt(barplot(table(gene_count_per_se)), 'number of genes in the target region per SE.pptx')
getwd()
```


# gene expression
```{r}
library(data.table)
fpkm = fread('tpm.csv')
#fpkm = fread('fpkm_remove_batch_effect.csv')
keep= rowSums(fpkm[,2:7]>1)>= 3
fpkm = fpkm[keep, ]

fpkm$p5_mean = (fpkm$rep1_young+fpkm$rep2_young+fpkm$rep3_young)/3
fpkm$p15_mean = (fpkm$rep1_old+fpkm$rep2_old+fpkm$rep3_old)/3

se_target_gene.p5 <- se_target_gene.p5.save
se_target_gene.p15 <- se_target_gene.p15.save

merge(se_target_gene.p5, fpkm, by='gene_id') -> se_target_gene.p5
merge(se_target_gene.p15, fpkm, by='gene_id') -> se_target_gene.p15

#hist(log2(se_target_gene.p5$p5_mean))
#hist(log2(se_target_gene.p15$p15_mean))

#hist(log2(fpkm$p5_mean))
#hist(log2(fpkm$p15_mean))

#se_target_gene.p15[!se_target_gene.p15$gene_id %in% se_target_gene.p5$gene_id,] -> new
# boxplot(log2(new$p5_mean))
# 
# boxplot(log2(se_target_gene.p5$p5_mean),
#         log2(se_target_gene.p15$p15_mean),
#         log2(fpkm$p5_mean),
#         log2(fpkm$p15_mean), outline=F,
#         names=c('P5_se_gene',
#                 'P15_se_gene',
#                 'P5_all_gene',
#                 'P15_all_gene'))
# 
# boxplot(log2(se_target_gene.p5$rep1_young),
#         log2(se_target_gene.p5$rep2_young),
#         log2(se_target_gene.p5$rep3_young),
#         log2(se_target_gene.p5$rep1_old),
#         log2(se_target_gene.p5$rep2_old),
#         log2(se_target_gene.p5$rep3_old),
#         log2(new$p15_mean),
#         log2(fpkm$p5_mean),
#         log2(fpkm$p15_mean), outline=F)


test <- function(fpkm, se_target_gene.p5, se_target_gene.p15){
  keep= rowSums(fpkm[,2:7]>1)>= 3
  fpkm = fpkm[keep, ]
  
  fpkm$p5_mean = (fpkm$rep1_young+fpkm$rep2_young+fpkm$rep3_young)/3
  fpkm$p15_mean = (fpkm$rep1_old+fpkm$rep2_old+fpkm$rep3_old)/3
  
  merge(se_target_gene.p5, fpkm, by='gene_id') -> se_target_gene.p5
  merge(se_target_gene.p15, fpkm, by='gene_id') -> se_target_gene.p15
  
  #hist(log2(se_target_gene.p5$p5_mean))
  #hist(log2(se_target_gene.p15$p15_mean))
  
  #hist(log2(fpkm$p5_mean))
  #hist(log2(fpkm$p15_mean))
  
  boxplot(log2(se_target_gene.p5$p5_mean),
          log2(se_target_gene.p15$p15_mean),
          log2(fpkm$p5_mean),
          log2(fpkm$p15_mean), outline=F,
          names=c('P5_se_gene',
                  'P15_se_gene',
                  'P5_all_gene',
                  'P15_all_gene'))
}


test(fread('tpm.csv'), se_target_gene.p5.save, se_target_gene.p15.save)
test(fread('fpkm.csv'), se_target_gene.p5.save, se_target_gene.p15.save)


pdf(NULL)
dev.control(displaylist="enable")
boxplot(log2(se_target_gene.p5$p5_mean),
        log2(se_target_gene.p15$p15_mean),
        log2(fpkm$p5_mean),
        log2(fpkm$p15_mean), outline=F,
        names=c('P5_se_gene',
                'P15_se_gene',
                'P5_all_gene',
                'P15_all_gene'),
        ylim=c(-5,15))

p.expression <- recordPlot()
invisible(dev.off())

output_ppt(p.expression, 'target gene expression mean.pptx')

library(ggplot2)

P5_se_gene <- data.frame(expression = log2(se_target_gene.p5$p5_mean), condition = 'P5_se_gene')
P15_se_gene <- data.frame(expression = log2(se_target_gene.p15$p15_mean), condition = 'P15_se_gene')
P5_all_gene <- data.frame(expression = log2(fpkm$p5_mean), condition = 'P5_all_gene')
P15_all_gene <- data.frame(expression = log2(fpkm$p15_mean), condition = 'P15_all_gene')

# Merge the four data frames into one long-format data frame
long_data <- bind_rows(P5_se_gene, P15_se_gene, P5_all_gene, P15_all_gene)
saveRDS(long_data, 'long_data.rds')

long_data = readRDS('long_data.rds')
p = ggplot(long_data, aes(x = condition, y = expression, fill = condition)) +
  geom_jitter(width = 0.2, aes(color = condition), size = 1.5, alpha = 0.1) +  # Add jitter with matching color
  # geom_violin(trim=T, adjust=1, width=0.8, color='black', alpha=0.7) +  # Hide the outliers in the boxplot
  geom_boxplot(outlier.shape = NA) +
  scale_colour_brewer(palette = "Set2") +
  scale_fill_brewer(palette = "Set2") +
  ylim(-5, 15) +
  labs(x = 'Condition', y = 'Log2(Expression)') +
  theme_classic()

output_ppt(p, 'boxplot_jitter-target gene expression mean.pptx')



wilcox.test(c(log2(se_target_gene.p5$p5_mean),log2(se_target_gene.p15$p15_mean)),
            c(log2(fpkm$p5_mean), log2(fpkm$p15_mean)),
            alternative = "two.sided")

wilcox.test(log2(se_target_gene.p5$p5_mean), log2(fpkm$p5_mean),
            alternative = "two.sided")

wilcox.test(log2(se_target_gene.p15$p15_mean),log2(fpkm$p15_mean),
            alternative = "two.sided")


# common genes and specific genes

se_target_gene.p5[!se_target_gene.p5$gene_id %in% se_target_gene.p15$gene_id,] -> se_target_gene.p5.sp  # n = 281
se_target_gene.p15[!se_target_gene.p15$gene_id %in% se_target_gene.p5$gene_id,] -> se_target_gene.p15.sp # n = 279
se_target_gene.p5[se_target_gene.p5$gene_id %in% se_target_gene.p15$gene_id,] -> se_target_gene.p5.common # n = 477

# new 20231123

P5_se_gene <- data.frame(expression = log2(se_target_gene.p5$p5_mean), condition = 'P5_se_gene')
P15_se_gene <- data.frame(expression = log2(se_target_gene.p15$p15_mean), condition = 'P15_se_gene')
P5_all_gene <- data.frame(expression = log2(fpkm$p5_mean), condition = 'P5_all_gene')
P15_all_gene <- data.frame(expression = log2(fpkm$p15_mean), condition = 'P15_all_gene')
P5_common_gene <- data.frame(expression = log2(se_target_gene.p5.common$p5_mean), condition = 'P5_common_gene')
P15_common_gene <- data.frame(expression = log2(se_target_gene.p5.common$p15_mean), condition = 'P15_common_gene')


# Merge the four data frames into one long-format data frame
plot.dat <- bind_rows(P5_se_gene, P15_se_gene, P5_all_gene, P15_all_gene,
                       P5_common_gene, P15_common_gene)

plot.dat$condition = factor(plot.dat$condition, levels = c('P5_common_gene', 'P5_se_gene', 'P5_all_gene',
                                                             'P15_common_gene', 'P15_se_gene', 'P15_all_gene'))

colnames(plot.dat) = c('value', 'key')
for(k in unique(plot.dat$key)){
  boxplot_filter(plot.dat[plot.dat$key==k,]) -> plot.dat[plot.dat$key==k,]
}

saveRDS(plot.dat, 'FigS2G.rds')

p = ggplot(plot.dat, aes(x = key, y = value, fill = key)) +
  geom_jitter(width = 0.2, aes(color = key), size = 1.5, alpha = 0.1) +  # Add jitter with matching color
  #geom_violin(trim=T, adjust=2, width=0.8, color='black', alpha=0.7) +  # Hide the outliers in the boxplot
  geom_boxplot(outlier.shape = NA, alpha=0.5) +
  scale_colour_brewer(palette = "Set2") +
  scale_fill_brewer(palette = "Set2") +
  theme_classic()


pdf(NULL)
dev.control(displaylist="enable")
boxplot(se_target_gene.p5.common$p5_mean %>% log2,
        se_target_gene.p5.common$p15_mean %>% log2,
  log2(se_target_gene.p5.sp$p5_mean),
        log2(se_target_gene.p15.sp$p15_mean),
        log2(fpkm$p5_mean),
        log2(fpkm$p15_mean), outline=F,
        names=c('common genes P5',
                'common genes P15',
          'P5_se_gene',
                'P15_se_gene',
                'P5_all_gene',
                'P15_all_gene'))

p.expression <- recordPlot()
invisible(dev.off())

output_ppt(p.expression, 'target gene expression mean specific genes.pptx')



wilcox.test(log2(se_target_gene.p5.common$p5_mean), log2(se_target_gene.p5.sp$p5_mean))
wilcox.test(log2(se_target_gene.p5.common$p5_mean), log2(fpkm$p5_mean))


wilcox.test(log2(se_target_gene.p5.common$p15_mean), log2(se_target_gene.p15.sp$p15_mean))
wilcox.test(log2(se_target_gene.p5.common$p15_mean), log2(fpkm$p15_mean))



# # violin plot
# plot.dat = rbind(
#   data.frame(value = se_target_gene.p5.common$p5_mean,name = 'common_p5'),
#   data.frame(value = se_target_gene.p5.common$p5_mean,name = 'common_p15'),
#   data.frame(value = se_target_gene.p5.sp$p5_mean,name = 'sp_p5'),
#   data.frame(value = se_target_gene.p15.sp$p15_mean,name = 'sp_p15'),
#   data.frame(value = fpkm$p5_mean,name = 'all_p5'),
#   data.frame(value = fpkm$p15_mean,name = 'all_p15')
#   )
# 
# ggplot(plot.dat, aes(x=name, y=log2(value), fill=name)) + geom_violin(position=position_dodge(width=0.9))+
#   geom_quasirandom()
# 




getGeneAnnoEn(se_target_gene.p5$gene_id) -> p5.anno
p5.anno[order(p5.anno$entrez_id),]-> p5.anno
p5.anno[!duplicated(p5.anno$gene_id),] -> p5.anno
merge(se_target_gene.p5, p5.anno, by='gene_id') -> se_target_gene.p5

getGeneAnnoEn(se_target_gene.p15$gene_id) -> p15.anno
p15.anno[order(p15.anno$entrez_id),]-> p15.anno
p15.anno[!duplicated(p15.anno$gene_id),] -> p15.anno
merge(se_target_gene.p15, p15.anno, by='gene_id') -> se_target_gene.p15

write.csv(se_target_gene.p5, 'p5_se_target_gene.csv', row.names = F)
write.csv(se_target_gene.p15, 'p15_se_target_gene.csv', row.names = F)


```


# pathway and go-term
```{r}
library(data.table)
library(rtracklayer)
library(GenomicRanges)
library(ChIPseeker)
library(clusterProfiler)
library(ReactomePA)
library(ggplot2)

# db='org.Mm.eg.db'
db='org.Hs.eg.db'
species = 'human'

enrichment_analysis <- function(entrez_id, db, species){
  go_res = enrichGO(
    gene = entrez_id,
    OrgDb=db,
    keyType = "ENTREZID",
    ont = "BP",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    readable = TRUE)
  
  
  pw_res = enrichPathway(gene=entrez_id, pvalueCutoff = 0.05, pAdjustMethod = "BH", organism = species)
  
  pw.plot <- dotplot(pw_res) + theme_classic()
  go.plot <- dotplot(go_res) + theme_classic()
  
  pw.plot
  go.plot
  
  return(list(pw=pw_res, go=go_res))
}


p5_pw = enrichment_analysis(se_target_gene.p5$entrez_id, db, species)
p15_pw = enrichment_analysis(se_target_gene.p15$entrez_id, db, species)

name='P5'

write.csv(data.frame(data.frame(p5_pw$pw)), paste0(name, ' reactome pathway.csv'))
write.csv(data.frame(data.frame(p5_pw$go)), paste0(name, ' GO enrichment pathway.csv'))

name='P15'

write.csv(data.frame(data.frame(p15_pw$pw)), paste0(name, ' reactome pathway.csv'))
write.csv(data.frame(data.frame(p15_pw$go)), paste0(name, ' GO enrichment pathway.csv'))


se_target_gene.p5 %>% dplyr::filter(grepl('common_se', region_id)) -> se_target_gene.common

common_pw = enrichment_analysis(se_target_gene.common$entrez_id, db, species)


name='common se target gene'

write.csv(data.frame(data.frame(common_pw$pw)), paste0(name, ' reactome pathway.csv'))
write.csv(data.frame(data.frame(common_pw$go)), paste0(name, ' GO enrichment pathway.csv'))
write.csv(se_target_gene.common, 'common_se_target_gene.csv', row.names = F)
```


# YY1 ChIP-seq, show target gene and SE both enriched in YY1

```{r}
# import YY1 peak
yy1_p15 = readbed('P15_peaks.narrowPeak')
yy1_p5 = readbed('P5_peaks.narrowPeak')
yy1 =reduce(union(yy1_p15, yy1_p5))

# import CTCF peak
ctcf_p5='P5_CTCF_macs2_peaks.narrowPeak'
ctcf_p15='P15_CTCF_macs2_peaks.narrowPeak'
ctcf_p5 = import(ctcf_p5)
ctcf_p15 = import(ctcf_p15)
ctcf = union(ctcf_p5, ctcf_p15)

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

```



# enhancer region - YY1 is more enriched than CTCF
```{r}
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


```



# TAD boundaries region - YY1 is more enriched than CTCF
```{r}
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


```


# bu - enhancer region and pol2 peak region
```{r}
library(regioneR)
pol2_p5 = readbed('C:\\Users\\Luyang_pc\\Documents\\BaiduNetdiskWorkspace\\OneDrive\\work\\3dGenome\\ChIP-seq_signal_tracks/Pol2_Y_peak.bed')
pol2_p15 = readbed('C:\\Users\\Luyang_pc\\Documents\\BaiduNetdiskWorkspace\\OneDrive\\work\\3dGenome\\ChIP-seq_signal_tracks/Pol2_O_peak.bed')
pol2 =reduce(union(pol2_p5, pol2_p15))

perm_ctcf.pol2 = permTest(A=ctcf, B=pol2, 
               randomize.function=randomizeRegions, 
               genome=univ, mask=g$mask, evaluate.function=numOverlaps, count.once=T, ntimes=1000)

perm_yy1.pol2 = permTest(A=yy1, B=pol2, 
               randomize.function=randomizeRegions, 
               genome=univ, mask=g$mask, evaluate.function=numOverlaps, count.once=T, ntimes=1000)

plot(perm_ctcf.pol2)


typical_enhancer = readbed('typical_enhancer.bed')

perm_ctcf.typical_enhancer = permTest(A=ctcf, B=typical_enhancer, 
               randomize.function=randomizeRegions, 
               genome=univ, mask=g$mask, evaluate.function=numOverlaps, count.once=T, ntimes=1000)


perm_yy1.typical_enhancer = permTest(A=yy1, B=typical_enhancer, 
               randomize.function=randomizeRegions, 
               genome=univ, mask=g$mask, evaluate.function=numOverlaps, count.once=T, ntimes=1000)
```



# YY1 KD RNA-seq
```{r}
yy1.fpkm = read.csv('YY1_KD_DEG/fpkm.csv')

yy1.fpkm$YY1_KD_mean = rowMeans(yy1.fpkm[,c('YY1_KD_rep1', 'YY1_KD_rep2', 'YY1_KD_rep3')])
yy1.fpkm$YY1_NT_mean = rowMeans(yy1.fpkm[,c('YY1_NT_rep1', 'YY1_NT_rep2', 'YY1_NT_rep3')])

se_target_gene.p5 = merge(se_target_gene.p5, yy1.fpkm, by='gene_id')
se_target_gene.p15 = merge(se_target_gene.p15, yy1.fpkm, by='gene_id')


se_target_gene.p5$yy1_kd_fc = log2((se_target_gene.p5$YY1_KD_mean+0.1)/(se_target_gene.p5$YY1_NT_mean+0.1))
se_target_gene.p15$yy1_kd_fc = log2((se_target_gene.p15$YY1_KD_mean+0.1)/(se_target_gene.p15$YY1_NT_mean+0.1))

yy1.fpkm$log2_fc = log2((yy1.fpkm$YY1_KD_mean+0.1)/(yy1.fpkm$YY1_NT_mean+0.1))
yy1.fpkm.filter = yy1.fpkm %>% dplyr::filter(gene_id %in% fpkm$gene_id)
yy1.fpkm.filter = yy1.fpkm.filter %>% dplyr::filter(!is.na(YY1_KD_mean))


quantile(abs(se_target_gene.p5$yy1_kd_fc), probs=seq(0,0.9, 0.25)) -> a
quantile(abs(se_target_gene.p15$yy1_kd_fc), probs=seq(0,0.9, 0.25)) -> b
quantile(abs(yy1.fpkm.filter$log2_fc), probs=seq(0,0.9, 0.25)) -> c

barplot(t(data.frame(a, b, c)%>% as.matrix()), beside=T)


boxplot(abs(se_target_gene.p5$yy1_kd_fc),
        abs(yy1.fpkm.filter$log2_fc), outline=F)

boxplot(abs(se_target_gene.p15$yy1_kd_fc),
        abs(yy1.fpkm.filter$log2_fc), outline=F)

t.test(abs(se_target_gene.p5$yy1_kd_fc), 
       abs(yy1.fpkm.filter$log2_fc))

t.test(abs(se_target_gene.p15$yy1_kd_fc), 
       abs(yy1.fpkm.filter$log2_fc))

# The expression of SE target genes were disregulated after YY1 KD. The gene
# expression changes are bigger than all expressed genes.


# combine target genes

rbind(se_target_gene.p5, se_target_gene.p15) -> combined.target_gene
combined.target_gene[!duplicated(combined.target_gene$gene_id),] -> combined.target_gene

quantile(abs(combined.target_gene$yy1_kd_fc), probs=seq(0,0.9, 0.25)) -> a
quantile(abs(yy1.fpkm.filter$log2_fc), probs=seq(0,0.9, 0.25)) -> b

barplot(t(data.frame(a, b)%>% as.matrix()), beside=T, ylim=c(0,0.8))

boxplot(abs(combined.target_gene$yy1_kd_fc),
        abs(yy1.fpkm.filter$log2_fc), outline=F)

t.test(abs(combined.target_gene$yy1_kd_fc), 
       abs(yy1.fpkm.filter$log2_fc))


#---------------------save plot----------
pdf(NULL)
dev.control(displaylist="enable")
barplot(t(data.frame(a, b)%>% as.matrix()), beside=T, ylim=c(0,0.8))
p.expression <- recordPlot()
invisible(dev.off())
output_ppt(p.expression, 'SE target genes are disregulated after YY1 knockdown - barplot.pptx')

pdf(NULL)
dev.control(displaylist="enable")
boxplot(abs(combined.target_gene$yy1_kd_fc),
        abs(yy1.fpkm.filter$log2_fc), outline=F)
p.expression <- recordPlot()
invisible(dev.off())
output_ppt(p.expression, 'SE target genes are disregulated after YY1 knockdown - boxplot.pptx')


# new - 2023-11-22
dat.t1 = data.frame(key='SE_target', value=abs(combined.target_gene$yy1_kd_fc))
dat.t2 = data.frame(key='all', value=abs(yy1.fpkm.filter$log2_fc))
plot.dat = rbind(dat.t1, dat.t2)


#plot.dat = gather(all_pro$dat)
plot.dat$key = factor(plot.dat$key, levels=c('SE_target', 'all'))
#plot.dat$value = log1p(plot.dat$value)
boxplot_filter(plot.dat[plot.dat$key=='SE_target',]) -> plot.dat[plot.dat$key=='SE_target',]
boxplot_filter(plot.dat[plot.dat$key=='all',]) -> plot.dat[plot.dat$key=='all',]

saveRDS(plot.dat, 'SE target genes are disregulated after YY1 knockdown.rds')

p = ggplot(plot.dat, aes(x = key, y = value, fill = key)) +
  geom_jitter(width = 0.2, aes(color = key), size = 1.5, alpha = 0.1) +  # Add jitter with matching color
  #geom_violin(trim=T, adjust=2, width=0.8, color='black', alpha=0.7) +  # Hide the outliers in the boxplot
  geom_boxplot(outlier.shape = NA, alpha=0.5) +
  scale_colour_brewer(palette = "Set2") +
  scale_fill_brewer(palette = "Set2") +
  theme_classic()
p







#---------------------


# is these genes enriched in DEGs of YY1KD?
se_target_gene.p5 %>% dplyr::filter(tx_biotype=='protein_coding') -> se_target_gene.p5.pc
se_target_gene.p15 %>% dplyr::filter(tx_biotype=='protein_coding') -> se_target_gene.p15.pc


YY1_sig_genes = 'YY1_KD_DEG/DESEQ2_sig_genes.csv'
YY1_sig_genes = fread(YY1_sig_genes)



sum(se_target_gene.p5.pc$gene_id %in% YY1_sig_genes$gene_id) / nrow(se_target_gene.p5.pc)
sum(se_target_gene.p15.pc$gene_id %in% YY1_sig_genes$gene_id) / nrow(se_target_gene.p15.pc)


# fisher's test or geometric test
a.fish = sum(se_target_gene.p5.pc$gene_id %in% YY1_sig_genes$gene_id)
b.fish = nrow(se_target_gene.p5.pc)

c.fish = nrow(YY1_sig_genes[YY1_sig_genes$biotype=='protein_coding'])

YY1_KD_allgene = fread('YY1_KD_DEG/DESEQ2_total_results.csv')
YY1_KD_allgene %>% dplyr::filter(!is.na(padj), biotype=='protein_coding') -> YY1_KD_allgene
d.fish = nrow(YY1_KD_allgene)

fisher.test.res = fisher.test(matrix(c(a.fish, b.fish, c.fish, d.fish),2,2))
fisher.test.res$p.value

####################################
# Yes, it is significantly enriched
####################################

# seperate YY1 up and YY1 down genes
YY1_sig_genes[YY1_sig_genes$biotype=='protein_coding'] %>% dplyr::filter(log2FoldChange<0) -> YY1_sig_genes.dn
YY1_sig_genes[YY1_sig_genes$biotype=='protein_coding'] %>% dplyr::filter(log2FoldChange>0) -> YY1_sig_genes.up

merge(se_target_gene.p5.pc, YY1_sig_genes, by='gene_id') -> se_target_gene.p5.pc
se_target_gene.p5.pc[se_target_gene.p5.pc$log2FoldChange > 0,] -> se_target_gene.p5.pc.up_YY1KD
se_target_gene.p5.pc[se_target_gene.p5.pc$log2FoldChange < 0,] -> se_target_gene.p5.pc.dn_YY1KD


barplot(c(nrow(se_target_gene.p5.pc.up_YY1KD),
        nrow(se_target_gene.p5.pc.dn_YY1KD)))

barplot(c(nrow(YY1_sig_genes.up),
        nrow(YY1_sig_genes.dn)))

fisher.test(matrix(c(nrow(se_target_gene.p5.pc.up_YY1KD),
                     nrow(se_target_gene.p5.pc.dn_YY1KD),
                     nrow(YY1_sig_genes.up),
                     nrow(YY1_sig_genes.dn)),2,2))


# Try with P15
merge(se_target_gene.p15.pc, YY1_sig_genes, by='gene_id') -> se_target_gene.p15.pc
se_target_gene.p15.pc[se_target_gene.p15.pc$log2FoldChange > 0,] -> se_target_gene.p15.pc.up_YY1KD
se_target_gene.p15.pc[se_target_gene.p15.pc$log2FoldChange < 0,] -> se_target_gene.p15.pc.dn_YY1KD

barplot(c(nrow(se_target_gene.p15.pc.up_YY1KD),
        nrow(se_target_gene.p15.pc.dn_YY1KD)))

barplot(c(nrow(YY1_sig_genes.up),
        nrow(YY1_sig_genes.dn)))

fisher.test(matrix(c(nrow(se_target_gene.p15.pc.up_YY1KD),
                     nrow(se_target_gene.p15.pc.dn_YY1KD),
                     nrow(YY1_sig_genes.up),
                     nrow(YY1_sig_genes.dn)),2,2))


# This is actually make sence? - separate common and specific ones. see if still can get the diff

```

# YY1 KD scatter plot
```{r}
yy1.deseq2 = read.csv('YY1_KD_DEG/DESEQ2_total_results.csv')
yy1.deseq2 = yy1.deseq2[!is.na(yy1.deseq2$padj),]


plot(log10(yy1.deseq2$baseMean), yy1.deseq2$log2FoldChange)
plot(yy1.deseq2$log2FoldChange, -log2(yy1.deseq2$padj))


dat.scater = yy1.deseq2 %>% dplyr::select(gene_id, baseMean,log2FoldChange, padj)

dat.scater[dat.scater$padj==0,]$padj = 1e-50



dat.scater %>% mutate(sig=cut(padj, breaks=c(0,0.05,1), labels=c('sig', 'notsig'))) -> dat.scater
dat.scater %>% mutate(direction=cut(log2FoldChange, breaks=c(-Inf, 0, Inf), labels=c('DN', 'UP'))) -> dat.scater
dat.scater %>% mutate(group=paste0(dat.scater$sig, dat.scater$direction))-> dat.scater


ggplot(dat.scater, aes(x=log2FoldChange, y=-log2(padj), color=group)) + geom_point(size=0.5) +
  theme_classic() + 
  scale_color_manual(values=c('gray', 'gray', 'blue','red'))

dat.scater$sig %>% table()

scatterplot.g = ggplot(dat.scater, aes(x=log10(baseMean), y=log2FoldChange, color=group)) + geom_point(size=0.1) +
  theme_classic() + 
  scale_color_manual(values=c('gray', 'gray', 'blue','red'))

output_ppt(scatterplot.g, 'YY1KD RNA-seq scatter plot.pptx')
scatterplot.g
```


# Fig 3B and C - YY1 KD of all SE target genes
```{r}
rbind(p5pair.genes.fpkm, p15pair.genes.fpkm) %>%
  dplyr::filter(!duplicated(pair)) -> all_se_pair.fpkm

all_se_pair.fpkm.uniq = all_se_pair.fpkm[!duplicated(all_se_pair.fpkm$gene_id),]

deseq2_scatter_plot_yy1kd <- function(yy1.deseq2){
  yy1.deseq2 = yy1.deseq2[!is.na(yy1.deseq2$yy1kd_padj),]
  dat.scater = yy1.deseq2 %>% dplyr::select(gene_id, yy1kd_baseMean,yy1kd_log2FoldChange, yy1kd_padj)
  if(sum(dat.scater$yy1kd_padj==0)!=0){
    dat.scater[dat.scater$yy1kd_padj==0,]$yy1kd_padj = 1e-50
  }
  dat.scater %>% mutate(sig=cut(yy1kd_padj, breaks=c(0,0.05,1), labels=c('sig', 'notsig'))) -> dat.scater
  dat.scater %>% mutate(direction=cut(yy1kd_log2FoldChange, breaks=c(-Inf, 0, Inf), labels=c('DN', 'UP'))) -> dat.scater
  dat.scater %>% mutate(group=paste0(dat.scater$sig, dat.scater$direction))-> dat.scater
  dat.scater %>% mutate(point_size = cut(yy1kd_padj, breaks=c(0,0.05,1), labels=c(0.5, 0.5))) -> dat.scater
  dat.scater %>% mutate(point_size=as.numeric(as.character(point_size))) -> dat.scater
  
  main <- ggplot(dat.scater, aes(x=log10(yy1kd_baseMean), y=yy1kd_log2FoldChange, color=group)) + 
    geom_point(aes(size=point_size)) +
    theme_classic() + 
    theme(legend.position = "none") +
    scale_color_manual(values=c('gray', 'gray', 'blue','red')) + 
    scale_size(range = c(0.5, 0.5))
  
  return(list(plot=main,dat=dat.scater))
}

deseq2_scatter_plot_yy1kd(all_se_pair.fpkm) -> se_pair.plot

se_pair.plot$plot + ylim(-5,5) -> se_pair.yy1kd.plot

se_pair.yy1kd.plot

se_pair.plot$dat$sig %>% table()

output_ppt_m(se_pair.yy1kd.plot, 'YY1 KD scatter plot of all SE targeted genes.pptx')


fisher.test(matrix(c(sum(se_pair.plot$dat$sig=='sig'),
                     nrow(se_pair.plot$dat),
                     sum(dat.scater$sig=='sig'),
                     length(dat.scater$sig)),2,2)) -> fisher.pvalue

fisher.pvalue$p.value
```


# YY1 peak decrease => SE-P pair broken => (gene expression change) --> no gene expression
```{r}
# regions with YY1 decrease
yy1_p5[!overlapsAny(yy1_p5, yy1_p15)] -> yy1.down

yy1_p15[!overlapsAny(yy1_p15, yy1_p5)] -> yy1.up

#TODO
# YY1 peak down != gene down-regulation: need to first find out genes that with YY1 peak down, their expression also down. 
# then check the gene list.


# # down-regulated genes with YY1 down
# YY1_sig_genes.dn
# 
# pros.df = data.frame(pros)
# 
# merge(YY1_sig_genes.dn, pros.df, by='gene_id') -> YY1_sig_genes.dn.Pro
# makeGRangesFromDataFrame(YY1_sig_genes.dn.Pro, keep.extra.columns = T) -> YY1_sig_genes.dn.Pro.gr
# 
# YY1_sig_genes.dn.Pro.gr[overlapsAny(YY1_sig_genes.dn.Pro.gr, yy1.down),] -> YY1_dn_genes.withYY1PeakDN_at_promoter
# YY1_dn_genes.withYY1PeakDN_at_promoter.genes = YY1_dn_genes.withYY1PeakDN_at_promoter[!duplicated(YY1_dn_genes.withYY1PeakDN_at_promoter$gene_id),]
# 
# # link with genes that dn reulated in YY1 KD & YY1 peak down at promoter
# merge(p5.pair.withYY1down.genes.pc, YY1_dn_genes.withYY1PeakDN_at_promoter.genes, by='gene_id') -> p5.pair.withYY1down.genes.pc
# merge(p5.pair.withYY1down.genes.pc, yy1.fpkm, by='gene_id') -> p5.pair.withYY1down.genes.pc
# merge(p5.pair.withYY1down.genes.pc, fpkm, by='gene_id') -> p5.pair.withYY1down.genes.pc


# make pair info
se_target_pro.p5$pair = paste(se_target_pro.p5$region_id, se_target_pro.p5$gene_id, sep='_')
se_target_pro.p15$pair = paste(se_target_pro.p15$region_id, se_target_pro.p15$gene_id, sep='_')

p5.pair = se_target_pro.p5[!duplicated(se_target_pro.p5$pair),]
p15.pair = se_target_pro.p15[!duplicated(se_target_pro.p15$pair),]

# only care about P5 now
# import SE region
se = import('super_enhancer.bed')
se.df = data.frame(se)
names(se.df)[1:3] = c('se_chr', 'se_start', 'se_end')

merge(p5.pair, se.df, by.x='region_id', by.y='name') -> p5.pair
p5.pair %>% dplyr::select(seqnames, start, end, gene_id, tx_biotype, region_id, se_chr, se_start, se_end, 
                          target_region_chr, target_region_start, target_region_end,pair, info) -> p5.pair

merge(p15.pair, se.df, by.x='region_id', by.y='name') -> p15.pair
p15.pair %>% dplyr::select(seqnames, start, end, gene_id, tx_biotype, region_id, se_chr, se_start, se_end, 
                          target_region_chr, target_region_start, target_region_end,pair, info) -> p15.pair

write.csv(p5.pair, 'P5_SE_and_its_target.csv',quote = F, row.names = F)
write.csv(p15.pair, 'P15_SE_and_its_target.csv',quote = F, row.names = F)
```
# continue with a separate file:
# /Volumes/Luyang1/onedriveMe/OneDrive/work/3dGenome/Find_enhancer_target/check_se_pro_pair_hic_reads_changes.R

# conclusion:
# p5.pair.withYY1down - higher hic contact at P5 than P15
# p15.pair.withYY1down - higher hic contact at P15 than P5








# try to define new lost YY1 with aging regions
```{r}




```




# expression related results is in gene_expression.Rmd


# Changes during aging

```{r}
YY1





```




















# GO-term of these genes
```{r}

getGeneAnnoEn(p5.pair.withYY1down.genes.pc$gene_id) -> p5.pair.withYY1down.genes.pc.info


p5.pair.withYY1down.genes.pc_pw = enrichment_analysis(p5.pair.withYY1down.genes.pc.info$entrez_id, db, species)


name='P5_pair_withYY1_down_protein_coding_genes'

write.csv(data.frame(data.frame(p5.pair.withYY1down.genes.pc_pw$pw)), paste0(name, ' reactome pathway.csv'))
write.csv(data.frame(data.frame(p5.pair.withYY1down.genes.pc_pw$go)), paste0(name, ' GO enrichment pathway.csv'))


getGeneAnnoEn(p15.pair.withYY1up.genes.pc$gene_id) -> p15.pair.withYY1up.genes.pc.info


p15.pair.withYY1up.genes.pc_pw = enrichment_analysis(p15.pair.withYY1up.genes.pc.info$entrez_id, db, species)


name='P15_pair_withYY1_up_protein_coding_genes'

write.csv(data.frame(data.frame(p15.pair.withYY1up.genes.pc_pw$pw)), paste0(name, ' reactome pathway.csv'))
write.csv(data.frame(data.frame(p15.pair.withYY1up.genes.pc_pw$go)), paste0(name, ' GO enrichment pathway.csv'))
```





# charecterize active super-enhacer: more CBP and Pol2 signal

```{r}

```

# TFs - important TFs
```{r}
# Key Transcription Factors in the Differentiation of Mesenchymal Stem Cells
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5010472/#:~:text=The%20major%20transcription%20factors%20that%20have%20been%20reported%20to%20have,during%20the%20differentiation%20of%20MSCs.

# Transcriptional Control of Mesenchymal Stem Cell Differentiation
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3083289/

# transcription-factor list
# https://www.rndsystems.com/research-area/mesenchymal-stem-cell-transcription-factors

# ***important***
# Osteogenesis depends on commissioning of a network of stem cell transcription factors that act as repressors of adipogenesis
# https://www.nature.com/articles/s41588-019-0359-1?proof=t

# MSC to NPC (neural progenitor-like cells)
#Significant transcriptomic changes are associated with differentiation of bone marrow-derived mesenchymal stem cells into neural progenitor-like cells in the presence of bFGF and EGF

# https://cellandbioscience.biomedcentral.com/articles/10.1186/s13578-020-00487-z
```

# Aging realted changes

## gene expression

What I am looking for is actually changes of SE-PRO pair


```{r}
se_target_pro.p5$pair = paste(se_target_pro.p5$region_id, se_target_pro.p5$gene_id, sep='_')
se_target_pro.p15$pair = paste(se_target_pro.p15$region_id, se_target_pro.p15$gene_id, sep='_')

p5.pair = se_target_pro.p5[!duplicated(se_target_pro.p5$pair),]
p15.pair = se_target_pro.p15[!duplicated(se_target_pro.p15$pair),]

# gene based list
p5.pair[, c('gene_id', 'pair')] -> p5.tmp
p15.pair[, c('gene_id', 'pair')] -> p15.tmp
gene_list = c(p5.tmp$gene_id, p15.tmp$gene_id) %>% unique()




# gene_based.df = data.frame(gene_id=gene_list)
# merge(gene_based.df, p5.tmp, by='gene_id', all.x=T) -> gene_based.df
# merge(gene_based.df, p15.tmp, by='gene_id', all.x=T) -> gene_based.df
# 
# # filter the list
# 
# gene_based.df %>% arrange(gene_id, pair.x) %>% dplyr::select(gene_id, pair.x) -> t1
# gene_based.df %>% arrange(gene_id, pair.y) %>% dplyr::select(gene_id, pair.y) -> t2
# gene_based.df <- cbind(t1, t2$pair.y)
# gene_based.df <- unique(gene_based.df)
# names(gene_based.df) <- c('gene_id', 'p5_pair', 'p15_pair')


```



```
```{r}
# YY1 ChIP-seq changes

# SE that swtich it's target region - and receiprocal - a single gene switch it's SE

# Hi-C data changes
# 1) Age-specific SE - and their target genes' expression
# 2) common SE but has significant difference in interactions - and their target genes' expression
# 3) And gene

```





# correlation between YY1 KD and aging - hypothesis
```{r}
se_target_gene.p5.exp = se_target_gene.p5 %>% dplyr::select(gene_id, biotype,
                                                         p5_mean, p15_mean, 
                                                         YY1_KD.FPKM, NT.FPKM, yy1_kd_fc)

se_target_gene.p15.exp = se_target_gene.p15 %>% dplyr::select(gene_id, biotype,
                                                         p5_mean, p15_mean, 
                                                         YY1_KD.FPKM, NT.FPKM, yy1_kd_fc)


se_target_gene.p5.exp$age_fc = log2((se_target_gene.p5.exp$p15_mean+0.1)/(se_target_gene.p5.exp$p5_mean+0.1))
se_target_gene.p15.exp$age_fc = log2((se_target_gene.p15.exp$p15_mean+0.1)/(se_target_gene.p15.exp$p5_mean+0.1))



plot(se_target_gene.p5.exp$yy1_kd_fc, se_target_gene.p5.exp$age_fc)
abline(lm(se_target_gene.p5.exp$age_fc ~ se_target_gene.p5.exp$yy1_kd_fc))

plot(se_target_gene.p15.exp$yy1_kd_fc, se_target_gene.p15.exp$age_fc)
abline(lm(se_target_gene.p15.exp$age_fc ~ se_target_gene.p15.exp$yy1_kd_fc))


# all gene
all_gene_exp = fpkm %>% dplyr::select('gene_id', 'gene_biotype', p5_mean, p15_mean)
all_gene_exp$aging_fc = log2((all_gene_exp$p15_mean+0.1)/(all_gene_exp$p5_mean+0.1))

merge(all_gene_exp, yy1.fpkm, by='gene_id') -> all_gene_exp

all_gene_exp$yy1_kd_fc = log2((all_gene_exp$YY1_KD.FPKM+0.1)/(all_gene_exp$NT.FPKM+0.1))

plot(all_gene_exp$yy1_kd_fc, all_gene_exp$aging_fc)
abline(lm(all_gene_exp$aging_fc ~ all_gene_exp$yy1_kd_fc))

```
























