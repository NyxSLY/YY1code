---
title: "super enhancer data analysis"
output: html_notebook
---

```{r}
library(rtracklayer)
library(ggplot2)
p5_se='/Volumes/luyang-1/work/3dGenome/enhancer/superEnhancer/P5/P5_superEnhancer_Gateway_SuperEnhancers.bed'
p15_se='/Volumes/luyang-1/work/3dGenome/enhancer/superEnhancer/P15/P15_superEnhancer_Gateway_SuperEnhancers.bed'
p5_se = import(p5_se)
p15_se = import(p15_se)
boxplot(width(p5_se), width(p15_se), outline=F)
library(ggplot2)

d1 = data.frame(group='1', value=width(p5_se))
d2 = data.frame(group='2', value=width(p15_se))

boxplot_filter = function(x){
    l = boxplot.stats(x$value)$stats[c(1, 5)]
    if(nrow(x[x$value>l[2],])>0){
        x[x$value>l[2],]$value=l[2]*1.2
    }
    if(nrow(x[x$value<l[1],])>0){
        x[x$value<l[1],]$value=l[1]*0.8
    }
    return(x)
}

dat = rbind(boxplot_filter(d1), boxplot_filter(d2))
pp = ggplot(dat, aes(x=group,y=value,fill=group))+
  geom_boxplot(outlier.shape=NA) + 
  scale_fill_manual(values=c("#F6402A", "#4F4DB3","#FDF6E4")) +
  theme_classic()

pdf('boxplot of SE length.pdf', w=4, h=5)
print(pp)
dev.off()
```

### overlaps between up/dn regions
```{R}
vennplot = function(dat1, dat2, name1, name2){
    o = findOverlaps(dat1, dat2)
    library(VennDiagram)
    library(dplyr)
    grid.newpage()
    draw.pairwise.venn(area1=dim(data.frame(dat1))[1], area2=dim(data.frame(dat2))[1], 
                       cross.area=subjectHits(o) %>% unique() %>% length(), 
                       category = c(name1, name2), 
                       lty = 1,
                       fill = c("red", "blue"), 
                       col = rep("black",2),
                       lwd = rep(1,2),
                       alpha = c(0.5, 0.5), 
                       cat.pos = c(0,0), 
                       cex=1.6,
                       cat.dist = rep(0.025, 2))
}
pdf('super enhancer overlap.pdf', w=5, h=5)
vennplot(p5_se, p15_se, 'P5', 'P15')
dev.off()

```

```{r}
p5_sp = p5_se[!overlapsAny(p5_se ,p15_se)]
p15_sp = p15_se[!overlapsAny(p15_se ,p5_se)]
common = p5_se[overlapsAny(p5_se ,p15_se)]
export(p5_sp,'p5_sp_SE.bed','BED')
export(p15_sp,'p15_sp_SE.bed','BED')
export(common,'common_SE.bed','BED')

boxplot(width(p5_sp), width(p15_sp), width(common), outline=F)
```

# super enhancer don't overlap with exons
```{r}
library(EnsDb.Hsapiens.v75)
hg19=EnsDb.Hsapiens.v75

exon = exons(hg19)
exon = data.frame(exon)
exon$seqnames = paste0('chr', exon$seqnames)
exon = makeGRangesFromDataFrame(exon)



sum(overlapsAny(p15_se, exon))
length(p15_se)

p15_se[!overlapsAny(p15_se, exon),] -> p15_donotOverlapWithExon
export(p15_donotOverlapWithExon,'p15_SE_donotOverlapWithExon.bed','BED')

p5_se[!overlapsAny(p5_se, exon),] -> p5_donotOverlapWithExon
export(p5_donotOverlapWithExon,'p5_SE_donotOverlapWithExon.bed','BED')



GenomicRanges::setdiff(p15_se[1], exon, ignore.strand=T)


sum(overlapsAny(p15_se, exon))

getwd()


# control
shift(p5_donotOverlapWithExon, 10000) -> p5_donotOverlapWithExon_control
p5_donotOverlapWithExon_control[!overlapsAny(p5__donotOverlapWithExon_control, exon)] -> p5_donotOverlapWithExon_control
write.table(p5_donotOverlapWithExon_control, 'p5_donotOverlapWithExon_control.bed', quote = F,
            row.names = F, col.names = F, sep = '\t')
```


# make non-exon region from SE
```{r}
gr = GRanges()
for(i in 1:length(common)){
  GenomicRanges::setdiff(common[i], exon, ignore.strand=T) -> new
  new$name = common[i]$name
  gr = c(gr, new)
}

# make GTF file
data.frame(gr) -> gr
gr$type = 'commonSE'
gr$source = 'SE'
gr$empty = '.'
gr$strand = '+'
gr$empty1 = '.'
gr$group = 'exon'
gr$info = paste0('gene_id "', gr$name, '"; transcript_id "', gr$name, '"; gene_type "',gr$type,'";')

gr %>% select(seqnames, source, group, start, end, empty, strand, empty1, info) -> gtf
write.table(gtf, 'common_SE_notNoExon.gtf', row.names = F, col.names = F, quote = F, sep='\t')

run <- function(common, name){
  gr = GRanges()
  for(i in 1:length(common)){
    GenomicRanges::setdiff(common[i], exon, ignore.strand=T) -> new
    new$name = common[i]$name
    gr = c(gr, new)
  }
  # make GTF file
  data.frame(gr) -> gr
  gr$type = 'commonSE'
  gr$source = 'SE'
  gr$empty = '.'
  gr$strand = '+'
  gr$empty1 = '.'
  gr$group = 'exon'
  gr$info = paste0('gene_id "', gr$name, '"; transcript_id "', gr$name, '"; gene_type "',gr$type,'";')
  
  gr %>% select(seqnames, source, group, start, end, empty, strand, empty1, info) -> gtf
  write.table(gtf, paste0(name,'_SE_notNoExon.gtf'), row.names = F, col.names = F, quote = F, sep='\t')
}

run(p5_sp, name='p5_se')
run(p15_sp, name='p15_se')




```







```{}
source('/Users/Luyang/scripts/R/commonAPI.R')
p5_sp_mat = readmat('/Volumes/luyang-1/work/3dGenome/enhancer/superEnhancer/heatmap/p5_sp_SE_matrixCryptic.mat.gz')
p5_sp_mat$rowsumDat

p15_sp_mat = readmat('/Volumes/luyang-1/work/3dGenome/enhancer/superEnhancer/heatmap/p15_sp_SE_matrixCryptic.mat.gz')
p15_sp_mat$rowsumDat

common_mat = readmat('/Volumes/luyang-1/work/3dGenome/enhancer/superEnhancer/heatmap/common_SE_matrixCryptic.mat.gz')
common_mat$rowsumDat

dim(p5_sp_mat$gene)
dim(p15_sp_mat$gene)
dim(common_mat$gene)

a = p5_sp_mat$rowsumDat[1:200,]
b = p15_sp_mat$rowsumDat[1:200,]
c = common_mat$rowsumDat[1:200,]

dat = data.frame(a[,c(1,2)],b[,c(1,2)],c[,c(1,2)])
library(pheatmap)

cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}
 
data_subset_norm <- t(apply(dat, 1, cal_z_score))
pheatmap(data_subset_norm)

```

```{R}
se = union(p5_se, p15_se)
hist(log10(width(se)))
```

```{r}
p5_enh = import('/Volumes/luyang-1/work/3dGenome/enhancer/superEnhancer/P5/P5_superEnhancer_Gateway_Enhancers.bed',
'BED')
p15_enh = import('/Volumes/luyang-1/work/3dGenome/enhancer/superEnhancer/P15/P15_superEnhancer_Gateway_Enhancers.bed',
'BED')

p5_enh = p5_enh[!overlapsAny(p5_enh, p5_se)]
p15_enh = p15_enh[!overlapsAny(p15_enh, p15_se)]

enh = p5_enh[overlapsAny(p5_enh,p15_enh)]
enh = intersect(p5_enh,p15_enh)

enh.notpro = enh[!overlapsAny(enh, pros)]

enh.notpro = data.frame(enh.notpro)
enh.notpro$id = paste(enh.notpro$seqnames, enh.notpro$start, enh.notpro$end, sep='_')
enh.notpro = enh.notpro[,c(1,2,3,6,4,5)]

write.table(enh.notpro, 'enhancer_notOnPromoter.bed', row.names = F, col.names = F, quote = F, sep = '\t')
```






