---
title: "DNA methylation methylKit"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_knit$set(dev.args = list(type = "cairo"))
library(ggplot2)
trace(grDevices::png, quote({
  if (missing(type) && missing(antialias)) {
    type <- "cairo-png"
    antialias <- "subpixel"
  }
}), print = FALSE)
plot=F

library(data.table)
library(stringr)
```



```{r}
library(methylKit)
library(dplyr)

file.list = list('P5_MSC_WGBS.deduplicated.bismark.cov',
                 'P15_MSC_WGBS.deduplicated.bismark.cov')

myobj = methRead(file.list,
                 sample.id=list('P5', 'P15'),
                 assembly='hg19',
                 context='CpG',
                 pipeline='bismarkCoverage',
                 treatment=c(0,1),
                 mincov = 10)

```

# statistics

## total methylated C vs unmethylated C
```{r}

p5.df = myobj[1] %>% data.frame
p15.df = myobj[2] %>% data.frame

totalC.p5 = sum(p5.df$coverage)
totalC.p15 = sum(p15.df$coverage)

methylC.p5 = sum(p5.df$numCs)
methylC.p15 = sum(p15.df$numCs)


if(plot){
  pdf('total methylC pct.pdf', w=4, h=6)
  barplot(c(methylC.p5/totalC.p5, methylC.p15/totalC.p15), ylim=c(0,0.8))
  dev.off()
}

barplot(c(methylC.p5/totalC.p5, methylC.p15/totalC.p15), ylim=c(0,0.8))

```





```{r fig.height=7}
if(plot){
  pdf('P5 CpG coverage.pdf', w=10, h=7)
  getMethylationStats(myobj[[1]],plot=TRUE,both.strands=FALSE)
  dev.off()
  pdf('P15 CpG coverage.pdf', w=10, h=7)
  getMethylationStats(myobj[[2]],plot=TRUE,both.strands=FALSE)
  dev.off()
  
}
getMethylationStats(myobj[[1]],plot=T,both.strands=FALSE)
getCoverageStats(myobj[[1]],plot=TRUE,both.strands=FALSE)
getMethylationStats(myobj[[2]],plot=T,both.strands=FALSE)
getCoverageStats(myobj[[2]],plot=TRUE,both.strands=FALSE)


```

# Filtering samples based on read coverage
```{}
filtered.myobj=filterByCoverage(myobj,lo.count=10,lo.perc=NULL,
                                      hi.count=NULL,hi.perc=99.9)
```



# Merge into 1kb windows and check sig diff region by Fisher
```{}
# because 1kb windows became the basic unit, single-base low coverage threshold can lower than 10. Use 5 here

myobj_lowCov = methRead(file.list,
           sample.id=list("P5","P15"),
           assembly="hg19",
           treatment=c(0,1),
           context="CpG",
           mincov = 5,
           pipeline='bismarkCoverage',
           )


filtered.myobj=filterByCoverage(myobj_lowCov,lo.count=5,lo.perc=NULL,
                                      hi.count=NULL,hi.perc=99.9)

noamlized.myobj = normalizeCoverage(filtered.myobj,method="median")

tiles = tileMethylCounts(noamlized.myobj,win.size=1000,step.size=1000,cov.bases = 10)
head(tiles[[1]],3)
meth=unite(tiles, destrand=FALSE)
meth

myDiff=calculateDiffMeth(meth)

myDiff.df = data.frame(myDiff)


myDiff25p.hyper=getMethylDiff(myDiff,difference=25,qvalue=0.01,type="hyper")
myDiff25p.hypo=getMethylDiff(myDiff,difference=25,qvalue=0.01,type="hypo")

# merge p-value and read count
meth.df = data.frame(meth)
meth.df$id = paste(meth.df$chr, meth.df$start, meth.df$end, sep='_')

myDiff.df = data.frame(myDiff)
myDiff.df$id = paste(myDiff.df$chr, myDiff.df$start, myDiff.df$end, sep='_')
dat.temp = myDiff.df[,c('id', 'pvalue', 'qvalue', 'meth.diff')]

data = merge(meth.df, dat.temp, by='id')

data %>% arrange(chr, start, end) -> data

# extract DMR - threshold defined as q-value < 0.01, absolute % difference > 25 

#sig = data %>% filter(qvalue<0.01, abs(meth.diff)>25)

sig = data %>% dplyr::filter(qvalue<0.01, abs(meth.diff)>25)
#sig = data %>% filter(qvalue<0.01)

#hist(abs(sig$meth.diff))
sig_dn = sig %>% dplyr::filter(meth.diff<0)
sig_up = sig %>% dplyr::filter(meth.diff>0)

hist(sig$meth.diff, xlab='CpG methylation % difference (P15 - P5)', main = '', freq=FALSE)

if(plot){
  pdf('nubmer of hyper and hypo DMR barplot.pdf', w=4, h=6)
  barplot(c(nrow(sig_up), nrow(sig_dn)), names=c('hyper', 'hypo'), ylim=c(0, 40000), main='number of hyper/hypo DMR region')
  dev.off()
  
  pdf('Length of hyper/hypo DMR region.pdf', w=4, h=6)
  barplot(c(sum(sig_up.gr@ranges@width)/1000000,sum(sig_dn.gr@ranges@width)/1000000), names=c('hyper', 'hypo'), 
        main='Length (Mb) of hyper/hypo DMR region', ylab='total length (Mb)')
  dev.off()
}

barplot(c(nrow(sig_up), nrow(sig_dn)), names=c('hyper', 'hypo'), ylim=c(0, 40000), main='number of hyper/hypo DMR region')


makeGRangesFromDataFrame(sig_up) %>% reduce -> sig_up.gr
makeGRangesFromDataFrame(sig_dn) %>% reduce -> sig_dn.gr

barplot(c(sum(sig_up.gr@ranges@width)/1000000,sum(sig_dn.gr@ranges@width)/1000000), names=c('hyper', 'hypo'), 
        main='Length (Mb) of hyper/hypo DMR region', ylab='total length (Mb)', ylim=c(0,35))

dim(sig)

```

# export
methylKit has more false positive
```{}
sig= sig %>% arrange(qvalue)
write.csv(sig, 'DMR.csv',quote=FALSE, row.names = FALSE)
head(sig)
```


# check sig diff at single-base level
```{r}
library(dplyr)
filtered.myobj=filterByCoverage(myobj,lo.count=10,lo.perc=NULL,
                                      hi.count=NULL,hi.perc=99.9)

noamlized.myobj = normalizeCoverage(filtered.myobj,method="median")


meth=unite(noamlized.myobj, destrand=FALSE)
meth

myDiff=calculateDiffMeth(meth)

myDiff.df = data.frame(myDiff)


myDiff25p.hyper=getMethylDiff(myDiff,difference=25,qvalue=0.01,type="hyper")
myDiff25p.hypo=getMethylDiff(myDiff,difference=25,qvalue=0.01,type="hypo")

# merge p-value and read count
meth.df = data.frame(meth)
meth.df$id = paste(meth.df$chr, meth.df$start, meth.df$end, sep='_')

myDiff.df = data.frame(myDiff)
myDiff.df$id = paste(myDiff.df$chr, myDiff.df$start, myDiff.df$end, sep='_')
dat.temp = myDiff.df[,c('id', 'pvalue', 'qvalue', 'meth.diff')]

data.single = merge(meth.df, dat.temp, by='id')

data.single %>% arrange(chr, start, end) -> data.single

write.csv(meth.df, 'see.csv')
# extract DMR - threshold defined as q-value < 0.01, absolute % difference > 25 

#sig = data %>% filter(qvalue<0.01, abs(meth.diff)>25)

sig.singleBase = data.single %>% dplyr::filter(qvalue<0.01, abs(data.single$meth.diff)>25)
sig_dn.singleBase = sig.singleBase %>% dplyr::filter(meth.diff<0)
sig_up.singleBase = sig.singleBase %>% dplyr::filter(meth.diff>0)

sig.singleBase.wide = data.single %>% dplyr::filter(qvalue<0.01)
sig_dn.singleBase.wide = sig.singleBase.wide %>% dplyr::filter(meth.diff<0)
sig_up.singleBase.wide = sig.singleBase.wide %>% dplyr::filter(meth.diff>0)

if(plot){
  pdf('CpG methylation % difference (P15 - P5).pdf', w=8, h=6)
  hist(sig.singleBase$meth.diff, xlab='CpG methylation % difference (P15 - P5)', main = '', freq=FALSE, col='skyblue')
  dev.off()
}
hist(sig.singleBase$meth.diff, xlab='CpG methylation % difference (P15 - P5)', main = '', freq=FALSE, col='skyblue')

sig.singleBase= sig.singleBase %>% arrange(qvalue)
write.csv(sig.singleBase, 'singleBase_DNAmethylation.csv',quote=FALSE, row.names = FALSE)


saveRDS(data.single, file='data.single.Rds')

```
## Obvious CpG DNA methylation loss 



# overlap with age-specific YY1
```{r}
data.single = readRDS('data.single.Rds')

fread('yy1_common.tsv') -> yy1.bindingsite.common
fread('yy1_p5.tsv') -> yy1.bindingsite.p5.sp
fread('yy1_p15.tsv') -> yy1.bindingsite.p15.sp

process = function(x){
  region = x$sequence_name

  chr = str_split(region, ':') %>%
  sapply(., '[[', 1)

  start = str_split(region, ':') %>%
    sapply(., '[[', 2) %>% str_split(., '-') %>% sapply(., '[[', 1)
  
  ss = x$start
  ee = x$stop
  
  start = as.numeric(start) + as.numeric(ss)
  end = as.numeric(start) + as.numeric(ee)
  
  
  df = data.frame(chr=chr, start=start, end=end)
  makeGRangesFromDataFrame(df) -> x.gr
  GenomicRanges::reduce(x.gr)
}

process(yy1.bindingsite.p5.sp) -> yy1.bs.p5.sp
process(yy1.bindingsite.p15.sp) -> yy1.bs.p15.sp
process(yy1.bindingsite.common) -> yy1.bs.common.sp


makeGRangesFromDataFrame(data.single, keep.extra.columns = T) -> data.single.gr

yy1_p5.sp.data = data.single.gr[overlapsAny(data.single.gr, yy1.bs.p5.sp)]
yy1_p15.sp.data = data.single.gr[overlapsAny(data.single.gr, yy1.bs.p15.sp)]
yy1_common.data = data.single.gr[overlapsAny(data.single.gr, yy1.bs.common.sp)]

boxplot(yy1_p5.sp.data$meth.diff,
        yy1_p15.sp.data$meth.diff,
        yy1_common.data$meth.diff)



# new 20231123
plot.dat1 = data.frame(key='p5sp', value=yy1_p5.sp.data$meth.diff)
plot.dat2 = data.frame(key='p15sp', value=yy1_p15.sp.data$meth.diff)
plot.dat3 = data.frame(key='common', value=yy1_common.data$meth.diff)

plot.dat = rbind(plot.dat1, plot.dat2, plot.dat3)
plot.dat$key = factor(plot.dat$key, levels = c('p5sp', 'p15sp', 'common'))

# no need to remove outlier, no big outlier
saveRDS(plot.dat, 'DNA methylation.rds')
p = ggplot(plot.dat, aes(x = key, y = value, fill = key)) +
  geom_jitter(width = 0.2, aes(color = key), size = 1.5, alpha = 0.1) +  # Add jitter with matching color
  #geom_violin(trim=T, adjust=2, width=0.8, color='black', alpha=0.7) +  # Hide the outliers in the boxplot
  geom_boxplot(outlier.shape = NA, alpha=0.5) +
  scale_colour_brewer(palette = "Set2") +
  scale_fill_brewer(palette = "Set2") +
  theme_classic()




```

# make better plot
```{r}
#source(file.path('D:/BaiduNetdiskWorkspace/OneDrive/scripts/R/rain_cloud/example_plot_code.r', '..',"halfViolinPlots.R"))
source(file.path('C:\\Users\\Luyang_pc\\Documents\\BaiduNetdiskWorkspace/OneDrive/scripts/R/rain_cloud/example_plot_code.r', '..',"halfViolinPlots.R"))
library("cowplot") 

# dat1 = data.frame(group = 'P5 sp', value=yy1_p5.sp.data$meth.diff)
# dat2 = data.frame(group = 'P15 sp', value=yy1_p15.sp.data$meth.diff)
# dat3 = data.frame(group = 'Common', value=yy1_common.data$meth.diff)

dat1 = data.frame(group = '1', value=yy1_p5.sp.data$meth.diff)
dat2 = data.frame(group = '2', value=yy1_p15.sp.data$meth.diff)
dat3 = data.frame(group = '3', value=yy1_common.data$meth.diff)

plotData = rbind(dat1, dat2, dat3)
colnames(plotData) = c("condition", 'value')

ggplot(plotData, aes(x = condition, y = value, fill = condition, color = condition)) +
  ggtitle("MethylC level diff at YY1 binding motif") +
  ylab("MethylC level diff") +
  scale_x_discrete(labels= c("P5 sp.", 'P15 sp.', 'Common')) + 
  theme_cowplot() +
  scale_shape_identity() +
  theme(legend.position = "none",
             plot.title = element_text(size = 20),
             axis.title = element_text(size = 15),
             axis.text = element_text(size = 15),
             axis.text.x = element_text(angle = 0, 
                hjust = 0,
                vjust = 0)) +
  scale_colour_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  geom_point(position = position_jitter(0.1), 
             size = 2, 
             alpha = 1, 
             aes(shape = 16)) +
  geom_flat_violin(position = position_nudge(x = -0.1, y = 0),
                   #aes(fill=NA),
             adjust = 3,
             alpha = 0.3, 
             trim = F, 
             scale = "area") +
  geom_boxplot(aes(x = as.numeric(condition) - 0.35, y = value), 
             notch = FALSE, 
             width = 0.3, 
             varwidth = FALSE, 
             outlier.shape = NA, 
             alpha = 0.6, 
             colour = "black", 
             show.legend = FALSE) -> p1;p1

output_ppt_m(p1, 'DNA methylation changes at YY1 binding motif.pptx')


```

# Age-specific YY1 on SE (and not on SE)
```{r}
#load('C:\\Users\\Luyang_pc\\Documents\\BaiduNetdiskWorkspace\\OneDrive\\work\\3dGenome\\Find_enhancer_target\\BASE.RData')
load('D:\\BaiduNetdiskWorkspace\\OneDrive\\work\\3dGenome\\Find_enhancer_target\\BASE.RData')

se$type = 'common'
se[grepl('p5_', se$name)]$type = 'p5'
se[grepl('p15_', se$name)]$type = 'p15'
se$type %>% table()

yy1.bs.p5.sp[overlapsAny(yy1.bs.p5.sp, se)] -> yy1.bs.p5.sp.onSE
yy1.bs.p15.sp[overlapsAny(yy1.bs.p15.sp, se)] -> yy1.bs.p15.sp.onSE
yy1.bs.common.sp[overlapsAny(yy1.bs.common.sp, se)] -> yy1.bs.common.sp.onSE

makeGRangesFromDataFrame(data.single, keep.extra.columns = T) -> data.single.gr

yy1_p5.sp.data1 = data.single.gr[overlapsAny(data.single.gr, yy1.bs.p5.sp.onSE)]
yy1_p15.sp.data1 = data.single.gr[overlapsAny(data.single.gr, yy1.bs.p15.sp.onSE)]
yy1_common.data1 = data.single.gr[overlapsAny(data.single.gr, yy1.bs.common.sp.onSE)]

boxplot(yy1_p5.sp.data1$meth.diff,
        yy1_p15.sp.data1$meth.diff,
        yy1_common.data1$meth.diff)



# No (almost) canornical YY1 binding site overlap with DNA methylation data


```




```{}
library(data.table)
library(commonAPI)
#ct = fread('D:\\work\\DNA_methylation\\bed\\final.preCallPeak.bed')
yy1_p5.sp = fread('young_sp_YY1.bed')
yy1_p15.sp = fread('old_sp_YY1.bed')

yy1_p5.sp = toGR(yy1_p5.sp)
yy1_p15.sp = toGR(yy1_p15.sp)


makeGRangesFromDataFrame(data.single, keep.extra.columns = T) -> data.single.gr

yy1_p5.sp.data = data.single.gr[overlapsAny(data.single.gr, yy1_p5.sp)]
yy1_p15.sp.data = data.single.gr[overlapsAny(data.single.gr, yy1_p15.sp)]


boxplot(yy1_p5.sp.data$meth.diff,
yy1_p15.sp.data$meth.diff)


boxplot(yy1_p5.sp.data$numCs1/yy1_p5.sp.data$coverage1-yy1_p5.sp.data$numCs2/yy1_p5.sp.data$coverage2)
abline(h=0)

boxplot(yy1_p15.sp.data$numCs1/yy1_p15.sp.data$coverage1-yy1_p15.sp.data$numCs2/yy1_p15.sp.data$coverage2)
abline(h=0)


# For each YY1 binding site as a unit, examine DNA methylation changes

run_def = function(yy1_p5.sp, data.single){
  data.single.gr = makeGRangesFromDataFrame(data.single, keep.extra.columns = T)
  findOverlapPairs(yy1_p5.sp, data.single.gr) -> x
  x = data.frame(x)
  x$first.id = paste(x$first.X.seqnames, x$first.X.start, x$first.X.end, sep='_')
  new.tab = data.frame(yy1=x$first.id, methy=x$second.X.id)
  merge(new.tab, data.single, by.x='methy', by.y='id') -> dat.yy1
  dat.yy1$pct_methyC1 = dat.yy1$numCs1/dat.yy1$coverage1
  dat.yy1$pct_methyC2 = dat.yy1$numCs2/dat.yy1$coverage2
  tapply(dat.yy1$pct_methyC1, dat.yy1$yy1, mean) -> dat.yy1.y
  tapply(dat.yy1$pct_methyC2, dat.yy1$yy1, mean) -> dat.yy1.o
  return(data.frame(yy1_id = names(dat.yy1.y),
                    ave_methy_young=dat.yy1.y,
                    ave_methy_old=dat.yy1.o))
}

run_def(yy1_p5.sp, data.single) -> dat.yy1.y
run_def(yy1_p15.sp, data.single) -> dat.yy1.o


boxplot(dat.yy1.y$ave_methy_old/dat.yy1.y$ave_methy_young, outline=F)
boxplot(dat.yy1.o$ave_methy_old/dat.yy1.o$ave_methy_young, outline=F)

boxplot(dat.yy1.y$ave_methy_young, dat.yy1.y$ave_methy_old, outline=F)
boxplot(dat.yy1.o$ave_methy_young, dat.yy1.o$ave_methy_old, outline=F)
boxplot(dat.yy1.y$ave_methy_young, dat.yy1.y$ave_methy_old, 
        dat.yy1.o$ave_methy_young, dat.yy1.o$ave_methy_old,
        sig.singleBase,outline=F)


hist(dat.yy1.y$ave_methy_young)
hist(dat.yy1.y$ave_methy_old)

hist(dat.yy1.o$ave_methy_young)
hist(dat.yy1.o$ave_methy_old)

hist(data.single$numCs1/data.single$coverage1)
hist(data.single$numCs2/data.single$coverage2)


data.frame(findOverlapPairs(yy1_p5.sp, data.single.gr))

findOverlaps(yy1_p5.sp, data.single.gr)



ct_flank_data.df = data.frame(ct_flank_data)
ct_flank_data.df$methyC_pct1 = ct_flank_data.df$numCs1/ct_flank_data.df$coverage1
ct_flank_data.df$methyC_pct2 = ct_flank_data.df$numCs2/ct_flank_data.df$coverage2


if(plot){
  pdf('methyC pct at in cTSS region (1kb region).pdf', w=4, h=6)
  boxplot(ct_flank_data.df$methyC_pct1, ct_flank_data.df$methyC_pct2, outline=F, names=c('P5', 'P15'), xlab='methyC %')
  dev.off()
}





```
