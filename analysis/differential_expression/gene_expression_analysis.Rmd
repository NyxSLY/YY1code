---
title: "gene expression changes and YY1 and gene"
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
#g = getGenomeAndMask(BSgenome.Hsapiens.UCSC.hg19.masked)
library(ggplot2)

#source('/Users/Luyang/scripts/R/commonAPI.R')
#source('/Users/Luyang/scripts/R/annotation.R')

#source('D:\\work\\script\\R\\commonAPI.R')
#source('D:\\work\\script\\R\\annotation.R')
library(commonAPI)

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

# functions
```{r message=F}

session::restore.session('../functions.RData')
```

# database preload

```{r}
library(EnsDb.Hsapiens.v75)
hg19=EnsDb.Hsapiens.v75
pros = promoters(hg19, upstream = 1000, downstream=1000)
pros = data.frame(pros)
pros$seqnames = as.character(pros$seqnames)
pros$seqnames = paste0('chr',pros$seqnames)
pros = makeGRangesFromDataFrame(pros, keep.extra.columns = T)
#pros = pros[pros$tx_biotype=='protein_coding',]


# import YY1 peak
yy1_p15 = readbed('../P15_peaks.narrowPeak')
yy1_p5 = readbed('../P5_peaks.narrowPeak')
yy1 =reduce(union(yy1_p15, yy1_p5))

# import CTCF peak
ctcf_p5='../P5_CTCF_macs2_peaks.narrowPeak'
ctcf_p15='../P15_CTCF_macs2_peaks.narrowPeak'
ctcf_p5 = import(ctcf_p5)
ctcf_p15 = import(ctcf_p15)
ctcf = union(ctcf_p5, ctcf_p15)


# fpkm
aging.fpkm = fread('../fpkm.csv')
# YY1 KD
yy1kd.fpkm = read.csv('../YY1_KD_DEG/fpkm.csv')


aging.fpkm$aging_mean = rowMeans(aging.fpkm[,2:7])
yy1kd.fpkm$yy1kd_mean = rowMeans(yy1kd.fpkm[,2:7])

#fc_data = '../DESEQ2_total_results.csv'
#fc_data = 'C:\\Users\\Luyang.000\\OneDrive\\work\\3dGenome\\Find_enhancer_target\\DESEQ_only_use_data2020/DESEQ2_total_results.csv'
fc_data = '../remove_batch_effect/DESEQ2_total_results_RUV.csv'
fc_data = fread(fc_data)
colnames(fc_data)[1] = 'gene_id'
fc_data.short = fc_data %>% dplyr::select(gene_id,
                                          baseMean,
                                          log2FoldChange,
                                          padj)


yy1kd_fc_data = '../YY1_KD_DEG/DESEQ2_total_results.csv'
yy1kd_fc_data = fread(yy1kd_fc_data)
yy1kd_fc_data.short = yy1kd_fc_data %>% dplyr::select(gene_id,
                                                      baseMean,
                                          log2FoldChange,
                                          padj)

yy1kd_fc_data.short = yy1kd_fc_data.short[!is.na(yy1kd_fc_data.short$padj),]
names(yy1kd_fc_data.short)[2:4] = c('yy1kd_baseMean','yy1kd_log2FoldChange', 'yy1kd_padj')

yy1kd.fpkm.expressed=yy1kd.fpkm[yy1kd.fpkm$yy1kd_mean>1, ]
yy1kd.fpkm.expressed %>% add_fc(., yy1kd_fc_data.short) -> yy1kd.fpkm.expressed
```

# Identify genes in the target region
```{r echo=T}
p5pair = fread('../P5_SE_and_its_target.csv')
p15pair = fread('../P15_SE_and_its_target.csv')

p5pair$target_region_id = paste(p5pair$target_region_chr,
                                p5pair$target_region_start,
                                p5pair$target_region_end,
                                sep='_')

p15pair$target_region_id = paste(p15pair$target_region_chr,
                                p15pair$target_region_start,
                                p15pair$target_region_end,
                                sep='_')

get_genes_in_target_region(p5pair, pros) -> p5pair.genes
get_genes_in_target_region(p15pair, pros) -> p15pair.genes

add_fpkm(p5pair.genes, aging.fpkm, yy1kd.fpkm) -> p5pair.genes.fpkm
add_fpkm(p15pair.genes, aging.fpkm, yy1kd.fpkm) -> p15pair.genes.fpkm

# only use genes with FPKM > 1
p5pair.genes.fpkm %>% dplyr::filter(aging_fpkm_mean > 1) -> p5pair.genes.fpkm
p15pair.genes.fpkm %>% dplyr::filter(aging_fpkm_mean > 1) -> p15pair.genes.fpkm


# add pair id
p5pair.genes.fpkm$target_start = sapply(strsplit(p5pair.genes.fpkm$target_region_id, "_"),'[[', 2)
p5pair.genes.fpkm$pair = paste(p5pair.genes.fpkm$se_id, p5pair.genes.fpkm$target_start, sep='_')


p15pair.genes.fpkm$target_start = sapply(strsplit(p15pair.genes.fpkm$target_region_id, "_"),'[[', 2)
p15pair.genes.fpkm$pair = paste(p15pair.genes.fpkm$se_id, p15pair.genes.fpkm$target_start, sep='_')


# add fold change and pvalue
p5pair.genes.fpkm %>% add_fc(., fc_data.short) -> p5pair.genes.fpkm
p15pair.genes.fpkm %>% add_fc(., fc_data.short) -> p15pair.genes.fpkm

p5pair.genes.fpkm %>% add_fc(., yy1kd_fc_data.short) -> p5pair.genes.fpkm
p15pair.genes.fpkm %>% add_fc(., yy1kd_fc_data.short) -> p15pair.genes.fpkm

```



# expression
## load data
```{}
P5_pair_withYY1down = fread('../P5_pair_withYY1down.csv')
P15_pair_withYY1up = fread('../P15_pair_withYY1up.csv')

P5_pair_withYY1down$target_region_id = paste(P5_pair_withYY1down$target_region_chr,
                                              P5_pair_withYY1down$target_region_start,
                                              P5_pair_withYY1down$target_region_end,
                                              sep='_')

P15_pair_withYY1up$target_region_id = paste(P15_pair_withYY1up$target_region_chr,
                                              P15_pair_withYY1up$target_region_start,
                                              P15_pair_withYY1up$target_region_end,
                                              sep='_')

# remove paradox gene
P5_pair_withYY1down$gene_id %in% P15_pair_withYY1up$gene_id %>% P5_pair_withYY1down$gene_id[.] -> paradox.gene

P5_pair_withYY1down %>% dplyr::filter(!gene_id %in% paradox.gene) -> P5_pair_withYY1down
P15_pair_withYY1up %>% dplyr::filter(!gene_id %in% paradox.gene) -> P15_pair_withYY1up

```




```{}
hic.dat = fread('../P5_pair_withYY1down_withHiCinfo.csv')
hic.dat.dn = hic.dat[hic.dat$P15<hic.dat$P5,]

P5_pair_withYY1down.hic = P5_pair_withYY1down[P5_pair_withYY1down$pair %in% hic.dat.dn$V1,]


merge(P5_pair_withYY1down.hic, fc_data.short, by='gene_id') -> P5_pair_withYY1down.fc
P5_pair_withYY1down.fc[!duplicated(P5_pair_withYY1down.fc$gene_id),] -> P5_pair_withYY1down.fc

P5_pair_withYY1down.fc$log2FoldChange %>% boxplot(outline=F)
abline(h=0)

P5_pair_withYY1down.fc %>% dplyr::filter(padj < 0.05) -> sig


nrow(sig)/nrow(P5_pair_withYY1down.fc)

nrow(fc_data.short %>% dplyr::filter(padj<0.05))/nrow(fc_data.short)


# matrix(c(11,192, 967, 57773),
#        nrow = 2,
#        dimnames = list(Guess = c("Milk", "Tea"),
#                        Truth = c("Milk", "Tea"))) -> tt
# 
# 
# fisher.test(tt)
```

```{}

deseq2_scatter_plot <- function(yy1.deseq2){
  yy1.deseq2 = yy1.deseq2[!is.na(yy1.deseq2$padj),]
  dat.scater = yy1.deseq2 %>% dplyr::select(gene_id, baseMean,log2FoldChange, padj)
  if(sum(dat.scater$padj==0)!=0){
    dat.scater[dat.scater$padj==0,]$padj = 1e-50
  }
  dat.scater %>% mutate(sig=cut(padj, breaks=c(0,0.05,1), labels=c('sig', 'notsig'))) -> dat.scater
  dat.scater %>% mutate(direction=cut(log2FoldChange, breaks=c(-Inf, 0, Inf), labels=c('DN', 'UP'))) -> dat.scater
  dat.scater %>% mutate(group=paste0(dat.scater$sig, dat.scater$direction))-> dat.scater
  dat.scater %>% mutate(point_size = cut(padj, breaks=c(0,0.05,1), labels=c(1.8, 1))) -> dat.scater
  dat.scater %>% mutate(point_size=as.numeric(as.character(point_size))) -> dat.scater
  
  main <- ggplot(dat.scater, aes(x=log10(baseMean), y=log2FoldChange, color=group)) + 
    geom_point(aes(size=point_size)) +
    theme_classic() + 
    theme(legend.position = "none") +
    scale_color_manual(values=c('gray', 'gray', 'blue','red')) + 
    scale_size(range = c(1, 1.8))
  
  density <- density_plot(dat.scater)
  return(list(main=main, density=density))
}

density_plot <- function(df){
  xdensity <- ggplot(df, aes(log2FoldChange, fill='#E69F00')) + 
  geom_density(alpha=.5, adjust = 1) + 
  scale_fill_manual(values = c('#E69F00')) + 
  theme(legend.position = "none") + 
    theme_classic()
}


deseq2_scatter_plot(P5_pair_withYY1down.fc) -> see


# scatterplot.g = ggplot(dat.scater, aes(x=log10(baseMean), y=log2FoldChange, color=group)) + geom_point(size=0.1) +
#   theme_classic() + 
#   scale_color_manual(values=c('gray', 'gray', 'blue','red'))
# 
# output_ppt(scatterplot.g, 'YY1KD RNA-seq scatter plot.pptx')
# scatterplot.g

output_ppt_m(see$main, 'Gain and Lost YY1 with age - scaterplot.pptx')

```


```{}
hic1.dat = fread('../P15_pair_withYY1up_withHiCinfo.csv')
hic1.dat.dn = hic1.dat[hic1.dat$P15<hic1.dat$P5,]

boxplot(hic1.dat$P5, hic1.dat$P15)

P15_pair_withYY1up.hic = P15_pair_withYY1up[P15_pair_withYY1up$pair %in% hic1.dat.dn$V1,]

merge(P15_pair_withYY1up.hic, fc_data.short, by='gene_id') -> P15_pair_withYY1up.fc
P15_pair_withYY1up.fc[!duplicated(P15_pair_withYY1up.fc$gene_id),] -> P15_pair_withYY1up.fc

P15_pair_withYY1up.fc$log2FoldChange %>% boxplot(outline=F)
abline(h=0)


deseq2_scatter_plot(P15_pair_withYY1up.fc) -> see1
see1$main + ylim(-1,1) -> new


deseq2_scatter_plot_yy1kd(P15_pair_withYY1up.fc %>% add_fc(yy1kd_fc_data.short))

#output_ppt_m(new, 'Gain and Lost YY1 with age - scaterplot.pptx')
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



# YY1 on enhancer
```{r}
YY1_on_se(p5pair.genes.fpkm, yy1_p5) -> p5pair.genes.fpkm.yy1
p5pair.genes.fpkm[!p5pair.genes.fpkm$pair %in% p5pair.genes.fpkm.yy1$pair,] -> p5pair.genes.fpkm.NOTyy1

boxplot(abs(p5pair.genes.fpkm.yy1$log2FoldChange), abs(p5pair.genes.fpkm.NOTyy1$log2FoldChange), outline=F)
boxplot(abs(p5pair.genes.fpkm.yy1$log2FoldChange), abs(p5pair.genes.fpkm$log2FoldChange), abs(fc_data.short$log2FoldChange), outline=F)

aging.fpkm.expressed=aging.fpkm[aging.fpkm$aging_mean>1, ]
aging.fpkm.expressed %>% add_fc(., fc_data.short) -> aging.fpkm.expressed
  
boxplot(p5pair.genes.fpkm.yy1$aging_fpkm_mean, p5pair.genes.fpkm$aging_fpkm_mean, aging.fpkm.expressed$aging_mean, outline=F)


p5pair.genes.fpkm %>% dplyr::filter(grepl('common', se_id)) -> p5pair.genes.fpkm.common
p5pair.genes.fpkm %>% dplyr::filter(!grepl('common', se_id)) -> p5pair.genes.fpkm.sp


whole_step <- function(p5pair.genes.fpkm, yy1_p5){
  YY1_on_se(p5pair.genes.fpkm, yy1_p5) -> p5pair.genes.fpkm.yy1
  p5pair.genes.fpkm[!p5pair.genes.fpkm$pair %in% p5pair.genes.fpkm.yy1$pair,] -> p5pair.genes.fpkm.NOTyy1
  
  boxplot(abs(p5pair.genes.fpkm.yy1$log2FoldChange), abs(p5pair.genes.fpkm.NOTyy1$log2FoldChange), outline=F, 
          names=c('with YY1', 'without YY1'), main='gene expression fold change')
  boxplot(abs(p5pair.genes.fpkm.yy1$log2FoldChange), abs(p5pair.genes.fpkm$log2FoldChange), abs(fc_data.short$log2FoldChange), outline=F,
          names=c('with YY1', 'without YY1', 'all expressed'), main='gene expression fold change')
  
  boxplot(abs(p5pair.genes.fpkm.yy1$log2FoldChange), abs(p5pair.genes.fpkm.yy1$yy1kd_log2FoldChange), outline=F,
          names=c('with YY1', 'with YY1 YY1KD'), main='YY1 KD fold change')

  aging.fpkm.expressed=aging.fpkm[aging.fpkm$aging_mean>1, ]
  boxplot(p5pair.genes.fpkm.yy1$aging_fpkm_mean, p5pair.genes.fpkm$aging_fpkm_mean, aging.fpkm.expressed$aging_mean, outline=F,
          names=c('with YY1', 'without YY1', 'all expressed'), main='gene expression (FPKM)')
  
  return(list(agine=get_sig_pct(p5pair.genes.fpkm.yy1, p5pair.genes.fpkm.NOTyy1),
              yy1kd=get_sig_pct1(p5pair.genes.fpkm.yy1, p5pair.genes.fpkm.NOTyy1)))
}


get_sig_pct <- function(p5pair.genes.fpkm.yy1, p5pair.genes.fpkm.NOTyy1){
  p5pair.genes.fpkm.yy1 %>% dplyr::filter(padj<0.05) -> p5pair.genes.fpkm.yy1.sig
  p5pair.genes.fpkm.NOTyy1 %>% dplyr::filter(padj<0.05) -> p5pair.genes.fpkm.NOTyy1.sig
  aging.fpkm.expressed %>% dplyr::filter(padj<0.05) -> aging.fpkm.expressed.sig
  
  return(list(yy1=p5pair.genes.fpkm.yy1.sig %>% nrow / nrow(p5pair.genes.fpkm.yy1),
              noYY1=p5pair.genes.fpkm.NOTyy1.sig %>% nrow / nrow(p5pair.genes.fpkm.NOTyy1),
              all=aging.fpkm.expressed.sig %>% nrow / nrow(aging.fpkm.expressed)))
}

get_sig_pct1 <- function(p5pair.genes.fpkm.yy1, p5pair.genes.fpkm.NOTyy1){
  print(names(p5pair.genes.fpkm.NOTyy1))
  p5pair.genes.fpkm.yy1 %>% dplyr::filter(yy1kd_padj<0.05) -> p5pair.genes.fpkm.yy1.sig
  p5pair.genes.fpkm.NOTyy1 %>% dplyr::filter(yy1kd_padj<0.05) -> p5pair.genes.fpkm.NOTyy1.sig
  yy1kd.fpkm.expressed %>% dplyr::filter(yy1kd_padj<0.05) -> yy1kd.fpkm.expressed
  
  return(list(yy1=p5pair.genes.fpkm.yy1.sig %>% nrow / nrow(p5pair.genes.fpkm.yy1),
              noYY1=p5pair.genes.fpkm.NOTyy1.sig %>% nrow / nrow(p5pair.genes.fpkm.NOTyy1),
              all=yy1kd.fpkm.expressed %>% nrow / nrow(aging.fpkm.expressed)))
}


whole_step(p5pair.genes.fpkm, yy1_p5)
whole_step(p15pair.genes.fpkm, yy1_p15)

aging.fpkm.expressed %>% 
  mutate(quantilegroup = ntile(aging_mean, 4)) -> aging.fpkm.expressed

aging.fpkm.expressed %>% dplyr::filter(quantilegroup==1) -> aging.fpkm.expressed.q1
aging.fpkm.expressed %>% dplyr::filter(quantilegroup==2) -> aging.fpkm.expressed.q2
aging.fpkm.expressed %>% dplyr::filter(quantilegroup==3) -> aging.fpkm.expressed.q3
aging.fpkm.expressed %>% dplyr::filter(quantilegroup==4) -> aging.fpkm.expressed.q4

boxplot(abs(p5pair.genes.fpkm.yy1$log2FoldChange), 
        abs(p5pair.genes.fpkm.NOTyy1$log2FoldChange)
        ,aging.fpkm.expressed.q1$log2FoldChange %>% abs,
        aging.fpkm.expressed.q2$log2FoldChange %>% abs,
        aging.fpkm.expressed.q3$log2FoldChange %>% abs,
        aging.fpkm.expressed.q4$log2FoldChange %>% abs,outline=F)
abline(h=median(abs(p5pair.genes.fpkm.yy1$log2FoldChange)))


tapply(aging.fpkm.expressed$aging_mean,
       aging.fpkm.expressed$quantilegroup,
       median)

median(p5pair.genes.fpkm.yy1$aging_fpkm_mean)
median(p5pair.genes.fpkm.NOTyy1$aging_fpkm_mean)
```
# YY1 on enhancer - finalize - just P5
genes targeted by SEs to which is bound in young cells (YY1 binding SEs) - more stable
genes targeted by SEs with no YY1 enrichment in young cells (YY1 non-binding SEs)  - less stable
only use common SE pair

```{r}
whole_step1 <- function(p5pair.genes.fpkm, yy1_p5, output_name){
  YY1_on_se(p5pair.genes.fpkm, yy1_p5) -> p5pair.genes.fpkm.yy1
  p5pair.genes.fpkm[!p5pair.genes.fpkm$pair %in% p5pair.genes.fpkm.yy1$pair,] -> p5pair.genes.fpkm.NOTyy1
  
  #---
  pdf(NULL)
  dev.control(displaylist="enable")
  
  boxplot(abs(p5pair.genes.fpkm.yy1$log2FoldChange), abs(p5pair.genes.fpkm.NOTyy1$log2FoldChange), outline=F, 
        names=c('with YY1', 'without YY1'), main='gene expression fold change')
  
  p.expression <- recordPlot()
  invisible(dev.off())
  output_ppt_m(p.expression, paste0(output_name, '.pptx'))
  #---
  pdf(NULL)
  dev.control(displaylist="enable")
  
  boxplot(abs(p5pair.genes.fpkm.yy1$log2FoldChange),
          abs(p5pair.genes.fpkm$log2FoldChange),
          abs(fc_data.short$log2FoldChange), outline=F,
          names=c('with YY1', 'all SE target', 'all expressed'), main='gene expression fold change')
  
  p.expression <- recordPlot()
  invisible(dev.off())
  output_ppt_m(p.expression, paste0(output_name, '.pptx'))
  
  #---
  pdf(NULL)
  dev.control(displaylist="enable")
  
  boxplot(abs(p5pair.genes.fpkm.yy1$log2FoldChange),
          abs(p5pair.genes.fpkm.NOTyy1$log2FoldChange),
          abs(p5pair.genes.fpkm.yy1$yy1kd_log2FoldChange), 
          abs(fc_data.short$log2FoldChange), outline=F,
          names=c('with YY1', 'without YY1', 'with YY1 YY1KD', 'All expressed'), main='YY1 KD fold change')
  
  p.expression <- recordPlot()
  invisible(dev.off())
  output_ppt_m(p.expression, paste0(output_name, '.pptx'))

  #---
  pdf(NULL)
  dev.control(displaylist="enable")
  
  aging.fpkm.expressed=aging.fpkm[aging.fpkm$aging_mean>1, ]
  boxplot(p5pair.genes.fpkm.yy1$aging_fpkm_mean, p5pair.genes.fpkm$aging_fpkm_mean, aging.fpkm.expressed$aging_mean, 
          outline=F,
          names=c('with YY1', 'without YY1', 'all expressed'), main='gene expression (FPKM)')
  
  p.expression <- recordPlot()
  invisible(dev.off())
  output_ppt_m(p.expression, paste0(output_name, '.pptx'))
  #
  return(list(agine=get_sig_pct(p5pair.genes.fpkm.yy1, p5pair.genes.fpkm.NOTyy1),
              yy1kd=get_sig_pct1(p5pair.genes.fpkm.yy1, p5pair.genes.fpkm.NOTyy1),
              withYY1=p5pair.genes.fpkm.yy1,
              noYY1=p5pair.genes.fpkm.NOTyy1))
}


whole_step(p5pair.genes.fpkm.common, yy1_p5)


# new - 2023-11-23
whole_step2 <- function(p5pair.genes.fpkm, yy1_p5){
  YY1_on_se(p5pair.genes.fpkm, yy1_p5) -> p5pair.genes.fpkm.yy1
  p5pair.genes.fpkm[!p5pair.genes.fpkm$pair %in% p5pair.genes.fpkm.yy1$pair,] -> p5pair.genes.fpkm.NOTyy1
  
  #---
  # pdf(NULL)
  # dev.control(displaylist="enable")
  
  # boxplot(abs(p5pair.genes.fpkm.yy1$log2FoldChange),
  #         abs(p5pair.genes.fpkm.NOTyy1$log2FoldChange),
  #         abs(fc_data.short$log2FoldChange), outline=F,
  #         names=c('with YY1', 'without YY1', 'All expressed'), main='Fig6C')
  # 
  # p.expression <- recordPlot()
  # invisible(dev.off())
  # output_ppt_m(p.expression, paste0(output_name, '.pptx'))
  
  
  plot.dat1 = data.frame(key='with YY1', value=abs(p5pair.genes.fpkm.yy1$log2FoldChange))
  plot.dat2 = data.frame(key='without YY1', value=abs(p5pair.genes.fpkm.NOTyy1$log2FoldChange))
  plot.dat3 = data.frame(key='All expressed', value=abs(fc_data.short$log2FoldChange))

  plot.dat = rbind(plot.dat1, plot.dat2, plot.dat3)
  plot.dat$key = factor(plot.dat$key, levels=c('with YY1',
                                               'without YY1',
                                               'All expressed'))


  for(k in unique(plot.dat$key)){
    boxplot_filter(plot.dat[plot.dat$key==k,]) -> plot.dat[plot.dat$key==k,]
  }

  saveRDS(plot.dat, 'Fig6C.rds')
}


```


# YY1 on enhancer - finalize - merge P5 and P15

genes linked by SE and SE is occupied by YY1 - most stable - highest gene expression
genes linked by SE but not occupied by YY1 - a little bit less stable
all expressed genes - least stable - low gene expression


highly expressed genes naturally tend to be stable - how to diff this.

```{r}

rbind(p5pair.genes.fpkm, p15pair.genes.fpkm) %>%
  dplyr::filter(!duplicated(pair)) -> all_se_pair.fpkm



sum(all_se_pair.fpkm$gene_id %>% duplicated())
all_se_pair.fpkm %>% dplyr::filter(duplicated(gene_id)) %>% dplyr::select(gene_id) -> dup.genes

# remove confusing gene - different isoforme located in different bin, or 2 adjscent SE regulate the same gene
all_se_pair.fpkm[!all_se_pair.fpkm$gene_id %in% dup.genes$gene_id,] -> all_se_pair.fpkm

whole_step1(all_se_pair.fpkm, yy1, 'YY1 stablize gene expression - all SE pair with and without YY1 - RUV4.pptx') -> stabilize.gene.out

whole_step2(all_se_pair.fpkm, yy1)

boxplot(abs(stabilize.gene.out$withYY1$log2FoldChange),
        abs(stabilize.gene.out$noYY1$log2FoldChange),outline=F)

# YY1 vs no YY1
wilcox.test(abs(stabilize.gene.out$withYY1$log2FoldChange),
        abs(stabilize.gene.out$noYY1$log2FoldChange))

# YY1 vs all
wilcox.test(abs(stabilize.gene.out$withYY1$log2FoldChange),
        abs(fc_data.short$log2FoldChange))

# no YY1 vs all
wilcox.test(abs(stabilize.gene.out$noYY1$log2FoldChange),
        abs(fc_data.short$log2FoldChange))


# YY1 vs YY1 KD
wilcox.test(abs(stabilize.gene.out$withYY1$log2FoldChange),
        abs(stabilize.gene.out$withYY1$yy1kd_log2FoldChange))


boxplot(abs(stabilize.gene.out$noYY1$log2FoldChange),
        abs(fc_data.short$log2FoldChange), outline=F)


wilcox.test(abs(stabilize.gene.out$withYY1$log2FoldChange),
        abs(fc_data.short$log2FoldChange))




```

# bu, se target gene vs all expressed gene
```{r}
# all SE target vs all gene
pdf(NULL)
dev.control(displaylist="enable")

aging.fpkm.expressed=aging.fpkm[aging.fpkm$aging_mean>1, ]

boxplot(abs(all_se_pair.fpkm$log2FoldChange),
        abs(fc_data.short$log2FoldChange), outline=F,
         names=c('SE target', 'all expressed'), main='gene expression (FPKM)')

p.expression <- recordPlot()
invisible(dev.off())
output_ppt(p.expression, 'SE target vs all expressed geen expression.pptx')

# SE gene vs all gene
wilcox.test(abs(all_se_pair.fpkm$log2FoldChange),
        abs(fc_data.short$log2FoldChange))
```


# with YY1 in Y, but lose YY1 in O -> should be less stable
```{r}
YY1_on_se(p5pair.genes.fpkm, yy1_p5) -> p5pair.genes.fpkm.yy1
YY1_on_se(p5pair.genes.fpkm.yy1, yy1_p15) -> p5pair.genes.fpkm.yy1.yy1p15


p5pair.genes.fpkm.yy1[!p5pair.genes.fpkm.yy1$pair %in% p5pair.genes.fpkm.yy1.yy1p15$pair, ] -> p5pair.genes.fpkm.yy1.yy1NOTp15


# gain of YY1
p5pair.genes.fpkm[!p5pair.genes.fpkm$pair %in% p5pair.genes.fpkm.yy1$pair,] -> p5pair.genes.fpkm.NOTyy1
YY1_on_se(p5pair.genes.fpkm.NOTyy1, yy1_p15) -> p5pair.genes.fpkm.gain_YY1


boxplot(abs(p5pair.genes.fpkm.yy1.yy1p15$log2FoldChange), abs(p5pair.genes.fpkm.yy1.yy1NOTp15$log2FoldChange),
        abs(p5pair.genes.fpkm.gain_YY1$log2FoldChange), outline=F, 
        names=c('always has YY1', 'lost YY1', 'gain YY1'))

boxplot(p5pair.genes.fpkm.yy1.yy1p15$log2FoldChange, p5pair.genes.fpkm.yy1.yy1NOTp15$log2FoldChange,
        p5pair.genes.fpkm.gain_YY1$log2FoldChange, outline=F, 
        names=c('always has YY1', 'lost YY1', 'gain YY1'))

boxplot(abs(p5pair.genes.fpkm.yy1.yy1p15$log2FoldChange), abs(p5pair.genes.fpkm.yy1.yy1NOTp15$log2FoldChange),
        abs(p5pair.genes.fpkm.yy1.yy1p15$yy1kd_log2FoldChange),
        abs(p5pair.genes.fpkm.yy1.yy1NOTp15$yy1kd_log2FoldChange), outline=F, 
        names=c('always has YY1', 'lost YY1 in O', 'always YY1 group KD','lost YY1 group KD'))

boxplot(abs(p5pair.genes.fpkm.yy1.yy1p15$log2FoldChange), abs(p5pair.genes.fpkm.yy1.yy1NOTp15$log2FoldChange),
        abs(p5pair.genes.fpkm.yy1.yy1NOTp15$yy1kd_log2FoldChange), outline=F, 
        names=c('always has YY1', 'lost YY1 in O', 'lost YY1 group KD'))


boxplot(abs(p5pair.genes.fpkm.yy1.yy1NOTp15$log2FoldChange),
        abs(p5pair.genes.fpkm.yy1.yy1NOTp15$yy1kd_log2FoldChange), outline=F, 
        names=c('lost YY1 in O', 'lost YY1 group KD'))


boxplot(abs(p5pair.genes.fpkm.yy1$log2FoldChange),
        abs(p5pair.genes.fpkm.yy1$yy1kd_log2FoldChange),
        abs(p5pair.genes.fpkm.NOTyy1$log2FoldChange),
        abs(p5pair.genes.fpkm.NOTyy1$yy1kd_log2FoldChange),
        outline=F,
        names=c('with YY1 in Y', 'with YY1 in Y YY1KD', 'wo YY1 in Y', 'wo YY1 in Y YY1KD'))


#---
output_name='aging data - always has YY1 vs lose YY1 in Old - RUV'
pdf(NULL)
dev.control(displaylist="enable")

boxplot(abs(p5pair.genes.fpkm.yy1.yy1p15$log2FoldChange), abs(p5pair.genes.fpkm.yy1.yy1NOTp15$log2FoldChange),
        abs(p5pair.genes.fpkm.gain_YY1$log2FoldChange), outline=F, 
        names=c('always has YY1', 'lost YY1', 'gain YY1'))

p.expression <- recordPlot()
invisible(dev.off())
output_ppt(p.expression, paste0(output_name, '.pptx'))

# new - 20231122
plot.dat1 = data.frame(key='always has YY1', value=abs(p5pair.genes.fpkm.yy1.yy1p15$log2FoldChange))
plot.dat2 = data.frame(key='lost YY1', value=abs(p5pair.genes.fpkm.yy1.yy1NOTp15$log2FoldChange))
plot.dat3 = data.frame(key='gain YY1', value=abs(p5pair.genes.fpkm.gain_YY1$log2FoldChange))

plot.dat = rbind(plot.dat1, plot.dat2, plot.dat3)
plot.dat$key = factor(plot.dat$key, levels=c('always has YY1',
                                             'lost YY1',
                                             'gain YY1'))

for(k in unique(plot.dat$key)){
  boxplot_filter(plot.dat[plot.dat$key==k,]) -> plot.dat[plot.dat$key==k,]
}

saveRDS(plot.dat, 'Fig6D.rds')

p = ggplot(plot.dat, aes(x = key, y = value, fill = key)) +
  geom_jitter(width = 0.2, aes(color = key), size = 1, alpha = 0.8) +  # Add jitter with matching color
  #geom_violin(trim=T, adjust=2, width=0.8, color='black', alpha=0.7) +  # Hide the outliers in the boxplot
  geom_boxplot(outlier.shape = NA, alpha=0.8) +
  scale_colour_brewer(palette = "Set2") +
  scale_fill_brewer(palette = "Set2") +
  theme_classic()


boxplot(abs(p5pair.genes.fpkm.yy1.yy1p15$log2FoldChange),
        abs(p5pair.genes.fpkm.yy1.yy1NOTp15$log2FoldChange), outline=F)

wilcox.test(abs(p5pair.genes.fpkm.yy1.yy1p15$log2FoldChange),
        abs(p5pair.genes.fpkm.yy1.yy1NOTp15$log2FoldChange))

t.test(abs(p5pair.genes.fpkm.yy1.yy1p15$log2FoldChange),
        abs(p5pair.genes.fpkm.yy1.yy1NOTp15$log2FoldChange))

#---
output_name='aging data vs YY1 KD data'
pdf(NULL)
dev.control(displaylist="enable")

boxplot(abs(p5pair.genes.fpkm.yy1.yy1NOTp15$log2FoldChange),
        abs(p5pair.genes.fpkm.yy1.yy1NOTp15$yy1kd_log2FoldChange), outline=F, 
        names=c('lost YY1 in O', 'lost YY1 group KD'))

p.expression <- recordPlot()
invisible(dev.off())
output_ppt(p.expression, paste0(output_name, '.pptx'))


wilcox.test(abs(p5pair.genes.fpkm.yy1.yy1NOTp15$log2FoldChange),
        abs(p5pair.genes.fpkm.yy1.yy1NOTp15$yy1kd_log2FoldChange))
t.test(abs(p5pair.genes.fpkm.yy1.yy1NOTp15$log2FoldChange),
        abs(p5pair.genes.fpkm.yy1.yy1NOTp15$yy1kd_log2FoldChange))

#--

boxplot(abs(p5pair.genes.fpkm.yy1.yy1p15$log2FoldChange), 
        abs(p5pair.genes.fpkm.yy1.yy1p15$yy1kd_log2FoldChange),
        abs(p5pair.genes.fpkm.yy1.yy1NOTp15$log2FoldChange),
        abs(p5pair.genes.fpkm.yy1.yy1NOTp15$yy1kd_log2FoldChange), outline=F, 
        names=c('always has YY1', 'always YY1 group KD', 'lost YY1 in O','lost YY1 group KD'))


boxplot(abs(p5pair.genes.fpkm.yy1.yy1p15$log2FoldChange)/abs(p5pair.genes.fpkm.yy1.yy1p15$yy1kd_log2FoldChange),
        abs(p5pair.genes.fpkm.yy1.yy1NOTp15$log2FoldChange)/abs(p5pair.genes.fpkm.yy1.yy1NOTp15$yy1kd_log2FoldChange), outline=F, 
        names=c('always has YY1 ratio', 'lost YY1 in O ratio'))


wilcox.test(abs(p5pair.genes.fpkm.yy1.yy1NOTp15$log2FoldChange),
        abs(p5pair.genes.fpkm.yy1.yy1NOTp15$yy1kd_log2FoldChange))
t.test(abs(p5pair.genes.fpkm.yy1.yy1NOTp15$log2FoldChange),
        abs(p5pair.genes.fpkm.yy1.yy1NOTp15$yy1kd_log2FoldChange))




wilcox.test(abs(p5pair.genes.fpkm.yy1$log2FoldChange),
        abs(p5pair.genes.fpkm.yy1$yy1kd_log2FoldChange))

wilcox.test(abs(p5pair.genes.fpkm.NOTyy1$log2FoldChange),
        abs(p5pair.genes.fpkm.NOTyy1$yy1kd_log2FoldChange))

boxplot(abs(p5pair.genes.fpkm.yy1$log2FoldChange),
        abs(p5pair.genes.fpkm.yy1$yy1kd_log2FoldChange),
        abs(p5pair.genes.fpkm.NOTyy1$log2FoldChange),
        abs(p5pair.genes.fpkm.NOTyy1$yy1kd_log2FoldChange),
        outline=F,
        names=c('with YY1 in Y', 'with YY1 in Y YY1KD', 'wo YY1 in Y', 'wo YY1 in Y YY1KD'))
# ------------------------------------------------------------------------------------------
# 
# YY1_on_se(p5pair.genes.fpkm, yy1_p5) -> p5pair.genes.fpkm.yy1
# YY1_on_se(p5pair.genes.fpkm.yy1, yy1_p15) -> p5pair.genes.fpkm.yy1.yy1p15
# 
# p5pair.genes.fpkm.yy1[!p5pair.genes.fpkm.yy1$pair %in% p5pair.genes.fpkm.yy1.yy1p15, ] -> p5pair.genes.fpkm.yy1.yy1NOTp15
# 
# boxplot(abs(p5pair.genes.fpkm.yy1.yy1p15$log2FoldChange), abs(p5pair.genes.fpkm.yy1.yy1NOTp15$log2FoldChange), outline=F, 
#         names=c('always has YY1', 'lost YY1 in O'))
# 
# 
# # ------------------------------------------------------------------------------------------
# 
# 
```

# compare between YY1 KD and aging data - Fig5 E
```{r}
boxplot(abs(p5pair.genes.fpkm.yy1.yy1p15$log2FoldChange),
        abs(p5pair.genes.fpkm.yy1.yy1p15$yy1kd_log2FoldChange),
        abs(p5pair.genes.fpkm.yy1.yy1NOTp15$log2FoldChange),
        abs(p5pair.genes.fpkm.yy1.yy1NOTp15$yy1kd_log2FoldChange),
        outline=F,
        names=c('maintain YY1', 'maintain YY1 YY1KD', 'lost YY1', 'lost YY1 YY1KD'))


boxplot(p5pair.genes.fpkm.yy1.yy1p15$yy1kd_log2FoldChange-p5pair.genes.fpkm.yy1.yy1p15$log2FoldChange,
        p5pair.genes.fpkm.yy1.yy1NOTp15$yy1kd_log2FoldChange-p5pair.genes.fpkm.yy1.yy1NOTp15$log2FoldChange,
        outline=F)


# all SE target vs all gene
pdf(NULL)
dev.control(displaylist="enable")

boxplot(abs(p5pair.genes.fpkm.yy1.yy1p15$log2FoldChange),
        abs(p5pair.genes.fpkm.yy1.yy1p15$yy1kd_log2FoldChange),
        abs(p5pair.genes.fpkm.yy1.yy1NOTp15$log2FoldChange),
        abs(p5pair.genes.fpkm.yy1.yy1NOTp15$yy1kd_log2FoldChange),
        outline=F,
        names=c('maintain YY1', 'maintain YY1 YY1KD', 'lost YY1', 'lost YY1 YY1KD'))

p.expression <- recordPlot()
invisible(dev.off())
output_ppt(p.expression, 'compare between YY1 KD and aging data.pptx')

# new - 20231122
plot.dat1 = data.frame(key='maintainYY1', value=abs(p5pair.genes.fpkm.yy1.yy1p15$log2FoldChange))
plot.dat2 = data.frame(key='maintain YY1 YY1KD', value=abs(p5pair.genes.fpkm.yy1.yy1p15$yy1kd_log2FoldChange))
plot.dat3 = data.frame(key='lost YY1', value=abs(p5pair.genes.fpkm.yy1.yy1NOTp15$log2FoldChange))
plot.dat4 = data.frame(key='lost YY1 YY1KD', value=abs(p5pair.genes.fpkm.yy1.yy1NOTp15$yy1kd_log2FoldChange))

plot.dat = rbind(plot.dat1, plot.dat2, plot.dat3, plot.dat4)
plot.dat$key = factor(plot.dat$key, levels=c('P5',
                                             'P15',
                                             'lost YY1',
                                             'lost YY1 YY1KD'))


for(k in unique(plot.dat$key)){
  boxplot_filter(plot.dat[plot.dat$key==k,]) -> plot.dat[plot.dat$key==k,]
}

saveRDS(plot.dat, 'Fig6E.rds')

p = ggplot(plot.dat, aes(x = key, y = value, fill = key)) +
  geom_jitter(width = 0.2, aes(color = key), size = 1.5, alpha = 0.1) +  # Add jitter with matching color
  #geom_violin(trim=T, adjust=2, width=0.8, color='black', alpha=0.7) +  # Hide the outliers in the boxplot
  geom_boxplot(outlier.shape = NA, alpha=0.5) +
  scale_colour_brewer(palette = "Set2") +
  scale_fill_brewer(palette = "Set2") +
  theme_classic()





getwd()
# all SE target vs all gene
pdf(NULL)
dev.control(displaylist="enable")

boxplot(p5pair.genes.fpkm.yy1.yy1p15$yy1kd_log2FoldChange-p5pair.genes.fpkm.yy1.yy1p15$log2FoldChange,
        p5pair.genes.fpkm.yy1.yy1NOTp15$yy1kd_log2FoldChange-p5pair.genes.fpkm.yy1.yy1NOTp15$log2FoldChange,
        outline=F)

p.expression <- recordPlot()
invisible(dev.off())
output_ppt_m(p.expression, 'compare between YY1 KD and aging data.pptx')


# new - 20231124
plot.dat1 = data.frame(key='Maintained YY1', value=p5pair.genes.fpkm.yy1.yy1p15$yy1kd_log2FoldChange-p5pair.genes.fpkm.yy1.yy1p15$log2FoldChange)
plot.dat2 = data.frame(key='Lost YY1', value=p5pair.genes.fpkm.yy1.yy1NOTp15$yy1kd_log2FoldChange-p5pair.genes.fpkm.yy1.yy1NOTp15$log2FoldChange)


plot.dat = rbind(plot.dat1, plot.dat2)
plot.dat$key = factor(plot.dat$key, levels=c('Maintained YY1', 'Lost YY1'))


for(k in unique(plot.dat$key)){
  boxplot_filter(plot.dat[plot.dat$key==k,]) -> plot.dat[plot.dat$key==k,]
}

p = ggplot(plot.dat, aes(x = key, y = value, fill = key)) +
  geom_jitter(width = 0.2, aes(color = key), size = 1.5, alpha = 0.1) +  # Add jitter with matching color
  #geom_violin(trim=T, adjust=2, width=0.8, color='black', alpha=0.7) +  # Hide the outliers in the boxplot
  geom_boxplot(outlier.shape = NA, alpha=0.5) +
  scale_colour_brewer(palette = "Set2") +
  scale_fill_brewer(palette = "Set2") +
  theme_classic()


saveRDS(plot.dat, 'FigS6E.rds')




# SE gene vs all gene
wilcox.test(abs(all_se_pair.fpkm$log2FoldChange),
        abs(fc_data.short$log2FoldChange))



wilcox.test(abs(p5pair.genes.fpkm.yy1.yy1p15$log2FoldChange),
            abs(p5pair.genes.fpkm.yy1.yy1p15$yy1kd_log2FoldChange))

wilcox.test(abs(p5pair.genes.fpkm.yy1.yy1p15$yy1kd_log2FoldChange),
            abs(p5pair.genes.fpkm.yy1.yy1NOTp15$log2FoldChange))

wilcox.test(abs(p5pair.genes.fpkm.yy1.yy1p15$yy1kd_log2FoldChange),
            abs(p5pair.genes.fpkm.yy1.yy1NOTp15$yy1kd_log2FoldChange))

wilcox.test(abs(p5pair.genes.fpkm.yy1.yy1NOTp15$log2FoldChange),
            abs(p5pair.genes.fpkm.yy1.yy1NOTp15$yy1kd_log2FoldChange))

wilcox.test(p5pair.genes.fpkm.yy1.yy1p15$yy1kd_log2FoldChange-p5pair.genes.fpkm.yy1.yy1p15$log2FoldChange,
            p5pair.genes.fpkm.yy1.yy1NOTp15$yy1kd_log2FoldChange-p5pair.genes.fpkm.yy1.yy1NOTp15$log2FoldChange)

sum(p5pair.genes.fpkm.yy1.yy1p15$padj<0.05)
sum(p5pair.genes.fpkm.yy1.yy1p15$yy1kd_padj<0.05)

sum(p5pair.genes.fpkm.yy1.yy1NOTp15$padj<0.05)
sum(p5pair.genes.fpkm.yy1.yy1NOTp15$yy1kd_padj<0.05)
```






# check YY1 on promoter
# with YY1 in Y, but lose YY1 in O -> should be less stable
```{r}



```


# 
```{r}
yy1_p5

yy1_p15



vennplot(yy1_p5, yy1_p15, 'Y', 'O')
```


# all SE pair with hic data

- significantly more hic contact between the pair => gene expression no change
- significantly less hic contact between the pair => gene expression no change

- overall, hic contact between the pair do not correlate with gene expression


```{r}
file='../all_se_pair_withHiCinfo.csv'
hic = fread(file)

#oission test
apply(hic, 1, function(x) {
  p_value = poisson.test(as.integer(x[6]), as.numeric(x[3]))
  p_value$p.value
  }) -> hic$p_values

# apply(hic, 1, function(x) {
#   p_value = matrix(c(as.integer(x[3]), 
#                      as.integer(x[4]),
#                      as.integer(x[6]), 
#                      as.integer(x[7])), 2, 2) %>% fisher.test()
#   p_value$p.value
#   }) -> hic$p_values


hic %>% mutate(q_value=p.adjust(hic$p_values)) -> hic
hic %>% dplyr::filter(q_value<0.05) -> sig

sig$V1 = sapply(strsplit(sig$V1, "\\."),'[[', 2)

sig.up = sig %>% dplyr::filter(p5_readcount<p15_readcount)
sig.dn = sig %>% dplyr::filter(p5_readcount>p15_readcount)


p5pair.genes.fpkm$target_start = sapply(strsplit(p5pair.genes.fpkm$target_region_id, "_"),'[[', 2)
p5pair.genes.fpkm$pair = paste(p5pair.genes.fpkm$se_id, p5pair.genes.fpkm$target_start, sep='_')
p5pair.genes.fpkm[p5pair.genes.fpkm$pair %in% sig.dn$V1,] -> p5pair.genes.fpkm.dn

merge(p5pair.genes.fpkm.dn, fc_data.short, by='gene_id') -> p5pair.genes.fpkm.dn

boxplot(p5pair.genes.fpkm.dn$log2FoldChange.y)
abline(h=0)

p5pair.genes.fpkm.dn %>% dplyr::filter(padj < 0.1)


#-------------------
p15pair.genes.fpkm$target_start = sapply(strsplit(p15pair.genes.fpkm$target_region_id, "_"),'[[', 2)
p15pair.genes.fpkm$pair = paste(p15pair.genes.fpkm$se_id, p15pair.genes.fpkm$target_start, sep='_')

p15pair.genes.fpkm[p15pair.genes.fpkm$pair %in% sig.up$V1,] -> p15pair.genes.fpkm.up
merge(p15pair.genes.fpkm.up, fc_data.short, by='gene_id') -> p15pair.genes.fpkm.up
boxplot(p15pair.genes.fpkm.up$log2FoldChange.y)
abline(h=0)

p15pair.genes.fpkm.up %>% dplyr::filter(padj < 0.1)
```

# repeat my previous analysis
- find pair with YY1 down or up
- identify genes in the target region
- make sure they have less or more HiC contact
- check gene expression
















```{}
cal_fc(P5_pair_withYY1down.fpkm) -> P5_pair_withYY1down.fpkm.fc
# only use genes with FPKM > 1
P5_pair_withYY1down.fpkm.fc %>% dplyr::filter(aging_fpkm_mean > 1) -> P5_pair_withYY1down.fpkm.fc


with(P5_pair_withYY1down.fpkm.fc,
     boxplot(rep1_fc,
             rep2_fc,
             rep3_fc,
             YY1_KD_rep1_fc,
             YY1_KD_rep2_fc,
             YY1_KD_rep3_fc, outline=F))
abline(h=1)


write.csv(P5_pair_withYY1down.fpkm.fc, 'P5_pair_withYY1down.fpkm.fc.csv')


P5_pair_withYY1down.fpkm.fc$target_region_id %>% table

P5_pair_withYY1down.fpkm.fc$fc_mean = (
                                         P5_pair_withYY1down.fpkm.fc$rep2_fc + 
                                         P5_pair_withYY1down.fpkm.fc$rep3_fc)/2

P5_pair_withYY1down.fpkm.fc$YY1_KD_fc_mean = (
                                         P5_pair_withYY1down.fpkm.fc$YY1_KD_rep2_fc+
                                         P5_pair_withYY1down.fpkm.fc$YY1_KD_rep3_fc)/2

boxplot(P5_pair_withYY1down.fpkm.fc$fc_mean,
        P5_pair_withYY1down.fpkm.fc$YY1_KD_fc_mean, outline=F)
abline(h=1)
```


```{r}
add_fpkm(P15_pair_withYY1up, aging.fpkm, yy1kd.fpkm) -> P15_pair_withYY1up.fpkm

hic.dat1 = fread('../P15_pair_withYY1up_withHiCinfo.csv')
hic.dat1.dn = hic.dat1[hic.dat1$P15>hic.dat1$P5,]

P15_pair_withYY1up.fpkm = P15_pair_withYY1up.fpkm[P15_pair_withYY1up.fpkm$pair %in% hic.dat1.dn$V1,]

# only use genes with FPKM > 1
P15_pair_withYY1up.fpkm %>% dplyr::filter(aging_fpkm_mean > 1) -> P15_pair_withYY1up.fpkm
cal_fc(P15_pair_withYY1up.fpkm) -> P15_pair_withYY1up.fpkm.fc
#P15_pair_withYY1up.fpkm.fc$aging_fpkm_mean %>% quantile() -> p15pair.genes.quantile

#P15_pair_withYY1up.fpkm %>% dplyr::filter(aging_fpkm_mean > p15pair.genes.quantile[1]) -> P15_pair_withYY1up.fpkm

with(P15_pair_withYY1up.fpkm.fc,
     boxplot(rep1_fc,
             rep2_fc,
             rep3_fc,
             YY1_KD_rep1_fc,
             YY1_KD_rep2_fc,
             YY1_KD_rep3_fc, outline=F))
abline(h=1)


write.csv(P15_pair_withYY1up.fpkm.fc, 'P15_pair_withYY1up.fpkm.fc.csv')


P15_pair_withYY1up.fpkm.fc$target_region_id %>% table


P15_pair_withYY1up.fpkm.fc$gene_biotype %>% table
```

```{r}
p5pair.genes.fpkm

pdplot <- function(p5pair.genes){
  boxplot(p5pair.genes.fpkm$rep1_old/p5pair.genes.fpkm$rep1_young,
        p5pair.genes.fpkm$rep2_old/p5pair.genes.fpkm$rep2_young,
        p5pair.genes.fpkm$rep3_old/p5pair.genes.fpkm$rep3_young, outline=F)
}



boxplot(p15pair.genes.fpkm$rep1_old/p15pair.genes.fpkm$rep1_young,
        p15pair.genes.fpkm$rep2_old/p15pair.genes.fpkm$rep2_young,
        p15pair.genes.fpkm$rep3_old/p15pair.genes.fpkm$rep3_young, outline=F)


p5pair.genes.fpkm %>% dplyr::filter(grepl('p5', se_id)) -> p5pair.genes.fpkm.sp
p15pair.genes.fpkm %>% dplyr::filter(grepl('p15', se_id)) -> p15pair.genes.fpkm.sp

pdplot(p5pair.genes.fpkm.sp)
pdplot(p15pair.genes.fpkm.sp)



rowMeans(aging.fpkm[, c(2:7)]) > 1 -> keep
aging.fpkm[keep,] -> aging.fpkm.expressed

pdplot(aging.fpkm.expressed)
abline(h = 1)


boxplot(aging.fpkm.expressed$rep1_old/aging.fpkm.expressed$rep1_young,
        p5pair.genes.fpkm$rep1_old/p5pair.genes.fpkm$rep1_young,
        p5pair.genes.fpkm.sp$rep1_old/p5pair.genes.fpkm.sp$rep1_young,
        P5_pair_withYY1down.fpkm.fc.new$rep1_old/P5_pair_withYY1down.fpkm.fc.new$rep1_young,
        P15_pair_withYY1up.fpkm.fc.new$rep1_old/P15_pair_withYY1up.fpkm.fc.new$rep1_young,outline=F)
abline(h=median(aging.fpkm.expressed$rep1_old/aging.fpkm.expressed$rep1_young, na.rm=T))

boxplot(aging.fpkm.expressed$rep2_old/aging.fpkm.expressed$rep2_young,
        p5pair.genes.fpkm$rep2_old/p5pair.genes.fpkm$rep2_young,
        p5pair.genes.fpkm.sp$rep2_old/p5pair.genes.fpkm.sp$rep2_young,
        P5_pair_withYY1down.fpkm.fc.new$rep2_old/P5_pair_withYY1down.fpkm.fc.new$rep2_young,
        P15_pair_withYY1up.fpkm.fc.new$rep2_old/P15_pair_withYY1up.fpkm.fc.new$rep2_young,outline=F)
abline(h=median(aging.fpkm.expressed$rep2_old/aging.fpkm.expressed$rep2_young, na.rm=T))


boxplot(aging.fpkm.expressed$rep3_old/aging.fpkm.expressed$rep3_young,
        p5pair.genes.fpkm$rep3_old/p5pair.genes.fpkm$rep3_young,
        p5pair.genes.fpkm.sp$rep3_old/p5pair.genes.fpkm.sp$rep3_young,
        P5_pair_withYY1down.fpkm.fc.new$rep3_old/P5_pair_withYY1down.fpkm.fc.new$rep3_young,
        P15_pair_withYY1up.fpkm.fc.new$rep3_old/P15_pair_withYY1up.fpkm.fc.new$rep3_young,outline=F)
abline(h=median(aging.fpkm.expressed$rep3_old/aging.fpkm.expressed$rep3_young, na.rm=T))




rbind(P5_pair_withYY1down.fpkm.fc %>% dplyr::filter(rep1_fc>1,
                                              rep2_fc>1,
                                              rep3_fc>1), 

P5_pair_withYY1down.fpkm.fc %>% dplyr::filter(rep1_fc<1,
                                              rep2_fc<1,
                                              rep3_fc<1)) -> P5_pair_withYY1down.fpkm.fc.new


rbind(P15_pair_withYY1up.fpkm.fc %>% dplyr::filter(rep1_fc>1,
                                              rep2_fc>1,
                                              rep3_fc>1), 

P15_pair_withYY1up.fpkm.fc %>% dplyr::filter(rep1_fc<1,
                                              rep2_fc<1,
                                              rep3_fc<1)) -> P15_pair_withYY1up.fpkm.fc.new



```




# absolute changes of genes pairs

# YY1 on SE and target gene expression
```{r}
boxplot(abs(P5_pair_withYY1down.fpkm.fc$rep1_fc), abs(aging.fpkm$rep1_old/aging.fpkm$rep1_young),outline=F)


p5pair.genes.fpkm

YY1_on_se(p5pair.genes.fpkm.high, yy1_p5) -> p5pair.genes.seOnYY1
p5pair.genes.fpkm.high[!p5pair.genes.fpkm.high$se_id %in% p5pair.genes.seOnYY1$se_id,] -> p5pair.genes.seNOTonYY1


boxplot(p5pair.genes.seOnYY1$rep1_young, 
        p5pair.genes.seOnYY1$rep2_young, 
        p5pair.genes.seOnYY1$rep3_young,
        p5pair.genes.seNOTonYY1$rep1_young,
        p5pair.genes.seNOTonYY1$rep2_young,
        p5pair.genes.seNOTonYY1$rep3_young, outline=F)


boxplot(p5pair.genes.seOnYY1$YY1_KD_rep1, 
        p5pair.genes.seOnYY1$YY1_KD_rep2, 
        p5pair.genes.seOnYY1$YY1_KD_rep3,
        p5pair.genes.seNOTonYY1$YY1_KD_rep1,
        p5pair.genes.seNOTonYY1$YY1_KD_rep2,
        p5pair.genes.seNOTonYY1$YY1_KD_rep3, outline=F)

```







































