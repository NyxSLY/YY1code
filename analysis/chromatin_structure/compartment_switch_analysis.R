---
title: "R Notebook"
output: html_notebook
---

```{r message=F}
#session::restore.session('../BASE.RData')
#session::restore.session('../functions.RData')
source('D:\\BaiduNetdiskWorkspace\\OneDrive\\scripts\\R\\commonAPI.R')
source('D:\\BaiduNetdiskWorkspace\\OneDrive\\scripts\\R\\annotation.R')
library(data.table)
library(rtracklayer)
library(GenomicRanges)
library(clusterProfiler)
library(ChIPseeker)
library(ReactomePA)
library(ggplot2)
```


```{r}
aging.deseq2 = read.csv('D:\\BaiduNetdiskWorkspace\\OneDrive\\work\\3dGenome\\Find_enhancer_target\\remove_batch_effect\\DESEQ2_total_results_RUV.csv')
aging.deseq2 = aging.deseq2[!is.na(aging.deseq2$padj),]
colnames(aging.deseq2)[1] = 'gene_id'
```


```{r}
library(EnsDb.Hsapiens.v75)
hg19=EnsDb.Hsapiens.v75
pros = promoters(hg19, upstream = 1000, downstream=1000)
pros = data.frame(pros)
pros = pros[pros$tx_biotype=='protein_coding',]


pros = pros[,c("seqnames", "start", "end", "gene_id")]

merged_pros <- pros %>%
  group_by(gene_id) %>%
  slice(1) %>%
  mutate(start = min(start),
         end = max(end))


dim(merged_pros)
```


```{r}
import('compartmentDiff_A2B.bed') -> a2b
import('compartmentDiff_B2A.bed') -> b2a
merged_pros %>% mutate(seqnames = paste0('chr', seqnames)) -> merged_pros


merge(merged_pros, aging.deseq2, by='gene_id') -> data.expression


data.expression = makeGRangesFromDataFrame(data.expression, keep.extra.columns = T)
```


```{r, warning=F}
data.expression[overlapsAny(data.expression, a2b)] -> genes.a2b
data.expression[overlapsAny(data.expression, b2a)] -> genes.b2a



boxplot(genes.a2b$log2FoldChange, genes.b2a$log2FoldChange, 
        data.expression$log2FoldChange,outline=F, names=c('a2b(22 genes)', 'b2a(only3gene)', 'all data'))


# try with only significant genes.
# BAD - 2 and 0 significant genes found in the datasets.

hist(genes.a2b$log2FoldChange)
hist(genes.b2a$log2FoldChange)

genes.a2b[genes.a2b$padj<0.05,]

genes.b2a[genes.b2a$padj<0.05,]

```























