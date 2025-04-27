---
title: "R Notebook"
output: html_notebook
---

```{r message=F}
session::restore.session('../BASE.RData')
session::restore.session('../functions.RData')
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

# Aging RNA-seq
```{r}
aging.deseq2 = read.csv('../remove_batch_effect/DESEQ2_total_results_RUV.csv')
aging.deseq2 = aging.deseq2[!is.na(aging.deseq2$padj),]
colnames(aging.deseq2)[1] = 'gene_id'

deseq2_scatter_plot(aging.deseq2, 0.5, 0.5, TRUE) + ylim(-5,5) -> aging.scater.plot
output_ppt(aging.scater.plot, 'Aging RNA-seq scatter plot.pptx')


# just plot axies
deseq2_scatter_plot_empty(aging.deseq2, 0.5, 0.5) + ylim(-5,5) -> axis_plot
output_ppt_m(axis_plot, 'Aging RNA-seq scatter plot.pptx')


aging.deseq2 %>% dplyr::filter(padj<0.05, log2FoldChange < 0) %>% nrow()
aging.deseq2 %>% dplyr::filter(padj<0.05, log2FoldChange < 0) %>% nrow()/nrow(aging.deseq2)


aging.deseq2 %>% dplyr::filter(padj<0.05, log2FoldChange > 0) %>% nrow()
aging.deseq2 %>% dplyr::filter(padj<0.05, log2FoldChange > 0) %>% nrow()/nrow(aging.deseq2)
```
# pathway and go-term
```{r}


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


aging.deseq2 %>% dplyr::filter(padj<0.05) -> sig 
getGeneAnnoEn(sig$gene_id) %>% .$entrez_id -> enterz_id
pw = enrichment_analysis(enterz_id, db, species)


name='Aging degs'

write.csv(data.frame(data.frame(pw$pw)), paste0(name, ' reactome pathway.csv'))
write.csv(data.frame(data.frame(pw$go)), paste0(name, ' GO enrichment pathway.csv'))


aging.deseq2 %>% dplyr::filter(padj<0.05) -> sig 


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




# sig genes overlap with SE
```{r}
# Identify genes in the target region
p5pair = p5.pair
p15pair = p15.pair

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

# add pair id
p5pair.genes$target_start = sapply(strsplit(p5pair.genes$target_region_id, "_"),'[[', 2)
p5pair.genes$pair = paste(p5pair.genes$se_id, p5pair.genes$target_start, sep='_')


p15pair.genes$target_start = sapply(strsplit(p15pair.genes$target_region_id, "_"),'[[', 2)
p15pair.genes$pair = paste(p15pair.genes$se_id, p15pair.genes$target_start, sep='_')


# add fold change and pvalue
p5pair.genes %>% add_fc(., fc_data.short) -> p5pair.genes.fc
p15pair.genes %>% add_fc(., fc_data.short) -> p15pair.genes.fc


p5pair.genes.fc %>% add_fc(., yy1kd_fc_data.short) -> p5pair.genes.fc
p15pair.genes.fc %>% add_fc(., yy1kd_fc_data.short) -> p15pair.genes.fc

rbind(p5pair.genes.fc, p15pair.genes.fc) %>%
  dplyr::filter(!duplicated(pair)) -> all_se_pair.fc

all_se_pair.fc.uniq = all_se_pair.fc[!duplicated(all_se_pair.fc$gene_id),]



# start
aging.deseq2 %>% dplyr::filter(padj<0.05) -> aging.sig
aging.sig %in% all_se_pair.fc.uniq$gene_id %>% sum()


all_se_pair.fc.uniq$gene_id %in% aging.sig$gene_id %>% sum()

all_se_pair.fc.uniq$gene_id %in% aging.sig$gene_id %>% all_se_pair.fc.uniq[.,] -> all_se_pair.fc.uniq.sig

# Non SE target DEGs

sig.genes$se_target = 'N'
sig.genes[sig.genes$gene_id %in% all_se_pair.fc.uniq$gene_id,]$se_target = 'Y'
sig.genes%>% dplyr::select(gene_id, baseMean,log2FoldChange, padj, se_target) -> new.dat

output_ppt_m(boxplot(abs(new.dat[new.dat$se_target=='Y',]$log2FoldChange),
  abs(new.dat[new.dat$se_target=='N',]$log2FoldChange),
        outline=F),
        'boxplot.pptx')

boxplot(abs(new.dat[new.dat$se_target=='Y',]$log2FoldChange),abs(new.dat[new.dat$se_target=='N',]$log2FoldChange),outline=F)

wilcox.test(abs(new.dat[new.dat$se_target=='N',]$log2FoldChange),
            abs(new.dat[new.dat$se_target=='Y',]$log2FoldChange),
            paired=F)

# TRY violin plot
se_target = abs(new.dat[new.dat$se_target=='Y',]$log2FoldChange)
non_se_target = abs(new.dat[new.dat$se_target=='N',]$log2FoldChange)

plotv.dat <- rbind(data.frame(value = se_target, label = rep("se_target", length(se_target))),
                data.frame(value = non_se_target, label = rep("non_se_target", length(non_se_target))))

boxplot_filter(plotv.dat) -> plotv.dat

ggplot(plotv.dat, aes(x = label, y = value)) +
  geom_violin(fill = "lightblue", color = "black", outlier.shape=NA) + 
  theme_classic() + ylim(0,5) -> pp; pp


make_plot(plotv.dat)









deseq2_scatter_plot(all_se_pair.fc.uniq, 1.8, 1) -> se_genes.aging.scatterplot
se_genes.aging.scatterplot

output_ppt(se_genes.aging.scatterplot,
           'aging RNA-seq of SE target genes.pptx')


deseq2_scatter_plot_empty(all_se_pair.fc.uniq, 1.8, 1) -> axis.plot
output_ppt_m(axis.plot,
           'aging RNA-seq of SE target genes.pptx')



all_se_pair.fc.uniq$gene_id %in% aging.sig$gene_id %>% sum()/nrow(all_se_pair.fc.uniq)

nrow(aging.sig)/nrow(aging.deseq2)


all_se_pair.fc.uniq %>% dplyr::filter(padj < 0.05, log2FoldChange<0) %>% nrow()
all_se_pair.fc.uniq %>% dplyr::filter(padj < 0.05, log2FoldChange<0) %>% nrow()/nrow(all_se_pair.fc.uniq)
all_se_pair.fc.uniq %>% dplyr::filter(padj < 0.05, log2FoldChange>0) %>% nrow()
all_se_pair.fc.uniq %>% dplyr::filter(padj < 0.05, log2FoldChange>0) %>% nrow()/nrow(all_se_pair.fc.uniq)


all_se_pair.fc.uniq %>% dplyr::filter(padj < 0.05) %>% nrow() -> a
nrow(all_se_pair.fc.uniq) - a -> b
aging.deseq2 %>% dplyr::filter(padj<0.05) %>% nrow() -> c
nrow(aging.deseq2) - c -> d

fisher.test(matrix(c(a,
                     b,
                     c,
                     d),2,2)) %>% {.$p.value}

# fisher.test(matrix(c(all_se_pair.fc.uniq %>% dplyr::filter(padj < 0.05) %>% nrow(),
#                      nrow(all_se_pair.fc.uniq),
#                      aging.deseq2 %>% dplyr::filter(padj<0.05) %>% nrow(),
#                      nrow(aging.deseq2)),2,2)) %>% {.$p.value}

fisher.pvalue$p.value




getGeneAnnoEn(all_se_pair.fc.uniq.sig$gene_id) %>% .$entrez_id -> enterz_id
pw = enrichment_analysis(enterz_id, db, species)


name='SE targeted Aging degs'

write.csv(data.frame(data.frame(pw$pw)), paste0(name, ' reactome pathway.csv'))
write.csv(data.frame(data.frame(pw$go)), paste0(name, ' GO enrichment pathway.csv'))
```


# GO-term enrichment
```{r}
# group and classification is done using revigo, with enriched go-term with pvalue as input

files = list.files('./GO-term/')
files %>% grepl('Revigo-.*csv', .) %>% files[.] %>% 
  lapply(function(x) fread(file.path('./GO-term',x))) -> main.go


lapply(main.go, function(x) dplyr::filter(x, Eliminated==FALSE)) -> main.go


files %>% grepl('Revigo-.*csv', .) %>% files[.]
names(main.go) = gsub('.csv', '', files %>% grepl('Revigo-.*csv', .) %>% files[.])

lapply(main.go, function(x) dplyr::select(x, TermID, Name, Value)) -> main.go.short

main.go.short$`Revigo-DEGs-with-pvalue` %>% 
  dplyr::mutate(alsoInSE=if_else(TermID %in% main.go.short$`Revigo-DEGs-of-SE-with-pvalue`$TermID,'Yes','No')) %>%
  arrange(Value) -> aging.degs.go


main.go.short$`Revigo-DEGs-of-SE-with-pvalue` %>% 
  dplyr::mutate(alsoInSE=if_else(TermID %in% main.go.short$`Revigo-DEGs-with-pvalue`$TermID,'Yes','No')) %>%
  arrange(Value) -> se.degs.go


# import group info
files = list.files('./GO-term/')
files %>% grepl('TreeMap.*csv', .) %>% files[.] %>% 
  lapply(function(x) fread(file.path('./GO-term',x),skip=4)) -> treemap.go

lapply(treemap.go, function(x){
  x %>% mutate(Representative=if_else(Representative=='null', Name, Representative))
}) %>% lapply(., function(x) dplyr::select(x,TermID, Representative)) -> treemap.go.short

names(treemap.go.short) = gsub('.csv', '', files %>% grepl('Revigo-.*csv', .) %>% files[.])

aging.degs.go %>% merge(., treemap.go.short$`Revigo-DEGs-with-pvalue`) -> aging.degs.go.group
se.degs.go %>% merge(., treemap.go.short$`Revigo-DEGs-of-SE-with-pvalue`) -> se.degs.go.group


aging.degs.go.group %>% arrange(Value)
tapply(aging.degs.go.group$Value, aging.degs.go.group$Representative, min) %>% sort() %>% names() -> levels
aging.degs.go.group %>% mutate(Representative=factor(Representative, levels=levels)) %>%
  arrange(Representative, Value) -> aging.degs.go.group.sorted


se.degs.go.group %>% arrange(Value)
tapply(se.degs.go.group$Value, se.degs.go.group$Representative, min) %>% sort() %>% names() -> levels
se.degs.go.group %>% mutate(Representative=factor(Representative, levels=levels)) %>%
  arrange(Representative, Value) -> se.degs.go.group.sorted


aging.degs.go.group.sorted %>% mutate(Name=gsub(',', '-', Name)) %>% 
  write.csv('./GO-term/aging_degs_GO-term_processed.csv', row.names = F,
                                         quote = F)

se.degs.go.group.sorted %>% mutate(Name=gsub(',', '-', Name)) %>% 
  write.csv('./GO-term/se_degs_GO-term_processed.csv', row.names = F,
                                         quote = F)


```


# after manual check and classification
```{r}
list.files('./GO-term') %>% grepl('.txt',.) %>% list.files('./GO-term')[.] %>%
  lapply(function(x) fread(file.path('./GO-term',x))) -> final.go

names(final.go) = gsub('.csv', '', list.files('./GO-term') %>% grepl('.txt',.) %>% list.files('./GO-term')[.])

list.files('.') %>% grepl('GO .*.csv',.) %>% list.files('.')[.] %>%
  lapply(function(x) fread(file.path('.',x))) %>% 
  lapply(function(x) dplyr::select(x,ID, pvalue, GeneRatio))-> initial.go


names(initial.go) = gsub('.csv', '', list.files('.') %>% grepl('GO .*.csv',.)  %>% list.files('.')[.])

rbind(final.go$`final-manual checked go-term - aging.txt` %>% dplyr::select(TermID, Name),
      final.go$`final-manual checked go-term - se.txt` %>% dplyr::select(TermID, Name)) -> base
base[!duplicated(base$TermID),] -> base


merge(base, initial.go$`Aging degs GO enrichment pathway`, by.x='TermID', by.y='ID', all.x=TRUE) -> base
merge(base, initial.go$`SE targeted Aging degs GO enrichment pathway`, by.x='TermID', by.y='ID', all.x=TRUE) -> merged.go


merged.go %>% write.csv('go-term plot data.csv')

# manual check and set up order
go_order= read.table('./GO-term/final_dot_plot_manual_factor_order.txt', sep = '\t')

merged.go.final = merged.go[merged.go$Name %in% go_order$V1,]
merged.go.final[merged.go.final$TermID!='GO:2000377',] -> merged.go.final

merged.go.final %>% dplyr::select(TermID, Name, pvalue.x, GeneRatio.x) -> merged.go.aging
names(merged.go.aging) = c('TermID', 'Name', 'Value', 'GeneRatio')

merged.go.final %>% dplyr::select(TermID, Name, pvalue.y, GeneRatio.y) -> merged.go.se
names(merged.go.se) = c('TermID', 'Name', 'Value', 'GeneRatio')
merged.go.aging$group='aging'
merged.go.se$group='se'

rbind(merged.go.aging, merged.go.se) -> plot.dat


plot.dat$GeneRatio %>% str_split('/') %>% sapply('[[', 1) %>% as.numeric() -> gene_number
dd = c(rep(1057, nrow(merged.go.aging)),
         rep(90, nrow(merged.go.se)))


plot.dat$GeneRatio=gene_number/dd


plot.dat$Value=-log2(plot.dat$Value)

plot.dat[is.na(plot.dat$Value),]$Value=0
plot.dat[is.na(plot.dat$GeneRatio),]$GeneRatio=0


```



```{r, fig.height=10, fig.width=15}



plot.dat[plot.dat$Value>15,]$Value=15
plot.dat$group = factor(plot.dat$group, levels = c('se', 'aging'))
plot.dat$Name = factor(plot.dat$Name, levels = go_order$V1)


spot.theme <- list(
  theme_classic(),
  theme(axis.ticks.x=element_blank(), axis.text.x=element_text(size = 19, angle = 90, hjust = 0)),
  theme(axis.ticks.y=element_blank(), axis.text.y=element_text(size = 19)),
  theme(axis.line=element_blank()),
  theme(text = element_text(size = 22)),
  # theme(legend.position = "none"),
  theme(plot.margin = unit(c(10,10,10,10), "mm")),
  scale_size_continuous(range = c(-0.3, 15)),
  scale_x_discrete(position = "top"))


ggplot(data=plot.dat, aes(Name, group)) + geom_point(aes(color = Value, size= GeneRatio)) + spot.theme +
  scale_color_gradient(low="white", high="red")

pdf('GO-term enrichment dot.pdf', height=8, width=15)
ggplot(data=plot.dat, aes(Name, group)) + geom_point(aes(color = Value, size= GeneRatio)) + spot.theme +
  scale_color_gradient(low="white", high="red")
dev.off()

go.plot <- ggplot(data=plot.dat, aes(Name, group)) + geom_point(aes(color = Value, size= GeneRatio)) + spot.theme +
  scale_color_gradient(low="white", high="red")


output_ppt(go.plot, 'GO-term enrichment dot.pptx')

```



# se DGEs hic contact - no failed
```{r}
wd='../total_normalized'
dirs = list.dirs(file.path(wd))
p5_dirs = dirs[grepl('P5', dirs)]
p15_dirs = dirs[grepl('P15', dirs)]


rbind(p5.pair %>% dplyr::select(pair, target_region_end, info),
      p15.pair %>% dplyr::select(pair, target_region_end, info)) -> pair.all
pair.all[!duplicated(pair.all$pair),] -> pair.all
merge(all_se_pair.fc.uniq, pair.all, by='pair') -> all_se_pair.fc.uniq.hicinput

result.p5.Y = run_hic(p5_dirs, all_se_pair.fc.uniq.hicinput)
result.p5.O = run_hic(p15_dirs, all_se_pair.fc.uniq.hicinput)

lapply(result.p5.Y, unlist) %>% unlist %>% data.frame() -> result.p5.Y.n
lapply(result.p5.O, unlist) %>% unlist %>% data.frame() -> result.p5.O.n

cbind(result.p5.Y.n, result.p5.O.n) -> result
names(result) = c('P5', 'P15')

boxplot((result$P5),(result$P15), outline=F)


boxplot((result$P15-result$P5), outline=F)
abline(h=0)




```



# overlap between YY1-KD and Aging 
```{r}
yy1kd_fc_data.short %>% dplyr::filter(yy1kd_padj < 0.05) -> yy1.sig

vennplot_simple(nrow(aging.sig),
                nrow(yy1.sig),
                sum(aging.sig$gene_id %in% yy1.sig$gene_id),
                'Aging',
                'YY1KD')


one_step <- function(com_data.sig, com_data2.sig, com_data, com_data2, name1, name2){
  vennplot_simple(nrow(com_data.sig),
                  nrow(com_data2.sig),
                  com_data.sig$gene_id %in% com_data2.sig$gene_id %>% sum %>% sum,
                  name1, name2)
  
  Total=c(com_data$gene_id, com_data2$gene_id) %>% unique() %>% length()
  Overlap = com_data.sig$gene_id %in% com_data2.sig$gene_id %>% sum
  group1=nrow(com_data.sig)
  group2=nrow(com_data2.sig)
  phyper(Overlap-1, group2, Total-group2, group1,lower.tail= FALSE)
}


one_step(aging.sig, yy1.sig, fc_data.short, yy1kd_fc_data.short, 'aging', 'YY1KD')

output_ppt(one_step(aging.sig, yy1.sig, fc_data.short, yy1kd_fc_data.short, 'aging', 'YY1KD'),
           name='aging DEGs overlaps with YY1 KD DEGs.pptx')

  

yy1kd.up = yy1.sig %>% dplyr::filter(yy1kd_log2FoldChange>0)
yy1kd.dn = yy1.sig %>% dplyr::filter(yy1kd_log2FoldChange<0)

vennplot_simple(nrow(aging.sig),
                nrow(yy1kd.up),
                sum(aging.sig$gene_id %in% yy1kd.up$gene_id),
                'Aging',
                'YY1KD')

vennplot_simple(nrow(aging.sig),
                nrow(yy1kd.dn),
                sum(aging.sig$gene_id %in% yy1kd.dn$gene_id),
                'Aging',
                'YY1KD')



# ======-----------------------
all_se_pair.fc.uniq %>% dplyr::filter(padj < 0.05, log2FoldChange<0) -> all_se_pair.fc.sig.dn
all_se_pair.fc.uniq %>% dplyr::filter(padj < 0.05, log2FoldChange>0) -> all_se_pair.fc.sig.up
all_se_pair.fc.uniq %>% dplyr::filter(padj < 0.05) -> all_se_pair.fc.sig

sum(all_se_pair.fc.sig$gene_id %in% yy1.sig$gene_id)


vennplot_simple(nrow(all_se_pair.fc.sig),
                nrow(yy1kd_fc_data.short),
                sum(all_se_pair.fc.sig$gene_id %in% yy1.sig$gene_id),
                'Aging SE',
                'YY1KD')



#================================
  

```

```{r}
library('pheatmap')
library(RColorBrewer)


aging.sig
yy1kd_fc_data.short
merge(aging.sig, yy1kd_fc_data.short, by='gene_id') -> aging.sig.yy1kd


aging.sig.yy1kd %>% dplyr::select(log2FoldChange, yy1kd_log2FoldChange) -> aging.sig.fc.merged
as.matrix(aging.sig.fc.merged) -> plot.data
row.names(plot.data) = aging.sig.yy1kd$gene_id

my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 1000)

aging.sig.fc.merged %>% arrange(-log2FoldChange) -> aging.sig.fc.merged

pheatmap(aging.sig.fc.merged, cluster_rows=F, show_rownames=F,
         cluster_cols=F, scale='none', treeheight_row=0, border_color = 'gray',
         cellheight=2)



plot(aging.sig.fc.merged$log2FoldChange,aging.sig.fc.merged$yy1kd_log2FoldChange, xlim=c(-4,4))
abline(lm(aging.sig.fc.merged$yy1kd_log2FoldChange~aging.sig.fc.merged$log2FoldChange))


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


aging.sig.fc.merged %>% mutate(group=if_else(log2FoldChange>0 & yy1kd_log2FoldChange >0, 'red', 
                                             if_else(log2FoldChange < 0 & yy1kd_log2FoldChange <0, 'blue', 'gray'))) -> aging.sig.fc.merged.scaterplot

ggplot(aging.sig.fc.merged.scaterplot, aes(x=log2FoldChange, y=yy1kd_log2FoldChange)) + geom_point(aes(color=group, size=0.5)) + theme_classic() +
  geom_smooth(method=lm, se=FALSE) +
  scale_color_manual(values=c('blue', 'gray','red')) + 
  scale_size(range = c(0.7, 0.7)) +
  xlim(-4,4) + ylim(-4,4) + 
  geom_hline(yintercept=0, linetype="dashed", 
                color = "black", size=0.5) +
  geom_vline(xintercept=0, linetype="dashed", 
                color = "black", size=0.5)

overlap.scatter <- ggplot(aging.sig.fc.merged.scaterplot, aes(x=log2FoldChange, y=yy1kd_log2FoldChange)) + geom_point(aes(color=group, size=0.5)) + theme_classic() +
  geom_smooth(method=lm, se=FALSE) +
  scale_color_manual(values=c('blue', 'gray','red')) + 
  scale_size(range = c(0.7, 0.7)) +
  xlim(-4,4) + ylim(-4,4) + 
  geom_hline(yintercept=0, linetype="dashed", 
                color = "black", size=0.5) +
  geom_vline(xintercept=0, linetype="dashed", 
                color = "black", size=0.5)


cor(aging.sig.fc.merged.scaterplot$log2FoldChange,aging.sig.fc.merged.scaterplot$yy1kd_log2FoldChange)


sum(aging.sig.fc.merged.scaterplot$group=='blue')/sum(aging.sig.fc.merged.scaterplot$log2FoldChange<0)
sum(aging.sig.fc.merged.scaterplot$group=='red')/sum(aging.sig.fc.merged.scaterplot$log2FoldChange>0)


output_ppt_m(overlap.scatter, 'aging DEGs overlaps with YY1 KD DEGs.pptx')


#---------------------------------------

aging.sig.se_target = aging.sig[aging.sig$gene_id %in% all_se_pair.fc.uniq$gene_id,]
yy1.sig.se_target = yy1.sig[yy1.sig$gene_id %in% all_se_pair.fc.uniq$gene_id, ]

vennplot_simple(nrow(aging.sig.se_target),
                nrow(yy1.sig.se_target),
                sum(aging.sig.se_target$gene_id %in% yy1.sig.se_target$gene_id),
                'Aging',
                'YY1KD')

output_ppt_m(vennplot_simple(nrow(aging.sig.se_target),
                nrow(yy1.sig.se_target),
                sum(aging.sig.se_target$gene_id %in% yy1.sig.se_target$gene_id),
                'Aging',
                'YY1KD'),
             'aging DEGs overlaps with YY1 KD DEGs.pptx')


all_se_pair.fc.uniq

aging.sig.yy1kd %>% dplyr::filter(gene_id %in% all_se_pair.fc.uniq$gene_id) %>% dplyr::select(log2FoldChange, yy1kd_log2FoldChange) ->
  aging.sig_se_target_gene.fc.merged


aging.sig_se_target_gene.fc.merged %>% arrange(-log2FoldChange) -> aging.sig_se_target_gene.fc.merged

aging.sig_se_target_gene.fc.merged %>% mutate(group=if_else(log2FoldChange>0 & yy1kd_log2FoldChange >0, 'red', 
                                             if_else(log2FoldChange < 0 & yy1kd_log2FoldChange <0, 'blue', 'gray'))) -> aging.sig_se_target_gene.fc.merged.scaterplot


overlap.scatter <- ggplot(aging.sig_se_target_gene.fc.merged.scaterplot, aes(x=log2FoldChange, y=yy1kd_log2FoldChange)) + geom_point(aes(color=group, size=0.5)) + theme_classic() +
  geom_smooth(method=lm, se=FALSE) +
  scale_color_manual(values=c('blue', 'gray','red')) + 
  scale_size(range = c(0.7, 0.7)) +
  xlim(-4,4) + ylim(-4,4) + 
  geom_hline(yintercept=0, linetype="dashed", 
                color = "black", size=0.5) +
  geom_vline(xintercept=0, linetype="dashed", 
                color = "black", size=0.5)




cor(aging.sig_se_target_gene.fc.merged.scaterplot$log2FoldChange,aging.sig_se_target_gene.fc.merged.scaterplot$yy1kd_log2FoldChange)


sum(aging.sig_se_target_gene.fc.merged.scaterplot$group=='blue')/sum(aging.sig_se_target_gene.fc.merged.scaterplot$log2FoldChange<0)
sum(aging.sig_se_target_gene.fc.merged.scaterplot$group=='red')/sum(aging.sig_se_target_gene.fc.merged.scaterplot$log2FoldChange>0)


output_ppt_m(overlap.scatter, 'aging DEGs overlaps with YY1 KD DEGs.pptx')

#-------------------------------------------------
{nrow(aging.sig.yy1kd %>% dplyr::filter(log2FoldChange < 0) %>% dplyr::filter(yy1kd_log2FoldChange < 0))} / {nrow(aging.sig.yy1kd %>% dplyr::filter(log2FoldChange < 0))}
{nrow(aging.sig.yy1kd %>% dplyr::filter(log2FoldChange > 0) %>% dplyr::filter(yy1kd_log2FoldChange > 0))} / {nrow(aging.sig.yy1kd %>% dplyr::filter(log2FoldChange > 0))}

fpkm %>% dplyr::select(gene_id, p5_mean, p15_mean) -> fpkm.mean
yy1.fpkm %>% dplyr::select(gene_id, YY1_KD_mean, YY1_NT_mean) -> yy1.fpkm.mean

merge(fpkm.mean, yy1.fpkm.mean, by='gene_id') -> fpkm.mean.merged

fpkm.mean.merged[fpkm.mean.merged$gene_id %in% aging.sig$gene_id][,c(2:5)] %>% as.matrix -> plot.data


getwd()
```














down-regulated genes - YY1 changes
```{r}
all_se_pair.fc.uniq %>% dplyr::filter(padj < 0.05, log2FoldChange<0) -> all_se_pair.fc.sig.dn
all_se_pair.fc.uniq %>% dplyr::filter(padj < 0.05, log2FoldChange>0) -> all_se_pair.fc.sig.up
all_se_pair.fc.uniq %>% dplyr::filter(padj < 0.05) -> all_se_pair.fc.sig

# these genes have more FC changes in aging data than in YY1 KD data.......
boxplot(all_se_pair.fc.sig.dn$log2FoldChange, all_se_pair.fc.sig.dn$yy1kd_log2FoldChange)

# generate region file

# promoter
library(rtracklayer)

generante_region <- function(all_se_pair.fc.sig.dn, name){
  all_se_pair.fc.sig.dn %>% makeGRangesFromDataFrame() %>% export(paste0(name,'_pro.bed'), format='BED')
  # super enhancer
  all_se_pair.fc.sig.dn %>% dplyr::select(se_chr, se_start, se_end) %>%
    makeGRangesFromDataFrame() -> se_region.gr
  enh_mid[overlapsAny(enh_mid, se_region.gr)] %>% export(paste0(name,'_SE_insideEnh.bed'), format='BED')
  #
  se_region.gr %>% export(paste0(name,'_SE.bed'), format='BED')
}

generante_region(all_se_pair.fc.sig.dn, 'se_target_gene_sigDN')
generante_region(all_se_pair.fc.sig.up, 'se_target_gene_sigUP')
generante_region(all_se_pair.fc.sig, 'se_target_gene_SIG')
generante_region(all_se_pair.fc.uniq, 'All_se_target_gene')


quantile(abs(all_se_pair.fc.uniq$log2FoldChange))

all_se_pair.fc.uniq %>% dplyr::filter(padj > 0.05, abs(log2FoldChange)<0.1) -> all_se_pair.fc.nochange
generante_region(all_se_pair.fc.nochange, 'se_target_gene_noChange')


# for all se target genes (scatter plot)
all_se_pair.fc.uniq %>% makeGRangesFromDataFrame() %>% export('All_se_target_gene_pro.bed', format='BED')
all_se_pair.fc.uniq %>% dplyr::select(se_chr, se_start, se_end) %>% makeGRangesFromDataFrame() %>% export('All_se_target_gene_SE.bed', format='BED')


# You need


enh_mid %>% export('All_enhancer.bed', format='BED')




data.frame(enh_mid) %>% View()



all_se_pair.fc.sig %>% makeGRangesFromDataFrame() %>%
  overlapsAny(yy1) %>% sum() / nrow(all_se_pair.fc.sig)

all_se_pair.fc.nochange %>% makeGRangesFromDataFrame() %>%
  overlapsAny(yy1) %>% sum() / nrow(all_se_pair.fc.nochange)


all_se_pair.fc.sig %>% makeGRangesFromDataFrame(seqnames.field = 'se_chr',
                                                 start.field = 'se_start',
                                                 end.field = 'se_end') %>%
  overlapsAny(yy1) %>% sum() / nrow(all_se_pair.fc.sig)



all_se_pair.fc.nochange %>% makeGRangesFromDataFrame(seqnames.field = 'se_chr',
                                                 start.field = 'se_start',
                                                 end.field = 'se_end') %>%
  overlapsAny(yy1) %>% sum() / nrow(all_se_pair.fc.nochange)

```







































