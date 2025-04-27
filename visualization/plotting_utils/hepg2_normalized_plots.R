# consider both SE-promoter loops from normal tissue and cancer

library(data.table)
library(commonAPI)
library(annotate)
library(dplyr)
library(EnsDb.Hsapiens.v86)
source('D:/BaiduNetdiskWorkspace/OneDrive/scripts/R/annotation.R')

library("rstudioapi")
curr_dir = dirname(getSourceEditorContext()$path)
setwd(curr_dir)


find_se_pro_loop = function(se_file, loop_file){
  se = fread(file.path(se_file))
  se = se[se$isSuper==1,]
  se = makeGRangesFromDataFrame(se, keep.extra.columns = T)

  # load loop
  loop = fread(file.path('contact_loop', loop_file))
  loop$id = paste0('ID', seq(1:nrow(loop)))
  colnames(loop)[1] = 'V1'

  makeGRangesFromDataFrame(loop, seqnames.field = 'V4',
                           start.field = 'V5',
                           end.field = 'V6',
                           keep.extra.columns = T) -> loop.right

  b_gr = toGR(loop)
  overlaps <- findOverlaps(se, b_gr)

  out = data.frame()
  for (i in 1:length(overlaps)) {
    temp <- cbind(data.frame(se[queryHits(overlaps[i])]), data.frame(b_gr[subjectHits(overlaps[i])]))
    out = rbind(out, temp)
  }

  b_gr = loop.right
  overlaps <- findOverlaps(se, b_gr)

  out1 = data.frame()
  for (i in 1:length(overlaps)) {
    temp <- cbind(data.frame(se[queryHits(overlaps[i])]), data.frame(b_gr[subjectHits(overlaps[i])]))
    out1 = rbind(out1, temp)
  }

  colnames(out1) = colnames(out)
  out2 = rbind(out, out1)


  out2$index = paste0(out2$REGION_ID, out2$id)

  # Analysis showed that duplicate indices occur because SE is too large, and both ends of the loop fall within one SE.
  # This situation does not reflect the direct relationship between SE and promoter, it's essentially an internal SE loop, and should be removed entirely.

  out2[duplicated(out2$index),]$index -> index.remove
  out2[!out2$index %in% index.remove,] -> out2
  colnames(out2)[1:3] = paste0('se_', colnames(out2)[1:3])
  out2 %>% dplyr::select("se_seqnames", 'se_start', 'se_end', REGION_ID, seqnames,
                         start, end, V4, V5, V6, id, index) -> out3
  return(out3)
}


# for cancer cell
se_file='HepG2_H3K27ac_SE/hepG2_H3K27ac_AllStitched_REGION_TO_GENE.txt'
loop_file='ENCFF264RQT_HepG2.bedpe'
out3_hepg2 = find_se_pro_loop(se_file, loop_file)

# for normal cell
se_file='HepG2_H3K27ac_normal_SE/hepG2_H3K27ac_AllStitched_REGION_TO_GENE.txt'
loop_file='ENCFF498XDH_normal_liver.bedpe'
out3_normal = find_se_pro_loop(se_file, loop_file)


out3_normal -> out3

target_region = out3 %>% dplyr::select(V4, V5, V6, id)
hg38=EnsDb.Hsapiens.v86
pros = promoters(hg38, upstream = 1000, downstream=1000)
pros = data.frame(pros)
pros$seqnames = as.character(pros$seqnames)
pros$seqnames = paste0('chr',pros$seqnames)
pros = makeGRangesFromDataFrame(pros, keep.extra.columns = T)

get_target_tss <- function(target_region){
  target_region = makeGRangesFromDataFrame(target_region, seqnames.field = 'V4',
                                           start.field = 'V5', end.field = 'V6',
                                           keep.extra.columns = T)
  target_region.df = data.frame(target_region)
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

get_target_tss(target_region) -> se_target_pro
se_target_pro[!duplicated(se_target_pro$gene_id), ] -> se_target_pro


merge(out3, se_target_pro, by.x='id',by.y='region_id') -> out4

colnames(out4) = c('loop_id', 'se_chr','se_start','se_end','se_id','loop_anchor_chr','loop_anchor_start',
                   'loop_anchor_end', 'loop_anchor_se_target_chr', 'loop_anchor_se_target_start', 'loop_anchor_se_target_end',
                   'se_loop_id', 'target_gene_chr','target_gene_start','target_gene_end',"width",
                   "strand", "tx_id", "tx_biotype", "tx_cds_seq_start", "tx_cds_seq_end", "gene_id",
                   "tx_name", "target_region_chr", "target_region_start", "target_region_end")


out5 = out4 %>% dplyr::select('loop_id', 'se_chr','se_start','se_end','se_id','loop_anchor_chr','loop_anchor_start',
                              'loop_anchor_end', 'loop_anchor_se_target_chr', 'loop_anchor_se_target_start', 'loop_anchor_se_target_end',
                              'se_loop_id','target_gene_chr','target_gene_start','target_gene_end',"strand", "tx_biotype","gene_id",
                              "target_region_chr", "target_region_start", "target_region_end")

# load RNA-seq
dir='RNA-seq'
fread(file.path(dir, 'edgeR_nom_HepG2_1.csv')) -> deg1
deg1[!is.na(deg1$FDR),] -> deg1

#deg1[!is.na(deg1$padj),] -> deg1
#
#
# #fread(file.path(dir, 's_p.csv')) -> deg
# #deg[!is.na(deg$padj),] -> deg
#
# library(biomaRt)
# ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# gene_info <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
#                    filters = "ensembl_gene_id",
#                    values = out5$gene_id,
#                    mart = ensembl)
#
# gene_info = gene_info[!duplicated(gene_info$hgnc_symbol),]
#
# merge(out5, gene_info, by.x='gene_id', by.y='ensembl_gene_id') -> out5
#
# # only keep protein coding gene
# out6 = out5[out5$tx_biotype =='protein_coding',]


# filter genes based on average FPKM
data.frame(deg1) -> deg1
deg1[,grepl('FPKM', colnames(deg1))] %>% rowMeans() -> deg1$mean_fpkm
deg1[deg1$mean_fpkm>1,] -> deg2

out7 = merge(out5, deg2, by.x='gene_id', by.y='Row.names')

out7.unique_gene = out7[!duplicated(out7$gene_id),]
out7.unique_gene[!is.na(out7.unique_gene$FDR),] -> out7.unique_gene

sum(out7.unique_gene$FDR<0.05)/nrow(out7.unique_gene)
sum(deg2$FDR<0.05)/nrow(deg2)

write.csv(out7.unique_gene, 'normal_as_base_processed_results.csv')


out7.unique_gene = fread('normal_as_base_processed_results.csv')
out7.unique_gene = data.frame(out7.unique_gene)
out7.unique_gene = out7.unique_gene[,-1]


# continue to test with YY1 info
library(rtracklayer)
yy1.hepg2 = import('YY1/cancer_sample/HepG2_YY1_common_peak.bed', format='BED')

# check promoter along
# load RNA-seq
dir='RNA-seq'
fread(file.path(dir, 'edgeR_nom_HepG2_1.csv')) -> deg1
deg1[!is.na(deg1$FDR),] -> deg1

pros = data.frame(pros)
merge(deg1, pros, by.x='Row.names', by.y='gene_id') -> deg2

makeGRangesFromDataFrame(deg2,
                         keep.extra.columns = T,
                         seqnames.field = 'seqnames',
                         start.field = 'start',
                         end.field = 'end') -> deg2.gr

deg2.gr[overlapsAny(deg2.gr, yy1.hepg2)] -> deg2.gr.YY1onPro
deg2.gr[!overlapsAny(deg2.gr, yy1.hepg2)] -> deg2.gr.YY1NOTonPro

deg2.gr.YY1onPro[!duplicated(deg2.gr.YY1onPro$Row.names)] -> deg2.gr.YY1onPro
deg2.gr.YY1NOTonPro[!duplicated(deg2.gr.YY1NOTonPro$Row.names)] -> deg2.gr.YY1NOTonPro

deg2.gr.YY1onPro = deg2.gr.YY1onPro[deg2.gr.YY1onPro$FDR<1]
deg2.gr.YY1NOTonPro = deg2.gr.YY1NOTonPro[deg2.gr.YY1NOTonPro$FDR<1]

boxplot(abs(deg2.gr.YY1onPro$logFC), abs(deg2.gr.YY1NOTonPro$logFC), outline=F,
        names=c('YY1 bound promoter', 'YY1 non-binding'))

wilcox.test(abs(deg2.gr.YY1onPro$logFC), abs(deg2.gr.YY1NOTonPro$logFC))


# Q3
yy1.normal = fread('YY1/normal_sample/ENCFF304RLZ.bed')
yy1.normal = toGR(yy1.normal)


yy1.common = yy1.normal[overlapsAny(yy1.normal, yy1.hepg2)]

yy1.normal.sp = yy1.normal[!overlapsAny(yy1.normal, yy1.hepg2)]
yy1.hepg2.sp = yy1.hepg2[!overlapsAny(yy1.hepg2, yy1.normal)]


# Actually starting Q3

run_func1 = function(yy1.hepg2){
  makeGRangesFromDataFrame(out7.unique_gene,
                           keep.extra.columns = T,
                           seqnames.field = 'se_chr',
                           start.field = 'se_start',
                           end.field = 'se_end') -> out7.unique_gene.seGR
  out7.unique_gene.seGR = out7.unique_gene.seGR[!is.na(out7.unique_gene.seGR$gene_biotype)]
  out7.unique_gene.seGR = out7.unique_gene.seGR[out7.unique_gene.seGR$FDR<1,]
  out7.unique_gene.seGR[overlapsAny(out7.unique_gene.seGR, yy1.hepg2)] -> data.seboundYY1
  out7.unique_gene.seGR[!overlapsAny(out7.unique_gene.seGR, yy1.hepg2)] -> data.seNOTboundYY1

  data.frame(data.seboundYY1) -> data.seboundYY1
  data.frame(data.seNOTboundYY1) -> data.seNOTboundYY1

  boxplot(abs(data.seboundYY1$logFC), abs(data.seNOTboundYY1$logFC), outline=F, names=c('YY1 bound SE', 'YY1 non-binding SE'))
  wilcox.test(abs(data.seboundYY1$logFC), abs(data.seNOTboundYY1$logFC)) -> wilcox_result
  print(wilcox_result)
  print(c(nrow(data.seboundYY1), nrow(data.seNOTboundYY1)))
  return(list(wilcox_test = wilcox_result, YY1_bound_se = data.seboundYY1, non_binding_se = data.seNOTboundYY1))
}


run_func1(yy1.common) -> yy1.common.res
run_func1(yy1.normal.sp) -> yy1.normal.res
run_func1(yy1.hepg2.sp) -> yy1.hepg2.res


# 2023-11-29

data1 = data.frame(key = 'YY1 bound SE', value=abs(yy1.common.res$YY1_bound_se$logFC))
data2 = data.frame(key = 'YY1 non-binding SE', value=abs(yy1.common.res$non_binding_se$logFC))

plot.dat = rbind(data1, data2)
plot.dat$key = factor(plot.dat$key, levels=c('YY1 bound SE',
                                             'YY1 non-binding SE'))


for(k in unique(plot.dat$key)){
  boxplot_filter(plot.dat[plot.dat$key==k,]) -> plot.dat[plot.dat$key==k,]
}

saveRDS(plot.dat, 'HepG2 common FigS5F.rds')

p = ggplot(plot.dat, aes(x = key, y = value, fill = key)) +
  geom_jitter(width = 0.2, aes(color = key), size = 1, alpha = 0.8) +  # Add jitter with matching color
  #geom_violin(trim=T, adjust=2, width=0.8, color='black', alpha=0.7) +  # Hide the outliers in the boxplot
  geom_boxplot(outlier.shape = NA, alpha=0.8) +
  scale_colour_brewer(palette = "Set2") +
  scale_fill_brewer(palette = "Set2") +
  theme_classic()

p
# promoter/SE OR strategy








