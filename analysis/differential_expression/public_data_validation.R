library(data.table)
library(commonAPI)
library(annotate)
library(dplyr)
source('D:/BaiduNetdiskWorkspace/OneDrive/scripts/R/annotation.R')
setwd('D:\\BaiduNetdiskWorkspace\\OneDrive\\work\\3dGenome\\manuscript\\writing\\Weiwei_comments_based_on_GSA\\public_datasets\\IMR90')

# Set the working directory to the current script's location
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# load RNA-seq
dir='D:\\BaiduNetdiskWorkspace\\OneDrive\\work\\3dGenome\\manuscript\\writing\\Weiwei_comments_based_on_GSA\\public_datasets\\IMR90\\RNA-seq'
fread(file.path(dir, 'all_gene_s_p_with_FPKM.csv')) -> deg1

deg1[!is.na(deg1$padj),] -> deg1

# filter genes based on average FPKM
data.frame(deg1) -> deg1
deg1[,grepl('FPKM', colnames(deg1))] %>% rowMeans() -> deg1$mean_fpkm
deg1[deg1$mean_fpkm>1,] -> deg2


# run this if not the first time - to save some time
out7.unique_gene = fread('processed_results.csv')

sum(out7.unique_gene$padj<0.05)/nrow(out7.unique_gene)
sum(deg2$padj<0.05)/nrow(deg2)

write.csv(out7.unique_gene, 'processed_results.csv')



# Figure 4

# sig gene
deg2.sig = deg2[deg2$padj<0.05,]
boxplot(log2(out7.unique_gene$mean_fpkm), log2(deg2$mean_fpkm), outline=F)

# new 20231123 - boxplot-jitter
plot.dat1 = data.frame(key = 'se_target_gene', value = log2(out7.unique_gene$mean_fpkm))
plot.dat2 = data.frame(key = 'all_gene', value = log2(deg2$mean_fpkm))
plot.dat = rbind(plot.dat1, plot.dat2)

plot.dat$key = factor(plot.dat$key, levels = c('se_target_gene', 'all_gene'))

for(k in unique(plot.dat$key)){
  boxplot_filter(plot.dat[plot.dat$key==k,]) -> plot.dat[plot.dat$key==k,]
}

saveRDS(plot.dat, 'FigX4-IMR90.rds')

p = ggplot(plot.dat, aes(x = key, y = value, fill = key)) +
  geom_jitter(width = 0.2, aes(color = key), size = 1.5, alpha = 0.1) +  # Add jitter with matching color
  #geom_violin(trim=T, adjust=2, width=0.8, color='black', alpha=0.7) +  # Hide the outliers in the boxplot
  geom_boxplot(outlier.shape = NA, alpha=0.5) +
  scale_colour_brewer(palette = "Set2") +
  scale_fill_brewer(palette = "Set2") +
  theme_classic()
p

wilcox.test(log2(out7.unique_gene$mean_fpkm), log2(deg2$mean_fpkm))

getwd()
# thoughts
# 1. SE target genes have significant higher gene expression value
# 2. If there is a trend that high expressed genes tend to unchange their expression
# 3. Then, I can only use the genes that have equal expressed value, and check if DEG is enriched.

# Extra, use enhancer except super-enhancer to run the analysis.
# OR, filter the SE region, do not use the whole SE region, but use the enhancer region that located in SE. (dont use the linker region)


# Figure 1

# make plot

## 1. scater plot

# YY1 KD scatter plot
dat.scater = deg2 %>% dplyr::select(ensembl_gene_id, baseMean,log2FoldChange, padj)

try(dat.scater[dat.scater$padj==0,]$padj <- 1e-50, silent=T)

dat.scater %>% mutate(sig=cut(padj, breaks=c(0,0.05,1), labels=c('sig', 'notsig'))) -> dat.scater
dat.scater %>% mutate(direction=cut(log2FoldChange, breaks=c(-Inf, 0, Inf), labels=c('DN', 'UP'))) -> dat.scater
dat.scater %>% mutate(group=paste0(dat.scater$sig, dat.scater$direction))-> dat.scater


scatterplot.g = ggplot(dat.scater, aes(x=log10(baseMean), y=log2FoldChange, color=group)) + geom_point(size=0.1) +
  theme_classic() +
  scale_color_manual(values=c('gray', 'gray', 'blue','red'))

scatterplot.g

output_ppt(scatterplot.g, 'All gene RNA-seq scatter plot.pptx')


# Figure 2

# se target gene

out7.unique_gene$hgnc_symbol %>% unique() -> se.target.gene
deg2[deg2$hgnc_symbol %in% se.target.gene,] -> deg2.setarget

dat.scater = deg2.setarget %>% dplyr::select(ensembl_gene_id, baseMean,log2FoldChange, padj)

try(dat.scater[dat.scater$padj==0,]$padj <- 1e-50, silent=T)

dat.scater %>% mutate(sig=cut(padj, breaks=c(0,0.05,1), labels=c('sig', 'notsig'))) -> dat.scater
dat.scater %>% mutate(direction=cut(log2FoldChange, breaks=c(-Inf, 0, Inf), labels=c('DN', 'UP'))) -> dat.scater
dat.scater %>% mutate(group=paste0(dat.scater$sig, dat.scater$direction))-> dat.scater


scatterplot.g = ggplot(dat.scater, aes(x=log10(baseMean), y=log2FoldChange, color=group)) + geom_point(size=0.1) +
  theme_classic() +
  scale_color_manual(values=c('gray', 'gray', 'blue','red'))

scatterplot.g



