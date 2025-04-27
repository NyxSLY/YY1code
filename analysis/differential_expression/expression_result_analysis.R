library(data.table)
library(dplyr)
library(stringr)
#library(commonAPI)
source('C:\\Users\\Luyang_pc\\Documents\\BaiduNetdiskWorkspace\\OneDrive\\scripts/R/commonAPI.R')
library(ggplot2)


#setwd('D:\\BaiduNetdiskWorkspace\\OneDrive\\work\\3dGenome\\Find_enhancer_target\\total_normalized')
setwd('C:\\Users\\Luyang_pc\\Documents\\BaiduNetdiskWorkspace\\OneDrive\\work\\3dGenome\\Find_enhancer_target\\total_normalized')


dirs = list.dirs()
p5_dirs = dirs[grepl('P5', dirs)]
p15_dirs = dirs[grepl('P15', dirs)]
sig_se = 'D:\\BaiduNetdiskWorkspace\\OneDrive\\work\\3dGenome\\Find_enhancer_target\\super_enhancer.bed'
sig_se = 'C:\\Users\\Luyang_pc\\Documents\\BaiduNetdiskWorkspace\\OneDrive\\work\\3dGenome\\Find_enhancer_target\\super_enhancer.bed'

p5_dat = 'D:\\BaiduNetdiskWorkspace\\OneDrive\\work\\3dGenome\\Find_enhancer_target\\P5_SE_and_its_target.csv'
p15_dat = 'D:\\BaiduNetdiskWorkspace\\OneDrive\\work\\3dGenome\\Find_enhancer_target\\P15_SE_and_its_target.csv'

p5_dat='C:\\Users\\Luyang_pc\\Documents\\BaiduNetdiskWorkspace\\OneDrive\\work\\3dGenome\\Find_enhancer_target\\P5_SE_and_its_target.csv'
p15_dat = 'C:\\Users\\Luyang_pc\\Documents\\BaiduNetdiskWorkspace\\OneDrive\\work\\3dGenome\\Find_enhancer_target\\P15_SE_and_its_target.csv'



process = function(file, model, chr, sig_se, enh.psy, binsize=40000, sp){
  model = fread(model)
  model = as.matrix(model)
  file = fread(file)
  file = file[,-1]
  file = as.matrix(file)
  ratio = (file+1)/(model+1)
  ratio = as.matrix(ratio)
  #dat = fread(sig_se)
  dat = sig_se
  dat.chr = dplyr::filter(dat, seqnames==chr)
  dat.chr$anchor_bin = dat.chr$se_end/binsize
  temp =dat.chr$info %>% str_split(':',simplify = T)
  dat.chr$distance = gsub('Kb', '', temp[,2]) %>% as.integer() * 1000

  dat.chr.up = dat.chr

  out = list()
  flank_out = list()
  if(nrow(dat.chr.up)==0){
    return(list(out=out, flank_out = flank_out, pct_overlap=NA))
  }
  
  for(i in 1:nrow(dat.chr.up)){
    line = dat.chr.up[i,]
    a = ratio[line$anchor_bin, (line$anchor_bin+line$distance/40000)]
    if((line$anchor_bin-3)<1 || (line$anchor_bin+3)> nrow(ratio)){
      return(list(out=out, flank_out = flank_out, pct_overlap=NA))
    }
    flank = ratio[(line$anchor_bin-3):(line$anchor_bin+3), (line$anchor_bin+line$distance/40000)]
    out[[i]] = a
    flank_out[[i]] = flank
  }
  
  return(list(out=out, flank_out = flank_out, pct_overlap=NA))
}


main = function(dir, sp, sig_se){
    print(dir)
    chr = str_extract(dir, "chr[0-9]+|chr[X|Y]")
    files = list.files(dir)
    model = file.path(dir,files[grepl('model.estimated.matrix.txt', files)])
    file = file.path(dir, files[grepl(paste(sub('./','',dir),chr,'matrix.txt', sep='.'), files)])
    enh.psy = file.path(dir,files[grepl('5e-2.bed', files)])
    return(process(file, model, chr, sig_se, enh.psy, binsize=40000, sp=sp))
}


run = function(dirs, sp, sig_se){
    result = list()
    for(dir in dirs){
        result[[sub('./','',dir)]] = main(dir, sp=sp, sig_se)
    }
    x = lapply(result, '[[', 2)
    xx = as.data.frame(t(matrix(unlist(x), nrow=length(x[[1]][[1]]))))
    barplot(colMeans(xx, na.rm=T))
    return(result)
}


# Use SE as a starting point to plot the corresponding gene promoter and surrounding HiC map

# The main thing to look at is the comparison of HiC plots for common, P5SP, P15SP regions between P15HiC and P5HiC

p15_dat = fread(p15_dat)
p5_dat = fread(p5_dat)

p15_dat.sp = dplyr::filter(p15_dat, !grepl('common', p15_dat$region_id))
p5_dat.sp = dplyr::filter(p5_dat, !grepl('common', p5_dat$region_id))

p15_dat.com = dplyr::filter(p15_dat, grepl('common', p15_dat$region_id))
p5_dat.com = dplyr::filter(p5_dat, grepl('common', p5_dat$region_id))



result.p15.O = run(p15_dirs, 'O', p15_dat)
result.p5.O = run(p5_dirs, 'O', p15_dat)

result.p5.Y = run(p5_dirs, 'O', p5_dat)
result.p15.Y = run(p15_dirs, 'O', p5_dat)


# result.p15.O = run(p15_dirs, 'O', p15_dat.sp)
# result.p5.O = run(p5_dirs, 'O', p15_dat.sp)
# 
# result.p5.Y = run(p5_dirs, 'O', p5_dat.sp)
# result.p15.Y = run(p15_dirs, 'O', p5_dat.sp)
# 
# 
# result.p15.O = run(p15_dirs, 'O', p15_dat.com)
# result.p5.O = run(p5_dirs, 'O', p15_dat.com)
# 
# result.p5.Y = run(p5_dirs, 'O', p5_dat.com)
# result.p15.Y = run(p15_dirs, 'O', p5_dat.com)


result.p15.O.x = lapply(result.p15.O, '[[', 2)
result.p15.O.xx = as.data.frame(t(matrix(unlist(result.p15.O.x), nrow=length(result.p15.O.x[[1]][[1]]))))
result.p15.Y.x = lapply(result.p15.Y, '[[', 2)
result.p15.Y.xx = as.data.frame(t(matrix(unlist(result.p15.Y.x), nrow=length(result.p15.Y.x[[1]][[1]]))))

result.p5.Y.x = lapply(result.p5.Y, '[[', 2)
result.p5.Y.xx = as.data.frame(t(matrix(unlist(result.p5.Y.x), nrow=length(result.p5.Y.x[[1]][[1]]))))
result.p5.O.x = lapply(result.p5.O, '[[', 2)
result.p5.O.xx = as.data.frame(t(matrix(unlist(result.p5.O.x), nrow=length(result.p5.O.x[[1]][[1]]))))

# ratio = result.p15.O.xx / result.p5.O.xx
# barplot(colMeans(ratio, na.rm=T))
# 
# ratio = result.p5.Y.xx / result.p15.Y.xx
# barplot(colMeans(ratio, na.rm=T))

ggplot2_barplot <- function(x){
    x.colmean = data.frame(name=seq(1:length(x)),value=colMeans(x, na.rm=T))
    p = ggplot(data=data.frame(x.colmean), aes(x=name, y=value)) + 
        geom_bar(stat='identity', fill='blue', color='black') + theme_minimal() +
      coord_cartesian(ylim = c(1, 2.5))+
        theme(
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            panel.border = element_blank(),
            panel.grid=element_blank(),
            axis.ticks.x = element_blank(),
            axis.ticks.y = element_line(color="black"),
            axis.text.x = element_blank(),
            axis.text.y = element_text(size=20, color='black'),
            axis.ticks.length=unit(.25, "cm"),
            axis.line.x = element_line(color="black", size=0.5),
            axis.line.y = element_line(color="black",, size=0.5),
            #legend.title = element_blank(),
            legend.position = "none"
        )
    return(p)
}

boxplot_filter(test)


output_ppt_m(ggplot2_barplot(result.p15.O.xx),
             'barplot of fig2b.pptx')

output_ppt_m(ggplot2_barplot(result.p5.Y.xx),
             'barplot of fig2b.pptx')

boxplot(result.p15.O.xx$V1,
        result.p15.O.xx$V2,
        result.p15.O.xx$V3,
        result.p15.O.xx$V4,
        result.p15.O.xx$V5,
        result.p15.O.xx$V6,
        result.p15.O.xx$V7,
        outline=F)


#test = melt(setDT(result.p15.O.xx))
#ggplot(boxplot_filter(test), aes(x=variable, y=value)) + geom_boxplot() + geom_jitter(shape = 15,
#                                                                                      color = "steelblue",
#                                                                                      position = position_jitter(width = 0.21))



setwd('/Volumes/luyang/work/3dGenome/enhancer/PSYCHIC/result')
pdf('P5 SP SE @ P5 matrix.pdf', w=5, h=3.5)
ggplot2_barplot(result.p5.Y.xx)
dev.off()

pdf('P15 SP SE @ P15 matrix.pdf', w=5, h=3.5)
ggplot2_barplot(result.p15.O.xx)
dev.off()



barplot(colMeans(result.p15.O.xx, na.rm=T))
barplot(colMeans(result.p5.O.xx, na.rm=T))
barplot(colMeans(result.p15.Y.xx, na.rm=T))
barplot(colMeans(result.p5.Y.xx, na.rm=T))

dim(result.p15.O.xx)
dim(result.p15.Y.xx)

# overlap with psy.enhancer

result.p15.O.overlap = lapply(result.p15.O, '[[', 3)
result.p15.O.overlap = as.data.frame(unlist(result.p15.O.overlap))

result.p15.Y.overlap = lapply(result.p15.Y, '[[', 3)
result.p15.Y.overlap = as.data.frame(unlist(result.p15.Y.overlap))

result.p5.O.overlap = lapply(result.p5.O, '[[', 3)
result.p5.O.overlap = as.data.frame(unlist(result.p5.O.overlap))

result.p5.Y.overlap = lapply(result.p5.Y, '[[', 3)
result.p5.Y.overlap = as.data.frame(unlist(result.p5.Y.overlap))

overlap.df = data.frame(p5.y=result.p5.Y.overlap[,1],
                        p5.o=result.p5.O.overlap[,1],
                        p15.y=result.p15.Y.overlap[,1],
                        p15.o=result.p15.O.overlap[,1])

barplot(colMeans(overlap.df, na.rm=T))
#-----------------------------



library(data.table)

#setwd('/Volumes/luyang-1/work/3dGenome/enhancer/PSYCHIC/examples/output')
model = '/Volumes/luyang-1/work/3dGenome/enhancer/PSYCHIC/result/evaluation/IC_test/output/P15_IC.chr21.model.estimated.matrix.txt'
model = fread(model)
model = as.matrix(model)

file = '/Volumes/luyang-1/work/3dGenome/enhancer/PSYCHIC/result/evaluation/IC_test/P5.chr21.matrix.txt'
file = fread(file)
file = file[,-1]
file = as.matrix(file)
dim(model)
dim(file)

ratio = (file+1)/(model+1)
ratio = as.matrix(ratio)

dat = '/Volumes/luyang-1/work/3dGenome/enhancer/4c/changed_SE/all_SE/40kb_total_norm/sig_regions.txt'
dat = fread(dat)

dat.chr = filter(dat, seqnames=='chr21')
dat.chr$anchor_bin = dat.chr$anchor_end/40000
dat.chr.up = filter(dat.chr, P15_readcount<P5_readcount)

out = list()
flank_out = list()
for(i in 1:nrow(dat.chr.up)){
    line = dat.chr.up[i,]
    a = ratio[line$anchor_bin, (line$anchor_bin+line$distance/40000)]
    flank = ratio[(line$anchor_bin-3):(line$anchor_bin+3), (line$anchor_bin+line$distance/40000)]
    out[[i]] = a
    flank_out[[i]] = flank
}

for(i in 1:16){
    barplot(flank_out[[i]], main = i)
}

# overlap between SE by H3K27ac and enhancer by HIC data PSYCHIC
enh.psy = '/Volumes/luyang-1/work/3dGenome/enhancer/PSYCHIC/result/evaluation/IC_test/output/P15_IC.chr21.enh_5e-2.bed'
enh.psy = readbed(enh.psy)

# enh.psy overlaps with H3K27ac peak
k27ac.o = readbed('/Volumes/Data1/ChIP_seq_11_2018/macs2/peak_bed/O3_H3K27ac_peaks.bed')
k27ac.y = readbed('/Volumes/Data1/ChIP_seq_11_2018/macs2/peak_bed/Y3_H3K27ac_peaks.bed')

enh.psy.k27ac = enh.psy[overlapsAny(enh.psy, k27ac.y)]
length(enh.psy.k27ac)/length(enh.psy)

sig.se.up = makeGRangesFromDataFrame(dat.chr.up[, c('seqnames', "anchor_start", 'anchor_end')], keep.extra.columns = T,
                                     start.field = 'anchor_start',
                                     end.field = 'anchor_end')
o = overlapsAny(sig.se.up, enh.psy)
sum(o)/length(o)  # no need to use enh.psy.k27ac, as all sig se overlap with k27ac peak for sure.


out.model = ''
for(i in 1:nrow(dat.chr.up)){
    line = dat.chr.up[i,]
    a = model[line$anchor_bin, (line$anchor_bin+line$distance/40000)]
    out.model[[i]] = as.numeric(a)
}


dat.chr.up$model_count = unlist(out.model)
dat.chr.up$if_psy = o

# so, half of the super-enhancer in P15 chr1 do not have sig more read counts than
# expected. But, has more reads than P5. So what is looks like in P5


# filter super-enhancer list before find sig changed SE?

# how did the paper define P-E pair? what if a enhancer regulate multiple promoter?





