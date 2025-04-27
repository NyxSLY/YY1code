library(data.table)
library(dplyr)
library(stringr)
source('/Users/Luyang/scripts/R/commonAPI.R')

script_path=dirname(rstudioapi::getActiveDocumentContext()$path)
#setwd('/Volumes/Luyang1/onedriveMe/OneDrive/work/3dGenome/Find_enhancer_target/total_normalized')
wd='/Volumes/Luyang1/onedriveMe/OneDrive/work/3dGenome/Find_enhancer_target'
#wd='C://Users//Luyang.000//OneDrive//work//3dGenome//Find_enhancer_target'
setwd(wd)

dirs = list.dirs(file.path(wd, 'total_normalized'))
p5_dirs = dirs[grepl('P5', dirs)]
p15_dirs = dirs[grepl('P15', dirs)]
pair_region_file1 = file.path(script_path,'P5_SE_and_its_target.csv')
pair_region_file2 = file.path(script_path,'P15_SE_and_its_target.csv')

pair_region1 = fread(pair_region_file1)
pair_region2 = fread(pair_region_file2)


rbind(pair_region1, pair_region2[grepl('p15',pair_region2$region_id),]) -> pair_region



process = function(file, model, chr, pair_region_file, binsize=40000, control=T){
    model = fread(model)
    model = as.matrix(model)
    file = fread(file)
    file = file[,-1]
    file = as.matrix(file)
    ratio = (file+1)/(model+1)
    ratio = as.matrix(ratio)
    dat = fread(pair_region_file)
    dat.chr = dplyr::filter(dat, seqnames==chr)
    out = list()
    flank_out = list()
    if(nrow(dat.chr)==0){
        return(NA)
    }
    
    dat.chr$target_bin = ceiling(dat.chr$target_region_end/binsize)
    distance = as.numeric(gsub('Kb', '', sapply(str_split(dat.chr$info, ':'), '[', 2)))/binsize*1000*-1
    dat.chr$anchor_bin = dat.chr$target_bin + distance
    dat.chr$control_bin = dat.chr$target_bin+3
    
    for(i in 1:nrow(dat.chr)){
        line = dat.chr[i,]
        readcount = file[line$anchor_bin, line$target_bin]
        ratio_count = ratio[line$anchor_bin, line$target_bin]
        #flank = ratio[(line$anchor_bin-3):(line$anchor_bin+3), line$target_bin]
        out[[line$pair]] = ratio_count
        flank_out[[line$pair]] = flank
    }
    
    # overlap between SE by H3K27ac and enhancer by HIC data PSYCHIC
    # enh.psy = readbed(enh.psy)
    # sig.se.up = makeGRangesFromDataFrame(dat.chr[, c('seqnames', "anchor_start", 'anchor_end')], keep.extra.columns = T,
    #                                      start.field = 'anchor_start',
    #                                      end.field = 'anchor_end')
    # sig.se.up = reduce(sig.se.up)
    # o = overlapsAny(sig.se.up, enh.psy)
    # pct_overlap = sum(o)/length(o)
    
    return(out)
}

process1 = function(file, model, chr, pair_region, binsize=40000, control=T){
    model = fread(model)
    model = as.matrix(model)
    file = fread(file)
    file = file[,-1]
    file = as.matrix(file)
    ratio = (file+1)/(model+1)
    ratio = as.matrix(ratio)
    dat = pair_region
    dat.chr = dplyr::filter(dat, seqnames==chr)
    out = list()
    flank_out = list()
    if(nrow(dat.chr)==0){
        return(NA)
    }
    
    dat.chr$target_bin = ceiling(dat.chr$target_region_end/binsize)
    distance = as.numeric(gsub('Kb', '', sapply(str_split(dat.chr$info, ':'), '[', 2)))/binsize*1000*-1
    dat.chr$anchor_bin = dat.chr$target_bin + distance
    dat.chr$control_bin = dat.chr$target_bin+3
    
    for(i in 1:nrow(dat.chr)){
        line = dat.chr[i,]
        readcount = file[line$anchor_bin, line$target_bin]
        ratio_count = ratio[line$anchor_bin, line$target_bin]
        bg_readcount = model[line$anchor_bin, line$target_bin]
        #flank = ratio[(line$anchor_bin-3):(line$anchor_bin+3), line$target_bin]
        out[[line$pair]] = list(ratio=ratio_count, rc=readcount, bg_rc=bg_readcount)
        flank_out[[line$pair]] = flank
    }
    
    # overlap between SE by H3K27ac and enhancer by HIC data PSYCHIC
    # enh.psy = readbed(enh.psy)
    # sig.se.up = makeGRangesFromDataFrame(dat.chr[, c('seqnames', "anchor_start", 'anchor_end')], keep.extra.columns = T,
    #                                      start.field = 'anchor_start',
    #                                      end.field = 'anchor_end')
    # sig.se.up = reduce(sig.se.up)
    # o = overlapsAny(sig.se.up, enh.psy)
    # pct_overlap = sum(o)/length(o)
    
    return(out)
}

main1 = function(dir, pair_region_file){
    print(dir)
    chr = str_extract(dir, "chr[0-9]+|chr[X|Y]")
    files = list.files(dir)
    model = file.path(dir,files[grepl('model.estimated.matrix.txt', files)])
    file = file.path(dir, files[grepl(paste(basename(dir),chr,'matrix.txt', sep='.'), files)])
    return(process1(file, model, chr, pair_region_file, binsize=40000))
}

run1 = function(dirs, pair_region_file){
    result = list()
    for(dir in dirs){
        result[[basename(dir)]] = main1(dir, pair_region_file)
    }
    #x = lapply(result, '[[', 2)
    #xx = as.data.frame(t(matrix(unlist(x), nrow=length(x[[1]][[1]]))))
    #barplot(colMeans(xx, na.rm=T))
    return(result)
}


main = function(dir, pair_region_file){
    print(dir)
    chr = str_extract(dir, "chr[0-9]+|chr[X|Y]")
    files = list.files(dir)
    model = file.path(dir,files[grepl('model.estimated.matrix.txt', files)])
    file = file.path(dir, files[grepl(paste(basename(dir),chr,'matrix.txt', sep='.'), files)])
    return(process(file, model, chr, pair_region_file, binsize=40000))
}

run = function(dirs, pair_region_file){
    result = list()
    for(dir in dirs){
        result[[basename(dir)]] = main(dir, pair_region_file)
    }
    #x = lapply(result, '[[', 2)
    #xx = as.data.frame(t(matrix(unlist(x), nrow=length(x[[1]][[1]]))))
    #barplot(colMeans(xx, na.rm=T))
    return(result)
}

result.p5.Y = run1(p5_dirs, pair_region)
result.p5.O = run1(p15_dirs, pair_region)

lapply(result.p5.Y, function(x){lapply(x, unlist)}) %>% unlist(recursive = F) %>% do.call('rbind', .) %>% 
    data.frame -> result.p5.Y.n
lapply(result.p5.O, function(x){lapply(x, unlist)}) %>% unlist(recursive = F) %>% do.call('rbind', .) %>% 
    data.frame -> result.p5.O.n

names(result.p5.Y.n) = c('p5_ratio', 'p5_readcount', 'p5_bg_readcount')
names(result.p5.O.n) = c('p15_ratio', 'p15_readcount', 'p15_bg_readcount')

# cbind(result.p5.Y.n, result.p5.O.n) -> result
# names(result) = c('P5', 'P15')
# 
# boxplot((result$P5),(result$P15), outline=F)
# grepl('common',rownames(result)) -> common_keep
# 
# result[common_keep,] -> result.com
# grepl('p5',rownames(result)) -> p5sp_keep
# grepl('p15',rownames(result)) -> p15sp_keep
# 
# result[p5sp_keep,] -> result.p5sp
# result[p15sp_keep,] -> result.p15sp
# 
# boxplot((result$P5),(result$P15), outline=F, main='YY1 DN')
# boxplot(result.com$P5,result.com$P15, outline=F, main='YY1 DN common')
# boxplot(result.p5sp$P5,result.p5sp$P15, outline=F, main='YY1 DN sp')
# boxplot(result.p15sp$P5,result.p15sp$P15, outline=F, main='YY1 DN sp')
# 
# wilcox.test(result$P5, result$P15, paired=T)
# wilcox.test(result.com$P5, result.com$P15, paired=T)
# 
# 
# result[!is.na(result$P5),] -> result
# rownames(result) <- sapply(strsplit(rownames(result), "\\."),'[[', 2)

getwd()
write.csv(result, 'all_se_pair_withHiCinfo.csv')


#----------------------
pdf(NULL)
dev.control(displaylist="enable")

boxplot((result$P5),(result$P15), outline=F, main='YY1 DN')

p.expression <- recordPlot()
invisible(dev.off())
output_ppt(p.expression, 'hic contact in all YY1 DN S-P pare.pptx')

#----------------------

pdf(NULL)
dev.control(displaylist="enable")
boxplot(result.com$P5,result.com$P15, outline=F, main='YY1 DN common')
p.expression <- recordPlot()
invisible(dev.off())
output_ppt(p.expression, 'hic contact in common YY1 DN S-P pare.pptx')

