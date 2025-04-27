# pathway and go-term
library(data.table)
library(rtracklayer)
library(GenomicRanges)
library(ChIPseeker)
library(clusterProfiler)
library(ReactomePA)
library(ggplot2)
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
load('BASE.RData')

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