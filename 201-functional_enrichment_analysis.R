rm(list = ls())
library(clusterProfiler)
library(annotables)
library(enrichplot)
load("result/intersect_Gene.RData")
data("grch37")

downgenlist <- dplyr::filter(grch37, ensgene %in% intersectdown)

upgenlist <- dplyr::filter(grch37, ensgene %in% intersetup)

upeg = bitr(unique(upgenlist$symbol), 
          fromType="SYMBOL", toType=c("ENTREZID", "SYMBOL"), OrgDb="org.Hs.eg.db")

downeg <- bitr(unique(downgenlist$symbol), 
               fromType="SYMBOL", toType=c("ENTREZID", "SYMBOL"), OrgDb="org.Hs.eg.db")


gmtfile <- system.file("extdata", "c2.cp.kegg.v6.1.entrez.gmt", package="clusterProfiler")

kegg <- read.gmt(gmtfile)

################### UP Gene

egmtkegg <- enricher(upeg$ENTREZID, 
                     TERM2GENE=kegg, 
                     pvalueCutoff = 0.9,
                     qvalueCutoff = 0.9,
                     minGSSize = 3)

dot5 <- dotplot(egmtkegg, showCategory=15, 
                color = "pvalue",
                x="GeneRatio",
                font.size = 12)

ggsave(plot = dot5, 
       filename = "result/up_kegg_enrich.tiff", 
       device = "tiff", 
       dpi = 300, 
       width = 24, 
       height = 15,
       units = "cm")

go <- enrichGO(upeg$ENTREZID, 
               OrgDb = "org.Hs.eg.db", 
               ont="bp",
               pvalueCutoff = 0.9,
               qvalueCutoff = 0.9, 
               minGSSize = 3)


dott <- barplot(go, showCategory=15, color = "pvalue") 

ggsave(plot = dott, 
       filename = "result/up_GO_enrich.tiff", 
       device = "tiff", 
       dpi = 300, 
       width = 24, 
       height = 15,
       units = "cm")

############## down Gene
egmtkegg <- enricher(downeg$ENTREZID, 
                     TERM2GENE=kegg, 
                     pvalueCutoff = 0.9,
                     qvalueCutoff = 0.9,
                     minGSSize = 3) 



dot5 <- dotplot(egmtkegg, showCategory=15, 
                color = "pvalue",
                x="GeneRatio",
                font.size = 10)

ggsave(plot = dot5, 
       filename = "result/down_kegg_enrich.tiff", 
       device = "tiff", 
       dpi = 300, 
       width = 24, 
       height = 15,
       units = "cm")

go <- enrichGO(downeg$ENTREZID, 
               OrgDb = "org.Hs.eg.db", 
               ont="bp",
               pvalueCutoff = 0.9,
               qvalueCutoff = 0.9, 
               minGSSize = 3)


dott <- barplot(go, showCategory=15, color = "pvalue") 

ggsave(plot = dott, 
       filename = "result/down_GO_enrich.tiff", 
       device = "tiff", 
       dpi = 300, 
       width = 24, 
       height = 15,
       units = "cm")




goread <- as.data.frame(go)
keggread <- egmtkegg@result

write.csv(goread, "result/enrich_tcga_DEG_go.csv", row.names = F, quote = F)
write.csv(keggread, "result/enrich_tcga_DEG_kegg.csv", row.names = F, quote = F)