rm(list = ls())
library(UpSetR)
load("result/GEO-GSE104131_mRNA_DEG.RData")
mrna_DEG$names <- rownames(mrna_DEG)

GSE104131up <- dplyr::filter(mrna_DEG, pvalue <= 0.01, log2FoldChange > 0.263) %>% 
              dplyr::select(names) %>% 
              unlist() %>% 
              unname()
GSE104131down <- dplyr::filter(mrna_DEG, pvalue <= 0.01, log2FoldChange < -0.263)%>% 
               dplyr::select(names) %>% 
               unlist() %>% 
                unname()
#############
load("result/GEO-GSE89223_mRNA_DEG.RData")
mrna_DEG$names <- rownames(mrna_DEG)

GSE89223up <- dplyr::filter(mrna_DEG, pvalue <= 0.01, log2FoldChange > 0.263) %>% 
  dplyr::select(names) %>% 
  unlist() %>% 
  unname()
GSE89223down <- dplyr::filter(mrna_DEG, pvalue <= 0.01, log2FoldChange < -0.263)%>% 
  dplyr::select(names) %>% 
  unlist() %>% 
  unname()

#############
load("result/TCGA-PRAD_mRNA_DEG.RData")
mrna_DEG$names <- rownames(mrna_DEG)

TCGA_PRADup <- dplyr::filter(mrna_DEG, pvalue <= 0.01, log2FoldChange > 0.263) %>% 
  dplyr::select(names) %>% 
  unlist() %>% 
  unname()
TCGA_PRADdown <- dplyr::filter(mrna_DEG, pvalue <= 0.01, log2FoldChange < -0.263)%>% 
  dplyr::select(names) %>% 
  unlist() %>% 
  unname()
########### upsetR

uiongeneup <- union(union(GSE104131up, GSE89223up), TCGA_PRADup)
uiongenedown <- union(union(GSE104131down, GSE89223down), TCGA_PRADdown)

indexup <- lapply(list(GSE104131up, GSE89223up, TCGA_PRADup), function(x){as.numeric(uiongeneup %in% x)}) %>% 
          as.data.frame()

names(indexup) <- c("GSE104131up", "GSE89223up", "TCGA_PRADup")
rownames(indexup) <- uiongeneup

indexdown <- lapply(list(GSE104131down, GSE89223down, TCGA_PRADdown), function(x){as.numeric(uiongenedown %in% x)}) %>% 
  as.data.frame()

names(indexdown) <- c("GSE104131down", "GSE89223down", "TCGA_PRADdown")


####
tiff("Plot00.tiff", width = 2600, height = 1500, units = 'px', res = 300)
 upset(indexup,
      main.bar.color = "blue", 
      matrix.color = "red",
      sets.bar.color = "blue",
      point.size = 3,
      line.size = 1,
      mb.ratio = c(0.7, 0.3),
      text.scale = 2) 
 dev.off()

 tiff("Plot1.tiff", width = 2600, height = 1500, units = 'px', res = 300)
 upset(indexdown,
       main.bar.color = "blue", 
       matrix.color = "red",
       sets.bar.color = "blue",
       point.size = 3,
       line.size = 1,
       mb.ratio = c(0.7, 0.3),
       text.scale = 2) 
 dev.off()

 intersetup <- intersect(intersect(GSE104131up, GSE89223up), TCGA_PRADup)
 intersectdown <- intersect(intersect(GSE104131down, GSE89223down), TCGA_PRADdown)
 

save(intersectdown, intersetup, file = "result/intersect_Gene.RData")


###########################3
## miRNA 画图
rm(list = ls())

mirbase <- read_csv("rawdata/mirbase_aliases.csv", col_names = F)

maturename <- str_split(mirbase$X2, ";")

mibaseindex <- data.frame(ids = rep(mirbase$X1, sapply(maturename, length)),
                          names = unlist(maturename)) %>% 
               dplyr::filter(names != "")

save(mibaseindex, file = "rawdata/mibaseIndex.RData")
###

load("result/GEO-GSE75415_miRNA_DEG.RData")
mrna_DEG$names <- rownames(mrna_DEG)

mrna_DEG <- inner_join(mrna_DEG, mibaseindex, by = c("names" = "names"))

GSE75415up <- dplyr::filter(mrna_DEG, pvalue <= 0.05, log2FoldChange > 0.263) %>% 
  dplyr::select(ids) %>% 
  unlist() %>% 
  as.character() %>% 
  unname()
GSE75415down <- dplyr::filter(mrna_DEG, pvalue <= 0.05, log2FoldChange < -0.263)%>% 
  dplyr::select(ids) %>% 
  unlist() %>% 
  as.character() %>% 
  unname()


#######

load("result/GEO-GSE76260_miRNA_DEG.RData")
mrna_DEG$names <- rownames(mrna_DEG)

paltform_index <- read_tsv("rawdata/miRNA/GSE76260/GPL8179_humanMI_V2_R0_XS0000124-MAP.txt/GPL8179_humanMI_V2_R0_XS0000124-MAP.txt")

targetname <- str_split(paltform_index$TargetMatureName, ",")

GPL8179index <- data.frame(ids = rep(paltform_index$SYMBOL, sapply(targetname, length)),
                          names = unlist(targetname)) %>% 
  dplyr::filter(names != "") %>% 
  distinct() %>% 
  inner_join(mibaseindex, by = c("names" = "names"))



mrna_DEG <- inner_join(mrna_DEG, GPL8179index, by = c("names" = "ids.x"))

GSE76260up <- dplyr::filter(mrna_DEG, pvalue <= 0.05, log2FoldChange > 0.263) %>% 
  dplyr::select(ids.y) %>% 
  unlist() %>% 
  as.character() %>% 
  unname()
GSE76260down <- dplyr::filter(mrna_DEG, pvalue <= 0.05, log2FoldChange < -0.263)%>% 
  dplyr::select(ids.y) %>% 
  unlist() %>% 
  as.character() %>% 
  unname()

###############

load("result/TCGA-PRAD_miRNA_DEG.RData")

mrna_DEG$names <- rownames(mrna_DEG)

mrna_DEG <- inner_join(mrna_DEG, mibaseindex, by = c("names" = "names"))

TCGA_PRADup <- dplyr::filter(mrna_DEG, pvalue <= 0.05, log2FoldChange > 0.263, baseMean > 0) %>% 
  dplyr::select(ids) %>% 
  unlist() %>% 
  as.character() %>% 
  unname()
TCGA_PRADdown <- dplyr::filter(mrna_DEG, pvalue <= 0.05, log2FoldChange < -0.263, baseMean > 0)%>% 
  dplyr::select(ids) %>% 
  unlist() %>% 
  as.character() %>% 
  unname()

####################

uiongeneup <- union(union(GSE75415up, GSE76260up), TCGA_PRADup)
uiongenedown <- union(union(GSE75415down, GSE76260down), TCGA_PRADdown)

indexup <- lapply(list(GSE75415up, GSE76260up, TCGA_PRADup), function(x){as.numeric(uiongeneup %in% x)}) %>% 
  as.data.frame()

names(indexup) <- c("GSE75415up", "GSE76260up", "TCGA_PRADup")
rownames(indexup) <- uiongeneup

indexdown <- lapply(list(GSE75415down, GSE76260down, TCGA_PRADdown), function(x){as.numeric(uiongenedown %in% x)}) %>% 
  as.data.frame()

names(indexdown) <- c("GSE75415down", "GSE76260down", "TCGA_PRADdown")


####
tiff("Plot00_mirna.tiff", width = 2600, height = 1500, units = 'px', res = 300)
upset(indexup,
      main.bar.color = "blue", 
      matrix.color = "red",
      sets.bar.color = "blue",
      point.size = 3,
      line.size = 1,
      mb.ratio = c(0.7, 0.3),
      text.scale = 2) 
dev.off()

tiff("Plot1_miRNA.tiff", width = 2600, height = 1500, units = 'px', res = 300)
upset(indexdown,
      main.bar.color = "blue", 
      matrix.color = "red",
      sets.bar.color = "blue",
      point.size = 3,
      line.size = 1,
      mb.ratio = c(0.7, 0.3),
      text.scale = 2) 
dev.off()

intersetup <- intersect(intersect(GSE75415up, GSE76260up), TCGA_PRADup)
intersectdown <- intersect(intersect(GSE75415down, GSE76260down), TCGA_PRADdown)


save(intersectdown, intersetup, file = "result/intersect_miRNA.RData")


