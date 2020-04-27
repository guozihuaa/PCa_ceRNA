rm(list = ls())
library(limma)
library(annotables)
library(EnhancedVolcano)
library(tibble)
load("PRAD_RNA_and_clinical.RData")

## 这里要注意表达矩阵的列名的顺序是否与临床样本分类，样本ID顺序是否一致

raw_clinic_data <- filter(clinical2, barcode %in% colnames(lihc_ranexp2))

DEGAll <- gdcDEAnalysis(counts     = lihc_ranexp2, 
                        group      = raw_clinic_data$shortLetterCode, 
                        comparison = 'TP-NT', 
                        method     = 'limma')

mrna_DEG <- DEGAll %>% 
  dplyr::select(baseMean = AveExpr,
                log2FoldChange = logFC,
                lfcSE = B,
                stat = t,
                pvalue = PValue,
                padj = FDR)


save(mrna_DEG, file = "result/TCGA-PRAD_mRNA_DEG.RData")

###

pfit1 <- EnhancedVolcano(mrna_DEG,
                
               lab = rownames(mrna_DEG),
               
               selectLab = c("a"),
                x = "log2FoldChange",
                
                y = "pvalue",
               subtitle = "Differential expression",
               caption = "FC cutoff, 1.333; p-value cutoff, 10e-4",
                xlab = bquote(~Log[2]~ "fold change"),
                
                ylab = bquote(~-Log[10]~adjusted~italic(P)),
                
                pCutoff = 0.0001,
                
                FCcutoff = 1.0,
                
                xlim=c(-4,4),
                
                transcriptLabSize = 0,
                
                colAlpha = 1,
                
                legend=c("NS","Log2 FC","P-value",
                         "P-value & Log2 FC"),
                
                legendPosition = "top",
                
                legendLabSize = 12,
                
                legendIconSize = 3.0,
               title = "TCGA-PRAD Tumor versus Normal")

ggpubr::ggexport(plot = pfit1, 
                 filename = paste("TCGA-PRAD", ".tiff", sep = ""), 
                 device = "tiff", 
                 res = 300,
                 width = 1800, 
                 height = 2200)


########## GEO 数据火山图

rm(list = ls())

load("rawdata/mRNA_GSE89223_clinin.RData")

raw_clinic_data <- filter(sampleldescript, title %in% colnames(exprdata))

DEGAll <- gdcDEAnalysis(counts     = exprdata, 
                        group      = raw_clinic_data$type, 
                        comparison = 'Tumor-Normal', 
                        method     = 'limma')

mrna_DEG <- DEGAll %>% 
  dplyr::select(baseMean = AveExpr,
                log2FoldChange = logFC,
                lfcSE = B,
                stat = t,
                pvalue = PValue,
                padj = FDR)

save(mrna_DEG, file = "result/GEO-GSE89223_mRNA_DEG.RData")
##

pfit1 <- EnhancedVolcano(mrna_DEG,
                         
                         lab = rownames(mrna_DEG),
                         
                         selectLab = c("a"),
                         x = "log2FoldChange",
                         
                         y = "pvalue",
                         
                         xlab = bquote(~Log[2]~ "fold change"),
                         
                         ylab = bquote(~-Log[10]~"P-value"),
                         
                         pCutoff = 0.05,
                         
                         FCcutoff = 1,
                         
                         xlim=c(-4,4),
                         ylim = c(0,10),
                         
                         transcriptLabSize = 0,
                         
                         colAlpha = 1,
                         
                         legend=c("NS","Log2 FC","P-value",
                                  "P-value & Log2 FC"),
                         
                         legendPosition = "top",
                         
                         legendLabSize = 10,
                         
                         legendIconSize = 3.0,
                         title = "GSE89223 Tumor versus Normal")

ggpubr::ggexport(plot = pfit1, 
                 filename = paste("GSE89223", ".tiff", sep = ""), 
                 device = "tiff", 
                 res = 300,
                 width = 1800, 
                 height = 2200)

###### GSE104131

rm(list = ls())
load("rawdata/mRNA_GSE104131_clinin.RData")
raw_clinic_data <- filter(sampleldescript, title %in% colnames(fpkm1))

fpkm1 <- as.matrix(fpkm1)[, raw_clinic_data$title]

fpkm1 <- log2(fpkm1 + 1)
raw_clinic_data <- sampleldescript
group_list=as.character(raw_clinic_data$type)
design <- model.matrix(~0+factor(group_list))
colnames(design)=levels(factor(group_list))
rownames(design)=colnames(fpkm1)


contrast.matrix<-makeContrasts(paste0(unique(group_list),
                                      collapse = "-"),levels = design)
fit <- lmFit(fpkm1,design)

fit2 <- contrasts.fit(fit, contrast.matrix) 
fit2 <- eBayes(fit2, trend = T)

mrna_DEG = topTable(fit2, coef=1, n=Inf) %>%
  rownames_to_column('gene')  %>%
  dplyr::filter(AveExpr>0.7315545) %>% 
  dplyr::select(baseMean = AveExpr,
                log2FoldChange = logFC,
                lfcSE = B,
                stat = t,
                pvalue = P.Value,
                padj = adj.P.Val,
                gene)%>%
  column_to_rownames('gene')

 

save(mrna_DEG, file = "result/GEO-GSE104131_mRNA_DEG.RData")

####

pfit1 <- EnhancedVolcano(mrna_DEG,
                         
                         lab = rownames(mrna_DEG),
                         
                         selectLab = c("a"),
                         x = "log2FoldChange",
                         
                         y = "pvalue",
                         
                         xlab = bquote(~Log[2]~ "fold change"),
                         
                         ylab = bquote(~-Log[10]~"P-value"),
                         
                         pCutoff = 0.05,
                         
                         FCcutoff = 0.5849625,
                         
                         xlim=c(-3,3),
                         ylim = c(0,8),
                         
                         transcriptLabSize = 0,
                         
                         colAlpha = 1,
                         
                         legend=c("NS","Log2 FC","P-value",
                                  "P-value & Log2 FC"),
                         
                         legendPosition = "top",
                         
                         legendLabSize = 10,
                         
                         legendIconSize = 3.0,
                         title = "GSE104131 Tumor versus Normal")

ggpubr::ggexport(plot = pfit1, 
                 filename = paste("GSE104131", ".tiff", sep = ""), 
                 device = "tiff", 
                 res = 300,
                 width = 1800, 
                 height = 2200)

#######################################################################################################################################################################################    
#############################################################
##miRNA 数据
##############################

### TCGA
rm(list = ls())
load("rawdata/mirexpr.RData")
load("PRAD_RNA_and_clinical.RData")

raw_clinic_data <- clinical2 %>% 
               mutate(title = str_sub(sample, 1, -2)) %>% 
           filter(title %in% colnames(mirExpr))

mirMatrix <- as.matrix(mirMatrix)[, raw_clinic_data$title]

DEGAll <- gdcDEAnalysis(counts     = mirMatrix, 
                        group      = raw_clinic_data$shortLetterCode, 
                        comparison = 'TP-NT', 
                        method     = 'limma',
                        filter = F)

mrna_DEG <- DEGAll %>% 
  dplyr::select(baseMean = AveExpr,
                log2FoldChange = logFC,
                lfcSE = B,
                stat = t,
                pvalue = PValue,
                padj = FDR)


save(mrna_DEG, file = "result/TCGA-PRAD_miRNA_DEG.RData")

pfit1 <- EnhancedVolcano(mrna_DEG,
                         
                         lab = rownames(mrna_DEG),
                         
                         selectLab = c("a"),
                         x = "log2FoldChange",
                         
                         y = "pvalue",
                         
                         xlab = bquote(~Log[2]~ "fold change"),
                         
                         ylab = bquote(~-Log[10]~"P-value"),
                         
                         pCutoff = 0.05,
                         
                         FCcutoff = 0.26,
                         
                         xlim=c(-3,3),
                         #ylim = c(0,8),
                         
                         transcriptLabSize = 0,
                         
                         colAlpha = 1,
                         
                         legend=c("NS","Log2 FC","P-value",
                                  "P-value & Log2 FC"),
                         
                         legendPosition = "top",
                         transcriptPointSize = 1.1,
                         
                         legendLabSize = 10,
                         
                         legendIconSize = 3.0,
                         title = "TCGA-PRAD Tumor versus Normal")

ggpubr::ggexport(plot = pfit1, 
                 filename = paste("TCGA-miRNA", ".tiff", sep = ""), 
                 device = "tiff", 
                 res = 300,
                 width = 1800, 
                 height = 2200)


########### mirna GSE21036
rm(list = ls())
load("rawdata/miRNA_GSE21036_cline.RData")

exprdata <- as.matrix(exprdata)[, sampleldescript$title]
raw_clinic_data <- sampleldescript

group_list=as.character(raw_clinic_data$type)
design <- model.matrix(~0+factor(group_list))
colnames(design)=levels(factor(group_list))
rownames(design)=colnames(exprdata)


contrast.matrix<-makeContrasts(paste0(unique(group_list),
                                      collapse = "-"),levels = design)
fit <- lmFit(exprdata,design)

fit2 <- contrasts.fit(fit, contrast.matrix) 
fit2 <- eBayes(fit2, trend = T)

mrna_DEG = topTable(fit2, coef=1, n=Inf) %>% 
  dplyr::select(baseMean = AveExpr,
                log2FoldChange = logFC,
                lfcSE = B,
                stat = t,
                pvalue = P.Value,
                padj = adj.P.Val)

save(mrna_DEG, file = "result/GEO-GSE21036_miRNA_DEG.RData")

####

pfit1 <- EnhancedVolcano(mrna_DEG,
                         
                         lab = rownames(mrna_DEG),
                         
                         selectLab = c("a"),
                         x = "log2FoldChange",
                         
                         y = "pvalue",
                         
                         xlab = bquote(~Log[2]~ "fold change"),
                         
                         ylab = bquote(~-Log[10]~"P-value"),
                         
                         pCutoff = 0.05,
                         
                         FCcutoff = 0.26,
                         
                         xlim=c(-3,3),
                         ylim = c(0,15),
                         
                         transcriptLabSize = 0,
                         
                         colAlpha = 1,
                         
                         legend=c("NS","Log2 FC","P-value",
                                  "P-value & Log2 FC"),
                         
                         legendPosition = "top",
                         transcriptPointSize = 1.2,
                         
                         legendLabSize = 10,
                         
                         legendIconSize = 3.0,
                         title = "GSE21036 Tumor versus Normal")

ggpubr::ggexport(plot = pfit1, 
                 filename = paste("GSE21036", ".tiff", sep = ""), 
                 device = "tiff", 
                 res = 300,
                 width = 1800, 
                 height = 2200)


########### mirna GSE76260
rm(list = ls())
load("rawdata/miRNA_GSE76260_clinic.RData")

exprdata <- exprdata[, sampleldescript$title]
raw_clinic_data <- sampleldescript

group_list=as.character(raw_clinic_data$type)
design <- model.matrix(~0+factor(group_list))
colnames(design)=levels(factor(group_list))
rownames(design)=colnames(exprdata)


contrast.matrix<-makeContrasts(paste0(unique(group_list),
                                      collapse = "-"),levels = design)
fit <- lmFit(exprdata,design)

fit2 <- contrasts.fit(fit, contrast.matrix) 
fit2 <- eBayes(fit2, trend = T)

mrna_DEG = topTable(fit2, coef=1, n=Inf) %>% 
  dplyr::select(baseMean = AveExpr,
                log2FoldChange = logFC,
                lfcSE = B,
                stat = t,
                pvalue = P.Value,
                padj = adj.P.Val)
save(mrna_DEG, file = "result/GEO-GSE76260_miRNA_DEG.RData")
####

pfit1 <- EnhancedVolcano(mrna_DEG,
                         
                         lab = rownames(mrna_DEG),
                         
                         selectLab = c("a"),
                         x = "log2FoldChange",
                         
                         y = "pvalue",
                         
                         xlab = bquote(~Log[2]~ "fold change"),
                         
                         ylab = bquote(~-Log[10]~"P-value"),
                         
                         pCutoff = 0.05,
                         
                         FCcutoff = 0.26,
                         
                         xlim=c(-2,2),
                         ylim = c(0,10),
                         
                         transcriptLabSize = 0,
                         
                         colAlpha = 1,
                         
                         legend=c("NS","Log2 FC","P-value",
                                  "P-value & Log2 FC"),
                         
                         legendPosition = "top",
                         transcriptPointSize = 1.2,
                         
                         legendLabSize = 10,
                         
                         legendIconSize = 3.0,
                         title = "GSE76260 Tumor versus Normal")

ggpubr::ggexport(plot = pfit1, 
                 filename = paste("GSE76260", ".tiff", sep = ""), 
                 device = "tiff", 
                 res = 300,
                 width = 1800, 
                 height = 2200)



