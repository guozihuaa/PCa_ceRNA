rm(list = ls())

library(GDCRNATools)
library(stringr)
library(dplyr)
library(annotables)
library(matrixStats)
####### ceRNA network analysis #######
load("rawdata/miRNA_lncRNA_inter.RData")
load("rawdata/miRNA_mRNA_inte.RData")
load("PRAD_RNA_and_clinical.RData")
load("rawdata/mirexpr.RData")

load("result/intersect_Gene.RData")
deggene <- union(intersetup, intersectdown)

load("result/intersect_miRNA.RData")
degmirna <- union(intersetup, intersectdown)

rnaExpr <- gdcVoomNormalization(counts = lihc_ranexp2, filter = FALSE)

#################
load("rawdata/mibaseIndex.RData")

newmirbase <- dplyr::filter(mibaseindex, names %in% rownames(mirExpr))

######
data("grch37")

annotiondata <- dplyr::select(grch37, ensgene, symbol)

rnaExpr2 <- rnaExpr[deggene, ] %>% 
            as.data.frame() %>% 
            rownames_to_column("IDs") %>% 
            inner_join(annotiondata, by =  c("IDs" = "ensgene"))

rnaExpr2 <- dplyr::select(rnaExpr2, -IDs) %>% 
            group_by(symbol) %>% 
            summarise_all(mean)

degmirna2 <- dplyr::filter(newmirbase, ids %in% degmirna)

miRNA3 <- mirExpr[as.character(degmirna2$names), ]

miRNA.mrna.interactions <- filter(mirna_mran, miRNA %in% rownames(miRNA3))
miRNA.lncRNA.interactions <- filter(all_mir_lncRNA, mirna %in% rownames(miRNA3))

#############


candidate_nodelnc <- dplyr::filter(grch37, 
                                              ensgene %in% deggene, 
                                              biotype != "protein_coding")

candidate_nodelncup <- dplyr::filter(grch37, 
                                   ensgene %in% intersetup, 
                                   biotype != "protein_coding")

candidate_nodelnc_do <- dplyr::filter(grch37, 
                                   ensgene %in% intersectdown, 
                                   biotype != "protein_coding")

candidate_neodegene <- dplyr::filter(grch37, 
                                     ensgene %in% deggene, 
                                     biotype == "protein_coding")

LNCtarget <- filter(miRNA.lncRNA.interactions, 
                                     target %in% candidate_nodelnc$symbol)

mRNAtarget <- filter(miRNA.mrna.interactions, 
                     Genesymbol %in% candidate_neodegene$symbol)

X<-split(LNCtarget, LNCtarget$target)
miRNA.target_lnc <- lapply(X, function(x){as.character(x$mirna)})

Y<-split(mRNAtarget, mRNAtarget$Genesymbol)
miRNA.target_mrna <- lapply(Y, function(x){as.character(x$miRNA)})

finalids <- intersect(c(names(miRNA.target_lnc), names(miRNA.target_mrna)), rnaExpr2$symbol)

rnaexp3 <- dplyr::filter(rnaExpr2, symbol %in% finalids) %>% 
           column_to_rownames("symbol") %>% 
           as.matrix()



rnaexp4 <- rnaexp3[rowSds(rnaexp3)!=0, ]

deLNC <- intersect(names(miRNA.target_lnc), rownames(rnaexp4))
dePC <- intersect(names(miRNA.target_mrna), rownames(rnaexp4))

colnames(rnaexp4) <- str_sub(colnames(rnaexp4), 1, -14)

sampleids <- intersect(colnames(rnaexp4), colnames(miRNA3))

rnaexp4 <- rnaexp4[,  sampleids]
miRNA3 <- miRNA3[,  sampleids]

exprofile <- t(rbind(miRNA3, rnaexp4)) %>% 
             data.frame(check.names = F) %>% 
             rownames_to_column("sampleID")

save(exprofile, file = "result/exprofile.RData")

ceOutput <- gdcCEAnalysis(lnc       = deLNC, 
                          pc          = dePC, 
                          lnc.targets = miRNA.target_lnc, 
                          pc.targets  = miRNA.target_mrna, 
                          rna.expr    = rnaexp4, 
                          mir.expr    = miRNA3)

ceOutput2 <- ceOutput[ceOutput$hyperPValue<0.05 & ceOutput$corPValue<0.01 & ceOutput$regSim > 0.5 & ceOutput$sppc > 0.01, ]

ceOutput2 <- ceOutput[ceOutput$hyperPValue<0.05 & ceOutput$corPValue<0.05 & ceOutput$regSim != 0,]

write.csv(ceOutput2, "result/gdcceRNAout.csv", row.names = F, quote = F)

edges <- gdcExportNetwork(ceNetwork = ceOutput2, net = 'edges')
nodes <- gdcExportNetwork(ceNetwork = ceOutput2, net = 'nodes')

write.csv(edges, "result/ceRNAedges.csv", row.names = F, quote = F)
write.csv(nodes, "result/ceRNAnodes.csv", row.names = F, quote = F)
save(ceOutput, edges, nodes, file = "result/gdcceRNAout.rds")
