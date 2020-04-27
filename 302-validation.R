library()
rm(list = ls())

mydata <- read_csv("result/network/pvalueDF.csv")
load("rawdata/for_figure_plot.RData")

t_exp$sample <- pheno$sampleID
resultnode <- dplyr::filter(mydata, pvalue < 0.05, auc > 0.8, p.value < 0.05)

all_edges <- readxl::read_xlsx("result/network/network.xlsx") %>% 
  filter(fromNode %in% resultnode$name, toNode %in% resultnode$name)

finalnoudes <- c(all_edges$fromNode, all_edges$toNode) %>% 
  unique()

forplot <- dplyr::filter(mydata, name %in% finalnoudes)

joindata <- mutate(joindata, 
                   meth = ifelse(shortLetterCode == "TP", 1, 0))

finaldata <- dplyr::filter(joindata, shortLetterCode %in% c("NT", "TP"))



####

library(glmnet)
x <- finaldata[, forplot$name]
y <- finaldata$meth



x$meth <- y

k <- dim(x)[1]
predictions <- c()
for (i in 1:k) {
  model <- glm(meth ~., family=binomial, data=x[-i,])
  predictions[i] <- predict(model, x[i,], type="resp")
}

library(pROC)
roc(y, predictions)


tiff("PlotAUC_ALL.tiff", width = 1800, height = 1500, units = 'px', res = 300)
roc1 <- roc(y, predictions)   
plot(roc1, col="red")  



legend(0.65, 0.3, 
       c("core-ceRNA (AUC=0.96)"), 
       col = c( "red"),
       text.col = c("red"), 
       lty = 1, 
       lwd = 2,
       pch = NA,
       merge = TRUE, bg = "gray90",
       cex = 1 )

dev.off()

model <- glm(meth ~., family=binomial, data=x)

result <- summary(model)

logistic_coeef <- result$coefficients

#### COXPH risk score
library(survival)
raw_clinic_data <- pheno
tmpsigure <- bind_cols(raw_clinic_data[, 3:4], t_exp[forplot$name]) %>% 
  na.omit()

res.cox <- coxph(Surv(DFDtime, DFSstatus) ~ ., 
                 data =  tmpsigure)


### risk score
pheno <- data.frame(pheno) 

riskscore <- c(as.matrix(t_exp[forplot$name]) %*% matrix(forplot$coef, ncol = 1))

pheno$riskscore <- (riskscore-min(riskscore)) / (max(riskscore)- min(riskscore))

pheno$class <- ifelse(pheno$riskscore > median(pheno$riskscore), "High", "Low")


pheno <- pheno %>% 
  dplyr::arrange(riskscore) %>% 
  mutate(ids = 1:length(riskscore))

p <- ggplot(pheno, aes(x = ids, y=riskscore, color=class)) +
  geom_point(size=1)+
  theme_bw()+
  #guides(fill = guide_legend(title = "LEFT"))+
  #guides(fill=FALSE)+
  scale_color_manual(values=c('red','blue'))+
  scale_x_continuous("Samples",limits=c(0,500))+
  ylab("Normalized risk score")+
  #
  geom_vline(aes(xintercept=245.5), colour="gray20", linetype="dashed",size=1)+
  theme(panel.border = element_rect(size = 1.5,fill = NA),
        # legend.position=c(0.65, 0.8),
        # legend.box.background = element_rect(),
        # #legend.key.size = unit(10,"point"),
        # legend.box.margin = margin(1, 0.3, 0.6, 0.6),
        # legend.text = element_text(size = 13, face="bold", colour="black"),
        # legend.title = element_text(face = "bold", size = 15),
        # legend.title.align = 0.5,
        #legend.spacing = 0.1,
        legend.position = "none",
        axis.text.x = element_text(face = "bold", size = 12, colour = "black"),
        axis.text.y = element_text(face = "bold", size = 12,colour = "black"),
        axis.title = element_text(face = "bold", size = 15)
  )

ggsave("result/Figures/risk.tiff", p, width = 7, height = 4, dpi = 300)

pheno2 <- dplyr::select(pheno, DFDtime, status = DFSstatus, ids) %>% 
          na.omit()

pheno2$status <- as.factor(pheno2$status)

p2 <- ggplot(pheno2, aes(x=ids, y=DFDtime)) + 
  geom_point(aes(col = status, shape=status, size = status)) +   # draw points
  xlim(c(0, 500)) +    # draw smoothing line
  scale_y_continuous("Disease free time (months)",limits=c(0,160))+
  xlab("Samples")+
  theme_bw()+
  #guides(fill=FALSE)+
  scale_color_manual(values=c('gray70','red'))+
  scale_size_manual(values=c(2,3))+
  geom_vline(aes(xintercept=245.5), colour="gray20", linetype="dashed",size=1)+
  theme(panel.border = element_rect(size = 1.5,fill = NA),
        legend.position="none",
        axis.text.x = element_text(face = "bold", size = 12, colour = "black"),
        axis.text.y = element_text(face = "bold", size = 12,colour = "black"),
        axis.title = element_text(face = "bold", size = 15)
  )

ggsave("result/Figures/sample.tiff", p2, width = 7, height = 4, dpi = 300)

######## 热图

pradsubstype <- TCGAquery_subtype(tumor = "PRAD") %>% 
                dplyr::select(sample, Subtype)

t_exp2 <- column_to_rownames(t_exp, "sample") 

t_exp3 <- t_exp2[pheno$sampleID, forplot$name] %>% 
          as.matrix() 


fulldata <- scale(t_exp3,
                    center=TRUE, 
                    scale=TRUE) %>% 
         t()


fulldata <- limitRange(fulldata)

library(pheatmap)


mycolor <- c(rgb(0,0,1,1), rgb(1,1,1), rgb(1,0,0,1))


annotation_col <- left_join(pheno, pradsubstype, by = c("sampleID" = "sample")) %>% 
                  column_to_rownames( "sampleID") %>% 
                  dplyr::select(class,  Tstage, Nstage, Subtype)

annotation_col$Subtypes <- as.character(annotation_col$Subtypes)
names(annotation_col) <- c("Risk score", "T-stage", "N-stage", "Subtypes")

annotation_col$`T-stage`[is.na(annotation_col$`T-stage`)] <- "T1"
annotation_col$`N-stage`[is.na(annotation_col$`N-stage`)] <- "Unknown"
annotation_col$Subtypes[is.na(annotation_col$Subtypes)] <- "Unknown"


tiff(file = "result/heatmap3.tiff",
     width = 3500,
     height = 4000,
     units = "px",
     res = 600)

pheatmap(fulldata, 
         cluster_rows = F, 
         cluster_cols = F, 
         show_colnames = F,
         annotation_col = annotation_col,
         border_color=NA,
         color = colorRampPalette(mycolor)(20),
         gaps_col = 246)

dev.off()
######

count <- as.matrix(table(annotation_col$`Risk score`, annotation_col$`T-stage`))

matdata <- as.matrix(unclass(count))

fisher.test(matdata, simulate.p.value = T)

count <- as.matrix(table(annotation_col$`Risk score`, annotation_col$`N-stage`))

matdata <- as.matrix(unclass(count))

fisher.test(matdata, simulate.p.value = T)


count <- as.matrix(table(annotation_col$`Risk score`, annotation_col$Subtypes))

matdata <- as.matrix(unclass(count))

fisher.test(matdata, simulate.p.value = T)


######### 

library(ConsensusClusterPlus)

dt = as.dist(1-cor(t(t_exp3), method="pearson"))



result <- ConsensusClusterPlus(dt,
                               maxK=6,
                               reps=100,
                               pItem=0.8,
                               pFeature=1,
                               title="example2",
                               distance="pearson",
                               clusterAlg="pam",
                               seed=1234,
                               plot="png")

library(cluster)
SIL = silhouette(result[[7]]$consensusClass, dt)


newclustr <- data.frame(ids = names(result[[4]]$consensusClass),
                        cluster = unlist(result[[4]]$consensusClass))


pheno2 <- inner_join(pheno, newclustr, by = c("sampleID" = "ids")) 

library(survival)
library(survminer)
fit <- survfit(Surv(DFDtime, DFSstatus) ~ cluster, 
               data = pheno2)

pfit1 <- ggsurvplot(fit,
                    palette = rainbow(4),
                    data = pheno2, 
                    censor.shape="|", 
                    pval = TRUE,
                    censor.size = 3,
                    break.time.by = 30,
                    
                    #legend = "left",
                    #legend.labs = c("Advanced", "Early"),
                    xlab = "Time in months", 
                    ylab = "Disease free probability",
                    
                    legend.labs = c("cluster1", "cluster2", "cluster3", "cluster4"),
                    pval.size = 6,
                    pval.coord = c(0,0.1),
                    ggtheme = theme_bw(),
                    
                    font.legend = 16,
                    font.x =  c(17, "bold", "black"),
                    font.y = c(17, "bold", "black"),
                    font.tickslab = c(15, "bold", "black"),
                    risk.table = TRUE,
                    risk.table.col = "strata",
                    risk.table.height = 0.25
                    #fontsize = 10
)

pfit1$table <- ggpar(
  pfit1$table,
  font.title    = c(13, "bold", "black"),
  #font.subtitle = c(15, "bold", "pink"),
  #font.caption  = c(11, "plain", "darkgreen"),
  font.x        = c(17, "bold", "black"),
  font.y        = c(17, "bold", "black"),
  font.xtickslab = c(15, "bold", "black"),
  font.ytickslab = 13)


ggpubr::ggexport(plot = pfit1, 
                 filename = "result/Figures/substypeos.tiff", 
                 device = "tiff", 
                 res = 300,
                 width = 2000, 
                 height = 2800)



######### PCA& tsne

###### 模式


library(ggplot2)
library(readr)
library(Rtsne)

load("result/data/2-pro/minrna_mrna_lnc_exp.rds")

bind_exp <- t_exp3

bind_exp <- t(scale(t(bind_exp)))

pca.norm = prcomp(bind_exp)

summary(pca.norm,loadings=TRUE)


mydata <- as.data.frame(pca.norm$x) %>% 
  dplyr::select(PC1, PC2) %>% 
  rownames_to_column("sample") %>% 
  inner_join(newclustr, by = c("sample" = "ids"))

mydata$cluster <- as.factor(mydata$cluster)


pstage <- ggplot(mydata, aes(x=PC1, y=PC2, color=cluster)) +
  geom_point(size=3) +
  labs(colour = "Cluster")+
  guides(colour = guide_legend(override.aes = list(size=6))) +
  xlab("PC1 - var: 44%") + ylab("PC2 - var: 23%") +
  ggtitle("Principal component analysis") +
  scale_color_manual(values = c("#FF0000FF","#80FF00FF", "#00FFFFFF", "#8000FFFF"))+
  theme_bw()+
  theme(panel.border = element_rect(size = 1.5,fill = NA),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 15, colour = "black"),
        legend.title = element_text(face = "bold", size = 12, colour = "black"),
        
        #legend.position="none",
        axis.text.x = element_text(face = "bold", size = 12, colour = "black"),
        axis.text.y = element_text(face = "bold", size = 12,colour = "black"),
        axis.title = element_text(face = "bold", size = 15)
  )

ggsave(plot = pstage, 
       filename = "result/Figures/00.tiff", 
       device = "tiff", 
       dpi = 300, 
       width = 20, 
       height = 19,
       units = "cm")


##### tsne
library(Rtsne)

plot(mydata)

tsne <- Rtsne(as.matrix(bind_exp), check_duplicates = FALSE, pca = TRUE, 
              perplexity=30, theta=0.1, dims=2,
              max_iter = 2000)

embedding <- as.data.frame(tsne$Y)
embedding$sample <- rownames(bind_exp)

embedding <- inner_join(embedding, newclustr, by = c("sample" = "ids"))

embedding$cluster <- as.factor(embedding$cluster)


p <- ggplot(embedding, aes(x=V1, y=V2, color=cluster)) +
  geom_point(size=1.25) +
  guides(colour = guide_legend(override.aes = list(size=6))) +
  xlab("") + ylab("") +
  ggtitle("t-SNE 2D Embedding of Products Data") +
  theme_light(base_size=20) +
  theme(strip.background = element_blank(),
        strip.text.x     = element_blank(),
        axis.text.x      = element_blank(),
        axis.text.y      = element_blank(),
        axis.ticks       = element_blank(),
        axis.line        = element_blank(),
        panel.border     = element_blank())

p

ggsave("tsne.png", p, width=8, height=6, units="in")



