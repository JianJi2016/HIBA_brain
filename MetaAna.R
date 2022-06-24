# library(readr)
# library(dplyr)
# library(textshape)
# library(DMwR)
# library(ropls)
# library(ggpubr)
# library(ggplot2)
# library(ggrepel)
# library(tibble)
# library(stats)
# library(ggsignif)
# library(textshape)
# library(pROC)
# library(plotROC)
# library(ComplexHeatmap)
# library(circlize)
pacman::p_load(readr,dplyr,textshape,DMwR,ropls,ggpubr,ggplot2,ggrepel,tibble,
               stats,ggsignif,pROC,plotROC,ComplexHeatmap,circlize)

file_name <- c("2020SAMPLE")
sample_import <- read_csv(paste0(file_name,".csv"))

# sample_import %>% subset(time == "24h") %>%
#   select(-time) %>%
#   subset(dose != "0 ppm") ->sample_import




colnames(sample_import)[1:2] <- c("Sample","Treatment")
sample_import$Treatment <- factor(sample_import$Treatment,
                                  levels = unique(sample_import$Treatment))


sampleinfo <- sample_import %>% 
  select(.,Sample, Treatment) %>%
  textshape::column_to_rownames(.)


metainfo <- colnames(sample_import)[-c(1:2)] %>% 
  as.data.frame() %>%
  rename(., Compounds = .)

compound_name <- colnames(sample_import)[-c(1:2)]

matrix <- sample_import %>% 
  select(.,-Sample, -Treatment)

write.csv(sample_import, "[1] data_original.csv", row.names = F)
write.csv(sampleinfo, "[2] sampleinfo.csv", row.names = F)
write.csv(metainfo, "[3] metainfo.csv", row.names = F)
write.csv(matrix, "[4] matrix.csv", row.names = F)

matrix_numeric <- do.call("cbind",lapply(matrix,as.numeric))

#  0 finding and knn filled
for (i in 1:length(compound_name)) {
  if(grep("TRUE",matrix_numeric[,i] == 0) %>% length() > 0){
    matrix_numeric[,i][grep("TRUE",matrix_numeric[,i] == 0)] <- NA
  }
}
matrix_fill <- knnImputation(matrix_numeric, k = 3)

# 
# # outlier finding and knn filled
# outlier.IQR <- function(x, multiple = 10, replace = TRUE, revalue = NA) { 
#   q <- quantile(x, na.rm = TRUE) 
#   IQR <- q[4] - q[2]
#   x1 <- which(x < q[2] - multiple * IQR | x > q[4] + multiple * IQR)
#   x2 <- x[x1]
#   if(length(x2) > 0){
#       outlier <- data.frame(location = x1, value = x2)}else{
#       outlier <- data.frame(location = 0, value = 0)}
#   if (replace == TRUE){
#       x[x1] <- revalue
#   }
#   return(x)
# }
# 
# sample_import$Treatment %>% unique() %>% as.character() -> group_info
# Group1_order <- grep(group_info[1],sample_import$Treatment)
# Group2_order <- grep(group_info[2],sample_import$Treatment)
# 
# for (i in 1:ncol(matrix_fill)) {
#   matrix_fill[,i] <- c(
#     outlier.IQR(matrix_fill[Group1_order,i]),
#     outlier.IQR(matrix_fill[Group2_order,i])
#   )
# }
#   
# matrix_fill <- knnImputation(matrix_fill, k = 3)


data_gapfilled <- cbind(sample_import[,c(1:2)],matrix_fill)
write.csv(data_gapfilled, "[5] data_gapfilled.csv", row.names = F)

# log transformation
matrix_log <- matrix_fill %>% log10

data_log <- matrix_log %>%
  cbind(sample_import[,c(1:2)],.)
write.csv(data_log, "[6] data_log.csv", row.names = F)

pareto_scale <- function(x){(x-mean(x, na.rm = T))/sqrt(sd(x, na.rm = T))}
auto_scale <- function(x){(x-mean(x, na.rm = T))/sd(x, na.rm = T)}

marix_pareto  <- apply(matrix_log, 2, pareto_scale)
data_pareto <- cbind(sample_import[,c(1:2)],marix_pareto)
write.csv(data_pareto, "[7] data_pareto.csv", row.names = F)


rownames(marix_pareto) <- rownames(sampleinfo)
PCA_file <- list(sampleinfo,metainfo,marix_pareto)
names(PCA_file) <- c("sampleinfo", "metainfo", "marix_pareto")


dir.create("PCA Results")
pca.path <- paste0(getwd(),"/PCA Results")


for(i in 1:ncol(PCA_file$marix_pareto)){
  a=mean(na.omit(PCA_file$marix_pareto[i,]))
  PCA_file$marix_pareto[i,is.na(PCA_file$marix_pareto[i,])] = a
}


pca_result <- opls(x = PCA_file$marix_pareto, predI = 2)
pac_modelscore <- pca_result@modelDF
pac_modelscore$PC = rownames(pac_modelscore)
write.csv(pac_modelscore, paste0(pca.path,"/PCA Model Score.csv"),row.names = F)

# PCA model score
PCA_modelscore <- ggbarplot(pac_modelscore, x = "PC", y = "R2X",
          fill = "PC",palette = "aaas",
          width = 0.5,
          font.main = c(16, "plain", "black"),            
          font.x = c(16, "plain", "black"),                  
          font.y = c(16, "plain", "black"),                
          font.legend = c(16, "plain", "black"),
          label = "R2X",
          xlab = "Components",
          ylab = "R2X Value",
          legend = "top",
          legend.title = "Components") +
  theme(axis.text = element_text(size = 16))
ggsave("PCA_modelscore.png", height = 4.5, width = 5,path = pca.path)

pdf(paste0(pca.path,"/PCA Model Score.pdf"),height = 4.5, width = 5)
PCA_modelscore
dev.off()


# PCA score plot
scoreMN <- pca_result@scoreMN
scoreMN <- cbind(scoreMN, PCA_file$sampleinfo)
scoreMN$samples <- rownames(scoreMN)
write.csv(scoreMN,paste0(pca.path,"/PCA ScoreMN.csv"))

# scoreMN$Treatment = factor(scoreMN$Treatment, 
#                                  levels = c(scoreMN$Treatment %>% 
#                                               unique() %>% .[1] %>% as.character(),
#                                             scoreMN$Treatment %>% 
#                                               unique() %>% .[2] %>% as.character()))

PCA_score <- ggscatter(scoreMN, x = "p1", y = "p2",
          color = "black", 
          fill = "Treatment",
          shape = 21,
          ellipse = T, 
          size = 5,
          ellipse.level = 0.95,
          mean.point = F,
          alpha = 0.8,
          font.label = 10, repel = TRUE,
          xlab = paste0("PC1 (",pca_result@modelDF$R2X[1]*100,"%)"),
          ylab = paste0("PC2 (",pca_result@modelDF$R2X[2]*100,"%)"),
          font.main = c(16, "plain", "black"),            
          font.x = c(16, "plain", "black"),                  
          font.y = c(16, "plain", "black"),                
          legend = "top",  
          legend.title = "Treatment",  
          font.legend = c(16, "plain", "black"),               
          rotate = F,                                         
          ticks = T,   
          label = "samples",
          tickslab = T,  
          palette = "aaas") +
  theme(axis.text = element_text(size = 16))
ggsave("PCA score.png",height =4.5, width = 5,path = pca.path)

pdf(paste0(pca.path,"/PCA score.pdf"),height = 4.5, width = 5)
PCA_score
dev.off()

# PCA loading plot
loadingMN <- pca_result@loadingMN
loadingMN <- as.data.frame(loadingMN)
loadingMN$compound <- as.vector(rownames(loadingMN))
colnames(loadingMN)[1:2] = c("PC1","PC2")

loadingMN_sort <- loadingMN %>% 
  mutate(distance = sqrt(PC1^2+PC2^2)) %>%
  arrange(desc(distance))

loadingMN_label <- loadingMN_sort %>%
  mutate(label = ifelse(distance > loadingMN_sort$distance[9], compound, ""))
write.csv(loadingMN_label, paste0(pca.path,"/PCA Loading.csv"))

PCA_loading <- ggscatter(loadingMN_label, x = "PC1", y = "PC2",
          color = "black",
          fill = "distance",
          size =  4,
          alpha = 0.8,
          shape = 21,
          # label = "label",
          legend = " ",
          font.main = c(16, "plain", "black"),            
          font.x = c(16, "plain", "black"),                  
          font.y = c(16, "plain", "black"),                
          font.legend = c(16, "plain", "black"),
          xlab = "Loading 1",
          ylab = "Loading 2") +
  scale_fill_gradient(low = "grey", high = "red",na.value = NA) +
  # scale_colour_gradient2() +
  geom_text_repel(data = loadingMN_label, size = 4,aes(label = label))+
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(axis.text = element_text(size = 16)) 

ggsave("PCA loading.png",height =4.5, width = 5,path = pca.path)

pdf(paste0(pca.path,"/PCA loading.pdf"),height = 4.5, width = 5)
PCA_loading
dev.off()



PCA_loading2 <- loadingMN_label %>% mutate(label2 <- ifelse(label == "","grey","red")) %>% 
  rename(.,lable2 = "label2 <- ifelse(label == \"\", \"grey\", \"red\")") %>% 
  mutate(marks = as.factor(lable2)) %>%
  ggscatter(., x = "PC1", y = "PC2",
            color = "marks",
            # fill = "distance",
            size =  4,
            alpha = 0.5,
            # shape = 21,
            # label = "label",
            legend = "",
            font.main = c(16, "plain", "black"),            
            font.x = c(16, "plain", "black"),                  
            font.y = c(16, "plain", "black"),                
            font.legend = c(16, "plain", "black"),
            xlab = "Loading 1",
            ylab = "Loading 2") +
  scale_color_manual(values= c("#E1E2E4", "red")) +
  geom_text_repel(data = loadingMN_label, size = 4,aes(label = label))+
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(axis.text = element_text(size = 16))

ggsave("PCA loading2.png",height =4.5, width = 5,path = pca.path)

pdf(paste0(pca.path,"/PCA loading2.pdf"),height = 4.5, width = 5)
PCA_loading2
dev.off()

# PLS-DA
dir.create("PLS-DA Results")
plsda.path <- paste0(getwd(),"/PLS-DA Results")

plsda_result <- opls(x = PCA_file$marix_pareto, 
                     y = PCA_file$sampleinfo[, 'Treatment'], 
                     predI = 2,
                     orthoI = NA,
                     crossvalI = 7,
                     permI = 1000
                     )

write.csv(plsda_result@modelDF,paste0(plsda.path,"/PLS-DA.ModelScore.csv"))

plsda_modelscore <- plsda_result@modelDF
plsda_plot <- plsda_modelscore %>% 
  tibble::rownames_to_column(var = "group") %>%
  select(group,`R2Y(cum)`,`Q2(cum)`) %>%
  reshape2::melt(.,1)

PLSDA_modelscore <- ggbarplot(plsda_plot,x= "group",y="value",
          fill = "variable",palette = "aaas",
          color = "black",
          font.main = c(16, "plain", "black"),            
          font.x = c(16, "plain", "black"),                  
          font.y = c(16, "plain", "black"),                
          font.legend = c(16, "plain", "black"),
          label = " ",
          xlab = "Components",
          ylab = "Score",
          # legend = "right",
          legend.title = " ",
          position = position_dodge(0.7))+
  theme(axis.text = element_text(size = 16))

ggsave("PLS-DA Model Score.png",height =4.5, width = 5,path = plsda.path)

pdf(paste0(plsda.path,"/PLS-DA Model Score.pdf"),height = 4.5, width = 5)
PLSDA_modelscore
dev.off()

# PLSDA score plot
plsda_scoreMN <- plsda_result@scoreMN
plsda_scoreMN <- cbind(plsda_scoreMN, PCA_file$sampleinfo)
plsda_scoreMN$samples <- rownames(plsda_scoreMN)
write.csv(plsda_scoreMN,paste0(plsda.path,"/PLS-DA score.csv"))

plsda_scoreMN$Treatment = factor(plsda_scoreMN$Treatment, 
                                 levels = c(plsda_scoreMN$Treatment %>% 
                                              unique() %>% .[1] %>% as.character(),
                                            plsda_scoreMN$Treatment %>% 
                                              unique() %>% .[2] %>% as.character()))

PLSDA_score <- ggscatter(plsda_scoreMN, x = "p1", y = "p2",
          color = "black", 
          fill = "Treatment",
          shape = 21,
          ellipse = T, 
          size = 5,
          ellipse.level = 0.95,
          mean.point = F,
          alpha = 0.8,
          font.label = 16, repel = TRUE,
          xlab = paste0("PC1 (",plsda_result@modelDF$R2X[1]*100,"%)"),
          ylab = paste0("PC2 (",plsda_result@modelDF$R2X[2]*100,"%)"),
          font.main = c(16, "plain", "black"),            
          font.x = c(16, "plain", "black"),                  
          font.y = c(16, "plain", "black"),                
          legend = "top",  
          legend.title = "Treatment",  
          font.legend = c(16, "plain", "black"),               
          rotate = F,                                         
          ticks = T,   
          label = "samples",
          tickslab = T,  
          palette = "aaas") +
  theme(axis.text = element_text(size = 16))

ggsave("PLS-DA score.png",height = 4.5, width = 5,path = plsda.path)

pdf(paste0(plsda.path,"/PLS-DA score.pdf"),height = 4.5, width = 5)
PLSDA_score
dev.off()

plsda_loadingMN <- plsda_result@loadingMN
plsda_loadingMN <- as.data.frame(plsda_loadingMN)
plsda_loadingMN$compound <- as.vector(rownames(plsda_loadingMN))

colnames(plsda_loadingMN)[1:2] = c("PC1","PC2")


plsda_loadingMN_sort <- plsda_loadingMN %>% 
  mutate(distance = PC1^2+PC2^2) %>%
  arrange(desc(distance))

plsda_loadingMN_label <- plsda_loadingMN_sort %>%
  mutate(label = ifelse(distance > plsda_loadingMN_sort$distance[9], compound, ""))
write.csv(plsda_loadingMN_label, paste0(plsda.path,"/PLSDA Loading.csv"))

PLSDA_loading <- ggscatter(plsda_loadingMN_label, x = "PC1", y = "PC2",
                         color = "black",
                         fill = "distance",
                         size =  4,
                         alpha = 0.5,
                         shape = 21,
                         # label = "label",
                         legend = " ",
                         font.main = c(16, "plain", "black"),            
                         font.x = c(16, "plain", "black"),                  
                         font.y = c(16, "plain", "black"),                
                         font.legend = c(16, "plain", "black"),
                         xlab = "Loading 1",
                         ylab = "Loading 2") +
  scale_fill_gradient(low = "grey", high = "red") +
  geom_text_repel(data = loadingMN_label, size = 4,aes(label = label))+
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(axis.text = element_text(size = 16))

ggsave("PLSDA loading.png",height =4.5, width = 5,path = plsda.path)

pdf(paste0(plsda.path,"/PLSDA loading.pdf"),height = 4.5, width = 5)
PLSDA_loading
dev.off()



PLSDA_loading2 <- plsda_loadingMN_label %>% mutate(label2 <- ifelse(label == "","grey","red")) %>% 
  rename(.,lable2 = "label2 <- ifelse(label == \"\", \"grey\", \"red\")") %>% 
  mutate(marks = as.factor(lable2)) %>%
  ggscatter(., x = "PC1", y = "PC2",
            color = "marks",
            # fill = "distance",
            size =  4,
            alpha = 0.8,
            # shape = 21,
            # label = "label",
            legend = "",
            font.main = c(16, "plain", "black"),            
            font.x = c(16, "plain", "black"),                  
            font.y = c(16, "plain", "black"),                
            font.legend = c(16, "plain", "black"),
            xlab = "Loading 1",
            ylab = "Loading 2") +
  scale_color_manual(values= c("#E1E2E4", "red")) +
  geom_text_repel(data = plsda_loadingMN_label, size = 4,aes(label = label))+
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(axis.text = element_text(size = 16))

ggsave("PLSDA loading2.png",height =4.5, width = 5,path = plsda.path)

pdf(paste0(plsda.path,"/PLSDA loading2.pdf"),height = 4.5, width = 5)
PLSDA_loading2
dev.off()

# VIP
plsda_VIP <- getVipVn(plsda_result) %>% 
  as.data.frame() %>%
  tibble::rownames_to_column(var = "name") %>%
  rename(VIP =".")
write.csv(plsda_VIP,paste0(plsda.path,"/PLS-DA VIP.csv"))

plsda_VIPplot <- plsda_VIP %>% arrange(desc(VIP)) %>%
  slice(1:10)

PLSDA_VIP <- ggbarplot(plsda_VIPplot, x = "name", y = "VIP",
            xlab = "",
            ylab = "VIP Value",
            fill = "VIP",
            legend = "",
            font.main = c(16, "plain", "black"),            
            font.x = c(16, "plain", "black"),                  
            font.y = c(16, "plain", "black"),                
            font.legend = c(16, "plain", "black"),
            orientation = "horiz",
            order = rev(plsda_VIPplot$name)) +
  scale_fill_gradient(low = "grey", high = "red") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_y_continuous(expand = c(0,0)) +
  theme(axis.text = element_text(size = 16))
  
ggsave("PLSDA VIP.png",height =4.5, width = 6,path = plsda.path)

pdf(paste0(plsda.path,"/PLSDA VIP.pdf"),height = 4.5, width = 6)
PLSDA_VIP
dev.off()

# permutation
permutation <- plsda_result@suppLs$permMN
permutation <- as.data.frame(permutation)
colnames(permutation)[2:3] = c("R2Y","Q2Y")
permutation <- permutation[,c(2,3,7)]
permutation2 <- reshape2::melt(permutation,3)

R2Y.value <- permutation2$value[permutation2$variable == "R2Y"]
Q2Y.value <- permutation2$value[permutation2$variable == "Q2Y"]
XR1 = 1
YR1 = max(R2Y.value)
XR2 = mean(permutation2$sim[permutation2$variable == "R2Y"])
YR2 = mean(R2Y.value)
slop1 <- (YR1-YR2)/(XR1-XR2)
Intecep1 <- YR1-slop1*XR1

XQ1 = 1
YQ1 = max(Q2Y.value)
XQ2 = mean(permutation2$sim[permutation2$variable == "Q2Y"])
YQ2 =mean(Q2Y.value)
slop2 <- (YQ1-YQ2)/(XQ1-XQ2)
Intecep2 <- YQ1-slop2*XQ1

PLSDA_permutation <- ggscatter(permutation2,x= "sim", y = "value",
          color = "black",
          fill = "variable",
          xlab = "Permutation correlation",
          ylab = "R2 and Q2",
          shape = 21,
          size = 5,
          palette = "aaas",
          alpha = 0.8,
          legend.title=" ",
          font.main = c(16, "plain", "black"),            
          font.x = c(16, "plain", "black"),                  
          font.y = c(16, "plain", "black"),                
          font.legend = c(16, "plain", "black")
) +
  geom_hline(yintercept = max(permutation$R2Y), linetype = "dashed") +
  geom_hline(yintercept = max(permutation$Q2Y), linetype = "dashed") +
  geom_abline(intercept = Intecep1, slope = slop1, color="#3B4992",
              size=0.5) +
  geom_abline(intercept = Intecep2, slope = slop2, color="#EE0000",
              size=0.5)+
  theme(axis.text = element_text(size = 16))

ggsave("PLS-DA permutation.png",height =4.5, width = 5,path = plsda.path)

pdf(paste0(plsda.path,"/PLS-DA permutation.pdf"),height = 4.5, width = 5)
PLSDA_permutation
dev.off()



# Splot
s <- as.matrix(PCA_file$marix_pareto)
T <- as.matrix(plsda_result@scoreMN)
p1 <- NULL
for (i in 1:ncol(s)) {
  p1[i] <- as.matrix(cov(s[,i], T),ncol =1)
}

pcorr1 <-  NULL
for (i in 1:length(p1)) {
  den <- apply(T, 2, sd)*sd(s[,i])
  pcorr1[i] <- p1[i]/den
}

Splotdata <- data.frame(p1 = p1,
                        pcorr1 =pcorr1,
                        name = compound_name)

Splotdata %>% mutate(distance = p1^2+pcorr1^2) %>% 
  arrange(desc(distance)) -> Splot_distance

Splot_distance$label[1:10] <- as.character(Splot_distance$name[1:10])
Splot_distance$label[11:nrow(Splot_distance)] <- ""

Splot_distance$color_label <- abs(Splot_distance$pcorr1)
Splot_distance$size_label <- abs(Splot_distance$p1)

PLSDA_Splot <- ggscatter(Splot_distance,x= "p1", y = "pcorr1",
          fill = "color_label",
          color = "black",
          shape = 21,
          legend = "top",
          size = "size_label",
          alpha = 0.5,
          xlab = "p[1]",
          ylab = "p(corr)[1]",
          font.main = c(16, "plain", "black"),
          font.x = c(16, "plain", "black"),
          font.y = c(16, "plain", "black"),
          font.legend = c(16, "plain", "black"))+
  scale_fill_gradient(low = "grey", high = "red") +
  geom_text_repel(data = Splot_distance, size = 3, aes(label = label)) +
  guides(fill = F, size = F) +
  theme(axis.text = element_text(size = 16)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed")

ggsave("PLS-DA Splot.png",height =4.5, width = 5, path = plsda.path)

pdf(paste0(plsda.path,"/PLS-DA Splot.pdf"),height = 4.5, width = 5)
PLSDA_Splot
dev.off()



Vplot <- cbind(plsda_result@loadingMN[,1],
               plsda_result@coefficientMN,
               plsda_result@vipVn,
               rownames(plsda_result@loadingMN))
colnames(Vplot) <- c("p1","corr","VIP","Sample")
Vplot <- as.data.frame(Vplot)
Vplot$VIP <- as.numeric(as.vector(Vplot$VIP))
Vplot$p1 <- as.numeric(as.vector(Vplot$p1))
Vplot$corr <- as.numeric(as.vector(Vplot$corr))
Vplot$label = ""
Vplot$label[match(plsda_VIPplot$name,as.character(Vplot$Sample))] =  plsda_VIPplot$name
write.csv(Vplot,paste0(plsda.path,"/PLS-DA Vplot.csv"))


PLSDA_Vplot <- ggscatter(Vplot,x= "p1", y = "VIP",
          fill = "VIP",
          color = "black",
          shape = 21,
          legend = "top",
          size = "VIP",
          alpha = 0.5,
          xlab = "Coeff",
          ylab = "VIP value",
          font.main = c(16, "plain", "black"),
          font.x = c(16, "plain", "black"),
          font.y = c(16, "plain", "black"),
          font.legend = c(16, "plain", "black"))+
  scale_fill_gradient(low = "grey", high = "red") +
  geom_text_repel(data = Vplot, size = 3, aes(label = label)) +
  guides(fill=FALSE,size = F) +
  theme(axis.text = element_text(size = 16))+
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 1, linetype = "dashed")

ggsave("PLS-DA Vplot.png",height =4.5, width = 5, path = plsda.path)

pdf(paste0(plsda.path,"/PLS-DA Vplot.pdf"),height = 4.5, width = 5)
PLSDA_Vplot
dev.off()



# volcano
dir.create("Volcano Results")
volcano.path <- paste0(getwd(),"/Volcano Results")

data_gapfilled$Treatment <- factor(data_gapfilled$Treatment,
                                   levels = data_gapfilled %>% 
                                     select(Treatment) %>% 
                                     unique() %>% 
                                     pull(1)
                                  )
factorNumber <- length(levels(data_gapfilled$Treatment))


data_original <- data_gapfilled
colnames(data_original)[-c(1,2)] <- compound_name
factors <-  data_original$Treatment %>% unique() %>% as.character()
treatments <- data_original$Treatment %>% unique() %>% length()
TreatmentA <- which(data_original$Treatment == factors[1]) ## non-cancer tissues
TreatmentB <- which(data_original$Treatment == factors[2]) ## tumor tissues
runttest <- function(data1, data2) {
  results <- t.test(data1,data2)
  foldchange <- median(data2)/median(data1)
  unlist(c(results$p.value,foldchange))
}
ttest_value <- NULL
FC_value <- NULL
for(i in 3:ncol(data_original)){
  ttest_value[i-2] = runttest(data_original[TreatmentA,i],data_original[TreatmentB,i])[1]
  FC_value[i-2] = runttest(data_original[TreatmentA,i],data_original[TreatmentB,i])[2]
}

data_original_tt_FC_result <- data.frame(Sample =colnames(data_original)[-c(1,2)], 
                                         pvalue =ttest_value,
                                         FC =FC_value, 
                                         p.adjust = p.adjust(ttest_value,method = "fdr")
                                         )
vacanno <- data_original_tt_FC_result
colnames(vacanno) = c("compound","Pvalue","FC","adjust.Pvalue")

# log the FC value
FC_log <- NULL
for (i in 1 : nrow(vacanno)){
  FC_log[i]  =  log2(vacanno$FC[i])
}
vacanno$FC_log = FC_log

# log the pvalue
pvalue_log <- NULL
for (i in 1 : nrow(vacanno)){
  pvalue_log[i]  =  -log10(vacanno$Pvalue[i])
}
vacanno$pvalue_log = pvalue_log

# find the significant p value, level set as p value < 0.1
sig_p <- NULL
for (i in 1 : nrow(vacanno)){
  if (vacanno$Pvalue [i] < 0.05){
    sig_p [i] = as.character(vacanno$compound)[i]
  }else{
    sig_p [i] = NA
  }
}
vacanno$sig_p = sig_p

# metabolites pass the FC selection, as for FC >1.5 or FC < 0.66
sig_f <- NULL
for (i in 1 : nrow(vacanno)){
  if (vacanno$FC [i] > 2 | vacanno$FC [i] < 0.5){
    sig_f [i] = as.character(vacanno$compound)[i]
  }else{
    sig_f [i] = NA
  }
}
vacanno$sig_f = sig_f

# metabolites pass the pvalue and FC
sig_both <- NULL
for (i in 1 : nrow(vacanno)){
  if ((is.na(vacanno$sig_p) | is.na(vacanno$sig_f))[i]){
    sig_both [i] = NA
  }else{
    sig_both [i] = as.character(vacanno$compound)[i]
  }
}
vacanno$sig_both = sig_both

#-- choose the first 20 to label if the number of significat metabolites over 20
sig_both_20 <- c(rep(NA, length(sig_both)))
num_sig <- length(grep("TRUE",!is.na(sig_both)))

if (num_sig >5) { 
  first_20 <- match(
    tail(
      sort(
        round(
          abs(
            log2(vacanno$FC[!is.na(sig_both)])
          ), digits=2
        )
      ),5
    ),
    round(abs(log2(vacanno$FC[!is.na(sig_both)])),digits=2)
  )
  
  first_20_FC <- log2(vacanno$FC[!is.na(sig_both)])[first_20]
  ordre_name <- match(first_20_FC, log2(vacanno$FC))
  
  for (i in ordre_name){
    sig_both_20[i] = sig_both[i]
  }
}else{sig_both_20 = sig_both}
vacanno$sig_both_20 = sig_both_20

#  set the colors, FC > 1 then blue if FC <-1 then red,just segregation of FC direction
FC_col <- NULL
for (i in 1 : nrow(vacanno)){
  if(log2(vacanno$FC[i]) > 0){  FC_col[i]  = "Red"}else{ FC_col[i] = "Blue"}
}
vacanno$FC_col = FC_col

# set the shape of metabolites whose p value over 0.1

p_shape <- NULL
for (i in 1 : nrow(vacanno)){
  if(-log10(vacanno$Pvalue[i]) > 1){ p_shape[i]  = 1}else{ p_shape[i] = 2}
}
vacanno$p_shape = p_shape
vacanno$p_shape = as.factor(vacanno$p_shape)

# set the transparancy, if metabolites did not pass the FC and p value test, 
# then their color are transparancy 
P_trans <- NULL
for (i in 1 : nrow(vacanno)){
  if((is.na(vacanno$sig_p) | is.na(vacanno$sig_f))[i]){
    P_trans[i]  = 0.9}else{ P_trans[i] = 1}
}
vacanno$P_trans = as.numeric(P_trans)
vacanno$FC_col <- as.factor(vacanno$FC_col)

write.csv(vacanno,paste0(volcano.path,"/Volcano Plot data.csv"),row.names = F)

cbind(
  plsda_VIP,
  vacanno %>%
    select(Pvalue:pvalue_log)
) %>%
  arrange(desc(VIP)) %>%
  mutate(sig_Metabolites = ifelse(VIP > 1 & Pvalue < 0.05,
                                  name,"")) -> 
  write.csv(.,paste0(volcano.path,"/Significant Metabolites.csv"), row.names = F)


Volcano <- ggpar(ggscatter(vacanno,x = "FC_log", 
                y = "pvalue_log",
                main = "Vacanno",
                size = "pvalue_log",
                ylab = "-log10(p)",
                xlab = "log2(FC)",
                fill = "FC_col",
                color = "black",
                shape = 21,
                font.x = c(16),
                font.y = c(16),
                font.label = c(16,"black"),
                alpha = "P_trans",
                repel = T,
                palette = "aaas") + 
        geom_vline(xintercept = c(a <- c(1,-1)), linetype = "dashed") +
        geom_hline(yintercept = 1.3, linetype = "dashed"), legend = "none")+
  ggtitle(" ") +
  geom_text_repel(size = 4,aes(FC_log, pvalue_log, label = sig_both_20)) +
  theme(axis.text = element_text(size = 16)) 

ggsave("Volcano Plot.png", height = 4.5, width = 5,path = volcano.path)

pdf(paste0(volcano.path,"/Volcano Plot.pdf"),height = 4.5, width = 5)
print(Volcano)
dev.off()

Volcano2 <- vacanno %>%
  mutate(darks = 
           ifelse(P_trans == 0.9 &
                    FC_log < 1 &
                    FC_log >-1 &
                    Pvalue > 0.05,
                    "grey", FC_col)) %>%
ggscatter(.,x = "FC_log", 
                  y = "pvalue_log",
                  main = "Vacanno",
                  size = "pvalue_log",
                  ylab = "-log10(p)",
                  xlab = "log2(FC)",
                  fill = "darks",
                  color = "black",
                  shape = 21,
                  font.x = c(16),
                  font.y = c(16),
                  font.label = c(16,"black"),
                  alpha = "P_trans",
                  repel = T) + 
      scale_fill_manual(values = c("#3B4992","#EE0000","grey")) +
      geom_vline(xintercept = c(a <- c(1,-1)), linetype = "dashed") +
      geom_hline(yintercept = 1.3, linetype = "dashed") +
      ggtitle(" ") +
      geom_text_repel(size = 4,aes(FC_log, pvalue_log, label = sig_both_20)) +
  guides(alpha = FALSE, fill = F, size = F) +
  theme(axis.text = element_text(size = 16)) 

ggsave("Volcano Plot2.png", height = 4.5, width = 5,path = volcano.path)

pdf(paste0(volcano.path,"/Volcano Plot2.pdf"),height = 4.5, width = 5)
print(Volcano2)
dev.off()

Volcano3 <- vacanno %>%
  mutate(darks = 
           ifelse(P_trans == 0.9 ,
                   "grey", FC_col)) %>%
  ggscatter(.,x = "FC_log", 
            y = "pvalue_log",
            main = "Vacanno",
            size = "pvalue_log",
            ylab = "-log10(p)",
            xlab = "log2(FC)",
            color = "darks",
            # color = "black",
            # shape = 21,
            font.x = c(16),
            font.y = c(16),
            font.label = c(16,"black"),
            # alpha = "P_trans",
            repel = T) + 
  scale_color_manual(values = c("#3B4992","grey","#EE0000")) +
  geom_vline(xintercept = c(a <- c(1,-1)), linetype = "dashed") +
  geom_hline(yintercept = 1.3, linetype = "dashed") +
  ggtitle(" ") +
  geom_text_repel(size = 4,aes(FC_log, pvalue_log, label = sig_both_20)) +
  guides(alpha = FALSE, color = F, size = F) +
  theme(axis.text = element_text(size = 16)) 

ggsave("Volcano Plot3.png", height = 4.5, width = 5,path = volcano.path)

pdf(paste0(volcano.path,"/Volcano Plot3.pdf"),height = 4.5, width = 5)
print(Volcano3)
dev.off()


# boxplot
dir.create("boxplot")
boxplot.path <- paste0(getwd(),"/boxplot")

Boxplot_App <- function(data, factor, i){
  ggpar(ggboxplot(data, 
                  x = colnames(data)[factor], 
                  y = paste0("`",colnames(data)[i],"`"), 
                  color = "black",
                  fill = colnames(data)[factor],
                  add = c("point"),
                  xlab = "Treatment",
                  width = 0.6,
                  ylab = "Relative Intensity",
                  add.params = list(size = 3, 
                                    jitter = 0.2)),
        legend.title = "Treatment",                         
        font.main = c(16, "plain", "black"),            
        font.x = c(16, "plain", "black"),                  
        font.y = c(16, "plain", "black"),                
        legend = " ",                                  
        font.legend = c(16, "plain", "black"),               
        rotate = F,                                         
        ticks = T,                                          
        tickslab = T, 
        palette = "aaas",
        xtickslab.rt = 0) +
    theme(axis.text = element_text(size = 16)) 
   }


Boxplot_App(data_gapfilled,2,3)



com_list <- NULL
for (i in 2 : length(levels(data_gapfilled$Treatment))) {
  com_list [[i-1]] = c(unique(data_gapfilled$Treatment)[1], 
                       unique(data_gapfilled$Treatment)[i])
}

my_comparisons <- com_list
position_num  <- c("1.1","1.2","1.3","1.4","1.9","2.1","2.3","2.5","2.7","2.9")
position_num = as.numeric(position_num)
compound_box = colnames(data_gapfilled)


for(i in 3:ncol(data_gapfilled)){
  AA <- Boxplot_App(data_gapfilled,2,i)
  sig_position <-NULL
  for (j in 1:length(levels(data_gapfilled$Treatment))){
    sig_position[j] = max(data_gapfilled[,i]) * position_num[j]
  }
  comound <- AA + stat_signif(comparisons = my_comparisons,
                              test = "t.test",
                              step_increase = 1,
                              y_position = sig_position
                            )
  ggsave(paste0(compound_box[i], ".png"),
         height = 4, width = 4,path = boxplot.path)
  
  pdf(paste0(boxplot.path,"/",compound_box[i], ".pdf"),height = 4, width = 4)
  print(comound)
  dev.off()
}

#  Zscore
dir.create("Zscore")
zscore.path <- paste0(getwd(),"/Zscore")
plsda_VIP %>% 
data_gapfilled

zscoreF <-function(x) {(x - mean(x))/sd(x)}
data_zscore <- apply(data_gapfilled[,-c(1,2)],2,zscoreF)
data_zscored <- cbind(data.frame(Treatment = data_gapfilled$Treatment),data_zscore)

c <-NULL
for(i in 2:ncol(data_zscored)){
  a <-NULL
  a <- matrix(c(data_zscored[,i], rep(colnames(data_zscored)[i], nrow(data_zscored))), nrow = nrow(data_zscored), ncol = 2)
  b <- cbind(data_zscored[,1],a)
  rownames(b) = rownames(data_zscored)
  c <-rbind(c, b)
}
colnames(c) = c("Treatment", "Z_score","Metabolites")
data_zscoreReady = c

for (i in 1:length(levels(data_gapfilled$Treatment))){
  print(paste0(i))
  data_zscoreReady[,1][data_zscoreReady[,1] == paste0(i)] = levels(data_gapfilled$Treatment)[i]
}
data_zscoreReady = data.frame(data_zscoreReady, stringsAsFactors = F)
data_zscoreReady$Treatment = factor(data_zscoreReady$Treatment,
                                    levels = c(data_zscoreReady$Treatment %>% 
                                                 unique() %>% .[1] %>% as.character(),
                                               data_zscoreReady$Treatment %>% 
                                                 unique() %>% .[2] %>% as.character())
                                    )
data_zscoreReady$Z_score = as.numeric(data_zscoreReady$Z_score)

VIP10 <- plsda_VIP %>% filter(VIP > 1) %>%
  arrange(desc(VIP)) %>%
  select(name) %>%
  as.data.frame() %>%
  as.vector() %>%
  pull() %>%
  .[1:10]

Zscore <- data_zscoreReady %>% filter(Metabolites %in% VIP10) %>%
  ggscatter(.,
            x = "Metabolites", 
            y = "Z_score",
            size = 5,
            xlab = "",
            ylab = "Z-score",
            ylim = c(-3, 3),
            combine = T,
            fill = "Treatment",
            shape = 21,
            color = "black", 
            palette = "aaas", 
            rotate = T,
            font.main = c(16, "plain", "black"), 
            font.x = c(16, "plain", "black"),                  
            font.y = c(16, "plain", "black"),                
            legend = "top", 
            font.legend = c(16, "plain", "black"), 
            repel = T) + 
  geom_vline(xintercept = c(a <- c(1:length(VIP10))), linetype = "dashed") +
  theme(axis.text = element_text(size = 16)) +
  guides(fill = guide_legend(title = ''))
ggsave("Z-score.png", height = 4, width = 6,path = zscore.path)
pdf(file = paste0(zscore.path,"/Zscore.pdf"),height =4, width =6)
print(Zscore)
dev.off()

# ------- ROC -----

dir.create("ROC")
roc.path <- paste0(getwd(),"/ROC")

ROC_group1 <- data_gapfilled %>% 
  filter(Treatment == data_gapfilled$Treatment %>%
           unique() %>%
           .[1] %>% 
           as.character()
         ) 

ROC_group2 <- data_gapfilled %>%
  filter(Treatment == data_gapfilled$Treatment %>%
           unique() %>%
           .[2] %>% 
           as.character()
         ) 

ROC_group1$Treatment = "1"
ROC_group2$Treatment = "0"
ROC_ready <- rbind(ROC_group1,ROC_group2)

MetaROC <- NULL
for (i in 3:ncol(ROC_ready)) {
  ROCvalue <- roc(ROC_ready[,2],
                  ROC_ready[,i],
                  plot=T,
                  levels=c("0", "1"),
                  direction = ">")
  
  ROC_value <- ROCvalue$auc[1]
  MetaROC[i-2] = round(ROCvalue$auc[1],2)
  cat("Show the metabolites AUC over 0.7")
  cat("\n")
  if(ROC_value > 0.5) {
    cat(paste0(colnames(ROC_ready)[i]," AUC: ",round(ROCvalue$auc[1],2)))
    cat("\n")
    ROCplot <- pROC::ggroc(ROCvalue,
                           legacy.axes = TRUE,
                           color = "#3B4992",
                           size = 0.5) +
      labs(x = "1-Specificity", y = "Sensitivity") +
      ggpubr::theme_pubr() +
      geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1),
                   color = "darkgrey", linetype = "dashed") +
      annotate(geom = "text",
               x = 0.7, y = 0.1,
               size = 5,
               color = "darkgrey",
               label = paste0("AUC = ",round(ROC_value,2))) +
      theme(axis.title.y = element_text(size = 16, angle = 90)) +
      theme(axis.title.x = element_text(size = 16, angle = 0)) +
      theme(axis.text = element_text(size = 16))
    print(ROCplot)
    ggsave(paste0(colnames(ROC_ready)[i],".png"),
           height = 4.5, width = 5, path =roc.path )
    pdf(paste0(roc.path,"/",colnames(ROC_ready)[i],".pdf"),height = 4.5, width = 5)
    print(ROCplot)
    dev.off()
  }
}

MetaROC <- as.data.frame(MetaROC)
rownames(MetaROC) = colnames(ROC_ready)[-(1:2)]
colnames(MetaROC)[1] = "AUC"

MetaROC$order <- c(1:nrow(MetaROC))
MetaROC$name <- rownames(MetaROC)
MetaROC$label <- ""
MetaROC$fill <- " "
for(i in 1:nrow(MetaROC)){
  if(MetaROC$AUC[i] > 0.9){
    MetaROC$label[i] = MetaROC$name[i]
    MetaROC$fill[i] = "red"
  }else{
    MetaROC$label[i] = NA
    MetaROC$fill[i] = "grey"
  }
}
MetaROC$fill <- as.factor(MetaROC$fill)

ROC_plot <- ggscatter(MetaROC,"order","AUC",
          shape = 21,fill = "fill",
          color = "black",size = "AUC",
          palette = "aaas",
          alpha = "fill",
          xlab = "Metabolites",
          ylab = "AUC score",
          font.main = c(16, "plain", "black"), 
          font.x = c(16, "plain", "black"),                  
          font.y = c(16, "plain", "black"),                
          # legend = "top", 
          # font.legend = c(16, "plain", "black"),
          legend = "") +
  geom_text_repel(
    data = MetaROC,
    aes(x=order,y=AUC,label = label),
    size = 3,
    segment.color = "black", show.legend = FALSE ) +
  geom_hline(yintercept = 0.9, linetype = "dashed") +
  theme(axis.text = element_text(size = 16))
ggsave("AUC plot.png",height = 4, width = 5, path =roc.path)

pdf(paste0(roc.path,"/AUC plot.pdf"),height = 4, width = 5)
ROC_plot
dev.off()

MetaROC %>% arrange(AUC) %>%
write.csv(.,paste0(roc.path,"/MetaROC.csv"))

# Heatmap

dir.create("Heatmap")
Heatmap.path <- paste0(getwd(),"/Heatmap")

VIP_10 <- plsda_VIP %>% arrange(desc(VIP)) %>%
  select(name) %>%
  head(10) %>%
  pull()
heatmap_data_process <- data_gapfilled[,match(VIP_10, colnames(data_gapfilled))]
heatmap_data_process_scale <- scale(heatmap_data_process)

color_table = c("#3e478d","#c8361f", "#3C894B", "#5a1d74", "#007f7e","#816800", "#003082","#8A0038", "#5B7E0F", "#9400D3")

myheatmap_color_table <- color_table[1:length(unique(data_gapfilled$Treatment))]
names(myheatmap_color_table) = unique(data_gapfilled$Treatment) %>% as.character()
heatmap_col = list(Treatment = myheatmap_color_table)

annot_df <- data_gapfilled %>% select(Sample, Treatment) %>%
  textshape :: column_to_rownames()

ha1 <- HeatmapAnnotation(df = annot_df, col = heatmap_col,
                     annotation_legend_param = list(Treatment = list(
                                                        title = "Treatment",
                                                        title_gp = gpar(fontsize = 12),
                                                        title_position = "topcenter",
                                                        labels_gp = gpar(fontsize = 12)
                                                        )
                                                    )
                       )


mycol <- colorRamp2(c(min(heatmap_data_process_scale), 
                      median(heatmap_data_process_scale), 
                      max(heatmap_data_process_scale)), 
                    c("blue", "white", "red"))

heatmapresult <- Heatmap(t(heatmap_data_process_scale),
                         heatmap_legend_param = list(title = "Metabolites", 
                                                     title_gp = gpar(fontsize = 12), 
                                                     labels_gp = gpar(fontsize = 12),
                                                     title_position = "topleft",
                                                     legend_direction = "vertical"),
                         col = mycol,
                         top_annotation = ha1)

Heatmap_plot <- ggplotify::as.ggplot(heatmapresult)

pdf(paste0(Heatmap.path,"/heatmap.pdf"),height = 4, width = 8)
Heatmap_plot
dev.off()


