
# Draw heatmaps
pacman::p_load(pheatmap, ggpubr, ggplot2,ggplotify)

test <- read.csv("test_a.csv", header = T, stringsAsFactors = F,row.names = "Name")
pheatmap(test, scale = "row", color = colorRampPalette(c("blue", "white", "red"))(50),clustering_distance_rows = "correlation")


annotation_col <- read.csv("annotation_col_a.csv", header = T,stringsAsFactors = F, row.names = "Name")


annotation.row <- read.csv("annotation_row_a.csv", header = T,stringsAsFactors = F, row.names = "Name")
labels_row <- annotation.row$Label
annotation.row$Label <- NULL
annotation_row <- annotation.row


# factor(annotation_col$CellType)
# factor(annotation_col$Time)
# factor(annotation_row$GeneClass)
# # Specify colors

ann_colors = list(
  # Time = c('12h' = "#7570B3", '24h' = "#66A61E"),
  Treatment = c('Con' = "#1B9E77", 'HIBA' = "#D95F02" ),
  Class = c('Bacteroidetes' = '#7570B3', 'Firmicutes' = "#E7298A",
            'Verrucomicrobia' = "#66A61E",
            'Actinobacteria' = "#1B9E77",
            'Deferribacteres' = "#D95F02",
            'Proteobacteria' = "grey")
)

test_b<- pheatmap(test, 
         scale = "row",
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
         # color = colorRampPalette(c("blue", "white", "red"))(50),
         annotation_col = annotation_col, 
         annotation_row = annotation_row,
         annotation_colors = ann_colors,
         labels_row = labels_row,
         border=FALSE,
         border_color = "black",
         treeheight_row = 30, 
         treeheight_col = 30,
         cellwidth = 15, 
         cellheight = 12, 
         fontsize = 8,
         clustering_method = "ward.D",
         # "ward.D","single", "complete", "average", "mcquitty", "median", "centroid", "ward.D2"
         # cluster_rows = FALSE, 
         cluster_cols = FALSE,
         gaps_col = c(5),
         cutree_col = 2
         # display_numbers = TRUE,
         # display_numbers = matrix(ifelse(test > 5, "*", ""), nrow(test))
         )


test_c = ggplotify::as.ggplot(test_b)

ggsave("test_c.pdf", width = 10,height = 10)



