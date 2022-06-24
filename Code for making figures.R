# Figure 1b 1c and 1d
# the data_B.csv is the raw data for Fig 1c
pacman::p_load(tidyverse,ggsci)
data_1 <- read.csv("data_B.csv",header = T)
data <- data_1
data$group = factor(data$group, levels = c("Movement speed",
                                           "Movement time",
                                           "Central area residence time",
                                           "12 h sleep duration",
                                           "incubation period",
                                           "sleep time"))
empty_bar=1
to_add = data.frame( matrix(NA, empty_bar*nlevels(data$group), ncol(data)) )
colnames(to_add) = colnames(data)
to_add$group=rep(levels(data$group), each=empty_bar)
data=rbind(data, to_add)
data=data %>% arrange(group)
data$id=seq(1, nrow(data))

# Get the name and the y position of each label
label_data=data
number_of_bar=nrow(label_data)
angle= 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_data$hjust<-ifelse( angle < -90, 1, 0)
label_data$angle<-ifelse(angle < -90, angle+180, angle)

# prepare a data frame for base lines
base_data=data %>% 
  group_by(group) %>% 
  summarize(start=min(id), end=max(id) - empty_bar) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))

# prepare a data frame for grid (scales)
grid_data = base_data
grid_data$end = grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start = grid_data$start - 1
grid_data=grid_data[-1,]

p = ggplot(data, aes(x=as.factor(id), y=value, fill=group)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  
  geom_bar(aes(x=as.factor(id), y=value, fill=group), 
           # width = "barwith",
           stat="identity"
  ) +
  geom_errorbar(aes(ymin = value- barwith, ymax = value + barwith), 
                width = 0.2, 
                color = "black",
                position = position_dodge(0.9))+
  
  
  # Add a val=100/75/50/25 lines. I do it at the beginning to make sur barplots are OVER it.
  geom_segment(data=grid_data, aes(x = end, y = 60, xend = start, yend = 60),
               colour = "red",  size=1 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 40, xend = start, yend = 40), 
               colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 20, xend = start, yend = 20), 
               colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  
  # Add text showing the value of each 100/75/50/25 lines
  annotate("text", x = rep(max(data$id),4), y = c(20, 40, 60, 80), label = c("20", "40", "60", "80") , 
           color="grey", size=3 , angle=0, fontface="bold", hjust=1) +
  
  
  ylim(-30,70) +
  theme_minimal() +
  theme(
    legend.position = "right",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm") 
  ) +
  coord_polar() + 
  guides(fill=F)+
  geom_text(data=label_data, aes(x=id, y=value+10, label=individual, hjust=hjust), 
            color="black", alpha=0.8, size=4, angle= label_data$angle, inherit.aes = FALSE ) +
  scale_fill_npg()+
  # Add base line information
  geom_segment(data=base_data, aes(x = start, y = -5, xend = end, yend = -5), 
               colour = "black", alpha=0.8, size=0.6 , inherit.aes = FALSE )  

print(p)
ggsave("B.pdf",
       width = 160,
       height = 160,
       units = "mm")


# Fig.2 Levels of neurotransmitters in brain and serum
neurotransmitters <- read.csv("neurotransmitters.csv")
library(ggpubr)
ggerrorplot(neurotransmitters,
            x = "Groups",
            y = "X5.HT",
            # add = "jitter",
            color = "Groups",
            palette = "npg") +
  coord_flip() 
ggsave("5-HT.pdf",height = 2,width = 5)




# Fig.3b
metabolites <- read.csv("metabolites number.csv")

library(ggpubr)
ggbarplot(metabolites,
          x = "type",
          y = "number",
          fill = "type",
          palette = "npg",
          label = T) +
  coord_flip() +
  theme(legend.position = "null")

ggsave("metabolites number.pdf")


# Fig.3c
class <- read.csv("ClassyFire categories.csv")
library(dplyr)

class %>% 
  group_by(class) %>%
  summarise(n = n()) %>%
  arrange(desc(n)) -> data

colnames(data) = c("category","count")

data$fraction <- data$count / sum(data$count)
data$ymax <- cumsum(data$fraction)
data$ymin <- c(0, head(data$ymax, n=-1))
data$labelPosition <- (data$ymax + data$ymin) / 2
data$label <- paste0(data$category, "\n value: ", data$count)
ggplot(data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=category)) +
  geom_rect() +
  coord_polar(theta="y") +
  xlim(c(2, 4)) +
  theme_void() +
  theme(legend.position = "none")

ggsave("ClassyFire categories.pdf")


# Fig. 5I the LDA analysis
df3 <- read.csv("LDA.csv")
library(dplyr)
df3$mpg_grp <- factor(ifelse(df3$LDAscore<0, "low", "high"), levels = c("low", "high"))

ggbarplot(df3 , x="Name", y= "LDAscore", 
          fill = "mpg_grp", color = "white", 
          palette = "npg", sort.val = "desc", 
          sort.by.groups = FALSE, 
          xlab = "", ylab = "LDA score (log10)", 
          legend.title="", rotate=TRUE)
ggsave("LDA barplot.pdf",width = 5, height = 4)

# 

pacman::p_load(ggpubr,ggplot2,magrittr)











# metabolic pathway based on metabolites in BG
pathway <- read.csv("msea_ora_result.csv")

library(ggpubr)
library(dplyr)

ggscatter(pathway[1:20,] %>%
            arrange(logP),
          x = "Name",
          y = "logP",
          shape = 21,
          fill = "logP",
          size = "Ratio") +
  coord_flip() +
  labs(y = "-logP",
       x = "")

ggsave("metabolic pathway.pdf", height = 6,width = 10)



# The correlation between the gut microbiota and the metabolome
# data for cytoscape to make Fig 6b-h

library(readr)
corr <- read.csv("correlation data.csv")
rownames(corr) = corr[,1]
corr = corr[,-1]
library(corrplot)
result_pair = cor(corr)
rows = rownames(result_pair)  # DD????
cols = colnames(result_pair)  # DD????
write.csv(result_pair,"result_pair.csv")
datap <- data.frame()
for (j in 19:34) {
    for (i in 1:18) {
    datap[i + 35*(j-1),1] = rows[i]
    datap[i + 35*(j-1),2] = result_pair[i,j]
    datap[i + 35*(j-1),3] = cols[j]
  }
}

library(dplyr)
colnames(datap) = c("Target", "value", "Source")
datap %>% 
  filter(
    value > 0.7| value < -0.7
  ) %>%
  filter(value != 1) -> datap4
write.csv(datap4,"edgeCytoscape.csv")






