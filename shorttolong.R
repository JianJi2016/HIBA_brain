data <- read_csv("brainshort.csv")


install.packages("reshape2")
library(reshape2)

data2 <- melt(data,1:3) %>% 
  select(2:5)

colnames(data2) = c("Region","Treatment","Metabolite","Intensity")

write.csv(data2,"longdata3.csv")
