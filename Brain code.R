library(readr)
library(sf)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
# library(magrittr)
library(scales)
library(stats)
library(ggpubr)

file_pre <- "./" # same directory as shiny file

# loading data for geomap
brain_sf <- st_read(dsn = paste(file_pre, "brain.geojson", sep = "/"), 
                    stringsAsFactors = FALSE, quiet = TRUE)

# loading Metabolites data
demo_df <- read_csv(file = paste(file_pre, "longdata3.csv", sep = "/"), 
                    col_names = TRUE) 

fun1 <- function(df){
  mypal <- colorRampPalette(c("#91D1C1","#FFEB84","#E74B35"))(101) 
  df1 <- df %>% select(Region, Treatment, Intensity) %>% 
    group_by(Region, Treatment) %>% 
    summarise(meanInten = mean(Intensity)) %>% ungroup() 
  minInten <- min(df1$meanInten)
  maxInten <- max(df1$meanInten)
  df_A <- df1 %>% mutate(normInten = round((meanInten - minInten)/(maxInten-minInten)*100)) %>% 
    mutate(color = mypal[normInten + 1]) %>% 
    select(Region, Treatment, color) %>% 
    replace_na(list(color = "#91D1C1"))
  
  
  df_B <- df %>% select(Region, Treatment, Intensity) %>% 
    right_join(Region_df, by = "Region") %>% 
    mutate(A2 = factor(x = Treatment, levels = c("Control", "C4", "Est"))) %>% 
    select(x, Intensity, Treatment = A2) 
  
  y_p <- df %>% select(Region, Treatment, Intensity) %>% 
    group_by(Region, Treatment) %>% 
    summarise(y_mean = mean(Intensity),
              y_diff = sd(Intensity, na.rm = TRUE)/sqrt(sum(!is.na(Intensity)))) %>% 
    ungroup() %>% 
    transmute(y_up = y_mean + y_diff) %>% pull(y_up) %>% 
    max() %>% `*`(1.1) 
  
  
  return(list(df_A, df_B, y_p))
  
}


Region_df <- data.frame(x = 1:3, stringsAsFactors = FALSE,
                        Region = c("Hippocampus",
                                   "Hypothalamus", 
                                   "Basal ganglia"))

meta <- demo_df$Metabolite %>% unique()
length(meta)


for (i in 1:364) {
  


reaction1 <- demo_df %>%
    filter(Region %in% Region_df$Region, 
           Metabolite == meta[i]) %>% 
    fun1()
  

sf_w3 <- brain_sf %>% left_join(filter(reaction1[[1]], Treatment == "Control"), by = "Region") %>% 
    select(-Treatment)

sf_w3$color[c(2,3,6:10)] = rep("white",7)


sf_w16 <- brain_sf %>% left_join(filter(reaction1[[1]], Treatment == "C4"), by = "Region") %>% 
  select(-Treatment)

sf_w16$color[c(2,3,6:10)] = rep("white",7)

sf_w59 <- brain_sf %>% left_join(filter(reaction1[[1]], Treatment == "Est"), by = "Region") %>% 
  select(-Treatment)
sf_w59$color[c(2,3,6:10)] = rep("white",7)

library(patchwork)
w3 <- ggplot(sf_w3) + 
    geom_sf(aes(fill = Region), color = "black", show.legend = FALSE) + 
    coord_sf() + 
    scale_fill_manual(values = sf_w3$color, 
                      labels = sf_w3$Region) + 
    theme_void()

w16 <- ggplot(sf_w16) + 
  geom_sf(aes(fill = Region), color = "black", show.legend = FALSE) + 
  coord_sf() + 
  scale_fill_manual(values = sf_w16$color, 
                    labels = sf_w16$Region) + 
  theme_void()

w59 <- ggplot(sf_w59) + 
  geom_sf(aes(fill = Region), color = "black", show.legend = FALSE) + 
  coord_sf() + 
  scale_fill_manual(values = sf_w59$color, 
                    labels = sf_w59$Region) + 
  theme_void()

w3/w16/w59

ggsave(paste0(meta[i],".pdf"))
}




