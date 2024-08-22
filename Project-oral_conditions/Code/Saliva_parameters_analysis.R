## ---------------------------
##
## Script name: Saliva parameters analysis
##
## Purpose of script: Statistical analysis to the salivary basic parameters
##
## Author: Eduardo Esteves
##
## Date Created: 2023-03-07
##
## Copyright (c) Eduardo Esteves, 2023
## Email: eduardojmesteves@gmail.com
##
## ---------------------------
##
## Notes:
##   
##
## ---------------------------

## set working directory for Mac and PC

setwd("/Users/ee/git/r")    # Edd's working directory (mac)

## ---------------------------

## load up the packages we will need:  (uncomment as required)

library(tidyverse)
library(hrbrthemes)
library(viridis)
library(scales)
library(BiocManager)
library(rstatix)
# source("functions/packages.R")       # loads up all the packages we need

## ---------------------------

## load up our functions into memory

# source("functions/summarise_data.R") 

## import the dataframe---------------------------

df = read.table(file = '/Users/ee/git/r/Project_periodontitis/Data/data.csv', 
                header = TRUE, sep = ',', quote = "", row.names = NULL, 
                stringsAsFactors = FALSE, check.names = FALSE)

data = df[,c('Patologia_Oral','Protein_concentration','Volume',"Flux",'pH')]

# Rename the oral conditions to English
data[data == "Saudavel"] <- "Healthy"
data[data == "Gengivite"] <- "Gingivitis"
data[data == "Mucosite"] <- "Mucositis"
data[data == "Peri-implantite"] <- "Peri-implantitis"
data[data == "Gengivite e Mucosite"] <- "Gingivitis and Mucositis"
data[data == "Gengivite e Peri-implantite"] <- "Gingivitis and Peri-implantitis"
data[data == "Mucosite e Periodontite"] <- "Mucositis and Periodontitis"
data[data == "Mucosite e Peri-implantite"] <- "Mucositis and Peri-implantitis"
data[data == "Periodontite e Peri-implantite"] <- "Periodontitis and Peri-implantitis"
data[data == "Gengivite e Mucosite Peri-implantite"] <- "Gingivitis, Mucositis and Peri-implantitis"

data$Patologia_Oral <- factor(data$Patologia_Oral,
                              levels=c("Healthy",
                                       "Gingivitis",
                                       "Mucositis",
                                       "Peri-implantitis",
                                       "Gingivitis and Mucositis",
                                       "Gingivitis and Peri-implantitis",
                                       "Mucositis and Periodontitis",
                                       "Mucositis and Peri-implantitis",
                                       "Periodontitis and Peri-implantitis",
                                       "Gingivitis, Mucositis and Peri-implantitis"))


## Oral Pathologies Individual Visualization---------------------------

data = data %>%
  mutate(group_conditions_short = recode(Patologia_Oral,"Healthy" = "H",
                                         "Gingivitis" = "G",
                                         "Mucositis" = "M",
                                         "Peri-implantitis" = "PI",
                                         "Gingivitis and Mucositis" = "G+M",
                                         "Gingivitis and Peri-implantitis" = "G+PI",
                                         "Mucositis and Periodontitis" = "M+P",
                                         "Mucositis and Peri-implantitis" = "M + PI",
                                         "Periodontitis and Peri-implantitis" = "P+PI",
                                         "Gingivitis, Mucositis and Peri-implantitis" = "G+M+PI"))

# BoxPlot - Total protein concentration
data %>%
  ggplot( aes(x=group_conditions_short, y=Protein_concentration, fill=group_conditions_short)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11),
    axis.text.x = element_text(angle = 30, hjust = 0.5, vjust = 0.5)
  ) +
  ggtitle("Protein concentration distribution visualisation using BoxPlot chart") +
  xlab("")

# kt <- kruskal.test(Protein_concentration ~ Patologia_Oral, data = data)
# if(kt$p.value > 0.05){
#   pw <-pairwise.wilcox.test(df$Protein_concentration, g=data$group_conditions,
#                             p.adjust.method = "BH")
# }

# BoxPlot - Saliva volume
data %>%
  ggplot( aes(x=group_conditions_short, y=Volume, fill=group_conditions_short)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11),
    axis.text.x = element_text(angle = 30, hjust = 0.5, vjust = 0.5)
  ) +
  ggtitle("Saliva Volume distribution visualisation using BoxPlot chart") +
  xlab("")


# BoxPlot - Saliva flux
data %>%
  ggplot( aes(x=group_conditions_short, y=Flux, fill=group_conditions_short)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11),
    axis.text.x = element_text(angle = 30, hjust = 0.5, vjust = 0.5)
  ) +
  ggtitle("Saliva Flux distribution visualisation using BoxPlot chart") +
  xlab("")


# BoxPlot - Saliva pH
data %>%
  ggplot( aes(x=group_conditions_short, y=pH, fill=group_conditions_short)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11),
    axis.text.x = element_text(angle = 30, hjust = 0.5, vjust = 0.5)
  ) +
  ggtitle("Saliva pH distribution visualisation using BoxPlot chart") +
  xlab("")

## Oral Pathologies Individual Statistical analysis---------------------------
# Summary statistics - Total Protein concentration
summary_statistics <- data %>% 
  group_by(Patologia_Oral) %>%
  get_summary_stats(Protein_concentration,Volume,Flux,pH, type = "common")

# Computation analysis
kruskal_protein_concentration <- data %>% kruskal_test(Protein_concentration ~ group_conditions_short)
kruskal_volume <- data %>% kruskal_test(Volume ~ group_conditions_short)
kruskal_flux <- data %>% kruskal_test(Flux ~ group_conditions_short)
kruskal_ph <- data %>% kruskal_test(pH ~ group_conditions_short)

kruskal_total = rbind(kruskal_protein_concentration,
                      kruskal_volume,
                      kruskal_flux,
                      kruskal_ph)
#write.csv(kruskal_total, file = "kruskal_total.csv")

pwc_ind_pc <- data %>% 
  dunn_test(Protein_concentration ~ Patologia_Oral, p.adjust.method = "bonferroni") 
pwc_ind_v <- data %>% 
  dunn_test(Volume ~ Patologia_Oral, p.adjust.method = "bonferroni")
pwc_ind_f <- data %>% 
  dunn_test(Flux ~ Patologia_Oral, p.adjust.method = "bonferroni") 
pwc_ind_ph <- data %>% 
  dunn_test(pH ~ Patologia_Oral, p.adjust.method = "bonferroni") 

pwc_ind = rbind(pwc_ind_pc,
                pwc_ind_v,
                pwc_ind_f,
                pwc_ind_ph)
                    
#write.csv(pwc_ind, file = "pwc_ind.csv")


## Oral Pathologies Grouped Visualization--------------------------- 
# Combine the oral conditions
data = data %>%
  mutate(group_conditions = recode(Patologia_Oral, 'Gingivitis' = 'Soft tissue',
                                   'Mucositis' = 'Soft tissue',
                                   'Gingivitis and Mucositis' = 'Soft tissue',
                                   'Peri-implantitis' = 'Hard tissue',
                                   'Periodontitis and Peri-implantitis' = 'Hard tissue',
                                   'Gingivitis and Peri-implantitis' = 'Mix',
                                   'Mucositis and Peri-implantitis' = 'Mix',
                                   'Gingivitis, Mucositis and Peri-implantitis' = 'Mix',
                                   'Mucositis and Periodontitis' = 'Mix'))

data$group_conditions <- factor(data$group_conditions , levels=c("Healthy", "Soft tissue", "Hard tissue", "Mix"))


# Violin basic - Protein concentration
data %>%
  ggplot( aes(x=group_conditions, y=Protein_concentration, fill=group_conditions)) +
  geom_violin() +
  geom_boxplot(width=0.1, color="grey", alpha=0.2) +
  scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") +
  theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("Protein concentration distribution visualisation using Violin chart") +
  xlab("")

# Violin basic - Saliva volume
data %>%
  ggplot( aes(x=group_conditions, y=Volume, fill=group_conditions)) +
  geom_violin() +
  geom_boxplot(width=0.1, color="grey", alpha=0.2) +
  scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") +
  theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("Saliva volume distribution visualisation using Violin chart") +
  xlab("")

# Violin basic - Flux
data %>%
  ggplot( aes(x=group_conditions, y=Flux, fill=group_conditions)) +
  geom_violin() +
  geom_boxplot(width=0.1, color="grey", alpha=0.2) +
  scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") +
  theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("Saliva flux distribution visualisation using Violin chart") +
  xlab("")


# Violin basic - pH
data %>%
  ggplot( aes(x=group_conditions, y=pH, fill=group_conditions)) +
  geom_violin() +
  geom_boxplot(width=0.1, color="grey", alpha=0.2) +
  scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") +
  theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("pH distribution visualisation using Violin chart") +
  xlab("")

## Oral Pathologies Grouped Statistical analysis---------------------------
# Summary statistics - Total Protein concentration
summary_statistics_grouped <- data %>% 
  group_by(group_conditions) %>%
  get_summary_stats(Protein_concentration,Volume,Flux,pH, type = "common")
#write.csv(summary_statistics_grouped, file = "summary_statistics_grouped.csv")


# Computation analysis
kruskal_protein_concentration_grouped <- data %>% kruskal_test(Protein_concentration ~ group_conditions)
kruskal_volume_grouped <- data %>% kruskal_test(Volume ~ group_conditions)
kruskal_flux_grouped <- data %>% kruskal_test(Flux ~ group_conditions)
kruskal_ph_grouped <- data %>% kruskal_test(pH ~ group_conditions)

kruskal_total_grouped = rbind(kruskal_protein_concentration_grouped,
                              kruskal_volume_grouped,
                              kruskal_flux_grouped,
                              kruskal_ph_grouped)
#write.csv(kruskal_total_grouped, file = "kruskal_total_grouped.csv")

pwc_grouped_pc <- data %>% 
  dunn_test(Protein_concentration ~ group_conditions, p.adjust.method = "bonferroni") 
pwc_grouped_v <- data %>% 
  dunn_test(Volume ~ group_conditions, p.adjust.method = "bonferroni")
pwc_grouped_f <- data %>% 
  dunn_test(Flux ~ group_conditions, p.adjust.method = "bonferroni") 
pwc_grouped_ph <- data %>% 
  dunn_test(pH ~ group_conditions, p.adjust.method = "bonferroni") 

pwc_grouped = rbind(pwc_grouped_pc,
                    pwc_grouped_v,
                    pwc_grouped_f,
                    pwc_grouped_ph)
#write.csv(pwc_grouped, file = "pwc_grouped.csv")
