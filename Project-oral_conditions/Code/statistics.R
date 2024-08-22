# Load the package and prepare the input data frame
#install.packages('pathfindR')
#install.packages("BiocManager")
#install.packages('tidyverse')
#install.packages('extrafont')
#install.packages("ggplot2")

library(pathfindR)
library(BiocManager)
library(tidyverse)
library(extrafont)
library(ggplot2)


        
df = read.table(file = '/Users/ee/git/r/Project_periodontitis/data.csv', 
                    header = TRUE, sep = ',', quote = "", row.names = NULL, 
                    stringsAsFactors = FALSE, check.names = FALSE)
colnames(df) <- make.unique(names(df))

col_names = names(df)
sapply(df, class)

colnames(df)[1] <- "id"
df[df == "Pacientes"] <- "id"
df[df == "Gengivite"] <- "Gingivitis"
df[df == "Mucosite"] <- "Mucositis"
df[df == "Gengivite e Peri-implantite"] <- "Gingivitis and Peri-implantitis"
df[df == "Gengivite e Mucosite"] <- "Gingivitis and Mucositis"
df[df == "Peri-implantite"] <- "Peri-implantitis"
df[df == "Mucosite e Peri-implantite"] <- "Mucositis and Peri-implantitis"
df[df == "Saudavel"] <- "Healthy"
df[df == "Periodontite e Peri-implantite"] <- "Periodontitis and Peri-implantitis"
df[df == "Gengivite e Mucosite Peri-implantite"] <- "Gingivitis; Peri-implantitis and Mucositis"
df[df == "Mucosite e Periodontite"] <- "Mucositis and Periodontitis"


##new df with id and patologia oral columns:
#df1 <- df[,c('id','Patologia_Oral')]

df = df %>%
  mutate(group_conditions = recode(Patologia_Oral, 'Gingivitis' = 'Soft tissue',
                                 'Mucositis' = 'Soft tissue',
                                 'Gingivitis and Mucositis' = 'Soft tissue',
                                 'Peri-implantitis' = 'Hard tissue',
                                 'Periodontitis and Peri-implantitis' = 'Hard tissue',
                                 'Gingivitis and Peri-implantitis' = 'mix',
                                 'Mucositis and Peri-implantitis' = 'mix',
                                 'Gingivitis; Peri-implantitis and Mucositis' = 'mix',
                                 'Mucositis and Periodontitis' = 'mix'))


disease_stat <- df %>%
  count(group_conditions)

healthy_n <- disease_stat %>%
  filter(group_conditions == "Healthy")

soft_tissue_n <- disease_stat %>%
  filter(group_conditions == "Soft tissue")

hardtissue_n <- disease_stat %>%
  filter(group_conditions == "Hard tissue")

mix_n <- disease_stat %>%
  filter(group_conditions == "mix")


healthy_color <- "#BEBEBE"
soft_color <- "#0000FF"
hard_color <- "#FF0000"
mix_color <- "#AE4371"


ggplot(data = df)+
  geom_point(mapping = aes(x=Protein_concentration, y=pH)) +
  facet_wrap(~ group_conditions, nrow = 2) +
  labs(
    title = "Total Protein Concentration VS pH",
    x = "Total Protein Concentration",
    y = "pH"
  )

kt <- kruskal.test(Protein_concentration ~ group_conditions, data = df)
if(kt$p.value < 0.05){
  pw <-pairwise.wilcox.test(df$Protein_concentration, g=df$group_conditions,
                       p.adjust.method = "BH")
}

strip_chart <- df %>%
  ggplot(aes(x=group_conditions, y=Protein_concentration, fill=factor(group_conditions), alpha = 0.7)) +
  stat_summary(fun = median, geom = "crossbar", show.legend = FALSE) +
  geom_jitter(show.legend = FALSE, color = "black", fill='pink') +
  labs(
    x = NULL,
    y = "Total Protein Concentration") +
  theme_classic() +
  scale_fill_manual(values = c("red","blue","green","purple"))
  # scale_x_discrete(breaks = c("Healthy","hard tissue","soft tissue","mix"),
  #                  labels = c(glue("Healthy<br>(N={healthy_n})"),
  #                             glue("hard tissue<br>(N={hardtissue_n})"),
  #                             glue("soft tissue<br>(N={soft_tissue_n})"),
  #                             glue("mix<br>(N={mix_n})")))
strip_chart


