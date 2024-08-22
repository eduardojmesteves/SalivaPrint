## ---------------------------
##
## Script name: Protein capillary electrophoresis (CE) analysis
##
## Purpose of script: Statistical analysis of protein CE results
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
## set working directory for Mac and PC---------------------------

setwd("/Users/eduardoesteves/dev/R/Project_periodontitis")

## load up the packages we will need---------------------------

library(BiocManager)
library(biosigner)
library(caret) #Feature plot of Predictor Variables https://www.rdocumentation.org/packages/caret/versions/6.0-93/topics/featurePlot
require(data.table)
library(dplyr) #Data manipulator
library(DataExplorer) #Exploratory Data Analysis (EDA) https://cran.r-project.org/web/packages/DataExplorer/readme/README.html 
library(extrafont)
library(ggplot2)
library(ggthemes)
library(hrbrthemes)
library(mixOmics)
#library(pathfindR)
library(reshape2)
library(readxl)
library(vtable)
library(showtext)
library(skimr) # Summary statistics https://cran.r-project.org/web/packages/skimr/vignettes/skimr.html
require(tidyverse)
library(tidyr)
library(viridisLite)
library(viridis)


# source("functions/packages.R")       # loads up all the packages we need

## ---------------------------

## load up our functions into memory

# source("functions/summarise_data.R") 

## ---------------------------
## import the data frame and prepare the data---------------------------

## Import the initial data frame
profiles <- read.table(file = '/Users/eduardoesteves/dev/R/Project_periodontitis/Data/profiles.csv', 
                       header = TRUE, sep = ',', quote = "", row.names = NULL, 
                       stringsAsFactors = FALSE, check.names = FALSE)

# Rename the column in profiles
colnames(profiles)[1] <- "id"

clinical <- read.table(file = '/Users/eduardoesteves/dev/R/Project_periodontitis/Data/clinical.csv', 
                       header = TRUE, sep = ',', quote = "", row.names = NULL, 
                       stringsAsFactors = FALSE, check.names = FALSE)

# Remove the duplicate columns
clinical <- clinical[, !duplicated(names(clinical))]
profiles <- profiles[, !duplicated(names(profiles))]

# Clinical plot
plot_str(clinical)
names(clinical)

#Remove the clutter information (columns with lower interest to the analysis)
summary_clinical <- dplyr::select(clinical, -c('Smoking',
                                               'Alcool',
                                               'Medication',
                                               'Medication_1',
                                               'Systemic_diseases',
                                               'psr_1q',
                                               'psr_2q',
                                               'psr_3q',
                                               'psr_4q',
                                               'psr_5q',
                                               'psr_6q'))

#Rename column names
names(summary_clinical)
setnames(summary_clinical, 
         old = c("Protein_concentration",
                 "n_tooth",
                 "n_implant",
                 "bop"),
         new = c("Protein concentration",
                 "Number of teeth",
                 "Number of implants",
                 "Bleeding on probing"))
names(summary_clinical)

#Rename oral pathology column values
unique(summary_clinical[,"oral_pathology_grouped", drop=FALSE])
table(summary_clinical$oral_pathology_grouped, useNA = 'always')


# Reorder the groups based on name
unique(profiles_melt$group)
summary_clinical$oral_pathology_grouped <- 
  factor(summary_clinical$oral_pathology_grouped,
         levels = c("Healthy",
                    "Soft_Tissues",
                    "Hard_Tissues"))

unique(summary_clinical[,"oral_pathology_grouped", drop=FALSE])

# Summary statistics table for the clinical data with grouped oral pathologies 
summary_clinical_table <- sumtable(summary_clinical,
                                   group = 'oral_pathology_grouped',
                                   align = 'center',
                                   group.long = FALSE,
                                   out = 'browser')


# library(data.table)
# fwrite(clinical, "clinical_data.csv")


## Exploratory Data Analysis---------------------------
### skim() package ###
#https://cran.r-project.org/web/packages/skimr/vignettes/skimr.html

skim(summary_clinical)
clinical_summary <- skim(summary_clinical)
write_csv(summary_clinical, file = "clinical_summary.csv")
# library(data.table)
# fwrite(df, "df.csv")
summary_clinical %>%
  skim() %>%
  partition()

skim(summary_clinical) %>%
  dplyr::filter(skim_variable == "oral_pathology_grouped")

skim(summary_clinical) %>%
  dplyr::select(skim_type, skim_variable, n_missing)
  
skim(clinical_sumtable) %>%
  dplyr::select(skim_type, skim_variable, numeric.mean)

# sum_table_by_pathology <- clinical_sumtable %>%
#   dplyr::group_by(oral_pathology) %>%
#   dplyr::select_if(is.numeric) %>%
#   skim()

write_csv(sum_table_by_pathology, file = "sum_table_by_pathology.csv")
### DataExplorer package ###
#https://boxuancui.github.io/DataExplorer/articles/dataexplorer-intro.html
introduce(summary_clinical)

plot_intro(summary_clinical)
plot_missing(summary_clinical)

#pdf("plot.pdf", width = 8, height = 5) # set the size of the plot

numeric_cols <- sapply(summary_clinical, is.numeric)
numeric_data <- summary_clinical[, numeric_cols]

plot_bar(numeric_data, by = "oral_pathology_grouped")

plot_bar(summary_clinical, by = "oral_pathology_grouped")
#dev.off() # close the PDF device

anyNA(summary_clinical$pH)
is.numeric(summary_clinical$pH)

plot_histogram(summary_clinical)

featurePlot(x = clinical$Protein_concentration, 
            y = clinical$oral_pathology_grouped, 
            plot = "scatter")

ggplot(clinical, aes(x = Protein_concentration, y = oral_pathology_grouped)) +
  geom_point() +
  labs(title = "Scatter Plot of Protein Concentration by Oral Pathology Group",
       x = "Protein Concentration",
       y = "Oral Pathology Group") +
  theme_minimal()

plot(summary_clinical)

## Profile data frame preparation for Visualisation and Analysis---------------------------

# Replace the group column in Profiles df with oral_pathology_grouped data from Clinical df
profiles_grouped <- profiles %>%
  left_join(summary_clinical %>% dplyr::select(id, 
                                       oral_pathology_grouped), by = "id") %>% # Merge the oral_pathology_grouped column into the profiles dataframe
  mutate(group = oral_pathology_grouped) %>% # Replace the group column in profiles with the merged oral_pathology_grouped column
  dplyr::select(-oral_pathology_grouped) %>% 
  na.omit() # Remove rows where group is NA

# Merge the data frame
clinical_profiles <- merge(clinical, profiles_grouped, by = "id", all = TRUE)

# Remove column "sample"
profiles_hold <- dplyr::select(profiles_grouped, -c("id"))

# Combine multiple columns of df into a single column
profiles_melt <- melt(profiles_hold)

# Reorder the groups based on name
unique(profiles_melt$group)
profiles_melt$group <- factor(profiles_melt$group, levels = c("Healthy",
                                                              "Soft_Tissues",
                                                              "Hard_Tissues"))
                              
unique(profiles_melt$group)
# Convert the 'value' column to a numeric data type
profiles_melt$variable <- as.numeric(as.character(profiles_melt$variable))

# Round the 'value' column to integers
profiles_melt$rounded_variable <- as.integer(round(profiles_melt$variable))

#create a summary table with mean and sd to be used on the geom_ribbon chart
summary_table <- profiles_melt %>% 
  group_by(group, rounded_variable) %>% 
  summarize(value.mean = mean(value),
            value.sd = sd(value))


## Profiles visualization---------------------------

p_plot_1 <- ggplot(summary_table, aes(x = rounded_variable,
                               y = rounded_variable,
                               fill = group)) +
  geom_ribbon(aes(ymin = value.mean - value.sd,ymax = value.mean + value.sd), alpha = 0.3) +
  geom_line(aes(y = value.mean, color = group), linetype = "solid") +
  scale_fill_manual(values = c("darkgreen", "orange","brown")) +
  scale_color_manual(values = c("darkgreen", "orange","brown")) +
  coord_cartesian(xlim = c(10, 140)) +
  scale_x_continuous(breaks = seq(10, 140, 10)) +
  labs(x = "Molecular Weight", y = "Fluorescence", title = "Protein Capillary Electrophoresis") +
  theme(panel.grid = element_blank()) # Remove the x and y lines
        #text = element_text(family = "Manrope")) # Plot the ribbon and mean lines
ggsave(path = "~/dev/R/Project_periodontitis/Results", "plot.png", p_plot_1, width = 10, height = 5, dpi = 300)



## Profiles data normalization---------------------------
#"Normalizing data with negative values can be a bit tricky."
# A common normalization factor for CE data is the total ion current (TIC), which represents the sum of all peak heights in a sample. 
tic <- rowSums(profiles_grouped[,3:ncol(profiles_grouped)]) # Calculate the normalization factor: To calculate the TIC, use the rowSums function in R to sum the peak heights across all variables for each sample
data_norm <- profiles_grouped[,3:ncol(profiles_grouped)] / tic # Divide each variable in the data frame by the TIC to obtain normalized values.
data_norm <- cbind(profiles_grouped[,1:2], data_norm) # Save the normalized data: add the sample ID column back
colnames(data_norm) <- colnames(profiles_grouped) # Save the normalized data: copy the column names

# Melt the data_norm dataframe
profiles_melt_tic <- melt(data_norm)

# Reorder the groups based on name
profiles_melt_tic$group <- factor(profiles_melt_tic$group, levels = c("Healthy",
                                                                      "Soft_Tissues",
                                                                      "Hard_Tissues"))

# Convert the 'value' column to a numeric data type
profiles_melt_tic$variable <- as.numeric(as.character(profiles_melt_tic$variable))

# Round the 'value' column to integers
profiles_melt_tic$rounded_variable <- as.integer(round(profiles_melt_tic$variable))

#create a summary table with mean and sd to be used on the geom_ribbon chart
summary_table_tic <- profiles_melt_tic %>%
  group_by(group, rounded_variable) %>% 
  summarize(value.mean = mean(value),
            value.sd = sd(value))

# Profiles visualization
p_plot_tic <- ggplot(summary_table_tic, aes(x = rounded_variable,
                                            y = rounded_variable,
                                            fill = group)) +
  geom_ribbon(aes(ymin = value.mean - value.sd,ymax = value.mean + value.sd), alpha = 0.3) +
  geom_line(aes(y = value.mean, color = group), linetype = "solid") +
  scale_fill_manual(values = c("darkgreen", "orange","brown")) +
  scale_color_manual(values = c("darkgreen", "orange","brown")) +
  coord_cartesian(xlim = c(10, 140)) +
  scale_x_continuous(breaks = seq(10, 140, 10)) +
  labs(x = "Molecular Weight", y = "Normalized Scale Intensity", title = "Protein Capillary Electrophoresis") +
  theme(panel.grid = element_blank()) # Remove the x and y lines
#text = element_text(family = "Manrope")) # Plot the ribbon and mean lines
ggsave(path = "~/dev/R/Project_periodontitis/Results", "plot_tic.png", p_plot_tic, width = 10, height = 5, dpi = 300)

#####
# Optional: Log transformation. 
# data_log <- log(data_norm + 1) 
# data_log <- cbind(profiles_1[,1], data_log)
# colnames(data_log) <- colnames(profiles_1)


# profiles_melt_tic1 <- melt(data_log)



# profiles_melt_tic1$group <- factor(profiles_melt_tic1$group, levels = c("Healthy",
#                                                                         "Gingivitis",
#                                                                         "Mucositis",
#                                                                         "Periimplantitis",
#                                                                         "Gingivitis_Mucositis_Periimplantitis",
#                                                                         "Periodontitis_Periimplantitis",
#                                                                         "Mucositis_Periimplantitis",
#                                                                         "MucositisePeriimplantitis",
#                                                                         "Gingivitis_Periimplantitis"))

# Convert the 'value' column to a numeric data type
# profiles_melt_tic1$variable <- as.numeric(as.character(profiles_melt_tic1$variable))

# Round the 'value' column to integers
# profiles_melt_tic1$rounded_variable <- as.integer(round(profiles_melt_tic1$variable))


#create a summary table with mean and sd to be used on the geom_ribbon chart
# summary_table_tic1 <- profiles_melt_tic1 %>%
#   group_by(group, rounded_variable) %>% 
#   summarize(value.mean = mean(value),
#             value.sd = sd(value))
# 
# 
# p_plot_tic1 <- ggplot(summary_table_tic1, aes(x = rounded_variable,
#                                             y = rounded_variable,
#                                             fill = group)) +
#   geom_ribbon(aes(ymin = value.mean - value.sd,ymax = value.mean + value.sd), alpha = 0.3) +
#   geom_line(aes(y = value.mean, color = group), linetype = "solid") +
#   scale_fill_manual(values = c("green", "blue", "orange","red")) +
#   scale_color_manual(values = c("green", "blue", "orange","red")) +
#   coord_cartesian(xlim = c(10, 140)) +
#   scale_x_continuous(breaks = seq(10, 140, 10)) +
#   labs(x = "Molecular Weight", y = "Normalized Scale Intensity", title = "Protein Capillary Electrophoresis") +
#   theme(panel.grid = element_blank()) # Remove the x and y lines
#         #text = element_text(family = "Manrope")) # Plot the ribbon and mean lines
# ggsave("plot_tic1.png", p_plot_tic1, width = 10, height = 5, dpi = 300)

#---

# One option is to use the scale() function to scale the y-values to a specified range. 
## This fuction below normalize each y-column to a range of 0 to 1:

normalize_ce_data <- function(data) {
  # Normalize each y-axis to range of 0 to 1
  y_cols <- colnames(data)[-1]
  y_norm <- lapply(data[y_cols], function(y) {
    scale(y, center = FALSE, scale = max(abs(y), na.rm = TRUE)) + 1
  })
  
  # Merge normalized data by x-axis
  data_norm <- cbind(data[1], y_norm)
  
  # Set column names
  colnames(data_norm)[-1] <- paste0("y_", seq_len(ncol(data_norm)-1), "_norm")
  
  # Return normalized data
  return(data_norm)
}

# Normalize CE data
ce_data_norm <- normalize_ce_data(profiles_1)

# Replace column names of normalized data frame with original column names
names(ce_data_norm) <- names(profiles_1)


profiles_melt_norm <- melt(ce_data_norm) 

# Reorder the groups based on name
profiles_melt_norm$group <- factor(profiles_melt$group, levels = c("Healthy",
                                                                   "Gingivitis",
                                                                   "Mucositis",
                                                                   "Periimplantitis"))

# Convert the 'value' column to a numeric data type
profiles_melt_norm$variable <- as.numeric(as.character(profiles_melt_norm$variable))

# Round the 'value' column to integers
profiles_melt_norm$rounded_variable <- as.integer(round(profiles_melt_norm$variable))

summary_table_norm <- profiles_melt_norm %>%
  group_by(group, rounded_variable) %>% 
  summarize(value.mean = mean(value),
            value.sd = sd(value))


p_plot_norm <- ggplot(summary_table_norm, aes(x = rounded_variable,
                                      y = rounded_variable,
                                      fill = group)) +
  geom_ribbon(aes(ymin = value.mean - value.sd,ymax = value.mean + value.sd), alpha = 0.3) +
  geom_line(aes(y = value.mean, color = group), linetype = "solid") +
  scale_fill_manual(values = c("green", "blue", "orange","red")) +
  scale_color_manual(values = c("green", "blue", "orange","red")) +
  coord_cartesian(xlim = c(10, 140)) +
  scale_x_continuous(breaks = seq(10, 140, 10)) +
  labs(x = "Molecular Weight", y = "Normalized Scale Intensity", title = "Protein Capillary Electrophoresis") +
  theme(panel.grid = element_blank(), # Remove the x and y lines
        text = element_text(family = "Manrope")) # Plot the ribbon and mean lines
ggsave("plot_norm.png", p_plot_norm, width = 10, height = 5, dpi = 300)


# Statistical analysis ####

## The data - Define predictor (X) and response (Y) variables ----
plsda <- data_norm[,3:399] # Select the MW columns
colnames(plsda) <- paste("mw", colnames(plsda), sep="_") # Add "mw_" to the column names 
X <- plsda # Assign plsda df to X
data_norm$group <- factor(data_norm$group, levels = c("Healthy",
                                                      "Soft_Tissues",
                                                      "Hard_Tissues")) # Make the group column a factor
Y <- data_norm$group # Assign plsda df to Y

dim(X) # check the dimensions of the X dataframe
summary(Y) # check the distribution of class labels


## Initial Analysis - Preliminary Analysis with PCA ----

# run pca method on data
pca.srbct <- pca(X, ncomp = 10, center = TRUE, scale = TRUE)
plot(pca.srbct)  # barplot of the eigenvalues (explained variance per component)

plotIndiv(pca.srbct, group = data_norm$group, ind.names = FALSE, # plot the samples projected
          legend = TRUE, title = 'PCA on SRBCT, comp 1 - 2') # onto the PCA subspace


## Initial sPLS-DA model ----
srbct.splsda <- splsda(X, Y, ncomp = 10)

# plot the samples projected onto the first two components of the PLS-DA subspace
plotIndiv(srbct.splsda , comp = 1:2, 
          group = data_norm$group, ind.names = FALSE,  # colour points by class
          ellipse = TRUE, # include 95% confidence ellipse for each class
          legend = TRUE, title = '(a) PLSDA with confidence ellipses')
# use the max.dist measure to form decision boundaries between classes based on PLS-DA data
background = background.predict(srbct.splsda, comp.predicted=2, dist = "max.dist")
# plot the samples projected onto the first two components of the PLS-DA subspace
plotIndiv(srbct.splsda, comp = 1:2,
          group = data_norm$group, ind.names = FALSE, # colour points by class
          background = background, # include prediction background for each class
          legend = TRUE, title = " (b) PLSDA with prediction background")


## Tuning sPLS-DA ----
# undergo performance evaluation in order to tune the number of components to use

perf.splsda.srbct <- perf(srbct.splsda, validation = "loo", 
                          folds = 5, nrepeat = 100, # use repeated cross-validation
                          progressBar = FALSE, auc = TRUE) # include AUC values

# plot the outcome of performance evaluation across all ten components
plot(perf.splsda.srbct, col = color.mixo(5:7), sd = TRUE,
     legend.position = "vertical")
#legend('topright', 
#       legend = c("overall", "BER", "max.dist", "centroids.dist", "mahalanobis.dist"), 
#       fill = color.mixo(5:7))


perf.splsda.srbct$choice.ncomp # what is the optimal value of components according to perf()



## Selecting the number of variables ----
### The keepX Parameter ====

# grid of possible keepX values that will be tested for each component
list.keepX <- c(1:10,  seq(20, 300, 10))

# undergo the tuning process to determine the optimal number of variables
tune.splsda.srbct <- tune.splsda(X, Y, ncomp = 2, # calculate for first 4 components
                                 validation = 'Mfold',
                                 folds = 5, nrepeat = 3, # use repeated cross-validation
                                 dist = 'max.dist', # use max.dist measure
                                 measure = "BER", # use balanced error rate of dist measure
                                 test.keepX = list.keepX,
                                 cpus = 2) # allow for paralleliation to decrease runtime

plot(tune.splsda.srbct, col = color.jet(2)) # plot output of variable number tuning

tune.splsda.srbct$choice.ncomp$ncomp # what is the optimal value of components according to tune.splsda()

tune.splsda.srbct$choice.keepX # what are the optimal values of variables according to tune.splsda()
optimal.ncomp <- tune.splsda.srbct$choice.ncomp$ncomp
optimal.keepX <- tune.splsda.srbct$choice.keepX[1:optimal.ncomp]

## Final Model ----

# form final model with optimised values for component and variable count
final.splsda <- splsda(X, Y, 
                       ncomp = 2, # The "optimal.ncomp" identify only 1 component ideal according to the tune.splda
                       keepX = optimal.keepX)


## Plots ----

plotIndiv(final.splsda, comp = c(1,2), # plot samples from final model
          group = data_norm$group, ind.names = FALSE, # colour by class label
          ellipse = TRUE, legend = TRUE, # include 95% confidence ellipse
          title = ' (a) sPLS-DA on SRBCT, comp 1 & 2')

plotIndiv(final.splsda, comp = c(1,3), # plot samples from final model
          group = profiles_1$group, ind.names = FALSE,  # colour by class label
          ellipse = TRUE, legend = TRUE, # include 95% confidence ellipse
          title = '(b) sPLS-DA on SRBCT, comp 1 & 3')


# set the styling of the legend to be homogeneous with previous plots
legend=list(legend = levels(Y), # set of classes
            col = unique(color.mixo(Y)), # set of colours
            title = "Profile grouping", # legend title
            cex = 0.7) # legend size

# generate the CIM, using the legend and colouring rows by each sample's class
cim <- cim(final.splsda, row.sideColors = color.mixo(Y), 
           legend = legend)
ggsave("cim.png", cim, width = 10, height = 5, dpi = 300)


## Variable Plots ----

# form new perf() object which utilises the final model
perf.splsda.srbct <- perf(final.splsda, 
                          folds = 5, nrepeat = 3, # use repeated cross-validation
                          validation = "Mfold", dist = "max.dist",  # use max.dist measure
                          progressBar = FALSE)

# plot the stability of each feature for the first three components, 'h' type refers to histogram
par(mfrow=c(1,3))
plot(perf.splsda.srbct$features$stable[[1]], type = 'h', 
     ylab = 'Stability', 
     xlab = 'Features', 
     main = '(a) Comp 1', las =2)
plot(perf.splsda.srbct$features$stable[[2]], type = 'h', 
     ylab = 'Stability', 
     xlab = 'Features', 
     main = '(b) Comp 2', las =2)
plot(perf.splsda.srbct$features$stable[[3]], type = 'h', 
     ylab = 'Stability', 
     xlab = 'Features',
     main = '(c) Comp 3', las =2)

perf.splsda.srbct$features$stable[[1]]
perf.splsda.srbct$features$stable[[2]]
perf.splsda.srbct$features$stable[[3]]

variable_plot_feature <- perf.splsda.srbct$features$stable[[1]]
variable_plot_feature <- t(data.frame(rbind(variable_plot_feature)))
write_csv(variable_plot_feature, file = "perf_splsda_srbct_features.csv")

var.name.short <- substr(profiles_1[,2:398], 1, 10) # form simplified gene names

plotVar(final.splsda, comp = c(1,2), var.names = list(var.name.short), cex = 3) # generate correlation circle plot


## Prediction ----

train <- sample(1:nrow(X), 10) # randomly select 10 samples in training
test <- setdiff(1:nrow(X), train) # rest is part of the test set

# store matrices into training and test set:
X.train <- X[train, ]
X.test <- X[test,]
Y.train <- Y[train]
Y.test <- Y[test]

# train the model
train.splsda.srbct <- splsda(X.train, Y.train, ncomp = optimal.ncomp, keepX = optimal.keepX)

# use the model on the Xtest set
predict.splsda.srbct <- predict(train.splsda.srbct, X.test, 
                                dist = "max.dist")

# evaluate the prediction accuracy for the first two components
predict.comp2 <- predict.splsda.srbct$class$max.dist[,1]
table(factor(predict.comp2, levels = levels(Y)), Y.test)

predict.comp2_table <- t(data.frame(rbind(table(factor(predict.comp2, levels = levels(Y)), Y.test))))


## Performance Plots ----
par(mfrow=c(1,2))
auc.splsda = auroc(final.splsda, roc.comp = 1, print = TRUE) # AUROC for the first component

auc.splsda = auroc(final.splsda, roc.comp = 3, print = TRUE) # AUROC for all three components









###############################################################

profiles_group <- read.table(file = '/Users/eduardoesteves/git/R/Project_periodontitis/Data/profiles_group.csv',
                             header = TRUE, sep = ',', quote = "", row.names = NULL, 
                             stringsAsFactors = FALSE, check.names = FALSE)
profiles_group <- dplyr::select(profiles_group, -c("id"))

## Statistical analysis with healthy/Hard/Soft/Mix groups####


plsda_group <- profiles_group[,2:398] # Select the MW columns
colnames(plsda) <- paste("mw", colnames(plsda), sep="_") # Add "mw_" to the column names 
X_group <- plsda_group # Assign plsda df to X
profiles_group$group <- factor(profiles_group$group) # Make the group column a factor
Y_group <- profiles_group$group # Assign plsda df to Y

dim(X_group) # check the dimensions of the X dataframe
summary(Y_group) # check the distribution of class labels


## Initial Analysis - Preliminary Analysis with PCA ----

# run pca method on data
pca.srbct_grouped = pca(X_group, ncomp = 10, center = TRUE, scale = TRUE) 
plot(pca.srbct_grouped)  # barplot of the eigenvalues (explained variance per component)

plotIndiv(pca.srbct_grouped, group = profiles_group$group, ind.names = FALSE, # plot the samples projected
          legend = TRUE, title = 'PCA on SRBCT, comp 1 - 2') # onto the PCA subspace


## Initial sPLS-DA model ----
srbct.splsda_grouped <- splsda(X_group, Y_group, ncomp = 10)

# plot the samples projected onto the first two components of the PLS-DA subspace
plotIndiv(srbct.splsda_grouped , comp = 1:2, 
          group = profiles_group$group, ind.names = FALSE,  # colour points by class
          ellipse = TRUE, # include 95% confidence ellipse for each class
          legend = TRUE, title = '(a) PLSDA with confidence ellipses')
# use the max.dist measure to form decision boundaries between classes based on PLS-DA data
background_grouped = background.predict(srbct.splsda_grouped, comp.predicted=2, dist = "max.dist")
# plot the samples projected onto the first two components of the PLS-DA subspace
plotIndiv(srbct.splsda_grouped, comp = 1:2,
          group = profiles_group$group, ind.names = FALSE, # colour points by class
          background = background_grouped, # include prediction background for each class
          legend = TRUE, title = " (b) PLSDA with prediction background")


## Tuning sPLS-DA ----
# undergo performance evaluation in order to tune the number of components to use

perf.splsda.srbct_group_loo <- perf(srbct.splsda_grouped, validation = "loo", 
                                folds = 1, nrepeat = 100, # use repeated cross-validation
                                progressBar = FALSE, auc = TRUE) # include AUC values

perf.splsda.srbct_group_mfold <- perf(srbct.splsda_grouped, validation = "Mfold", 
                                    folds = 5, nrepeat = 100, # use repeated cross-validation
                                    progressBar = FALSE, auc = TRUE) # include AUC values                          

# plot the outcome of performance evaluation across all ten components
plot(perf.splsda.srbct_group_loo, col = color.mixo(5:7), sd = TRUE,
     legend.position = "vertical")

plot(perf.splsda.srbct_group_mfold, col = color.mixo(5:7), sd = TRUE,
     legend.position = "vertical")

perf.splsda.srbct_group_loo$choice.ncomp # what is the optimal value of components according to perf()
perf.splsda.srbct_group_mfold$choice.ncomp


## Selecting the number of variables ----
### The keepX Parameter ====

# grid of possible keepX values that will be tested for each component
list.keepX_group <- c(1:10,  seq(20, 300, 10))

# undergo the tuning process to determine the optimal number of variables
tune.splsda.srbct_group_mfold <- tune.splsda(X_group, Y_group, ncomp = 3, # calculate for first 4 components
                                 validation = 'Mfold',
                                 folds = 5, nrepeat = 100, # use repeated cross-validation
                                 dist = 'max.dist', # use max.dist measure
                                 measure = "BER", # use balanced error rate of dist measure
                                 test.keepX = list.keepX_group,
                                 cpus = 2) # allow for paralleliation to decrease runtime

plot(tune.splsda.srbct_group_mfold, col = color.jet(3)) # plot output of variable number tuning

tune.splsda.srbct_group_mfold$choice.ncomp$ncomp # what is the optimal value of components according to tune.splsda()

tune.splsda.srbct_group_mfold$choice.keepX # what are the optimal values of variables according to tune.splsda()
optimal.ncomp_group <- tune.splsda.srbct_group_mfold$choice.ncomp$ncomp
optimal.keepX_group <- tune.splsda.srbct_group_mfold$choice.keepX[1:optimal.ncomp_group]



#tune.splsda.srbct_group_loo <- tune.splsda(X_group, Y_group, ncomp = 2, # calculate for first 4 components
                                             #validation = 'loo',
                                             #folds = 5, nrepeat = 100, # use repeated cross-validation
                                             #dist = 'max.dist', # use max.dist measure
                                             #measure = "BER", # use balanced error rate of dist measure
                                             #test.keepX = list.keepX_group,
                                             #cpus = 2) # allow for paralleliation to decrease runtime


#plot(tune.splsda.srbct_group_loo, col = color.jet(2)) # plot output of variable number tuning

#tune.splsda.srbct_group_loo$choice.ncomp$ncomp # what is the optimal value of components acco

## Final Model ----

# form final model with optimised values for component and variable count
final.splsda_group <- splsda(X_group, Y_group, 
                       ncomp = 2, # The "optimal.ncomp" identify only 1 component ideal according to the tune.splda
                       keepX = optimal.keepX_group)


## Plots ----

plotIndiv(final.splsda_group, comp = c(1,2), # plot samples from final model
          group = profiles_group$group, ind.names = FALSE, # colour by class label
          ellipse = TRUE, legend = TRUE, # include 95% confidence ellipse
          title = ' (a) sPLS-DA on SRBCT, comp 1 & 2')

plotIndiv(final.splsda, comp = c(1,3), # plot samples from final model
          group = profiles_1$group, ind.names = FALSE,  # colour by class label
          ellipse = TRUE, legend = TRUE, # include 95% confidence ellipse
          title = '(b) sPLS-DA on SRBCT, comp 1 & 3')


# set the styling of the legend to be homogeneous with previous plots
legend_group <- list(legend = levels(Y_group), # set of classes
            col = unique(color.mixo(Y_group)), # set of colours
            title = "Profile grouping", # legend title
            cex = 0.7) # legend size

# generate the CIM, using the legend and colouring rows by each sample's class
cim <- cim(final.splsda, row.sideColors = color.mixo(Y), 
           legend = legend)
ggsave("cim.png", cim, width = 10, height = 5, dpi = 300)


## Variable Plots ----

# form new perf() object which utilises the final model
perf.splsda.srbct_varplots_grouped <- perf(final.splsda_group, 
                          folds = 5, nrepeat = 100, # use repeated cross-validation
                          validation = "Mfold", dist = "max.dist",  # use max.dist measure
                          progressBar = TRUE)

# plot the stability of each feature for the first three components, 'h' type refers to histogram
par(mfrow=c(1,3))
plot(perf.splsda.srbct_varplots_grouped$features$stable[[1]], type = 'h', 
     ylab = 'Stability', 
     xlab = 'Features', 
     main = '(a) Comp 1', las =2)
plot(perf.splsda.srbct$features$stable[[2]], type = 'h', 
     ylab = 'Stability', 
     xlab = 'Features', 
     main = '(b) Comp 2', las =2)
plot(perf.splsda.srbct$features$stable[[3]], type = 'h', 
     ylab = 'Stability', 
     xlab = 'Features',
     main = '(c) Comp 3', las =2)

perf.splsda.srbct$features$stable[[1]]
perf.splsda.srbct$features$stable[[2]]
perf.splsda.srbct$features$stable[[3]]

variable_plot_feature <- perf.splsda.srbct$features$stable[[1]]
variable_plot_feature <- t(data.frame(rbind(variable_plot_feature)))
write_csv(variable_plot_feature, file = "perf_splsda_srbct_features.csv")

var.name.short <- substr(profiles_1[,2:398], 1, 10) # form simplified gene names

plotVar(final.splsda, comp = c(1,2), var.names = list(var.name.short), cex = 3) # generate correlation circle plot


## Prediction ----

train <- sample(1:nrow(X), 10) # randomly select 10 samples in training
test <- setdiff(1:nrow(X), train) # rest is part of the test set

# store matrices into training and test set:
X.train <- X[train, ]
X.test <- X[test,]
Y.train <- Y[train]
Y.test <- Y[test]

# train the model
train.splsda.srbct <- splsda(X.train, Y.train, ncomp = optimal.ncomp, keepX = optimal.keepX)

# use the model on the Xtest set
predict.splsda.srbct <- predict(train.splsda.srbct, X.test, 
                                dist = "max.dist")

# evaluate the prediction accuracy for the first two components
predict.comp2 <- predict.splsda.srbct$class$max.dist[,1]
table(factor(predict.comp2, levels = levels(Y)), Y.test)

predict.comp2_table <- t(data.frame(rbind(table(factor(predict.comp2, levels = levels(Y)), Y.test))))


## Performance Plots ----
par(mfrow=c(1,2))
auc.splsda = auroc(final.splsda, roc.comp = 1, print = TRUE) # AUROC for the first component

auc.splsda = auroc(final.splsda, roc.comp = 3, print = TRUE) # AUROC for all three components

