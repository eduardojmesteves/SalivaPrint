library(pathfindR)
library(BiocManager)
library(tidyverse)
library(extrafont)
library(ggplot2)
library(mixOmics)

df <-read.table(file = '/Users/ee/Downloads/proteins_dataset.csv', 
                header = TRUE, sep = ',', quote = "", row.names = NULL, 
                stringsAsFactors = FALSE, check.names = FALSE)

df2 <-read.csv(file = '/Users/ee/Downloads/Y2.csv', header = F)

df$oral_diagnosis = as.factor(df$oral_diagnosis)
str(df2)
xx <- as.factor(df2$V1)

X1 = subset(df, select = -c(oral_diagnosis,sample) )
Y1 = subset(df, select = -c(oral_diagnosis) )

Y1 <-read.table(file = '/Users/ee/Downloads/vecto.csv', 
                header = TRUE, sep = ',', quote = "", row.names = NULL, 
                stringsAsFactors = FALSE, check.names = FALSE)

Y11 <- as.data.frame(unclass(Y1),stringsAsFactors=TRUE)

Y1t <- setNames(as.data.frame(t(Y1[-1])), Y1[[1]])
vec1 <- Y1$healthy                
df$D <- as.factor(Y1$healthy)    




MyResult.splsda <- splsda(X1, xx, keepX = c(50,50)) # 1 Run the method
plotIndiv(MyResult.splsda)                          # 2 Plot the samples (coloured by classes automatically)
plotVar(MyResult.splsda)                            # 3 Plot the variables

plotIndiv(MyResult.splsda, ind.names = FALSE, legend=TRUE,
          ellipse = TRUE, star = TRUE, title = 'sPLS-DA on SRBCT',
          X.label = 'PLS-DA 1', Y.label = 'PLS-DA 2')



background <- background.predict(MyResult.splsda, comp.predicted=2,
                                 dist = "max.dist") 
plotIndiv(MyResult.splsda, comp = 1:2, group = srbct$class,
          ind.names = FALSE, title = "Maximum distance",
          legend = TRUE,  background = background)

