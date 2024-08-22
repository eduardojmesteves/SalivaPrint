library(pathfindR)
library(BiocManager)
library(tidyverse)
library(extrafont)
library(ggplot2)
library(mixOmics)
options(rgl.useNULL = TRUE)
library(rgl)

df <-read.table(file = '/Users/ee/git/r/Project_profile/proteins_dataset.csv', 
                header = TRUE, sep = ',', quote = "", row.names = NULL, 
                stringsAsFactors = FALSE, check.names = FALSE)

df2 <-read.csv(file = '/Users/ee/git/r/Project_profile/factor.csv', header = F)

#df$oral_diagnosis = as.factor(df$oral_diagnosis)
str(df2)
Y1 <- as.factor(df2$V1)

X1 = subset(df, select = -c(oral_diagnosis,sample) )


MyResult.splsda <- splsda(X1, Y1, keepX = c(394,394)) # 1 Run the method
plotIndiv(MyResult.splsda)                          # 2 Plot the samples (coloured by classes automatically)

plotVar(MyResult.splsda)                            # 3 Plot the variables

plotIndiv(MyResult.splsda, ind.names = FALSE, legend=TRUE,
          ellipse = TRUE, star = TRUE, title = 'sPLS-DA on SRBCT',
          X.label = 'PLS-DA 1', Y.label = 'PLS-DA 2')

MyResult.plsda <- plsda(X1, Y1) # 1 Run the method
plotIndiv(MyResult.plsda)       # 2 Plot the samples (coloured by classes automatically)

plotVar(MyResult.plsda, cutoff = 0.7)        # 3 Plot the variables

plotIndiv(MyResult.plsda, ind.names = FALSE, legend=TRUE,
          ellipse = TRUE, star = TRUE, title = 'sPLS-DA on SRBCT',
          X.label = 'PLS-DA 1', Y.label = 'PLS-DA 2')


# background <- background.predict(MyResult.splsda, comp.predicted=2,
#                                  dist = "max.dist") 
# plotIndiv(MyResult.splsda, comp = 1:2, group = srbct$class,
#           ind.names = FALSE, title = "Maximum distance",
#           legend = TRUE,  background = background)

MyResult.splsda2 <- splsda(X1,Y1, ncomp=3, keepX=c(15,10,5))
selectVar(MyResult.splsda2, comp=1)$value
plotLoadings(MyResult.splsda2, contrib = 'max', method = 'mean')
plotIndiv(MyResult.splsda2, style="3d")

MyResult.plsda2 <- plsda(X1,Y1, ncomp=10)
set.seed(30) # for reproducbility in this vignette, otherwise increase nrepeat
MyPerf.plsda <- perf(MyResult.plsda2, validation = "Mfold", folds = 3, 
                     progressBar = FALSE, nrepeat = 10) # we suggest nrepeat = 50
plot(MyPerf.plsda, col = color.mixo(5:7), sd = TRUE, legend.position = "horizontal")
