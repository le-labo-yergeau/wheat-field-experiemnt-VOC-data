setwd("~/Dropbox/Work/INRS/Data/VOC/Field_exp_wheat")

library(MetaboAnalystR)
library(mixOmics)
library(data.table)
# 1. normalize data using metaboanalyst
#laod data exported from aligned list to metaboanalyst (file containing groups defined in mzMine and unqiue sample names)
mSet<-InitDataObjects("conc", "stat", FALSE)
mSet<-Read.TextData(mSet, "Field_exp_wheat_alignment1_metabo.csv", "colu", "disc")
mSet<-SanityCheckData(mSet)
mSet<-RemoveMissingPercent(mSet, percent=0.5)
mSet<-ImputeVar(mSet, method="min")
mSet<-FilterVariable(mSet, "iqr", "F", 25)
mSet<-PreparePrenormData(mSet)
#normalizaition (sample norm: by median, data transform: log transf, data scaling: auto)
mSet<-Normalization(mSet, "MedianNorm", "LogNorm", "AutoNorm", ratio=FALSE, ratioNum=20)
#save transformed data in wd
mSet<-PlotNormSummary(mSet, "norm_0_", "png", 72, width=NA)
mSet<-PlotSampleNormSummary(mSet, "snorm_0_", "png", 72, width=NA)
mSet<-SaveTransformedData(mSet)
mSet<-ANOVA.Anal(mSet, F, 0.05, "fisher", FALSE)
mSet<-PlotANOVA(mSet, "aov_0_", "png", 72, width=NA)
mSet<-PCA.Anal(mSet)
mSet<-PlotPCAScree(mSet, "pca_scree_0_", "png", 72, width=NA, 5)
mSet<-PlotPCA2DScore(mSet, "pca_score2d_0_", "png", 72, width=NA, 1,2,0.95,0,0)
mSet<-PLSR.Anal(mSet, reg=TRUE)
mSet<-PlotPLS2DScore(mSet, "pls_score2d_0_", "png", 72, width=NA, 1,2,0.95,0,0)

#load normalised data set (if necessary remove first row and add sample_id, e.g. 100_DS_5)
mset=data.frame(fread("data_normalized_labeled.csv", header=T))
#conversion to log 
#mset=log10(mset[c(2:49)])
head(mset)
row.names(mset)=mset[,1]
mset=mset[,-1]
mset= as.data.frame(t(mset))
mset=mset[order(rownames(mset)),]

#get the mapping file
treat=data.frame(fread("mapping_file.tsv", sep="\t", header=T))
row.names(treat)=treat[,1]
treat=treat[order(treat$SampleID), ]
treat=transform(treat,perc_SWHC=as.factor(perc_SWHC), block=as.factor(Block))
head(treat)

#sanity check
row.names(mset)==paste("X",treat$SampleID, sep="")

#pca 
data.pca = pca(mset)  # 1 Run the method
data.pca

#The amount of variance explained
plot(data.pca)

plotIndiv(data.pca) # 2 Plot the samples

plotVar(data.pca) # 3 Plot the variables
plotIndiv(data.pca, group = treat$perc_SWHC, 
          legend = TRUE, legend.title = 'perc_SWHC', title = 'PCA volatiles')
plotIndiv(data.pca, ind.names = FALSE,
          group = treat$perc_SWHC,
          pch = as.factor(treat$Cultivar),
          legend = TRUE, title = 'Volatiles field wheat',
          legend.title = 'perc_SWHC', legend.title.pch = 'Cultivar')

#plot loadings 
#plotLoadings(data.pca)

data.pca2
#The amount of variance explained
plot(data.pca2)

#plot loadings 
#plotLoadings(data.pca2)

#3d plot 
library(rgl)
plotIndiv(data.pca2,
          group = treat$perc_SWHC, style="3d",
          legend = TRUE, title = 'Volatiles-all, PCA comp 1 - 2 - 3')

#plsda SWHC
X = mset
Y = as.factor(treat$perc_SWHC)
summary(Y)
length(Y)
data.splsda = splsda(X, Y) # 1 Run the method
plotIndiv(data.splsda) # 2 Plot the samples
plotVar(data.splsda, var.names=FALSE) # 3 Plot the variables

plotIndiv(data.splsda, ind.names = FALSE, legend=TRUE,
          ellipse = TRUE, star = FALSE, title = 'Volatiles-perc_SWHC',
          X.label = 'PLS-DA 1', Y.label = 'PLS-DA 2')

data.splsda2 = splsda(X,Y, ncomp=3)
selectVar(data.splsda2, comp=1)$value
plotIndiv(data.splsda2, style="3d")
plotLoadings(data.splsda2, contrib = 'max', method = 'mean')

#plsda Cultivar
X = mset
Y = treat$Cultivar
summary(Y)
length(Y)
data.splsda = splsda(X, Y) # 1 Run the method
plotIndiv(data.splsda) # 2 Plot the samples
plotVar(data.splsda, var.names=FALSE) # 3 Plot the variables

plotIndiv(data.splsda, ind.names = FALSE, legend=TRUE,
          ellipse = TRUE, star = FALSE, title = 'Volatiles-cultivar',
          X.label = 'PLS-DA 1', Y.label = 'PLS-DA 2')

data.splsda2 = splsda(X,Y, ncomp=3)
selectVar(data.splsda2, comp=1)$value
plotIndiv(data.splsda2, style="3d")
plotLoadings(data.splsda2, contrib = 'max', method = 'mean')
```

