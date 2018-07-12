# File: 05_countMatrixEDA.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: quality checks on the count matrix with covariates
# Date: 11/07/2018

source('header.R')

## load the data
library(RMySQL)

db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
dbListTables(db)

# get the query
g_did
q = paste0('select MetaFile.* from MetaFile
           where (MetaFile.idData = 28) AND (MetaFile.comment like "%count%")')
dfSample = dbGetQuery(db, q)
dfSample
n = paste0(dfSample$location, dfSample$name)
load(n)

## load the metadata i.e. covariates
q = paste0('select Sample.* from Sample where Sample.idData = 28')
dfSample = dbGetQuery(db, q)
dim(dfSample)
dfSample
# close connection after getting data
dbDisconnect(db)

## make count matrix
names(lCounts)
mCounts = do.call(cbind, lCounts)
colnames(mCounts) = gsub('.bam_\\w+.bam', '', colnames(mCounts))

# reorder the count matrix columns according to the order in samples table
i = match(dfSample$title, colnames(mCounts))
mCounts = mCounts[,i]
# sanity check
identical(dfSample$title, colnames(mCounts))

mData = mCounts
dim(mData)

# drop the samples where average across rows is less than 3
i = rowMeans(mData)
table( i < 3)
mData = mData[!(i< 3),]
dim(mData)

library(DESeq2)
sf = estimateSizeFactorsForMatrix(mData)
mData.norm = sweep(mData, 2, sf, '/')

## load CDiagnostics and test
## compare the normalised and raw data
library(downloader)
url = 'https://raw.githubusercontent.com/uhkniazi/CDiagnosticPlots/master/CDiagnosticPlots.R'
download(url, 'CDiagnosticPlots.R')

# load the required packages
source('CDiagnosticPlots.R')
# delete the file after source
unlink('CDiagnosticPlots.R')

oDiag.1 = CDiagnosticPlots(log(mData.norm+1), 'Normalised')
oDiag.2 = CDiagnosticPlots(log(mData+1), 'Original')

# the batch variable we wish to colour by, 
# this can be any grouping/clustering in the data capture process
fBatch = factor(dfSample$group1)
fBatch = cut(as.numeric(dfSample$group2), 3, labels = paste0('time', 1:3))
fBatch = factor(dfSample$group3)
# well labels
lab = sapply(strsplit(dfSample$description, ';'), function(x) return(x[4]))

pdf('results/matrixClustering.pdf')

## compare the 2 methods using various plots
par(mfrow=c(1,2))
boxplot.median.summary(oDiag.1, fBatch, legend.pos = 'topright', axis.label.cex = 0.7)
boxplot.median.summary(oDiag.2, fBatch, legend.pos = 'topright', axis.label.cex = 0.7)

plot.mean.summary(oDiag.1, fBatch, axis.label.cex = 0.7)
plot.mean.summary(oDiag.2, fBatch, axis.label.cex = 0.7)

plot.sigma.summary(oDiag.1, fBatch, axis.label.cex = 0.7)
plot.sigma.summary(oDiag.2, fBatch, axis.label.cex = 0.7)

plot.missing.summary(oDiag.1, fBatch, axis.label.cex = 0.7, cex.main=1)
plot.missing.summary(oDiag.2, fBatch, axis.label.cex = 0.7, cex.main=1)

plot.PCA(oDiag.1, fBatch, cex.main=1)
plot.PCA(oDiag.2, fBatch, cex.main=1)
# 
plot.dendogram(oDiag.1, fBatch, labels_cex = 0.8, cex.main=0.7)
plot.dendogram(oDiag.2, fBatch, labels_cex = 0.8, cex.main=0.7)

## change parameters 
l = CDiagnosticPlotsGetParameters(oDiag.1)
l$PCA.jitter = F
l$HC.jitter = F

oDiag.1 = CDiagnosticPlotsSetParameters(oDiag.1, l)
oDiag.2 = CDiagnosticPlotsSetParameters(oDiag.2, l)
plot.PCA(oDiag.1, fBatch)
plot.PCA(oDiag.2, fBatch)
plot.PCA(oDiag.1, fBatch, csLabels = dfSample$group2)
plot.dendogram(oDiag.1, fBatch, labels_cex = 0.7)
plot.dendogram(oDiag.2, fBatch, labels_cex = 0.7)

## extreme values
oDiag.1 = calculateExtremeValues(oDiag.1)
oDiag.2 = calculateExtremeValues(oDiag.2)
m1 = mGetExtremeValues(oDiag.1)
m2 = mGetExtremeValues(oDiag.2)

## samples with most extreme values
apply(m1, 2, function(x) sum(x > 0))
apply(m2, 2, function(x) sum(x > 0))

## variables that are contributing to this
v1 = apply(m1, 1, function(x) sum(x > 0))
v2 = apply(m2, 1, function(x) sum(x > 0))

which(v1 > 0)
which(v2 > 0)

############ make a lattice plot
library(lattice)
## spikein normalisation
sf2 = estimateSizeFactorsForMatrix(mData, controlGenes = grep('ERCC', x = rownames(mData)))
mData.ercc = sweep(mData, 2, sf2, '/')

ivMean.spike = colMeans(mData.ercc)
ivMean.gene = colMeans(mData.norm)
ivMean.raw = colMeans(mData)

dfData = data.frame(Spike=ivMean.spike, Gene=ivMean.gene, Raw=ivMean.raw, treatment=dfSample$group1, Title=dfSample$title)
dfData.st = stack(dfData[,1:3])
dfData.st$treatment = dfData$treatment
dfData.st$Title = factor(dfData$Title)


xyplot(log(values) ~ Title, data=dfData.st, auto.key=list(columns=3), main='Average Gene Expression in Each Sample', type='b', groups=ind,
        par.settings=list(superpose.line = list(lwd=0.5, lty=1)),
        xlab='', ylab='Mean Expression', pch=20, scales=list(x=list(cex=0.7, rot=45)))

xyplot((values) ~ Title, data=dfData.st, auto.key=list(columns=3), main='Average Gene Expression in Each Sample', type='b', groups=ind,
       par.settings=list(superpose.line = list(lwd=0.5, lty=1)),
       xlab='', ylab='Mean Expression', pch=20, scales=list(x=list(cex=0.7, rot=45)))
