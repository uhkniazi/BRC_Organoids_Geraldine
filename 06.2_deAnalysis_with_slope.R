# File: 06_deAnalysis.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: Modelling and selecting DE genes
# Date: 12/07/2018

## load the data
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
tail(rownames(mData))
## remove the ercc genes
i = grep('ERCC', x = rownames(mData))
mData = mData[-i,]

library(DESeq2)
sf = estimateSizeFactorsForMatrix(mData)
mData.norm = sweep(mData, 2, sf, '/')
dim(mData.norm)
str(mData.norm)
## perform DE analysis
## delete sample section after testing
mData.norm = round(mData.norm, 0)
#set.seed(123)
# i = sample(1:nrow(mData.norm), 100, replace = F)
# dfData = data.frame(t(mData.norm[i,]))

dfData = data.frame(t(mData.norm))
dfData = stack(dfData)
dim(dfData)
dfData$fBatch = factor(dfSample$group1)
dfData$fAdjust1 = factor(dfSample$group3)
dfData$time = as.numeric(dfSample$group2)
dfData$Coef = factor(dfData$fBatch:dfData$ind)
dfData$Coef.adj1 = factor(dfData$fAdjust1:dfData$ind)
## adjust the time as a continuous variable in each ind/gene
dim(dfData)
dfData = droplevels.data.frame(dfData)
dfData = dfData[order(dfData$Coef, dfData$Coef.adj1), ]
str(dfData)

# ## setup the model
# library(lme4)
# fit.lme1 = glmer.nb(values ~ 1  + time + (1 | Coef) + (1 | Coef.adj1) + (0 + time | ind), data=dfData)
# summary(fit.lme1)
# fit.lme2 = glmer.nb(values ~ 1 + time + (1 | Coef) + (1 | Coef.adj1), data=dfData)
# summary(fit.lme2)
# fit.lme3 = glmer.nb(values ~ 1 + (1 | Coef) + (1 | Coef.adj1) + (0 + time | ind), data=dfData)
# summary(fit.lme3)
# 
# anova(fit.lme1, fit.lme2, fit.lme3)
# 
# ran = ranef(fit.lme1, condVar=F)
# 
# plot(log(fitted(fit.lme1)), resid(fit.lme1), pch=20, cex=0.7)
# lines(lowess(log(fitted(fit.lme1)), resid(fit.lme1)), col=2)
# plot(log(fitted(fit.lme1)), resid(fit.lme1, type='pearson'), pch=20, cex=0.7)

## setup the stan model
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

stanDso = rstan::stan_model(file='nbinomResp2RandomIntercepts1RandomSlope.stan')

# library(MASS)
# s = max(1, fitdistr(dfData$values, 'negative binomial')$estimate['size'])
## calculate hyperparameters for variance of coefficients
l = gammaShRaFromModeSD(sd(log(dfData$values+0.5)), 2*sd(log(dfData$values+0.5)))
# # ## set initial values
# ran = ranef(fit.lme1)
# r1 = ran$Coef
# r2 = ran$Coef.adj1
# r3 = ran$Coef.adj2
# 
r1 = rep(0, times=nlevels(dfData$Coef))
r2 = rep(0, times=nlevels(dfData$Coef.adj1))
r3 = rep(0, times=nlevels(dfData$ind))
initf = function(chain_id = 1) {
  list(sigmaRan1 = 1, sigmaRan2=1, sigmaRanSlope1=1, rGroupsJitter1=r1, rGroupsJitter2=r2,
       rGroupsSlope3=r3, iSize=8, intercept=5, slope=0)
}

### try a t model without mixture
lStanData = list(Ntotal=nrow(dfData), Nclusters1=nlevels(dfData$Coef),
                 Nclusters2=nlevels(dfData$Coef.adj1),
                 Nclusters3=nlevels(dfData$ind),
                 NgroupMap1=as.numeric(dfData$Coef),
                 NgroupMap2=as.numeric(dfData$Coef.adj1),
                 NgroupMap3=as.numeric(dfData$ind),
                 X = dfData$time,
                 y=dfData$values, 
                 gammaShape=l$shape, gammaRate=l$rate,
                 intercept_mean = mean(log(dfData$values+0.5)), intercept_sd= sd(log(dfData$values+0.5))*1,
                 slope_mean = 0, slope_sd=1)

fit.stan = sampling(stanDso, data=lStanData, iter=1000, chains=6, 
                    pars=c('sigmaRan1', 'sigmaRan2', 'sigmaRanSlope1', 'intercept', 'slope',
                           'iSize',  
                           'rGroupsJitter1', 
                           'rGroupsJitter2',
                           'rGroupsSlope3'),
                    cores=6, init=initf)#, control=list(adapt_delta=0.99, max_treedepth = 15))
save(fit.stan, file='temp/fit.stan.nb_17july.rds')

print(fit.stan, c('intercept', 'slope', 'sigmaRan1', 'sigmaRan2', 'sigmaRanSlope1', 'iSize'), digits=3)
traceplot(fit.stan, c('intercept', 'slope'))
traceplot(fit.stan, 'sigmaRan1')
traceplot(fit.stan, 'sigmaRanSlope1')
print(fit.stan, 'rGroupsJitter1')

## get the coefficient of interest - Modules in our case from the random coefficients section
mCoef = extract(fit.stan)$rGroupsJitter1
dim(mCoef)
# ## get the intercept at population level
iIntercept = as.numeric(extract(fit.stan)$intercept)
## add the intercept to each random effect variable, to get the full coefficient
mCoef = sweep(mCoef, 1, iIntercept, '+')

## function to calculate statistics for differences between coefficients
getDifference = function(ivData, ivBaseline){
  stopifnot(length(ivData) == length(ivBaseline))
  # get the difference vector
  d = ivData - ivBaseline
  # get the z value
  z = mean(d)/sd(d)
  # get 2 sided p-value
  p = pnorm(-abs(mean(d)/sd(d)))*2
  return(list(z=z, p=p))
}

## split the data into the comparisons required
d = data.frame(cols=1:ncol(mCoef), mods=levels(dfData$Coef))
## split this factor into sub factors
f = strsplit(as.character(d$mods), ':')
d = cbind(d, do.call(rbind, f))
head(d)
colnames(d) = c(colnames(d)[1:2], c('fBatch', 'ind'))
d$split = factor(d$ind)

levels(d$fBatch)
## repeat this for each comparison

## get a p-value for each comparison
l = tapply(d$cols, d$split, FUN = function(x, base='SIO', deflection='ILC3') {
  c = x
  names(c) = as.character(d$fBatch[c])
  dif = getDifference(ivData = mCoef[,c[deflection]], ivBaseline = mCoef[,c[base]])
  r = data.frame(ind= as.character(d$ind[c[base]]), coef.base=mean(mCoef[,c[base]]), 
                 coef.deflection=mean(mCoef[,c[deflection]]), zscore=dif$z, pvalue=dif$p)
  r$difference = r$coef.deflection - r$coef.base
  #return(format(r, digi=3))
  return(r)
})

dfResults = do.call(rbind, l)
dfResults$adj.P.Val = p.adjust(dfResults$pvalue, method='BH')

### plot the results
dfResults$logFC = dfResults$difference
dfResults$P.Value = dfResults$pvalue
head(rownames(dfResults))
library(org.Mm.eg.db)
## remove X from annotation names
dfResults$ind = gsub('X', '', as.character(dfResults$ind))

df = AnnotationDbi::select(org.Mm.eg.db, keys = as.character(dfResults$ind), columns = 'SYMBOL', keytype = 'ENTREZID')
i = match(dfResults$ind, df$ENTREZID)
df = df[i,]
dfResults$SYMBOL = df$SYMBOL
identical(dfResults$ind, df$ENTREZID)
## produce the plots 
f_plotVolcano(dfResults, 'ILC2 vs SIO', fc.lim=c(-2.5, 2.5))

m = tapply(dfData$values, dfData$ind, mean)
i = match(rownames(dfResults), names(m))
m = m[i]
identical(names(m), rownames(dfResults))
plotMeanFC(log(m), dfResults, 0.01, 'TNFa vs Control')
table(dfResults$pvalue < 0.05)
## save the results 
write.csv(dfResults, file='results/DEAnalysisTNFaVsControl.xls')



######### do a comparison with deseq2
dfDesign = data.frame(Treatment = dfSample$group1, fAdjust1 = dfSample$group3, fAdjust2 = as.numeric(dfSample$group2),
                      row.names=colnames(mData))

oDseq = DESeqDataSetFromMatrix(mData, dfDesign, design = ~ Treatment + fAdjust1 )
oDseq = DESeq(oDseq)
plotDispEsts(oDseq)
oRes = results(oDseq, contrast=c('Treatment', 'ILC3', 'SIO'))
temp = as.data.frame(oRes)
i = match((dfResults$ind), rownames(temp))
temp = temp[i,]
identical((dfResults$ind), rownames(temp))
plot(dfResults$logFC, temp$log2FoldChange, pch=20)
table(oRes$padj < 0.01)
