# File: 02_EDA_Proteomics.R
# Auth: umar.niazi@kcl.ac.uk
# Date: 24/03/2022
# Desc: exploratory analysis for the proteomics data from stanford dataset

source('header.R')

# load the data
lData = f_LoadObject('dataUpload/stanford_data_rlist.rds')

dfMeta = lData$covariates
mCounts = lData$proteomics$data
dim(mCounts)
mCounts = t(mCounts)
mCounts[1:10, 1:10]
range(mCounts)
# remove NA values
f = apply(mCounts, 2, is.na)
f2 = rowSums(f)
table(f2)
# it appears some rows are NA i.e. proteins not detected? or data absent
f = mCounts[which(f2 > 0), ]
f[,1:10]
dim(mCounts); dim(na.omit(mCounts))
mCounts = na.omit(mCounts)
range(mCounts)
plot(density(mCounts))
## some EDA diagnostic plots on the data matrix
library(downloader)
url = 'https://raw.githubusercontent.com/uhkniazi/CDiagnosticPlots/master/CDiagnosticPlots.R'
download(url, 'CDiagnosticPlots.R')

# load the required packages
source('CDiagnosticPlots.R')
# delete the file after source
unlink('CDiagnosticPlots.R')

colnames(mCounts) = as.character(dfMeta$fGroups)
oDiag.1 = CDiagnosticPlots(mCounts, 'Proteomics')

# the batch variable we wish to colour by, 
# this can be any grouping/clustering in the data capture process
str(dfMeta)
fBatch = factor(dfMeta$Csite)
levels(fBatch)

# choose a different grouping variable
summary(dfMeta$CGAdelivery)
fBatch = cut(dfMeta$CGAdelivery, 4, include.lowest = T)
table(fBatch)

summary(dfMeta$CGAsampling)
fBatch = cut(dfMeta$CGAsampling, 4, include.lowest = T)
table(fBatch)

summary(dfMeta$CGSamDel)
fBatch = cut(dfMeta$CGSamDel, 5, include.lowest = T)
table(fBatch)

boxplot.median.summary(oDiag.1, fBatch, legend.pos = 'topright', axis.label.cex = 0.5)
plot.mean.summary(oDiag.1, fBatch, axis.label.cex = 0.5)
plot.sigma.summary(oDiag.1, fBatch, axis.label.cex = 0.5)
plot.missing.summary(oDiag.1, fBatch, axis.label.cex = 0.5, cex.main=1)
plot.PCA(oDiag.1, fBatch, cex.main=1)#, csLabels = as.character(dfMeta$fGroups))
plot.dendogram(oDiag.1, fBatch, labels_cex = 0.8, cex.main=0.7)
## change parameters 
l = CDiagnosticPlotsGetParameters(oDiag.1)
l$PCA.jitter = F
l$HC.jitter = F
oDiag.1 = CDiagnosticPlotsSetParameters(oDiag.1, l)
plot.PCA(oDiag.1, fBatch, legend.pos = 'bottom')
plot.dendogram(oDiag.1, fBatch, labels_cex = 0.7)

## extract the PCA components and model the variation
######## modelling of PCA components to assign sources of variance to covariates in the design
par(p.old)
plot(oDiag.1@lData$PCA$sdev)
# use the first 3 principal components
mPC = oDiag.1@lData$PCA$x[,1:3]

## try a linear mixed effect model to account for varince
library(lme4)
# prepare data for input
dfData = data.frame(mPC)
dfData = stack(dfData)
str(dfData)
dfData$values = as.numeric(scale(dfData$values))

library(lattice)
densityplot(~ values, data=dfData)
densityplot(~ values | ind, data=dfData, scales=list(relation='free'))

# add covariates of interest to the data frame
str(lData$covariates)
dfData$fTreatment = lData$covariates$fGroups
dfData$fSite = lData$covariates$Csite
dfData$sampling = cut(lData$covariates$CGAsampling, 5, include.lowest = T)
dfData$samDel = cut(lData$covariates$CGSamDel, 4, include.lowest = T)

densityplot(~ values | ind, groups=fTreatment, data=dfData, auto.key = list(columns=3), scales=list(relation='free'))
densityplot(~ values | ind, groups=fSite, data=dfData, auto.key = list(columns=3), scales=list(relation='free'))
densityplot(~ values | ind, groups=sampling, data=dfData, auto.key = list(columns=3), scales=list(relation='free'))
densityplot(~ values | ind, groups=samDel, data=dfData, auto.key = list(columns=4))
plot(lData$covariates$CGAsampling ~ lData$covariates$CGSamDel)
coplot(lData$covariates$CGAsampling ~ lData$covariates$CGSamDel | lData$covariates$fGroups)
coplot(lData$covariates$CGAsampling ~ lData$covariates$CGSamDel | lData$covariates$Csite)

# format data for modelling, i.e. create coefficients to estimate
str(dfData)
dfData$Coef.1 = factor(dfData$fTreatment:dfData$ind)
dfData$Coef.2 = factor(dfData$fSite:dfData$ind)
dfData$Coef.3 = factor(dfData$sampling:dfData$ind)
dfData$Coef.4 = factor(dfData$samDel:dfData$ind)
str(dfData)

fit.lme1 = lmer(values ~ 1  + (1 | Coef.1), data=dfData)
fit.lme2 = lmer(values ~ 1  + (1 | Coef.1) + (1 | Coef.2), data=dfData)
fit.lme3 = lmer(values ~ 1  + (1 | Coef.1) + (1 | Coef.2) + (1 | Coef.3), data=dfData)
fit.lme4 = lmer(values ~ 1  + (1 | Coef.1) + (1 | Coef.2) + (1 | Coef.3) + (1 | Coef.4), data=dfData)

anova(fit.lme1, fit.lme2, fit.lme3, fit.lme4)

summary(fit.lme1)
summary(fit.lme2)

plot((fitted(fit.lme2)), resid(fit.lme2), pch=20, cex=0.7)
lines(lowess((fitted(fit.lme2)), resid(fit.lme2)), col=2)
hist(dfData$values, prob=T)
lines(density(fitted(fit.lme2)))

## fit model with stan with various model sizes
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(rethinking)

stanDso = rstan::stan_model(file='tResponsePartialPooling.stan')

######## models of various sizes using stan
str(dfData)
m1 = model.matrix(values ~ Coef.1 - 1, data=dfData)
m2 = model.matrix(values ~ Coef.2 - 1, data=dfData)
m3 = model.matrix(values ~ Coef.3 - 1, data=dfData)
m4 = model.matrix(values ~ Coef.4 - 1, data=dfData)
m = cbind(m1, m2, m3, m4)

lStanData = list(Ntotal=nrow(dfData), Ncol=ncol(m), X=m,
                 NscaleBatches=4, NBatchMap=c(rep(1, times=nlevels(dfData$Coef.1)),
                                              rep(2, times=nlevels(dfData$Coef.2)),
                                              rep(3, times=nlevels(dfData$Coef.3)),
                                              rep(4, times=nlevels(dfData$Coef.4))),
                 y=dfData$values)

fit.stan.4 = sampling(stanDso, data=lStanData, iter=1000, chains=2, pars=c('betas', 'populationMean', 'sigmaPop', 'sigmaRan',
                                                                           'nu', 'mu', 'log_lik'),
                      cores=2, control=list(adapt_delta=0.99, max_treedepth = 12))
print(fit.stan.4, c('populationMean', 'sigmaPop', 'sigmaRan', 'nu', 'betas'), digits=3)

traceplot(fit.stan.4, 'populationMean')
traceplot(fit.stan.4, 'sigmaPop')
traceplot(fit.stan.4, 'sigmaRan')

### model without the sampling related covariates
m = cbind(m1, m2)
lStanData = list(Ntotal=nrow(dfData), Ncol=ncol(m), X=m,
                 NscaleBatches=2, NBatchMap=c(rep(1, times=nlevels(dfData$Coef.1)),
                                              rep(2, times=nlevels(dfData$Coef.2))),
                 y=dfData$values)

fit.stan.2 = sampling(stanDso, data=lStanData, iter=1000, chains=2, pars=c('betas', 'populationMean', 'sigmaPop', 'sigmaRan',
                                                                           'nu', 'mu', 'log_lik'),
                      cores=2, control=list(adapt_delta=0.99, max_treedepth = 12))
print(fit.stan.2, c('populationMean', 'sigmaPop', 'sigmaRan', 'nu', 'betas'), digits=3)

traceplot(fit.stan.2, 'populationMean')
traceplot(fit.stan.2, 'sigmaPop')
traceplot(fit.stan.2, 'sigmaRan')

## some model scores and comparisons
compare(fit.stan.4, fit.stan.2)
compare(fit.stan.4, fit.stan.2, func = LOO)
plot(compare(fit.stan.4, fit.stan.2))

############### new simulated data
###############
### generate some posterior predictive data
## generate random samples from alternative t-distribution parameterization
## see https://grollchristian.wordpress.com/2013/04/30/students-t-location-scale/
rt_ls <- function(n, df, mu, a) rt(n,df)*a + mu
## follow the algorithm in section 14.3 page 363 in Gelman 2013
simulateOne = function(mu, sigma, nu){
  yrep = rt_ls(length(mu), nu, mu,  sigma)
  return(yrep)
}

## sample n values, numerous times
mDraws.sim = matrix(NA, nrow = nrow(dfData), ncol=300)
l = extract(fit.stan.2)
for (i in 1:300){
  p = sample(1:nrow(l$mu), 1)
  mDraws.sim[,i] = simulateOne(l$mu[p,], 
                               l$sigmaPop[p],
                               l$nu[p])
}

dim(mDraws.sim)
plot(density(dfData$values), main='posterior predictive density plots, model 4')
apply(mDraws.sim, 2, function(x) lines(density(x), lwd=0.5, col='lightgrey'))
lines(density(dfData$values))

## plot residuals
plot(dfData$values - colMeans(l$mu) ~ colMeans(l$mu))
lines(lowess(colMeans(l$mu), dfData$values - colMeans(l$mu)))
apply(l$mu[sample(1:nrow(l$mu), 100),], 1, function(x) {
  lines(lowess(x, dfData$values - x), lwd=0.5, col=2)
})

## plot the original PCA and replicated data
plot(dfData$values[dfData$ind == 'PC1'], dfData$values[dfData$ind == 'PC2'], 
     col=c(1,2)[as.numeric(dfData$fTreatment[dfData$ind == 'PC1'])], main='PCA Components - original and simulated',
     xlab='PC1', ylab='PC2')
points(rowMeans(mDraws.sim)[dfData$ind == 'PC1'], rowMeans(mDraws.sim)[dfData$ind == 'PC2'],
       col=c(1,2)[as.numeric(dfData$fTreatment[dfData$ind == 'PC1'])], pch='1')

plot(dfData$values[dfData$ind == 'PC1'], dfData$values[dfData$ind == 'PC2'], 
     col=c(1,2)[as.numeric(dfData$fTreatment[dfData$ind == 'PC1'])], main='PCA Components - original and model 2',
     xlab='PC1', ylab='PC2', xlim=c(-4, 4), ylim=c(-3, 4), pch=15)

apply(mDraws.sim, 2, function(x) {
  points(x[dfData$ind == 'PC1'], x[dfData$ind == 'PC2'],
         col=c(1,2)[as.numeric(dfData$fTreatment[dfData$ind == 'PC1'])], pch=20)
})

## try a different colour
plot(dfData$values[dfData$ind == 'PC1'], dfData$values[dfData$ind == 'PC2'], 
     col=c(1:5)[as.numeric(dfData$fSite[dfData$ind == 'PC1'])], main='PCA Components - original and model 2',
     xlab='PC1', ylab='PC2', xlim=c(-4.5, 4.5), ylim=c(-4, 4), pch=15, cex=1.5)

apply(mDraws.sim, 2, function(x) {
  points(x[dfData$ind == 'PC1'], x[dfData$ind == 'PC2'],
         col=c(1:5)[as.numeric(dfData$fSite[dfData$ind == 'PC1'])], pch=20, cex=0.3)
})

points(dfData$values[dfData$ind == 'PC1'], dfData$values[dfData$ind == 'PC2'], 
       col=c(1:5)[as.numeric(dfData$fSite[dfData$ind == 'PC1'])], pch=15, cex=1.5)

legend('topleft', legend = levels(dfData$fSite), fill=c(1:5))

# colour the data with the regions covariates
plot(dfData$values[dfData$ind == 'PC1'], dfData$values[dfData$ind == 'PC2'], 
     col=c(1:5)[as.numeric(dfData$fSite[dfData$ind == 'PC1'])], main='PCA Components - original and simulated',
     xlab='PC1', ylab='PC2', pch=c(10,4)[as.numeric(dfData$fTreatment)[dfData$ind == 'PC1']])
legend('topleft', legend = c('norm', 'pre'), pch=c(10, 4))
points(rowMeans(mDraws.sim)[dfData$ind == 'PC1'], rowMeans(mDraws.sim)[dfData$ind == 'PC2'],
       col=c(1:5)[as.numeric(dfData$fSite[dfData$ind == 'PC1'])], pch='1')
legend('topright', legend = levels(dfData$fSite), fill=c(1:5))


m = cbind(extract(fit.stan.2)$sigmaRan, extract(fit.stan.2)$sigmaPop) 
dim(m)
m = log(m)
colnames(m) = c('Treatment', 'Region', 'Residual')
pairs(m, pch=20, cex=0.5, col='grey')

df = stack(data.frame(m[,-3]))
histogram(~ values | ind, data=df, xlab='Log SD', scales=list(relation='free'))

## calculate bayesian p-value for this test statistic
getPValue = function(Trep, Tobs){
  left = sum(Trep <= Tobs)/length(Trep)
  right = sum(Trep >= Tobs)/length(Trep)
  return(min(left, right))
}
## define some test quantities to measure the lack of fit
## define a test quantity T(y, theta)
## variance
T1_var = function(Y) return(var(Y))

## min quantity
T1_min = function(Y){
  return(min(Y))
} 

## max quantity
T1_max = function(Y){
  return(max(Y))
} 

## mean quantity
T1_mean = function(Y){
  return(mean(Y))
} 

## mChecks
ivResp = dfData$values
mChecks = matrix(NA, nrow=4, ncol=1)
rownames(mChecks) = c('Variance', 'Max', 'Min', 'Mean')
colnames(mChecks) = c('model 1')

t1 = apply(mDraws.sim, 2, T1_var)
mChecks['Variance', 1] = getPValue(t1, var(ivResp))

## testing for outlier detection i.e. the minimum value show in the histograms earlier
t1 = apply(mDraws.sim, 2, T1_min)
t2 = T1_min(ivResp)
mChecks['Min',1] = getPValue(t1, t2)

## maximum value
t1 = apply(mDraws.sim, 2, T1_max)
t2 = T1_max(ivResp)
mChecks['Max', 1] = getPValue(t1, t2)

## mean value
t1 = apply(mDraws.sim, 2, T1_mean)
t2 = T1_mean(ivResp)
mChecks['Mean', 1] = getPValue(t1, t2)

mChecks
