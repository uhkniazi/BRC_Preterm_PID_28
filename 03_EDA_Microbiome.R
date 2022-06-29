# File: 03_EDA_Microbiome.R
# Auth: umar.niazi@kcl.ac.uk
# Date: 29/6/2022
# Desc: exploratory analysis for the microbiome data

source('header.R')

# load the data
dfMeta = read.csv(file.choose(), header=T)
dfData = read.csv(file.choose(), header=T)
rownames(dfMeta) = 1:nrow(dfMeta)
rownames(dfData) = 1:nrow(dfMeta)

# one sample is missing a metadata info, remove that
which(dfMeta$BMI_group == '')
dfMeta = dfMeta[-52,]
dfData = dfData[-52,]
identical(rownames(dfMeta), rownames(dfData))
dfMeta = droplevels.data.frame(dfMeta)

## structure of the metadata
xtabs( ~ dfMeta$outcome_numeric + dfMeta$Age_group)
xtabs( ~ dfMeta$outcome_numeric + dfMeta$BMI_group)
f = factor(dfMeta$BMI_group):factor(dfMeta$Age_group)
levels(f)
o = dfMeta$outcome_numeric
xtabs( ~ o + f)
table(dfMeta$outcome_numeric)
t = as.matrix(xtabs( ~ o + f))
t = t/sum(rowSums(t))
# estimate P(outcome | age, bmi)
p_age.bmi = colSums(t)
p_outcome_given_age.bmi = round(sweep(t, 2, p_age.bmi, FUN = '/'), 3)

b = barplot(p_outcome_given_age.bmi, beside = F, xaxt='n')
axis(1, at = b, labels = colnames(p_outcome_given_age.bmi), las=2, cex.axis=0.7,
     tick=F)
legend('top', legend = c('0', '1'), fill=c('black', 'grey'), xpd=T,
       horiz=T, inset=c(0, -0.2))
## some acrobatics to check data matrix
mCounts = as.matrix(dfData)
dim(mCounts)
dim(dfMeta)
mCounts = t(mCounts)
mCounts[1:10, 1:10]
range(mCounts)
# remove NA values
f = apply(mCounts, 2, is.na)
dim(f)
f2 = rowSums(f)
table(f2)
# # it appears some rows are NA i.e. proteins not detected? or data absent
# f = mCounts[which(f2 > 0), ]
# f[,1:10]
dim(mCounts); dim(na.omit(mCounts))
# mCounts = na.omit(mCounts)
range(mCounts)
plot(density(mCounts))
## try a log transformation of the data to control spread and sparsness
plot(density(log(mCounts+1)))
table(mCounts == 0)
## most of the data has zeros
summary(log(as.vector(mCounts+1)))
table(log(mCounts+1) <= 1.4) # 3rd quantile
# data has a lot of zeros - work with binary matrix instead?
mCounts.orig = mCounts
mCounts = log(mCounts+1)
mCounts[mCounts > 0] = 1
range(mCounts)
## some EDA diagnostic plots on the data matrix
library(downloader)
url = 'https://raw.githubusercontent.com/uhkniazi/CDiagnosticPlots/experimental/CDiagnosticPlots.R'
download(url, 'CDiagnosticPlots.R')

# load the required packages
source('CDiagnosticPlots.R')
# delete the file after source
unlink('CDiagnosticPlots.R')
#colnames(mCounts) = as.character(dfMeta$outcome_numeric)
oDiag.1 = CDiagnosticPlots(mCounts, 'Microbiome Binary')

# the batch variable we wish to colour by, 
# this can be any grouping/clustering in the data capture process
str(dfMeta)
fBatch = factor(dfMeta$outcome_numeric)
fBatch = factor(dfMeta$Age_group)
fBatch = factor(dfMeta$BMI_group)
fBatch = factor(dfMeta$Age_group):factor(dfMeta$BMI_group)
levels(fBatch)
table(fBatch)
xtabs( ~ dfMeta$outcome_numeric + dfMeta$Age_group)
iPch = c(17, 20)[as.numeric(dfMeta$outcome_numeric)+1]

boxplot.median.summary(oDiag.1, fBatch, legend.pos = 'topright', axis.label.cex = 0.5)
plot.mean.summary(oDiag.1, fBatch, axis.label.cex = 0.5)
plot.sigma.summary(oDiag.1, fBatch, axis.label.cex = 0.5)
plot.missing.summary(oDiag.1, fBatch, axis.label.cex = 0.5, cex.main=1)
## change parameters 
l = CDiagnosticPlotsGetParameters(oDiag.1)
l$PCA.scaleSubjects = F
l$PCA.scaleVariables = F
l$HC.scaleSubjects = F
l$HC.scaleVaribles = F
l$PCA.jitter = F
l$HC.jitter = F
oDiag.1 = CDiagnosticPlotsSetParameters(oDiag.1, l)
plot.PCA(oDiag.1, fBatch, cex.main=1, pch = iPch, pch.cex = 1, legend.pos = 'topright')#, csLabels = as.character(dfMeta$fGroups))
legend('top', legend = c('0', '1'), pch=c(17, 20))
plot.dendogram(oDiag.1, fBatch, labels_cex = 0.8, cex.main=0.7)
#ob = calculateExtremeValues(oDiag.1)
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
str(dfMeta)
dfData$fTreatment = as.factor(dfMeta$outcome_numeric)
dfData$fAge = as.factor(dfMeta$Age_group)
dfData$fBmi = as.factor(dfMeta$BMI_group)

densityplot(~ values | ind, groups=fTreatment, data=dfData, auto.key = list(columns=3), scales=list(relation='free'))
densityplot(~ values | ind, groups=fAge, data=dfData, auto.key = list(columns=3), scales=list(relation='free'))
densityplot(~ values | ind, groups=fBmi, data=dfData, auto.key = list(columns=3), scales=list(relation='free'))
densityplot(~ values | ind+fAge, groups=fTreatment, data=dfData, auto.key = list(columns=4), scales=list(relation='free'))
densityplot(~ values | ind+fBmi, groups=fTreatment, data=dfData, auto.key = list(columns=4), scales=list(relation='free'))
densityplot(~ values | ind+fBmi+fAge, groups=fTreatment, data=dfData, auto.key = list(columns=4), scales=list(relation='free'))

# format data for modelling, i.e. create coefficients to estimate
str(dfData)
dfData$Coef.1 = factor(dfData$fTreatment:dfData$ind)
dfData$Coef.2 = factor(dfData$fAge:dfData$ind)
dfData$Coef.3 = factor(dfData$fBmi:dfData$ind)
str(dfData)

fit.lme1 = lmer(values ~ 1  + (1 | Coef.1), data=dfData)
fit.lme2 = lmer(values ~ 1  + (1 | Coef.1) + (1 | Coef.2), data=dfData)
fit.lme3 = lmer(values ~ 1  + (1 | Coef.1) + (1 | Coef.2) + (1 | Coef.3), data=dfData)

anova(fit.lme1, fit.lme2, fit.lme3)

summary(fit.lme1)
summary(fit.lme2)

# plot((fitted(fit.lme2)), resid(fit.lme2), pch=20, cex=0.7)
# lines(lowess((fitted(fit.lme2)), resid(fit.lme2)), col=2)
# hist(dfData$values, prob=T)
# lines(density(fitted(fit.lme2)))

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
m = cbind(m1, m2, m3)

lStanData = list(Ntotal=nrow(dfData), Ncol=ncol(m), X=m,
                 NscaleBatches=3, NBatchMap=c(rep(1, times=nlevels(dfData$Coef.1)),
                                              rep(2, times=nlevels(dfData$Coef.2)),
                                              rep(3, times=nlevels(dfData$Coef.3))),
                 y=dfData$values)

fit.stan.3 = sampling(stanDso, data=lStanData, iter=1000, chains=2, pars=c('betas', 'populationMean', 'sigmaPop', 'sigmaRan',
                                                                           'nu', 'mu', 'log_lik'),
                      cores=2, control=list(adapt_delta=0.99, max_treedepth = 12))
print(fit.stan.3, c('populationMean', 'sigmaPop', 'sigmaRan', 'nu', 'betas'), digits=3)

traceplot(fit.stan.3, 'populationMean')
traceplot(fit.stan.3, 'sigmaPop')
traceplot(fit.stan.3, 'sigmaRan')

### model without the sampling related covariates
m = m1 # cbind(m1, m2)
lStanData = list(Ntotal=nrow(dfData), Ncol=ncol(m), X=m,
                 NscaleBatches=1, NBatchMap=c(rep(1, times=nlevels(dfData$Coef.1))),
                                             # rep(2, times=nlevels(dfData$Coef.2))),
                 y=dfData$values)

fit.stan.1 = sampling(stanDso, data=lStanData, iter=1000, chains=2, pars=c('betas', 'populationMean', 'sigmaPop', 'sigmaRan',
                                                                           'nu', 'mu', 'log_lik'),
                      cores=2, control=list(adapt_delta=0.99, max_treedepth = 12))
print(fit.stan.1, c('populationMean', 'sigmaPop', 'sigmaRan', 'nu', 'betas'), digits=3)

traceplot(fit.stan.1, 'populationMean')
traceplot(fit.stan.1, 'sigmaPop')
traceplot(fit.stan.1, 'sigmaRan')

## some model scores and comparisons
compare(fit.stan.3, fit.stan.1)
#compare(fit.stan.4, fit.stan.2, func = LOO)
plot(compare(fit.stan.3, fit.stan.1))

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
l = extract(fit.stan.3)
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

plot.PCA(oDiag.1, factor(dfMeta$outcome_numeric))
## plot the original PCA and replicated data
plot(dfData$values[dfData$ind == 'PC1'], dfData$values[dfData$ind == 'PC2'], 
     col=rainbow(2)[as.numeric(dfData$fTreatment[dfData$ind == 'PC1'])], main='PCA Components - original and simulated',
     xlab='PC1', ylab='PC2')
points(rowMeans(mDraws.sim)[dfData$ind == 'PC1'], rowMeans(mDraws.sim)[dfData$ind == 'PC2'],
       col=rainbow(2)[as.numeric(dfData$fTreatment[dfData$ind == 'PC1'])], pch='1')

plot(dfData$values[dfData$ind == 'PC1'], dfData$values[dfData$ind == 'PC2'], 
     col=rainbow(2)[as.numeric(dfData$fTreatment[dfData$ind == 'PC1'])], main='PCA Components - original and model 2',
     xlab='PC1', ylab='PC2', xlim=c(-4, 4), ylim=c(-3, 3), pch=15)

apply(mDraws.sim, 2, function(x) {
  points(x[dfData$ind == 'PC1'], x[dfData$ind == 'PC2'], cex=0.5,
         col=rainbow(2)[as.numeric(dfData$fTreatment[dfData$ind == 'PC1'])], pch=20)
})
points(dfData$values[dfData$ind == 'PC1'], dfData$values[dfData$ind == 'PC2'], 
       col=rainbow(2)[as.numeric(dfData$fTreatment[dfData$ind == 'PC1'])], pch=15, cex=1.5)


## try a different colour
plot(dfData$values[dfData$ind == 'PC1'], dfData$values[dfData$ind == 'PC2'], 
     col=rainbow(3)[as.numeric(dfData$fBmi[dfData$ind == 'PC1'])], main='PCA Components - original and simulated',
     xlab='PC1', ylab='PC2')
points(rowMeans(mDraws.sim)[dfData$ind == 'PC1'], rowMeans(mDraws.sim)[dfData$ind == 'PC2'],
       col=rainbow(3)[as.numeric(dfData$fBmi[dfData$ind == 'PC1'])], pch='1')


plot(dfData$values[dfData$ind == 'PC1'], dfData$values[dfData$ind == 'PC2'], 
     col=rainbow(3)[as.numeric(dfData$fBmi[dfData$ind == 'PC1'])], main='PCA Components - original and model 2',
     xlab='PC1', ylab='PC2', xlim=c(-4, 4), ylim=c(-3, 3), pch=15)

apply(mDraws.sim, 2, function(x) {
  points(x[dfData$ind == 'PC1'], x[dfData$ind == 'PC2'], cex=0.5,
         col=rainbow(3)[as.numeric(dfData$fBmi[dfData$ind == 'PC1'])], pch=20)
})

# points(dfData$values[dfData$ind == 'PC1'], dfData$values[dfData$ind == 'PC2'], 
#        col=c(1:5)[as.numeric(dfData$fSite[dfData$ind == 'PC1'])], pch=15, cex=1.5)
# 
# legend('topleft', legend = levels(dfData$fSite), fill=c(1:5))

m = cbind(extract(fit.stan.3)$sigmaRan, extract(fit.stan.3)$sigmaPop) 
dim(m)
m = log(m)
colnames(m) = c('Preterm', 'Age', 'BMI', 'Residual')
pairs(m, pch=20, cex=0.5, col='grey')

df = stack(data.frame(m[,-4]))
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
