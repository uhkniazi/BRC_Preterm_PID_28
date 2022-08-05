# File: 02_VarSelection_Metabolomics.R
# Auth: umar.niazi@kcl.ac.uk
# Date: 3/8/2022
# Desc: variable selection for updated metabolomics data

#### data loading and formatting
source('header.R')

# load the data
dfMeta = read.csv(file.choose(), header=T, stringsAsFactors = T)
dfData = read.csv(file.choose(), header=T)
rownames(dfMeta) = 1:nrow(dfMeta)
rownames(dfData) = 1:nrow(dfMeta)
table(is.na(dfMeta))
table(is.na(dfData))
str(dfMeta)
# one sample is missing a metadata info, remove that
which(dfMeta$BMI_group == '')
dfMeta = dfMeta[-52,]
dfData = dfData[-52,]
identical(rownames(dfMeta), rownames(dfData))
dfMeta = droplevels.data.frame(dfMeta)
dfMeta$outcome_numeric = ifelse(dfMeta$outcome == 'TERM', 0, 1)
colnames(dfMeta)
dfMeta = dfMeta[,c(6, 3, 5)]
str(dfMeta)

## some acrobatics to check data matrix
mCounts = as.matrix(dfData)
dim(mCounts)
mCounts = t(mCounts)

## variables with most 0s
ivProb = apply(mCounts, 1, function(inData) {
  # inData[inData < 1] = 0
  inData = as.logical(inData)
  lData = list('success'=sum(inData), fail=sum(!inData))
  return(mean(rbeta(1000, lData$success + 0.5, lData$fail + 0.5)))
})
names(ivProb) = rownames(mCounts)
i = which(ivProb >= 0.3)
length(i)
ivProb[-i]
dim(mCounts)
mCounts = mCounts[i,]

## transformation
## https://homepages.inf.ed.ac.uk/rbf/HIPR2/pixexp.htm
fTransform = function(x, c=1, r = 0.5){
  return(c * (x^r))
}

m = apply(mCounts, 2, fTransform, c=1, r=0.4)
dim(m); dim(mCounts)
mCounts = m
plot(density(mCounts))

### reduce the number of variables with various steps
## remove variables with 0 sd i.e. not changing 
dim(mCounts)
s = apply(mCounts, 1, sd)
summary(s)
s = which(s == 0)
length(s)

lData.train = list(data=t(mCounts), covariates=dfMeta)

# additional variable to map the sample numbers as 
# we may drop some after matching and overlap
lData.train$covariates$index = rownames(lData.train$covariates)
lData.train$covariates$fGroups = factor(lData.train$covariates$outcome_numeric)
## propensity score calculation to match samples
colnames(lData.train$covariates)
# set up data for binomial model
df = data.frame(lData.train$covariates[, c(5, 2, 3)])
str(df)
f = glm(fGroups ~ Age_group+BMI_group+Age_group:BMI_group , data=df, family=binomial(link='logit'))
summary(f)
## some plots see here for more details
## https://github.com/uhkniazi/BRC_Allergy_Oliver_PID_18/blob/5c5b9ee10f0b97a6bd8479f67140653df2c4e543/04_responseVsDiversity.R#L201
## create the plots for regression with each predictor (not input) fixed at its average
## see Data Analysis ... Regression & Multilevel M [Gelman] for jitter.binary function
jitter.binary = function(a, jitt=.05){
  ifelse (a==0, runif (length(a), 0, jitt), runif (length(a), 1-jitt, 1))
}

# propensity score
ivPropensity = fitted(f)
plot(density(ivPropensity))
lData.train$propensity = ivPropensity
fG = lData.train$covariates$fGroups
hist2(ivPropensity[fG=='1'], 
      ivPropensity[fG=='0'], legends = c('pr', 'n'))

colnames(lData.train$covariates)
z = lData.train$covariates$outcome_numeric
df = data.frame(fGroups = z, a=lData.train$covariates$Age_group, b=lData.train$covariates$BMI_group)
str(df)
f = glm(fGroups ~ a + b + a:b, data = df,
        family=binomial(link='logit'))
summary(f)

library(Matching)
colnames(df)
# making sure no duplications as due to smaller sample sizes
# duplications create over training and correlations between variables
#set.seed(123)
matches = Match(Tr=df$fGroups, X = as.matrix(fitted(f)), estimand = 'ATT', replace = T, ties = T)
summary(matches)
iIndex = unique(c(matches$index.treated, matches$index.control))
table(duplicated(iIndex))
# save data before unmatched data is dropped
lData.train.full = lData.train
f = rep('dropped', times=length(df$fGroups))
f[iIndex] = 'matched' 
lData.train.full$covariates$matched = factor(f)
lData.train$data = lData.train$data[iIndex,]
lData.train$covariates = lData.train$covariates[iIndex,]

## calculate propensity score after matching
colnames(lData.train$covariates)
df = data.frame(lData.train$covariates[, c(5, 2, 3)])
str(df)
f = glm(fGroups ~ Age_group+BMI_group+Age_group:BMI_group , data=df, family=binomial(link='logit'))
summary(f)
ivPropensity = fitted(f)
fG = lData.train$covariates$fGroups
hist2(ivPropensity[fG=='1'], 
      ivPropensity[fG=='0'], legends = c('pr', 'n'), main='After matching', xlab='Propensity')

# ## structure of the metadata after matching
# dfMeta = lData.train$covariates
# dfMeta = droplevels.data.frame(dfMeta)
# xtabs( ~ dfMeta$outcome_numeric + dfMeta$Age_group)
# xtabs( ~ dfMeta$outcome_numeric + dfMeta$BMI_group)
# f = dfMeta$BMI_group:dfMeta$Age_group
# o = dfMeta$outcome_numeric
# xtabs( ~ o + f)
# table(dfMeta$outcome_numeric)
# t = as.matrix(xtabs( ~ o + f))
# t = t/sum(rowSums(t))
# # estimate P(outcome | age, bmi)
# p_age.bmi = colSums(t)
# p_outcome_given_age.bmi = round(sweep(t, 2, p_age.bmi, FUN = '/'), 3)
# 
# b = barplot(p_outcome_given_age.bmi, beside = F, xaxt='n')
# axis(1, at = b, labels = colnames(p_outcome_given_age.bmi), las=2, cex.axis=0.7,
#      tick=F)
# legend('top', legend = c('0', '1'), fill=c('black', 'grey'), xpd=T,
#        horiz=T, inset=c(0, -0.2))
# df = data.frame(d=lData.train$data[,1], y=lData.train$covariates$fGroups, c=lData.train$covariates$Csite, p=ivPropensity)
# 
# ## try some model fits with covariates
# fit.1 = glm(y ~ d, data=df, family='binomial')
# fit.2 = glm(y ~ d + c, data=df, family='binomial')
# fit.3 = glm(y ~ d + p, data=df, family='binomial')
# 
# summary(fit.1); summary(fit.2); summary(fit.3)
# 
# ## regress out the predictor variable
# fit.reg = lm(d ~ p, data=df)
# p.resid = resid(fit.reg)
# fit.4 = glm(y ~ p.resid, data=df, family='binomial')
# summary(fit.1); summary(fit.2); summary(fit.3); summary(fit.4)
# df = data.frame(c=lData.train$covariates$CGSamDel,
#                 f=lData.train$covariates$fGroups)
# hist2(df$c[df$f=='pre'], df$c[df$f=='norm'], main = 'CGSamDel', 
#       legends = c('pre', 'norm'), legend.pos = 'topleft' )

## select variables showing average difference
p.vals = lapply(1:ncol(lData.train$data), function(x){
  df = data.frame(y=lData.train$data[,x], d=lData.train$covariates$fGroups, a=lData.train$covariates$Age_group,
                  b=lData.train$covariates$BMI_group)
  f = lm(y ~ d + a + b, data=df)
  s = summary(f)$coefficients
  return(s['d1', 4])
})

names(p.vals) = colnames(lData.train$data)
dfPvals = do.call(rbind, p.vals)
dfPvals = cbind(dfPvals, p.adjust(dfPvals[,1], method = 'BH'))
colnames(dfPvals) = c('pvalue', 'p.adj')
dfPvals = data.frame(dfPvals)
f = which(dfPvals$pvalue < 0.1)
length(f)
cvTopVariables.lm = rownames(dfPvals)#[f]

if(!require(downloader) || !require(methods)) stop('Library downloader and methods required')

url = 'https://raw.githubusercontent.com/uhkniazi/CCrossValidation/experimental/CCrossValidation.R'
download(url, 'CCrossValidation.R')

# load the required packages
source('CCrossValidation.R')
# delete the file after source
unlink('CCrossValidation.R')

########################## perform a random forest step
dfData = data.frame(lData.train$data[,cvTopVariables.lm])
fGroups = lData.train$covariates$fGroups

oVar.r = CVariableSelection.RandomForest(dfData, fGroups, boot.num = 100)

plot.var.selection(oVar.r)

# ########## adjust the data first
# m = apply(lData.train$data[,cvTopVariables.lm], 2, function(x){
#   return(resid(lm(x ~ ivPropensity)))
# })
# 
# lData.train$data_adjusted = m
# dfData = data.frame(m)
# fGroups = lData.train$covariates$fGroups
# 
# oVar.rAdj = CVariableSelection.RandomForest(dfData, fGroups, boot.num = 100)
# 
# plot.var.selection(oVar.rAdj)

## there does seem to be a difference in the top variables after adjustment
######################## Stan section for binomial regression approach
dfData = data.frame((lData.train$data[,cvTopVariables.lm]))
dim(dfData)
dfData$fGroups = lData.train$covariates$fGroups
table(dfData$fGroups)
rm(fGroups)
levels(dfData$fGroups)
lData = list(resp=ifelse(dfData$fGroups == '1', 1, 0), mModMatrix=model.matrix(fGroups ~ 1 + ., data=dfData))

library(rethinking)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
stanDso = rstan::stan_model(file='binomialRegressionSharedCoeffVariance.stan')

lStanData = list(Ntotal=length(lData$resp), Ncol=ncol(lData$mModMatrix), X=lData$mModMatrix,
                 y=lData$resp)

## give initial values
# initf = function(chain_id = 1) {
#   list(betas=rep(0, times=ncol(lStanData$X)), tau=0.5)
# }


fit.stan = sampling(stanDso, data=lStanData, iter=2000, chains=4, pars=c('tau', 'betas2', 'log_lik'), cores=4,# init=initf,
                    control=list(adapt_delta=0.99, max_treedepth = 13))

# save(fit.stan, file='temp/fit.stan.binom_preterm_met.rds')

print(fit.stan, c('betas2', 'tau'))
print(fit.stan, 'tau')
traceplot(fit.stan, 'tau')

## get the coefficient of interest
mCoef = extract(fit.stan)$betas2
dim(mCoef)
## get the intercept
iIntercept = mCoef[,1]
mCoef = mCoef[,-1]
colnames(mCoef) = colnames(lData$mModMatrix)[2:ncol(lData$mModMatrix)]

## coeftab object 
ct.1 = coeftab(fit.stan)
rn = rownames(ct.1@coefs)
i = grep('betas', rn)
rownames(ct.1@coefs)[i[-1]] = colnames(mCoef)
rownames(ct.1@se)[i[-1]] = colnames(mCoef)
plot(ct.1, pars=colnames(mCoef))

m = abs(colMeans(mCoef))
m = sort(m, decreasing = T)

l2 = barplot(m[1:20], 
             las=2, xaxt='n', col='grey', main='Top Variables - binomial')
axis(1, at = l2, labels = names(m)[1:20], tick = F, las=2, cex.axis=0.9 )

### predictive performance of the model
## binomial prediction
mypred = function(theta, data){
  betas = theta # vector of betas i.e. regression coefficients for population
  ## data
  mModMatrix = data$mModMatrix
  # calculate fitted value
  iFitted = mModMatrix %*% betas
  # using logit link so use inverse logit
  #iFitted = plogis(iFitted)
  return(iFitted)
}


mCoef = extract(fit.stan)$betas2
dim(mCoef)
colnames(mCoef) = c('Intercept', colnames(lData$mModMatrix)[2:ncol(lData$mModMatrix)])
library(lattice)
## get the predicted values
## create model matrix
X = as.matrix(cbind(rep(1, times=nrow(dfData)), dfData[,colnames(mCoef)[-1]]))
colnames(X) = colnames(mCoef)
head(X)
ivPredict = plogis(mypred(colMeans(mCoef), list(mModMatrix=X))[,1])
xyplot(ivPredict ~ fGroups, xlab='Actual Group', main= 'Binomial',
       ylab='Predicted Probability of Pre-term',
       data=dfData)


# ## find correlated variables
dfData.bk = dfData
dim(dfData)
length(cvTopVariables.lm)
mData = as.matrix(dfData[,-27])
length(as.vector(mData))
mCor = cor(mData, use="na.or.complete")
library(caret)
image(mCor)
### find the columns that are correlated and should be removed
n = findCorrelation((mCor), cutoff = 0.8, names=T)
# data.frame(n)
# # sapply(n, function(x) {
# #   (abs(mCor[,x]) >= 0.7)
# # })
# # 
# n = sapply(n, function(x) {
#   rownames(mCor)[(abs(mCor[,x]) >= 0.8)]
# })

#cvDrop.colinear = names(n[(sapply(n, length) > 2)])
cvDrop.colinear = n

cvTopVariables.rf = rownames(CVariableSelection.RandomForest.getVariables(oVar.r))[1:20]
cvTopVariables.bin = names(m)[1:20]
table(cvTopVariables.bin %in% cvTopVariables.rf)
cvTopVariables = unique(c(cvTopVariables.rf, cvTopVariables.bin))
length(cvTopVariables)
cvTopVariables = cvTopVariables[!(cvTopVariables %in% cvDrop.colinear)]
## subset selection
dfData = data.frame(lData.train$data[,cvTopVariables])
dim(dfData)
fGroups = lData.train$covariates$fGroups
oVar.sub = CVariableSelection.ReduceModel(dfData, fGroups, boot.num = 100)
plot.var.selection(oVar.sub)
table(fGroups)
log(25)

cvVar = CVariableSelection.ReduceModel.getMinModel(oVar.sub, size = 3)

## there is a small bug as the 2 proteins 590 and 59 have similar prefix
## that is why this additional protein is added to the results
## it is due to this line 
## https://github.com/uhkniazi/CCrossValidation/blob/55e7ab5ee27403a190fe7daf507bd57731493ae4/CCrossValidation.R#L931
#cvVar = cvVar[1:3]
# cvVar = cvVar[1:6]
table(cvVar %in% cvTopVariables.bin)
table(cvVar %in% cvTopVariables.rf)

## prepare input data
dfData = data.frame(lData.train$data[,cvVar])
fGroups = lData.train$covariates$fGroups
levels(fGroups)
table(fGroups)
# cross validation
oCV.lda = CCrossValidation.LDA(dfData[,cvVar], dfData[,cvVar], fGroups, fGroups, level.predict = '1',
                               boot.num = 100, k.fold = 10) 

plot.cv.performance(oCV.lda)

## use the full data as test data
dfData.test = data.frame(lData.train.full$data[,cvVar])
fGroups.test = lData.train.full$covariates$fGroups
dim(dfData.test)
table(fGroups.test)

oCV.lda = CCrossValidation.LDA(test.dat = dfData.test,
                               train.dat = dfData,
                               test.groups = fGroups.test,
                               train.groups = fGroups,
                               level.predict = '1',
                               boot.num = 100, k.fold = 10) 

plot.cv.performance(oCV.lda)

## names of these proteins
cvVar.names = cvVar

########################################################################
## refit the binomial model and make some figures
########################################################################
dfData = cbind(lData.train$data[,cvVar])#, lData.train$covariates[,c(4, 5)])
head(dfData)
dim(dfData)
dfData = data.frame((dfData))

dfData$fGroups = lData.train$covariates$fGroups
table(dfData$fGroups)
rm(fGroups)
levels(dfData$fGroups)
head(dfData)
lData = list(resp=ifelse(dfData$fGroups == '1', 1, 0), mModMatrix=model.matrix(fGroups ~ 1 + ., data=dfData))

stanDso = rstan::stan_model(file='binomialRegressionGuessMixture.stan')

lStanData = list(Ntotal=length(lData$resp), Ncol=ncol(lData$mModMatrix), X=lData$mModMatrix,
                 y=lData$resp)

## give initial values
# initf = function(chain_id = 1) {
#   list(betas=rep(0, times=ncol(lStanData$X)), tau=0.5)
# }


fit.stan = sampling(stanDso, data=lStanData, iter=2000, chains=4, pars=c('betas', 'log_lik', 'guess'), cores=4,# init=initf,
                    control=list(adapt_delta=0.99, max_treedepth = 13))

# save(fit.stan, file='temp/fit.stan.binom_preterm_met.rds')

print(fit.stan, c('betas', 'guess'))
# print(fit.stan, 'tau')
# traceplot(fit.stan, 'tau')

## get the coefficient of interest
mCoef = extract(fit.stan)$betas
dim(mCoef)
## get the intercept
iIntercept = mCoef[,1]
mCoef = mCoef[,-1]
colnames(mCoef) = colnames(lData$mModMatrix)[2:ncol(lData$mModMatrix)]

## correlation of coefficients
mCor = cor(mCoef)
dim(mCor)
pairs(mCoef, pch=20, cex=0.6)
heatmap(abs(mCor), Rowv = NA, Colv = NA, symm = T, scale = 'none', cexRow = 1, cexCol = 1)

## coeftab object 
ct.1 = coeftab(fit.stan)
rn = rownames(ct.1@coefs)
i = grep('betas', rn)
rownames(ct.1@coefs)[i[-1]] = colnames(mCoef)
rownames(ct.1@se)[i[-1]] = colnames(mCoef)
plot(ct.1, pars=colnames(mCoef))

### predictive performance of the model
mCoef = extract(fit.stan)$betas
dim(mCoef)
colnames(mCoef) = c('Intercept', colnames(lData$mModMatrix)[2:ncol(lData$mModMatrix)])
## get the predicted values
## create model matrix
head(dfData)
X = as.matrix(cbind(rep(1, times=nrow(dfData)), dfData[,-(ncol(dfData))]))
colnames(X) = colnames(mCoef)
head(X)
ivPredict = plogis(mypred(colMeans(mCoef), list(mModMatrix=X))[,1])
xyplot(ivPredict ~ fGroups, xlab='Actual Group', main= 'Binomial 3 Variables',
       ylab='Predicted Probability of Pre-term',
       data=dfData)

f1 = fit.stan # 3 var and without guess
f2 = fit.stan # 3 var with guess
## go back to refit with less number of variables
plot(compare(f1, f2))
compare(f1, f2)
# model without guess parameter is slightly better
## data unmatched data is also being used here
df = data.frame((lData.train.full$data[,colnames(mCoef)[-1]]))
df$fGroups = lData.train.full$covariates$fGroups:lData.train.full$covariates$matched
head(df)
X = as.matrix(cbind(rep(1, times=nrow(df)), df[,-(ncol(df))]))
colnames(X) = colnames(mCoef)
head(X)
ivPredict = plogis(mypred(colMeans(mCoef), list(mModMatrix=X))[,1])
xyplot(ivPredict ~ fGroups, xlab='Actual Group', main= 'Binomial 3 Variables with Dropped',
       ylab='Predicted Probability of Pre-term',
       data=df)
fPredict = ifelse(ivPredict > 0.5, 'P', 'N')
fOriginal = lData.train.full$covariates$fGroups
table(fPredict, fOriginal)

fPredict = ifelse(ivPredict > 0.3, 'P', 'N')
table(fPredict, fOriginal)

# ## which samples from the control group are not predicted correctly
# iOutliers = which(ivPredict > 0.6 & df$fGroups == '0')
# f = as.numeric(df$fGroups) - 1
# f[iOutliers] = 'outlier'
# f = factor(f, levels = c('0', 'outlier', '1'))
# levels(f)
# 
# xyplot(ivPredict ~ f, xlab='Actual Group', main= 'including unmatched outliers',
#        ylab='Predicted Probability of Pre-term',
#        data=df)
# lData.train.full$outlier_index = iOutliers

### figures
fGroups.jitt = jitter.binary(as.numeric(dfData$fGroups)-1)

plot(dfData$Threonine, fGroups.jitt, pch=20, xlab='Threonine', ylab='Probability of preterm',
     main='Prediction of Preterm')
x = seq(min(dfData$Threonine), max(dfData$Threonine), length.out = 100)
colnames(dfData)
m = cbind(1, matrix(colMeans(dfData[,-4]), nrow = length(x), ncol = 3, byrow = T))
m[,2] = x
c = colMeans(extract(fit.stan)$betas)
lines(x, plogis(m %*% c), col='black')
# m = cbind(1, x, min(dfData$TRANCE), min(dfData$CGAsampling), min(dfData$CGSamDel))
# lines(x, plogis(m %*% c), col='red')
# m = cbind(1, x, max(dfData$TRANCE), max(dfData$CGAsampling), max(dfData$CGSamDel))
# lines(x, plogis(m %*% c), col='green')
# legend('left', legend = c('Min', 'Average', 'Max'), fill=c('red', 'black', 'green'))

plot(dfData$Pyruvate, fGroups.jitt, pch=20, xlab='Pyruvate', ylab='Probability of preterm',
     main='Prediction of Preterm')
x = seq(min(dfData$Pyruvate), max(dfData$Pyruvate), length.out = 100)
colnames(dfData)
m = cbind(1, matrix(colMeans(dfData[,-4]), nrow = length(x), ncol = 3, byrow = T))
m[,3] = x
c = colMeans(extract(fit.stan)$betas)
lines(x, plogis(m %*% c), col='black')

plot(dfData$Taurine, fGroups.jitt, pch=20, xlab='Taurine', ylab='Probability of preterm',
     main='Prediction of Preterm')
x = seq(min(dfData$Taurine), max(dfData$Taurine), length.out = 100)
colnames(dfData)
m = cbind(1, matrix(colMeans(dfData[,-4]), nrow = length(x), ncol = 3, byrow = T))
m[,4] = x
c = colMeans(extract(fit.stan)$betas)
lines(x, plogis(m %*% c), col='black')

m = apply(dfData[,-4], 2, min)
m = cbind(1, matrix(m, nrow=length(x), ncol=3, byrow = T))
m[,4] = x
c = colMeans(extract(fit.stan)$betas)
lines(x, plogis(m %*% c), col='red')

m = apply(dfData[,-4], 2, max)
m = cbind(1, matrix(m, nrow=length(x), ncol=3, byrow = T))
m[,4] = x
c = colMeans(extract(fit.stan)$betas)
lines(x, plogis(m %*% c), col='green')

# #########################################################################
# ### repeat model fitting on unadjusted data with adjustment covariate
# #########################################################################
# dfData = data.frame(lData.train$data[,cvVar])
# dim(dfData)
# ## some acrobatics to rename variables
# df = dfKey
# df = df[df$key %in% cvVar,]
# identical(df$key, cvVar)
# i = match(cvVar, df$key)
# identical(cvVar, df$key[i])
# df = df[i,]
# identical(colnames(dfData), df$key)
# colnames(dfData) = df$original_names
# dfData$site = dfMeta$Csite
# dfData$fGroups = lData.train$covariates$fGroups
# table(dfData$fGroups)
# rm(fGroups)
# levels(dfData$fGroups)
# lData = list(resp=ifelse(dfData$fGroups == 'pre', 1, 0), mModMatrix=model.matrix(fGroups ~ 1 + ., data=dfData))
# 
# # stanDso = rstan::stan_model(file='binomialRegressionSharedCoeffVariance.stan')
# 
# lStanData = list(Ntotal=length(lData$resp), Ncol=ncol(lData$mModMatrix), X=lData$mModMatrix,
#                  y=lData$resp)
# 
# ## give initial values
# # initf = function(chain_id = 1) {
# #   list(betas=rep(0, times=ncol(lStanData$X)), tau=0.5)
# # }
# 
# 
# fit.stan.2 = sampling(stanDso, data=lStanData, iter=2000, chains=4, pars=c('tau', 'betas2', 'log_lik'), cores=4,# init=initf,
#                     control=list(adapt_delta=0.99, max_treedepth = 13))
# 
# # save(fit.stan, file='temp/fit.stan.binom_preterm_met.rds')
# 
# print(fit.stan.2, c('betas2', 'tau'))
# # print(fit.stan, 'tau')
# # traceplot(fit.stan, 'tau')
# 
# ## get the coefficient of interest
# mCoef = extract(fit.stan.2)$betas2
# dim(mCoef)
# ## get the intercept
# iIntercept = mCoef[,1]
# mCoef = mCoef[,-1]
# colnames(mCoef) = colnames(lData$mModMatrix)[2:ncol(lData$mModMatrix)]
# 
# ## coeftab object 
# p = paste0('betas2[', 2:7, ']')
# plot(coeftab(fit.stan, fit.stan.2), pars=p)
# ct.1 = coeftab(fit.stan.2)
# rn = rownames(ct.1@coefs)
# i = grep('betas', rn)
# rownames(ct.1@coefs)[i[-1]] = colnames(mCoef)
# rownames(ct.1@se)[i[-1]] = colnames(mCoef)
# plot(ct.1, pars=colnames(mCoef))
