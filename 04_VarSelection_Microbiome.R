# File: 04_VarSelection_Microbiome.R
# Auth: umar.niazi@kcl.ac.uk
# Date: 6/7/2022
# Desc: variable selection for microbial data

#### data loading and formatting
source('header.R')

dfMeta = read.csv(file.choose(), header=T, stringsAsFactors = T, row.names = 1)
dfData = read.csv(file.choose(), header=T, row.names = 1)
identical(rownames(dfMeta), rownames(dfData))
table(is.na(dfMeta))
table(is.na(dfData))
str(dfMeta)
# one sample is missing a metadata info, remove that
which(dfMeta$BMI_group == '')
dfMeta = dfMeta[-52,]
dfData = dfData[-52,]
identical(rownames(dfMeta), rownames(dfData))
dfMeta = droplevels.data.frame(dfMeta)
# dfMeta$outcome_numeric = ifelse(dfMeta$outcome == 'TERM', 0, 1)
colnames(dfMeta)
dfMeta = dfMeta[,c(3, 5, 8)]
str(dfMeta)

## trim and transform the count matrix
mCounts = as.matrix(dfData)
range(mCounts)
dim(mCounts)
dim(dfMeta)
mCounts = t(mCounts)
mCounts.orig = mCounts
mCounts[mCounts > 0] = 1
range(mCounts)

## variables with most 0s
ivProb = apply(mCounts, 1, function(inData) {
  inData[inData < 1] = 0
  inData = as.logical(inData)
  lData = list('success'=sum(inData), fail=sum(!inData))
  return(mean(rbeta(1000, lData$success + 0.5, lData$fail + 0.5)))
})

hist(ivProb, main='Proportion of non-zeros in each variable', xlab='')
summary(ivProb)
quantile(ivProb, 0:10/10)
table(dfMeta$outcome_numeric)
25/(58+25)
# about 30% of the data is from minority class
# drop variables with lower probability of expression
## use a binary matrix for lower expressed variables
i = which(ivProb >= 0.3)
length(i)
dim(mCounts)
mCounts = mCounts[i,]

## use original values for the highly expressed variables
## in order to preserve more information 
i = which(ivProb >= 0.9)
length(i)
i = names(i)
dim(mCounts.orig)
m = mCounts.orig[i,]
head(m); dim(m)
m = log(m+1)
range(m)
dim(mCounts)
range(mCounts)
mCounts[i,] = m[i,]
range(mCounts)

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
# f = glm(fGroups ~ fCluster+Age_group+BMI_group+Age_group:BMI_group , data=df, family=binomial(link='logit'))
# summary(f)
f2 = glm(fGroups ~ Age_group+BMI_group+Age_group:BMI_group , data=df, family=binomial(link='logit'))
summary(f2)
# anova(f, f2)

## some plots see here for more details
## https://github.com/uhkniazi/BRC_Allergy_Oliver_PID_18/blob/5c5b9ee10f0b97a6bd8479f67140653df2c4e543/04_responseVsDiversity.R#L201
## create the plots for regression with each predictor (not input) fixed at its average
## see Data Analysis ... Regression & Multilevel M [Gelman] for jitter.binary function
jitter.binary = function(a, jitt=.05){
  ifelse (a==0, runif (length(a), 0, jitt), runif (length(a), 1-jitt, 1))
}
f = f2
# propensity score
ivPropensity = fitted(f)
plot(density(ivPropensity))
lData.train$propensity = ivPropensity
fG = lData.train$covariates$fGroups
hist2(ivPropensity[fG=='1'], 
      ivPropensity[fG=='0'], legends = c('pr', 'n'), main = 'Before Matching', xlab='Propensity')

## as there is no overlap, use the covariates only in overlapping areas
## use the matching library

# perform matching using propensity score
colnames(lData.train$covariates)
z = lData.train$covariates$outcome_numeric
df = data.frame(fGroups = z, a=lData.train$covariates$Age_group, b=lData.train$covariates$BMI_group) 
                # c=lData.train$covariates$fCluster)
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


str(lData.train$covariates)

## skipping section on reducing variable count on average differences
## as it doesn't work very well with binary response variables
# ## select variables showing average difference
# p.vals = lapply(1:ncol(lData.train$data), function(x){
#   df = data.frame(y=jitter.binary(lData.train$data[,'Prevotella_oris']), d=lData.train$covariates$fGroups, a=lData.train$covariates$Age_group,
#                   b=lData.train$covariates$BMI_group, 
#                   c = lData.train$covariates$fCluster)
#   f = lm(y ~ d + a + b + c, data=df)
#   s = summary(f)$coefficients
#   return(s['d1', 4])
# })
# 
# names(p.vals) = colnames(lData.train$data)
# dfPvals = do.call(rbind, p.vals)
# dfPvals = cbind(dfPvals, p.adjust(dfPvals[,1], method = 'BH'))
# colnames(dfPvals) = c('pvalue', 'p.adj')
# dfPvals = data.frame(dfPvals)
# f = which(dfPvals$pvalue < 0.1)
# length(f)
# 
# cvTopVariables.lm = rownames(dfPvals)[f]
cvTopVariables.lm = colnames(lData.train$data)

if(!require(downloader) || !require(methods)) stop('Library downloader and methods required')

url = 'https://raw.githubusercontent.com/uhkniazi/CCrossValidation/experimental/CCrossValidation.R'
download(url, 'CCrossValidation.R')

# load the required packages
source('CCrossValidation.R')
# delete the file after source
unlink('CCrossValidation.R')

########################## perform a random forest step
# fCluster = lData.train$covariates$fCluster
dfData = data.frame(lData.train$data[,cvTopVariables.lm])
fGroups = lData.train$covariates$fGroups
dim(dfData)
levels(fGroups)
# xtabs(~ fGroups + fCluster)

oVar.r = CVariableSelection.RandomForest(dfData, fGroups, boot.num = 100, big.warn = F)
plot.var.selection(oVar.r)

## skip the sections 
# oVar.r.c1 = CVariableSelection.RandomForest(dfData[fCluster == 'C1',], fGroups[fCluster == 'C1'], boot.num = 100, big.warn = F)
# plot.var.selection(oVar.r.c1)
# 
# oVar.r.c2 = CVariableSelection.RandomForest(dfData[fCluster == 'C2',], fGroups[fCluster == 'C2'], boot.num = 100, big.warn = F)
# plot.var.selection(oVar.r.c2)


### not using the adjustment approach
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
dfData = data.frame(lData.train$data[,cvTopVariables.lm])
dim(dfData)
head(dfData)
dfData$fGroups = lData.train$covariates$fGroups
# dfData$fCluster = as.numeric(lData.train$covariates$fCluster)-1
#fCluster = lData.train$covariates$fCluster
#dfData = dfData[fCluster == 'C2',]
table(dfData$fGroups)
dim(dfData)
rm(fGroups)
rm(fCluster)
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
axis(1, at = l2, labels = names(m)[1:20], tick = F, las=2, cex.axis=0.8 )

plot(ct.1, pars=names(m)[1:20])
cvTopVariables.bin.full = names(m)[1:20]
## not using cluster anymore as coefficient is zero
# cvTopVariables.bin.C1 = names(m)[1:10]
# cvTopVariables.bin.C2 = names(m)[1:10]

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

fit.stan.full = fit.stan
# fit.stan.C1 = fit.stan
# fit.stan.C2 = fit.stan 

# ## find correlated variables
dfData.bk = dfData
dim(dfData)
length(cvTopVariables.lm)
mData = t(mCounts[cvTopVariables.lm,])#as.matrix(dfData[,-230])
dim(mData)
s = apply(mData, 2, sd)
table( s == 0)
#mData = mData[,-(which(s == 0))]
length(as.vector(mData))
mCor = cor(mData, use="na.or.complete")
library(caret)
image(mCor)
### find the columns that are correlated and should be removed
n = findCorrelation((mCor), cutoff = 0.7, names=T)
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
#cvTopVariables.rf = c(cvTopVariables.rf, rownames(CVariableSelection.RandomForest.getVariables(oVar.r.c1))[1:10])
#cvTopVariables.rf = c(cvTopVariables.rf, rownames(CVariableSelection.RandomForest.getVariables(oVar.r.c2))[1:10])
#cvTopVariables.rf = unique(cvTopVariables.rf)
length(cvTopVariables.rf)

# cvTopVariables.bin = names(m)[1:20]
cvTopVariables.bin = cvTopVariables.bin.full# unique(c(cvTopVariables.bin.full, cvTopVariables.bin.C1, cvTopVariables.bin.C2))
length(cvTopVariables.bin)
table(cvTopVariables.bin %in% cvTopVariables.rf)
cvTopVariables = unique(c(cvTopVariables.rf, cvTopVariables.bin))
length(cvTopVariables)
cvTopVariables = cvTopVariables[!(cvTopVariables %in% cvDrop.colinear)]
## subset selection
dfData = data.frame(lData.train$data[,cvTopVariables])
dim(dfData)
fGroups = lData.train$covariates$fGroups
# fCluster = lData.train$covariates$fCluster
oVar.sub = CVariableSelection.ReduceModel(dfData, fGroups, boot.num = 100)
plot.var.selection(oVar.sub)
table(fGroups)
log(25)
cvVar = CVariableSelection.ReduceModel.getMinModel(oVar.sub, size = 4)
length(cvVar)
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

# ### repeat on subset clusters
# # dfData = data.frame(lData.train$data[,cvVar])
# # fGroups = lData.train$covariates$fGroups
# # levels(fGroups)
# table(fGroups[fCluster == 'C1'])
# head(dfData)
# # cross validation
# oCV.lda.C1 = CCrossValidation.LDA(dfData[fCluster == 'C1',], dfData[fCluster == 'C1',],
#                                   fGroups[fCluster == 'C1'], fGroups[fCluster == 'C1'], level.predict = '1',
#                                boot.num = 100, k.fold = 2) 
# 
# plot.cv.performance(oCV.lda.C1)
# 
# dfData = data.frame(v1 = lData.train$data[,cvVar.C2])
# fGroups = lData.train$covariates$fGroups
# levels(fGroups)
# table(fGroups[fCluster == 'C2'])
# head(dfData)
# # cross validation
# oCV.lda.C2 = CCrossValidation.LDA(data.frame(v1=dfData[fCluster == 'C2',]), data.frame(v1=dfData[fCluster == 'C2',]),
#                                   fGroups[fCluster == 'C2'], fGroups[fCluster == 'C2'], level.predict = '1',
#                                   boot.num = 100, k.fold = 10) 
# 
# plot.cv.performance(oCV.lda.C2)


## use the full data as test data
dfData = data.frame(lData.train$data[,cvVar])
fGroups = lData.train$covariates$fGroups
levels(fGroups)
table(fGroups)

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

######################################################################
## binomial regression and cross-validation using stan instead of LDA
######################################################################
rm(stanDso)
url = 'https://raw.githubusercontent.com/uhkniazi/CCrossValidation/experimental/bernoulli.stan'
download(url, 'bernoulli.stan')

dfData = data.frame(lData.train$data[,cvVar])
fGroups = lData.train$covariates$fGroups
levels(fGroups)
table(fGroups)

dfData.test = data.frame(lData.train.full$data[,cvVar])
fGroups.test = lData.train.full$covariates$fGroups
dim(dfData.test)
table(fGroups.test)


oCV.s = CCrossValidation.StanBern(train.dat = dfData, 
                                  test.dat = dfData.test, 
                                  test.groups = fGroups.test, 
                                  train.groups = fGroups,
                                  level.predict = '1',
                                  boot.num = 10, k.fold = 10, 
                                  ncores = 2, nchains = 2) 

save(oCV.s, file='temp/oCV.s.rds')

plot.cv.performance(oCV.s)
unlink('bernoulli.stan')
unlink('bernoulli.rds')

########################################################################
## refit the binomial model and make some figures
########################################################################
dfData = cbind(lData.train$data[,cvVar.names])#, lData.train$covariates[,c(4, 5)])
head(dfData)
dim(dfData)
dfData = data.frame((dfData))
dfData = stack(dfData)
head(dfData)
dfData$fGroups = lData.train$covariates$fGroups
# dfData$fCluster = lData.train$covariates$fCluster
str(dfData)
xyplot(values ~ fGroups | ind, data=dfData, auto.key = list(columns=2), scales=list(relation='free'),
       type='p',
       par.strip.text=list(cex=0.7), varwidth=T, main='Gardnerella_vaginalis')

# bwplot(Streptococcus_pseudoporcinus ~ fGroups | fCluster, data=dfData, panel=function(x, y, ...) panel.bwplot(x, y, pch='|',...), type='b',
#        par.strip.text=list(cex=0.7), varwidth=T, main='Streptococcus_pseudoporcinus')

dfData = cbind(lData.train$data[,cvVar.names])#, lData.train$covariates[,c(4, 5)])
head(dfData)
dim(dfData)
dfData = data.frame((dfData))
dfData$fGroups = lData.train$covariates$fGroups
#dfData$fCluster = as.numeric(lData.train$covariates$fCluster)-1
str(dfData)

table(dfData$fGroups)
rm(fGroups)
rm(fCluster)
levels(dfData$fGroups)
head(dfData)
lData = list(resp=ifelse(dfData$fGroups == '1', 1, 0), mModMatrix=model.matrix(fGroups ~ 1 + ., data=dfData))

stanDso = rstan::stan_model(file='binomialRegression.stan')

lStanData = list(Ntotal=length(lData$resp), Ncol=ncol(lData$mModMatrix), X=lData$mModMatrix,
                 y=lData$resp)

## give initial values
# initf = function(chain_id = 1) {
#   list(betas=rep(0, times=ncol(lStanData$X)), tau=0.5)
# }


fit.stan = sampling(stanDso, data=lStanData, iter=2000, chains=4, pars=c('betas', 'log_lik'), cores=4,# init=initf,
                    control=list(adapt_delta=0.99, max_treedepth = 13))

# save(fit.stan, file='temp/fit.stan.binom_preterm_met.rds')

print(fit.stan, c('betas'))#, 'guess'))
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
ncol(dfData)
X = as.matrix(cbind(rep(1, times=nrow(dfData)), dfData[,-5]))
colnames(X) = colnames(mCoef)
head(X)
ivPredict = plogis(mypred(colMeans(mCoef), list(mModMatrix=X))[,1])
xyplot(ivPredict ~ fGroups, xlab='Actual Group', main= 'Binomial 10 Variables',
       ylab='Predicted Probability of Pre-term',
       data=dfData)

f1 = fit.stan # without guess
# f2 = fit.stan # with guess
## go back to refit with less number of variables
plot(compare(f1, f2))
compare(f1, f2)
## model without guess parameter appears to be better as it has minimal effect in this case
## remake figure with full data
## this step requires to refit the model without scaling, as 
## data unmatched data is also being used here
fit.stan = f1
df = data.frame((lData.train.full$data[,colnames(mCoef)[-1]]))
df$fGroups = lData.train.full$covariates$fGroups:lData.train.full$covariates$matched
head(df)
ncol(df)
X = as.matrix(cbind(rep(1, times=nrow(df)), df[,-5]))
colnames(X) = colnames(mCoef)
head(X)
ivPredict = plogis(mypred(colMeans(mCoef), list(mModMatrix=X))[,1])
xyplot(ivPredict ~ fGroups, xlab='Actual Group', main= 'Binomial 10 Variables with Unmatched',
       ylab='Predicted Probability of Pre-term',
       data=df)
## which samples from the control group are not predicted correctly
iOutliers = which(ivPredict > 0.6 & df$fGroups == '0:dropped')
df = lData.train.full$covariates[iOutliers,]

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
colnames(dfData)

plot(dfData$Gardnerella_vaginalis, fGroups.jitt, pch=20, xlab='Gardnerella_vaginalis', ylab='Probability of preterm',
     main='Prediction of Preterm')
x = seq(min(dfData$Gardnerella_vaginalis), max(dfData$Gardnerella_vaginalis), length.out = 100)
m = cbind(1, matrix(colMeans(dfData[,-5]), nrow = length(x), ncol = 4, byrow = T))
m[,2] = x
c = colMeans(extract(fit.stan)$betas)
lines(x, plogis(m %*% c), col='black')
# m = cbind(1, x, min(dfData$TRANCE), min(dfData$CGAsampling), min(dfData$CGSamDel))
# lines(x, plogis(m %*% c), col='red')
# m = cbind(1, x, max(dfData$TRANCE), max(dfData$CGAsampling), max(dfData$CGSamDel))
# lines(x, plogis(m %*% c), col='green')
# legend('left', legend = c('Min', 'Average', 'Max'), fill=c('red', 'black', 'green'))

colnames(dfData)
plot(dfData$Prevotella_oris, fGroups.jitt, pch=20, xlab='Prevotella_oris', ylab='Probability of preterm',
     main='Prediction of Preterm')
x = seq(min(dfData$Prevotella_oris), max(dfData$Prevotella_oris), length.out = 100)
m = cbind(1, matrix(colMeans(dfData[,-5]), nrow = length(x), ncol = 4, byrow = T))
m[,3] = x
head(m)
c = colMeans(extract(fit.stan)$betas)
lines(x, plogis(m %*% c), col='black')

colnames(dfData)
plot(dfData$Megasphaera_genomosp., fGroups.jitt, pch=20, xlab='Megasphaera_genomosp.', ylab='Probability of preterm',
     main='Prediction of Preterm')
x = seq(min(dfData$Megasphaera_genomosp.), max(dfData$Megasphaera_genomosp.), length.out = 100)
m = cbind(1, matrix(colMeans(dfData[,-5]), nrow = length(x), ncol = 4, byrow = T))
m[,4] = x
c = colMeans(extract(fit.stan)$betas)
lines(x, plogis(m %*% c), col='black')



# colnames(dfData)
# plot(dfData$Dialister_micraerophilus, fGroups.jitt, pch=20, xlab='Dialister_micraerophilus', ylab='Probability of preterm',
#      main='Prediction of Preterm')
# x = seq(min(dfData$Dialister_micraerophilus), max(dfData$Dialister_micraerophilus), length.out = 100)
# m = cbind(1, matrix(colMeans(dfData[,-5]), nrow = length(x), ncol = 4, byrow = T))
# m[,4] = x
# c = colMeans(extract(fit.stan)$betas)
# lines(x, plogis(m %*% c), col='black')
# 
# colnames(dfData)
# plot(dfData$Propionibacterium_sp., fGroups.jitt, pch=20, xlab='Propionibacterium_sp.', ylab='Probability of preterm',
#      main='Prediction of Preterm')
# x = seq(min(dfData$Propionibacterium_sp.), max(dfData$Propionibacterium_sp.), length.out = 100)
# m = cbind(1, matrix(colMeans(dfData[,-5]), nrow = length(x), ncol = 4, byrow = T))
# m[,5] = x
# c = colMeans(extract(fit.stan)$betas)
# lines(x, plogis(m %*% c), col='black')


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
