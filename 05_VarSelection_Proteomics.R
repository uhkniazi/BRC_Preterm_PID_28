# File: 05_VarSelection_Proteomics.R
# Auth: umar.niazi@kcl.ac.uk
# Date: 12/04/2022
# Desc: variable selection for proteomics data from stanford dataset

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

### reduce the number of variables with various steps
## remove variables with 0 sd i.e. not changing 
s = apply(mCounts, 1, sd)
summary(s)
s = which(s == 0)
length(s)

lData.train = list(data=t(mCounts), covariates=dfMeta)

rm(lData)
## select variables showing average difference
p.vals = lapply(1:ncol(lData.train$data), function(x){
  df = data.frame(y=lData.train$data[,x], d=lData.train$covariates$fGroups, c=lData.train$covariates$Csite,
                  e=lData.train$covariates$CGAsampling, f=lData.train$covariates$CGSamDel)
  f = lm(y ~ d + c + e + f, data=df)
  s = summary(f)$coefficients
  return(s['dpre', 4])
})

names(p.vals) = colnames(lData.train$data)
dfPvals = do.call(rbind, p.vals)
dfPvals = cbind(dfPvals, p.adjust(dfPvals[,1], method = 'BH'))
colnames(dfPvals) = c('pvalue', 'p.adj')
dfPvals = data.frame(dfPvals)
f = which(dfPvals$pvalue < 0.01)
length(f)
cvTopVariables.lm = rownames(dfPvals)[f]

# additional variable to map the sample numbers as 
# we may drop some after matching and overlap
lData.train$covariates$index = 1:nrow(lData.train$covariates)
## propensity score calculation to match samples
cor(lData.train$covariates$CGAsampling, lData.train$covariates$CGSamDel)
colnames(lData.train$covariates)
df = data.frame(lData.train$covariates[, c(7 ,4, 5)])
df$CGAsampling = scale(df$CGAsampling); df$CGSamDel = scale(df$CGSamDel)
f = glm(fGroups ~ CGAsampling + CGSamDel, data=df, family=binomial(link='logit'))
summary(f)
# propensity score
ivPropensity = fitted(f)
lData.train$propensity = ivPropensity
fG = lData.train$covariates$fGroups
hist2(ivPropensity[fG=='pre'], 
      ivPropensity[fG=='norm'], legends = c('pr', 'n'))

## as there is no overlap, use the covariates only in overlapping areas
# ## bin the propensity score vector
# ## to find overlapping bins from both pre and norm groups
# b = cut(ivPropensity, breaks = 5)
# iNorm = b[fG=='norm']
# iIndex = which(b %in% unique(iNorm))
# ivPropensity.sub = ivPropensity[iIndex]
# fG.sub = fG[iIndex]
# hist2(ivPropensity.sub[fG.sub=='pre'], 
#       ivPropensity.sub[fG.sub=='norm'], legends = c('pr', 'n'))
# table(fG)
# table(fG.sub)
# ## this seems to leave almost no pre observations behind
## repeat on individual covariates to choose common support regions
ivPropensity = lData.train$covariates$CGSamDel

hist2(ivPropensity[fG=='pre'], 
      ivPropensity[fG=='norm'], legends = c('pr', 'n'))

## bin the propensity score vector
## to find overlapping bins from both pre and norm groups
b = cut(ivPropensity, breaks = 10)
levels(b)
iNorm = b[fG=='norm']
iPre = b[fG == 'pre']
iInter = sort(intersect(iPre, iNorm))
iIndex = which(b %in% iInter)
ivPropensity.sub = ivPropensity[iIndex]
fG.sub = fG[iIndex]
hist2(ivPropensity.sub[fG.sub=='pre'], 
      ivPropensity.sub[fG.sub=='norm'], legends = c('pr', 'n'))
table(fG)
table(fG.sub)

## this match is more useful overlap - use this index to subset the data
lData.train.original = lData.train
lData.train$data = lData.train$data[iIndex,]
lData.train$covariates = lData.train$covariates[iIndex,]
# perform matching using propensity score
# library(MatchIt)
z = as.numeric(lData.train$covariates$fGroups) - 1
df = data.frame(fGroups = z, Csite=lData.train$covariates$Csite,
                CGAsampling=scale(lData.train$covariates$CGAsampling),
                CGSamDel = scale(lData.train$covariates$CGSamDel))
f = glm(fGroups ~ CGAsampling + CGSamDel, data = df,
        family=binomial(link='logit'))
summary(f)
df$propensity=fitted(f)
hist(df$propensity)
hist2(df$propensity[df$fGroups==1], df$propensity[df$fGroups==0])
library(Matching)
# matches = matchit(fGroups ~ CGAsampling + CGSamDel, data = df, 
#                   replace = T)
# summary(matches)
colnames(df)
#matches = Match(Tr=df$fGroups, X = df$propensity, estimand = 'ATE')
matches = Match(Tr=df$fGroups, X = as.matrix(df[,c(3, 4)]), estimand = 'ATE', replace = T)
summary(matches)
iIndex = c(matches$index.treated, matches$index.control)

lData.train$data = lData.train$data[iIndex,]
lData.train$covariates = lData.train$covariates[iIndex,]

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
df = data.frame(c=lData.train$covariates$CGSamDel,
                f=lData.train$covariates$fGroups)
hist2(df$c[df$f=='pre'], df$c[df$f=='norm'], main = 'CGSamDel', 
      legends = c('pre', 'norm'), legend.pos = 'topleft' )

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
dfData = data.frame(lData.train$data[,cvTopVariables.lm])
dim(dfData)
dfData$fGroups = lData.train$covariates$fGroups
table(dfData$fGroups)
rm(fGroups)
levels(dfData$fGroups)
lData = list(resp=ifelse(dfData$fGroups == 'pre', 1, 0), mModMatrix=model.matrix(fGroups ~ 1 + ., data=dfData))

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
             las=2, xaxt='n', col='grey', main='Top Variables - binomial', sub='adjusted')
axis(1, at = l2, labels = names(m)[1:20], tick = F, las=2, cex.axis=0.7 )

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
xyplot(ivPredict ~ fGroups, xlab='Actual Group', main= 'Binomial Adjusted',
       ylab='Predicted Probability of Pre-term',
       data=dfData)


# ## find correlated variables
dim(dfData)
mData = as.matrix(dfData[,-31])
length(as.vector(mData))
mCor = cor(mData, use="na.or.complete")
library(caret)
image(mCor)
### find the columns that are correlated and should be removed
n = findCorrelation((mCor), cutoff = 0.7, names=T)
data.frame(n)
# sapply(n, function(x) {
#   (abs(mCor[,x]) >= 0.7)
# })
# 
n = sapply(n, function(x) {
  rownames(mCor)[(abs(mCor[,x]) >= 0.7)]
})

cvTopVariables.rf = rownames(CVariableSelection.RandomForest.getVariables(oVar.r))[1:20]
cvTopVariables.bin = names(m)[1:20]
table(cvTopVariables.bin %in% cvTopVariables.rf)
cvTopVariables = unique(c(cvTopVariables.rf, cvTopVariables.bin))
length(cvTopVariables)
## subset selection
dfData = data.frame(lData.train$data[,cvTopVariables])
fGroups = lData.train$covariates$fGroups
oVar.sub = CVariableSelection.ReduceModel(dfData, fGroups, boot.num = 100)
plot.var.selection(oVar.sub)
table(fGroups)
log(40)
# select 3 to 4 variables
cvVar = CVariableSelection.ReduceModel.getMinModel(oVar.sub, size = 4)

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
# cross validation
oCV.lda = CCrossValidation.LDA(dfData[,cvVar], dfData[,cvVar], fGroups, fGroups, level.predict = 'pre',
                               boot.num = 100, k.fold = 10) 

plot.cv.performance(oCV.lda)

## names of these proteins
lData = f_LoadObject('dataUpload/stanford_data_rlist.rds')
dfKey = lData$proteomics$key
table(dfKey$key %in% cvVar)
cvVar.names = dfKey[dfKey$key %in% cvVar, 'original_names'] 

########################################################################
## refit the binomial model and make some figures
########################################################################
dfData = data.frame(lData.train$data_adjusted[,cvVar])
dim(dfData)
## some acrobatics to rename variables
df = dfKey
df = df[df$key %in% cvVar,]
identical(df$key, cvVar)
i = match(cvVar, df$key)
identical(cvVar, df$key[i])
df = df[i,]
identical(colnames(dfData), df$key)
colnames(dfData) = df$original_names

dfData$fGroups = lData.train$covariates$fGroups
table(dfData$fGroups)
rm(fGroups)
levels(dfData$fGroups)
lData = list(resp=ifelse(dfData$fGroups == 'pre', 1, 0), mModMatrix=model.matrix(fGroups ~ 1 + ., data=dfData))

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

### predictive performance of the model
mCoef = extract(fit.stan)$betas2
dim(mCoef)
colnames(mCoef) = c('Intercept', colnames(lData$mModMatrix)[2:ncol(lData$mModMatrix)])
## get the predicted values
## create model matrix
X = as.matrix(cbind(rep(1, times=nrow(dfData)), dfData[,colnames(mCoef)[-1]]))
colnames(X) = colnames(mCoef)
head(X)
ivPredict = plogis(mypred(colMeans(mCoef), list(mModMatrix=X))[,1])
xyplot(ivPredict ~ fGroups, xlab='Actual Group', main= 'Binomial Adjusted, 6 Variables',
       ylab='Predicted Probability of Pre-term',
       data=dfData)

#########################################################################
### repeat model fitting on unadjusted data with adjustment covariate
#########################################################################
dfData = data.frame(lData.train$data[,cvVar])
dim(dfData)
## some acrobatics to rename variables
df = dfKey
df = df[df$key %in% cvVar,]
identical(df$key, cvVar)
i = match(cvVar, df$key)
identical(cvVar, df$key[i])
df = df[i,]
identical(colnames(dfData), df$key)
colnames(dfData) = df$original_names
dfData$site = dfMeta$Csite
dfData$fGroups = lData.train$covariates$fGroups
table(dfData$fGroups)
rm(fGroups)
levels(dfData$fGroups)
lData = list(resp=ifelse(dfData$fGroups == 'pre', 1, 0), mModMatrix=model.matrix(fGroups ~ 1 + ., data=dfData))

# stanDso = rstan::stan_model(file='binomialRegressionSharedCoeffVariance.stan')

lStanData = list(Ntotal=length(lData$resp), Ncol=ncol(lData$mModMatrix), X=lData$mModMatrix,
                 y=lData$resp)

## give initial values
# initf = function(chain_id = 1) {
#   list(betas=rep(0, times=ncol(lStanData$X)), tau=0.5)
# }


fit.stan.2 = sampling(stanDso, data=lStanData, iter=2000, chains=4, pars=c('tau', 'betas2', 'log_lik'), cores=4,# init=initf,
                    control=list(adapt_delta=0.99, max_treedepth = 13))

# save(fit.stan, file='temp/fit.stan.binom_preterm_met.rds')

print(fit.stan.2, c('betas2', 'tau'))
# print(fit.stan, 'tau')
# traceplot(fit.stan, 'tau')

## get the coefficient of interest
mCoef = extract(fit.stan.2)$betas2
dim(mCoef)
## get the intercept
iIntercept = mCoef[,1]
mCoef = mCoef[,-1]
colnames(mCoef) = colnames(lData$mModMatrix)[2:ncol(lData$mModMatrix)]

## coeftab object 
p = paste0('betas2[', 2:7, ']')
plot(coeftab(fit.stan, fit.stan.2), pars=p)
ct.1 = coeftab(fit.stan.2)
rn = rownames(ct.1@coefs)
i = grep('betas', rn)
rownames(ct.1@coefs)[i[-1]] = colnames(mCoef)
rownames(ct.1@se)[i[-1]] = colnames(mCoef)
plot(ct.1, pars=colnames(mCoef))
