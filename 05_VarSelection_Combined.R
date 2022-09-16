# File: 05_VarSelection_Combined.R
# Auth: umar.niazi@kcl.ac.uk
# Date: 16/9/2022
# Desc: variable selection for microbial and metabolite data

source('header.R')

dfMeta = read.csv(file.choose(), header=T, stringsAsFactors = T, row.names = 1)
dfData.met = read.csv(file.choose(), header=T, row.names = 1)
dfData.mic = read.csv(file.choose(), header=T, row.names = 1)

identical(rownames(dfMeta), rownames(dfData.met))
identical(rownames(dfMeta), rownames(dfData.mic))
# one sample is missing a metadata info, remove that
which(dfMeta$BMI_group == '')
dfMeta = dfMeta[-52,]
dfData.met = dfData.met[-52,]
dfData.mic = dfData.mic[-52,]

dfMeta = droplevels.data.frame(dfMeta)
# dfMeta$outcome_numeric = ifelse(dfMeta$outcome == 'TERM', 0, 1)
colnames(dfMeta)
dfMeta = dfMeta[,c(3, 5, 8)]
str(dfMeta)

########## set up the microbiome data
mCounts = as.matrix(dfData.mic)
range(mCounts)
dim(mCounts)
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
mCounts.mic = mCounts

###########################
#### set up the metabolomics data
mCounts = as.matrix(dfData.met)
dim(mCounts)
mCounts = t(mCounts)
mCounts[1:10, 1:10]
range(mCounts)
dim(mCounts); dim(na.omit(mCounts))
plot(density(mCounts))
table(mCounts == 0)
table(mCounts < 1)

## transformation
## https://homepages.inf.ed.ac.uk/rbf/HIPR2/pixexp.htm
fTransform = function(x, c=1, r = 0.5){
  return(c * (x^r))
}

x = as.vector(mCounts)
plot(density(fTransform(x, c=1, r = 0.1)))
m = apply(mCounts, 2, fTransform, c=1, r=0.1)
range(m)
dim(m); dim(mCounts)
mCounts.met = m
#########################################################
## select variables of interest from previous analyses
m1 = mCounts.met[c(), ]
m2 = mCounts.mic[c(), ]
identical(colnames(m1), colnames(m2))
mCounts = rbind(m1, m2)

lData.train = list(data=t(mCounts), covariates=dfMeta)
lData.train$covariates$fGroups = factor(lData.train$covariates$outcome_numeric)
## skipping matching section

######### ML steps
if(!require(downloader) || !require(methods)) stop('Library downloader and methods required')

url = 'https://raw.githubusercontent.com/uhkniazi/CCrossValidation/experimental/CCrossValidation.R'
download(url, 'CCrossValidation.R')

# load the required packages
source('CCrossValidation.R')
# delete the file after source
unlink('CCrossValidation.R')

## subset selection
dfData = data.frame(lData.train$data)
dim(dfData)
fGroups = lData.train$covariates$fGroups
oVar.sub = CVariableSelection.ReduceModel(dfData, fGroups, boot.num = 100)
plot.var.selection(oVar.sub)
table(fGroups)
log(25)
cvVar = CVariableSelection.ReduceModel.getMinModel(oVar.sub, size = 6)
length(cvVar)

### CV with LDA
dfData = data.frame(lData.train$data[,cvVar])
fGroups = lData.train$covariates$fGroups
levels(fGroups)
table(fGroups)
# cross validation
oCV.lda = CCrossValidation.LDA(dfData[,cvVar], dfData[,cvVar], fGroups, fGroups, level.predict = '1',
                               boot.num = 100, k.fold = 10) 

plot.cv.performance(oCV.lda)

########################################################################
## fit the binomial model and make some figures
########################################################################
library(rethinking)
library(rstan)

dfData = cbind(lData.train$data[,cvVar])
head(dfData)
dim(dfData)
dfData = data.frame((dfData))

dfData$fGroups = lData.train$covariates$fGroups
table(dfData$fGroups)
rm(fGroups)
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

print(fit.stan, c('betas'))
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

library(lattice)
xyplot(ivPredict ~ fGroups, xlab='Actual Group', main= 'Metabolites & Microbiome',
       ylab='Predicted Probability of Pre-term',
       data=dfData)

# f1 = fit.stan # without guess
# f2 = fit.stan # with guess
## go back to refit with less number of variables
# plot(compare(f1, f2))
# compare(f1, f2)
# fit.stan = f1 ## without guess has better WAIC

fPredict = ifelse(ivPredict > 0.5, 'P', 'N')
fOriginal = lData.train$covariates$fGroups
table(fPredict, fOriginal)

fPredict = ifelse(ivPredict > 0.3, 'P', 'N')
table(fPredict, fOriginal)

