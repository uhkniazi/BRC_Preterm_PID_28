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

## select variables showing average difference
p.vals = lapply(1:ncol(lData.train$data), function(x){
  df = data.frame(y=lData.train$data[,x], d=lData.train$covariates$fGroups, c=lData.train$covariates$Csite)
  f = lm(y ~ d + c, data=df)
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

## propensity score calculation
f = glm(fGroups ~ Csite, data=lData.train$covariates, family=binomial(link='logit'))
ivPropensity = fitted(f)

df = data.frame(d=lData.train$data[,1], y=lData.train$covariates$fGroups, c=lData.train$covariates$Csite, p=ivPropensity)

## try some model fits with covariates
fit.1 = glm(y ~ d, data=df, family='binomial')
fit.2 = glm(y ~ d + c, data=df, family='binomial')
fit.3 = glm(y ~ d + p, data=df, family='binomial')

summary(fit.1); summary(fit.2); summary(fit.3)

## regress out the predictor variable
fit.reg = lm(d ~ p, data=df)
p.resid = resid(fit.reg)
fit.4 = glm(y ~ p.resid, data=df, family='binomial')
summary(fit.1); summary(fit.2); summary(fit.3); summary(fit.4)

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

########## adjust the data first
m = apply(lData.train$data[,cvTopVariables.lm], 2, function(x){
  return(resid(lm(x ~ ivPropensity)))
})

dfData = data.frame(m)
fGroups = lData.train$covariates$fGroups

oVar.rAdj = CVariableSelection.RandomForest(dfData, fGroups, boot.num = 100)

plot.var.selection(oVar.rAdj)

## there does seem to be a difference in the top variables after adjustment
