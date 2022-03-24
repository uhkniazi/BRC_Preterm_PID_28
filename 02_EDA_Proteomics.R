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
f = apply(mCounts, 2, is.na)
f2 = rowSums(f)
table(f2)
# it appears some rows are NA i.e. proteins not detected? or data absent
f = mCounts[which(f2 > 0), ]
f[,1:10]
dim(mCounts); dim(na.omit(mCounts))
mCounts = na.omit(mCounts)
