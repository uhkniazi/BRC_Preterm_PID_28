# File: preterm.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: preliminary analysis of preterm data
# Date: 21/01/2022

source('header.R')

load('dataExternal/PretermData.rda')

# gestational age - surrogate for pre-term yes no treatment variable
tapply(CGAdelivery, factor(Ccohort), summary)
par(mfrow=c(2,2))
hist2(CGAdelivery[Ccohort==1], CGAdelivery[Ccohort==0], main = 'CGAdelivery', 
      legends = c('pre', 'norm'), legend.pos = 'topleft' )

tapply(CGAsampling, factor(Ccohort), summary)
hist2(CGAsampling[Ccohort==1], CGAsampling[Ccohort==0], main = 'CGAsampling',
      legends = c('pre', 'norm'))

tapply(CGSamDel, factor(Ccohort), summary)
hist2(CGSamDel[Ccohort==1], CGSamDel[Ccohort==0], main = 'CGSamDel',
      legends = c('pre', 'norm'), legend.pos = 'topleft')

t = as.table(xtabs(~ Ccohort + Csite))
barplot(t, beside = T, las=2, col=c('black', 'lightgrey'), xlim=c(0, 20))
legend('topright', legend = c('normal', 'preterm'), 
       fill=c('black', 'lightgrey'))

length(unique(CpatientID))
length(CpatientID)

dfData = data.frame(CpatientID, Ccohort, CGAdelivery, CGAsampling, 
                    CGSamDel, Csite=factor(Csite), fGroups=factor(Ccohort, labels = c('norm', 'pre')))

dim(CMetabolomics)
CMetabolomics[1:10, 1:10]
# metabolomics matrix has very long column names
# replace with shorter names
dfKey = data.frame(key=paste('Met', 1:ncol(CMetabolomics), sep='_'),
                   original_names=colnames(CMetabolomics), stringsAsFactors = F)

# # there are also duplicated names/versions of some metabolites
# # e.g. looking at Adenine, see columns
# i = which(colnames(CMetabolomics) %in% "Adenine")
# > i
# [1]   45  127  154 1881 2911 3055 4472
table(duplicated(dfKey$original_names))
# FALSE  TRUE 
# 2256  4374 
colnames(CMetabolomics) = dfKey$key
lMetabolomics = list(data=CMetabolomics, key=dfKey, 
                     desc=c('column names of original matrix replaced with shorter names in key'))

# similarly replace the names in the other two data matrices 
# to create unique column names
dfKey = data.frame(key=paste('Trans', 1:ncol(CTranscriptomics), sep='_'),
                   original_names=colnames(CTranscriptomics), stringsAsFactors = F)
colnames(CTranscriptomics) = dfKey$key
lTranscriptomics = list(data=CTranscriptomics, key=dfKey,
                        desc='column names replaced with unique names, see key')

dfKey = data.frame(key=paste('Prot', 1:ncol(CProteomics), sep='_'),
                   original_names=colnames(CProteomics), stringsAsFactors = F)
identical(colnames(CProteomics), dfKey$original_names)
colnames(CProteomics) = dfKey$key
lProteomics = list(data=CProteomics,
                   key = dfKey,
                   desc='column names replaced with unique names, see key')

lData = list(covariates=dfData,
             transcriptomics=lTranscriptomics,
             proteomics=lProteomics,
             metabolomics=lMetabolomics)
lData$source='data downloaded from https://nalab.stanford.edu/multiomicsmulticohortpreterm/'
lData$publication = 'https://pubmed.ncbi.nlm.nih.gov/33337494/'

dir.create('dataUpload')
save(lData, file='dataUpload/stanford_data_rlist.rds')
