## transformation
## https://homepages.inf.ed.ac.uk/rbf/HIPR2/pixexp.htm
fTransform = function(x, c=1, r = 0.5){
  return(c * (x^r))
}

df = read.csv(file.choose(), header=T)
m = as.matrix(df)
x = as.vector(m)
