ordsample <-
function(n, marginal, Sigma, support=list(), Spearman=FALSE, cormat="discrete", tol = 1e-6)
{
    if (!all(unlist(lapply(marginal, function(x) (sort(x)==x & min(x)>0 & max(x)<1))))) stop("Error in assigning marginal distributions!")
    if(!isSymmetric(Sigma)) stop("Correlation matrix not symmetric")
    if(min(eigen(Sigma)$values)<0) stop("Correlation matrix contains negative eigenvalues")
    if(!all.equal(diag(unname(Sigma)), rep(1, nrow(Sigma)))) stop("Correlation matrix diagonal not equal to 1")
# k=number of variables
k <- length(marginal)
# kj=number of categories for the k variables (vector of k integer numbers)
kj <- numeric(k)
len <- length(support)
for(i in 1:k)
{
kj[i] <- length(marginal[[i]]) + 1
if(len==0)
{
support[[i]] <- 1:kj[i]
}
}
if(cormat=="discrete")
{
Sigmac <- ordcont(marginal=marginal, Sigma=Sigma, support=support, Spearman=Spearman)[[1]]
Sigma <- Sigmac
}
# sample of size n from k-dimensional normal with vector of zero means and correlation matrix Sigma
valori <- mvrnorm(n, rep(0, k), Sigma, tol = tol)
if(n==1) valori <- matrix(valori,nrow=1)
# discretization according to the marginal distributions
for(i in 1:k)
{
valori[,i] <- as.integer(cut(valori[,i], breaks=c(min(valori[,i])-1, qnorm(marginal[[i]]), max(valori[,i])+1)))
valori[,i] <- support[[i]][valori[,i]]
}
return(valori)
}

