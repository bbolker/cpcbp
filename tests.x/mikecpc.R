library(cpcbp)
source('test_examples.R')


hmean <- function(L) {
  ## harmonic mean of a list of vectors
  inv <- lapply(L,function(x) 1/x)
  ## compress lists into a matrix, compute means, invert ...
  1/rowMeans(do.call(cbind,inv))
}

reorder.ee <- function(x) {
  oo <- order(hmean(x$evals),decreasing=TRUE)
  ordfun <- function(x) {
    if (is.matrix(x)) x[,oo] else x[oo]
  }
  ## reorder all elements of all elements
  ## (evecs, evals) -> (by group)
  ## in decreasing order of harmonic mean
  for (e in c("evecs","evals")) {
    x[[e]] <- lapply(x[[e]], ordfun)
  }
  return(x)
}

cpc <- function(cov,n){
  common <- nrow(cov[[1]])
  return(partial_cpc(cov,n,q=common -1,B=(diag(common))))
}

mikecpc <- function(covs,npts,ncp=NULL){
  if(is.null(ncp))(ncp=nrow(covs[[1]])-1)
  mlist <- list()
  p <- nrow(covs[[1]])
  k <- length(covs)
  full <- cpc(covs,npts)
  if(ncp == (p-1)){
    mlist$evecs <- full
    mlist$par <- p*(p-1)/2 + k*p}
  else{mlist$evecs <- partial_cpc(covs,npts,q=ncp,B=mikecpc(covs,npts)$evecs[[1]])
  mlist$par <- p*(p-1)/2 + k*p + (k-1)*(p-ncp)*(p-ncp-1)/2}
  mlist$evals <- list()
  mlist$cov <- list()
  ncov <- length(covs)
  for(j in 1:ncov){ 
    ## Flury eq. 1.20, p. 70; cov[[j]]==S_j
    mlist$evals[[j]] <- diag(t(mlist$evecs[[j]]) %*% covs[[j]] %*% mlist$evecs[[
j]])
    mlist$cov[[j]] <- mlist$evecs[[j]] %*% diag(mlist$evals[[j]])  %*%
        t(mlist$evecs[[j]])}
  return(reorder.ee(mlist))
}

cpcchisq <- function(cov1,cov2,npts){
  num <- sapply(cov1,det)
  dem <- sapply(cov2,det)
  ll <- sum((npts-1)*log(num/dem))  ## n or n-1?
  ## ll <- Reduce("+",mapply("*",lapply(mapply("/",num,dem),log),(npts-1)))
  return(ll)
}

fullcpc <- function(covs,npts,ncp=NULL){
  if(is.null(ncp))(ncp = nrow(covs[[1]]) -1)
  cpclist <- list()
  k <-length(covs)
  p <- nrow(covs[[1]])
  arbpar <- k*( p*(p-1)/2 + p )
  n <- nrow(covs[[1]])
  for(i in 1:(n-1)){
    cpclist[[i]] <- mikecpc(covs,npts,ncp = i)
  }
  
  mlist <- list()
  mlist$evecs <- cpclist[[ncp]]$evecs
  mlist$evals <- cpclist[[ncp]]$evals
  mlist$cov <- cpclist[[ncp]]$cov
#  mlist$chisqtest <- matrix(1,nrow=ncp,ncol=3)
  cpclist[[n]] <- NA
  cpclist[[n]]$cov <- covs
  if(ncp == 1){
    mlist$chisqtest <- data.frame(statistic = cpcchisq(cpclist[[ncp]]$cov,cpclist[[n]]$cov,npts), 
  df = arbpar - cpclist[[ncp]]$par, pval = NA)}
  if(ncp > 1){
    mlist$chisqtest <- data.frame(statistic=rep(1,ncp),df=NA,pval=NA)

    for(i in 1:(ncp-1)){
      mlist$chisqtest[i,1] <-  cpcchisq(cpclist[[ncp]]$cov,cpclist[[ncp-i]]$cov,npts)
      mlist$chisqtest[i,2] <- cpclist[[ncp-i]]$par - cpclist[[ncp]]$par
    }
    mlist$chisqtest[ncp,1] <- cpcchisq(cpclist[[ncp]]$cov,cpclist[[n]]$cov,npts)
    mlist$chisqtest[ncp,2] <- arbpar - cpclist[[ncp]]$par} 
                            
  
  return(mlist)
}

mikepoolcpc <- function(covs,npts){
  mlist <- list()
  poolvar <- Reduce('+',mapply('*',covs,(npts-1),SIMPLIFY = FALSE))/(sum(npts)-length(npts))
  mlist$evals <- diag(t(eigen(poolvar)$vectors) %*% poolvar %*% eigen(poolvar)$vectors)
  mlist$evecs <- eigen(poolvar)$vectors
  mlist$cov <- poolvar
  return(mlist)}


mikebp <- function(x,f,evec, center = TRUE)if (missing(evec)) {
  if (missing(f)) {
    if (is.list(x)) {
      f = x[[2]]
      x = x[[1]]
    }
    else stop("must specify either CPC1 or a grouping factor")
  }
  
  covs <- lapply(split(data.frame(x),f),cov)
  evecs = fullcpc(covs,table(f))$evecs[[1]]
  if (length(evecs) == 1 && is.na(evecs)) 
    return(NA)
  evec = evecs[, 1]
if (center) 
  x = scale(x, center = TRUE, scale = FALSE)
t(bpmat(evec) %*% t(x))
  
}
