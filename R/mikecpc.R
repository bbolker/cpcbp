hmean <- function(L){
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
  return(partial_cpc(cov,n,q=common-1,B=(diag(common))))
}


cpcchisq <- function(cov1,cov2,npts){
  num <- sapply(cov1,det)
  dem <- sapply(cov2,det)
  ll <- sum((npts-1)*log(num/dem))  ## n or n-1?
  return(ll)
}

cpc_object <- function(covs,npts,ncp=NULL){
  if(is.null(ncp))(ncp=nrow(covs[[1]])-1)
  mlist <- list()
  p <- nrow(covs[[1]])
  k <- length(covs)
  full <- cpc(covs,npts)
  if(ncp == (p-1)){
    mlist$evecs <- full
    mlist$par <- p*(p-1)/2 + k*p}
  else{mlist$evecs <- partial_cpc(covs,npts,q=ncp,B=full[[1]])
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

cpc_LRT <- function(cpc1,cpc2=NULL,covs,npts){
    k <- length(covs)
    p <- nrow(covs[[1]])
    arbpar <- k*( p*(p-1)/2 + p)
    chisqtest <- data.frame(statistic = NA,df = NA,pval = NA)
    if(is.null(cpc2)){
        chisqtest[1,1] <- cpcchisq(cpc1$cov,covs,npts)
        chisqtest[1,2] <- arbpar - cpc1$par
        chisqtest[1,3] <- 1 - pchisq(chisqtest[1,1],chisqtest[1,2])
    }
    else{
        chisqtest[1,1] <- cpcchisq(cpc1$cov,cpc2$cov,npts)
        chisqtest[1,2] <- cpc2$par - cpc1$par
        chisqtest[1,3] <- 1 - pchisq(chisqtest[1,1],chisqtest[1,2])
    }
    return(chisqtest)
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
