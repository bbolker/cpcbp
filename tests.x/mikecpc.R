library(cpcbp)
source('test_examples.R')

mikecpc <- function(covs,npts,ncp=(ncol(covs[[1]])-1)){
  mlist <- list()
  mlist$evecs <- partial_cpc(covs,npts,ncp)
  mlist$evals <- list()
  mlist$cov <- list()
  ncov <- length(covs)
  for(j in 1:ncov){ 
    ## Flury eq. 1.20, p. 70; cov[[j]]==S_j
    mlist$evals[[j]] <- diag(t(mlist$evecs[[j]]) %*% covs[[j]] %*% mlist$evecs[[
j]])
    mlist$cov[[j]] <- mlist$evecs[[j]] %*% diag(mlist$evals[[j]])  %*% t(mlist$e
vecs[[j]])}
  return(mlist)
}

cpcchisq <- function(cov1,cov2,npts){
  num <- sapply(cov1,det)
  dem <- sapply(cov2,det)
  ll <- sum((npts-1)*log(num/dem))  ## n or n-1?
  ## ll <- Reduce("+",mapply("*",lapply(mapply("/",num,dem),log),(npts-1)))
  return(ll)
}

fullcpc <- function(covs,npts,ncp=(ncol(covs[[1]])-1)){
  cpclist <- list()
  n <- ncol(covs[[1]])
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
    mlist$chisqtest <- matrix(c(cpcchisq(cpclist[[ncp]]$cov,cpclist[[n]]$cov,npts), 
  "df","pval"),nrow=1,ncol=3)}
  if(ncp > 1){
    mlist$chisqtest <- matrix(1,nrow=ncp,ncol=3)
    mlist$chisqtest[,2] <- "df"
    mlist$chisqtest[,3] <- "pval"
    for(i in 1:(ncp-1)){
      mlist$chisqtest[i,1] <-  cpcchisq(cpclist[[ncp]]$cov,cpclist[[ncp-i]]$cov,npts)
    }
    mlist$chisqtest[ncp,1] <- cpcchisq(cpclist[[ncp]]$cov,cpclist[[n]]$cov,npts)} 
                            
  
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
