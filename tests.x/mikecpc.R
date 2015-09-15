mikecpc <- function(covs,npts,ncp=(ncol(covs[[1]])-1)){
  mlist <- list()
  mlist$evecs <- partial_cpc(covs,npts,ncp)
  mlist$evals <- list()
  mlist$cov <- list()
  ncov <- length(covs)
  for(j in 1:ncov){ mlist$evals[[j]] <- diag(t(mlist$evecs[[j]]) %*% covs[[j]] %*% mlist$evecs[[j]])
  mlist$cov[[j]] <- mlist$evecs[[j]] %*% diag(mlist$evals[[j]])  %*% t(mlist$evecs[[j]])
  }
  mlist$chisqtest <- list()
  mlist$chisqtest$df <- (ncp+1)*(ncp)*(ncov-1)/2
  mlist$chisqtest$chisq <- cpcchisq(mlist$cov,covs,npts)
  return(mlist)
}


cpcchisq <- function(cov1,cov2,npts){
  num <- lapply(cov1,det)
  dem <- lapply(cov2,det)
  ll <- Reduce("+",mapply("*",lapply(mapply("/",num,dem),log),(npts-1)))
  return(ll)
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
  evecs = mikecpc(covs,table(f))$evecs[[1]]
  if (length(evecs) == 1 && is.na(evecs)) 
    return(NA)
  evec = evecs[, 1]
if (center) 
  x = scale(x, center = TRUE, scale = FALSE)
t(bpmat(evec) %*% t(x))
  
}
