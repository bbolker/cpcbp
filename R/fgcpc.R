##fgwrapper 
fgcpc <- function(covs,n){
  p <- nrow(covs[[1]])
  ncov <- length(covs)
  mlist <- list()
  for(i in 1:ncov){
  mlist$evecs[[i]] <- FGalgorithm(1e-6,1e-6,p,n,covs)
  }
  for(j in 1:ncov){ 
    ## Flury eq. 1.20, p. 70; cov[[j]]==S_j
    mlist$evals[[j]] <- diag(t(mlist$evecs[[j]]) %*% covs[[j]] %*% mlist$evecs[[j]])
  }
  return(reorder.ee(mlist))
}
