source("mikecpc.R")
ff <- fullcpc(test_marten,n_marten)
ff2 <- fullcpc(test_marten,n_marten,2)

evals <- function(mat){
  aa <- lapply(mat,eigen)
  ## a <- lapply(aa,"[[","values")
  a <- list()
  for(i in 1:length(mat)){
    a[[i]] <- aa[[i]]$values
  }
  return(a)
  
}

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

ff <- reorder.ee(ff)
ff2 <- reorder.ee(ff2)

## compute pairwise cosine distances among (C)PCs as a way
## of testing this?
    
## part 2: try to figure out what's up with ordering
## in the vole example?

## Vole partial CPC (1 CPC)
v1 <- fullcpc(test_vole,n_vole,1)

## http://stackoverflow.com/questions/2535234/find-cosine-similarity-in-r
cosineDist <- function(x){
  as.dist(1 - x%*%t(x)/(sqrt(rowSums(x^2) %*% t(rowSums(x^2))))) 
}

cdist2 <- function(x,y) x %*% y / sqrt(x%*%x * y%*%y)
cdistmat <- function(x,y) {
    n1 <- ncol(x)
    n2 <- ncol(y)
    res <- matrix(NA,n1,n2)
    for (i in 1:n1) {
        for (j in 1:n2) {
            res[i,j] <- cdist2(x[,i],y[,j])
        }
    }
    res
}
tmpf <- function(i,j) {
    round(abs(zapsmall(cdistmat(v1$evecs[[i]],v1$evecs[[j]]))),3)
}
tmpf(1,2)  ## very similar
tmpf(1,3)  ## pretty similar
tmpf(1,4)  ## flipped


## Q: 
