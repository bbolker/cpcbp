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

