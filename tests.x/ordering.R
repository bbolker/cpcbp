source("mikecpc.R")
ff <- fullcpc(test_marten,n_marten)
ff2 <- fullcpc(test_marten,n_marten,2)
hmean <- function(L) {
    ## harmonic mean of a list of vectors
    inv <- lapply(L,function(x) 1/x)
    ## compress lists into a matrix, compute means, invert ...
    1/rowMeans(do.call(cbind,inv))
}
ff2$evals
(hh <- hmean(ff2$evals))
oo <- order(hh,decreasing=TRUE)
lapply(ff2$evecs,
       function(y) y[,oo])


(hh <- hmean(ff$evals))
oo <- order(hh,decreasing=TRUE)
lapply(ff$evecs,
       function(y) y[,oo])


reorder.ee <- function(x) {
    oo <- order(hmean(x$evals,decreasing=TRUE))
    ## reorder all elements of all elements
    ## (evecs, evals) -> (by group)
    ## in decreasing order of harmonic mean
    x <- lapply(x,
                function(z) lapply(z,function(y) y[,oo]))
    return(x)
}
    
    
