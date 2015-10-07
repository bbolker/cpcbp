library(cpcbp)

##phillips cpc without dataset... only takes in the list of covs 
newphil <- function (x=x,covs, npts, progdir, progname, ansfn, datfn, 
                     outfn, unlink.temp = TRUE, use = "complete.obs", verbose = FALSE){ 
  
  tmpdirf = function(x, base) {
    if (missing(x)) 
      switch(.Platform$OS.type, unix = tempfile(base), 
             windows = tempfile(base, tmpdir = getwd()))
    else x
  }
  ansfn = tmpdirf(ansfn, "cpcans")
  datfn = tmpdirf(datfn, "cpcdat")
  outfn = tmpdirf(outfn, "cpcout")
  
  run.cpc(covs, npts, progdir, progname, ansfn, datfn, outfn)
  ngrp = length(covs)
  nvar = ncol(covs[[1]])
  if (length(grep("ill-formed", readLines(outfn))) > 0) {
    r = NA
  }
  else r = read.cpc(outfn, ngrp, nvar, verbose = verbose)
  if (unlink.temp) {
    unlink(ansfn)
    unlink(outfn)
    unlink(datfn)
  }
  return(r)
}

##test1 iris data...change cm to mm measurements to compare with Flury



