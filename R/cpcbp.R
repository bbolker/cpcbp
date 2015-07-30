require(MASS)

makeSymmetric <- function(m,upper=TRUE) {
  if (upper) {
    m[lower.tri(m)] <- t(m)[lower.tri(m)]
  } else {
    m[upper.tri(m)] <- t(m)[upper.tri(m)]
  }
  m
}

covmat <- function(vars,cors,tol=1e-7) {
  n <- length(vars)
  cors <- rep(cors,length.out=n*(n-1)/2)
  m <- outer(sqrt(vars),sqrt(vars))
  ## is there a more efficient way to create a symmetric matrix?
  m[lower.tri(m)] <- m[lower.tri(m)]*cors
  m = makeSymmetric(m,upper=FALSE)
  if (any(eigen(m)$values<(-tol)))
    warning("non-positive definite matrix")
  m
}

strip.blanks <- function(str) {
 sub(' +$', '', sub('^ +', '', str))
}

skip.blank.lines <- function(con) {
  junk <- readLines(con,1)
  n <- 0
  while(nchar(junk)==0 || length(grep("[^ ]",junk))==0) {
    junk <- readLines(con,1)
    n <- n+1
  }
  pushBack(junk,con)
  n
}

skip.pattern.lines <- function(con,pattern,blank.ok=TRUE) {
  junk <- readLines(con,1)
  n <- 0
  while(length(grep(pattern,junk))>0 ||
        (blank.ok && (nchar(junk)==0 || length(grep("[^ ]",junk))==0))) {
    junk <- readLines(con,1)
    n <- n+1
  }
  pushBack(junk,con)
  n
}

phillips.getpmat <- function(con,n,hdrlines=1) {
  readLines(con,hdrlines) ## strip header
  lines <- readLines(con,n)
  s1 <- strsplit(lines,":")
  labs <- sapply(s1, function(z)strip.blanks(z[1]))
  vals <- t(sapply(sapply(s1,"[",2),
                   function(z)
                   as.numeric(strsplit(z,"[:(),=]")[[1]][c(1,3,5)])))
  dimnames(vals) <- list(labs,c("chi-square","df","p"))
  vals
}

phillips.getvec <- function(con,n=3,hdrlines=1) {
  readLines(con,hdrlines)
  vec = numeric()
  while(length(vec)<n) {
    readLines(con,1) ## skip vector names
    s1 = try(scan(con,n=n, nlines=1,quiet = TRUE))
    if (class(s1)=="try-error") browser()
    vec = c(vec,s1)
  }
  vec
}

phillips.getsubmat <- function(con,n=3) {
  readLines(con,1)
  mstr <- readLines(con,n=n)
  t(sapply(strsplit(mstr," +"),
         function(z)as.numeric(z[-(1:2)])))
}

phillips.getmat <- function(con,n=3,hdrlines=1) {
  readLines(con,hdrlines)
  m = matrix(nrow=n,ncol=0)
  while(ncol(m)<n) {
    m = cbind(m,phillips.getsubmat(con,n))
    skip.blank.lines(con)
  }
  m
}

phillips.getcolonval <- function(con) {
 str <- readLines(con,1)
 as.numeric(strsplit(str,":")[[1]][2])
}


run.cpc <- function(covs,npts,
                    progdir,
                    progname,
                    ansfn,datfn,
                    outfn) {
  if (missing(progname))
      progname=paste("phillips-cpc-",.Platform$OS.type,".exe",sep="")
  if (!missing(progdir)) {
    progpath = file.path(progdir,progname)
  } else {
    progpath <- system.file("exec",progname,package="cpcbp")
  }
  ngrp = length(covs)
  nvar = ncol(covs[[1]])
  npts = rep(npts,length.out=ngrp)
  outf <- file(description=datfn,open="w") ## open for  writing
  on.exit(close(outf),add=TRUE)
  writeLines(c(paste("group",1:ngrp,sep="",collapse=" "),
               paste("var",1:nvar,sep="",collapse=" ")),con=outf)
  for (i in 1:ngrp) {
    writeLines(as.character(npts[i]),outf)
    write.table(covs[[i]],file=outf,
              row.names=FALSE,col.names=FALSE,quote=FALSE)
  }
  flush(outf) ## needed now that we don't close until exit!!
  if (.Platform$OS.type=="unix") {
    ans <- file(description=ansfn,open="w") ## open for writing
    on.exit(close(ans),add=TRUE)
    writeLines(c(datfn,outfn,"y","n","n",""),con=ans)
    flush(ans)
    system(paste(progpath,"<",ansfn,"1> /dev/null 2>&1"))
  } else {
    ## triggers codetools error (no input argument in Linux) but oh well ...
    ## ? assume shell() has same arguments as system() ?
    print(progpath)
    system(paste("\"",progpath,"\"",sep=""),
           input=c(basename(datfn),basename(outfn),
                      "y","n","n",""))
  }
}

phillips.cpcvec <- function(x,f,covs,
                     npts,progdir,
                     progname,ansfn,datfn,outfn,
                     unlink.temp=TRUE,
                     use="complete.obs") {
  if (missing(f)) {
    f <- x[[2]]
    x <- x[[1]]
  }
  P1 = phillips.cpc(x=x,f=f,covs=covs,
    npts=npts,
    progdir=progdir,progname=progname,
    ansfn=ansfn,
    datfn=datfn,
    outfn=outfn,unlink.temp=unlink.temp,
    use=use)
  if (length(P1)==1 && is.na(P1)) NA else P1$evecs.CPC
}

phillips.cpc <- function(x,f,covs,
                         npts,
                         progdir,
                         progname,
                         ansfn,
                         datfn,
                         outfn,
                         unlink.temp=TRUE,
                         use="complete.obs",
                         verbose=FALSE) {
  if (missing(f)) {
    f <- x[[2]]
    x <- x[[1]]
  }
  tmpdirf = function(x,base) {
    if(missing(x)) switch(.Platform$OS.type,
                          unix=tempfile(base),
                          windows=tempfile(base,tmpdir=getwd())) else x
  }
  ansfn = tmpdirf(ansfn,"cpcans")
  datfn = tmpdirf(datfn,"cpcdat")
  outfn = tmpdirf(outfn,"cpcout")
  if (missing(covs)) {
    if (missing(x) && missing(f))
      stop("must specify either covariances or data and grouping factor")
    datalist = split.data.frame(x,f)
    covs <- lapply(datalist,cov,use=use)
    npts <- sapply(datalist,nrow)
  } else {
    if (missing(npts)) {
      if (missing(x) || missing(f)) {
        stop("must specify either number of points per group or data and grouping factor")
        datalist = split.data.frame(x,f)
        npts <- sapply(datalist,nrow)
      }
    }
  }
  run.cpc(covs,npts,progdir,progname,ansfn,datfn,outfn)
  ngrp = length(covs)
  nvar = ncol(covs[[1]])
  if (length(grep("ill-formed",readLines(outfn)))>0) {
    r = NA
  } else r = read.cpc(outfn,ngrp,nvar,verbose=verbose)
  if (unlink.temp) {
    unlink(ansfn)
    unlink(outfn)
    unlink(datfn)
  }
  return(r)
}

phillips.getsection <- function(con,hdrlines,ntests,
                                 nvar,
                                 ngrp,
                                 matvals=FALSE, ## additional value with each matrix?
                                 grpevals=FALSE,## evals for each group?
                                 grps=TRUE,
                                 nevecs=1,
                                verbose=FALSE) {   ## multiple evals/mats?
    readLines(con,hdrlines)
    if (!grps) ngrp <- 1
    if (!grpevals) nevals <- 1 else nevals <- ngrp
    skip.pattern.lines(con,"^\\*\\*\\*")
    if (verbose) cat("getting critical value and parameter\n")
    crit <- phillips.getcolonval(con)
    par <- phillips.getcolonval(con)
    skip.blank.lines(con)
    if (verbose) cat("getting test p-value matrix\n")
    testmat <- phillips.getpmat(con,ntests)
    skip.blank.lines(con)
    ## readLines(con,1)
    if (verbose) cat("getting eigenvalues\n")
    if (grpevals) readLines(con,1)
    evals <- lapply(1:nevals,
                    function(x) {
                      m <- phillips.getvec(con,n=nvar)
                      skip.blank.lines(con)
                      m
                    })
    if (!grpevals) evals <- evals[[1]]
    evecs <- list()
    if (verbose) cat("getting eigenvector matrices\n")
    for (i in 1:nevecs) {
      skip.blank.lines(con)
      evecs <- append(evecs,list(phillips.getmat(con,n=nvar)))
    }
    if (nevecs==1) evecs <- evecs[[1]]
    if (verbose) cat("getting covariance matrices\n")
    skip.blank.lines(con)
    cov <- lapply(1:ngrp,
                  function(x) {
                    if (matvals) { readLines(con,1) } ## skip prop const
                    m <- phillips.getmat(con,n=nvar)
                    skip.blank.lines(con)
                    m
                  })
    if (!grps) cov <- cov[[1]]
    list(crit=crit,par=par,testmat=testmat,evals=evals,
         evecs=evecs,cov=cov)
  }    

read.cpc <- function(outfn,ngrp,nvar,verbose=FALSE) {
  ##
  outf <- file(outfn,open="r")
  on.exit(close(outf),add=TRUE)
  ## get data for pooled/equal
  if (verbose) cat("getting pooled information\n")
  pool <- phillips.getsection(outf,hdrlines=2,nvar+1,nvar,ngrp,
                              matvals=FALSE,grpevals=FALSE,grps=FALSE,
                              verbose=verbose)
  ## Extract data for proportionality
  if (verbose) cat("getting proportionality information\n")
  propdat <- phillips.getsection(outf,1,nvar,nvar,ngrp,
                                  matvals=TRUE,grpevals=FALSE,grps=TRUE,
                                 verbose=verbose)
  ## Extract data for common principal components
  if (verbose) cat("getting CPC information\n")
  cpcdatlist <- list()
  cpcdatlist <- append(cpcdatlist,
                       list(phillips.getsection(outf,1,nvar-1,nvar,ngrp,
                                            matvals=FALSE,grpevals=TRUE,
                                            grps=TRUE,verbose=verbose)))
  for (i in (nvar-2):1) {
    if (verbose) cat("getting",i,"CPC information\n")
    cpcdatlist <- append(cpcdatlist,
                         list(phillips.getsection(outf,1,i,nvar,ngrp,
                                              matvals=FALSE,grpevals=TRUE,
                                              grps=TRUE,nevecs=ngrp,
                                                  verbose=verbose)))
  }
  rest <- readLines(outf)
  lastline <- rest[length(rest)]
  s = strsplit(lastline," ")[[1]]
  cpc1.pval <- as.numeric(s[length(s)])
  return(list(evecs.CPC=cpcdatlist[[1]]$evecs,cpc1.pval=cpc1.pval,
              datlist=cpcdatlist,pool=pool,propdat=propdat))
}

bpmat <- function(evec) {
  return(diag(length(evec)) - evec %*% t(evec))
}

bpfun <- function(x,f,evec,center=TRUE) {
  if (missing(evec)) {
    if (missing(f)) {
      if (is.list(x)) {
        f = x[[2]]
        x = x[[1]]
      } else stop("must specify either CPC1 or a grouping factor")
    }
    evecs = cpc.options()$cpcvecfun(x,f)
    if (length(evecs)==1 && is.na(evecs)) return(NA)
    evec = evecs[,1]
  }
  if (center) x=scale(x,center=TRUE,scale=FALSE)
  ## should do this by group, not overall mean -- for unbalanced data
  t(bpmat(evec) %*% t(x))
}

bp.contrasts <- function(formula) {
  to = terms(formula)
  interc = as.logical(attr(to,"intercept"))
  x = eval.parent(attr(to,"variables"))[[attr(to,"response")]]
  f = eval.parent(attr(to,"variables"))[[2]]  ## mild hack
  x = bpfun(x,f)
  bperr <- bp.error(x,f)
  n = c(table(f))
  r = list()
  for (i in 1:ncol(x)) {
    tmpform = formula
    tmpform[[2]]=quote(x[,i])
    ## switch environment
    environment(tmpform) <- new.env()
    lm1 = lm(tmpform)
    v = coef(summary(lm1))
    cm = contrasts(f)
    if (interc) {
      cm = cbind(1,cm)
    } else {
      cm = cbind(c(1,rep(0,nrow(cm)-1)),cm)
    }
    ## is this right??
    v[,"Std. Error"] = sqrt((v[,"Std. Error"]^2*sum(n)+
       cm %*% bperr[,i]*n^2)/sum(n))
    v[,"t value"] = v[,"Estimate"]/v[,"Std. Error"]
    v[,"Pr(>|t|)"] = 2*pt(abs(v[,"t value"]),df=sum(n)-1,lower.tail=FALSE)
    r[[i]] = v
  }
  r
}
  
bp.means <- function(x,f,center=FALSE) {
  if (missing(f)) {
    f <- x[[2]]
    x <- x[[1]]
  }
  f <- f[complete.cases(x)]
  x <- x[complete.cases(x),]
  cpcvecfun = cpc.options()$cpcvecfun
  cpc.evecs = cpcvecfun(x,f)
  cpc1 = cpc.evecs[,1]
  if (center) {
      ctr = colMeans(x)
      x = meancorrect(x)
    } else {
      ctr = rep(0,ncol(x))
    }
  bpx = bpfun(x,cpc1)
  datalist = split.data.frame(bpx,f)
  m.bp = t(sapply(datalist,colMeans))
  bpe <- bp.error(x,f,cpcmat=cpc.evecs)
  varm <- t(sapply(datalist,apply,2,function(x){var(x)/length(x)}))
  varm=colSums(varm)
  m.bp = m.bp[2,]-m.bp[1,]
  label1="meandiffs"
  tote = sqrt(varm+bpe)
  L = list(m.bp,sqrt(varm),tote)
  names(L) = c(label1,"sd.raw","sd.corr")
  L
}

pooled.cpc <- function (x, f, use="complete.obs") 
{
  if (missing(f)) {
    f <- x[[2]]
    x <- x[[1]]
  }
  datalist = split.data.frame(as.data.frame(x),f)
  Xsc <- do.call("rbind", lapply(datalist, scale, scale = FALSE))
  eigen(cov(Xsc,use=use))$vectors
}

calc.cpcerr <- function(x,f,
                        cpcmat=NULL,
                        calc.cov=TRUE,
                        use="complete.obs") {
  if (missing(f)) {
    f <- x[[2]]
    x <- x[[1]]
  }
  cpcvecfun = cpc.options()$cpcvecfun
  ngrp = length(levels(f))
  ngrpn = table(f)
  nvar = ncol(x)
  N = nrow(x)
  datalist = split.data.frame(x,f)
  e <- lapply(datalist,function(z) {eigen(cov(z,use=use))$values})
  theta.hat <- matrix(0,nrow=nvar,ncol=nvar)
  for (i in 1:ngrp) { ## fancy: uses "outer"
    theta.i = outer(e[[i]],e[[i]],
      function(x,y) { (x*y)/(x-y)^2})*(N/ngrpn[i])
    theta.hat = theta.hat + 1/theta.i
  }
  theta.hat = 1/theta.hat ## calc. harmonic mean
  if (missing(cpcmat))
    cpcmat <- cpcvecfun(x,f) ## get common principal components
  ## fancy: sum with j!=h is equivalent to matrix mult. with
  ## diag(theta)==0
  diag(theta.hat) <- 0
  if (!calc.cov) {
    ## VARIANCES of all elements with themselves
    return(1/N*(cpcmat^2 %*% theta.hat))
  } else {
    ## VAR-COV matrix of eig. 1 only
    covarr <- matrix(0,nrow=nvar,ncol=nvar)
    for (j in 1:nvar) {
      covarr <- covarr + theta.hat[1,j]*(cpcmat[,j] %*% t(cpcmat[,j]))
    }
    return(covarr/N)
  }
}

bp.vcov <- function(x,f,cpcmat,eigvar) {
  if (missing(f)) {
    f <- x[[2]]
    x <- x[[1]]
  }
  nvar = ncol(x)
  cpcvecfun = cpc.options()$cpcvecfun
  negtol=cpc.options()$neg.bpvar.tol
  if (missing(cpcmat)) cpcmat <- cpcvecfun(x,f)
  if (missing(eigvar)) eigvar <- calc.cpcerr(x,f,calc.cov=TRUE,
                                             cpcmat=cpcmat)
  cpcv1 <- cpcmat[,1] ## principal eigenvector
  bpvar <- matrix(nrow=nvar,ncol=nvar)
  cpcv1 <- cpcmat[,1] ## principal eigenvector
  bpcov <- array(dim=rep(nvar,3))
  for (i in 1:nvar) {
    for (j in 1:nvar) {
      ## fixed bug: square product
      ## added cov. term
      ## fixed bug: transposed cpcmat correctly
      bpcov[i,j,j] <- (cpcv1[i]*cpcv1[j])^2*
        (eigvar[i,i]/cpcv1[i]^2+eigvar[j,j]/cpcv1[j]^2)+
          2*cpcv1[i]*cpcv1[j]*eigvar[i,j]
      if (j<nvar) {
        for (k in (j+1):nvar) {
          bpcov[i,j,k] = 2*cpcv1[i]*cpcv1[j]*eigvar[i,k]+
            2*cpcv1[i]*cpcv1[k]*eigvar[i,j]+
              cpcv1[i]^2*eigvar[j,k]+
                cpcv1[j]*cpcv1[k]*eigvar[i,i]
          bpcov[i,j,k] = bpcov[i,j,k] -
            cpcv1[i]*cpcv1[j]*eigvar[i,k] -
              cpcv1[i]*cpcv1[k]*eigvar[i,j] -
                eigvar[i,j]*eigvar[i,k]
        }
      }
    } ##
    bpcov[i,,] <- makeSymmetric(bpcov[i,,])
  }
  return(bpcov)
}

bp.error <- function(x,f,
                     cpcmat,m,eigvar,
                     use="complete.obs",
                     debug=FALSE,center=TRUE) {
  if (missing(f)) {
    f <- x[[2]]
    x <- x[[1]]
  }
  crossterms=TRUE
  cpcvecfun = cpc.options()$cpcvecfun
  negtol=cpc.options()$neg.bpvar.tol
  npts = table(f)
  ngrp = length(levels(f))
  nvar = ncol(x)
  if (missing(cpcmat)) cpcmat <- cpcvecfun(x,f)
  if (missing(m)) m <- do.call("rbind",by(x,f,colMeans))
  if (missing(eigvar)) eigvar <- calc.cpcerr(x,f,calc.cov=TRUE,
                                             cpcmat=cpcmat)
  bpvar <- matrix(nrow=nvar,ncol=nvar)
  cpcv1 <- cpcmat[,1] ## principal eigenvector
  if (debug) print(bpcov)
  ## for each variable, for each group:
  ## if (debug) print(bperr) ## ??
  if (center) m = scale(m,center=TRUE,scale=FALSE)
  ## mdev = abs(mdev) ## ??? with this, bperr is exactly
                      ## TWICE the sum of bperrmat ... ?
  ## i.e., (2x)^2 rather than x^2+x^2 ...
  ## need to reconstruct why I did this in the first place!
  ## i.e., because we know that the error is the SAME for
  ## each group, not added independently
  ##  so we really DO need totdev to be added --- but how
  ##  do we weight number of points per group?
  ##     (could test with totdev vs. dev to see which one actually
  ##      works better)
  ##      ... I'm worried about the effect on the covariances ...
  bperrmat = matrix(nrow=ngrp,ncol=nvar)
  bpcov = bp.vcov(x,f,cpcmat,eigvar)
  for (j in 1:ngrp) {
    for (i in 1:nvar) {
      bperrmat[j,i] = crossprod(crossprod(bpcov[i,,],m[j,]),m[j,])
    }
  }
  if (debug) print(bperrmat)
  ##
  if (any(bperrmat<(-negtol))) warning("estimated BP variance <0:",
                                       min(bperrmat))
  varnames = colnames(x)
  grpnames = levels(f)
  dimnames(bperrmat) = list(grpnames,varnames)
  return(bperrmat)
}

bp.anova <- function(x,f,verbose=FALSE,
                     use="complete.obs",
                     debug=FALSE,
                     center=TRUE,
                     nsim=10000) {
  if (missing(f)) {
    f <- x[[2]]
    x <- x[[1]]
  }
  cpcvecfun = cpc.options()$cpcvecfun
  nvar <- ncol(x)
  N <- nrow(x)
  ngrpn <- as.numeric(table(f))
  cpc.evecs = cpcvecfun(x,f)
  e11 = cpc.evecs[,1]
  if (center) x = scale(x,center=TRUE,scale=FALSE)
  ## should scale by group means rather than overall
  ##   means (no difference if balanced)
  bperr <- bp.error(x,f,cpcmat=cpc.evecs,
                    use=use,debug=debug)
  bpssq <- colSums(sweep(bperr,1,ngrpn,"*"))
  bpx = bpfun(x,evec=e11)
  alist <- list()
  for (i in 1:nvar) {
    a <- anova(lm(bpx[,i] ~ f))
    a$"Pr(>F)"[1] <- pbperr(a$"F value"[1],bpms=bpssq[i]/N,
                            atab=a,nsim=nsim,lower.tail=FALSE)
    alist[[i]] <- a
  }
  names(alist) <- colnames(x)
  if (!verbose) alist else
  list(alist=alist,bp=bpx,cpc.evecs=cpc.evecs,bperr=bperr,bpssq=bpssq)
}


plot_multigrp <- function(x,f,vars=1:ncol(x),
                          cols,
                          pchs,
                          eqsc=FALSE,xlim=NULL,ylim=NULL,...) {
  if (missing(f)) {
    f <- x[[2]]
    x <- x[[1]]
  }
  ngrp=length(levels(f))
  if (missing(cols)) cols <- 1:ngrp
  if (missing(pchs)) pchs <- 1:ngrp
  ## xlist = split.data.frame(as.data.frame(x),f)
  x = x[,vars]
  if (length(vars)==2) {
    plot(x,col=cols[f],pch=pchs[f],xlim=xlim,ylim=ylim,...)
  } else {
    pairs(x,col=cols[f],pch=pchs[f],xlim=xlim,ylim=ylim,...)
  }
}

plot_dat.theor <- function(x,f,vars=1:2,cols=1:length(levels(f)),theor,
                           lines=TRUE,ellipses=TRUE,...) {
  if (missing(f)) {
    f <- x[[2]]
    x <- x[[1]]
  }
  if (ellipses) require("ellipse")
  plot(x[,vars[1]],x[,vars[2]],
       xlab=colnames(x)[vars[1]],
       ylab=colnames(x)[vars[2]],col=cols[f],...)
  points(theor$mean[,vars[1]],theor$mean[,vars[2]],pch=16,col=cols)
  slope=theor$eigs[vars[2],1]/theor$eigs[vars[1],1]
  int = theor$mean[,vars[2]]-slope*theor$mean[,vars[1]]
  if (lines) mapply(function(int,col) abline(a=int,b=slope,col=col), int,cols)
  if (ellipses) invisible(mapply(function(ctr,col) {
    lines(ellipse::ellipse(theor$varcov,
                  centre=as.numeric(ctr)),
                  col=col)},
                                 split.data.frame(theor$mean[,vars],1:2),
                                 as.list(cols)))
}

sim.theor <- function(vars=c(10,10,10),cors=.8,
                    m1a=rep(0,length(vars)),
                    offset=1, offset2=0) {
  ## ndim <- length(vars)
  ## ngrp = length(offset)
  VC0 <- covmat(vars,cors)
  e0 <- eigen(VC0)
  offset = c(0,offset)
  offset2 = c(0,offset2)
  meanvals = t(mapply(function(O1,O2) {
    m1a+O1*e0$vectors[,1]+O2*e0$vectors[,2]
  },offset,offset2))
  list(mean=meanvals,varcov=VC0,eigs=e0$vec)
}
  
simdata <- function(vars=c(10,10,10),cors=.8,npts=200,
                    seed,m1a=rep(0,length(vars)),
                    offset=1, offset2=0) {
  ## ndim <- length(vars)
  VC0 <- covmat(vars,cors)
  e0 <- eigen(VC0)
  offset = c(0,offset)
  if (missing(offset2)) {
    offset2=rep(0,length(offset))
  } else offset2 = c(0,offset2)
  lens = c(length(offset),length(offset2),length(npts))
  ngrp = max(lens)
  if (!all(lens==1 | lens==ngrp)) stop ("all params must have same length")
  if (length(npts)==1) npts = rep(npts,ngrp)
  if (!missing(seed)) set.seed(seed)
  x = do.call("rbind",mapply(function(n,O1,O2) {
    mvrnorm(n,mu=m1a+O1*e0$vectors[,1]+O2*e0$vectors[,2],Sigma=VC0)
  },npts,offset,offset2,SIMPLIFY=FALSE))
  list(data=x,f=factor(rep(1:length(npts),npts)))
}

meancorrect <- function(x) {
  scale(x,scale=FALSE,center=colMeans(x,na.rm=TRUE))
}

coverfun <- function(p,alpha=0.05,one.tailed=FALSE,na.rm=TRUE) {
  if (na.rm) p <- p[!is.na(p)]
  if (!one.tailed) {
    return(sum(p>alpha/2 & p<(1-alpha/2))/length(p))
  } else {
    return(sum(p>alpha)/length(p))
  }
}

.cpc.options = list(cpcvecfun=phillips.cpcvec,
  neg.bpvar.tol = 5e-3)

cpc.options <- function(cpcvecfun,neg.bpvar.tol) {
  if (!missing(cpcvecfun)) { .cpc.options$cpcvecfun = cpcvecfun }
  if (!missing(neg.bpvar.tol)) { .cpc.options$neg.bpvar.tol = neg.bpvar.tol }
  invisible(.cpc.options)
}


## simulate data from the null hypothesis: means collapsed onto
##   CPC1 axis, covariances equivalent to the 1-CPC estimate from
##   Phillips/Flury
nullsim <- function(x,f,P1,plot.it=FALSE) {
  if (missing(f)) {
    f <- x[[2]]
    x <- x[[1]]
  }
   ngrp = length(levels(f))
   npts = table(f)
   if (missing(P1)) P1 = phillips.cpc(x,f)  
   covs = P1$datlist[[length(P1$datlist)]]$cov  ## extract 1-CPC covs
   cpc1 = P1$evecs.CPC[,1]                      ## CPC1 vector
   x2 = split.data.frame(x,f)
   means = t(sapply(x2,colMeans))
   nullmat = cpc1%*%t(cpc1)
   nullmeans = nullmat%*% t(means)              ## project onto CPC1
   M =  mapply(mvrnorm,npts,as.data.frame(nullmeans),covs,SIMPLIFY=FALSE)
   do.call("rbind",M)                      
}


## extract F statistics
bp.fstats = function(obj) {
  sapply(obj,function(x) x$"F value"[1])
}

## generate null F statistics
nullFstats <- function(x,f,n=1,bpanova=FALSE) {
  if (missing(f)) {
    f <- x[[2]]
    x <- x[[1]]
  }
  if (n==1) {
    bpx = nullsim(x,f)
    anovafun = if (bpanova) {
      bp.anova
    } else {
      function(x,f) {
        lapply(x,function(z)
               anova(lm(z ~ f)))
      }
    }
    bpn = bpfun(bpx,f)  ## back-project
    if (length(bpn)==1 && is.na(bpn)) return(rep(NA,nrow(bpx)))
    bp.fstats(anovafun(as.data.frame(bpn),f))
  } else t(replicate(n,nullFstats(x,f,1)))
}

bp.pvals <- function(obj) sapply(obj,function(x)x$"Pr(>F)"[1])

rbperr = function(n,df1,df2,bpms,atab,scaled,dvar=FALSE) {
  if (!missing(atab)) {
    df1 = atab$Df[1]
    df2 = atab$Df[2]
    if (missing(scaled) || scaled) {
      s2.err = atab$"Mean Sq"[2]
      npts = df2+2
      bpms = bpms*(npts/df1)/s2.err
    }
  }
  bperr = rchisq(n,1)*bpms
  if (!dvar) {
    num = rchisq(n,df1)/df1
    denom = rchisq(n,df2)/df2
    return((num+bperr)/denom)
  } else {
    f = rf(n,df1,df2)
    return(f+bperr)
  }
}

qbperr = function(q,df1,df2,bpms,atab,nsim=20000,...) {
  quantile(rbperr(nsim,df1=df1,df2=df2,bpms=bpms,atab=atab,...),
           prob=q)
}

pbperr = function(x,df1,df2,bpms,atab,nsim=20000,lower.tail=TRUE,...) {
  r = rbperr(nsim,df1=df1,df2=df2,bpms=bpms,atab=atab,...)
  if (lower.tail) sum(r<x)/nsim else sum(r>x)/nsim
}

dbperr = function(x,df1,df2,bpms,atab,nsim=20000,...) {
  r = rbperr(nsim,df1=df1,df2=df2,bpms,atab=atab,...)
  d = density(r,from=0)
  approx(d$x,d$y,xout=x)$y
}
