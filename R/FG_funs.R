##  utility: define a nearly-diagonal matrix
##  (alpha off-diagonal plus diagonal values)
alphaMat <- function(diagVal,alpha) {
    d <- matrix(alpha,nrow=length(diagVal),ncol=length(diagVal))
    diag(d) <- diagVal
    d
}

## distance from diagonality
phi_fun <- function(F) prod(diag(F))/det(F)

## combined distance of a list of matrices
Phi_fun <- function(mL,B=diag(nrow(mL[[1]]))) {
    phivec <- sapply(mL,function(M) phi_fun(t(B) %*% M %*% B))
    prod(phivec)
}

## is it pos def?
posDef <- function(M) {
 all(Re(eigen(M,only.values=TRUE)$values)>0)
}

new_cpc <- function(matList,n=rep(1,length(matList)),
                    verbose=FALSE,
                    eps.F=1e-6,eps.G=1e-6,
                    Bstart=diag(p)) {
    ## FIXME: test matList for condition number (see FG remark 5.1)
    dimmat <- sapply(matList,dim)
    if (any(dimmat[1,]!=dimmat[2,])) stop("all matrices must be square")
    if (length(unique(dimmat[1,]))>1) stop("all matrices must have the same dimensions")
    if (!all(sapply(matList,posDef))) stop("all matrices must be positive definite")
    ## FIXME: check for symmetry
    ## FIXME: check for already-diagonal matrices
    phi_start <- sapply(matList,phi_fun)
    p <- nrow(matList[[1]])  ## matrix dimensions
    if (p==1) return(matrix(1,nrow=1))  ## trivial case
    k <- length(matList)     ## number of matrices
    
    
    ##step F0: We create a pxp matrix as our initial approximation
    B=Bstart
    
    ##step F2: Create all pairs of columns
    Phi=Phi_fun(matList,B)
    Phi_previous=Inf

    iter <- 0
    
    while ((Phi_diff <- (Phi_previous-Phi))>eps.F || 
             (iter==1 && Phi_diff==0)) {
      if (iter==1 && Phi_diff==0) {
             ## FG remark 5.3, 5.4: sometimes diag(p) is a bad starting point
             combMat <- Reduce("+",mapply("*",matList,n,SIMPLIFY=FALSE))
             B <- eigen(combMat)$vec
         }
        iter <- iter+1
        
        Bprevious=B
        Phi_previous=Phi
        
        for (j in 2:p) {
            for (l in 1:(j-1)) {

                if (verbose) cat("F-step: ",j,l,"\n")
                ## j <- 2; l <- 1
                H <- cbind(B[,l],B[,j]) ## or H <- B[,c(l,j)](
                
                ##T is the matrix for holding the T matrices in Step F21
                
                ## or we could organize this as an array ... tensor package ???
                ## T=array(1,dim=c(2,2,k))
                T=matrix(1,nrow=4,ncol=k)
                
                ## Step F21: Create k T matrices by extracting from "matList"
                
                for (i in 1:k) {
                    ## matrix(unlist(matList[a])) == matList[[a]]
                    ## matList[[1:2]] == 0 ??? WHY ???
                    T1.1=drop(t((B[,l]))%*%matList[[i]]%*%(B[,l]))
                    T2.1=drop(t((B[,j]))%*%matList[[i]]%*%(B[,l]))
                    T1.2=drop(t((B[,l]))%*%matList[[i]]%*%(B[,j]))
                    T2.2=drop(t((B[,j]))%*%matList[[i]]%*%(B[,j]))
                    
                    ## for the future: is there a tensor product that will do this all at once?
                    
                    ##Store 2x2 matrices in column form in matrix E for easy extraction in later steps
                    ## maybe call this matrix (or array) T instead of E
                    T[,i]=matrix(c(T1.1,T2.1,T1.2,T2.2),nrow=4,ncol=1)
                }
                
                ##Step F22: Perform G algorithm
                
                ##Step G0
                Q=diag(2)
                
                ##Q previous is used in determining when to stop G-algorithm
                d.w=matrix(1,ncol=2)
                ## G2: Compute the deltas for each T matrix
                if (verbose) cat("G-step:\n")
                ## V=1
                V <- Inf
                eps.G <- 1e-6
                T2 <- T
                ctr <- 0
                while(V>eps.G){
                    if (verbose) cat(" current V: ",V,"\n")
                    ctr <- ctr+1
                    Qprevious=Q
                    for (w in 1:k) {
                        for(x in 1:2) {
                            d.w[x]=drop(t(Q[,x])%*%matrix(c(T[,w]),nrow=2,ncol=2)%*%Q[,x])
                        }
                        alpha=(d.w[1]-d.w[2])/(d.w[1]*d.w[2])
                        T2[,w]=n[w]*alpha*T[,w]
                    }
                    
                    ## can take a shortcut here: rowSums(E) -> vector of the sums of each row
                    Tf=matrix(rowSums(T2),2)
                    Teigen=eigen(Tf)
                    ## want the arrangement that gives the minimum distance from Qprev
                    ## we can swap columns or multiply one or both by -1
                    ## browser()
                    if(Teigen$vec[1,1]<0){Teigen$vec[,1]<-(-Teigen$vec[,1])}
                    if(Teigen$vec[1,2]<0){Teigen$vec[,2]<-(-Teigen$vec[,2])}
                    Q1=Teigen$vec[,1:2]
                    Q2=Teigen$vec[,2:1]
                    Qf1=Qprevious-Q1
                    Qf2=Qprevious-Q2
                    V1=norm(Qf1,type=c("f"))
                    V2=norm(Qf2,type=c("f"))
                    V=min(V1,V2)
                    if (V==V1) {
                      Q=Q1
                    }else if (V==V2) {
                      Q=Q2
                    }
                    ## or: Q <- Teigen$vec[,1:2]
                }  ## while V<eps.G
                
                if ((-pi/2)<=acos(Q1[1,1]) & acos(Q1[1,1]) <=(pi/2)) {
                    Q=Q1
                } else {Q=Q2}
                
                Hstar=H%*%Q
                B[,l]=as.matrix(Hstar[,1])
                B[,j]=as.matrix(Hstar[,2])
                
            } ## loop over i
        }  ## loop over j
        Phi=Phi_fun(matList,B)
        if (verbose) cat(" current Phi: ",Phi,"\n")
    
    }  ## while diff(phi) > eps.F

    if (verbose) cat("yay!\n")
    return(B)
}