
h_cmpfun <- function(sigmaList,matList,n) {
  sigmadetvec <- sapply(sigmaList,det)
  matdetvec <- sapply(matList,det)
  k=length(matList)
  v=0
  for (i in 1:k) {
    v=v+sum(n[i]*log(sigmadetvec[i]/matdetvec[i]))
  }
  v
}

##' Test whether all CPC are identical
##' 
##' Tests whether the matrices share CPC based on equation 2.18 (Flury, 1984).
##' The null hypothesis is that all groups have the same PC.
##' @param matList list of matrices
##' @param n numeric vector of sample sizes
##' @param B CPC matrix
##' @references Flury, B.N.  Common principal components in K groups.  J.Amer.Statist.Assoc.  1984, 79, 892-898.
##' @examples source("test_examples.R")
##' h_cpc(test_turtle,n_turtle)
##' @export
h_cpc <- function(matList,n,B=new_cpc(matList,n)) {
  k=length(matList)   ## number of groups
  p=nrow(matList[[1]]) ## number of variables (matrix dimension)
  ## n= weights (samples per group)
  chi=matrix(1,k)
  rotate <- function(X,B) t(B) %*% X %*% B
  ## compute similarity score ...
  sigmaList <- list()
  for (i in 1:k) {
    D_hat=diag(diag(rotate(matList[[i]],B)))
    sigmaList[[i]]=rotate(D_hat,t(B))
  }
  chiSum <- h_cmpfun(sigmaList,matList,n)
  p.val <- pchisq(chiSum,df=(k-1)*p*(p-1)/2,lower.tail=FALSE)
  datname <- deparse(substitute(matList))  ## recover _name_ of data argument
  vv <- list(data.name=datname,
             method="Flury test of common principal components",
             statistic=c("chi-sq"=chiSum),p.value=p.val)
  class(vv) <- "htest"
  vv
}




##' Comparing CPC to set of eigenvectors
##' 
##' Tests whether to reject the null hypothesis that a set of q  hypothetical
##' eigenvectors are equal to the first q likelihood eigenvectors in CPC, based
##' on equation 3.3 (Flury, 1986).
##' The value q must be a value of 2 or greater.  For q=1, see the function h_one.
##' @param matList list of matrices
##' @param n numeric vector of sample sizes
##' @param betaList list of eigenvectors (organized from largest to smallest eigenvalue) to be compared to7
##' @references Flury B.N.  Asymptotic theory for common principal component analysis.  Annals Statist.  1986.  14, 418-430.
##' @export
h_many=function(matList,n,betaList) {
  B=new_cpc(matList,n)
  k=length(matList)
  q=length(betaList)
  T=sum(n)
  counter=0
  part1=matrix(1,(q*(q-1)/2))
  for (l in 2:q) {
    for (j in 1:(l-1)) {
      x=matrix(1,k)
      counter=counter+1
      for (i in 1:k){
        D_hat=matrix(diag(t(B)%*%matList[[i]]%*%B))
        x[i]=(n[i]/T) *D_hat[j]-D_hat[l]/(D_hat[j]*D_hat[l])
      }
      g_hatinv=sum(x)
      part1[counter]=g_hatinv*(t(B[,l])%*%betaList[[j]]-t(B[,j])%*%betaList[[l]])^2
    }
  }
  
  counter=0
  part2=matrix(1,(q*(p-q)))
  for(l in (q+1):p) {
    for (j in 1:q) {
      counter=counter+1
      for (i in 1:k) {
        D_hat=matrix(diag(t(B)%*%matList[[i]]%*%B))
        x[i]=(n[i]/T) *D_hat[j]-D_hat[l]/(D_hat[j]*D_hat[l])
      }
      g_hatinv=sum(x)
      part2[counter]=g_hatinv*(t(B[,l])%*%betaList[[j]])^2
    }
  }
  chiSum=(T*(0.25*sum(part1)+sum(part2)))
  p.val <- pchisq(chiSum,df=q*(p-(q+1)/2),lower.tail=FALSE)
  datname <- deparse(substitute(matList))  ## recover _name_ of data argument
  vv <- list(data.name=datname,
             method="Flury test of common principal components",
             statistic=c("chi-sq"=chiSum),p.value=p.val)
  class(vv) <- "htest"
  vv
}




##' Comparing first CPC to one hypothetical eigenvector
##' 
##' Tests whether to reject the null hypothesis that the CPC with the 
##' largest eigenvalue is similar to a hypothetical eigenvector, based on 
##' equation 3.4 (Flury, 1986).  If comparing 2 or more eigenvectors,
##' see the function h_many.
##' @param matList list of matrices
##' @param n numeric vector of sample sizes
##' @param v eigenvector to be compared to
##' @references Flury B.N.  Asymptotic theory for common principal component analysis.  Annals Statist.  1986.  14, 418-430.
##' @examples source("test_examples.R")
##' v=rep((1/sqrt(3)),3)
##' h_one(test_turtle,n_turtle,v)
##' @export
h_one=function(matList,n,v) {
  B=new_cpc(matList,n)
  k=length(matList)
  p=nrow(matList[[1]])
  X=matrix(1,k)
  for (i in 1:k) {
    D=diag(diag(t(B)%*%matList[[i]]%*%B))
    sigma_hat=B%*%D%*%t(B)
    X[i]=n[i]%*%(D[1,1]%*%t(v)%*%solve(sigma_hat)%*%v+(1/D[1,1])%*%t(v)%*%sigma_hat%*%v-2)
  }
  chiSum=sum(X)
  p.val <- pchisq(chiSum,df=p-1,lower.tail=FALSE)
  datname <- deparse(substitute(matList))  ## recover _name_ of data argument
  vv <- list(data.name=datname,
             method="Flury test of common principal components",
             statistic=c("chi-sq"=chiSum),p.value=p.val)
  class(vv) <- "htest"
  vv
}




##' Reducing Parameter Space
##' 
##' Tests whether to reject the NH that the last p-q eigenvalues are very small
##' and approximately zero variance, based on approximation 4.3 (Flury, 1986).
##' @param matList list of covariance matrices
##' @param n numeric vector of sample sizes
##' @param q number of eigenvectors to be kept
##' @param f a fraction value between 0 and 1 which determines the relative
##' contribution of an eigenvalue to CPC
##' @param alpha confidence level
##' @references Flury B.N.  Asymptotic theory for common principal component analysis.  Annals Statist.  1986.  14, 418-430.
##' @export
h_reduce <- function(matList,n,q,f,alpha=0.05) {  
  p=nrow(matList[[1]])
  B=new_cpc(matList,n)
  k=length(matList)
  z=matrix(1,k)
  for (i in 1:k) {
    denom=1
    lambda_hat=diag(t(B)%*%matList[[i]]%*%B)
    D=diag(diag(t(B)%*%matList[[i]]%*%B))
    sigma_hat=B%*%D%*%t(B)
    c_hat=sum(D[1:q])
    d_hat=sum(diag(sigma_hat))-c_hat
    
    num=(sqrt(n[i])*((1-f)*d_hat-f*c_hat))
    denom1=(f^2)*sum((lambda_hat[1:q])^2)
    denom2=((1-f)^2)*sum((lambda_hat[(q+1):p])^2)
    z[i]=num/sqrt(2*(denom1+denom2))
  }
  beta=1-(1-alpha)^(1/k)
  z.beta <- pnorm(beta, mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE)
  datname <- deparse(substitute(matList))  ## recover _name_ of data argument
  vv <- list(data.name=datname,
             method="Flury test of common principal components: warning this
             is a comparison of the value to the z.beta value, NOT a pvalue!!!
             The output value should be greater than z.beta to reject the hypothesis
             I don't know how to change the title...",
             statistic=c("value"=max(z)),p.value=z.beta)
  class(vv) <- "htest"
  vv
  # z beta upper beta quantile of the standard normal distribution
}






##' Hypothesis of Sphericity
##' 
##' Tests whether to reject null hypothesis that the 
##' last p-q eigenvalues are equal, based on approximation 4.9 (Flury,1986).
##' The value p is the dimension of the matrices
##' and the value q is the number of eigenvectors to be kept.
##' @param matList list of covariance matrices
##' @param n numeric vector of sample sizes
##' @param q number of eigenvectors to be kept
##' @references Flury B.N.  Asymptotic theory for common principal component analysis.  Annals Statist.  1986.  14, 418-430.
##' @examples source("test_examples.R")
##' h_sphere(test_turtle,n_turtle,q=1)
##' @export
## approximation 4.9 determines whether to reject NH that the
## last p-q eigenvalues are equal (hypothesis of sphereficity)
h_sphere <- function(matList,n,q) {
  B=new_cpc(matList,n)
  k=length(matList)
  p=nrow(matList[[1]])
  chi=matrix(1,k)
  for (i in 1:k) {
    denom=1
    lambda_hat=matrix(diag(t(B)%*%matList[[i]]%*%B),p)
    lambda_star=(sum(lambda_hat[(q+1):p])/(p-q))
    for (a in (q+1):p) {
      denom=denom*lambda_hat[a]
    }
    chi[i]=n[i]*log((lambda_star^(p-q))/(denom))
  }
  chiSum=sum(chi)
  p.val <- pchisq(chiSum,df=(p-q-1)*(p-q+2*k)/2,lower.tail=FALSE)
  datname <- deparse(substitute(matList))  ## recover _name_ of data argument
  vv <- list(data.name=datname,
             method="Flury test of common principal components",
             statistic=c("chi-sq"=chiSum),p.value=p.val)
  class(vv) <- "htest"
  vv
}
