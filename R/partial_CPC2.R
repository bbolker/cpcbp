##' Partial CPC Algorithm
##' 
##' Performs partial CPC Algorithm where the groups
##'  share 1 through p-2 CPCs, from Flury (1988) Book.  It will return
##'  k matrices in a form of a list, showing the common CPC to all matrices
##'  and also the uncommon CPC specific to each of the groups.  Please ensure
##'  that you order the columns of B=new_cpc(matList) according to
##'  the first q columns to be kept in this model.
##' @param matList list of matrices
##' @param n numeric vector of sample sizes
##' @param q the first q components to be kept
##' @param B CPC matrix
##' @param alpha_stop stopping value for small adjustments by jacobi rotations
##' @references Flury Book 1988
##' @examples partial_cpc(test_vole,n_vole,q=1,B=new_cpc(matList,n),alpha_stop=0.001)
##' @export
partial_cpc <- function(matList,n,q,B=NULL,alpha_stop=0.001) {
  debug <- FALSE
  
  if(q == ncol(matList[[1]]) - 1){
    return(full_cpc(matlist,n))}
  else{B=full_cpc(matList,n)[[1]]
  p <- nrow(matList[[1]])  ## matrix dimensions
  k <- length(matList)     ## number of matrices
  
  
  ##kinda looks like a quadratic form-ish
  quad <- function(X,B) t(B) %*% X %*% B
  
  ## Jacobi-rotate matrix B by an angle theta in the m-j plane
  jacobi <- function(B,theta,m,j) {
    H <- cbind(B[,m],B[,j])
    J=matrix(c(cos(theta),sin(theta),
               -sin(theta),cos(theta)),2)
    Hstar=H%*%J
    B[,m]=as.matrix(Hstar[,1])
    B[,j]=as.matrix(Hstar[,2])
    B
  }
  
  ##sweeper adjusts the uncommon parts of each matrix
  sweeper <- function(matList,B_star,q,n) {
    p <- nrow(B_star)
    k=length(matList)
    
    ##Partition the matrix B_star into the common and uncommon parts
    B1_star=B_star[,1:q,drop=FALSE]
    B2_star=B_star[,(q+1):p,drop=FALSE]
    spareList=matList
    for (i in 1:k) {
      spareList[[i]]=quad(matList[[i]],B1_star)
    }
    Q=new_cpc(spareList)
    B1_plus=B1_star%*%Q
    
    gvalue=matrix(0,k)
    for (i in 1:k) {
      spareList=matList
      B2_plus=B2_star%*%(eigen(quad(matList[[i]],B2_star))$vec)
      B_plus=cbind(B1_plus,B2_plus)
      gvalue[i]=n[i]%*%sum(log(diag(quad(matList[[i]],B_plus))))
    }
    sum(gvalue)
  }
  
  
  
  ##THE FUNCTION REALLY STARTS HERE!
  
  ##initial gvalue calculator
  gvalue=numeric(k)
  for (i in 1:k) {
    gvalue[i]=(n[i]%*%sum(log(diag(quad(matList[[i]],B)))))
  }
  
  ##
  ## probably equivalent to
  ldqFun <- function(m) sum(log(diag(quad(m,B))))
  n*sapply(matList,ldqFun)
  
  g_initial=sum(gvalue)
  
  
  alpha_start=0.1
  alpha=2*alpha_start
  num=0
  while (alpha>alpha_stop) {
    alpha=alpha/2
    if (debug) cat("alpha ->",alpha,"\n")
    
    counter=1
    while (counter>0) {
      counter=0
      for (m in 1:q) {
        for (j in (q+1):p) {
          if (debug) cat(m,j,"\n")
          B_initial=B
          
          if (debug) cat("rotate by alpha=",alpha,"\n") 
          ## rotate by alpha in one direction ...
          B_star=jacobi(B,alpha,m,j)
          g_plus=sweeper(matList,B_star,q,n)
          
          if (g_plus<g_initial) {
            if (debug) cat("g_plus=",g_plus,"<g_initial=",g_initial,"\n") 
            ## improvement
            B=B_star
            g_initial=g_plus
            counter=counter+1
          } else {
            if (debug) cat("rotate by -alpha=",alpha, "\n") 
            ## alpha rotation failed, try -alpha
            B_star=jacobi(B,-alpha,m,j)
            g_plus=sweeper(matList,B_star,q,n)
            
            if (g_plus<g_initial) {
              ## worked
              if (debug) cat("g_plus=",g_plus,"<g_initial=",g_initial,"\n") 
              B=B_star
              g_initial=g_plus
              counter=counter+1
            }
          }
        } ## loop over j
      }  ## loop over m
      num=num+1
    }
  }
  
  ##Taking the matrix B for a final run through to get specific B matrices
  ## which contain a part common to all matrices, and a part uncommon and specific
  ##to each matrix
  
  B1_star=B[,1:q,drop=FALSE]
  B2_star=B[,(q+1):p,drop=FALSE]
  spareList=matList
  for (i in 1:k) {
    spareList[[i]]=quad(matList[[i]],B1_star)
  }
  Q=new_cpc(spareList)
  B1_plus=B1_star%*%Q
  
  Blist=matList
  for (i in 1:k) {
    B2_plus=B2_star%*%(eigen(quad(matList[[i]],B2_star))$vec)
    B_plus=cbind(B1_plus,B2_plus)
    Blist[[i]]=B_plus
  }
  return(Blist)
} }
