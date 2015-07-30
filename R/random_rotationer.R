##' Random Rotation Matrix
##' 
##' Generates a random rotation matrix of dimension n, based on the subgroup
##' algorithm provided by Diaconis, P. and Shahshahani, M. (1987).
##' @param n dimension of matrix
##' @references Diaconis, P.; Shahshahani, M.  The subgroup algorithm for generating uniform random variables.  Probability in the Engineering and Informational Sciences.  1987, 1, 15-32.
##' @export
random_rotation_matrix <- function(n) {
  alpha=runif(1,min=-1,max=1)
  theta=asin(alpha)
  R=matrix(c(cos(theta),sin(theta),
             -sin(theta),cos(theta)),2)
  for (i in 3:n) {
    I=diag(i)
    I[2:i,2:i]=R
    ## pick a point at random on the i-sphere
    v=rnorm(i, mean = 0, sd = 1)
    d=sqrt(sum(v^2))
    newv=v/d
    newvector=(I[,1])-newv
    x=newvector/sqrt(t(newvector)%*%newvector)
    reflect=diag(i)-2*outer(x,x)
    R=reflect%*%I
  }
  R
}

##' Random Rotation
##' 
##' Applies a random rotation to a matrix based on the subgroup 
##' algorithm provided by Diaconis, P. and Shahshahani, M. (1987).
##' @param M matrix to be rotated
##' @references Diaconis, P.; Shahshahani, M.  The subgroup algorithm for generating uniform random variables.  Probability in the Engineering and Informational Sciences.  1987, 1, 15-32.
##' @export
random_rotation <- function(M) {
  R <- random_rotation_matrix(nrow(M))
  solve(R) %*% M %*% R
}