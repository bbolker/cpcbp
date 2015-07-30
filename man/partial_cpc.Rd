\name{partial_cpc}
\alias{partial_cpc}
\title{Partial CPC Algorithm}
\usage{
  partial_cpc(matList, n, q, B = new_cpc(matList, n),
    alpha_stop = 0.001)
}
\arguments{
  \item{matList}{list of matrices}

  \item{n}{numeric vector of sample sizes}

  \item{q}{the first q components to be kept}

  \item{B}{CPC matrix}

  \item{alpha_stop}{stopping value for small adjustments by
  jacobi rotations}
}
\description{
  Performs partial CPC Algorithm where the groups share 1
  through p-2 CPCs, from Flury (1988) Book.  It will return
  k matrices in a form of a list, showing the common CPC to
  all matrices and also the uncommon CPC specific to each
  of the groups.  Please ensure that you order the columns
  of B=new_cpc(matList) according to the first q columns to
  be kept in this model.
}
\examples{
partial_cpc(test_vole,n_vole,q=1,B=new_cpc(matList,n),alpha_stop=0.001)
}
\references{
  Flury Book 1988
}

