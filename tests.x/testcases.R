##testcase 
comp <- function(phil,mike){
  list(phileval=phil$evals, mikeeval = mike$evals, philevecs = phil$evecs,
       mikeevecs = mike$evecs)
}

fourbyfour <- function(var,n){
  phil <- newphil(covs=var,npts=n)
  mike1 <- fullcpc(var,n)
  mike2 <- fullcpc(var,n,2)
  mike3 <- fullcpc(var,n,1)
  comp1 <- comp(phil$datlist[[1]],mike1)
  comp2 <- comp(phil$datlist[[2]],mike2)
  comp3 <- comp(phil$datlist[[3]],mike3)
  return(list(comp1,comp2,comp3))
}

threebythree <- function(var,n){
  phil <- newphil(covs=var,npts=n)
  mike1 <- fullcpc(var,n)
  mike2 <- fullcpc(var,n,1)
  comp1 <- comp(phil$datlist[[1]],mike1)
  comp2 <- comp(phil$datlist[[2]],mike2)
  return(list(comp1,comp2))
}

voletest <- threebythree(test_vole,n_vole)
turtletest <- threebythree(test_turtle,n_turtle)

iristest <- fourbyfour(test_iris,n_iris)
martentest <- fourbyfour(test_marten,n_marten)
banktest <- fourbyfour(test_bank,n_bank)
