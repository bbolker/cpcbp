
## testing pooled.cpc vs call.cpc accuracy (assume they
## get the relative signs correct)
set.seed(1001)
e1 =   eigen(covmat(c(10,10,10),cor=0.8))$vectors

fixesign = function(x,y) {
  sweep(x,2,sign(y[2,]/x[2,]),"*")
  ## condition on second row to avoid zero in first column!
}
  
tmpf = function(n) {
  X = simdata(npts=10000)
  c1 = fixesign(call.cpc(X),e1)
  c2 = fixesign(pooled.cpc(X),e1)
  de1 =  sum((abs(e1)-abs(c1))^2)
  de2 = sum((abs(e1)-abs(c2))^2)
  d12 = sum((abs(c1)-abs(c2))^2)
  c(de1=de1,de2=de2,d12=d12)

}
r = t(sapply(1:1000,function(n){
  if (n%%100==0) cat(n,"\n"); tmpf(n)
}))

par(mfrow=c(2,2))
invisible(apply(log10(r),2,hist))
invisible(apply(r,2,hist))

