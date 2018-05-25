

# Internal functions ------------------------------------------------------

calculateMortality = function(N, add, w, access, dt, Ystar=3.5, delta=0.9, niter=7) {

  eta_crit = 0.57
  Yso = getYso(N=N, w=w, dt=dt, Ystar=Ystar)
  M   = getStartM(N=N, w=w, delta=delta, dt=dt, Ystar=Ystar) # columns are predators
  Mj  = M*access # init

  ml = iterativeMortality(rMj=as.matrix(Mj), radd=as.numeric(add),
                          rw=as.numeric(w), rN=as.numeric(N),
                          rYso=as.numeric(Yso), rniter=as.integer(niter))

  ratio = pmax(1 - (ml$Yj/Yso)/eta_crit, 0)
  ratio[is.nan(ratio)] = 0

  return(list(M=ml$M, starv=c(ratio)))
}


getYso = function(N, w, dt, Ystar=3.5) Ystar*w*N*dt

getStartM = function(N, w, delta, dt, Ystar=3.5) {
  Ys = Ystar*w*N*dt
  startC = function(i, Ys, N, w, delta) pmin(Ys[i], delta*w*N)/w
  C = sapply(seq_along(N), FUN=startC, Ys=Ys, N=N, w=w, delta=delta)
  N = floor(N)
  C = floor(C)
  M = log(N/(N-C))
  M[is.nan(M)] = 0
  M[!is.finite(M)] = 0
  return(M)
}

#' @export
mortality.senecence.spec = function(x, par, tiny=1e-6) {

  A = par$A
  x = x/A

  M0  = if(is.null(par$M0))  3*log(10) else par$M0
  M50 = if(is.null(par$M50)) 0.89 else par$M50
  M75 = if(is.null(par$M75)) 0.91 else par$M75

  s1 = (M50*log(3))/(M75-M50)
  s2 = s1/M50
  selec = M0/(1+exp(s1-(s2*x)))
  selec[selec<tiny] = 0

  return(selec)

}
