

# Internal functions ------------------------------------------------------

calculateMortality = function(N, add, w, access, dt, Ystar=3.5, delta=0.9, niter=7) {

  eta_crit = 0.57
  Yso = getYso(N=N, w=w, dt=dt, Ystar=Ystar)
  M   = getStartM(N=N, w=w, delta=delta, dt=dt, Ystar=Ystar) # columns are predators
  Mj = M*access # init

  # YYj = matrix(ncol=ncol(Mj), nrow=niter)
  # cf = matrix(ncol=ncol(Mj), nrow=niter)

  for(i in 1:niter) {
    Zj = rowSums(Mj) + add # mortality for each "prey" (1-exp(-Z))*N is total deads of each prey
    Cj = (Mj/Zj)*(1-exp(-Zj))*N # total deads of each prey by predator
    Cj[is.nan(Cj)] = 0
    Yj = colSums(w*Cj)
    cfj = matrix(pmin(Yso/Yj, 1), nrow=length(Yj), ncol=length(Yj), byrow=TRUE)
    cfj[is.nan(cfj)] = 0
    Mj = cfj*Mj
    # YYj[i, ] = Yj
    # cf[i, ] = cfj[1,]
  }
  ratio = pmax(1 - (Yj/Yso)/eta_crit, 0)
  ratio[is.nan(ratio)] = 0

  return(list(M=Mj, starv=c(ratio)))
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
