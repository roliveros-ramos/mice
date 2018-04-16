

# Internal functions ------------------------------------------------------

calculateMortality = function(N, F, w, access, dt, Ystar=0.4, delta=0.9, niter=7) {

  Yso = getYso(N=N, w=w, dt=dt, Ystar=Ystar)
  M   = getStartM(N=N, w=w, delta=delta, dt=dt, Ystar=Ystar) # columns are predators
  Mj = M*access # init

  # YYj = matrix(ncol=ncol(Mj), nrow=niter)
  # cf = matrix(ncol=ncol(Mj), nrow=niter)

  for(i in 1:niter) {
    Zj = rowSums(Mj) + F # mortality for each "prey" (1-exp(-Z))*N is total deads of each prey
    Cj = (Mj/Zj)*(1-exp(-Zj))*N # total deads of each prey by predator
    Cj[is.nan(Cj)] = 0
    Yj = colSums(w*Cj)
    cfj = matrix(pmin(Yso/Yj, 1), nrow=length(Yj), ncol=length(Yj), byrow=TRUE)
    cfj[is.nan(cfj)] = 0
    Mj = cfj*Mj
    # YYj[i, ] = Yj
    # cf[i, ] = cfj[1,]
  }
  ratio = pmax(1-Yj/Yso, 0)
  return(list(M=Mj, starv=c(ratio)))
}
