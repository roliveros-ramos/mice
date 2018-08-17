

# Internal functions ------------------------------------------------------

calculateMortality = function(N, F, add, Mstarv, w, access, dt,
                              Ystar=3.5, delta=0.9, eta_crit = 0.57, niter=7L) {

  # Ystar is given in 1/year.

  if(length(Ystar)==1) Ystar = rep(Ystar, length(N))
  if(length(delta)==1) delta = rep(delta, length(N))
  # Yso = getYso(N=N, w=w, dt=dt, Ystar=Ystar)
  # Mj   = getStartM(N=N, w=w, delta=delta, dt=dt, Ystar=Ystar, access=access) # columns are predators
  # Mj  = M*access # init in 1/dt units!

  # iterative mortality assumes rates in the same units
  # ml = iterativeMortality(Mj=Mj, add=add*dt, w=w, N=N, Yso=Yso, niter=niter)
  ml = iterativeMortalityF(N=N, w=w, delta=delta, dt=dt, Ystar=Ystar, access=access,
                           F=F, add=add, Mstarv=Mstarv, eta_crit=eta_crit, niter=niter)

  # ratio = pmax(1 - (ml$Yj/Yso)/eta_crit, 0)
  # ratio[is.nan(ratio)] = 0

  # outputs in 1/year
  # return(list(M=ml$M/dt, starv=ratio))
  return(ml)
}

getYso = function(N, w, dt, Ystar=3.5) Ystar*w*N*dt


#' @export
mortality.senecence.spec = function(x, par, tiny=1e-6) {

  A = par$A
  x = x/A

  M0  = if(is.null(par$M0))  3*log(10) else par$M0
  M50 = if(is.null(par$M50)) 0.89 else par$M50
  M75 = if(is.null(par$M75)) (M50+0.02) else par$M75

  if(is.na(M50)) return(rep(0, length(x)))

  if(M75<M50) stop("M75 must be greater than M50")

  s1 = (M50*log(3))/(M75-M50)
  s2 = s1/M50
  selec = M0/(1+exp(s1-(s2*x)))
  selec[selec<tiny] = 0

  return(selec)

}

#' @export
mortality.constant.spec = function(x, par, tiny=1e-6) {

  Mb  = if(is.null(par$Mb))  1e-3 else par$Mb

  return(rep(Mb, length(x)))

}
