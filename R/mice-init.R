
initGroups = function(groups, dt) {
  out = lapply(groups, FUN=.initSpecies, dt=dt)
  return(out)
}

.initSpecies = function(par, dt) {
  A = par$A
  Linf = par$Linf
  k = par$k
  t0 = par$t0
  a = par$a
  b = par$b
  am = par$am
  M0 = par$M0
  B0 = par$B0*1e6 # tonnes -> g
  psmin = par$predRange[1] # to be updated with other criteria
  psmax = par$predRange[2]
  M = if(is.null(M0)) -log(0.05)/A else M0
  out = data.frame(age=seq(0, A, by=dt))
  out$size = VB(age=out$age, Linf=Linf, k=k, t0=t0)
  out$w = a*out$size^b
  mw = exp(-M*out$age)*out$w
  out$N = B0*(mw/sum(mw))/out$w
  out$mature = (out$age >= am) + 0
  out$ssb = out$w*out$mature
  out$psize_min = out$size*psmin
  out$psize_max = out$size*psmax
  # out$group = par$group
  out$name = par$name
  return(out)
}





