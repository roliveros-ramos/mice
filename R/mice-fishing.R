
getSelectivity = function(L, fleets) {

  out = matrix(0, nrow=length(L), ncol=length(fleets))

  for(i in seq_along(fleets)) {
    out[, i] = fleets[[i]]$selectivity(L)
  }

  return(out)

}

getFishingMortality = function(fsh, pop) {

  grs = .getVar(pop, "name")
  ndt = length(fsh[[1]]$F)

  f = lapply(sapply(fsh, FUN="[", i="F"), FUN=rep, each=length(grs))
  f = array(do.call(c, f), dim=c(length(grs), ndt, length(fsh)))
  f = aperm(f, perm = c(1,3,2)) # group x fishery x time
  colnames(f) = names(fsh)
  accessMatrix = array(sapply(sapply(fsh, FUN="[", i="target"),
                              FUN="[", i=grs), dim=dim(f))
  f = f*accessMatrix
  return(f)

}


# Selectivity models ------------------------------------------------------

# TO_DO: set to zero for L<0
# TO_DO: parameter check

#' @export
selectivity.logistic.spec = function(x, par, tiny=1e-6) {

  L50 = par$L50
  L75 = par$L75

  s1 = (L50*log(3))/(L75-L50)
  s2 = s1/L50
  selec = 1/(1+exp(s1-(s2*x)))
  selec[selec<tiny] = 0

  return(selec)

}

#' @export
selectivity.gaussian.spec = function(x, par, tiny=1e-5) {

  L50 = par$L50
  L75 = par$L75

  sd = (L75-L50)/qnorm(0.75)
  mean = L50
  selec = dnorm(x, mean=mean, sd=sd)
  selec = selec/max(selec, na.rm=TRUE)
  selec[selec<tiny] = 0

  return(selec)

}

#' @export
selectivity.uniform.spec = function(x, par, tiny=1e-5) {

  selec = rep(1, length(x))

  return(selec)

}

#' @export
'selectivity.knife-edge.spec' = function(x, par, tiny=1e-5) {

  L50 = par$L50

  selec = numeric(length(x))
  selec[x >= L50] = 1

  return(selec)

}
