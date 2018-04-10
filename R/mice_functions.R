

getYso = function(N, w, dt, Ystar=0.4) Ystar*w*N*dt

getStartM = function(N, w, delta, dt, Ystar=0.4) {
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


VB = function(age, Linf, k, t0) Linf*(1-exp(-k*(age-t0)))





.getVar = function(pop, what) {
  out = list()
  for(i in seq_along(pop)) {
    out[[i]] = pop[[i]][[what]]
  }
  out = c(unlist(out))
  return(out)
}

.getSkeleton = function(pop, what=1) {
  out = list()
  for(i in seq_along(pop)) {
    out[[i]] = pop[[i]][[what]]
    out[[i]][] = NA
  }
  out = as.relistable(out)
  return(out)
}

getAccessibility2 = function(pop, predRange) {
  size = .getVar(pop, "size")
  .inside = function(x, y, limits) ((x >= y*limits[1]) & (x <= y*limits[2])) + 0
  out = outer(X = size, Y = size, limits=predRange, FUN=.inside)
  return(out)
}

getAccessibility = function(pop) {
  size = .getVar(pop, "size")
  psize = cbind(.getVar(pop, "psize_min"), .getVar(pop, "psize_max"))
  .inside = function(limits, x) ((x >= limits[1]) & (x <= limits[2])) + 0
  out = apply(psize, 1, FUN=.inside, x=size)
  return(out)
}


updateN = function(N, skeleton, plus=FALSE) {

  # skeleton[] = N
  N = relist(N, skeleton)

  .updateN = function(x, plus) {
    n = length(x)
    nTmp = c(0, head(x, -1))
    if(isTRUE(plus)) nTmp[n] = nTmp[n] + tail(x, 1)
    return(nTmp)
  }

  out = c(unlist(lapply(N, FUN=.updateN, plus=plus)))

  return(out)
}

updateR = function(N, skeleton, R) {

  # skeleton[] = N
  N = relist(N, skeleton)

  for(i in seq_along(N)) N[[i]][1] = R[i]

  out = c(unlist(N))

  return(out)
}


getRecruitment = function(SSB, SSN, skeleton, groups) {

  SSB = 1e-6*sapply(relist(SSB, skeleton), FUN=sum) # grams to tonnes
  SSN = sapply(relist(SSN, skeleton), FUN=sum) # total individuals, male + female

  R = numeric(length(groups))

  for(i in seq_along(groups)) {

    thisGroup = groups[[i]]
    recType = thisGroup$recType
    recBy   = thisGroup$recBy
    if(is.null(recType)) recType = Ricker
    if(is.null(recBy)) recBy = "biomass"
    recModel = match.fun(recType)
    x = switch(recBy, biomass=SSB[i], number=SSN[i])
    R[i] = recModel(x, par=thisGroup)
  }

  return(R)

}

Ricker = function(x, par) par$alpha*x*exp(-par$beta*x)

getBiomass = function(N, w, skeleton) {

  B = relist(N*w, skeleton)
  B = 1e-6*unlist(lapply(B, FUN=sum))
  return(B)
}
