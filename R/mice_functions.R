



#' Run a simulation of the MICE model
#'
#' @param LHT List with the life history parameters
#' @param B0 Initial biomass for all species
#' @param predRange Minimum and maximum predator-prey ratios
#' @param Mstarv Maximum starvation mortality
#' @param Ystar Optimal annual food ration per gram of biomass
#' @param delta Fraction of prey population available to predators
#' @param dt Time step for simulation, in fractions of a year
#' @param T Time horizon for the simulation (years).
#'
#' @return A list with the abundance (N), length (L) and biomass (B) of
#' all the species modeled.
#' @export
#'
#' @examples
runMICE = function(LHT, B0, F=0, predRange = c(0.1, 0.25), Mstarv = 0.3,
                   Ystar = 0.4, delta = 0.9, dt = 1/4, T = 100/dt) {
  
  pop = initGroups(LHT=LHT, B0=B0, dt=dt)
  access = getAccessibility(pop=pop, predRange = predRange)
  
  skeleton = .getSkeleton(pop)
  
  N0 = .getVar(pop, "N")
  L0 = .getVar(pop, "size")
  ssb = .getVar(pop, "ssb")
  w = .getVar(pop, "w")
  
  N = matrix(nrow=length(N0), ncol=T+2)
  L = matrix(nrow=length(L0), ncol=T+2)
  B = matrix(NA, nrow=T+1, ncol=length(skeleton))
  
  N[, 1] = N0
  L[, 1] = L0
  B[1, ] = getBiomass(N[, 1], w, skeleton)
  
  for(i in 1:T) {
    mort = calculateMortality(N=N[, i], w=w, F=F*dt, access=access, dt)
    Z = rowSums(mort$M) + F*dt + Mstarv*mort$starv 
    N[, i+1] = updateN(N[, i]*exp(-Z), skeleton=skeleton, plus=FALSE)
    SSB = ssb*N[, i+1]
    R = getRecruitment(SSB, skeleton=skeleton, LHT) 
    N[, i+1] = updateR(N[, i+1], skeleton=skeleton, R)
    B[i+1, ] = getBiomass(N[, i+1], w, skeleton=skeleton)
  }
  
  return(list(N=N, L=L, B=B))
  
}



# Internal functions ------------------------------------------------------



calculateMortality = function(N, w, access, F, dt, Ystar=0.4, delta=0.9, niter=7) {
  
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


initGroups = function(LHT, B0, dt) {
  par = c(LHT, B0=list(B0))
  spp = seq_along(LHT$A)
  # out = do.call(rbind, lapply(spp, FUN=.initSpecies, dt=dt, par=par))
  out = lapply(spp, FUN=.initSpecies, dt=dt, par=par)
  return(out)
}


.initSpecies = function(i, dt, par) {
  A = par$A[i]
  Linf = par$Linf[i]
  k = par$k[i]
  t0 = par$t0[i]
  a = par$a[i]
  b = par$b[i]
  am = par$am[i]
  M0 = par$M0[i]
  B0 = par$B0[i]*1e6 # tonnes -> g
  M = if(is.null(M0)) -log(0.05)/A else M0
  out = data.frame(age=seq(0, A, by=dt))
  out$size = VB(age=out$age, Linf=Linf, k=k, t0=t0)
  out$w = a*out$size^b
  mw = exp(-M*out$age)*out$w
  out$N = B0*(mw/sum(mw))/out$w
  out$mature = (out$age >= am) + 0
  out$ssb = out$w*out$mature
  out$group = i
  return(out)
}


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

getAccessibility = function(pop, predRange) {
  size = .getVar(pop, "size")
  .inside = function(x, y, limits) ((x >= y*limits[1]) & (x <= y*limits[2])) + 0
  out = outer(X = size, Y = size, limits=predRange, FUN=.inside)
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


getRecruitment = function(SSB, skeleton, LHT) {
  
  # skeleton[] = SSB
  SSB = relist(SSB, skeleton)
  
  SSB = 1e-6*sapply(SSB, FUN=sum)
  
  spp = seq_along(SSB)
  
  iRicker = function(i, x, alpha, beta) Ricker(x[i], alpha[i], beta[i])
  R = sapply(spp, FUN=iRicker, x=SSB, alpha=LHT$alpha, beta=LHT$beta)
  return(R)
}

Ricker = function(x, alpha, beta) alpha*x*exp(-beta*x)

getBiomass = function(N, w, skeleton) {
  
  B = relist(N*w, skeleton)
  B = 1e-6*unlist(lapply(B, FUN=sum))
  return(B)
}
