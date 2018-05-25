#' @export
runMICE2 = function(groups, fleets, T, ndtPerYear=4,
                   Mstarv = 0.3, Ystar = 3.5, delta = 0.9, par=NULL,
                   Fmult=1, niter=7, verbose=TRUE) {

  CPU.time = proc.time()

  ndt = ndtPerYear*T # total number of simulated time steps
  dt = 1/ndtPerYear # time step

  time = seq(from=0, to=T, length=ndt+1)

  groups = checkGroups(groups, ndt=ndt) # TO_DO: check, set defaults
  if(!is.null(par)) groups = updateParameters(groups, par)

  fleets = checkFleets(fleets) # TO_DO: check, set defaults
  if(!is.null(par)) fleets = updateParameters(fleets, par)

  types = sapply(groups, FUN="[", i="type")

  pop = initGroups(groups=groups, dt=dt)

  fsh = initFleets(fleets=fleets, groups=groups, ndtPerYear=ndtPerYear, T=T) # on annual scale or already on dt? F$total?


  Ft = getFishingMortality(fsh, pop)

  access = getAccessibility(pop=pop, groups=groups)

  skeleton = .getSkeleton(pop, what=2)

  N0 = .getVar(pop, "N")
  L0 = .getVar(pop, "size")
  w_ssb = .getVar(pop, "ssb")
  isMature = .getVar(pop, "mature")
  w = .getVar(pop, "w_mean") # make it dynamic
  Mold = .getVar(pop, "Mold")*dt

  N = matrix(nrow=length(N0), ncol=ndt+1) # abundance
  C = matrix(nrow=length(N0), ncol=ndt) # catch (number)
  L = matrix(nrow=length(L0), ncol=ndt+1) # length
  B = matrix(NA, nrow=ndt+1, ncol=length(skeleton))

  CatchbyFleet = array(dim=c(length(N0), length(fleets), ndt)) # catch (number)
  CatchbyFleetByGroup = array(dim=c(ndt, length(groups), length(fleets))) # catch (number)

  N[, 1] = N0
  L[, 1] = L0
  B[1, ] = getBiomass(N[, 1], w, skeleton)

  for(t in seq_len(ndt)) {

    if(isTRUE(verbose)) {
      pb = txtProgressBar(style=3)
      setTxtProgressBar(pb, (t-1)/ndt)
    }

    Fst = getSelectivity(L=L[, t], fleets=fsh)*Ft[, , t]*dt
    F = rowSums(Fst)*Fmult
    mort = calculateMortality2(N=N[, t], add=F+Mold, w=w, access=access,
                              dt=dt, Ystar=Ystar, delta=delta, niter=niter)
    M  = rowSums(mort$M)
    Ms = Mstarv*mort$starv
    Z = M + Ms + Mold + F
    deaths = N[, t]*(1-exp(-Z))
    C[, t]   = updateC(deaths*F/Z)
    CatchbyFleet[, , t] = deaths*Fst/Z
    CatchbyFleetByGroup[t,,] = rowsum(CatchbyFleet[, , t],
                                      group=.getVar(pop, "name"))
    N[, t+1] = updateN(N[, t]*exp(-Z), skeleton=skeleton, plus=FALSE)
    SSB = w_ssb*N[, t+1]
    SSN = isMature*N[, t+1]
    R = getRecruitment(SSB=SSB, SSN=SSN, t=t, skeleton=skeleton, groups=groups)
    N[, t+1] = updateR(N[, t+1], skeleton=skeleton, R)
    L[, t+1] = updateL(L0)
    B[t+1, ] = getBiomass(N[, t+1], w, skeleton=skeleton)

  }

  if(isTRUE(verbose)) {
    pb = txtProgressBar(style=3)
    setTxtProgressBar(pb, (t-1)/ndt)
  }

  # post-processing outputs

  Y = t(rowsum(1e-6*C*w, group=.getVar(pop, "name"))) # grams to tonnes

  colnames(Y) = colnames(B) = sapply(groups, FUN="[[", i="name")

  CPU.time = proc.time() - CPU.time

  if(isTRUE(verbose)) message(sprintf("\nModel run time: %0.1f seconds.", CPU.time[3]))

  raw = list(N=N, L=L, C=C, CatchbyFleet=CatchbyFleet,
             CatchbyFleetByGroup=CatchbyFleetByGroup)
  out = list(B=B, Y=Y, raw=raw, CPU.time=CPU.time, time=time, groupTypes=types)

  class(out) = "mice.model"

  return(out)

}



calculateMortality2 = function(N, add, w, access, dt, Ystar=3.5, delta=0.9, niter=7) {

  eta_crit = 0.57
  Yso = getYso(N=N, w=w, dt=dt, Ystar=Ystar)
  M   = getStartM(N=N, w=w, delta=delta, dt=dt, Ystar=Ystar) # columns are predators
  Mj = M*access # init

  ml = iterativeMortality2(Mj=Mj, add=add, w=w, N=N, Yso=Yso, niter=niter)


  ratio = pmax(1 - (ml$Yj/Yso)/eta_crit, 0)
  ratio[is.nan(ratio)] = 0

  return(list(M=ml$M, starv=c(ratio)))
}



iterativeMortality2 = function(Mj, add, w, N, Yso, niter) {
  for(k in 1:niter) {

    # Zj = rowSums(Mj) + add # mortality for each "prey" (1-exp(-Z))*N is total deads of each prey
    Zj = numeric(nrow(Mj))
    Dj = numeric(nrow(Mj))
    for(i in 1:nrow(Mj)) {
      for(j in 1:ncol(Mj)) {
        Zj[i] = Zj[i] + Mj[i,j] # Zj[] == 0
      }
      Zj[i] = Zj[i] + add[i]
      Dj[i] = 1 - exp(-Zj[i])
    }

    # Cj = (Mj/Zj)*(1-exp(-Zj))*N # total deads of each prey by predator
    # Cj[is.nan(Cj)] = 0
    # Yj = colSums(w*Cj)
    # cfj = matrix(pmin(Yso/Yj, 1), nrow=length(Yj), ncol=length(Yj), byrow=TRUE)
    # cfj[is.nan(cfj)] = 0
    Cj = matrix(0, ncol=ncol(Mj), nrow=nrow(Mj))
    Yj = numeric(ncol(Mj)) # init on 0
    cfj = numeric(ncol(Mj))
    cfj[] = 1
    ans = 0
    for(i in 1:nrow(Mj)) {
      for(j in 1:ncol(Mj)) {
        ans = (Mj[i,j]/Zj[i])*Dj[i]*N[i]
        if(is.nan(ans)) {
          Cj[i,j] = 0
        } else {
          Cj[i,j] = ans
        }
        Yj[j] = Yj[j] + w[i]*Cj[i,j]
      }
    }

    cfj_tmp = 0
    for(j in 1:ncol(Mj)) {
      cfj_tmp = Yso[j]/Yj[j]
      if(is.nan(cfj_tmp)) cfj_tmp = 0
      if(cfj_tmp < 1) cfj[j] = cfj_tmp
    }

    # Mj = cfj*Mj
    for(i in 1:nrow(Mj)) {
      for(j in 1:ncol(Mj)) {
        Mj[i,j] = cfj[j]*Mj[i,j]
      }
    }


  }

  return(list(M=Mj, Yj=Yj))
}


# RCPP version ------------------------------------------------------------

#' @export
runMICE3 = function(groups, fleets, T, ndtPerYear=4,
                    Mstarv = 0.3, Ystar = 3.5, delta = 0.9, par=NULL,
                    Fmult=1, niter=7, verbose=TRUE) {

  CPU.time = proc.time()

  ndt = ndtPerYear*T # total number of simulated time steps
  dt = 1/ndtPerYear # time step

  time = seq(from=0, to=T, length=ndt+1)

  groups = checkGroups(groups, ndt=ndt) # TO_DO: check, set defaults
  if(!is.null(par)) groups = updateParameters(groups, par)

  fleets = checkFleets(fleets) # TO_DO: check, set defaults
  if(!is.null(par)) fleets = updateParameters(fleets, par)

  types = sapply(groups, FUN="[", i="type")

  pop = initGroups(groups=groups, dt=dt)

  fsh = initFleets(fleets=fleets, groups=groups, ndtPerYear=ndtPerYear, T=T) # on annual scale or already on dt? F$total?


  Ft = getFishingMortality(fsh, pop)

  access = getAccessibility(pop=pop, groups=groups)

  skeleton = .getSkeleton(pop, what=2)

  N0 = .getVar(pop, "N")
  L0 = .getVar(pop, "size")
  w_ssb = .getVar(pop, "ssb")
  isMature = .getVar(pop, "mature")
  w = .getVar(pop, "w_mean") # make it dynamic
  Mold = .getVar(pop, "Mold")*dt

  N = matrix(nrow=length(N0), ncol=ndt+1) # abundance
  C = matrix(nrow=length(N0), ncol=ndt) # catch (number)
  L = matrix(nrow=length(L0), ncol=ndt+1) # length
  B = matrix(NA, nrow=ndt+1, ncol=length(skeleton))

  CatchbyFleet = array(dim=c(length(N0), length(fleets), ndt)) # catch (number)
  CatchbyFleetByGroup = array(dim=c(ndt, length(groups), length(fleets))) # catch (number)

  N[, 1] = N0
  L[, 1] = L0
  B[1, ] = getBiomass(N[, 1], w, skeleton)

  for(t in seq_len(ndt)) {

    if(isTRUE(verbose)) {
      pb = txtProgressBar(style=3)
      setTxtProgressBar(pb, (t-1)/ndt)
    }

    Fst = getSelectivity(L=L[, t], fleets=fsh)*Ft[, , t]*dt
    F = rowSums(Fst)*Fmult
    mort = calculateMortality3(N=N[, t], add=F+Mold, w=w, access=access,
                               dt=dt, Ystar=Ystar, delta=delta, niter=niter)
    M  = rowSums(mort$M)
    Ms = Mstarv*mort$starv
    Z = M + Ms + Mold + F
    deaths = N[, t]*(1-exp(-Z))
    C[, t]   = updateC(deaths*F/Z)
    CatchbyFleet[, , t] = deaths*Fst/Z
    CatchbyFleetByGroup[t,,] = rowsum(CatchbyFleet[, , t],
                                      group=.getVar(pop, "name"))
    N[, t+1] = updateN(N[, t]*exp(-Z), skeleton=skeleton, plus=FALSE)
    SSB = w_ssb*N[, t+1]
    SSN = isMature*N[, t+1]
    R = getRecruitment(SSB=SSB, SSN=SSN, t=t, skeleton=skeleton, groups=groups)
    N[, t+1] = updateR(N[, t+1], skeleton=skeleton, R)
    L[, t+1] = updateL(L0)
    B[t+1, ] = getBiomass(N[, t+1], w, skeleton=skeleton)

  }

  if(isTRUE(verbose)) {
    pb = txtProgressBar(style=3)
    setTxtProgressBar(pb, (t-1)/ndt)
  }

  # post-processing outputs

  Y = t(rowsum(1e-6*C*w, group=.getVar(pop, "name"))) # grams to tonnes

  colnames(Y) = colnames(B) = sapply(groups, FUN="[[", i="name")

  CPU.time = proc.time() - CPU.time

  if(isTRUE(verbose)) message(sprintf("\nModel run time: %0.1f seconds.", CPU.time[3]))

  raw = list(N=N, L=L, C=C, CatchbyFleet=CatchbyFleet,
             CatchbyFleetByGroup=CatchbyFleetByGroup)
  out = list(B=B, Y=Y, raw=raw, CPU.time=CPU.time, time=time, groupTypes=types)

  class(out) = "mice.model"

  return(out)

}


# original one
calculateMortality3 = function(N, add, w, access, dt, Ystar=3.5, delta=0.9, niter=7) {

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
