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
runMICE = function(groups, fleets, T, ndtPerYear=4,
                   Mstarv = 0.3, Ystar = 3.5, delta = 0.9, Fmult=1, niter=7,
                   verbose=TRUE) {

  CPU.time = proc.time()

  ndt = ndtPerYear*T # total number of simulated time steps
  dt = 1/ndtPerYear # time step

  time = seq(from=0, to=T, length=ndt+1)

  groups = checkGroups(groups, ndt=ndt) # TO_DO: check, set defaults
  fleets = checkFleets(fleets) # TO_DO: check, set defaults

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

    Fst = getSelectivity(L=L[, t], fleets=fsh)*Ft[, , t]*dt
    F = rowSums(Fst)*Fmult
    mort = calculateMortality(N=N[, t], add=F+Mold, w=w, access=access,
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

  # post-processing outputs

  Y = t(rowsum(1e-6*C*w, group=.getVar(pop, "name"))) # grams to tonnes

  colnames(Y) = colnames(B) = sapply(groups, FUN="[[", i="name")

  CPU.time = proc.time() - CPU.time

  if(isTRUE(verbose)) message(sprintf("Model run time: %0.1f seconds.", CPU.time[3]))

  raw = list(N=N, L=L, C=C, CatchbyFleet=CatchbyFleet,
             CatchbyFleetByGroup=CatchbyFleetByGroup)
  out = list(B=B, Y=Y, raw=raw, CPU.time=CPU.time, time=time, groupTypes=types)

  class(out) = "mice.model"

  return(out)

}



# Methods for main class --------------------------------------------------

#' @export
plot.mice.model = function(x, what="B", ...) {

  xtime = if(what=="Y") head(x$time, -1) else x$time
  iPlot = which(x$groupTypes == "functional_group")
  matplot(xtime, x[[what]][, iPlot], lty=1, type="l", las=1, lwd=2, xlab="", ylab="", ...)
  title(main = "Biomass by functional group")
  return(invisible())

}


