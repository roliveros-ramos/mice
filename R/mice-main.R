#' Run a simulation of the MICE model
#'
#' @param groups A list containing the information to create the functional groups. See details.
#' @param fleets A list containing the information to create the fleets. See details.
#' @param environment A list containing the environmental information. One value per time step is required.
#' @param ndtPerYear Number of time steps per year.
#' @param Mstarv Maximum starvation mortality.
#' @param Ystar Optimal annual food ration per gram of biomass.
#' @param delta Fraction of prey population available to predators.
#' @param par A list with parameters values. These parameters take precedence over 'groups' and 'fleets'.
#' @param Fmult Fishing multiplier. All the fishing mortalities are multiplied by this value.
#' @param niter Number of iterations for the calculation of predation mortality.
#' @param verbose Logical, should running messages be produced?
#'
#' @return A list with the abundance (N), length (L) and biomass (B) of
#' all the species modeled.
#' @export
#'
#' @examples
runMICE = function(groups, fleets, environment=NULL, T, ndtPerYear=4,
                   Mstarv = 0.3, Ystar = 3.5, delta = 0.9, par=NULL,
                   Fmult=1, prices=NULL, niter=7, verbose=TRUE) {

  CPU.time = proc.time()

  ndt = ndtPerYear*T # total number of simulated time steps
  dt = 1/ndtPerYear # time step

  time = seq(from=0, to=T, length=ndt+1)

  groups = checkGroups(groups, ndt=ndt) # TO_DO: check, set defaults
  if(!is.null(par)) groups = updateParameters(groups, par)

  fleets = checkFleets(fleets) # TO_DO: check, set defaults
  if(!is.null(par)) fleets = updateParameters(fleets, par)

  if(length(Fmult)==1) Fmult = rep(Fmult, length(fleets))
  if(length(Fmult)!=length(fleets)) stop("Wrong number of Fishing multipliers.")

  environment = checkEnvironment(environment, ndt=ndt) # TO_DO: check, set defaults
  # if(!is.null(par)) environment = updateParameters(environment, par)

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

  Fmult = matrix(Fmult, nrow=length(L0), ncol=length(Fmult), byrow=TRUE)

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

    Fst = Fmult*getSelectivity(L=L[, t], fleets=fsh)*Ft[, , t]*dt
    F = rowSums(Fst)
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
    if(any(N[,t]<0)) {
      print(t)
    }
    SSB = w_ssb*N[, t+1]
    SSN = isMature*N[, t+1]
    R = getRecruitment(SSB=SSB, SSN=SSN, t=t, skeleton=skeleton,
                       groups=groups, environment=environment)
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

  tcb = rowSums(B[, types=="functional_group"])
  if(!is.null(prices)) {
    prix = matrix(prices, nrow=nrow(B), ncol=sum(types=="functional_group"), byrow=TRUE)
    value = rowSums(B[, types=="functional_group"]*prix)
  } else {
    value=NA
  }

  CPU.time = proc.time() - CPU.time

  if(isTRUE(verbose)) message(sprintf("\nModel run time: %0.1f seconds.", CPU.time[3]))

  raw = list(N=N, L=L, C=C, CatchbyFleet=CatchbyFleet,
             CatchbyFleetByGroup=CatchbyFleetByGroup)
  out = list(B=B, Y=Y, tcb=tcb, value=value, raw=raw, CPU.time=CPU.time, time=time, groupTypes=types)

  class(out) = "mice.model"

  return(out)

}


#' Calculate Prey-Predator sizes ratios based on a regression model
#'
#' @param predType The type of predation interaction (piscivorous, planktivorous, predacious).
#' @param Linf The maximum length of the species.
#' @param thr Lengths where a change in predation behaviour is seen.
#' @param smallest The smallest size to calculate the relationship. Default to 0.1 cm.
#'
#' @return A table with the prey-predator sizes ratios.
#' @export
#'
#' @examples
predPreyRatio = function(predType, Linf, thr=NULL, smallest=1e-1) {

  if(!is.null(thr)) {
    if(any(thr>=Linf)) stop("Threshold (thr) values must be lower than Linf.")
    if(any(thr<=smallest)) stop("Threshold (thr) values must be greater than 'smallest'.")
    if(any(diff(thr)<=0)) warning("Threshold values are not in order, sorting them...")
    thr = sort(thr)
  }

  validTypes = c("piscivorous", "planktivorous", "predacious")
  if(!all(predType %in% validType)) stop("Invalid predType.")

  nthr = c(smallest, thr, Linf)
  nstages = length(nthr) - 1

  msg = sprintf("You must provide %d values for 'predType'", nstages)
  if(length(predType)!=nstages) stop(msg)

  output = NULL

  for(i in seq_len(nstages)) {

    x = data.frame(PredLength=seq(from=nthr[i], to=nthr[i+1], by=0.1),
                   feedType=predType[i])
    x$logPredLength = log10(x$PredLength)
    x$pMin = 10^predict(preySizeModel$min, newdata=x)
    x$pMax = 10^predict(preySizeModel$max, newdata=x)
    modMin = lm(pMin ~ PredLength + 0, data=x)
    modMax = lm(pMax ~ PredLength + 0, data=x)

    out = c(1/coef(modMin), 1/coef(modMax))
    names(out) = c("predPrey.sizeRatio.min", "predPrey.sizeRatio.max")
    output = cbind(output, out)

  }

  return(output)

}



# Methods for main class --------------------------------------------------

#' @export
plot.mice.model = function(x, skip=0, ...) {

  sim = x
  nr = sum(sim$groupTypes == "resource")
  ns = ncol(sim$B) - nr
  ind = sim$time > skip
  par(mfrow=c(ns,1), mar=c(0,0,0,0), oma=c(3,4,1,4))
  for(i in seq_len(ns)) {
    icol =if(i!=7) i else 10
    plot(sim$time[ind], 1e-3*sim$B[ind,i+nr], col=icol, axes=FALSE, type="l")
    axis(2*((i+1)%%2+1), las=1); box()
    mtext(toupper(colnames(sim$B)[i+nr]), 3, adj=0.95, line=-1, cex=0.5, col=icol)
  }
  axis(1)
  return(invisible())

}


