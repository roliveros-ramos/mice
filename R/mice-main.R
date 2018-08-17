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
runMICE = function(groups, fleets, environment=NULL, T, ndtPerYear=4,
                   par=NULL, Fmult=1, prices=NULL, eta_crit=0.57, niter=7,
                   subDtPerYear=24, verbose=TRUE) {

  CPU.time = proc.time()

  ndt  = ndtPerYear*T # total number of simulated time steps
  dt   = 1/ndtPerYear # time step
  xndt = ceiling(subDtPerYear/ndtPerYear)
  subDtPerYear = xndt*ndtPerYear
  xdt  = 1/subDtPerYear

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
  dietMatrix = getDietMatrix(pop=pop, groups=groups)

  skeleton = .getSkeleton(pop, what=2)

  isResource = which(.getVar(pop, "type") == "resource")

  N0 = .getVar(pop, "N")
  L0 = .getVar(pop, "size")
  dL = .getVar(pop, "dL")
  # length-weight
  a  = .getVar(pop, "a")
  b  = .getVar(pop, "b")

  egg_size = getEggSize(pop)
  egg_tl   = getEggTL(pop)

  Ystar = .getVar(pop, "Ystar")
  delta = .getVar(pop, "delta")

  w_ssb    = .getVar(pop, "ssb")
  isMature = .getVar(pop, "mature")
  w_mean   = .getVar(pop, "w_mean") # make it dynamic
  Madd     = .getVar(pop, "Mold") + .getVar(pop, "Mb")
  Mstarv   = .getVar(pop, "Mstarv")

  tl       = .getVar(pop, "TL") # initial TL for resources only

  No = matrix(nrow=length(N0), ncol=ndt+1)                   # abundance (start)
  Lo = matrix(nrow=length(L0), ncol=ndt+1)                   # length (start)

  N = matrix(nrow=length(N0), ncol=ndt+1)                   # abundance
  L = matrix(nrow=length(L0), ncol=ndt+1)                   # length

  C = matrix(nrow=length(N0), ncol=ndt)                     # catch (number)
  Bage = matrix(nrow=length(N0), ncol=ndt)                  # Biomass by age

  TL = matrix(nrow=length(N0), ncol=ndt)                    # Trophic level

  B = matrix(NA_real_, nrow=ndt, ncol=length(skeleton))   # biomass
  Y = matrix(NA_real_, nrow=ndt, ncol=length(skeleton))   # biomass

  CatchbyFleet        = array(dim=c(length(N0), length(fleets), ndt)) # catch (number)
  CatchbyFleetByGroup = array(dim=c(length(groups), length(fleets), ndt)) # catch (number)

  Fmult = matrix(Fmult, nrow=length(L0), ncol=length(Fmult), byrow=TRUE)

  No[, 1] = N0 # init numbers at the beginning of time step
  Lo[, 1] = L0 # init size at the beginning of the time step

  # B[1, ] = getBiomass(N[, 1], w, skeleton)

  access = array(dim = c(length(L0), length(L0), xndt))
  for(s in seq_len(xndt)) {
    psize = L0 + cbind((s-1)*dL, s*dL)/xndt
    access[,, s] = getAccessibility(size=psize, pop=pop, dietMatrix=dietMatrix)
  }

  if(isTRUE(verbose)) pb = txtProgressBar(style=3)

  for(t in seq_len(ndt)) {

    if(isTRUE(verbose)) setTxtProgressBar(pb, (t-1)/ndt)

    Nx = matrix(nrow=length(N0), ncol=xndt+1) # abundance for subdt
    Lx = matrix(nrow=length(L0), ncol=xndt+1) # length for subdt
    Bx = matrix(nrow=length(L0), ncol=xndt+1) # length for subdt

    Cx = matrix(nrow=length(N0), ncol=xndt)
    Yx = matrix(nrow=length(N0), ncol=xndt)

    preyed = matrix(0, nrow=length(N0), ncol=length(N0))

    Nx[, 1] = No[, t]
    Lx[, 1] = Lo[, t]
    Bx[, 1] = No[, t]*a*Lo[, t]^b

    for(s in seq_len(xndt)) { # for predation and growth

      Fst = Fmult*getSelectivity(L=Lx[, s], fleets=fsh)*Ft[, , t]
      F   = rowSums(Fst) # annual rate

      size05 = Lx[, s] + 0.5*dL/xndt # mid-step length
      w05    = a*size05^b # mid-step weight
      mort = calculateMortality(N=Nx[, s], F=F, add=Madd, Mstarv=Mstarv, w=w05,
                                access=access[, , s],
                                dt=xdt, Ystar=Ystar, delta=delta,
                                eta_crit=eta_crit, niter=niter)
      # M  = rowSums(mort$M)
      # Ms = Mstarv*mort$starv
      # Z  = M + Ms + Madd + F
      # deaths = Nx[, s]*(1-exp(-Z*xdt))

      Z      = mort$Z
      deaths = mort$deaths

      # preyed = preyed + ((w05*deaths/Z)*mort$M)
      preyed = preyed + mort$preyed

      Nx[, s+1] = Nx[, s] - deaths # abundance at the end of sub-dt
      Lx[, s+1] = Lx[, s] + dL/xndt # length at the end of sub-dt
      Bx[, s+1] = Nx[, s+1]*a*Lx[, s+1]^b # biomass at the end of sub-dt

      Cx[, s]   = deaths*F/Z # catch during sub-dt
      Yx[, s]   = Cx[, s]*w05 # catch during sub-dt

      # CatchbyFleet[, , t] = deaths*Fst/Z
      # CatchbyFleetByGroup[, , t] = rowsum(CatchbyFleet[, , t], group=.getVar(pop, "name"), reorder=FALSE)

    }

    TL[, t] = calculateTL(preyed, tl, isResource, t)
    tl      = updateTL(TL=TL[, t], skeleton=skeleton, egg_tl=egg_tl)

    C[, t] = rowSums(Cx)

    bio = getBiomass(Bx, skeleton=skeleton)
    Bage[, t] = bio$Bage

    B[t, ] = bio$B  # biomass by functional group
    Y[t, ] = 1e-6*rowsum(rowSums(Yx), group=.getVar(pop, "name"), reorder = FALSE)

    # reproduction
    SSB = Nx[, xndt]*w_ssb     # abundance at the end (instantaneous reproduction)
    SSN = Nx[, xndt]*isMature
    R   = getRecruitment(SSB=SSB, SSN=SSN, t=t, skeleton=skeleton, groups=groups,
                         environment=environment)

    # save state variables
    No[, t+1] = updateN(Nx[, xndt], skeleton=skeleton, recruits=R, plus=FALSE) # reproduction
    Lo[, t+1] = updateL(Lx[, xndt], skeleton=skeleton, egg_size=egg_size)

  }

  if(isTRUE(verbose)) {
    pb = txtProgressBar(style=3)
    setTxtProgressBar(pb, (t-1)/ndt)
  }

  # post-processing outputs


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

  # temp
  N = No
  L = Lo

  raw = list(N=N, L=L, C=C, TL=TL, Bage=Bage, CatchbyFleet=CatchbyFleet,
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
#' @examples predPreyRatio(predType = c("planktivorous", "predacious"), Linf = 30, thr = 15)
predPreyRatio = function(predType, Linf, thr=NULL, smallest=1e-1) {

  if(!is.null(thr)) {
    if(any(thr>=Linf)) stop("Threshold (thr) values must be lower than Linf.")
    if(any(thr<=smallest)) stop("Threshold (thr) values must be greater than 'smallest'.")
    if(any(diff(thr)<=0)) warning("Threshold values are not in order, sorting them...")
    thr = sort(thr)
  }

  validTypes = c("piscivorous", "planktivorous", "predacious")
  if(!all(predType %in% validTypes)) stop("Invalid predType.")

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

  colnames(output) = levels(cut(seq(0, Linf), breaks=c(0, thr, Linf)))

  return(output)

}


# Methods for main class --------------------------------------------------




