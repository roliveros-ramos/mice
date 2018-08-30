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
                   subDtPerYear=24, verbose=TRUE, pop=NULL, resources=NULL) {

  CPU.time = proc.time()

  # time step calculations
  ndt  = ndtPerYear*T # total number of simulated time steps
  dt   = 1/ndtPerYear # time step
  xndt = ceiling(subDtPerYear/ndtPerYear)
  subDtPerYear = xndt*ndtPerYear
  xdt  = 1/subDtPerYear

  time = seq(from=0, to=T, length=ndt+1)

  groups = checkGroups(groups, ndtPerYear=ndtPerYear,
                       subDtPerYear=subDtPerYear, T=T, resources=resources)
  # TO_DO: check, set defaults
  if(!is.null(par)) groups = updateParameters(groups, par)

  fleets = checkFleets(fleets) # TO_DO: check, set defaults
  if(!is.null(par)) fleets = updateParameters(fleets, par)

  if(length(Fmult)==1) Fmult = rep(Fmult, length(fleets))
  if(length(Fmult)!=length(fleets)) stop("Wrong number of Fishing multipliers.")

  environment = checkEnvironment(environment, ndt=ndt) # TO_DO: check, set defaults
  # if(!is.null(par)) environment = updateParameters(environment, par)

  types = sapply(groups, FUN="[", i="type")
  isResourceGroup = (types == "resource")

  pop = if(is.null(pop)) initGroups(groups=groups, dt=dt) else pop
  fsh = initFleets(fleets=fleets, groups=groups, ndtPerYear=ndtPerYear, T=T) # on annual scale or already on dt? F$total?
  resources = getResources(groups=groups)

  Ft = getFishingMortality(fsh, pop)
  dietMatrix = getDietMatrix(pop=pop, groups=groups)

  skeleton = .getSkeleton(pop, what=2)

  isResource = which(.getVar(pop, "type") == "resource")

  gNames = .getVar(pop, "name")

  N0 = .getVar(pop, "N")
  L0 = .getVar(pop, "size")
  dL = .getVar(pop, "dL")
  # length-weight
  a  = .getVar(pop, "a")
  b  = .getVar(pop, "b")

  egg_size = getEggSize(pop)
  egg_tl   = getEggTL(pop)

  Ystar = .getVar(pop, "Ystar")
  delta = getDeltaBar(.getVar(pop, "delta"), xndt)

  w_ssb    = .getVar(pop, "ssb")
  isMature = .getVar(pop, "mature")
  w_mean   = .getVar(pop, "w_mean") # make it dynamic
  Madd     = .getVar(pop, "Mold") + .getVar(pop, "Mb")
  Mstarv   = .getVar(pop, "Mstarv")

  tl       = .getVar(pop, "TL") # initial TL for resources only

  nGroups = length(skeleton)
  nSizeGroups = length(N0)

  No = matrix(nrow=nSizeGroups, ncol=ndt+1)                  # abundance (start)
  Lo = matrix(nrow=nSizeGroups, ncol=ndt+1)                  # length (start)

  N = matrix(nrow=nSizeGroups, ncol=ndt+1)                   # abundance
  L = matrix(nrow=nSizeGroups, ncol=ndt+1)                   # length

  C    = matrix(nrow=nSizeGroups, ncol=ndt)                  # catch (number)
  Bage = matrix(nrow=nSizeGroups, ncol=ndt)                  # Biomass by age

  TL = matrix(nrow=nSizeGroups, ncol=ndt)                    # Trophic level

  B = matrix(NA_real_, nrow=ndt, ncol=nGroups)     # biomass
  Y = matrix(NA_real_, nrow=ndt, ncol=nGroups)     # yield

  CatchbyFleet        = array(dim=c(nSizeGroups, length(fleets), ndt)) # catch (number)
  YieldbyFleet        = array(dim=c(nSizeGroups, length(fleets), ndt)) # catch (number)
  YieldbyFleetByGroup = array(dim=c(nGroups,     length(fleets), ndt)) # catch (number)

  Fmult = matrix(Fmult, nrow=nSizeGroups, ncol=length(Fmult), byrow=TRUE)

  No[, 1] = N0 # init numbers at the beginning of time step
  Lo[, 1] = L0 # init size at the beginning of the time step

  # B[1, ] = getBiomass(N[, 1], w, skeleton)

  access = array(dim = c(nSizeGroups, nSizeGroups, xndt))
  for(s in seq_len(xndt)) {
    psize = L0 + cbind((s-1)*dL, s*dL)/xndt
    access[,, s] = getAccessibility(size=psize, pop=pop, dietMatrix=dietMatrix)
  }

  if(isTRUE(verbose)) pb = txtProgressBar(style=3)

  for(t in seq_len(ndt)) {

    if(isTRUE(verbose)) setTxtProgressBar(pb, (t-1)/ndt)

    Nx = matrix(nrow=nSizeGroups, ncol=xndt+1) # abundance for subdt
    Lx = matrix(nrow=nSizeGroups, ncol=xndt+1) # length for subdt
    Bx = matrix(nrow=nSizeGroups, ncol=xndt+1) # length for subdt

    Cx = matrix(nrow=nSizeGroups, ncol=xndt)
    Yx = matrix(nrow=nSizeGroups, ncol=xndt)

    CatchbyFleetx        = array(dim=c(nSizeGroups, length(fleets), xndt)) # catch (number)
    CatchbyFleetByGroupx = array(dim=c(nGroups,     length(fleets), xndt)) # catch (number)

    preyed = matrix(0, nrow=nSizeGroups, ncol=nSizeGroups)

    Nx[, 1] = No[, t]
    Lx[, 1] = Lo[, t]
    Bx[, 1] = No[, t]*a*Lo[, t]^b

    for(s in seq_len(xndt)) { # for predation and growth

      step = (t-1)*xndt + s

      size05 = Lx[, s] + 0.5*dL/xndt  # mid-step length
      w05    = a*size05^b             # mid-step weight, improve formula

      Fst = Fmult*getSelectivity(L=size05, fleets=fsh)*Ft[, , t]
      F   = rowSums(Fst) # annual rate

      mort = calculateMortality(N=Nx[, s], F=F, add=Madd, Mstarv=Mstarv, w=w05,
                                access=access[, , s], dt=xdt,
                                Ystar=Ystar, delta=delta, eta_crit=eta_crit,
                                niter=niter)

      Z      = mort$Z
      deaths = mort$deaths
      preyed = preyed + mort$preyed

      # abundance at the end of sub-dt
      Nx[, s+1] = Nx[, s] - deaths
      Nx[, s+1] = updateResources(Nx[, s+1], resources[, step], isResource)
      # length at the end of sub-dt
      Lx[, s+1] = Lx[, s] + dL/xndt
      # biomass at the end of sub-dt
      Bx[, s+1] = Nx[, s+1]*a*Lx[, s+1]^b
      # catch during sub-dt
      Cx[, s]   = deaths*F/Z
      # yield during sub-dt
      Yx[, s]   = Cx[, s]*w05

      CatchbyFleetx[, , s] = deaths*Fst/Z
      YieldbyFleetx[, , s] = CatchbyFleetx[, , s]*w05
      YieldbyFleetByGroupx[, , s] = rowsum(YieldbyFleetx[, , x], group=gNames, reorder=FALSE)

    } # end of sub-dt

    TL[, t] = calculateTL(preyed, tl, isResource, t)
    tl      = updateTL(TL=TL[, t], skeleton=skeleton, egg_tl=egg_tl)

    C[, t] = rowSums(Cx)

    CatchbyFleet[, , t] = rowSums(CatchbyFleetx, dims=2)
    YieldbyFleet[, , t] = 1e-6*rowSums(YieldbyFleetx, dims=2)
    YieldbyFleetByGroup[, , t] = 1e-6*rowSums(YieldbyFleetByGroupx, dims=2)

    bio = getBiomass(Bx, skeleton=skeleton)
    Bage[, t] = bio$Bage

    B[t, ] = bio$B  # biomass by functional group
    Y[t, ] = 1e-6*rowsum(rowSums(Yx), group=.getVar(pop, "name"), reorder = FALSE)

    # reproduction
    SSB = Nx[, xndt+1]*w_ssb     # abundance at the end (instantaneous reproduction)
    SSN = Nx[, xndt+1]*isMature
    R   = getRecruitment(SSB=SSB, SSN=SSN, t=t, skeleton=skeleton, groups=groups,
                         environment=environment)

    # save state variables
    No[, t+1] = updateN(Nx[, xndt+1], skeleton=skeleton, recruits=R,
                        plus=FALSE, isResource=isResourceGroup) # reproduction
    Lo[, t+1] = updateL(Lx[, xndt+1], skeleton=skeleton, egg_size=egg_size, isResource=isResourceGroup)

  } # end of t timestep

  if(isTRUE(verbose)) {
    pb = txtProgressBar(style=3)
    setTxtProgressBar(pb, (t-1)/ndt)
  }

  # post-processing outputs

  pop = updatePopulation(pop, N=No[, t+1], L=Lo[, t+1], w=a*Lo[, t+1]^b, TL=tl)

  colnames(Y) = colnames(B) = sapply(groups, FUN="[[", i="name")

  tcb = rowSums(B[, types=="functional_group"])

  if(!is.null(prices)) { # needs to be updated!
    # prices is matrix nGroups x nFleets
    prix = matrix(prices, nrow=nrow(B), ncol=sum(types=="functional_group"), byrow=TRUE)
    value = rowSums(Y[, types=="functional_group"]*prix)
  } else {
    value=NA
  }

  CPU.time = proc.time() - CPU.time

  if(isTRUE(verbose)) message(sprintf("\nModel run time: %0.1f seconds.", CPU.time[3]))

  # temp
  N = No
  L = Lo

  raw = list(N=N, L=L, C=C, TL=TL, Bage=Bage, CatchbyFleet=CatchbyFleet,
             YieldbyFleet=YieldbyFleet, YieldbyFleetByGroup=YieldbyFleetByGroup)

  out = list(B=B, Y=Y, tcb=tcb, value=value, raw=raw, pop=pop, CPU.time=CPU.time, time=time, groupTypes=types)

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




