
# Initialization of Functional Groups -------------------------------------

checkGroups = function(groups, ndtPerYear, subDtPerYear, T, resources) {

  groupNames = sapply(groups, FUN = "[[", i="name")
  names(groups) = groupNames
  ndt = ndtPerYear*T

  if(!is.null(resources)) groups = updateResourceBiomass(groups, resources)

  for(i in seq_along(groups)) {

    if(is.null(groups[[i]]$type)) groups[[i]]$type = "functional_group"

    groups[[i]]$target = .initExclusion(preyExclude = groups[[i]]$preyExclude,
                                        target      = groups[[i]]$target,
                                        groupNames  = groupNames)

    # check for resources
    if(groups[[i]]$type=="resource") {

      groups[[i]]$biomass = .checkResourceBiomass(groups[[i]]$biomass,
                                                  ndtPerYear, subDtPerYear, T)
      groups[[i]]$recruitment = "resource"

      ngroups = groups[[i]]$ngroups
      if(is.null(ngroups)) groups[[i]]$ngroups = 5

      groups[[i]]$TL = .checkTL(groups[[i]]$TL)

    } # end of check for resources

    ## get recruitment function
    recruitmentType = groups[[i]][["recruitment", exact=TRUE]]

    if(is.null(recruitmentType)) recruitmentType = "ricker"

    # recruitment seasonality
    recruitmentSeasonality = groups[[i]]$recruitmentSeasonality
    if(is.null(recruitmentSeasonality))
      recruitmentSeasonality = rep(1, ndtPerYear)/ndtPerYear
    if(length(recruitmentSeasonality)==ndtPerYear) {
      recruitmentSeasonality = recruitmentSeasonality/sum(recruitmentSeasonality)
      recruitmentSeasonality = rep(recruitmentSeasonality, length=ndt)
    }
    if(length(recruitmentSeasonality)!=ndt)
      stop(sprintf("Recruitment seasonality must be a vector of length %s or %s.", ndtPerYear, ndt))
    if(anyNA(recruitmentSeasonality)) stop("Recruitment seasonality must not include NAs.")
    groups[[i]]$recruitmentSeasonality = recruitmentSeasonality

    recruitmentModel = match.fun(paste("recruitment", recruitmentType, "spec", sep="."))

    src = attr(recruitmentModel, "srcref")

    attr(recruitmentModel, "type") = recruitmentType
    attr(recruitmentModel, "source") = src
    class(recruitmentModel) = "recruitment.function"

    groups[[i]]$recruitment = recruitmentModel

    # end of recruitment

  }

  return(groups)

}

initGroups = function(groups, dt) {
  out = lapply(groups, FUN=.initGroups, dt=dt)
  names(out) = sapply(groups, FUN = "[[", i="name")
  return(out)
}

.initGroups = function(par, dt) {
  if(is.null(par$type)) return(.initSpecies(par=par, dt=dt))
  if(par$type=="resource") return(.initResources(par=par, dt=dt))
  return(.initSpecies(par=par, dt=dt))
}

.initSpecies = function(par, dt) {
  ns = 20 # number of bins for growth, constant for now...
  A = par$A
  Linf = par$Linf
  k = par$k
  t0 = par$t0
  a = par$a
  b = par$b
  am = par$am
  B0 = par$B0*1e6 # tonnes -> g
  psmin = par$predRange[1] # to be updated with other criteria
  psmax = par$predRange[2]
  M = if(is.null(par$M)) -log(0.01)/A else par$M
  egg_size = if(is.null(par$egg_size)) 1e-4 else par$egg_size

  out = data.frame(name=par$name, age=seq(0, A, by=dt), stringsAsFactors=FALSE)

  L1 = VB(age=out$age, Linf=Linf, k=k, t0=t0, egg_size=egg_size)    # t=t
  L2 = VB(age=out$age+dt, Linf=Linf, k=k, t0=t0, egg_size=egg_size) # t=t+dt

  out$size = L1 # size at t
  out$w    = a*out$size^b # weight at t

  out$dL = (L2-L1)
  out$a  = a
  out$b  = b

  out$s_mean = (L1+L2)/2 # mean length in [t, t+dt]
  out$w_mean = (a*(L2^(b+1) - L1^(b+1))/(L2-L1))/(b+1) # mean weight in [t, t+dt]

  out$s_min = L1 + 1*(L2-L1)/ns # lower threshold for prey size
  out$s_max = L1 + (ns-1)*(L2-L1)/ns # upper threshold for prey size

  mw = exp(-M*out$age)*out$w
  out$N = B0*(mw/sum(mw))/out$w
  out$mature = (out$age+dt >= am) + 0
  out$ssb = out$w_mean*out$mature

  out$feedType = .getFeedingType(par)
  # temporal
  out$feedType[out$s_mean<5] = "planktivorous"
  out$logPredLength = log10(out$s_min)
  out$psize_min = 10^predict(preySizeModel$min, newdata=out) # prey upon with mean length
  out$logPredLength = log10(out$s_max)
  out$psize_max = 10^predict(preySizeModel$max, newdata=out) # prey upon with mean length
  out$logPredLength = NULL

  out$Mold   = mortality.senecence.spec(out$age + 0.5*dt, par=par)
  out$Mb     = mortality.constant.spec(out$age, par=par)
  out$Mstarv = if(is.null(par$Mstarv)) 1 else par$Mstarv

  out$Ystar  = if(is.null(par$Ystar)) 3.5 else par$Ystar
  out$delta  = if(is.null(par$delta)) 0.9 else par$delta

  out$TL = NA # test is 2 is a better initialization
  out$egg_tl = if(is.null(par$egg_tl)) 2 else par$egg_tl

  # out$group = par$group
  out$type = par$type

  class(out) = c("mice_species", class(out))

  return(out)
}

.initResources = function(par, dt) {

  B0 = par$biomass[1]*1e6 # tonnes -> g

  ngroups = par$ngroups

  age = rep(0, ngroups)

  out = data.frame(name=par$name, age=age, stringsAsFactors=FALSE)

  logSize = log10(c(par$size_min, par$size_max))

  L = 10^seq(from=logSize[1], to=logSize[2], length=ngroups+1)

  L1 = head(L, -1)
  L2 = tail(L, -1)

  out$size = (L1+L2)/2 # size at t
  out$w    = 1 # weight at t

  out$dL = 0
  out$a  = 1
  out$b  = 0

  out$s_mean = (L1+L2)/2 # mean length in [t, t+dt]
  out$w_mean = 1 # mean weight in [t, t+dt]

  out$s_min = L1 # lower threshold for prey size
  out$s_max = L2 # upper threshold for prey size

  out$N = B0/ngroups # biomass is constant in log-space bin size (Sheldon, 1972).
  out$mature = 0
  out$ssb = 0
  out$psize_min = -1 # prey upon with mean length
  out$psize_max = -1 # prey upon with mean length
  # out$group = par$group

  out$Mold = 0
  out$Mb   = 1e-8
  out$Mstarv = 0

  out$Ystar  = 0
  out$delta  = if(is.null(par$delta)) 0.99999999 else par$delta

  TL = par$TL
  Dn = diff(TL)
  eta = 0.1
  i = seq_len(ngroups+1)-1
  Di = (log((eta^Dn -1)*i + ngroups) - log(ngroups))/log(eta)
  TLi = TL[1] + Di
  out$TL = TLi[-1] - diff(TLi)/2
  out$egg_tl = out$TL[1]

  out$type = par$type

  class(out) = c("mice_resources", class(out))

  return(out)
}





# Initialization of Fleets ------------------------------------------------

checkFleets = function(fleets) {

  fleetNames = sapply(fleets, FUN = "[[", i="name")
  names(fleets) = fleetNames

  return(fleets)

}

initFleets = function(fleets, groups, ndtPerYear, T) {
  groupNames = sapply(groups, FUN = "[[", i="name")
  out = lapply(fleets, FUN=.initFleet, ndtPerYear=ndtPerYear, T=T, groupNames=groupNames)
  names(out) = sapply(fleets, FUN = "[[", i="name")
  return(out)
}


.initFleet = function(par, ndtPerYear, T, groupNames, tiny=1e-3) {

  # add seasonality

  # move to checkFleets
  ndt = ndtPerYear*T
  F = par$F

  if(length(F)==1) F = rep(F, ndt) # constant F
  if(length(F)==T) F = rep(F, each=ndtPerYear) # F by Year

  if(length(F)!=ndt) stop("Fishing mortality is not right.") # improve: message

  par$F = F

  selectivityType = par$selectivity
  selectivityModel = match.fun(paste("selectivity", selectivityType, "spec", sep="."))

  par$selectivity = function(L) selectivityModel(x=L, par=par, tiny=tiny)
  attr(par$selectivity, "type") = selectivityType
  class(par$selectivity) = "selectivity.function"

  par$target = .initTargets(target=par$target, groupNames=groupNames)

  return(par)

}

.initTargets = function(target, groupNames) {

  out = setNames(numeric(length(groupNames)), nm = groupNames)

  if(is.null(target)) {
    out[] = 1
    return(out)
  }
  if(is.character(target)) {
    if(!all(target %in% groupNames)) {
      ind = which(!(target %in% groupNames))
      msg = paste(target[ind], collapse=",")
      msg = sprintf("target names (%s) do not match species group names.", msg)
      stop(msg)
    }
    out[target] = 1
    return(out)
  }
  if(!is.null(names(target))) {
    if(!all(names(target) %in% groupNames))
      stop("Target names do not match species group names.")
    out[names(target)] = target
    return(out)
  }
}


.initExclusion = function(preyExclude, target, groupNames) {

  out = setNames(numeric(length(groupNames)), nm = groupNames)
  out[] = 1

  if(is.null(preyExclude) & is.null(target)) {
    return(out)
  }

  if(!is.null(target)) {
    if(!is.numeric(target)) stop("Target must be numeric.")
    if(is.null(names(target))) stop("Target must be a named vector.")
    if(!all(names(target) %in% groupNames))
      stop("Target names do not match species group names.")
    out[names(target)] = target
  }

  if(!is.null(preyExclude)) {

    if(!is.character(preyExclude)) stop("preyExclude must be a character vector.")

    if(is.character(preyExclude)) {
      if(!all(preyExclude %in% groupNames)) {
        ind = which(!(preyExclude %in% groupNames))
        msg = paste(preyExclude[ind], collapse=",")
        msg = sprintf("preyExclude names (%s) do not match species group names.", msg)
        stop(msg)
      }
      out[preyExclude] = 0
    }
  }

  return(out)
}

# Initialization of the environmental variables

checkEnvironment = function(environment, ndt) {

  if(is.null(environment)) return(NULL)
  if(!is.list(environment)) stop("The 'environment' argument must be a list.")
  ll = sapply(environment, length)
  if(any(ll!=ndt)) stop("You must provide one value per time step.")
  return(environment)

}

# update Parameters -------------------------------------------------------

updateParameters = function(target, par) {

  if(!is.list(par)) stop("par must be a list.")

  parNames = names(par)
  if(any(parNames=="")) stop("all 'par' elements must be named.")

  validNames = sapply(target, FUN = "[[", i="name")

  for(parName in parNames) {
    gNames = names(par[[parName]])
    gNames = gNames[which(gNames %in% validNames)]
    if(length(gNames)==0) next
    for(iName in gNames) {
      target[[iName]][parName] = par[[parName]][iName]
    }
  }

  return(target)

}

updateResourceBiomass = function(target, par) {

  if(is.null(par)) return(target)

  if(!is.list(par)) stop("par must be a list.")

  gNames = names(par)
  if(any(gNames=="")) stop("all 'resources' elements must be named.")

  validNames = sapply(target, FUN = "[[", i="name")
  gNames = gNames[which(gNames %in% validNames)]

  if(length(gNames)==0) return(target)

  for(iName in gNames) {
    target[[iName]][["biomass"]] = par[[iName]]
  }

  return(target)

}


# Trophic level -----------------------------------------------------------

.checkTL = function(TL) {

  if(is.null(TL)) stop("TL must be provided for resources.")
  if(length(TL)>2) stop("TL must be a vector of length 2 (min and max TL).")
  if(length(TL)==1) TL = rep(TL, 2)
  if(TL[1]>TL[2]) stop("min TL is greater than max TL.")
  if(TL[1]==1 & TL[2]!=1) {
    TL[2]=1
    warning("Autotroph group, ignoring max TL.")
    }
  return(TL)

}

initialTL = function(preyed, tl) {

  maxit = 100
  preyed = t(t(preyed)/colSums(preyed))
  preyed[is.nan(preyed)] = 0

  TL = tl

  natl = sum(is.na(TL))
  mtl = 1
  niter = 1

  while(((natl>0) | (mtl>1e-5)) & (niter<maxit)) {
    TL0 = TL
    tlp = preyed*TL
    tlp[tlp==0] = NA
    omitTL = apply(tlp, 2, FUN=function(x) all(is.na(x)))
    TL = colSums(tlp, na.rm=TRUE) + 1
    TL[omitTL] = NA
    TL[!is.na(tl)] = tl[!is.na(tl)]
    natl = sum(is.na(TL))
    natl0 = sum(is.na(TL0))
    if(natl==0 & natl0==0) mtl = abs(mean(TL0-TL))
    niter = niter + 1
  }

  return(TL)

}


calculateTL = function(preyed, tl, isResource, t) {

  if(t==1) {
    tl  = initialTL(preyed, tl)
    return(tl)
  }

  # otl = TL[, t-1]

  preyed = t(t(preyed)/colSums(preyed))
  preyed[is.nan(preyed)] = 0
  tlp = tl*preyed
  TL = colSums(tlp, na.rm=TRUE) + 1
  TL[TL==1] = tl[TL==1] # if no prey, keep last TL.
  TL[isResource] = tl[isResource] # keep constant resources' TL

  return(TL)

}

