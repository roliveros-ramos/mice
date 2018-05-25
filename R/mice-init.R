
# Initialization of Functional Groups -------------------------------------

checkGroups = function(groups, ndt) {

  groupNames = sapply(groups, FUN = "[[", i="name")
  names(groups) = groupNames

  for(i in seq_along(groups)) {

    if(is.null(groups[[i]]$type)) groups[[i]]$type = "functional_group"

    groups[[i]]$target = .initExclusion(groups[[i]]$preyExclude, groupNames)

    # check for resources
    if(groups[[i]]$type=="resource") {
      if(is.null(groups[[i]]$biomass))
        stop("Biomass for each resource must be provided.")
      if(length(groups[[i]]$biomass)==1) {
        groups[[i]]$biomass = rep(groups[[i]]$biomass, ndt)
      }
      if(length(groups[[i]]$biomass)!=ndt)
        stop("One resource biomass value for every time step or a single one must be provided.")

      groups[[i]]$recruitment = "resource"

    } # end of check for resources

    ## get recruitment function
    recruitmentType = groups[[i]]$recruitment

    if(is.null(recruitmentType)) recruitmentType = "ricker"

    recruitmentModel = match.fun(paste("recruitment", recruitmentType, "spec", sep="."))

    src = attr(recruitmentModel, "srcref")

    # thisGroup = groups[[i]]
    #
    # recModel = function(biomass, abundance, t) {
    #   recruitmentModel(biomass=biomass, abundance=abundance, t=t, par=thisGroup)
    #  }

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
  out = data.frame(name=par$name, age=seq(0, A, by=dt), stringsAsFactors=FALSE)

  L1 = VB(age=out$age, Linf=Linf, k=k, t0=t0)    # t=t
  L2 = VB(age=out$age+dt, Linf=Linf, k=k, t0=t0) # t=t+dt

  out$size = L1 # size at t
  out$w    = a*out$size^b # weight at t

  out$s_mean = (L1+L2)/2 # mean length in [t, t+dt]
  out$w_mean = (a*(L2^(b+1) - L1^(b+1))/(L2-L1))/(b+1) # mean weight in [t, t+dt]

  out$s_min = L1 + 1*(L2-L1)/ns # lower threshold for prey size
  out$s_max = L1 + (ns-1)*(L2-L1)/ns # upper threshold for prey size

  mw = exp(-M*out$age)*out$w
  out$N = B0*(mw/sum(mw))/out$w
  out$mature = (out$age >= am) + 0
  out$ssb = out$w_mean*out$mature

  out$feedType = .getFeedingType(par)
  # temporal
  out$feedType[out$s_mean<5] = "planktivorous"
  out$logPredLength = log10(out$s_min)
  out$psize_min = 10^predict(preySizeModel$min, newdata=out) # prey upon with mean length
  out$logPredLength = log10(out$s_max)
  out$psize_max = 10^predict(preySizeModel$max, newdata=out) # prey upon with mean length
  out$logPredLength = NULL

  out$Mold = mortality.senecence.spec(out$age + 0.5*dt, par=par)
  # out$group = par$group
  out$type = par$type

  return(out)
}

.initResources = function(par, dt) {

  B0 = par$biomass[1]*1e6 # tonnes -> g

  out = data.frame(name=par$name, age=0, stringsAsFactors=FALSE)

  L1 = par$size_min    # t=t
  L2 = par$size_max # t=t+dt

  out$size = (L1+L2)/2 # size at t
  out$w    = 1 # weight at t

  out$s_mean = (L1+L2)/2 # mean length in [t, t+dt]
  out$w_mean = 1 # mean weight in [t, t+dt]

  out$s_min = L1 # lower threshold for prey size
  out$s_max = L2 # upper threshold for prey size

  out$N = B0
  out$mature = 0
  out$ssb = 0
  out$psize_min = -1 # prey upon with mean length
  out$psize_max = -1 # prey upon with mean length
  # out$group = par$group

  out$Mold = 0

  out$type = par$type
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


.initFleet = function(par, ndtPerYear, T, groupNames, tiny=1e-6) {

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


.initExclusion = function(preyExclude, groupNames) {

  out = setNames(numeric(length(groupNames)), nm = groupNames)
  out[] = 1

  if(is.null(preyExclude)) {
    return(out)
  }

  if(is.character(preyExclude)) {
    if(!all(preyExclude %in% groupNames)) {
      ind = which(!(preyExclude %in% groupNames))
      msg = paste(preyExclude[ind], collapse=",")
      msg = sprintf("preyExclude names (%s) do not match species group names.", msg)
      stop(msg)
    }
    out[preyExclude] = 0
    return(out)
  }

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

