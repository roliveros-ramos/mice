
getRecruitment = function(SSB, SSN, t, skeleton, groups, environment) {

  SSB = 1e-6*sapply(relist(SSB, skeleton), FUN=sum) # grams to tonnes
  SSN = sapply(relist(SSN, skeleton), FUN=sum) # total individuals, male + female

  R = numeric(length(groups))

  for(i in seq_along(groups)) {

    R[i] = groups[[i]]$recruitment(biomass=SSB[i], abundance=SSN[i], t=t,
                                   par=groups[[i]], environment=environment)

  } # end for groups

  return(R)

}

updateR = function(N, skeleton, R) {

  N = relist(N, skeleton)
  for(i in seq_along(N)) N[[i]][1] = R[i]
  out = c(unlist(N))

  return(out)
}

# Recruitment models ------------------------------------------------------


rho = function(T, par) {
  if(is.null(par[["Topt"]]))
    stop("Optimal temperature value must be provided.") else Topt = par[["Topt"]]
  if(is.null(par[["u"]])) u = 3 else u = par[["u"]]
  if(is.null(par[["Ts"]])) Ts = 10 else Ts = par[["Ts"]]

  out = pmax(1 - u^((T-Topt)-Ts/2) - u^(-(T-Topt)-Ts/2), 0)
  return(out)
}



#' @export
recruitment.ricker.spec = function(biomass, abundance, t, par, environment) {

  if(is.null(par[["alpha"]])) stop("ricker model: alpha has not been provided.")
  if(is.null(par[["beta"]])) stop("ricker model: beta has not been provided.")
  if(is.null(par[["delta"]])) par$delta = 1

  envVar = par$environment$recruitment
  xrho = if(is.null(envVar)) 1 else rho(environment[[envVar]][t], par)

  R = par$alpha*(biomass^par$delta)*exp(-par$beta*biomass)*xrho

  return(R)

}

#' @export
recruitment.resource.spec = function(biomass, abundance, t, par, environment) {

  R = 1e6*par$biomass[t] # tonnes to grams
  if(is.na(R)) stop("Recruitment forcing not available for time t.")
  return(R)

}

#' @export
recruitment.forcing.spec = function(biomass, abundance, t, par, environment) {

  R = par$recruits[t]
  if(is.na(R)) stop("Recruitment forcing not available for time t.")
  return(R)

}

#' @export
recruitment.pups.spec = function(biomass, abundance, t, par, environment) {

  if(is.null(par[["pups"]])) stop("pups model: pups number has not been provided.")

  envVar = par$environment$recruitment
  xrho = if(is.null(envVar)) 1 else rho(environment[[envVar]][t], par)

  R = pmax(0.5*abundance*par$pups*(1-0.5*biomass/par$K)*xrho, 0)
  return(R)

}


# all functions -----------------------------------------------------------

getRecruitment2 = function(SSB, SSN, t, skeleton, groups) {

  SSB = 1e-6*sapply(relist(SSB, skeleton), FUN=sum) # grams to tonnes
  SSN = sapply(relist(SSN, skeleton), FUN=sum) # total individuals, male + female

  R = numeric(length(groups))

  for(i in seq_along(groups)) {

    thisGroup = groups[[i]]

    if(thisGroup$type=="resource") {
      R[i] = 1e6*thisGroup$biomass[t] # tonnes to grams
    } else {
      recType = thisGroup$recType
      recBy   = thisGroup$recBy
      if(is.null(recType)) recType = "Ricker"
      if(is.null(recBy)) recBy = "biomass"
      recModel = match.fun(recType) # update names
      x = switch(recBy, biomass=SSB[i], number=SSN[i])
      R[i] = recModel(x, par=thisGroup)
    }

  } # end for groups

  return(R)

}

Ricker = function(x, par) par$alpha*x*exp(-par$beta*x)


