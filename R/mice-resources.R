

.checkResourceBiomass = function(biomass, ndtPerYear, subDtPerYear, T) {

  ndt  = ndtPerYear*T
  Ndt  = subDtPerYear*T
  xdt  = subDtPerYear/ndtPerYear

  if(is.null(biomass))
    stop("Biomass for each resource must be provided, check resources definition.")

  if(length(biomass)==1) return(rep(biomass, Ndt))

  if(length(biomass)==ndt) return(rep(biomass, each=xdt))

  stop("One resource biomass value for every time step or a single one must be provided.")

  return(out) # a vector of length T*subDtPerYear

}


getResources = function(groups) {

  isResource = which(sapply(groups, FUN="[", i="type") == "resource")
  groups = groups[isResource]
  ngroups = sapply(groups, FUN="[[", i="ngroups")
  biomass = sapply(groups, FUN="[[", i="biomass")
  biomass = apply(biomass, 1, rep.int, times=ngroups)
  biomass = biomass/rep.int(ngroups, times=ngroups)
  return(1e6*biomass)

}

updateResources = function(x, resources, isResource) {
  x[isResource] = resources
  return(x)
}



