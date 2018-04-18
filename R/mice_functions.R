
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

# getAccessibility2 = function(pop, predRange) {
#   size = .getVar(pop, "size")
#   .inside = function(x, y, limits) ((x >= y*limits[1]) & (x <= y*limits[2])) + 0
#   out = outer(X = size, Y = size, limits=predRange, FUN=.inside)
#   return(out)
# }
#
# getAccessibility3 = function(pop) {
#   size = .getVar(pop, "size")
#   psize = cbind(.getVar(pop, "psize_min"), .getVar(pop, "psize_max"))
#   .inside = function(limits, x) ((x >= limits[1]) & (x <= limits[2])) + 0
#   out = apply(psize, 1, FUN=.inside, x=size)
#   return(out)
# }

getAccessibility = function(pop, groups) {

  xi = .getVar(pop, "s_min") # minimum size of prey
  yi = .getVar(pop, "s_max") # maximum size of prey

  x = .getVar(pop, "psize_min") # minimum prey size (for predator)
  y = .getVar(pop, "psize_max") # maximum prey size (for predator)

  predator = cbind(x, y)
  prey     = cbind(xi, yi)

  .inside = function(predator, prey) {
    x = predator[1]
    y = predator[2]
    xi = prey[, 1]
    yi = prey[, 2]
    prey_range = yi - xi # size range of prey
    m = pmax.int(xi, x)
    M = pmin.int(yi, y)
    r = M - m
    out = numeric(nrow(prey))
    out[r==0] = 1
    ind = which(r>0)
    out[ind] = r[ind]/prey_range[ind]
    return(out)
  }

  grs = .getVar(pop, "name")
  dietMatrix = sapply(sapply(groups, FUN="[", i="target"), FUN="[", i=grs)
  colnames(dietMatrix) = sapply(groups, FUN="[[", i="name")
  dietMatrix = dietMatrix[, grs]

  out = apply(predator, 1, FUN=.inside, prey=prey)*dietMatrix

  return(out)
}

.getFeedingType = function(par) {
  return(par$feedType)
}

updateN = function(N, skeleton, plus=FALSE) {

  # skeleton[] = N
  N[N<1] = 0
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


updateL = function(L) {
  return(L)
}

updateC = function(C) {
  return(C)
}


getBiomass = function(N, w, skeleton) {

  B = relist(N*w, skeleton)
  B = 1e-6*unlist(lapply(B, FUN=sum))
  return(B)
}
