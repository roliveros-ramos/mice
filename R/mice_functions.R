
.getVar = function(pop, what) {
  out = list()
  for(i in seq_along(pop)) {
    tmp = pop[[i]][[what]]
    out[[i]] = if(!is.null(tmp)) tmp else rep(NA, nrow(pop[[i]]))
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


getDietMatrix = function(pop, groups) {

  grs = .getVar(pop, "name")
  dietMatrix = sapply(sapply(groups, FUN="[", i="target"), FUN="[", i=grs)
  colnames(dietMatrix) = sapply(groups, FUN="[[", i="name")
  dietMatrix = dietMatrix[, grs]
  return(dietMatrix)

}

getAccessibility = function(size, pop, dietMatrix) {

  prey_size = size
  predator_range = matrix(nrow=nrow(size), ncol=2)

  out = data.frame(feedType=.getVar(pop, "feedType"))
  out$feedType[rowMeans(prey_size)<5] = "planktivorous"
  out$logPredLength = log10(prey_size[,1])
  predator_range[, 1] = 10^predict(preySizeModel$min, newdata=out) # prey upon with "size"
  out$logPredLength = log10(prey_size[,2])
  predator_range[, 2] = 10^predict(preySizeModel$max, newdata=out) # prey upon with "size"

  .inside = function(predator_range, prey_size) {
    if(anyNA(predator_range)) return(numeric(nrow(prey_size)))
    x = predator_range[1]
    y = predator_range[2]
    xi = prey_size[, 1]
    yi = prey_size[, 2]
    prey_range = yi - xi # size range of prey
    m = pmax.int(xi, x)
    M = pmin.int(yi, y)
    r = M - m
    out = numeric(nrow(prey_size))
    out[r==0] = 1
    ind = which(r>0)
    out[ind] = r[ind]/prey_range[ind]
    return(out)
  }

  out = apply(predator_range, 1, FUN=.inside, prey_size=prey_size)*dietMatrix

  return(out)
} # for long time step

# getAccessibility = function(size, pop, dietMatrix) {
#
#   out = data.frame(feedType=.getVar(pop, "feedType"))
#   out$feedType[size<5] = "planktivorous"
#   out$logPredLength = log10(size)
#
#   predator = matrix(nrow=nrow(out), ncol=2)
#
#   predator[, 1] = 10^predict(preySizeModel$min, newdata=out) # prey upon with "size"
#   predator[, 2] = 10^predict(preySizeModel$max, newdata=out) # prey upon with "size"
#
#   .inside = function(predator, prey) {
#     if(anyNA(predator)) return(numeric(prey))
#     out = 0 + (prey >= predator[1]) & (prey <= predator[2])
#     return(out)
#   }
#
#   out = apply(predator, 1, FUN=.inside, prey=size)*dietMatrix
#
#   return(out)
#
# }

.getFeedingType = function(par) {
  return(par$feedType)
}

# updateN = function(N, skeleton, plus=FALSE) {
#
#   N[N<1] = 0
#   N = relist(N, skeleton)
#
#   .updateN = function(x, plus) {
#     n = length(x)
#     nTmp = c(0, head(x, -1))
#     if(isTRUE(plus)) nTmp[n] = nTmp[n] + tail(x, 1)
#     return(nTmp)
#   }
#
#   out = c(unlist(lapply(N, FUN=.updateN, plus=plus)))
#
#   return(out)
# }

updateN = function(N, skeleton, recruits, plus) {

  N[N<1] = 0
  N = relist(N, skeleton)
  if(length(recruits)!=length(N)) stop("Recruits vector doesn't match functional groups number.")

  .updateN = function(x, R, plus) {
    n = length(x)
    nTmp = c(R, head(x, -1))
    if(isTRUE(plus)) nTmp[n] = nTmp[n] + tail(x, 1)
    return(nTmp)
  }

  for(i in seq_along(N)) N[[i]] = .updateN(x=N[[i]], R=recruits[i], plus=plus)

  out = c(unlist(N))

  return(out)
}

updateL = function(L, skeleton, egg_size) {

  L = relist(L, skeleton)
  if(length(egg_size)!=length(L)) stop("Egg size vector doesn't match functional groups number.")

  .updateL = function(x, egg_size) {
    return(c(egg_size, head(x, -1)))
  }

  for(i in seq_along(L)) L[[i]] = .updateL(x=L[[i]], egg_size=egg_size[i])

  out = c(unlist(L))

  return(out)

}

updateTL = function(TL, skeleton, egg_tl) {

  TL = relist(TL, skeleton)
  if(length(egg_tl)!=length(TL)) stop("Egg TL vector doesn't match functional groups number.")

  .updateTL = function(x, egg_tl) {
    return(c(egg_tl, head(x, -1)))
  }

  for(i in seq_along(TL)) TL[[i]] = .updateTL(x=TL[[i]], egg_tl=egg_tl[i])

  out = c(unlist(TL))

  return(out)

}

updateC = function(C) {
  return(C)
}


# getBiomass = function(N, w, skeleton) {
#
#   B = relist(N*w, skeleton)
#   B = 1e-6*unlist(lapply(B, FUN=sum))
#   return(B)
#
# }

getBiomass = function(B, skeleton) {

  n = ncol(B)
  w = rep(2, n)
  w[c(1,n)] = 1
  w = w/sum(w)
  Bage = 1e-6*colSums(t(B)*w)
  # if(isTRUE(bySize)) return(1e-6*B)
  B = relist(Bage, skeleton)
  B = unlist(lapply(B, FUN=sum))
  return(list(B=B, Bage=Bage))

}


