
#' @export
plot.mice.model = function(x, skip=0, ...) {

  sim = x
  xtime = sim$time[-1] - 0.5*diff(sim$time)[1]
  year = ceiling(sim$time[-1])
  nr = sum(sim$groupTypes == "resource")
  ns = ncol(sim$B) - nr
  ind = xtime > skip
  par(mfrow=c(ns+1,1), mar=c(0,0,0,0), oma=c(3,4,1,4))

  TE = calculateTrophicEff(x)
  ylim = c(0, 1.3*max(TE, na.rm=TRUE))
  matplot(x=unique(year)+0.5, TE, type="l", lty=1, axes=FALSE, ylim=ylim)
  mtext(toupper("Trophic Efficiency"), 3, adj=0.05, line=-1, cex=0.5)
  abline(h=0.10, lty=3)
  axis(4, las=1)
  box()

  for(i in seq_len(ns)) {
    icol =if(i!=7) i else 10
    b = 1e-3*sim$B[ind,i+nr]
    y = 1e-3*sim$Y[ind,i+nr]
    ymax = max(b, na.rm=TRUE)
    plot(xtime[ind], b, col=icol, axes=FALSE, type="l",
         ylim=c(0,ymax))
    lines(xtime[ind], 4*y, type="h", col=icol)
    axis(2*((i+1)%%2+1), las=1); box()
    mtext(toupper(colnames(sim$B)[i+nr]), 3, adj=0.95, line=-1, cex=0.5, col=icol)
  }
  axis(1)

  return(invisible())

}


# Auxiliar functions ------------------------------------------------------



calculateTrophicEff = function(object) {

  .calculateTE = function(TL, biomass) {
    x = as.numeric(TL)
    y = as.numeric(biomass)
    xd = cut(x, breaks=1:6, labels = FALSE, include.lowest = TRUE)
    yx = tapply(y, xd, sum)
    ye = exp(diff(log(yx)))
    te = setNames(rep(NA, 4), 2:5)
    te[names(ye)] = ye
    return(te)
  }

  year = ceiling(object$time[-1])
  TE = matrix(nrow=length(unique(year)), ncol=4)
  for(i in unique(year)) {
    iYear = which(year==i)
    TE[i, ] = .calculateTE(TL=object$raw$TL[, iYear],
                           biomass=object$raw$Bage[, iYear])
  }

  return(TE)

}


# Other plots -------------------------------------------------------------


#' @export
plotAll = function(x, skip=0, what="biomass", ...) {

  sim = x
  xtime = sim$time[-1] - 0.5*diff(sim$time)[1]
  year = ceiling(sim$time[-1])
  nr = sum(sim$groupTypes == "resource")
  ns = ncol(sim$B) - nr
  ind = xtime > skip

  X = if(what=="biomass") sim$B else sim$Y

  # par(mfrow=c(ns,1), mar=c(0,0,0,0), oma=c(3,4,1,4))

  for(i in seq_len(ns)) {
    icol =if(i!=7) i else 10
    b = 1e-3*X[ind,i+nr]
    ymax = max(b, na.rm=TRUE)
    plot(xtime[ind], b, col=icol, axes=FALSE, type="l",
         ylim=c(0,ymax))
    axis(2*((i+1)%%2+1), las=1); box()
    mtext(toupper(colnames(sim$B)[i+nr]), 3, adj=0.95, line=-1, cex=0.5, col=icol)
  }
  axis(1)

  return(invisible())

}


