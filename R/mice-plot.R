
#' @export
plot.mice.model = function(x, skip=0, ...) {

  sim = x
  xtime = sim$time[-1] - 0.5*diff(sim$time)[1]
  nr = sum(sim$groupTypes == "resource")
  ns = ncol(sim$B) - nr
  ind = xtime > skip
  par(mfrow=c(ns+1,1), mar=c(0,0,0,0), oma=c(3,4,1,4))

  calculateTrophicEff(x, names.arg = NA)

  for(i in seq_len(ns)) {
    icol =if(i!=7) i else 10
    b = 1e-3*sim$B[ind,i+nr]
    y = 1e-3*sim$Y[ind,i+nr]
    ymax = max(b, na.rm=TRUE)
    plot(xtime[ind], b, col=icol, axes=FALSE, type="l",
         ylim=c(0,ymax))
    lines(xtime[ind], y, type="h", col=icol)
    axis(2*((i+1)%%2+1), las=1); box()
    mtext(toupper(colnames(sim$B)[i+nr]), 3, adj=0.95, line=-1, cex=0.5, col=icol)
  }
  axis(1)

  return(invisible())

}


# Auxiliar functions ------------------------------------------------------


calculateTrophicEff = function(object, plot=TRUE, names.arg=NULL) {

  # B(TL) = B(1)*eff^(TL-1)
  # log B(TL) = log B(TL=1) + (TL-1)*log(eff)
  # y = log B(TL)
  # x = TL - 1
  # y = log B(1) + log(eff)*x

  x = as.numeric(object$raw$TL)
  y = object$raw$Bage

  xd = cut(x, breaks=1:6, labels = FALSE, include.lowest = TRUE)

  yx = tapply(y, xd, sum)
  xm = as.numeric(names(yx)) - 1
  ym = as.numeric(log10(yx))
  off = ym[1]
  ym = ym - off
  # how to better approximate the efficiency?
  mod = lm(ym ~ xm + 0)
  eff = 10^(coef(mod)[1])
  msg = sprintf("TE = %02.02f%%", 100*eff)
  if(plot) {
    # plot(xm,ym, xlab="TL", ylab="log(biomass)", las=1,
    #      pch=19, axes=FALSE)
    # abline(mod, col="red", lwd=2)
    # abline(c(0, log(0.1)), col="blue", lty=2)
    # mtext(msg, 3, adj=0.95, line=-2)
    # axis(1, at=axTicks(1), labels = axTicks(1)+1)
    # t2 = round(axTicks(2)+off)
    # axis(2, las=1, at=t2-off, labels = t2)
    # box()

    barplot(yx/yx[1], ylim=c(0,1.2), names.arg=names.arg, las=1)
    mtext(msg, 3, adj=0.95, line=-2)

  }

  names(eff) = NULL

  cat(msg, "\n")

}

