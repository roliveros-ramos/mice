
#' @export
plot.mice.model = function(x, skip=0, ...) {

  sim = x
  xtime = sim$time[-1] - 0.5*diff(sim$time)[1]
  nr = sum(sim$groupTypes == "resource")
  ns = ncol(sim$B) - nr
  ind = xtime > skip
  par(mfrow=c(ns,1), mar=c(0,0,0,0), oma=c(3,4,1,4))

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


