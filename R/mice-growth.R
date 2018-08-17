
VB = function(age, Linf, k, t0, thr=0.5, egg_size=1e-4) {

  ref = Linf*(1-exp(-k*(thr-t0)))
  out1 = ((ref-egg_size)/thr)*age + egg_size
  out = Linf*(1-exp(-k*(age-t0)))
  out[age<thr] = out1[age<thr]
  return(out)

}

getEggSize = function(pop) {
  out = sapply(pop, FUN=function(x) x$size[1])
  return(unlist(out))
}

getEggTL = function(pop) {
  out = sapply(pop, FUN=function(x) x$egg_tl[1])
  return(unlist(out))
}
