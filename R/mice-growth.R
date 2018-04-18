
VB = function(age, Linf, k, t0, thr=0.5, egg=1e-4) {

  ref = Linf*(1-exp(-k*(thr-t0)))
  out1 = ((ref-egg)/thr)*age + egg
  out = Linf*(1-exp(-k*(age-t0)))
  out[age<thr] = out1[age<thr]
  return(out)

}
