#' Run a simulation of the MICE model
#'
#' @param LHT List with the life history parameters
#' @param B0 Initial biomass for all species
#' @param predRange Minimum and maximum predator-prey ratios
#' @param Mstarv Maximum starvation mortality
#' @param Ystar Optimal annual food ration per gram of biomass
#' @param delta Fraction of prey population available to predators
#' @param dt Time step for simulation, in fractions of a year
#' @param T Time horizon for the simulation (years).
#'
#' @return A list with the abundance (N), length (L) and biomass (B) of
#' all the species modeled.
#' @export
#'
#' @examples
runMICE = function(groups, fleets, F=0, Mstarv = 0.3,
                   Ystar = 0.4, delta = 0.9, dt = 1/4, T = 100/dt, verbose=TRUE) {

  CPU.time = proc.time()

  pop = initGroups(groups=groups, dt=dt)
  access = getAccessibility(pop=pop)

  skeleton = .getSkeleton(pop)

  N0 = .getVar(pop, "N")
  L0 = .getVar(pop, "size")
  w_ssb = .getVar(pop, "ssb")
  isMature = .getVar(pop, "mature")
  w = .getVar(pop, "w")

  N = matrix(nrow=length(N0), ncol=T+2)
  L = matrix(nrow=length(L0), ncol=T+2)
  B = matrix(NA, nrow=T+1, ncol=length(skeleton))

  N[, 1] = N0
  L[, 1] = L0
  B[1, ] = getBiomass(N[, 1], w, skeleton)

  for(i in 1:T) {
    mort = calculateMortality(N=N[, i], w=w, F=F*dt, access=access, dt)
    M  = rowSums(mort$M)
    Ms = Mstarv*mort$starv
    Z = M + F*dt + Ms
    N[, i+1] = updateN(N[, i]*exp(-Z), skeleton=skeleton, plus=FALSE)
    SSB = w_ssb*N[, i+1]
    SSN = isMature*N[, i+1]
    R = getRecruitment(SSB=SSB, SSN=SSN, skeleton=skeleton, groups=groups)
    N[, i+1] = updateR(N[, i+1], skeleton=skeleton, R)
    B[i+1, ] = getBiomass(N[, i+1], w, skeleton=skeleton)
  }

  CPU.time = proc.time() - CPU.time

  if(isTRUE(verbose)) message(sprintf("Model run time: %0.1f seconds.", CPU.time[3]))

  return(list(N=N, L=L, B=B, CPU.time=CPU.time))

}
