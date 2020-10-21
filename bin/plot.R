# Source: View(autoEstCont)
PlotContaminationFraction <- function(object, contaminationRange = c(0.01, 0.8), title, verbose = TRUE) {
  dd <- object$fit$dd
  priorRho <- object$fit$priorRho
  priorRhoStdDev <- object$fit$priorRhoStdDev
  p.L = function(x, alpha) {
    if (x == 0) {
      0
    }
    else {
      qgamma(alpha, x)
    }
  }
  p.U = function(x, alpha) {
    qgamma(1 - alpha, x + 1)
  }
  alpha = 0.95
  alpha = (1 - alpha)/2
  dd$rhoHigh = sapply(seq(nrow(dd)), function(e) p.U(dd$obsCnt[e], 
    alpha)/dd$expCnt[e])
  dd$rhoLow = sapply(seq(nrow(dd)), function(e) p.L(dd$obsCnt[e], 
    alpha)/dd$expCnt[e])
  rhoProbes = seq(0, 1, 0.001)
  v2 = (priorRhoStdDev/priorRho)^2
  k = 1 + v2^-2/2 * (1 + sqrt(1 + 4 * v2))
  theta = priorRho/(k - 1)
  tmp = sapply(rhoProbes, function(e) {
    tmp = dd[dd$useEst, ]
    mean(dgamma(e, k + tmp$obsCnt, scale = theta/(1 + theta * 
      tmp$expCnt)))
  })
  xx = dgamma(rhoProbes, k, scale = theta)
  w = which(rhoProbes >= contaminationRange[1] & rhoProbes <= 
    contaminationRange[2])
  rhoEst = (rhoProbes[w])[which.max(tmp[w])]
  rhoFWHM = range((rhoProbes[w])[which(tmp[w] >= (max(tmp[w])/2))])
  contEst = rhoEst
  if (verbose) 
    message(sprintf("Estimated global rho of %.2f", rhoEst))
  plot(rhoProbes, tmp, "l", xlim = c(0, 1), ylim = c(0, 
    max(c(xx, tmp))), frame.plot = FALSE, xlab = "Contamination Fraction", 
    ylab = "Probability Density")
  lines(rhoProbes, xx, lty = 2)
  abline(v = rhoProbes[which.max(tmp)], col = "red")
  legend(x = "topright", legend = c(sprintf("prior rho %g(+/-%g)", 
    priorRho, priorRhoStdDev), sprintf("post rho %g(%g,%g)", 
    rhoEst, rhoFWHM[1], rhoFWHM[2]), "rho max"), lty = c(2, 
    1, 1), col = c("black", "black", "red"), bty = "n")
  title(title)
}
