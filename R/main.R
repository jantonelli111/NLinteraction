#' Estimate nonlinear association and identify interactions
#' 
#' This function takes in the observed data (y, x, c) and estimates
#' a potentially nonlinear, interactive effect between x and y while
#' adjusting for c linearly
#'
#' @param y              The outcome to be analyzed
#' @param x              an n by p matrix of exposures to evaluate
#' @param c              An n by q matrix of additional covariates to adjust for
#'                                          
#'
#' @return An list of values that contain the treatment effect, confidence interval for the 
#'         treatment effect, and full posterior draws for the treatment effect.
#'
#' @export
#' @examples
#'
#' n = 200
#' p = 10
#' pc = 1
#' 
#' sigma = matrix(0.3, p, p)
#' diag(sigma) = 1
#' X = rmvnorm(n, mean=rep(0,p), sigma = sigma)
#' 
#' C = matrix(rnorm(n*pc), nrow=n)
#' 
#' TrueH = function(X) {
#'   return(0.5*(X[,2]*X[,3]) - 0.6*(X[,4]^2 * X[,5]))
#' }
#' 
#' Y = 5 + C + TrueH(X) + rnorm(n)
#' 
#' NLmod = NLint(Y=Y, X=X, C=C)
#' 
#' ## Print posterior inclusion probabilities
#' NLmod$MainPIP
#' 
#' ## Show the two way interaction probability matrix
#' NLmod$InteractionPIP
#' 



NLint = function(Y=Y, X=X, C=C, nChains = 2, nIter = 10000, 
                 nBurn = 2000, thin = 8, c = 0.001, d = 0.001,
                 sigB="EB", k = 15, ns = 3, alph=3, gamm=dim(X)[2], 
                 probSamp1 = 0.5, threshold = 0.25) {
  
  SigmaC = 1000*diag(dim(C)[2]+1)
  muC = rep(0, dim(C)[2]+1)
  
  designC = cbind(rep(1, dim(X)[1]), C)
  
  ## Prior mean vector
  muB = rep(0, ns)
  
  ## Design matrix with splines
  Xstar = array(NA, dim=c(n,p,ns+1))
  Xstar[,,1] = 1
  for (j in 1 : p) {
    Xstar[,j,2:(ns+1)] = scale(splines::ns(X[,j], df=ns))
  }
  

  ## Empirical Bayes if the user wants to use it
  if (sigB == "EB") {
    SigMin = MCMCmixtureMinSig(Y=Y, X=X, C=C, Xstar = Xstar, nPerms = 10, 
                               nIter = 500, c = c, d = d,
                               sigBstart=0.5, muB=muB,
                               SigmaC=SigmaC, muC=muC, 
                               k = k, ns = ns, threshold=0.25)
    
    SigEstEB = MCMCmixtureEB(Y=Y, X=X, C=C, Xstar = Xstar, nChains = nChains, nIter = nIter, 
                             nBurn = nBurn, thin = thin, c = c, d = d,
                             sigBstart=0.5, muB=muB,
                             SigmaC=SigmaC, muC=muC, 
                             k = k, ns = ns, alph=alph, gamm=gamm, probSamp1 = probSamp1)
    
    SigEst = max(SigEstEB, SigMin)
    
    posterior = MCMCmixture(Y=Y, X=X, C=C, Xstar = Xstar, nChains = nChains, nIter = nIter, 
                            nBurn = nBurn, thin = thin, c = c, d = d,
                            sigB=SigEst, muB=muB,
                            SigmaC=SigmaC, muC=muC, 
                            k = k, ns = ns, alph=alph, gamm=gamm, probSamp1 = probSamp1) 
  } else {
    posterior = MCMCmixture(Y=Y, X=X, C=C, Xstar = Xstar, nChains = nChains, nIter = nIter, 
                            nBurn = nBurn, thin = thin, c = c, d = d,
                            sigB=sigB, muB=muB,
                            SigmaC=SigmaC, muC=muC, 
                            k = k, ns = ns, alph=alph, gamm=gamm, probSamp1 = probSamp1) 
  }
  
  keep = nBurn + (1:floor((nIter - nBurn) / thin))*thin
  totalScans = length(keep)
  
  waic = WaicMixture(Xstar=Xstar, designC=designC, totalScans=totalScans, 
                     nChains=nChains, zetaPost = posterior$zeta, 
                     betaList = posterior$beta, betaCPost = posterior$betaC, 
                     sigmaPost = posterior$sigma, n=n, k=k, ns=ns)
  
  intMean = InteractionMatrix(zetaPost = posterior$zeta, totalScans = totalScans,
                                    nChains = 2, p = p, k = k)
  
  inclusions = InclusionVector(zetaPost = posterior$zeta, totalScans = totalScans,
                                    nChains = 2, p = p, k = k)
  
  
  l = list(posterior = posterior,
           waic = waic,
           InteractionPIP = intMean,
           MainPIP = inclusions)
  return(l)
}
