#' Estimate nonlinear association and identify interactions
#' 
#' This function takes in the observed data (y, x, c) and estimates
#' a potentially nonlinear, interactive effect between x and y while
#' adjusting for c linearly
#'
#' @param Y              The outcome to be analyzed.
#' @param X              an n by p matrix of exposures to evaluate. These should be continuous.
#' @param C              An n by q matrix of additional covariates to adjust for. If there are no
#'                       additional covariates (i.e q=0), then C should be set to NULL when specifying
#'                       the model
#' 
#' @param nScans         The number of MCMC scans to run.
#' @param nBurn          The number of MCMC scans that will be dropped as a burn-in.
#' @param thin           This number represents how many iterations between each scan.
#' 
#' @param c              The first hyperparameter of the inverse gamma prior on the residual variance
#' @param d              The second hyperparameter of the inverse gamma prior on the residual variance
#' 
#' @param sigB           Either a numeric value to be used for the value of sigma_beta,
#'                       the variance of the slab distribution,
#'                       or "EB" is specified to indicate that it will be estimated
#'                       via empirical Bayes. We recommend "EB" unless the user has
#'                       strong prior beliefs about the magnitude of the nonzero regression
#'                       coefficients
#'                       
#' @param k              The number of total components to allow in the model
#' @param ns             The degrees of freedom of the splines used to estimate nonlinear functions of the exposures
#' @param alph           The first hyperparameter for the beta prior on the probability of an exposure being included
#'                       in a given component
#' @param gamm           The second hyperparameter for the beta prior on the probability of an exposure being included
#'                       in a given component
#'                       
#' @param probSamp1      The proportion of the time the MCMC updates one variable at a time for each component as opposed to
#'                       two at a time. We have found that sampling one at a time is faster computationally, but we don't 
#'                       recommend probSamp1 = 1, because it is possible for some interactions to be significant without 
#'                       their main effects being significant, and these associations are much easier to pick up 2 at a time
#'                      
#' @param threshold      thresholding parameter when finding the lower bound for the slab variance. This parameter represents
#'                       the percentage of time a null association enters the model when tau_h = 0.5. Smaller values are more
#'                       conservative and prevent false discoveries. We recommend either 0.25 or 0.1. This lower bound is only
#'                       calculated if the "EB" option is selected, as it is used to make sure the empirical Bayes variance isn't
#'                       too small                                                                                                                                                                                             
#'
#' @return A list containing the full posterior draws of all parameters in the model, 
#'         the waic associated with the model, the posterior inclusion probabilities (PIPs)
#'         for each exposure entering into the model, and the matrix of 2-way interaction
#'         probabilities
#'         
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
                 probSamp1 = 0.9, threshold = 0.1) {
  
  designC = cbind(rep(1, dim(X)[1]), C)
  
  SigmaC = 1000*diag(dim(designC)[2])
  muC = rep(0, dim(designC)[2])
  
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
           MainPIP = inclusions,
           ns = ns,
           k = k)
  return(l)
}







#' Calculate posterior inclusion probabilities for an interaction
#' 
#' This function takes a fitted NLint object and can return the posterior
#' probability that any given set of variables are jointly interacting with each other.
#'
#' @param NLmod              The fitted NLint object
#' @param Xsub               The exposures for which you want to evaluate whether there is an interaction
#'
#' @return A number indicating the posterior probability that the exposures in Xsub are jointly interacting together.
#'         If Xsub is a scalar then the marginal posterior inclusion probability is returned, while if Xsub is a vector
#'         of length 2, the two way interaction is returned and so on.
#'         
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
#' ## Find 2 way interaction probability between exposure 4 and 5
#' InteractionProb(NLmod=NLmod, Xsub=c(4,5))


InteractionProb = function(NLmod, Xsub) {
  zeta = NLmod$posterior$zeta
  nChains = dim(zeta)[1]
  totalScans = dim(zeta)[2]
  p = dim(zeta)[3]
  k = dim(zeta)[4]
  
  intMat = matrix(NA, totalScans, nChains)
  for (ni in 1 : totalScans) {
    for (nc in 1 : nChains) {
      int_ind = 0
      for (h in 1 : k) {
        if (all(zeta[nc,ni,Xsub,h] == rep(1, length(Xsub)))) int_ind = 1
      }
      intMat[ni,nc] = int_ind
    }
  }
  return(mean(intMat))
}

