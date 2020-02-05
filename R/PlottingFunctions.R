#' Plot 2-way interaction matrix
#' 
#' This function takes in an object which is a model from the NLint function
#' and plots the two-way interaction matrix
#'
#' @param NLmod              The fitted model
#' @param userLabels         True or false indicating whether the user will specify labels. Default is FALSE
#' @param xLabels            A p-dimensional vector of labels for the X axis. This will only get used if userLabels=TRUE
#' @param yLabels            A p-dimensional vector of labels for the Y axis. This will only get used if userLabels=TRUE
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
#' plotInt(NLmod = NLmod)

plotInt = function(NLmod, userLabels = FALSE, 
                   xLabels = NULL, yLabels = NULL) {
  p = dim(NLmod$posterior$zeta)[3]
  
  if (userLabels) {
    mat = matrix(NLmod$InteractionPIP, nrow=p, ncol=p,
                 dimnames = list(xLabels, yLabels))
    mat.melted = reshape2::melt(mat)
    
    g = ggplot2::ggplot(mat.melted, ggplot2::aes(x = Var1, y = Var2, fill = value)) + 
      ggplot2::geom_tile() + ggplot2::coord_equal() +
      ggplot2::ylab("") + ggplot2::xlab("") + 
      ggplot2::ggtitle("Posterior inclusion probabilities") +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face="bold")) + 
      ggplot2::theme(legend.key.size = ggplot2::unit(1.2, "cm")) +
      ggplot2::theme(legend.text = ggplot2::element_text(size = 15)) + 
      ggplot2::theme(legend.title=ggplot2::element_text(size=15, face="bold")) + 
      ggplot2::theme(axis.text = ggplot2::element_text(size=14))
  } else {
    mat = matrix(NLmod$InteractionPIP, nrow=p, ncol=p,
                 dimnames = list(paste("X", 1:p, sep=''), paste("X", 1:p, sep='')))
    mat.melted = reshape2::melt(mat)
    
    g = ggplot2::ggplot(mat.melted, ggplot2::aes(x = Var1, y = Var2, fill = value)) + 
      ggplot2::geom_tile() + ggplot2::coord_equal() +
      ggplot2::ylab("") + ggplot2::xlab("") + 
      ggplot2::ggtitle("Posterior inclusion probabilities") +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face="bold")) + 
      ggplot2::theme(legend.key.size = ggplot2::unit(1.2, "cm")) +
      ggplot2::theme(legend.text = ggplot2::element_text(size = 15)) + 
      ggplot2::theme(legend.title=ggplot2::element_text(size=15, face="bold")) + 
      ggplot2::theme(axis.text = ggplot2::element_text(size=14))
  }
  
  return(g)
}


#' Plot posterior predictive exposure response surface for an exposure
#' 
#' This function takes in an object which is a model from the NLint function
#' and plots the exposure response curve and pointwise 95% credible interval
#' band for exposure j1, while fixing exposure j2 at quantile_j2, the
#' remaining exposures at quantile_rest, and the covariates at their mean
#'
#' @param NLmod              The fitted model
#' @param X                  an n by p matrix of exposures to evaluate. These should be continuous.
#' @param C                  An n by q matrix of additional covariates to adjust for.
#' @param j1                 The exposure for which we will plot the exposure response function
#' @param j2                 The second exposure which we will fix at a certain value. This should
#'                           be a secondary exposure of interest in the sense that we want to see
#'                           if the exposure response curve of j1 is different at different levels of
#'                           j2 (i.e interaction)
#' @param gridLength         The number of points to evaluate the curve at. More points make for a smoother curve,
#'                           but also increase the computation time associated with calculating the curve
#' @param quantile_j2        The quantile at which we will fix exposure j2
#' @param quantile_rest      The quantile at which we will fix the remaining exposures
#' @param ...                Plotting parameters to be passed on such as ylim, ylab, main, etc.
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
#' par(mfrow=c(1,2), pty='s')
#' plotSurface1d(NLmod = NLmod, X=X, C=C, j1=4, j2=5,
#'               gridLength=30, quantile_j2=0.2, quantile_rest = 0.5)
#' plotSurface1d(NLmod = NLmod, X=X, C=C, j1=4, j2=5,
#'               gridLength=30, quantile_j2=0.8, quantile_rest = 0.5)

plotSurface1d = function(NLmod, X, C, j1, j2, gridLength = 100, 
                         quantile_j2, quantile_rest = 0.5,...) {
  n = dim(X)[1]
  ns = NLmod$ns
  k = NLmod$k
  p = dim(X)[2]
  
  ## Design matrix with splines
  Xstar = array(NA, dim=c(n,p,ns+1))
  Xstar[,,1] = 1
  for (j in 1 : p) {
    Xstar[,j,2:(ns+1)] = scale(splines::ns(X[,j], df=ns))
  }
  
  if (NLmod$speed == TRUE) {
    
    PlotInteraction(j1=j1, j2=j2, X=X, Xstar=Xstar, 
                    C=C, quantile_j2=quantile_j2, 
                    quantile_rest=quantile_rest, 
                    ns=ns, gridLength=gridLength, p=p, 
                    zetaPost=NLmod$posterior$zeta, 
                    betaList=NLmod$posterior$beta, 
                    betaCPost=NLmod$posterior$betaC, 
                    totalScans=dim(NLmod$posterior$betaC)[2], 
                    nChains=dim(NLmod$posterior$betaC)[1], k=k,...)
  } else {
    PlotInteraction_MH(j1=j1, j2=j2, X=X, Xstar=Xstar, 
                    C=C, quantile_j2=quantile_j2, 
                    quantile_rest=quantile_rest, 
                    ns=ns, gridLength=gridLength, p=p, 
                    zetaPost=NLmod$posterior$zeta, 
                    betaList=NLmod$posterior$beta, 
                    betaCPost=NLmod$posterior$betaC, 
                    totalScans=dim(NLmod$posterior$betaC)[2], 
                    nChains=dim(NLmod$posterior$betaC)[1], k=k,...)
  }
  
}






#' Plot bivariate posterior predictive exposure response surface
#' 
#' This function takes in an object which is a model from the NLint function
#' and plots the bivariate exposure response curve and pointwise posterior standard deviations
#' for exposures j1 and j2, keeping the remaining exposures at quantile_rest and the covariates at their mean
#'
#' @param NLmod              The fitted model
#' @param X                  an n by p matrix of exposures to evaluate. These should be continuous.
#' @param C                  An n by q matrix of additional covariates to adjust for.
#' @param j1                 The exposure for which we will plot the exposure response function
#' @param j2                 The second exposure which we will fix at a certain value. This should
#'                           be a secondary exposure of interest in the sense that we want to see
#'                           if the exposure response curve of j1 is different at different levels of
#'                           j2 (i.e interaction)
#' @param gridLength_j1      The number of points to evaluate the curve at for the first exposure. More points make for a smoother curve,
#'                           but also increase the computation time associated with calculating the curve
#' @param gridLength_j2      The number of points to evaluate the curve at for the first exposure. More points make for a smoother curve,
#'                           but also increase the computation time associated with calculating the curve
#' @param quantile_rest      The quantile at which we will fix the remaining exposures
#' @param minDist            If the euclidean distance between all observed data points is
#'                           greater than this number for any point on the plot, we wont' plot
#'                           the surface there. The default is infinity so that all points on
#'                           the grid are plotted.
#' @param ...                Plotting parameters to be passed on such as ylim, ylab, main, etc.
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
#' par(mfrow=c(1,2), pty='s')
#' plotSurface2d(NLmod = NLmod, X=X, C=C, j1=4, j2=5,
#'               gridLength_j1=20, gridLength_j2 = 20, quantile_rest = 0.5)


plotSurface2d = function(NLmod, X, C, j1, j2, gridLength_j1 = 20,
                         gridLength_j2 = 20, minDist=Inf,
                         quantile_rest = 0.5,...) {
  n = dim(X)[1]
  ns = NLmod$ns
  k = NLmod$k
  p = dim(X)[2]
  
  ## Design matrix with splines
  Xstar = array(NA, dim=c(n,p,ns+1))
  Xstar[,,1] = 1
  for (j in 1 : p) {
    Xstar[,j,2:(ns+1)] = scale(splines::ns(X[,j], df=ns))
  }
  
  grid_j1 = seq(quantile(X[,j1], .025), quantile(X[,j1], .975), length=gridLength_j1)
  grid_j2 = seq(quantile(X[,j2], .025), quantile(X[,j2], .975), length=gridLength_j2)
  
  if (NLmod$speed == TRUE) {
    
    PlotInteractionHeatmapMeanSD(j1=j1, j2=j2, X=X, Xstar=Xstar,
                                 C=C, grid_j1=grid_j1, grid_j2=grid_j2,
                                 quantile_rest=quantile_rest,
                                 ns=ns, p=p, minDist=minDist,
                                 zetaPost=NLmod$posterior$zeta, 
                                 betaList=NLmod$posterior$beta, 
                                 betaCPost=NLmod$posterior$betaC, 
                                 totalScans=dim(NLmod$posterior$betaC)[2], 
                                 nChains=dim(NLmod$posterior$betaC)[1], k=k,...)
  } else {
    PlotInteractionHeatmapMeanSD_MH(j1=j1, j2=j2, X=X, Xstar=Xstar,
                                    C=C, grid_j1=grid_j1, grid_j2=grid_j2,
                                    quantile_rest=quantile_rest,
                                    ns=ns, p=p, minDist=minDist,
                                    zetaPost=NLmod$posterior$zeta, 
                                    betaList=NLmod$posterior$beta, 
                                    betaCPost=NLmod$posterior$betaC, 
                                    totalScans=dim(NLmod$posterior$betaC)[2], 
                                    nChains=dim(NLmod$posterior$betaC)[1], k=k,...)
  }
  
}






#' Plot mean bivariate posterior predictive exposure response surface
#' 
#' This function takes in an object which is a model from the NLint function
#' and plots the bivariate exposure response curve
#' for exposures j1 and j2, keeping the remaining exposures at quantile_rest and the covariates at their mean
#'
#' @param NLmod              The fitted model
#' @param X                  an n by p matrix of exposures to evaluate. These should be continuous.
#' @param C                  An n by q matrix of additional covariates to adjust for.
#' @param j1                 The exposure for which we will plot the exposure response function
#' @param j2                 The second exposure which we will fix at a certain value. This should
#'                           be a secondary exposure of interest in the sense that we want to see
#'                           if the exposure response curve of j1 is different at different levels of
#'                           j2 (i.e interaction)
#' @param gridLength_j1      The number of points to evaluate the curve at for the first exposure. More points make for a smoother curve,
#'                           but also increase the computation time associated with calculating the curve
#' @param gridLength_j2      The number of points to evaluate the curve at for the first exposure. More points make for a smoother curve,
#'                           but also increase the computation time associated with calculating the curve
#' @param quantile_rest      The quantile at which we will fix the remaining exposures
#' @param minDist            If the euclidean distance between all observed data points is
#'                           greater than this number for any point on the plot, we wont' plot
#'                           the surface there. The default is infinity so that all points on
#'                           the grid are plotted.
#' @param ...                Plotting parameters to be passed on such as ylim, ylab, main, etc.
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
#' par(mfrow=c(1,2), pty='s')
#' plotSurface2d(NLmod = NLmod, X=X, C=C, j1=4, j2=5,
#'               gridLength_j1=20, gridLength_j2 = 20, quantile_rest = 0.5)


plotSurface2dMean = function(NLmod, X, C, j1, j2, gridLength_j1 = 20,
                         gridLength_j2 = 20, minDist=Inf,
                         quantile_rest = 0.5,...) {
  n = dim(X)[1]
  ns = NLmod$ns
  k = NLmod$k
  p = dim(X)[2]
  
  ## Design matrix with splines
  Xstar = array(NA, dim=c(n,p,ns+1))
  Xstar[,,1] = 1
  for (j in 1 : p) {
    Xstar[,j,2:(ns+1)] = scale(splines::ns(X[,j], df=ns))
  }
  
  grid_j1 = seq(quantile(X[,j1], .025), quantile(X[,j1], .975), length=gridLength_j1)
  grid_j2 = seq(quantile(X[,j2], .025), quantile(X[,j2], .975), length=gridLength_j2)
  
  if (NLmod$speed == TRUE) {
    
    PlotInteractionHeatmapMean(j1=j1, j2=j2, X=X, Xstar=Xstar,
                                 C=C, grid_j1=grid_j1, grid_j2=grid_j2,
                                 quantile_rest=quantile_rest,
                                 ns=ns, p=p, minDist=minDist,
                                 zetaPost=NLmod$posterior$zeta, 
                                 betaList=NLmod$posterior$beta, 
                                 betaCPost=NLmod$posterior$betaC, 
                                 totalScans=dim(NLmod$posterior$betaC)[2], 
                                 nChains=dim(NLmod$posterior$betaC)[1], k=k,...)
  } else {
    PlotInteractionHeatmapMean_MH(j1=j1, j2=j2, X=X, Xstar=Xstar,
                                    C=C, grid_j1=grid_j1, grid_j2=grid_j2,
                                    quantile_rest=quantile_rest,
                                    ns=ns, p=p, minDist=minDist,
                                    zetaPost=NLmod$posterior$zeta, 
                                    betaList=NLmod$posterior$beta, 
                                    betaCPost=NLmod$posterior$betaC, 
                                    totalScans=dim(NLmod$posterior$betaC)[2], 
                                    nChains=dim(NLmod$posterior$betaC)[1], k=k,...)
  }
  
}







#' Plot standard deviation of bivariate posterior predictive exposure response surface
#' 
#' This function takes in an object which is a model from the NLint function
#' and plots the pointwise posterior standard deviations of the bivariate exposure response surface
#' for exposures j1 and j2, keeping the remaining exposures at quantile_rest and the covariates at their mean
#'
#' @param NLmod              The fitted model
#' @param X                  an n by p matrix of exposures to evaluate. These should be continuous.
#' @param C                  An n by q matrix of additional covariates to adjust for.
#' @param j1                 The exposure for which we will plot the exposure response function
#' @param j2                 The second exposure which we will fix at a certain value. This should
#'                           be a secondary exposure of interest in the sense that we want to see
#'                           if the exposure response curve of j1 is different at different levels of
#'                           j2 (i.e interaction)
#' @param gridLength_j1      The number of points to evaluate the curve at for the first exposure. More points make for a smoother curve,
#'                           but also increase the computation time associated with calculating the curve
#' @param gridLength_j2      The number of points to evaluate the curve at for the first exposure. More points make for a smoother curve,
#'                           but also increase the computation time associated with calculating the curve
#' @param quantile_rest      The quantile at which we will fix the remaining exposures
#' @param minDist            If the euclidean distance between all observed data points is
#'                           greater than this number for any point on the plot, we wont' plot
#'                           the surface there. The default is infinity so that all points on
#'                           the grid are plotted.
#' @param ...                Plotting parameters to be passed on such as ylim, ylab, main, etc.
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
#' par(mfrow=c(1,2), pty='s')
#' plotSurface2d(NLmod = NLmod, X=X, C=C, j1=4, j2=5,
#'               gridLength_j1=20, gridLength_j2 = 20, quantile_rest = 0.5)


plotSurface2dSD = function(NLmod, X, C, j1, j2, gridLength_j1 = 20,
                         gridLength_j2 = 20, minDist=Inf,
                         quantile_rest = 0.5,...) {
  n = dim(X)[1]
  ns = NLmod$ns
  k = NLmod$k
  p = dim(X)[2]
  
  ## Design matrix with splines
  Xstar = array(NA, dim=c(n,p,ns+1))
  Xstar[,,1] = 1
  for (j in 1 : p) {
    Xstar[,j,2:(ns+1)] = scale(splines::ns(X[,j], df=ns))
  }
  
  grid_j1 = seq(quantile(X[,j1], .025), quantile(X[,j1], .975), length=gridLength_j1)
  grid_j2 = seq(quantile(X[,j2], .025), quantile(X[,j2], .975), length=gridLength_j2)
  
  if (NLmod$speed == TRUE) {
    
    PlotInteractionHeatmapSD(j1=j1, j2=j2, X=X, Xstar=Xstar,
                                 C=C, grid_j1=grid_j1, grid_j2=grid_j2,
                                 quantile_rest=quantile_rest,
                                 ns=ns, p=p, minDist=minDist,
                                 zetaPost=NLmod$posterior$zeta, 
                                 betaList=NLmod$posterior$beta, 
                                 betaCPost=NLmod$posterior$betaC, 
                                 totalScans=dim(NLmod$posterior$betaC)[2], 
                                 nChains=dim(NLmod$posterior$betaC)[1], k=k,...)
  } else {
    PlotInteractionHeatmapSD_MH(j1=j1, j2=j2, X=X, Xstar=Xstar,
                                    C=C, grid_j1=grid_j1, grid_j2=grid_j2,
                                    quantile_rest=quantile_rest,
                                    ns=ns, p=p, minDist=minDist,
                                    zetaPost=NLmod$posterior$zeta, 
                                    betaList=NLmod$posterior$beta, 
                                    betaCPost=NLmod$posterior$betaC, 
                                    totalScans=dim(NLmod$posterior$betaC)[2], 
                                    nChains=dim(NLmod$posterior$betaC)[1], k=k,...)
  }
  
}