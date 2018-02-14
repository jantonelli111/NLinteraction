##########################################################
## Function to calculate predicted values         ########
##########################################################

PredictionsMixture = function(XstarOld, XstarNew, designC, totalScans, nChains, zetaPost, 
                              betaList, betaCPost, k, ns) {
  n = dim(XstarNew)[1]
  PredictedPost = array(NA, dim=c(nChains, totalScans, n))
  hHatPost = array(NA, dim=c(nChains, totalScans, n))
  counter = 1
  for (ni in 1 : totalScans) {
    for (nc in 1 : nChains) {
      
      f_jhi_temp = array(NA, dim=c(n,k))
      for (h in 1 : k) {
        wv3 = which(zetaPost[nc,ni,,h] == 1)
        
        if (length(wv3) == 0) {
          f_jhi_temp[,h] = rep(0,n)
        } else if (length(wv3) == 1) {
          tempXstar = XstarNew[,wv3,]
          
          f_jhi_temp[,h] = tempXstar %*% as.vector(betaList[[ni]][[nc]][[h]])
        } else {
          tempXstarOld = matrix(NA, dim(XstarOld)[1], (ns+1)^length(wv3))
          nsMat = expand.grid(rep(list(1 : (ns + 1)), length(wv3)))
          for (ii in 1 : nrow(nsMat)) {
            tempXstarOld[,ii] = XstarOld[,wv3[1], nsMat[ii,1]]
            for (jj in 2 : length(wv3)) {
              tempXstarOld[,ii] = tempXstarOld[,ii]*XstarOld[,wv3[jj],nsMat[ii,jj]]
            }
          }
          tempXstarOld = cbind(dim(XstarOld)[1], tempXstarOld[,-1])
          
          tempXstarNew = matrix(NA, dim(XstarNew)[1], (ns+1)^length(wv3))
          nsMat = expand.grid(rep(list(1 : (ns + 1)), length(wv3)))
          for (ii in 1 : nrow(nsMat)) {
            tempXstarNew[,ii] = XstarNew[,wv3[1], nsMat[ii,1]]
            for (jj in 2 : length(wv3)) {
              tempXstarNew[,ii] = tempXstarNew[,ii]*XstarNew[,wv3[jj],nsMat[ii,jj]]
            }
          }
          tempXstarNew = cbind(rep(1, n), tempXstarNew[,-1])
          
          for (jj in 2 : ncol(tempXstarNew)) {
            tempXstarNew[,jj] = (tempXstarNew[,jj] - mean(tempXstarOld[,jj])) / sd(tempXstarOld[,jj])
          }
          
          f_jhi_temp[,h] = tempXstarNew %*% as.vector(betaList[[ni]][[nc]][[h]])
        }
      }
      PredictedPost[nc,counter,] = apply(f_jhi_temp, 1, sum) + 
        (designC %*% betaCPost[nc,ni,])
      hHatPost[nc,counter,] = apply(f_jhi_temp, 1, sum)
    }
    counter = counter + 1
  }
  return(list(PredictedPost = PredictedPost, hHatPost = hHatPost))
}

WaicMixture = function(Xstar, designC, totalScans, nChains, zetaPost, 
                                     betaList, betaCPost, sigmaPost, n, k, ns) {
  LPost = array(NA, dim=c(nChains, totalScans, n))
  counter = 1
  for (ni in 1:totalScans) {
    for (nc in 1 : nChains) {
      
      f_jhi_temp = array(NA, dim=c(n,k))
      for (h in 1 : k) {
        wv3 = which(zetaPost[nc,ni,,h] == 1)
        
        if (length(wv3) == 0) {
          f_jhi_temp[,h] = rep(0,n)
        } else if (length(wv3) == 1) {
          tempXstar = Xstar[,wv3,]
          
          f_jhi_temp[,h] = tempXstar %*% as.vector(betaList[[ni]][[nc]][[h]])
        } else {
          tempXstar = matrix(NA, n, (ns+1)^length(wv3))
          nsMat = expand.grid(rep(list(1 : (ns + 1)), length(wv3)))
          for (ii in 1 : nrow(nsMat)) {
            tempXstar[,ii] = Xstar[,wv3[1], nsMat[ii,1]]
            for (jj in 2 : length(wv3)) {
              tempXstar[,ii] = tempXstar[,ii]*Xstar[,wv3[jj],nsMat[ii,jj]]
            }
          }
          
          tempXstar = cbind(rep(1, n), scale(cbind(tempXstar[,-1])))
          
          f_jhi_temp[,h] = tempXstar %*% as.vector(betaList[[ni]][[nc]][[h]])
        }
      }
      
      mean = apply(f_jhi_temp, 1, sum) + 
        (designC %*% betaCPost[nc,ni,])
      
      pdf = (1 / sqrt(2*pi*sigmaPost[nc,ni]))*exp(-((Y - mean)^2)/(2*sigmaPost[nc,ni]))
      LPost[nc,counter,] = pdf
      
    }
    counter = counter + 1
  }
  
  lppd = apply(log(apply(LPost, c(1,3), mean)), 1, sum)
  pwaic = apply(apply(log(LPost), c(1,3), var), 1, sum)
  
  waic = lppd - pwaic
  return(-2*mean(waic))
}

##########################################################
## Function to plot interactions between 2 variables    ##
##########################################################

PlotInteraction = function(j1, j2, X, Xstar, C, quantile_j2, quantile_rest, 
                                  ns, gridLength, p, zetaPost, betaList, 
                                  betaCPost, totalScans, nChains, k,...) {

  n = dim(X)[1]
  NewDesignMat = matrix(NA, gridLength, p)
  for (j in 1 : p) {
    NewDesignMat[,j] = quantile(X[,j], quantile_rest)
  }
  
  NewDesignMat[,j1] = seq(quantile(X[,j1], .025), quantile(X[,j1], .975), length=gridLength)
  NewDesignMat[,j2] = quantile(X[,j2], quantile_j2)
  
  
  NewDesign = array(NA, dim=c(gridLength,p,ns+1))
  NewDesign[,,1] = 1
  for (j in 1 : p) {
    temp_ns_object = splines::ns(X[,j], df=ns)
    temp_sds = apply(temp_ns_object, 2, sd)
    temp_means = apply(temp_ns_object, 2, mean)
    NewDesign[,j,2:(ns+1)] = t((t(predict(temp_ns_object, NewDesignMat[,j])) - temp_means) / temp_sds)
  }
  
  NewDesignC = matrix(NA, gridLength, pc+1)
  NewDesignC[,1] = 1
  for (jc in 1 : pc) {
    NewDesignC[,jc+1] = mean(C[,jc])
  }


  predictions = PredictionsMixture(XstarOld = Xstar, XstarNew = NewDesign, designC = NewDesignC, 
                                   totalScans = totalScans, nChains = nChains,
                                   zetaPost = zetaPost, 
                                   betaList = betaList, betaCPost = betaCPost, k=k, ns=ns)
  

  dots = list(...)
  if ("ylim" %in% ls(dots)) {
    plot(NewDesignMat[,j1], apply(predictions$PredictedPost, 3, mean), type='l', lwd=3,...)
  } else {
    plot(NewDesignMat[,j1], apply(predictions$PredictedPost, 3, mean), type='l', lwd=3,
         ylim=range(c(apply(predictions$PredictedPost, 3, quantile, c(.025, .975)))),...)
  }
  
  polygon(x = c(NewDesignMat[,j1], rev(NewDesignMat[,j1])), 
          y = c(apply(predictions$PredictedPost, 3, quantile, .025), 
                rev(apply(predictions$PredictedPost, 3, quantile, .975))),
          density = -1, col = "grey", border=NA)
  lines(NewDesignMat[,j1], apply(predictions$PredictedPost, 3, mean), lwd=2)
}

##########################################################################
## Function to plot heatmap between 2 variables conditional on third    ##
##########################################################################

PlotInteractionHeatmapMeanSD = function(j1, j2, grid_j1, grid_j2, Xstar, X, C, quantile_j3, quantile_rest, 
                                  ns, p, zetaPost, betaList, minDist=Inf,
                                  betaCPost, totalScans, nChains, k,...) {
  
  colors = colorRampPalette(c("blue", "green", "red"))
  
  NewDesignMat = matrix(NA, length(grid_j1)*length(grid_j2), p)
  for (j in 1 : p) {
    NewDesignMat[,j] = quantile(X[,j], quantile_rest)
  }
  
  ij_grid = expand.grid(1:length(grid_j1), 1:length(grid_j2))
  
  NewDesignMat[,j1] = grid_j1[ij_grid[,1]]
  NewDesignMat[,j2] = grid_j2[ij_grid[,2]]
  
  
  NewDesign = array(NA, dim=c(length(grid_j1)*length(grid_j2),p,ns+1))
  NewDesign[,,1] = 1
  for (j in 1 : p) {
    temp_ns_object = splines::ns(X[,j], df=ns)
    temp_sds = apply(temp_ns_object, 2, sd)
    temp_means = apply(temp_ns_object, 2, mean)
    NewDesign[,j,2:(ns+1)] = t((t(predict(temp_ns_object, NewDesignMat[,j])) - temp_means) / temp_sds)
  }
  
  NewDesignC = matrix(NA, length(grid_j1)*length(grid_j2), pc+1)
  NewDesignC[,1] = 1
  for (jc in 1 : pc) {
    NewDesignC[,jc+1] = mean(C[,jc])
  }
  
  predictions = PredictionsMixture(XstarOld = Xstar, XstarNew = NewDesign, designC = NewDesignC, 
                                   totalScans = totalScans, nChains = nChains,
                                   zetaPost = zetaPost, 
                                   betaList = betaList, betaCPost = betaCPost, k=k, ns=ns)
  
  ## for plotting we need to turn predicted values into matrix
  NewPredictedPostMat = array(NA, dim=c(nChains, totalScans, length(grid_j1), length(grid_j2)))
  for (i in 1 : length(grid_j1)) {
    for (j in 1 : length(grid_j2)) {
      dist1 = grid_j1[i] - X[,j1]
      dist2 = grid_j2[j] - X[,j2]
      dist = sqrt(dist1^2 + dist2^2)
      
      if (min(dist) > minDist) {
        NewPredictedPostMat[,,i,j] = NA
      } else {
        wTemp = which(ij_grid[,1] == i & ij_grid[,2] == j)
        NewPredictedPostMat[,,i,j] = predictions$PredictedPost[,,wTemp] 
      }
    }
  }

  par(mfrow=c(1,2), pty='s')
  fields::image.plot(grid_j1, grid_j2, apply(NewPredictedPostMat, 3:4, mean, na.rm=TRUE),...)  
  fields::image.plot(grid_j1, grid_j2, apply(NewPredictedPostMat, 3:4, sd, na.rm=TRUE),...)  
}

PlotInteractionHeatmapMean = function(j1, j2, grid_j1, grid_j2, Xstar, X, C, quantile_j3, quantile_rest, 
                                        ns, p, zetaPost, betaList, minDist=Inf,
                                        betaCPost, totalScans, nChains, k,...) {
  
  colors = colorRampPalette(c("blue", "green", "red"))
  
  NewDesignMat = matrix(NA, length(grid_j1)*length(grid_j2), p)
  for (j in 1 : p) {
    NewDesignMat[,j] = quantile(X[,j], quantile_rest)
  }
  
  ij_grid = expand.grid(1:length(grid_j1), 1:length(grid_j2))
  
  NewDesignMat[,j1] = grid_j1[ij_grid[,1]]
  NewDesignMat[,j2] = grid_j2[ij_grid[,2]]
  
  
  NewDesign = array(NA, dim=c(length(grid_j1)*length(grid_j2),p,ns+1))
  NewDesign[,,1] = 1
  for (j in 1 : p) {
    temp_ns_object = splines::ns(X[,j], df=ns)
    temp_sds = apply(temp_ns_object, 2, sd)
    temp_means = apply(temp_ns_object, 2, mean)
    NewDesign[,j,2:(ns+1)] = t((t(predict(temp_ns_object, NewDesignMat[,j])) - temp_means) / temp_sds)
  }
  
  NewDesignC = matrix(NA, length(grid_j1)*length(grid_j2), pc+1)
  NewDesignC[,1] = 1
  for (jc in 1 : pc) {
    NewDesignC[,jc+1] = mean(C[,jc])
  }
  
  predictions = PredictionsMixture(XstarOld = Xstar, XstarNew = NewDesign, designC = NewDesignC, 
                                   totalScans = totalScans, nChains = nChains,
                                   zetaPost = zetaPost, 
                                   betaList = betaList, betaCPost = betaCPost, k=k, ns=ns)
  
  ## for plotting we need to turn predicted values into matrix
  NewPredictedPostMat = array(NA, dim=c(nChains, totalScans, length(grid_j1), length(grid_j2)))
  for (i in 1 : length(grid_j1)) {
    for (j in 1 : length(grid_j2)) {
      dist1 = grid_j1[i] - X[,j1]
      dist2 = grid_j2[j] - X[,j2]
      dist = sqrt(dist1^2 + dist2^2)
      
      if (min(dist) > minDist) {
        NewPredictedPostMat[,,i,j] = NA
      } else {
        wTemp = which(ij_grid[,1] == i & ij_grid[,2] == j)
        NewPredictedPostMat[,,i,j] = predictions$PredictedPost[,,wTemp] 
      }
    }
  }
  
  fields::image.plot(grid_j1, grid_j2, apply(NewPredictedPostMat, 3:4, mean, na.rm=TRUE),...)  
}

PlotInteractionHeatmapSD = function(j1, j2, grid_j1, grid_j2, Xstar, X, C, quantile_j3, quantile_rest, 
                                        ns, p, zetaPost, betaList, minDist=Inf,
                                        betaCPost, totalScans, nChains, k,...) {
  
  colors = colorRampPalette(c("blue", "green", "red"))
  
  NewDesignMat = matrix(NA, length(grid_j1)*length(grid_j2), p)
  for (j in 1 : p) {
    NewDesignMat[,j] = quantile(X[,j], quantile_rest)
  }
  
  ij_grid = expand.grid(1:length(grid_j1), 1:length(grid_j2))
  
  NewDesignMat[,j1] = grid_j1[ij_grid[,1]]
  NewDesignMat[,j2] = grid_j2[ij_grid[,2]]
  
  
  NewDesign = array(NA, dim=c(length(grid_j1)*length(grid_j2),p,ns+1))
  NewDesign[,,1] = 1
  for (j in 1 : p) {
    temp_ns_object = splines::ns(X[,j], df=ns)
    temp_sds = apply(temp_ns_object, 2, sd)
    temp_means = apply(temp_ns_object, 2, mean)
    NewDesign[,j,2:(ns+1)] = t((t(predict(temp_ns_object, NewDesignMat[,j])) - temp_means) / temp_sds)
  }
  
  NewDesignC = matrix(NA, length(grid_j1)*length(grid_j2), pc+1)
  NewDesignC[,1] = 1
  for (jc in 1 : pc) {
    NewDesignC[,jc+1] = mean(C[,jc])
  }
  
  predictions = PredictionsMixture(XstarOld = Xstar, XstarNew = NewDesign, designC = NewDesignC, 
                                   totalScans = totalScans, nChains = nChains,
                                   zetaPost = zetaPost, 
                                   betaList = betaList, betaCPost = betaCPost, k=k, ns=ns)
  
  ## for plotting we need to turn predicted values into matrix
  NewPredictedPostMat = array(NA, dim=c(nChains, totalScans, length(grid_j1), length(grid_j2)))
  for (i in 1 : length(grid_j1)) {
    for (j in 1 : length(grid_j2)) {
      dist1 = grid_j1[i] - X[,j1]
      dist2 = grid_j2[j] - X[,j2]
      dist = sqrt(dist1^2 + dist2^2)
      
      if (min(dist) > minDist) {
        NewPredictedPostMat[,,i,j] = NA
      } else {
        wTemp = which(ij_grid[,1] == i & ij_grid[,2] == j)
        NewPredictedPostMat[,,i,j] = predictions$PredictedPost[,,wTemp] 
      }
    }
  }
  
  fields::image.plot(grid_j1, grid_j2, apply(NewPredictedPostMat, 3:4, sd, na.rm=TRUE),...)  
}



##########################################################
## Function to estimate probability of interaction between 2 variables    ##
##########################################################

InteractionMatrix = function(zetaPost, totalScans, nChains, p, k) {
  intArray = array(NA, c(nChains, totalScans, p, p))
  counter = 1
  for (ni in 1 : totalScans) {
    for (nc in 1 : nChains) {
      for (i in 1 : (p-1)) {
        for (j in (i + 1) : p) {
          INDij = 0
          for (kk in 1 : k) {
            if (all(zetaPost[nc,ni,c(i,j),kk] == c(1,1))) INDij = 1
          }
          intArray[nc,counter,i,j] = INDij
        }
      }
    }
    counter = counter + 1
  }
  
  return(apply(intArray, c(3,4), mean))
}


InclusionVector = function(zetaPost, totalScans, nChains, p, k) {
  inclusionArray = array(NA, c(nChains, totalScans, p))
  for (ni in 1 : totalScans) {
    for (nc in 1 : nChains) {
      inclusionArray[nc,ni,] = (apply(zetaPost[nc,ni,,], 1, sum) > 0)
    }
  }
  return(apply(inclusionArray, 3, mean))
}
