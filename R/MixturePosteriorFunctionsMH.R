##########################################################
## Function to calculate predicted values         ########
##########################################################

PredictionsMixture_MH = function(XstarOld, XstarNew, designC, totalScans, nChains, zetaPost, 
                              betaList, betaCPost, k, ns) {
  n = dim(XstarOld)[1]
  n2 = dim(XstarNew)[1]
  PredictedPost = array(NA, dim=c(nChains, totalScans, n2))
  hHatPost = array(NA, dim=c(nChains, totalScans, n2))
  
  for (ni in 1 : totalScans) {
    for (nc in 1 : nChains) {
      
      model = sort(as.character(unique(unlist(ActivePred_MH(zetaPost[nc,ni,,], k=k)))))
      ## remove intercept from model
      model = model[-1]
      
      f_jhi_nc = rep(NA, n2)
      
      if (length(model) == 0) {
        f_jhi_nc = rep(0,n2)
      } else {
        ## remove any unnecessary terms (main effects if interactions are later included)
        modelNumeric = list()
        for (j in 1 : length(model)) {
          modelNumeric[[j]] = as.numeric(numextract_MH(model[j]))
        }
        
        modelNumeric2 = list()
        for (j1 in 1 : length(model)) {
          OK=TRUE
          for (j2 in (1 : length(model))[-j1]) {
            if(all(modelNumeric[[j1]] %in% modelNumeric[[j2]])) OK=FALSE
          }
          if (OK == TRUE) modelNumeric2[[length(modelNumeric2)+1]] = modelNumeric[[j1]]
        }
        
        ## First term in the model
        tempTerms = modelNumeric2[[1]]
        stl = length(tempTerms)
        
        if (stl == 1) {
          tempXstarOld2 = as.matrix(XstarOld[,tempTerms,-1])
          tempXstarNew2 = as.matrix(XstarNew[,tempTerms,-1])
        } else {
          ## Old
          tempXstarOld2 = matrix(NA, n, (ns+1)^(stl))
          
          nsMat = expand.grid(rep(list(1 : (ns + 1)),stl))
          for (ii in 1 : nrow(nsMat)) {
            tempXstarOld2[,ii] = XstarOld[,tempTerms[1], nsMat[ii,1]]
            for (jj in 2 : ncol(nsMat)) {
              tempXstarOld2[,ii] = tempXstarOld2[,ii]*XstarOld[,tempTerms[jj], nsMat[ii,jj]]
            }
          }
          tempXstarOld2 = tempXstarOld2[,-1]
          
          ## New
          tempXstarNew2 = matrix(NA, n2, (ns+1)^(stl))
          
          nsMat = expand.grid(rep(list(1 : (ns + 1)),stl))
          for (ii in 1 : nrow(nsMat)) {
            tempXstarNew2[,ii] = XstarNew[,tempTerms[1], nsMat[ii,1]]
            for (jj in 2 : ncol(nsMat)) {
              tempXstarNew2[,ii] = tempXstarNew2[,ii]*XstarNew[,tempTerms[jj], nsMat[ii,jj]]
            }
          }
          tempXstarNew2 = tempXstarNew2[,-1]
        }
        
        ## add any remaining terms
        if (length(modelNumeric2) > 1) {
          for (j in 2 : length(modelNumeric2)) {
            tempTerms = modelNumeric2[[j]]
            stl = length(tempTerms)
            
            if (stl == 1) {
              tempXstarNew3 = XstarNew[,tempTerms,]
              tempXstarOld3 = XstarOld[,tempTerms,]
            } else {
              ## Old
              tempXstarOld3 = matrix(NA, n, (ns+1)^(stl))
              
              nsMat = expand.grid(rep(list(1 : (ns + 1)),stl))
              for (ii in 1 : nrow(nsMat)) {
                tempXstarOld3[,ii] = XstarOld[,tempTerms[1], nsMat[ii,1]]
                for (jj in 2 : ncol(nsMat)) {
                  tempXstarOld3[,ii] = tempXstarOld3[,ii]*XstarOld[,tempTerms[jj], nsMat[ii,jj]]
                }
              }
              
              ## New
              tempXstarNew3 = matrix(NA, n2, (ns+1)^(stl))
              
              nsMat = expand.grid(rep(list(1 : (ns + 1)),stl))
              for (ii in 1 : nrow(nsMat)) {
                tempXstarNew3[,ii] = XstarNew[,tempTerms[1], nsMat[ii,1]]
                for (jj in 2 : ncol(nsMat)) {
                  tempXstarNew3[,ii] = tempXstarNew3[,ii]*XstarNew[,tempTerms[jj], nsMat[ii,jj]]
                }
              }
            }
            tempXstarNew2 = cbind(tempXstarNew2, tempXstarNew3[,-1])
            tempXstarOld2 = cbind(tempXstarOld2, tempXstarOld3[,-1])
          }
        }
        
        for (jj in 1 : ncol(tempXstarNew2)) {
          tempXstarNew2[,jj] = (tempXstarNew2[,jj] - mean(tempXstarOld2[,jj])) / 
            sd(tempXstarOld2[,jj])
        }
        
        f_jhi_nc = as.vector(tempXstarNew2 %*% as.vector(betaList[[ni]][[nc]]))
      }
      
      if (dim(designC)[2] == 1) {
        PredictedPost[nc,ni,] = f_jhi_nc + 
          (designC * betaCPost[nc,ni])
        hHatPost[nc,ni,] = f_jhi_nc
      } else {
        PredictedPost[nc,ni,] = f_jhi_nc + 
          (designC %*% betaCPost[nc,ni,])
        hHatPost[nc,ni,] = f_jhi_nc 
      }
    }
  }
  return(list(PredictedPost = PredictedPost, hHatPost = hHatPost))
}

WaicMixture_MH = function(Y, Xstar, designC, totalScans, nChains, zetaPost, 
                       betaList, betaCPost, sigmaPost, n, k, ns) {
  LPost = array(NA, dim=c(nChains, totalScans, n))
  counter = 1
  
  for (ni in 1:totalScans) {
    for (nc in 1 : nChains) {

      model = sort(as.character(unique(unlist(ActivePred_MH(zetaPost[nc,ni,,], k=k)))))
      ## remove intercept from model
      model = model[-1]
      
      f_jhi_nc = rep(NA, n)

      if (length(model) == 0) {
        f_jhi_nc = rep(0,n)
      } else {
        ## remove any unnecessary terms (main effects if interactions are later included)
        modelNumeric = list()
        for (j in 1 : length(model)) {
          modelNumeric[[j]] = as.numeric(numextract_MH(model[j]))
        }
        
        modelNumeric2 = list()
        for (j1 in 1 : length(model)) {
          OK=TRUE
          for (j2 in (1 : length(model))[-j1]) {
            if(all(modelNumeric[[j1]] %in% modelNumeric[[j2]])) OK=FALSE
          }
          if (OK == TRUE) modelNumeric2[[length(modelNumeric2)+1]] = modelNumeric[[j1]]
        }
        
        ## First term in the model
        tempTerms = modelNumeric2[[1]]
        stl = length(tempTerms)
        
        if (stl == 1) {
          tempXstar2 = Xstar[,tempTerms,-1]
        } else {
          tempXstar2 = matrix(NA, n, (ns+1)^(stl))
          
          nsMat = expand.grid(rep(list(1 : (ns + 1)),stl))
          for (ii in 1 : nrow(nsMat)) {
            tempXstar2[,ii] = Xstar[,tempTerms[1], nsMat[ii,1]]
            for (jj in 2 : ncol(nsMat)) {
              tempXstar2[,ii] = tempXstar2[,ii]*Xstar[,tempTerms[jj], nsMat[ii,jj]]
            }
          }
          tempXstar2 = tempXstar2[,-1]
        }
        
        ## add any remaining terms
        if (length(modelNumeric2) > 1) {
          for (j in 2 : length(modelNumeric2)) {
            tempTerms = modelNumeric2[[j]]
            stl = length(tempTerms)
            
            if (stl == 1) {
              tempXstar3 = Xstar[,tempTerms,]
            } else {
              tempXstar3 = matrix(NA, n, (ns+1)^(stl))
              
              nsMat = expand.grid(rep(list(1 : (ns + 1)),stl))
              for (ii in 1 : nrow(nsMat)) {
                tempXstar3[,ii] = Xstar[,tempTerms[1], nsMat[ii,1]]
                for (jj in 2 : ncol(nsMat)) {
                  tempXstar3[,ii] = tempXstar3[,ii]*Xstar[,tempTerms[jj], nsMat[ii,jj]]
                }
              }
            }
            tempXstar2 = cbind(tempXstar2, tempXstar3[,-1])
          }
        }
        
        tempXstar2 = scale(tempXstar2)
        
        f_jhi_nc = as.vector(tempXstar2 %*% as.vector(betaList[[ni]][[nc]]))
        
      }
      
      if (dim(designC)[2] == 1) {
        mean = f_jhi_nc + 
          (designC * betaCPost[nc,ni])
      } else {
        mean = f_jhi_nc + 
          (designC %*% betaCPost[nc,ni,]) 
      }
      
      pdf = (1 / sqrt(2*pi*sigmaPost[nc,ni]))*exp(-((Y - mean)^2)/(2*sigmaPost[nc,ni]))
      LPost[nc,ni,] = pdf
    }
  }
  
  lppd = apply(log(apply(LPost, c(1,3), mean)), 1, sum)
  pwaic = apply(apply(log(LPost), c(1,3), var), 1, sum)
  
  waic = lppd - pwaic
  return(-2*mean(waic))
}

##########################################################
## Function to plot interactions between 2 variables    ##
##########################################################

PlotInteraction_MH = function(j1, j2, X, Xstar, C, quantile_j2, quantile_rest, 
                           ns, gridLength, p, zetaPost, betaList, 
                           betaCPost, totalScans, nChains, k,...) {
  
  if (is.null(C)) {
    pc = 0
    NewDesignC = matrix(NA, gridLength, pc+1)
    NewDesignC[,1] = 1
  } else {
    pc = dim(C)[2]
    NewDesignC = matrix(NA, gridLength, pc+1)
    NewDesignC[,1] = 1
    for (jc in 1 : pc) {
      NewDesignC[,jc+1] = mean(C[,jc])
    }
  }
  
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
  
  
  predictions = PredictionsMixture_MH(XstarOld = Xstar, XstarNew = NewDesign, designC = NewDesignC, 
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

PlotInteractionHeatmapMeanSD_MH = function(j1, j2, grid_j1, grid_j2, Xstar, X, C, quantile_j3, quantile_rest, 
                                        ns, p, zetaPost, betaList, minDist=Inf,
                                        betaCPost, totalScans, nChains, k,...) {
  
  n = dim(X)[1]
  if (is.null(C)) {
    pc = 0
    NewDesignC = matrix(NA, length(grid_j1)*length(grid_j2), pc+1)
    NewDesignC[,1] = 1
  } else {
    pc = dim(C)[2]
    NewDesignC = matrix(NA, length(grid_j1)*length(grid_j2), pc+1)
    NewDesignC[,1] = 1
    for (jc in 1 : pc) {
      NewDesignC[,jc+1] = mean(C[,jc])
    }
  }
  
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
  
  predictions = PredictionsMixture_MH(XstarOld = Xstar, XstarNew = NewDesign, designC = NewDesignC, 
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

PlotInteractionHeatmapMean_MH = function(j1, j2, grid_j1, grid_j2, Xstar, X, C, quantile_j3, quantile_rest, 
                                      ns, p, zetaPost, betaList, minDist=Inf,
                                      betaCPost, totalScans, nChains, k,...) {
  
  n = dim(X)[1]
  if (is.null(C)) {
    pc = 0
    NewDesignC = matrix(NA, length(grid_j1)*length(grid_j2), pc+1)
    NewDesignC[,1] = 1
  } else {
    pc = dim(C)[2]
    NewDesignC = matrix(NA, length(grid_j1)*length(grid_j2), pc+1)
    NewDesignC[,1] = 1
    for (jc in 1 : pc) {
      NewDesignC[,jc+1] = mean(C[,jc])
    }
  }
  
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
  
  predictions = PredictionsMixture_MH(XstarOld = Xstar, XstarNew = NewDesign, designC = NewDesignC, 
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

PlotInteractionHeatmapSD_MH = function(j1, j2, grid_j1, grid_j2, Xstar, X, C, quantile_j3, quantile_rest, 
                                    ns, p, zetaPost, betaList, minDist=Inf,
                                    betaCPost, totalScans, nChains, k,...) {
  
  n = dim(X)[1]
  if (is.null(C)) {
    pc = 0
    NewDesignC = matrix(NA, length(grid_j1)*length(grid_j2), pc+1)
    NewDesignC[,1] = 1
  } else {
    pc = dim(C)[2]
    NewDesignC = matrix(NA, length(grid_j1)*length(grid_j2), pc+1)
    NewDesignC[,1] = 1
    for (jc in 1 : pc) {
      NewDesignC[,jc+1] = mean(C[,jc])
    }
  }
  
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
  
  predictions = PredictionsMixture_MH(XstarOld = Xstar, XstarNew = NewDesign, designC = NewDesignC, 
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

InteractionMatrix_MH = function(zetaPost, totalScans, nChains, p, k) {
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


InclusionVector_MH = function(zetaPost, totalScans, nChains, p, k) {
  inclusionArray = array(NA, c(nChains, totalScans, p))
  for (ni in 1 : totalScans) {
    for (nc in 1 : nChains) {
      inclusionArray[nc,ni,] = (apply(zetaPost[nc,ni,,], 1, sum) > 0)
    }
  }
  return(apply(inclusionArray, 3, mean))
}


