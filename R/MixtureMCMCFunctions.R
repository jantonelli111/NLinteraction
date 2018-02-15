ActivePred = function(zeta, k) {
  active = list()
  
  h = 1
  wv = which(zeta[,h] == 1)
  active[[h]] = c(0)
  if (length(wv) > 0) {
    active[[h]] = c(active[[h]], wv)
    
  }
  if (length(wv) > 1) {
    for (l in 2 : length(wv)) {
      comb = combn(wv, l)
      for (j in 1 : ncol(comb)) {
        
        active[[h]] = c(active[[h]], paste(sort(comb[1:l,j]), sep='', collapse="-"))
      }
    } 
  }
  
  for (h in 2 : k) {
    wv = which(zeta[,h] == 1)
    active[[h]] = c(0)
    if (length(wv) > 0) {
      for (l in 1 : length(wv)) {
        active[[h]] = c(active[[h]], wv[l]*(sum(zeta[wv[l],1:(h-1)]) == 0)) 
      }
    }
    if (length(wv) > 1) {
      for (l in 2 : length(wv)) {
        comb = combn(wv, l)
        for (j in 1 : ncol(comb)) {
          if (!paste(sort(comb[1:l,j]), sep='', collapse="-") %in% 
              unlist(active[1:(h-1)])) {
            active[[h]] = c(active[[h]], paste(sort(comb[1:l,j]), sep='', collapse="-"))
          }
        }
      }
    }
  }
  return(active)
}

UpdateBetaOne = function(tempZeta, f_jhi_nc, betaC,
                         sigmaP, tau, k, sigB, Xstar, tempBeta, h,
                         designC, ns) {
  
  activeZ = ActivePred(tempZeta, k=k)
  numZero = length(which(apply(tempZeta[,-h], 2, sum) == 0))
  
  groups = sample(1:p, 1, replace=FALSE)
  
  tempY = Y - apply(f_jhi_nc[,-h], 1, sum) -
    (designC %*% betaC)
  
  tempZeta10 = tempZeta00 = tempZeta01 = tempZeta
  tempZeta10[groups,h] = 1
  tempZeta00[groups,h] = 0
  
  wv = which(tempZeta00[,h] == 1)
  
  p00 = p10 = -Inf
  p01 = rep(-Inf, 2^length(wv) - 1)
  
  ## First calculate probability for reduced model
  
  if (length(wv) == 0) {
    p00 = log(1 - tau[h])
    # } else if (length(wv) == 1 & sum(tempZeta[wv,-h]) == 0) {
  } else if (length(wv) == 1) {
    
    if (h > 1 & sum(tempZeta[wv,1:(h-1)]) > 0) {
      p00 = log(tau[h] * (1 - tau[h]))
    } else {
      tempXstar = Xstar[,wv,]
      SigmaB = diag(sigmaP*sigB, dim(tempXstar)[2])
      muB = rep(0, dim(tempXstar)[2])
      muBeta = solve((t(tempXstar) %*% tempXstar)/
                       sigmaP + solve(SigmaB)) %*%
        ((t(tempXstar) %*% tempY)/
           sigmaP + solve(SigmaB) %*% muB)
      covBeta = solve((t(tempXstar) %*% tempXstar)/
                        sigmaP + solve(SigmaB))
      
      p00 = log(tau[h] * (1 - tau[h])) + 
        mvtnorm::dmvnorm(rep(0, ns), mean=muB[-1], sigma=as.matrix(SigmaB[-1,-1]), log=TRUE) - 
        mvtnorm::dmvnorm(rep(0, ns), mean=muBeta[-1], sigma=as.matrix(covBeta[-1,-1]), log=TRUE)
    }
    
  } else {
    
    tempXstar = matrix(NA, n, (ns+1)^length(wv))
    nsMat = expand.grid(rep(list(1 : (ns + 1)), length(wv)))
    activeMat = expand.grid(rep(list(c(0,rep(1,ns))), length(wv)))
    activeInd = "0"
    for (ii in 1 : nrow(nsMat)) {
      if (sum(activeMat[ii,]) > 0) {
        activeInd = c(activeInd, paste(sort(wv[which(activeMat[ii,]==1)]), collapse="-"))
      }
      tempXstar[,ii] = Xstar[,wv[1], nsMat[ii,1]]
      for (jj in 2 : length(wv)) {
        tempXstar[,ii] = tempXstar[,ii]*Xstar[,wv[jj],nsMat[ii,jj]]
      }
    }
    
    tempXstar = cbind(rep(1, n), scale(tempXstar[,-1]))
    
    ## exclude 1 because we dont' want the intercept in the active list
    w = 2:dim(tempXstar)[2]
    if (h > 1) {
      w = which(!activeInd %in% unlist(activeZ[1:(h-1)]))
    }
    
    if (length(w) == 0) {
      p00 = log((tau[h]^length(wv)) * (1 - tau[h]))
    } else if (length(w) == 1) {
      tempXstar2 = tempXstar[,w]
      
      SigmaB = sigB
      muB = rep(0, 1)
      muBeta = solve((t(tempXstar2) %*% tempXstar2)/
                       sigmaP + solve(SigmaB)) %*%
        ((t(tempXstar2) %*% tempY)/
           sigmaP + solve(SigmaB) %*% muB)
      covBeta = solve((t(tempXstar2) %*% tempXstar2)/
                        sigmaP + solve(SigmaB))
      p00 = log((tau[h]^length(wv)) * (1 - tau[h])) + 
        mvtnorm::dmvnorm(rep(0, length(muB)), mean=muB, sigma=as.matrix(SigmaB), log=TRUE) - 
        mvtnorm::dmvnorm(rep(0, length(muB)), mean=muBeta, sigma=as.matrix(covBeta), log=TRUE)
    } else {
      tempXstar2 = tempXstar[,w]
      
      SigmaB = diag(sigmaP*sigB, dim(tempXstar2)[2])
      muB = rep(0, dim(tempXstar2)[2])
      muBeta = solve((t(tempXstar2) %*% tempXstar2)/
                       sigmaP + solve(SigmaB)) %*%
        ((t(tempXstar2) %*% tempY)/
           sigmaP + solve(SigmaB) %*% muB)
      covBeta = solve((t(tempXstar2) %*% tempXstar2)/
                        sigmaP + solve(SigmaB))
      p00 = log((tau[h]^length(wv)) * (1 - tau[h])) + 
        mvtnorm::dmvnorm(rep(0, length(muB)), mean=muB, sigma=as.matrix(SigmaB), log=TRUE) - 
        mvtnorm::dmvnorm(rep(0, length(muB)), mean=muBeta, sigma=as.matrix(covBeta), log=TRUE)
      
    }
  }
  
  ## Now calculate probability for full interactions
  if ((numZero > 0 & sum(tempZeta10[,h]) > 1) |
      sum(tempZeta10[,h]) == 1) {
    if (length(wv) == 0) {
      
      if (h > 1 & sum(tempZeta[groups,1:(h-1)]) > 0) {
        p10 = log(tau[h])
      } else {
        tempXstar = Xstar[,groups,]
        SigmaB = diag(sigmaP*sigB, dim(tempXstar)[2])
        muB = rep(0, dim(tempXstar)[2])
        muBeta = solve((t(tempXstar) %*% tempXstar)/
                         sigmaP + solve(SigmaB)) %*%
          ((t(tempXstar) %*% tempY)/
             sigmaP + solve(SigmaB) %*% muB)
        covBeta = solve((t(tempXstar) %*% tempXstar)/
                          sigmaP + solve(SigmaB))
        
        p10 = log(tau[h]) + 
          mvtnorm::dmvnorm(rep(0, ns), mean=muB[-1], sigma=as.matrix(SigmaB[-1,-1]), log=TRUE) - 
          mvtnorm::dmvnorm(rep(0, ns), mean=muBeta[-1], sigma=as.matrix(covBeta[-1,-1]), log=TRUE)
      }
      
    } else {
      
      
      wv2 = c(wv, groups)
      
      tempXstar = matrix(NA, n, (ns+1)^length(wv2))
      nsMat = expand.grid(rep(list(1 : (ns + 1)), length(wv2)))
      activeMat = expand.grid(rep(list(c(0,rep(1,ns))), length(wv2)))
      activeInd = "0"
      for (ii in 1 : nrow(nsMat)) {
        if (sum(activeMat[ii,]) > 0) {
          activeInd = c(activeInd, paste(sort(wv2[which(activeMat[ii,]==1)]), collapse="-"))
        }
        tempXstar[,ii] = Xstar[,wv2[1], nsMat[ii,1]]
        for (jj in 2 : length(wv2)) {
          tempXstar[,ii] = tempXstar[,ii]*Xstar[,wv2[jj],nsMat[ii,jj]]
        }
      }
      
      tempXstar = cbind(rep(1, n), scale(tempXstar[,-1]))
      
      ## exclude 1 because we dont' want the intercept in the active list
      w = 2:dim(tempXstar)[2]
      if (h > 1) {
        w = which(!activeInd %in% unlist(activeZ[1:(h-1)]))
      }
      
      if (length(w) == 0) {
        p10 = log((tau[h]^length(wv2)))
      } else if (length(w) == 1) {
        tempXstar2 = tempXstar[,w]
        
        SigmaB = sigB
        muB = rep(0, 1)
        muBeta = solve((t(tempXstar2) %*% tempXstar2)/
                         sigmaP + solve(SigmaB)) %*%
          ((t(tempXstar2) %*% tempY)/
             sigmaP + solve(SigmaB) %*% muB)
        covBeta = solve((t(tempXstar2) %*% tempXstar2)/
                          sigmaP + solve(SigmaB))
        p10 = log((tau[h]^length(wv2))) + 
          mvtnorm::dmvnorm(rep(0, length(muB)), mean=muB, sigma=as.matrix(SigmaB), log=TRUE) - 
          mvtnorm::dmvnorm(rep(0, length(muB)), mean=muBeta, sigma=as.matrix(covBeta), log=TRUE)
      } else {
        tempXstar2 = tempXstar[,w]
        
        SigmaB = diag(sigmaP*sigB, dim(tempXstar2)[2])
        muB = rep(0, dim(tempXstar2)[2])
        muBeta = solve((t(tempXstar2) %*% tempXstar2)/
                         sigmaP + solve(SigmaB)) %*%
          ((t(tempXstar2) %*% tempY)/
             sigmaP + solve(SigmaB) %*% muB)
        covBeta = solve((t(tempXstar2) %*% tempXstar2)/
                          sigmaP + solve(SigmaB))
        p10 = log((tau[h]^length(wv2))) + 
          mvtnorm::dmvnorm(rep(0, length(muB)), mean=muB, sigma=as.matrix(SigmaB), log=TRUE) - 
          mvtnorm::dmvnorm(rep(0, length(muB)), mean=muBeta, sigma=as.matrix(covBeta), log=TRUE)
        
      }
    }
  }
  
  ## Now calculate probability for reduced model + effect elsewhere
  if (numZero > 0 & sum(tempZeta10[,h]) > 1) {
    
    PossMatTemp = expand.grid(rep(list(c(0,1)), length(wv)))
    PossMat = PossMatTemp[-nrow(PossMatTemp),]
    
    wh01 = which(apply(tempZeta, 2, sum) == 0 & 1:k != h)[1]
    
    if (length(wv) == 1) {
      
      tempZeta01 = tempZeta
      tempZeta01[,h] = 0
      tempZeta01[wv,h] = 1
      tempZeta01[groups,wh01] = 1
      activeZ01 = ActivePred(tempZeta01, k=k)
      
      tempXstar = cbind(Xstar[,wv,], Xstar[,groups,-1])
      activeInd = as.character(c(0, rep(wv, ns), rep(groups, ns)))
      
      if (h == 1) {
        w = 2:(ns+1)
        if (as.character(groups) %in% unlist(activeZ01[1:(wh01-1)])) {
          2:dim(tempXstar)[2]
        }
      } else if (wh01 == 1) {
        w = (ns+2):dim(tempXstar)[2]
        if (as.character(wv) %in% unlist(activeZ01[1:(h-1)])) {
          2:dim(tempXstar)[2]
        }
      } else if (h < wh01) {
        w = which(!activeInd[1:(ns+1)] %in% unlist(activeZ01[1:(h-1)]))
        w = c(w, which(!activeInd %in% unlist(activeZ01[1:(wh01-1)])))
      } else {
        w = which(!activeInd[1:(ns+1)] %in% unlist(activeZ01[1:(h-1)]))
        w = c(w, which(1:dim(tempXstar)[2] > (ns+1) &
                         !activeInd %in% unlist(activeZ01[1:(wh01-1)])))
      }
      
      if (length(w) == 0) {
        p01[1] = log(tau[h]^(length(wv)+1) * (1 - tau[h]))
      } else if (length(w) == 1) {
        tempXstar2 = tempXstar[,w]
        
        SigmaB = sigB
        muB = rep(0, 1)
        muBeta = solve((t(tempXstar2) %*% tempXstar2)/
                         sigmaP + solve(SigmaB)) %*%
          ((t(tempXstar2) %*% tempY)/
             sigmaP + solve(SigmaB) %*% muB)
        covBeta = solve((t(tempXstar2) %*% tempXstar2)/
                          sigmaP + solve(SigmaB))
        p01[1] = log(tau[h]^(length(wv)+1) * (1 - tau[h])) + 
          mvtnorm::dmvnorm(rep(0, length(muB)), mean=muB, sigma=as.matrix(SigmaB), log=TRUE) - 
          mvtnorm::dmvnorm(rep(0, length(muB)), mean=muBeta, sigma=as.matrix(covBeta), log=TRUE)
      } else {
        tempXstar2 = tempXstar[,w]
        
        SigmaB = diag(sigmaP*sigB, dim(tempXstar2)[2])
        muB = rep(0, dim(tempXstar2)[2])
        muBeta = solve((t(tempXstar2) %*% tempXstar2)/
                         sigmaP + solve(SigmaB)) %*%
          ((t(tempXstar2) %*% tempY)/
             sigmaP + solve(SigmaB) %*% muB)
        covBeta = solve((t(tempXstar2) %*% tempXstar2)/
                          sigmaP + solve(SigmaB))
        p01[1] = log(tau[h]^(length(wv)+1) * (1 - tau[h])) + 
          mvtnorm::dmvnorm(rep(0, length(muB)), mean=muB, sigma=as.matrix(SigmaB), log=TRUE) - 
          mvtnorm::dmvnorm(rep(0, length(muB)), mean=muBeta, sigma=as.matrix(covBeta), log=TRUE)
        
      }
    } else {
      
      for (iii in 1 : nrow(PossMat)) {
        tempXstar = matrix(NA, n, (ns+1)^length(wv))
        nsMat = expand.grid(rep(list(1 : (ns + 1)), length(wv)))
        activeMat = expand.grid(rep(list(c(0,rep(1,ns))), length(wv)))
        activeInd = "0"
        for (ii in 1 : nrow(nsMat)) {
          if (sum(activeMat[ii,]) > 0) {
            activeInd = c(activeInd, paste(sort(wv[which(activeMat[ii,]==1)]), collapse="-"))
          }
          tempXstar[,ii] = Xstar[,wv[1], nsMat[ii,1]]
          for (jj in 2 : length(wv)) {
            tempXstar[,ii] = tempXstar[,ii]*Xstar[,wv[jj],nsMat[ii,jj]]
          }
        }
        
        if (iii == 1) {
          tempXstar2 = Xstar[,groups,]
          activeInd = c(activeInd, as.character(rep(groups, ns)))
          wv3Temp = groups
        } else {
          wvTemp = which(PossMat[iii,] == 1)
          wv3Temp = c(wv[wvTemp], groups)
          
          tempXstar2 = matrix(NA, n, (ns+1)^length(wv3Temp))
          nsMat = expand.grid(rep(list(1 : (ns + 1)), length(wv3Temp)))
          activeMat = expand.grid(rep(list(c(0,rep(1,ns))), length(wv3Temp)))
          for (ii in 1 : nrow(nsMat)) {
            if (sum(activeMat[ii,]) > 0) {
              activeInd = c(activeInd, paste(sort(wv3Temp[which(activeMat[ii,]==1)]), collapse="-"))
            }
            tempXstar2[,ii] = Xstar[,wv3Temp[1], nsMat[ii,1]]
            for (jj in 2 : length(wv3Temp)) {
              tempXstar2[,ii] = tempXstar2[,ii]*Xstar[,wv3Temp[jj],nsMat[ii,jj]]
            }
          }
        }
        
        tempXstar = cbind(rep(1, n), scale(cbind(tempXstar[,-1], tempXstar2[,-1])))
        
        tempZeta01 = tempZeta
        tempZeta01[,h] = 0
        tempZeta01[wv,h] = 1
        tempZeta01[wv3Temp,wh01] = 1
        activeZ01 = ActivePred(tempZeta01, k=k)
        
        if (h == 1) {
          w = 2:((ns+1)^length(wv))
          w = c(w, which(!activeInd %in% unlist(activeZ01[1:(wh01-1)]) &
                           1:dim(tempXstar)[2] > ((ns+1)^length(wv))))
        } else if (wh01 == 1) {
          w = (((ns+1)^length(wv)) + 1):dim(tempXstar)[2]
          w = c(w, which(!activeInd %in% unlist(activeZ01[1:(h-1)]) &
                           1:dim(tempXstar)[2] <= ((ns+1)^length(wv))))
        } else if (h < wh01) {
          w = which(!activeInd[1:(((ns+1)^length(wv)))] %in% unlist(activeZ01[1:(h-1)]))
          w = c(w, which(!activeInd %in% unlist(activeZ01[1:(wh01-1)])))
        } else {
          w = which(!activeInd[1:(((ns+1)^length(wv)))] %in% unlist(activeZ01[1:(h-1)]))
          w = c(w, which(1:dim(tempXstar)[2] > (((ns+1)^length(wv))) &
                           !activeInd %in% unlist(activeZ01[1:(wh01-1)])))
        }
        
        if (length(w) == 0) {
          p01[iii] = log((tau[h]^(length(wv) + length(wv3Temp))) * 
                           (1 - tau[h]))
        } else if (length(w) == 1) {
          tempXstar2 = tempXstar[,w]
          
          SigmaB = sigB
          muB = rep(0, 1)
          muBeta = solve((t(tempXstar2) %*% tempXstar2)/
                           sigmaP + solve(SigmaB)) %*%
            ((t(tempXstar2) %*% tempY)/
               sigmaP + solve(SigmaB) %*% muB)
          covBeta = solve((t(tempXstar2) %*% tempXstar2)/
                            sigmaP + solve(SigmaB))
          p01[iii] = log((tau[h]^(length(wv) + length(wv3Temp))) * 
                           (1 - tau[h])) + 
            mvtnorm::dmvnorm(rep(0, length(muB)), mean=muB, sigma=as.matrix(SigmaB), log=TRUE) - 
            mvtnorm::dmvnorm(rep(0, length(muB)), mean=muBeta, sigma=as.matrix(covBeta), log=TRUE)
        } else {
          tempXstar2 = tempXstar[,w]
          
          SigmaB = diag(sigmaP*sigB, dim(tempXstar2)[2])
          muB = rep(0, dim(tempXstar2)[2])
          muBeta = solve((t(tempXstar2) %*% tempXstar2)/
                           sigmaP + solve(SigmaB)) %*%
            ((t(tempXstar2) %*% tempY)/
               sigmaP + solve(SigmaB) %*% muB)
          covBeta = solve((t(tempXstar2) %*% tempXstar2)/
                            sigmaP + solve(SigmaB))
          p01[iii] = log((tau[h]^(length(wv) + length(wv3Temp))) * 
                           (1 - tau[h])) + 
            mvtnorm::dmvnorm(rep(0, length(muB)), mean=muB, sigma=as.matrix(SigmaB), log=TRUE) - 
            mvtnorm::dmvnorm(rep(0, length(muB)), mean=muBeta, sigma=as.matrix(covBeta), log=TRUE)
          
        }
        
      }
      
    }
  }
  
  if (length(wv) == 0) {
    maxlog = max(c(p00, p10))
    updated_p = exp(c(p00, p10) - maxlog)
    samp = sample(1:2, 1, p = updated_p)
  } else if (length(wv) == 1 | numZero == 0) {
    maxlog = max(c(p00, p10, p01))
    updated_p = exp(c(p00, p10, p01) - maxlog)
    samp = sample(1:(length(updated_p)), 1, p = updated_p)
  } else {
    maxlog = max(c(p00, p10, p01))
    updated_p = exp(c(p00, p10, p01) - maxlog)
    samp = sample(1:(nrow(PossMat) + 2), 1, p = updated_p) 
  }
  
  if (samp == 1) {
    tempZeta = tempZeta00
  } else if (samp == 2) {
    tempZeta = tempZeta10
  } else if (samp == 3) {
    tempZeta = tempZeta00
    tempZeta[groups,wh01] = 1
  } else {
    tempZeta = tempZeta00
    tempZeta[c(wv, groups),wh01] = c(as.numeric(PossMat[samp-2,]), 1)
  }
  
  activeZ = ActivePred(tempZeta, k=k)
  
  wv3 = which(tempZeta[,h] == 1)
  if (length(wv3) == 0) {
    tempBeta[[h]] = c(0)
    f_jhi_nc[,h] = rep(0, n)
  } else if (length(wv3) == 1) {
    
    if (h > 1 & sum(tempZeta[wv3,1:(h-1)]) > 0) {
      tempBeta[[h]] = rep(0, ns+1)
      f_jhi_nc[,h] = rep(0,n)
    } else {
      tempXstar = Xstar[,wv3,]
      SigmaB = diag(sigmaP*sigB, dim(tempXstar)[2])
      muB = rep(0, dim(tempXstar)[2])
      muBeta = solve((t(tempXstar) %*% tempXstar)/
                       sigmaP + solve(SigmaB)) %*%
        ((t(tempXstar) %*% tempY)/
           sigmaP + solve(SigmaB) %*% muB)
      covBeta = solve((t(tempXstar) %*% tempXstar)/
                        sigmaP + solve(SigmaB))
      
      tempBeta[[h]] = mvtnorm::rmvnorm(1, muBeta, as.matrix(covBeta))
      f_jhi_nc[,h] = tempXstar %*% as.vector(tempBeta[[h]])
    }
  } else {
    
    tempXstar = matrix(NA, n, (ns+1)^length(wv3))
    nsMat = expand.grid(rep(list(1 : (ns + 1)), length(wv3)))
    activeMat = expand.grid(rep(list(c(0,rep(1,ns))), length(wv3)))
    activeInd = "0"
    for (ii in 1 : nrow(nsMat)) {
      if (sum(activeMat[ii,]) > 0) {
        activeInd = c(activeInd, paste(sort(wv3[which(activeMat[ii,]==1)]), collapse="-"))
      }
      tempXstar[,ii] = Xstar[,wv3[1], nsMat[ii,1]]
      for (jj in 2 : length(wv3)) {
        tempXstar[,ii] = tempXstar[,ii]*Xstar[,wv3[jj],nsMat[ii,jj]]
      }
    }
    
    tempXstar = cbind(rep(1, n), scale(cbind(tempXstar[,-1])))
    
    w = 1:dim(tempXstar)[2]
    if (h > 1) {
      w = which(!activeInd %in% unlist(activeZ[1:(h-1)]))
    }
    
    
    if (length(w) == 0) {
      tempBeta[[h]] = rep(0, dim(tempXstar)[2])
      f_jhi_nc[,h] = tempXstar %*% as.vector(tempBeta[[h]])
    } else if (length(w) == 1) {
      tempXstar2 = tempXstar[,w]
      
      SigmaB = sigB
      muB = rep(0, 1)
      muBeta = solve((t(tempXstar2) %*% tempXstar2)/
                       sigmaP + solve(SigmaB)) %*%
        ((t(tempXstar2) %*% tempY)/
           sigmaP + solve(SigmaB) %*% muB)
      covBeta = solve((t(tempXstar2) %*% tempXstar2)/
                        sigmaP + solve(SigmaB))
      
      temp = rep(0, dim(tempXstar)[2])
      temp[w] = mvtnorm::rmvnorm(1, muBeta, as.matrix(covBeta))
      tempBeta[[h]] = temp
      f_jhi_nc[,h] = tempXstar %*% as.vector(tempBeta[[h]])
    } else {
      tempXstar2 = tempXstar[,w]
      
      SigmaB = diag(sigmaP*sigB, dim(tempXstar2)[2])
      muB = rep(0, dim(tempXstar2)[2])
      muBeta = solve((t(tempXstar2) %*% tempXstar2)/
                       sigmaP + solve(SigmaB)) %*%
        ((t(tempXstar2) %*% tempY)/
           sigmaP + solve(SigmaB) %*% muB)
      covBeta = solve((t(tempXstar2) %*% tempXstar2)/
                        sigmaP + solve(SigmaB))
      
      temp = rep(0, dim(tempXstar)[2])
      temp[w] = mvtnorm::rmvnorm(1, muBeta, as.matrix(covBeta))
      tempBeta[[h]] = temp
      f_jhi_nc[,h] = tempXstar %*% as.vector(tempBeta[[h]])
    }
    
  }
  
  if (samp > 2) {
    
    wv3_01 = which(tempZeta[,wh01] == 1)
    if (length(wv3_01) == 0) {
      tempBeta[[wh01]] = c(0)
      f_jhi_nc[,wh01] = rep(0, n)
    } else if (length(wv3_01) == 1) {
      
      if (wh01 > 1 & sum(tempZeta[wv3_01,1:(wh01-1)]) > 0) {
        tempBeta[[wh01]] = rep(0, ns+1)
        f_jhi_nc[,wh01] = rep(0, n)
      } else {
        tempXstar = Xstar[,wv3_01,]
        SigmaB = diag(sigmaP*sigB, dim(tempXstar)[2])
        muB = rep(0, dim(tempXstar)[2])
        muBeta = solve((t(tempXstar) %*% tempXstar)/
                         sigmaP + solve(SigmaB)) %*%
          ((t(tempXstar) %*% tempY)/
             sigmaP + solve(SigmaB) %*% muB)
        covBeta = solve((t(tempXstar) %*% tempXstar)/
                          sigmaP + solve(SigmaB))
        
        tempBeta[[wh01]] = mvtnorm::rmvnorm(1, muBeta, as.matrix(covBeta))
        f_jhi_nc[,wh01] = tempXstar %*% as.vector(tempBeta[[wh01]])
      }
    } else {
      
      tempXstar = matrix(NA, n, (ns+1)^length(wv3_01))
      nsMat = expand.grid(rep(list(1 : (ns + 1)), length(wv3_01)))
      activeMat = expand.grid(rep(list(c(0,rep(1,ns))), length(wv3_01)))
      activeInd = "0"
      for (ii in 1 : nrow(nsMat)) {
        if (sum(activeMat[ii,]) > 0) {
          activeInd = c(activeInd, paste(sort(wv3_01[which(activeMat[ii,]==1)]), collapse="-"))
        }
        tempXstar[,ii] = Xstar[,wv3_01[1], nsMat[ii,1]]
        for (jj in 2 : length(wv3_01)) {
          tempXstar[,ii] = tempXstar[,ii]*Xstar[,wv3_01[jj],nsMat[ii,jj]]
        }
      }
      
      tempXstar = cbind(rep(1, n), scale(cbind(tempXstar[,-1])))
      
      w = 1:dim(tempXstar)[2]
      if (h > 1) {
        w = which(!activeInd %in% unlist(activeZ[1:(h-1)]))
      }
      
      
      if (length(w) == 0) {
        tempBeta[[wh01]] = rep(0, dim(tempXstar)[2])
        f_jhi_nc[,wh01] = tempXstar %*% as.vector(tempBeta[[wh01]])
      } else if (length(w) == 1) {
        tempXstar2 = tempXstar[,w]
        
        SigmaB = sigB
        muB = rep(0, 1)
        muBeta = solve((t(tempXstar2) %*% tempXstar2)/
                         sigmaP + solve(SigmaB)) %*%
          ((t(tempXstar2) %*% tempY)/
             sigmaP + solve(SigmaB) %*% muB)
        covBeta = solve((t(tempXstar2) %*% tempXstar2)/
                          sigmaP + solve(SigmaB))
        
        temp = rep(0, dim(tempXstar)[2])
        temp[w] = mvtnorm::rmvnorm(1, muBeta, as.matrix(covBeta))
        tempBeta[[wh01]] = temp
        f_jhi_nc[,wh01] = tempXstar %*% as.vector(tempBeta[[wh01]])
      } else {
        tempXstar2 = tempXstar[,w]
        
        SigmaB = diag(sigmaP*sigB, dim(tempXstar2)[2])
        muB = rep(0, dim(tempXstar2)[2])
        muBeta = solve((t(tempXstar2) %*% tempXstar2)/
                         sigmaP + solve(SigmaB)) %*%
          ((t(tempXstar2) %*% tempY)/
             sigmaP + solve(SigmaB) %*% muB)
        covBeta = solve((t(tempXstar2) %*% tempXstar2)/
                          sigmaP + solve(SigmaB))
        
        temp = rep(0, dim(tempXstar)[2])
        temp[w] = mvtnorm::rmvnorm(1, muBeta, as.matrix(covBeta))
        tempBeta[[wh01]] = temp
        f_jhi_nc[,wh01] = tempXstar %*% as.vector(tempBeta[[wh01]])
      }
      
    }
  }
  
  return(list(f_jhi_nc = f_jhi_nc,
              beta = tempBeta,
              zeta = tempZeta))
}


UpdateBetaTwo = function(tempZeta, f_jhi_nc, betaC,
                         sigmaP, tau, k, sigB, Xstar, tempBeta, h,
                         designC, ns) {
  
  activeZ = ActivePred(tempZeta, k=k)
  numZero = length(which(apply(tempZeta[,-h], 2, sum) == 0))
  
  groups = sample(1:p, 2, replace=FALSE)
  
  tempY = Y - apply(f_jhi_nc[,-h], 1, sum) -
    (designC %*% betaC)
  
  
  
  tempZeta10 = tempZeta00 = tempZeta01 = tempZeta11 = tempZeta
  tempZeta11[groups,h] = 1
  tempZeta00[groups,h] = 0
  tempZeta10[groups[1],h] = 1
  tempZeta10[groups[2],h] = 0
  tempZeta01[groups[1],h] = 0
  tempZeta01[groups[2],h] = 1
  
  wv = which(tempZeta00[,h] == 1)
  
  p00 = p10 = p01 = p11 = -Inf
  p00_1 = p00_2 = p00_12 = rep(-Inf, 2^length(wv) - 1)
  if (length(wv) == 0) {
    p00_1 = p00_2 = p00_12 = -Inf
  }
  p01_1 = p10_2 = rep(-Inf, 2^(length(wv)+1) - 1)
  
  ## First calculate probability for reduced model
  
  if (length(wv) == 0) {
    p00 = log((1 - tau[h])^2)
  } else if (length(wv) == 1) {
    
    if (h > 1 & sum(tempZeta[wv,1:(h-1)]) > 0) {
      p00 = log(tau[h] * (1 - tau[h]))
    } else {
      SigmaB = diag(sigmaP*sigB, ns+1)
      muB = rep(0, ns+1)
      muBeta = solve((t(Xstar[,wv,]) %*% Xstar[,wv,])/
                       sigmaP + solve(SigmaB)) %*%
        ((t(Xstar[,wv,]) %*% tempY)/
           sigmaP + solve(SigmaB) %*% muB)
      covBeta = solve((t(Xstar[,wv,]) %*% Xstar[,wv,])/
                        sigmaP + solve(SigmaB))
      
      p00 = log(tau[h] * (1 - tau[h])^2) + 
        mvtnorm::dmvnorm(rep(0, ns), mean=muB[-1], sigma=as.matrix(SigmaB[-1,-1]), log=TRUE) - 
        mvtnorm::dmvnorm(rep(0, ns), mean=muBeta[-1], sigma=as.matrix(covBeta[-1,-1]), log=TRUE)
    }
  } else {
    tempXstar = matrix(NA, n, (ns+1)^length(wv))
    nsMat = expand.grid(rep(list(1 : (ns + 1)), length(wv)))
    activeMat = expand.grid(rep(list(c(0,rep(1,ns))), length(wv)))
    activeInd = "0"
    for (ii in 1 : nrow(nsMat)) {
      if (sum(activeMat[ii,]) > 0) {
        activeInd = c(activeInd, paste(sort(wv[which(activeMat[ii,]==1)]), collapse="-"))
      }
      tempXstar[,ii] = Xstar[,wv[1], nsMat[ii,1]]
      for (jj in 2 : length(wv)) {
        tempXstar[,ii] = tempXstar[,ii]*Xstar[,wv[jj],nsMat[ii,jj]]
      }
    }
    
    ## exclude 1 because we dont' want the intercept in the active list
    w = 2:dim(tempXstar)[2]
    if (h > 1) {
      w = which(!activeInd %in% unlist(activeZ[1:(h-1)]))
    }
    
    if (length(w) == 0) {
      p00 = log((tau[h]^length(wv)) * (1 - tau[h]))
    } else if (length(w) == 1) {
      tempXstar2 = tempXstar[,w]
      
      SigmaB = sigB
      muB = rep(0, 1)
      muBeta = solve((t(tempXstar2) %*% tempXstar2)/
                       sigmaP + solve(SigmaB)) %*%
        ((t(tempXstar2) %*% tempY)/
           sigmaP + solve(SigmaB) %*% muB)
      covBeta = solve((t(tempXstar2) %*% tempXstar2)/
                        sigmaP + solve(SigmaB))
      p00 = log((tau[h]^length(wv)) * (1 - tau[h])) + 
        mvtnorm::dmvnorm(rep(0, length(muB)), mean=muB, sigma=as.matrix(SigmaB), log=TRUE) - 
        mvtnorm::dmvnorm(rep(0, length(muB)), mean=muBeta, sigma=as.matrix(covBeta), log=TRUE)
    } else {
      tempXstar2 = tempXstar[,w]
      
      SigmaB = diag(sigmaP*sigB, dim(tempXstar2)[2])
      muB = rep(0, dim(tempXstar2)[2])
      muBeta = solve((t(tempXstar2) %*% tempXstar2)/
                       sigmaP + solve(SigmaB)) %*%
        ((t(tempXstar2) %*% tempY)/
           sigmaP + solve(SigmaB) %*% muB)
      covBeta = solve((t(tempXstar2) %*% tempXstar2)/
                        sigmaP + solve(SigmaB))
      p00 = log((tau[h]^length(wv)) * (1 - tau[h])) + 
        mvtnorm::dmvnorm(rep(0, length(muB)), mean=muB, sigma=as.matrix(SigmaB), log=TRUE) - 
        mvtnorm::dmvnorm(rep(0, length(muB)), mean=muBeta, sigma=as.matrix(covBeta), log=TRUE)
      
    }
  }
  
  ## Now calculate probability for one higher interaction (group 1)
  if ((numZero > 0 & sum(tempZeta10[,h]) > 1) |
      sum(tempZeta10[,h]) == 1) {
    if (length(wv) == 0) {
      if (h > 1 & sum(tempZeta[groups[1],1:(h-1)]) > 0) {
        p10 = log(tau[h])
      } else {
        tempXstar = Xstar[,groups[1],]
        SigmaB = diag(sigmaP*sigB, dim(tempXstar)[2])
        muB = rep(0, dim(tempXstar)[2])
        muBeta = solve((t(tempXstar) %*% tempXstar)/
                         sigmaP + solve(SigmaB)) %*%
          ((t(tempXstar) %*% tempY)/
             sigmaP + solve(SigmaB) %*% muB)
        covBeta = solve((t(tempXstar) %*% tempXstar)/
                          sigmaP + solve(SigmaB))
        
        p10 = log(tau[h]) + 
          mvtnorm::dmvnorm(rep(0, ns), mean=muB[-1], sigma=as.matrix(SigmaB[-1,-1]), log=TRUE) - 
          mvtnorm::dmvnorm(rep(0, ns), mean=muBeta[-1], sigma=as.matrix(covBeta[-1,-1]), log=TRUE)
      }
    } else {
      wv2 = c(wv, groups[1])
      
      tempXstar = matrix(NA, n, (ns+1)^length(wv2))
      nsMat = expand.grid(rep(list(1 : (ns + 1)), length(wv2)))
      activeMat = expand.grid(rep(list(c(0,rep(1,ns))), length(wv2)))
      activeInd = "0"
      for (ii in 1 : nrow(nsMat)) {
        if (sum(activeMat[ii,]) > 0) {
          activeInd = c(activeInd, paste(sort(wv2[which(activeMat[ii,]==1)]), collapse="-"))
        }
        tempXstar[,ii] = Xstar[,wv2[1], nsMat[ii,1]]
        for (jj in 2 : length(wv2)) {
          tempXstar[,ii] = tempXstar[,ii]*Xstar[,wv2[jj],nsMat[ii,jj]]
        }
      }
      
      tempXstar = cbind(rep(1, n), scale(tempXstar[,-1]))
      
      ## exclude 1 because we dont' want the intercept in the active list
      w = 2:dim(tempXstar)[2]
      if (h > 1) {
        w = which(!activeInd %in% unlist(activeZ[1:(h-1)]))
      }
      
      if (length(w) == 0) {
        p10 = log((tau[h]^length(wv2)))
      } else if (length(w) == 1) {
        tempXstar2 = tempXstar[,w]
        
        SigmaB = sigB
        muB = rep(0, 1)
        muBeta = solve((t(tempXstar2) %*% tempXstar2)/
                         sigmaP + solve(SigmaB)) %*%
          ((t(tempXstar2) %*% tempY)/
             sigmaP + solve(SigmaB) %*% muB)
        covBeta = solve((t(tempXstar2) %*% tempXstar2)/
                          sigmaP + solve(SigmaB))
        p10 = log((tau[h]^length(wv2))) + 
          mvtnorm::dmvnorm(rep(0, length(muB)), mean=muB, sigma=as.matrix(SigmaB), log=TRUE) - 
          mvtnorm::dmvnorm(rep(0, length(muB)), mean=muBeta, sigma=as.matrix(covBeta), log=TRUE)
      } else {
        tempXstar2 = tempXstar[,w]
        
        SigmaB = diag(sigmaP*sigB, dim(tempXstar2)[2])
        muB = rep(0, dim(tempXstar2)[2])
        muBeta = solve((t(tempXstar2) %*% tempXstar2)/
                         sigmaP + solve(SigmaB)) %*%
          ((t(tempXstar2) %*% tempY)/
             sigmaP + solve(SigmaB) %*% muB)
        covBeta = solve((t(tempXstar2) %*% tempXstar2)/
                          sigmaP + solve(SigmaB))
        p10 = log((tau[h]^length(wv2))) + 
          mvtnorm::dmvnorm(rep(0, length(muB)), mean=muB, sigma=as.matrix(SigmaB), log=TRUE) - 
          mvtnorm::dmvnorm(rep(0, length(muB)), mean=muBeta, sigma=as.matrix(covBeta), log=TRUE)
        
      }
    }
  }
  
  ## Now calculate probability for one higher interaction (group 2)
  if ((numZero > 0 & sum(tempZeta01[,h]) > 1) |
      sum(tempZeta01[,h]) == 1) {
    if (length(wv) == 0) {
      if (h > 1 & sum(tempZeta[groups[2],1:(h-1)]) > 0) {
        p01 = log(tau[h])
      } else {
        tempXstar = Xstar[,groups[2],]
        SigmaB = diag(sigmaP*sigB, dim(tempXstar)[2])
        muB = rep(0, dim(tempXstar)[2])
        muBeta = solve((t(tempXstar) %*% tempXstar)/
                         sigmaP + solve(SigmaB)) %*%
          ((t(tempXstar) %*% tempY)/
             sigmaP + solve(SigmaB) %*% muB)
        covBeta = solve((t(tempXstar) %*% tempXstar)/
                          sigmaP + solve(SigmaB))
        
        p01 = log(tau[h]) + 
          mvtnorm::dmvnorm(rep(0, ns), mean=muB[-1], sigma=as.matrix(SigmaB[-1,-1]), log=TRUE) - 
          mvtnorm::dmvnorm(rep(0, ns), mean=muBeta[-1], sigma=as.matrix(covBeta[-1,-1]), log=TRUE)
      }
    } else {
      wv2 = c(wv, groups[2])
      
      tempXstar = matrix(NA, n, (ns+1)^length(wv2))
      nsMat = expand.grid(rep(list(1 : (ns + 1)), length(wv2)))
      activeMat = expand.grid(rep(list(c(0,rep(1,ns))), length(wv2)))
      activeInd = "0"
      for (ii in 1 : nrow(nsMat)) {
        if (sum(activeMat[ii,]) > 0) {
          activeInd = c(activeInd, paste(sort(wv2[which(activeMat[ii,]==1)]), collapse="-"))
        }
        tempXstar[,ii] = Xstar[,wv2[1], nsMat[ii,1]]
        for (jj in 2 : length(wv2)) {
          tempXstar[,ii] = tempXstar[,ii]*Xstar[,wv2[jj],nsMat[ii,jj]]
        }
      }
      
      tempXstar = cbind(rep(1, n), scale(tempXstar[,-1]))
      
      ## exclude 1 because we dont' want the intercept in the active list
      w = 2:dim(tempXstar)[2]
      if (h > 1) {
        w = which(!activeInd %in% unlist(activeZ[1:(h-1)]))
      }
      
      if (length(w) == 0) {
        p01 = log((tau[h]^length(wv2)))
      } else if (length(w) == 1) {
        tempXstar2 = tempXstar[,w]
        
        SigmaB = sigB
        muB = rep(0, 1)
        muBeta = solve((t(tempXstar2) %*% tempXstar2)/
                         sigmaP + solve(SigmaB)) %*%
          ((t(tempXstar2) %*% tempY)/
             sigmaP + solve(SigmaB) %*% muB)
        covBeta = solve((t(tempXstar2) %*% tempXstar2)/
                          sigmaP + solve(SigmaB))
        p01 = log((tau[h]^length(wv2))) + 
          mvtnorm::dmvnorm(rep(0, length(muB)), mean=muB, sigma=as.matrix(SigmaB), log=TRUE) - 
          mvtnorm::dmvnorm(rep(0, length(muB)), mean=muBeta, sigma=as.matrix(covBeta), log=TRUE)
      } else {
        tempXstar2 = tempXstar[,w]
        
        SigmaB = diag(sigmaP*sigB, dim(tempXstar2)[2])
        muB = rep(0, dim(tempXstar2)[2])
        muBeta = solve((t(tempXstar2) %*% tempXstar2)/
                         sigmaP + solve(SigmaB)) %*%
          ((t(tempXstar2) %*% tempY)/
             sigmaP + solve(SigmaB) %*% muB)
        covBeta = solve((t(tempXstar2) %*% tempXstar2)/
                          sigmaP + solve(SigmaB))
        p01 = log((tau[h]^length(wv2))) + 
          mvtnorm::dmvnorm(rep(0, length(muB)), mean=muB, sigma=as.matrix(SigmaB), log=TRUE) - 
          mvtnorm::dmvnorm(rep(0, length(muB)), mean=muBeta, sigma=as.matrix(covBeta), log=TRUE)
        
      }
    }
  }
  
  ## Now calculate probability for two higher order interactions (group 1 + 2)
  if (numZero > 0) {
    
    wv2 = c(wv, groups)
    tempXstar = matrix(NA, n, (ns+1)^length(wv2))
    nsMat = expand.grid(rep(list(1 : (ns + 1)), length(wv2)))
    activeMat = expand.grid(rep(list(c(0,rep(1,ns))), length(wv2)))
    activeInd = "0"
    for (ii in 1 : nrow(nsMat)) {
      if (sum(activeMat[ii,]) > 0) {
        activeInd = c(activeInd, paste(sort(wv2[which(activeMat[ii,]==1)]), collapse="-"))
      }
      tempXstar[,ii] = Xstar[,wv2[1], nsMat[ii,1]]
      for (jj in 2 : length(wv2)) {
        tempXstar[,ii] = tempXstar[,ii]*Xstar[,wv2[jj],nsMat[ii,jj]]
      }
    }
    
    tempXstar = cbind(rep(1, n), scale(tempXstar[,-1]))
    
    ## exclude 1 because we dont' want the intercept in the active list
    w = 2:dim(tempXstar)[2]
    if (h > 1) {
      w = which(!activeInd %in% unlist(activeZ[1:(h-1)]))
    }
    
    if (length(w) == 0) {
      p11 = log((tau[h]^length(wv2)))
    } else if (length(w) == 1) {
      tempXstar2 = tempXstar[,w]
      
      SigmaB = sigB
      muB = rep(0, 1)
      muBeta = solve((t(tempXstar2) %*% tempXstar2)/
                       sigmaP + solve(SigmaB)) %*%
        ((t(tempXstar2) %*% tempY)/
           sigmaP + solve(SigmaB) %*% muB)
      covBeta = solve((t(tempXstar2) %*% tempXstar2)/
                        sigmaP + solve(SigmaB))
      p11 = log((tau[h]^length(wv2))) + 
        mvtnorm::dmvnorm(rep(0, length(muB)), mean=muB, sigma=as.matrix(SigmaB), log=TRUE) - 
        mvtnorm::dmvnorm(rep(0, length(muB)), mean=muBeta, sigma=as.matrix(covBeta), log=TRUE)
    } else {
      tempXstar2 = tempXstar[,w]
      
      SigmaB = diag(sigmaP*sigB, dim(tempXstar2)[2])
      muB = rep(0, dim(tempXstar2)[2])
      muBeta = solve((t(tempXstar2) %*% tempXstar2)/
                       sigmaP + solve(SigmaB)) %*%
        ((t(tempXstar2) %*% tempY)/
           sigmaP + solve(SigmaB) %*% muB)
      covBeta = solve((t(tempXstar2) %*% tempXstar2)/
                        sigmaP + solve(SigmaB))
      p11 = log((tau[h]^length(wv2))) + 
        mvtnorm::dmvnorm(rep(0, length(muB)), mean=muB, sigma=as.matrix(SigmaB), log=TRUE) - 
        mvtnorm::dmvnorm(rep(0, length(muB)), mean=muBeta, sigma=as.matrix(covBeta), log=TRUE)
      
    }
  }
  
  
  ## Now calculate probability for reduced model + effect with group 1
  if ((numZero > 0 & sum(tempZeta10[,h]) > 1)) {
    
    PossMatTemp = expand.grid(rep(list(c(0,1)), length(wv)))
    PossMat00_1 = PossMatTemp[-nrow(PossMatTemp),]
    
    wh01 = which(apply(tempZeta, 2, sum) == 0 & 1:k != h)[1]
    
    if (length(wv) == 1) {
      
      tempZeta0010 = tempZeta
      tempZeta0010[,h] = 0
      tempZeta0010[wv,h] = 1
      tempZeta0010[groups[1],wh01] = 1
      activeZ0010 = ActivePred(tempZeta0010, k=k)
      
      tempXstar = cbind(Xstar[,wv,], Xstar[,groups[1],-1])
      activeInd = as.character(c(0, rep(wv, ns), rep(groups[1], ns)))
      
      if (h == 1) {
        w = 2:(ns+1)
        if (as.character(groups[1]) %in% unlist(activeZ0010[1:(wh01-1)])) {
          2:dim(tempXstar)[2]
        }
      } else if (wh01 == 1) {
        w = (ns+2):dim(tempXstar)[2]
        if (as.character(wv) %in% unlist(activeZ0010[1:(h-1)])) {
          2:dim(tempXstar)[2]
        }
      } else if (h < wh01) {
        w = which(!activeInd[1:(ns+1)] %in% unlist(activeZ0010[1:(h-1)]))
        w = c(w, which(!activeInd %in% unlist(activeZ0010[1:(wh01-1)])))
      } else {
        w = which(!activeInd[1:(ns+1)] %in% unlist(activeZ0010[1:(h-1)]))
        w = c(w, which(1:dim(tempXstar)[2] > (ns+1) &
                         !activeInd %in% unlist(activeZ0010[1:(wh01-1)])))
      }
      
      if (length(w) == 0) {
        p00_1[1] = log(tau[h]^(length(wv)+1) * (1 - tau[h])^2)
      } else if (length(w) == 1) {
        tempXstar2 = tempXstar[,w]
        
        SigmaB = sigB
        muB = rep(0, 1)
        muBeta = solve((t(tempXstar2) %*% tempXstar2)/
                         sigmaP + solve(SigmaB)) %*%
          ((t(tempXstar2) %*% tempY)/
             sigmaP + solve(SigmaB) %*% muB)
        covBeta = solve((t(tempXstar2) %*% tempXstar2)/
                          sigmaP + solve(SigmaB))
        p00_1[1] = log(tau[h]^(length(wv)+1) * (1 - tau[h])^2) + 
          mvtnorm::dmvnorm(rep(0, length(muB)), mean=muB, sigma=as.matrix(SigmaB), log=TRUE) - 
          mvtnorm::dmvnorm(rep(0, length(muB)), mean=muBeta, sigma=as.matrix(covBeta), log=TRUE)
      } else {
        tempXstar2 = tempXstar[,w]
        
        SigmaB = diag(sigmaP*sigB, dim(tempXstar2)[2])
        muB = rep(0, dim(tempXstar2)[2])
        muBeta = solve((t(tempXstar2) %*% tempXstar2)/
                         sigmaP + solve(SigmaB)) %*%
          ((t(tempXstar2) %*% tempY)/
             sigmaP + solve(SigmaB) %*% muB)
        covBeta = solve((t(tempXstar2) %*% tempXstar2)/
                          sigmaP + solve(SigmaB))
        p00_1[1] = log(tau[h]^(length(wv)+1) * (1 - tau[h])^2) + 
          mvtnorm::dmvnorm(rep(0, length(muB)), mean=muB, sigma=as.matrix(SigmaB), log=TRUE) - 
          mvtnorm::dmvnorm(rep(0, length(muB)), mean=muBeta, sigma=as.matrix(covBeta), log=TRUE)
        
      }
    } else {
      
      for (iii in 1 : nrow(PossMat00_1)) {
        tempXstar = matrix(NA, n, (ns+1)^length(wv))
        nsMat = expand.grid(rep(list(1 : (ns + 1)), length(wv)))
        activeMat = expand.grid(rep(list(c(0,rep(1,ns))), length(wv)))
        activeInd = "0"
        for (ii in 1 : nrow(nsMat)) {
          if (sum(activeMat[ii,]) > 0) {
            activeInd = c(activeInd, paste(sort(wv[which(activeMat[ii,]==1)]), collapse="-"))
          }
          tempXstar[,ii] = Xstar[,wv[1], nsMat[ii,1]]
          for (jj in 2 : length(wv)) {
            tempXstar[,ii] = tempXstar[,ii]*Xstar[,wv[jj],nsMat[ii,jj]]
          }
        }
        
        if (iii == 1) {
          tempXstar2 = Xstar[,groups[1],]
          activeInd = c(activeInd, as.character(rep(groups[1], ns)))
          wv3Temp = groups[1]
        } else {
          wvTemp = which(PossMat00_1[iii,] == 1)
          wv3Temp = c(wv[wvTemp], groups[1])
          
          tempXstar2 = matrix(NA, n, (ns+1)^length(wv3Temp))
          nsMat = expand.grid(rep(list(1 : (ns + 1)), length(wv3Temp)))
          activeMat = expand.grid(rep(list(c(0,rep(1,ns))), length(wv3Temp)))
          for (ii in 1 : nrow(nsMat)) {
            if (sum(activeMat[ii,]) > 0) {
              activeInd = c(activeInd, paste(sort(wv3Temp[which(activeMat[ii,]==1)]), collapse="-"))
            }
            tempXstar2[,ii] = Xstar[,wv3Temp[1], nsMat[ii,1]]
            for (jj in 2 : length(wv3Temp)) {
              tempXstar2[,ii] = tempXstar2[,ii]*Xstar[,wv3Temp[jj],nsMat[ii,jj]]
            }
          }
        }
        
        tempXstar = cbind(rep(1, n), scale(cbind(tempXstar[,-1], tempXstar2[,-1])))
        
        tempZeta0010 = tempZeta
        tempZeta0010[,h] = 0
        tempZeta0010[wv,h] = 1
        tempZeta0010[wv3Temp,wh01] = 1
        activeZ0010 = ActivePred(tempZeta0010, k=k)
        
        if (h == 1) {
          w = 2:((ns+1)^length(wv))
          w = c(w, which(!activeInd %in% unlist(activeZ0010[1:(wh01-1)]) &
                           1:dim(tempXstar)[2] > ((ns+1)^length(wv))))
        } else if (wh01 == 1) {
          w = (((ns+1)^length(wv)) + 1):dim(tempXstar)[2]
          w = c(w, which(!activeInd %in% unlist(activeZ0010[1:(h-1)]) &
                           1:dim(tempXstar)[2] <= ((ns+1)^length(wv))))
        } else if (h < wh01) {
          w = which(!activeInd[1:(((ns+1)^length(wv)))] %in% unlist(activeZ0010[1:(h-1)]))
          w = c(w, which(!activeInd %in% unlist(activeZ0010[1:(wh01-1)])))
        } else {
          w = which(!activeInd[1:(((ns+1)^length(wv)))] %in% unlist(activeZ0010[1:(h-1)]))
          w = c(w, which(1:dim(tempXstar)[2] > (((ns+1)^length(wv))) &
                           !activeInd %in% unlist(activeZ0010[1:(wh01-1)])))
        }
        
        if (length(w) == 0) {
          p00_1[iii] = log((tau[h]^(length(wv) + length(wv3Temp))) * 
                             (1 - tau[h])^2)
        } else if (length(w) == 1) {
          tempXstar2 = tempXstar[,w]
          
          SigmaB = sigB
          muB = rep(0, 1)
          muBeta = solve((t(tempXstar2) %*% tempXstar2)/
                           sigmaP + solve(SigmaB)) %*%
            ((t(tempXstar2) %*% tempY)/
               sigmaP + solve(SigmaB) %*% muB)
          covBeta = solve((t(tempXstar2) %*% tempXstar2)/
                            sigmaP + solve(SigmaB))
          p00_1[iii] = log((tau[h]^(length(wv) + length(wv3Temp))) * 
                             (1 - tau[h])^2) + 
            mvtnorm::dmvnorm(rep(0, length(muB)), mean=muB, sigma=as.matrix(SigmaB), log=TRUE) - 
            mvtnorm::dmvnorm(rep(0, length(muB)), mean=muBeta, sigma=as.matrix(covBeta), log=TRUE)
        } else {
          tempXstar2 = tempXstar[,w]
          
          SigmaB = diag(sigmaP*sigB, dim(tempXstar2)[2])
          muB = rep(0, dim(tempXstar2)[2])
          muBeta = solve((t(tempXstar2) %*% tempXstar2)/
                           sigmaP + solve(SigmaB)) %*%
            ((t(tempXstar2) %*% tempY)/
               sigmaP + solve(SigmaB) %*% muB)
          covBeta = solve((t(tempXstar2) %*% tempXstar2)/
                            sigmaP + solve(SigmaB))
          p00_1[iii] = log((tau[h]^(length(wv) + length(wv3Temp))) * 
                             (1 - tau[h])^2) + 
            mvtnorm::dmvnorm(rep(0, length(muB)), mean=muB, sigma=as.matrix(SigmaB), log=TRUE) - 
            mvtnorm::dmvnorm(rep(0, length(muB)), mean=muBeta, sigma=as.matrix(covBeta), log=TRUE)
          
        }
        
      }
      
    }
    
    
  }
  
  ## Now calculate probability for reduced model + effect with group 2
  if ((numZero > 0 & sum(tempZeta01[,h]) > 1)) {
    
    PossMatTemp = expand.grid(rep(list(c(0,1)), length(wv)))
    PossMat00_2 = PossMatTemp[-nrow(PossMatTemp),]
    
    wh01 = which(apply(tempZeta, 2, sum) == 0 & 1:k != h)[1]
    
    if (length(wv) == 1) {
      
      tempZeta0001 = tempZeta
      tempZeta0001[,h] = 0
      tempZeta0001[wv,h] = 1
      tempZeta0001[groups[2],wh01] = 1
      activeZ0001 = ActivePred(tempZeta0001, k=k)
      
      tempXstar = cbind(Xstar[,wv,], Xstar[,groups[2],-1])
      activeInd = as.character(c(0, rep(wv, ns), rep(groups[2], ns)))
      
      if (h == 1) {
        w = 2:(ns+1)
        if (as.character(groups[2]) %in% unlist(activeZ0001[1:(wh01-1)])) {
          2:dim(tempXstar)[2]
        }
      } else if (wh01 == 1) {
        w = (ns+2):dim(tempXstar)[2]
        if (as.character(wv) %in% unlist(activeZ0001[1:(h-1)])) {
          2:dim(tempXstar)[2]
        }
      } else if (h < wh01) {
        w = which(!activeInd[1:(ns+1)] %in% unlist(activeZ0001[1:(h-1)]))
        w = c(w, which(!activeInd %in% unlist(activeZ0001[1:(wh01-1)])))
      } else {
        w = which(!activeInd[1:(ns+1)] %in% unlist(activeZ0001[1:(h-1)]))
        w = c(w, which(1:dim(tempXstar)[2] > (ns+1) &
                         !activeInd %in% unlist(activeZ0001[1:(wh01-1)])))
      }
      
      if (length(w) == 0) {
        p00_2[1] = log(tau[h]^(length(wv)+1) * (1 - tau[h])^2)
      } else if (length(w) == 1) {
        tempXstar2 = tempXstar[,w]
        
        SigmaB = sigB
        muB = rep(0, 1)
        muBeta = solve((t(tempXstar2) %*% tempXstar2)/
                         sigmaP + solve(SigmaB)) %*%
          ((t(tempXstar2) %*% tempY)/
             sigmaP + solve(SigmaB) %*% muB)
        covBeta = solve((t(tempXstar2) %*% tempXstar2)/
                          sigmaP + solve(SigmaB))
        p00_2[1] = log(tau[h]^(length(wv)+1) * (1 - tau[h])^2) + 
          mvtnorm::dmvnorm(rep(0, length(muB)), mean=muB, sigma=as.matrix(SigmaB), log=TRUE) - 
          mvtnorm::dmvnorm(rep(0, length(muB)), mean=muBeta, sigma=as.matrix(covBeta), log=TRUE)
      } else {
        tempXstar2 = tempXstar[,w]
        
        SigmaB = diag(sigmaP*sigB, dim(tempXstar2)[2])
        muB = rep(0, dim(tempXstar2)[2])
        muBeta = solve((t(tempXstar2) %*% tempXstar2)/
                         sigmaP + solve(SigmaB)) %*%
          ((t(tempXstar2) %*% tempY)/
             sigmaP + solve(SigmaB) %*% muB)
        covBeta = solve((t(tempXstar2) %*% tempXstar2)/
                          sigmaP + solve(SigmaB))
        p00_2[1] = log(tau[h]^(length(wv)+1) * (1 - tau[h])^2) + 
          mvtnorm::dmvnorm(rep(0, length(muB)), mean=muB, sigma=as.matrix(SigmaB), log=TRUE) - 
          mvtnorm::dmvnorm(rep(0, length(muB)), mean=muBeta, sigma=as.matrix(covBeta), log=TRUE)
        
      }
    } else {
      
      for (iii in 1 : nrow(PossMat00_2)) {
        tempXstar = matrix(NA, n, (ns+1)^length(wv))
        nsMat = expand.grid(rep(list(1 : (ns + 1)), length(wv)))
        activeMat = expand.grid(rep(list(c(0,rep(1,ns))), length(wv)))
        activeInd = "0"
        for (ii in 1 : nrow(nsMat)) {
          if (sum(activeMat[ii,]) > 0) {
            activeInd = c(activeInd, paste(sort(wv[which(activeMat[ii,]==1)]), collapse="-"))
          }
          tempXstar[,ii] = Xstar[,wv[1], nsMat[ii,1]]
          for (jj in 2 : length(wv)) {
            tempXstar[,ii] = tempXstar[,ii]*Xstar[,wv[jj],nsMat[ii,jj]]
          }
        }
        
        if (iii == 1) {
          tempXstar2 = Xstar[,groups[2],]
          activeInd = c(activeInd, as.character(rep(groups[2], ns)))
          wv3Temp = groups[2]
        } else {
          wvTemp = which(PossMat00_2[iii,] == 1)
          wv3Temp = c(wv[wvTemp], groups[2])
          
          tempXstar2 = matrix(NA, n, (ns+1)^length(wv3Temp))
          nsMat = expand.grid(rep(list(1 : (ns + 1)), length(wv3Temp)))
          activeMat = expand.grid(rep(list(c(0,rep(1,ns))), length(wv3Temp)))
          for (ii in 1 : nrow(nsMat)) {
            if (sum(activeMat[ii,]) > 0) {
              activeInd = c(activeInd, paste(sort(wv3Temp[which(activeMat[ii,]==1)]), collapse="-"))
            }
            tempXstar2[,ii] = Xstar[,wv3Temp[1], nsMat[ii,1]]
            for (jj in 2 : length(wv3Temp)) {
              tempXstar2[,ii] = tempXstar2[,ii]*Xstar[,wv3Temp[jj],nsMat[ii,jj]]
            }
          }
        }
        
        tempXstar = cbind(rep(1, n), scale(cbind(tempXstar[,-1], tempXstar2[,-1])))
        
        tempZeta0001 = tempZeta
        tempZeta0001[,h] = 0
        tempZeta0001[wv,h] = 1
        tempZeta0001[wv3Temp,wh01] = 1
        activeZ0001 = ActivePred(tempZeta0001, k=k)
        
        if (h == 1) {
          w = 2:((ns+1)^length(wv))
          w = c(w, which(!activeInd %in% unlist(activeZ0001[1:(wh01-1)]) &
                           1:dim(tempXstar)[2] > ((ns+1)^length(wv))))
        } else if (wh01 == 1) {
          w = (((ns+1)^length(wv)) + 1):dim(tempXstar)[2]
          w = c(w, which(!activeInd %in% unlist(activeZ0001[1:(h-1)]) &
                           1:dim(tempXstar)[2] <= ((ns+1)^length(wv))))
        } else if (h < wh01) {
          w = which(!activeInd[1:(((ns+1)^length(wv)))] %in% unlist(activeZ0001[1:(h-1)]))
          w = c(w, which(!activeInd %in% unlist(activeZ0001[1:(wh01-1)])))
        } else {
          w = which(!activeInd[1:(((ns+1)^length(wv)))] %in% unlist(activeZ0001[1:(h-1)]))
          w = c(w, which(1:dim(tempXstar)[2] > (((ns+1)^length(wv))) &
                           !activeInd %in% unlist(activeZ0001[1:(wh01-1)])))
        }
        
        if (length(w) == 0) {
          p00_2[iii] = log((tau[h]^(length(wv) + length(wv3Temp))) * 
                             (1 - tau[h])^2)
        } else if (length(w) == 1) {
          tempXstar2 = tempXstar[,w]
          
          SigmaB = sigB
          muB = rep(0, 1)
          muBeta = solve((t(tempXstar2) %*% tempXstar2)/
                           sigmaP + solve(SigmaB)) %*%
            ((t(tempXstar2) %*% tempY)/
               sigmaP + solve(SigmaB) %*% muB)
          covBeta = solve((t(tempXstar2) %*% tempXstar2)/
                            sigmaP + solve(SigmaB))
          p00_2[iii] = log((tau[h]^(length(wv) + length(wv3Temp))) * 
                             (1 - tau[h])^2) + 
            mvtnorm::dmvnorm(rep(0, length(muB)), mean=muB, sigma=as.matrix(SigmaB), log=TRUE) - 
            mvtnorm::dmvnorm(rep(0, length(muB)), mean=muBeta, sigma=as.matrix(covBeta), log=TRUE)
        } else {
          tempXstar2 = tempXstar[,w]
          
          SigmaB = diag(sigmaP*sigB, dim(tempXstar2)[2])
          muB = rep(0, dim(tempXstar2)[2])
          muBeta = solve((t(tempXstar2) %*% tempXstar2)/
                           sigmaP + solve(SigmaB)) %*%
            ((t(tempXstar2) %*% tempY)/
               sigmaP + solve(SigmaB) %*% muB)
          covBeta = solve((t(tempXstar2) %*% tempXstar2)/
                            sigmaP + solve(SigmaB))
          p00_2[iii] = log((tau[h]^(length(wv) + length(wv3Temp))) * 
                             (1 - tau[h])^2) + 
            mvtnorm::dmvnorm(rep(0, length(muB)), mean=muB, sigma=as.matrix(SigmaB), log=TRUE) - 
            mvtnorm::dmvnorm(rep(0, length(muB)), mean=muBeta, sigma=as.matrix(covBeta), log=TRUE)
          
        }
        
      }
      
    }
    
    
  }
  
  ## Now calculate probability for reduced model + int effect with group 1 + 2
  if ((numZero > 0 & sum(tempZeta11[,h]) > 2)) {
    PossMatTemp = expand.grid(rep(list(c(0,1)), length(wv)))
    PossMat00_12 = PossMatTemp[-nrow(PossMatTemp),]
    
    wh01 = which(apply(tempZeta, 2, sum) == 0 & 1:k != h)[1]
    
    if (length(wv) == 1) {
      
      tempZeta0011 = tempZeta
      tempZeta0011[,h] = 0
      tempZeta0011[wv,h] = 1
      tempZeta0011[groups,wh01] = 1
      activeZ0011 = ActivePred(tempZeta0011, k=k)
      activeInd = as.character(c(0, rep(wv, ns)))
      
      tempXstar = Xstar[,wv,]
      
      wv3Temp = c(groups)
      
      tempXstar2 = matrix(NA, n, (ns+1)^length(wv3Temp))
      nsMat = expand.grid(rep(list(1 : (ns + 1)), length(wv3Temp)))
      activeMat = expand.grid(rep(list(c(0,rep(1,ns))), length(wv3Temp)))
      for (ii in 1 : nrow(nsMat)) {
        if (sum(activeMat[ii,]) > 0) {
          activeInd = c(activeInd, paste(sort(wv3Temp[which(activeMat[ii,]==1)]), collapse="-"))
        }
        tempXstar2[,ii] = Xstar[,wv3Temp[1], nsMat[ii,1]]
        for (jj in 2 : length(wv3Temp)) {
          tempXstar2[,ii] = tempXstar2[,ii]*Xstar[,wv3Temp[jj],nsMat[ii,jj]]
        }
      }
      
      tempXstar = cbind(rep(1, n), scale(cbind(tempXstar[,-1], tempXstar2[,-1])))
      
      if (h == 1) {
        w = 2:(ns+1)
        w = c(w, which(!activeInd %in% unlist(activeZ0011[1:(wh01-1)])))
      } else if (wh01 == 1) {
        w = (ns+2):dim(tempXstar)[2]
        w = c(w, which(!activeInd %in% unlist(activeZ0011[1:(h-1)])))
      } else if (h < wh01) {
        w = which(!activeInd[1:(ns+1)] %in% unlist(activeZ0011[1:(h-1)]))
        w = c(w, which(!activeInd %in% unlist(activeZ0011[1:(wh01-1)])))
      } else {
        w = which(!activeInd[1:(ns+1)] %in% unlist(activeZ0011[1:(h-1)]))
        w = c(w, which(1:dim(tempXstar)[2] > (ns+1) &
                         !activeInd %in% unlist(activeZ0011[1:(wh01-1)])))
      }
      
      if (length(w) == 0) {
        p00_12[1] = log(tau[h]^(length(wv)+2) * (1 - tau[h])^2)
      } else if (length(w) == 1) {
        tempXstar2 = tempXstar[,w]
        
        SigmaB = sigB
        muB = rep(0, 1)
        muBeta = solve((t(tempXstar2) %*% tempXstar2)/
                         sigmaP + solve(SigmaB)) %*%
          ((t(tempXstar2) %*% tempY)/
             sigmaP + solve(SigmaB) %*% muB)
        covBeta = solve((t(tempXstar2) %*% tempXstar2)/
                          sigmaP + solve(SigmaB))
        p00_12[1] = log(tau[h]^(length(wv)+2) * (1 - tau[h])^2) + 
          mvtnorm::dmvnorm(rep(0, length(muB)), mean=muB, sigma=as.matrix(SigmaB), log=TRUE) - 
          mvtnorm::dmvnorm(rep(0, length(muB)), mean=muBeta, sigma=as.matrix(covBeta), log=TRUE)
      } else {
        tempXstar2 = tempXstar[,w]
        
        SigmaB = diag(sigmaP*sigB, dim(tempXstar2)[2])
        muB = rep(0, dim(tempXstar2)[2])
        muBeta = solve((t(tempXstar2) %*% tempXstar2)/
                         sigmaP + solve(SigmaB)) %*%
          ((t(tempXstar2) %*% tempY)/
             sigmaP + solve(SigmaB) %*% muB)
        covBeta = solve((t(tempXstar2) %*% tempXstar2)/
                          sigmaP + solve(SigmaB))
        p00_12[1] = log(tau[h]^(length(wv)+2) * (1 - tau[h])^2) + 
          mvtnorm::dmvnorm(rep(0, length(muB)), mean=muB, sigma=as.matrix(SigmaB), log=TRUE) - 
          mvtnorm::dmvnorm(rep(0, length(muB)), mean=muBeta, sigma=as.matrix(covBeta), log=TRUE)
        
      }
    } else {
      
      for (iii in 1 : nrow(PossMat00_12)) {
        tempXstar = matrix(NA, n, (ns+1)^length(wv))
        nsMat = expand.grid(rep(list(1 : (ns + 1)), length(wv)))
        activeMat = expand.grid(rep(list(c(0,rep(1,ns))), length(wv)))
        activeInd = "0"
        for (ii in 1 : nrow(nsMat)) {
          if (sum(activeMat[ii,]) > 0) {
            activeInd = c(activeInd, paste(sort(wv[which(activeMat[ii,]==1)]), collapse="-"))
          }
          tempXstar[,ii] = Xstar[,wv[1], nsMat[ii,1]]
          for (jj in 2 : length(wv)) {
            tempXstar[,ii] = tempXstar[,ii]*Xstar[,wv[jj],nsMat[ii,jj]]
          }
        }
        
        if (iii == 1) {
          wv3Temp = groups
          
          tempXstar2 = matrix(NA, n, (ns+1)^length(wv3Temp))
          nsMat = expand.grid(rep(list(1 : (ns + 1)), length(wv3Temp)))
          activeMat = expand.grid(rep(list(c(0,rep(1,ns))), length(wv3Temp)))
          for (ii in 1 : nrow(nsMat)) {
            if (sum(activeMat[ii,]) > 0) {
              activeInd = c(activeInd, paste(sort(wv3Temp[which(activeMat[ii,]==1)]), collapse="-"))
            }
            tempXstar2[,ii] = Xstar[,wv3Temp[1], nsMat[ii,1]]
            for (jj in 2 : length(wv3Temp)) {
              tempXstar2[,ii] = tempXstar2[,ii]*Xstar[,wv3Temp[jj],nsMat[ii,jj]]
            }
          }
        } else {
          wvTemp = which(PossMat00_12[iii,] == 1)
          wv3Temp = c(wv[wvTemp], groups)
          
          tempXstar2 = matrix(NA, n, (ns+1)^length(wv3Temp))
          nsMat = expand.grid(rep(list(1 : (ns + 1)), length(wv3Temp)))
          activeMat = expand.grid(rep(list(c(0,rep(1,ns))), length(wv3Temp)))
          for (ii in 1 : nrow(nsMat)) {
            if (sum(activeMat[ii,]) > 0) {
              activeInd = c(activeInd, paste(sort(wv3Temp[which(activeMat[ii,]==1)]), collapse="-"))
            }
            tempXstar2[,ii] = Xstar[,wv3Temp[1], nsMat[ii,1]]
            for (jj in 2 : length(wv3Temp)) {
              tempXstar2[,ii] = tempXstar2[,ii]*Xstar[,wv3Temp[jj],nsMat[ii,jj]]
            }
          }
        }
        
        tempXstar = cbind(rep(1, n), scale(cbind(tempXstar[,-1], tempXstar2[,-1])))
        
        tempZeta0011 = tempZeta
        tempZeta0011[,h] = 0
        tempZeta0011[wv,h] = 1
        tempZeta0011[wv3Temp,wh01] = 1
        activeZ0011 = ActivePred(tempZeta0011, k=k)
        
        if (h == 1) {
          w = 2:((ns+1)^length(wv))
          w = c(w, which(!activeInd %in% unlist(activeZ0011[1:(wh01-1)]) &
                           1:dim(tempXstar)[2] > ((ns+1)^length(wv))))
        } else if (wh01 == 1) {
          w = (((ns+1)^length(wv)) + 1):dim(tempXstar)[2]
          w = c(w, which(!activeInd %in% unlist(activeZ0011[1:(h-1)]) &
                           1:dim(tempXstar)[2] <= ((ns+1)^length(wv))))
        } else if (h < wh01) {
          w = which(!activeInd[1:(((ns+1)^length(wv)))] %in% unlist(activeZ0011[1:(h-1)]))
          w = c(w, which(!activeInd %in% unlist(activeZ0011[1:(wh01-1)])))
        } else {
          w = which(!activeInd[1:(((ns+1)^length(wv)))] %in% unlist(activeZ0011[1:(h-1)]))
          w = c(w, which(1:dim(tempXstar)[2] > (((ns+1)^length(wv))) &
                           !activeInd %in% unlist(activeZ0011[1:(wh01-1)])))
        }
        
        if (length(w) == 0) {
          p00_12[iii] = log((tau[h]^(length(wv) + length(wv3Temp))) * 
                              (1 - tau[h])^2)
        } else if (length(w) == 1) {
          tempXstar2 = tempXstar[,w]
          
          SigmaB = sigB
          muB = rep(0, 1)
          muBeta = solve((t(tempXstar2) %*% tempXstar2)/
                           sigmaP + solve(SigmaB)) %*%
            ((t(tempXstar2) %*% tempY)/
               sigmaP + solve(SigmaB) %*% muB)
          covBeta = solve((t(tempXstar2) %*% tempXstar2)/
                            sigmaP + solve(SigmaB))
          p00_12[iii] = log((tau[h]^(length(wv) + length(wv3Temp))) * 
                              (1 - tau[h])^2) + 
            mvtnorm::dmvnorm(rep(0, length(muB)), mean=muB, sigma=as.matrix(SigmaB), log=TRUE) - 
            mvtnorm::dmvnorm(rep(0, length(muB)), mean=muBeta, sigma=as.matrix(covBeta), log=TRUE)
        } else {
          tempXstar2 = tempXstar[,w]
          
          SigmaB = diag(sigmaP*sigB, dim(tempXstar2)[2])
          muB = rep(0, dim(tempXstar2)[2])
          muBeta = solve((t(tempXstar2) %*% tempXstar2)/
                           sigmaP + solve(SigmaB)) %*%
            ((t(tempXstar2) %*% tempY)/
               sigmaP + solve(SigmaB) %*% muB)
          covBeta = solve((t(tempXstar2) %*% tempXstar2)/
                            sigmaP + solve(SigmaB))
          p00_12[iii] = log((tau[h]^(length(wv) + length(wv3Temp))) * 
                              (1 - tau[h])^2) + 
            mvtnorm::dmvnorm(rep(0, length(muB)), mean=muB, sigma=as.matrix(SigmaB), log=TRUE) - 
            mvtnorm::dmvnorm(rep(0, length(muB)), mean=muBeta, sigma=as.matrix(covBeta), log=TRUE)
          
        }
        
      }
      
    }
    
    
    
  }
  
  
  ## Now calculate probability for reduced model with interaction of group 1
  ## + effect with group 2
  if (numZero > 0) {
    wv2 = c(wv, groups[1])
    PossMatTemp = expand.grid(rep(list(c(0,1)), length(wv2)))
    PossMat10_2 = PossMatTemp[-nrow(PossMatTemp),]
    
    wh01 = which(apply(tempZeta, 2, sum) == 0 & 1:k != h)[1]
    
    if (length(wv2) == 1) {
      
      tempZeta1001 = tempZeta
      tempZeta1001[,h] = 0
      tempZeta1001[wv2,h] = 1
      tempZeta1001[groups[2],wh01] = 1
      activeZ1001 = ActivePred(tempZeta1001, k=k)
      
      tempXstar = cbind(Xstar[,wv2,], Xstar[,groups[2],-1])
      activeInd = as.character(c(0, rep(wv2, ns), rep(groups[2], ns)))
      
      if (h == 1) {
        w = 2:(ns+1)
        if (as.character(groups[2]) %in% unlist(activeZ1001[1:(wh01-1)])) {
          2:dim(tempXstar)[2]
        }
      } else if (wh01 == 1) {
        w = (ns+2):dim(tempXstar)[2]
        if (as.character(wv2) %in% unlist(activeZ1001[1:(h-1)])) {
          2:dim(tempXstar)[2]
        }
      } else if (h < wh01) {
        w = which(!activeInd[1:(ns+1)] %in% unlist(activeZ1001[1:(h-1)]))
        w = c(w, which(!activeInd %in% unlist(activeZ1001[1:(wh01-1)])))
      } else {
        w = which(!activeInd[1:(ns+1)] %in% unlist(activeZ1001[1:(h-1)]))
        w = c(w, which(1:dim(tempXstar)[2] > (ns+1) &
                         !activeInd %in% unlist(activeZ1001[1:(wh01-1)])))
      }
      
      if (length(w) == 0) {
        p10_2[1] = log(tau[h]^(length(wv2)+1) * (1 - tau[h]))
      } else if (length(w) == 1) {
        tempXstar2 = tempXstar[,w]
        
        SigmaB = sigB
        muB = rep(0, 1)
        muBeta = solve((t(tempXstar2) %*% tempXstar2)/
                         sigmaP + solve(SigmaB)) %*%
          ((t(tempXstar2) %*% tempY)/
             sigmaP + solve(SigmaB) %*% muB)
        covBeta = solve((t(tempXstar2) %*% tempXstar2)/
                          sigmaP + solve(SigmaB))
        p10_2[1] = log(tau[h]^(length(wv2)+1) * (1 - tau[h])) + 
          mvtnorm::dmvnorm(rep(0, length(muB)), mean=muB, sigma=as.matrix(SigmaB), log=TRUE) - 
          mvtnorm::dmvnorm(rep(0, length(muB)), mean=muBeta, sigma=as.matrix(covBeta), log=TRUE)
      } else {
        tempXstar2 = tempXstar[,w]
        
        SigmaB = diag(sigmaP*sigB, dim(tempXstar2)[2])
        muB = rep(0, dim(tempXstar2)[2])
        muBeta = solve((t(tempXstar2) %*% tempXstar2)/
                         sigmaP + solve(SigmaB)) %*%
          ((t(tempXstar2) %*% tempY)/
             sigmaP + solve(SigmaB) %*% muB)
        covBeta = solve((t(tempXstar2) %*% tempXstar2)/
                          sigmaP + solve(SigmaB))
        p10_2[1] = log(tau[h]^(length(wv2)+1) * (1 - tau[h])) + 
          mvtnorm::dmvnorm(rep(0, length(muB)), mean=muB, sigma=as.matrix(SigmaB), log=TRUE) - 
          mvtnorm::dmvnorm(rep(0, length(muB)), mean=muBeta, sigma=as.matrix(covBeta), log=TRUE)
        
      }
    } else {
      
      for (iii in 1 : nrow(PossMat10_2)) {
        tempXstar = matrix(NA, n, (ns+1)^length(wv2))
        nsMat = expand.grid(rep(list(1 : (ns + 1)), length(wv2)))
        activeMat = expand.grid(rep(list(c(0,rep(1,ns))), length(wv2)))
        activeInd = "0"
        for (ii in 1 : nrow(nsMat)) {
          if (sum(activeMat[ii,]) > 0) {
            activeInd = c(activeInd, paste(sort(wv2[which(activeMat[ii,]==1)]), collapse="-"))
          }
          tempXstar[,ii] = Xstar[,wv2[1], nsMat[ii,1]]
          for (jj in 2 : length(wv2)) {
            tempXstar[,ii] = tempXstar[,ii]*Xstar[,wv2[jj],nsMat[ii,jj]]
          }
        }
        
        if (iii == 1) {
          tempXstar2 = Xstar[,groups[2],]
          activeInd = c(activeInd, as.character(rep(groups[2], ns)))
          wv3Temp = groups[2]
        } else {
          wvTemp = which(PossMat10_2[iii,] == 1)
          wv3Temp = c(wv2[wvTemp], groups[2])
          
          tempXstar2 = matrix(NA, n, (ns+1)^length(wv3Temp))
          nsMat = expand.grid(rep(list(1 : (ns + 1)), length(wv3Temp)))
          activeMat = expand.grid(rep(list(c(0,rep(1,ns))), length(wv3Temp)))
          for (ii in 1 : nrow(nsMat)) {
            if (sum(activeMat[ii,]) > 0) {
              activeInd = c(activeInd, paste(sort(wv3Temp[which(activeMat[ii,]==1)]), collapse="-"))
            }
            tempXstar2[,ii] = Xstar[,wv3Temp[1], nsMat[ii,1]]
            for (jj in 2 : length(wv3Temp)) {
              tempXstar2[,ii] = tempXstar2[,ii]*Xstar[,wv3Temp[jj],nsMat[ii,jj]]
            }
          }
        }
        
        tempXstar = cbind(rep(1, n), scale(cbind(tempXstar[,-1], tempXstar2[,-1])))
        
        tempZeta1001 = tempZeta
        tempZeta1001[,h] = 0
        tempZeta1001[wv2,h] = 1
        tempZeta1001[wv3Temp,wh01] = 1
        activeZ1001 = ActivePred(tempZeta1001, k=k)
        
        if (h == 1) {
          w = 2:((ns+1)^length(wv2))
          w = c(w, which(!activeInd %in% unlist(activeZ1001[1:(wh01-1)]) &
                           1:dim(tempXstar)[2] > ((ns+1)^length(wv2))))
        } else if (wh01 == 1) {
          w = (((ns+1)^length(wv2)) + 1):dim(tempXstar)[2]
          w = c(w, which(!activeInd %in% unlist(activeZ1001[1:(h-1)]) &
                           1:dim(tempXstar)[2] <= ((ns+1)^length(wv2))))
        } else if (h < wh01) {
          w = which(!activeInd[1:(((ns+1)^length(wv2)))] %in% unlist(activeZ1001[1:(h-1)]))
          w = c(w, which(!activeInd %in% unlist(activeZ1001[1:(wh01-1)])))
        } else {
          w = which(!activeInd[1:(((ns+1)^length(wv2)))] %in% unlist(activeZ1001[1:(h-1)]))
          w = c(w, which(1:dim(tempXstar)[2] > (((ns+1)^length(wv2))) &
                           !activeInd %in% unlist(activeZ1001[1:(wh01-1)])))
        }
        
        if (length(w) == 0) {
          p10_2[iii] = log((tau[h]^(length(wv2) + length(wv3Temp))) * 
                             (1 - tau[h]))
        } else if (length(w) == 1) {
          tempXstar2 = tempXstar[,w]
          
          SigmaB = sigB
          muB = rep(0, 1)
          muBeta = solve((t(tempXstar2) %*% tempXstar2)/
                           sigmaP + solve(SigmaB)) %*%
            ((t(tempXstar2) %*% tempY)/
               sigmaP + solve(SigmaB) %*% muB)
          covBeta = solve((t(tempXstar2) %*% tempXstar2)/
                            sigmaP + solve(SigmaB))
          p10_2[iii] = log((tau[h]^(length(wv2) + length(wv3Temp))) * 
                             (1 - tau[h])) + 
            mvtnorm::dmvnorm(rep(0, length(muB)), mean=muB, sigma=as.matrix(SigmaB), log=TRUE) - 
            mvtnorm::dmvnorm(rep(0, length(muB)), mean=muBeta, sigma=as.matrix(covBeta), log=TRUE)
        } else {
          tempXstar2 = tempXstar[,w]
          
          SigmaB = diag(sigmaP*sigB, dim(tempXstar2)[2])
          muB = rep(0, dim(tempXstar2)[2])
          muBeta = solve((t(tempXstar2) %*% tempXstar2)/
                           sigmaP + solve(SigmaB)) %*%
            ((t(tempXstar2) %*% tempY)/
               sigmaP + solve(SigmaB) %*% muB)
          covBeta = solve((t(tempXstar2) %*% tempXstar2)/
                            sigmaP + solve(SigmaB))
          p10_2[iii] = log((tau[h]^(length(wv2) + length(wv3Temp))) * 
                             (1 - tau[h])) + 
            mvtnorm::dmvnorm(rep(0, length(muB)), mean=muB, sigma=as.matrix(SigmaB), log=TRUE) - 
            mvtnorm::dmvnorm(rep(0, length(muB)), mean=muBeta, sigma=as.matrix(covBeta), log=TRUE)
          
        }
        
      }
      
    }
    
    
    
  }
  
  ## Now calculate probability for reduced model with interaction of group 2
  ## + effect with group 1
  if (numZero > 0) {
    
    wv2 = c(wv, groups[2])
    PossMatTemp = expand.grid(rep(list(c(0,1)), length(wv2)))
    PossMat01_1 = PossMatTemp[-nrow(PossMatTemp),]
    
    wh01 = which(apply(tempZeta, 2, sum) == 0 & 1:k != h)[1]
    
    if (length(wv2) == 1) {
      
      tempZeta0110 = tempZeta
      tempZeta0110[,h] = 0
      tempZeta0110[wv2,h] = 1
      tempZeta0110[groups[1],wh01] = 1
      activeZ0110 = ActivePred(tempZeta0110, k=k)
      
      tempXstar = cbind(Xstar[,wv2,], Xstar[,groups[1],-1])
      activeInd = as.character(c(0, rep(wv2, ns), rep(groups[1], ns)))
      
      if (h == 1) {
        w = 2:(ns+1)
        if (as.character(groups[1]) %in% unlist(activeZ0110[1:(wh01-1)])) {
          2:dim(tempXstar)[2]
        }
      } else if (wh01 == 1) {
        w = (ns+2):dim(tempXstar)[2]
        if (as.character(wv2) %in% unlist(activeZ0110[1:(h-1)])) {
          2:dim(tempXstar)[2]
        }
      } else if (h < wh01) {
        w = which(!activeInd[1:(ns+1)] %in% unlist(activeZ0110[1:(h-1)]))
        w = c(w, which(!activeInd %in% unlist(activeZ0110[1:(wh01-1)])))
      } else {
        w = which(!activeInd[1:(ns+1)] %in% unlist(activeZ0110[1:(h-1)]))
        w = c(w, which(1:dim(tempXstar)[2] > (ns+1) &
                         !activeInd %in% unlist(activeZ0110[1:(wh01-1)])))
      }
      
      if (length(w) == 0) {
        p01_1[1] = log(tau[h]^(length(wv2)+1) * (1 - tau[h]))
      } else if (length(w) == 1) {
        tempXstar2 = tempXstar[,w]
        
        SigmaB = sigB
        muB = rep(0, 1)
        muBeta = solve((t(tempXstar2) %*% tempXstar2)/
                         sigmaP + solve(SigmaB)) %*%
          ((t(tempXstar2) %*% tempY)/
             sigmaP + solve(SigmaB) %*% muB)
        covBeta = solve((t(tempXstar2) %*% tempXstar2)/
                          sigmaP + solve(SigmaB))
        p01_1[1] = log(tau[h]^(length(wv2)+1) * (1 - tau[h])) + 
          mvtnorm::dmvnorm(rep(0, length(muB)), mean=muB, sigma=as.matrix(SigmaB), log=TRUE) - 
          mvtnorm::dmvnorm(rep(0, length(muB)), mean=muBeta, sigma=as.matrix(covBeta), log=TRUE)
      } else {
        tempXstar2 = tempXstar[,w]
        
        SigmaB = diag(sigmaP*sigB, dim(tempXstar2)[2])
        muB = rep(0, dim(tempXstar2)[2])
        muBeta = solve((t(tempXstar2) %*% tempXstar2)/
                         sigmaP + solve(SigmaB)) %*%
          ((t(tempXstar2) %*% tempY)/
             sigmaP + solve(SigmaB) %*% muB)
        covBeta = solve((t(tempXstar2) %*% tempXstar2)/
                          sigmaP + solve(SigmaB))
        p01_1[1] = log(tau[h]^(length(wv2)+1) * (1 - tau[h])) + 
          mvtnorm::dmvnorm(rep(0, length(muB)), mean=muB, sigma=as.matrix(SigmaB), log=TRUE) - 
          mvtnorm::dmvnorm(rep(0, length(muB)), mean=muBeta, sigma=as.matrix(covBeta), log=TRUE)
        
      }
    } else {
      
      for (iii in 1 : nrow(PossMat01_1)) {
        tempXstar = matrix(NA, n, (ns+1)^length(wv2))
        nsMat = expand.grid(rep(list(1 : (ns + 1)), length(wv2)))
        activeMat = expand.grid(rep(list(c(0,rep(1,ns))), length(wv2)))
        activeInd = "0"
        for (ii in 1 : nrow(nsMat)) {
          if (sum(activeMat[ii,]) > 0) {
            activeInd = c(activeInd, paste(sort(wv2[which(activeMat[ii,]==1)]), collapse="-"))
          }
          tempXstar[,ii] = Xstar[,wv2[1], nsMat[ii,1]]
          for (jj in 2 : length(wv2)) {
            tempXstar[,ii] = tempXstar[,ii]*Xstar[,wv2[jj],nsMat[ii,jj]]
          }
        }
        
        if (iii == 1) {
          tempXstar2 = Xstar[,groups[1],]
          activeInd = c(activeInd, as.character(rep(groups[1], ns)))
          wv3Temp = groups[1]
        } else {
          wvTemp = which(PossMat01_1[iii,] == 1)
          wv3Temp = c(wv2[wvTemp], groups[1])
          
          tempXstar2 = matrix(NA, n, (ns+1)^length(wv3Temp))
          nsMat = expand.grid(rep(list(1 : (ns + 1)), length(wv3Temp)))
          activeMat = expand.grid(rep(list(c(0,rep(1,ns))), length(wv3Temp)))
          for (ii in 1 : nrow(nsMat)) {
            if (sum(activeMat[ii,]) > 0) {
              activeInd = c(activeInd, paste(sort(wv3Temp[which(activeMat[ii,]==1)]), collapse="-"))
            }
            tempXstar2[,ii] = Xstar[,wv3Temp[1], nsMat[ii,1]]
            for (jj in 2 : length(wv3Temp)) {
              tempXstar2[,ii] = tempXstar2[,ii]*Xstar[,wv3Temp[jj],nsMat[ii,jj]]
            }
          }
        }
        
        tempXstar = cbind(rep(1, n), scale(cbind(tempXstar[,-1], tempXstar2[,-1])))
        
        tempZeta0110 = tempZeta
        tempZeta0110[,h] = 0
        tempZeta0110[wv2,h] = 1
        tempZeta0110[wv3Temp,wh01] = 1
        activeZ0110 = ActivePred(tempZeta0110, k=k)
        
        if (h == 1) {
          w = 2:((ns+1)^length(wv2))
          w = c(w, which(!activeInd %in% unlist(activeZ0110[1:(wh01-1)]) &
                           1:dim(tempXstar)[2] > ((ns+1)^length(wv2))))
        } else if (wh01 == 1) {
          w = (((ns+1)^length(wv2)) + 1):dim(tempXstar)[2]
          w = c(w, which(!activeInd %in% unlist(activeZ0110[1:(h-1)]) &
                           1:dim(tempXstar)[2] <= ((ns+1)^length(wv2))))
        } else if (h < wh01) {
          w = which(!activeInd[1:(((ns+1)^length(wv2)))] %in% unlist(activeZ0110[1:(h-1)]))
          w = c(w, which(!activeInd %in% unlist(activeZ0110[1:(wh01-1)])))
        } else {
          w = which(!activeInd[1:(((ns+1)^length(wv2)))] %in% unlist(activeZ0110[1:(h-1)]))
          w = c(w, which(1:dim(tempXstar)[2] > (((ns+1)^length(wv2))) &
                           !activeInd %in% unlist(activeZ0110[1:(wh01-1)])))
        }
        
        if (length(w) == 0) {
          p01_1[iii] = log((tau[h]^(length(wv2) + length(wv3Temp))) * 
                             (1 - tau[h]))
        } else if (length(w) == 1) {
          tempXstar2 = tempXstar[,w]
          
          SigmaB = sigB
          muB = rep(0, 1)
          muBeta = solve((t(tempXstar2) %*% tempXstar2)/
                           sigmaP + solve(SigmaB)) %*%
            ((t(tempXstar2) %*% tempY)/
               sigmaP + solve(SigmaB) %*% muB)
          covBeta = solve((t(tempXstar2) %*% tempXstar2)/
                            sigmaP + solve(SigmaB))
          p01_1[iii] = log((tau[h]^(length(wv2) + length(wv3Temp))) * 
                             (1 - tau[h])) + 
            mvtnorm::dmvnorm(rep(0, length(muB)), mean=muB, sigma=as.matrix(SigmaB), log=TRUE) - 
            mvtnorm::dmvnorm(rep(0, length(muB)), mean=muBeta, sigma=as.matrix(covBeta), log=TRUE)
        } else {
          tempXstar2 = tempXstar[,w]
          
          SigmaB = diag(sigmaP*sigB, dim(tempXstar2)[2])
          muB = rep(0, dim(tempXstar2)[2])
          muBeta = solve((t(tempXstar2) %*% tempXstar2)/
                           sigmaP + solve(SigmaB)) %*%
            ((t(tempXstar2) %*% tempY)/
               sigmaP + solve(SigmaB) %*% muB)
          covBeta = solve((t(tempXstar2) %*% tempXstar2)/
                            sigmaP + solve(SigmaB))
          p01_1[iii] = log((tau[h]^(length(wv2) + length(wv3Temp))) * 
                             (1 - tau[h])) + 
            mvtnorm::dmvnorm(rep(0, length(muB)), mean=muB, sigma=as.matrix(SigmaB), log=TRUE) - 
            mvtnorm::dmvnorm(rep(0, length(muB)), mean=muBeta, sigma=as.matrix(covBeta), log=TRUE)
          
        }
        
      }
      
    }
    
  }
  
  
  
  
  
  maxlog = max(c(p00, p10, p01, p11, p00_1, 
                 p00_2, p00_12, p10_2, p01_1))
  updated_p = exp(c(p00, p10, p01, p11, p00_1, 
                    p00_2, p00_12, p10_2, p01_1) - maxlog)
  samp = sample(c("00", "10", "01", "11",
                  paste("00_10-", 1:length(p00_1), sep=''),
                  paste("00_20-", 1:length(p00_2), sep=''),
                  paste("00_12-", 1:length(p00_12), sep=''),
                  paste("10_20-", 1:length(p10_2), sep=''),
                  paste("01_10-", 1:length(p01_1), sep='')), 1, p = updated_p)
  
  if (!samp %in% c("00", "10", "01", "11")) {
    samp1 = substr(samp, 1, 5)
    samp2 = as.numeric(substr(samp, 7, nchar(samp)))
  }
  
  
  if (samp == "00") {
    tempZeta = tempZeta00
  } else if (samp == "10") {
    tempZeta = tempZeta10
  } else if (samp == "01") {
    tempZeta = tempZeta01
  } else if (samp == "11") {
    tempZeta = tempZeta11
  } else if (samp1 == "00_10" & samp2 == 1) {
    tempZeta = tempZeta00
    tempZeta[c(groups[1]),wh01] = 1
  } else if (samp1 == "00_10" & samp2 > 1) {
    tempZeta = tempZeta00
    tempZeta[c(wv, groups[1]),wh01] = c(as.numeric(PossMat00_1[samp2,]), 1)
  } else if (samp1 == "00_20" & samp2 == 1) {
    tempZeta = tempZeta00
    tempZeta[c(groups[2]),wh01] = 1
  } else if (samp1 == "00_20" & samp2 > 1) {
    tempZeta = tempZeta00
    tempZeta[c(wv, groups[2]),wh01] = c(as.numeric(PossMat00_2[samp2,]), 1)
  } else if (samp1 == "00_12" & samp2 == 1) {
    tempZeta = tempZeta00
    tempZeta[c(groups),wh01] = 1
  } else if (samp1 == "00_12" & samp2 > 1) {
    tempZeta = tempZeta00
    tempZeta[c(wv, groups),wh01] = c(as.numeric(PossMat00_12[samp2,]), 1, 1)
  } else if (samp1 == "10_20" & samp2 == 1) {
    tempZeta = tempZeta10
    tempZeta[c(groups[2]),wh01] = 1
  } else if (samp1 == "10_20" & samp2 > 1) {
    tempZeta = tempZeta10
    tempZeta[c(wv, groups),wh01] = c(as.numeric(PossMat10_2[samp2,]), 1)
  } else if (samp1 == "01_10" & samp2 == 1) {
    tempZeta = tempZeta01
    tempZeta[c(groups[1]),wh01] = 1
  } else if (samp1 == "01_10" & samp2 > 1) {
    tempZeta = tempZeta01
    tempZeta[c(wv, groups[2], groups[1]),wh01] = c(as.numeric(PossMat01_1[samp2,]), 1)
  }
  
  activeZ = ActivePred(tempZeta, k=k)
  
  wv3 = which(tempZeta[,h] == 1)
  if (length(wv3) == 0) {
    tempBeta[[h]] = c(0)
    f_jhi_nc[,h] = rep(0, n)
  } else if (length(wv3) == 1) {
    
    if (h > 1 & sum(tempZeta[wv3,1:(h-1)]) > 0) {
      tempBeta[[h]] = rep(0, ns+1)
      f_jhi_nc[,h] = rep(0,n)
    } else {
      tempXstar = Xstar[,wv3,]
      SigmaB = diag(sigmaP*sigB, dim(tempXstar)[2])
      muB = rep(0, dim(tempXstar)[2])
      muBeta = solve((t(tempXstar) %*% tempXstar)/
                       sigmaP + solve(SigmaB)) %*%
        ((t(tempXstar) %*% tempY)/
           sigmaP + solve(SigmaB) %*% muB)
      covBeta = solve((t(tempXstar) %*% tempXstar)/
                        sigmaP + solve(SigmaB))
      
      tempBeta[[h]] = mvtnorm::rmvnorm(1, muBeta, as.matrix(covBeta))
      f_jhi_nc[,h] = tempXstar %*% as.vector(tempBeta[[h]])
    }
  } else {
    
    tempXstar = matrix(NA, n, (ns+1)^length(wv3))
    nsMat = expand.grid(rep(list(1 : (ns + 1)), length(wv3)))
    activeMat = expand.grid(rep(list(c(0,rep(1,ns))), length(wv3)))
    activeInd = "0"
    for (ii in 1 : nrow(nsMat)) {
      if (sum(activeMat[ii,]) > 0) {
        activeInd = c(activeInd, paste(sort(wv3[which(activeMat[ii,]==1)]), collapse="-"))
      }
      tempXstar[,ii] = Xstar[,wv3[1], nsMat[ii,1]]
      for (jj in 2 : length(wv3)) {
        tempXstar[,ii] = tempXstar[,ii]*Xstar[,wv3[jj],nsMat[ii,jj]]
      }
    }
    
    tempXstar = cbind(rep(1, n), scale(cbind(tempXstar[,-1])))
    
    w = 1:dim(tempXstar)[2]
    if (h > 1) {
      w = which(!activeInd %in% unlist(activeZ[1:(h-1)]))
    }
    
    
    if (length(w) == 0) {
      tempBeta[[h]] = rep(0, dim(tempXstar)[2])
      f_jhi_nc[,h] = tempXstar %*% as.vector(tempBeta[[h]])
    } else if (length(w) == 1) {
      tempXstar2 = tempXstar[,w]
      
      SigmaB = sigB
      muB = rep(0, 1)
      muBeta = solve((t(tempXstar2) %*% tempXstar2)/
                       sigmaP + solve(SigmaB)) %*%
        ((t(tempXstar2) %*% tempY)/
           sigmaP + solve(SigmaB) %*% muB)
      covBeta = solve((t(tempXstar2) %*% tempXstar2)/
                        sigmaP + solve(SigmaB))
      
      temp = rep(0, dim(tempXstar)[2])
      temp[w] = mvtnorm::rmvnorm(1, muBeta, as.matrix(covBeta))
      tempBeta[[h]] = temp
      f_jhi_nc[,h] = tempXstar %*% as.vector(tempBeta[[h]])
    } else {
      tempXstar2 = tempXstar[,w]
      
      SigmaB = diag(sigmaP*sigB, dim(tempXstar2)[2])
      muB = rep(0, dim(tempXstar2)[2])
      muBeta = solve((t(tempXstar2) %*% tempXstar2)/
                       sigmaP + solve(SigmaB)) %*%
        ((t(tempXstar2) %*% tempY)/
           sigmaP + solve(SigmaB) %*% muB)
      covBeta = solve((t(tempXstar2) %*% tempXstar2)/
                        sigmaP + solve(SigmaB))
      
      temp = rep(0, dim(tempXstar)[2])
      temp[w] = mvtnorm::rmvnorm(1, muBeta, as.matrix(covBeta))
      tempBeta[[h]] = temp
      f_jhi_nc[,h] = tempXstar %*% as.vector(tempBeta[[h]])
    }
    
  }
  
  if (!samp %in% c("00", "10", "01", "11")) {
    
    wv3_01 = which(tempZeta[,wh01] == 1)
    if (length(wv3_01) == 0) {
      tempBeta[[wh01]] = c(0)
      f_jhi_nc[,wh01] = rep(0, n)
    } else if (length(wv3_01) == 1) {
      
      if (wh01 > 1 & sum(tempZeta[wv3_01,1:(wh01-1)]) > 0) {
        tempBeta[[wh01]] = rep(0, ns+1)
        f_jhi_nc[,wh01] = rep(0, n)
      } else {
        tempXstar = Xstar[,wv3_01,]
        SigmaB = diag(sigmaP*sigB, dim(tempXstar)[2])
        muB = rep(0, dim(tempXstar)[2])
        muBeta = solve((t(tempXstar) %*% tempXstar)/
                         sigmaP + solve(SigmaB)) %*%
          ((t(tempXstar) %*% tempY)/
             sigmaP + solve(SigmaB) %*% muB)
        covBeta = solve((t(tempXstar) %*% tempXstar)/
                          sigmaP + solve(SigmaB))
        
        tempBeta[[wh01]] = mvtnorm::rmvnorm(1, muBeta, as.matrix(covBeta))
        f_jhi_nc[,wh01] = tempXstar %*% as.vector(tempBeta[[wh01]])
      }
    } else {
      
      tempXstar = matrix(NA, n, (ns+1)^length(wv3_01))
      nsMat = expand.grid(rep(list(1 : (ns + 1)), length(wv3_01)))
      activeMat = expand.grid(rep(list(c(0,rep(1,ns))), length(wv3_01)))
      activeInd = "0"
      for (ii in 1 : nrow(nsMat)) {
        if (sum(activeMat[ii,]) > 0) {
          activeInd = c(activeInd, paste(sort(wv3_01[which(activeMat[ii,]==1)]), collapse="-"))
        }
        tempXstar[,ii] = Xstar[,wv3_01[1], nsMat[ii,1]]
        for (jj in 2 : length(wv3_01)) {
          tempXstar[,ii] = tempXstar[,ii]*Xstar[,wv3_01[jj],nsMat[ii,jj]]
        }
      }
      
      tempXstar = cbind(rep(1, n), scale(cbind(tempXstar[,-1])))
      
      w = 1:dim(tempXstar)[2]
      if (h > 1) {
        w = which(!activeInd %in% unlist(activeZ[1:(h-1)]))
      }
      
      
      if (length(w) == 0) {
        tempBeta[[wh01]] = rep(0, dim(tempXstar)[2])
        f_jhi_nc[,wh01] = tempXstar %*% as.vector(tempBeta[[wh01]])
      } else if (length(w) == 1) {
        tempXstar2 = tempXstar[,w]
        
        SigmaB = sigB
        muB = rep(0, 1)
        muBeta = solve((t(tempXstar2) %*% tempXstar2)/
                         sigmaP + solve(SigmaB)) %*%
          ((t(tempXstar2) %*% tempY)/
             sigmaP + solve(SigmaB) %*% muB)
        covBeta = solve((t(tempXstar2) %*% tempXstar2)/
                          sigmaP + solve(SigmaB))
        
        temp = rep(0, dim(tempXstar)[2])
        temp[w] = mvtnorm::rmvnorm(1, muBeta, as.matrix(covBeta))
        tempBeta[[wh01]] = temp
        f_jhi_nc[,wh01] = tempXstar %*% as.vector(tempBeta[[wh01]])
      } else {
        tempXstar2 = tempXstar[,w]
        
        SigmaB = diag(sigmaP*sigB, dim(tempXstar2)[2])
        muB = rep(0, dim(tempXstar2)[2])
        muBeta = solve((t(tempXstar2) %*% tempXstar2)/
                         sigmaP + solve(SigmaB)) %*%
          ((t(tempXstar2) %*% tempY)/
             sigmaP + solve(SigmaB) %*% muB)
        covBeta = solve((t(tempXstar2) %*% tempXstar2)/
                          sigmaP + solve(SigmaB))
        
        temp = rep(0, dim(tempXstar)[2])
        temp[w] = mvtnorm::rmvnorm(1, muBeta, as.matrix(covBeta))
        tempBeta[[wh01]] = temp
        f_jhi_nc[,wh01] = tempXstar %*% as.vector(tempBeta[[wh01]])
      }
      
    }
  }
  
  
  
  
  return(list(f_jhi_nc = f_jhi_nc,
              beta = tempBeta,
              zeta = tempZeta))
  
}



MCMCmixture = function(Y, X, C, Xstar, nChains = 2, nIter = 30000, nBurn = 10000, thin = 20,
                       c = 0.001, d = 0.001, sigB, muB, SigmaC, muC,
                       k = 10, ns = 4, alph, gamm, probSamp1=1) {
  
  designC = cbind(rep(1,n), C)
  
  n = dim(X)[1]
  p = dim(X)[2]
  pc = dim(designC)[2] - 1
  
  zetaPost = array(NA, dim=c(nChains, nIter, p, k))
  alphaPost = array(NA, dim=c(nChains, nIter))
  gammaPost = array(NA, dim=c(nChains, nIter))
  tauPost = array(NA, dim=c(nChains, nIter, k))
  sigmaPost = array(NA, dim=c(nChains, nIter))
  betaPost = array(NA, dim=c(nChains, nIter, p, k, ns+1))
  betaCPost = array(NA, dim=c(nChains, nIter, pc+1))
  
  ## set initial values
  
  for (nc in 1 : nChains) {
    zetaPost[nc,1,,] = 0
  }
  
  alphaPost[,1] = alph
  gammaPost[,1] = gamm
  
  tauPost[,1,] = runif(nChains*k)
  
  sigmaPost[,1] = rgamma(nChains, 1, 1)
  betaPost[,1,,,] = 0
  betaPost[,1,,,1] = 1
  betaCPost[,1,] = rnorm(nChains*(pc+1), sd=0.1)
  
  ## array that calculates value of functions for each subject for any cluster or covariate
  f_jhi = array(NA, dim=c(nChains,n,k))
  for (nc in 1 : nChains) {
    for (h in 1 : k) {
      f_jhi[nc,,h] = 0
    }
  }
  
  betaList = list()
  betaList[[1]] = list()
  for (nc in 1 : nChains) {
    betaList[[1]][[nc]] = list()
    for (h in 1 : k) {
      betaList[[1]][[nc]][[h]] = c(0)
    }
  }
  
  ## begin MCMC
  for (ni in 2 : nIter) {
    betaList[[ni]] = list()
    for (nc in 1 : nChains) {
      
      betaList[[ni]][[nc]] = list()
      for (h in 1 : k) {
        betaList[[ni]][[nc]][[h]] = c(0)
      }
      
      if(ni %% 100 == 0 & nc == 1) print(ni)
      
      alphaPost[nc,ni] = alph
      gammaPost[nc,ni] = gamm
      
      ## sample sigma squared
      tempSum0 = 0
      tempSum2 = 0
      for (h in 1 : k) {
        tempSum0 = tempSum0 + sum(unlist(betaList[[ni-1]][[nc]][[h]])[-1] != 0)
        tempSum2 = tempSum2 + sum(unlist(betaList[[ni-1]][[nc]][[h]])[-1]^2)
      }
      SampleDiffSq = (Y - apply(f_jhi[nc,,], 1, sum) - 
                        (designC %*% betaCPost[nc,ni-1,]))^2
      sigmaPost[nc,ni] = 1/rgamma(1, c + n/2 + tempSum0/2, d + sum(SampleDiffSq)/2 + tempSum2/(2*sigB))
      
      ## update tau_h
      for (h in 1 : k) {
        tauPost[nc,ni,h] = rbeta(1, (alphaPost[nc,ni]*gammaPost[nc,ni]/k) + sum(zetaPost[nc,ni-1,,h]), 
                                 gammaPost[nc,ni] + sum(1 - zetaPost[nc,ni-1,,h]))
      }
      
      tempZeta = zetaPost[nc,ni-1,,]
      tempBeta = betaList[[ni-1]][[nc]]
      for (h in 1 : k) {
        
        OneOrTwo = sample(1:2, 1, p=c(probSamp1, 1-probSamp1), replace=FALSE)
        
        if (OneOrTwo == 1) {
          UpdateBeta = UpdateBetaOne(tempZeta=tempZeta, f_jhi_nc=f_jhi[nc,,],
                                     betaC = betaCPost[nc,ni-1,], sigmaP=sigmaPost[nc,ni],
                                     tau=tauPost[nc,ni,], k=k, sigB=sigB, Xstar=Xstar,
                                     tempBeta = tempBeta, h=h, designC=designC, ns=ns)
        } else {
          UpdateBeta = UpdateBetaTwo(tempZeta=tempZeta, f_jhi_nc=f_jhi[nc,,],
                                     betaC = betaCPost[nc,ni-1,], sigmaP=sigmaPost[nc,ni],
                                     tau=tauPost[nc,ni,], k=k, sigB=sigB, Xstar=Xstar,
                                     tempBeta = tempBeta, h=h, designC=designC, ns=ns)
        }

        tempZeta = UpdateBeta$zeta
        tempBeta = UpdateBeta$beta
        f_jhi[nc,,] = UpdateBeta$f_jhi_nc
      }
      
      betaList[[ni]][[nc]] = tempBeta
      
      ## Remove redundant columns and combine coefficients
      for (h in 1 : k) {
        for (h2 in (1 : k)[-h]) {
          
          wh = which(tempZeta[,h] == 1)
          wh2 = which(tempZeta[,h2] == 1)
          if (all(wh2 %in% wh) == TRUE & length(wh2) > 0) {
            
            nsMat = expand.grid(rep(list(1 : (ns + 1)), length(wh)))
            nsMat2 = data.frame(matrix(1, (ns+1)^length(wh2), length(wh)))
            names(nsMat) = paste("X", wh, sep='')
            names(nsMat2) = paste("X", wh, sep='')
            
            nsMat2[,which(wh %in% wh2)] = expand.grid(rep(list(1 : (ns + 1)), length(wh2)))
            
            Names = do.call(paste, c(nsMat, sep="-"))
            Names2 = do.call(paste, c(nsMat2, sep="-"))
            
            wCombine = match(Names2, Names)
            
            betaList[[ni]][[nc]][[h]][wCombine] = betaList[[ni]][[nc]][[h]][wCombine] + 
              betaList[[ni]][[nc]][[h2]]
            
            
            if(length(wh) == 1) {
              tempXstar = Xstar[,wh,]
            } else {
              tempXstar = matrix(NA, n, (ns+1)^length(wh))
              nsMat = expand.grid(rep(list(1 : (ns + 1)), length(wh)))
              for (ii in 1 : nrow(nsMat)) {
                tempXstar[,ii] = Xstar[,wh[1], nsMat[ii,1]]
                for (jj in 2 : length(wh)) {
                  tempXstar[,ii] = tempXstar[,ii]*Xstar[,wh[jj],nsMat[ii,jj]]
                }
              }
            }
            
            f_jhi[nc,,h] = tempXstar %*%  as.vector(betaList[[ni]][[nc]][[h]])
            
            tempZeta[,h2] = rep(0,p)
            f_jhi[nc,,h2] = rep(0, n)
            
            betaList[[ni]][[nc]][[h2]] = c(0)
          }
        }
      }
      
      
      zetaPost[nc,ni,,] = tempZeta
      
      ## Update betaC
      tempY = Y - apply(f_jhi[nc,,], 1, sum)
      
      muBetaC = solve((t(designC) %*% designC)/sigmaPost[nc,ni-1] + solve(SigmaC)) %*% 
        ((t(designC) %*% tempY)/sigmaPost[nc,ni-1] + solve(SigmaC) %*% muC)
      covBetaC = solve((t(designC) %*% designC)/sigmaPost[nc,ni-1] + solve(SigmaC))
      
      betaCPost[nc,ni,] = mvtnorm::rmvnorm(1, mean=muBetaC, sigma = covBetaC)
    }
  }
  
  keep = nBurn + (1:floor((nIter - nBurn) / thin))*thin
  
  posterior = list(alpha = alphaPost[,keep],
                   gamma = gammaPost[,keep],
                   sigma = sigmaPost[,keep],
                   tau = tauPost[,keep,],
                   zeta = zetaPost[,keep,,],
                   beta = betaList[keep],
                   betaC = betaCPost[,keep,])
  
  return(posterior)
}


MCMCmixtureEB = function(Y, X, C, Xstar, nChains = 2, nIter = 30000, nBurn = 10000, thin = 20,
                         c = 0.001, d = 0.001, sigBstart, muB, SigmaC, muC,
                         k = 10, ns = 4, alph, gamm, probSamp1=1, SigMin = 0.1) {
  
  designC = cbind(rep(1,n), C)
  
  n = dim(X)[1]
  p = dim(X)[2]
  pc = dim(designC)[2] - 1
  
  zetaPost = array(NA, dim=c(nChains, nIter, p, k))
  alphaPost = array(NA, dim=c(nChains, nIter))
  gammaPost = array(NA, dim=c(nChains, nIter))
  tauPost = array(NA, dim=c(nChains, nIter, k))
  sigmaPost = array(NA, dim=c(nChains, nIter))
  betaPost = array(NA, dim=c(nChains, nIter, p, k, ns+1))
  betaCPost = array(NA, dim=c(nChains, nIter, pc+1))
  sigBPost = array(NA, dim=c(nChains, nIter))
  
  ## set initial values
  
  for (nc in 1 : nChains) {
    zetaPost[nc,1,,] = 0
  }
  
  alphaPost[,1] = alph
  gammaPost[,1] = gamm
  
  tauPost[,1,] = runif(nChains*k)
  
  sigmaPost[,1] = rgamma(nChains, 1, 1)
  betaPost[,1,,,] = 0
  betaPost[,1,,,1] = 1
  betaCPost[,1,] = rnorm(nChains*(pc+1), sd=0.1)
  sigBPost[,1] = sigBstart
  sigB = sigBPost[1,1]
  
  ## array that calculates value of functions for each subject for any cluster or covariate
  f_jhi = array(NA, dim=c(nChains,n,k))
  for (nc in 1 : nChains) {
    for (h in 1 : k) {
      f_jhi[nc,,h] = 0
    }
  }
  
  betaList = list()
  betaList[[1]] = list()
  for (nc in 1 : nChains) {
    betaList[[1]][[nc]] = list()
    for (h in 1 : k) {
      betaList[[1]][[nc]][[h]] = c(0)
    }
  }
  
  ## begin MCMC
  for (ni in 2 : nIter) {
    betaList[[ni]] = list()
    for (nc in 1 : nChains) {
      
      betaList[[ni]][[nc]] = list()
      for (h in 1 : k) {
        betaList[[ni]][[nc]][[h]] = c(0)
      }
      
      if(ni %% 100 == 0 & nc == 1) print(ni)
      
      alphaPost[nc,ni] = alph
      gammaPost[nc,ni] = gamm
      
      ## sample sigma squared
      tempSum0 = 0
      tempSum2 = 0
      for (h in 1 : k) {
        tempSum0 = tempSum0 + sum(unlist(betaList[[ni-1]][[nc]][[h]])[-1] != 0)
        tempSum2 = tempSum2 + sum(unlist(betaList[[ni-1]][[nc]][[h]])[-1]^2)
      }
      SampleDiffSq = (Y - apply(f_jhi[nc,,], 1, sum) - 
                        (designC %*% betaCPost[nc,ni-1,]))^2
      sigmaPost[nc,ni] = 1/rgamma(1, c + n/2 + tempSum0/2, d + sum(SampleDiffSq)/2 + tempSum2/(2*sigB))
      
      ## update tau_h
      for (h in 1 : k) {
        tauPost[nc,ni,h] = rbeta(1, (alphaPost[nc,ni]*gammaPost[nc,ni]/k) + sum(zetaPost[nc,ni-1,,h]), 
                                 gammaPost[nc,ni] + sum(1 - zetaPost[nc,ni-1,,h]))
      }
      
      tempZeta = zetaPost[nc,ni-1,,]
      tempBeta = betaList[[ni-1]][[nc]]
      for (h in 1 : k) {
        
        OneOrTwo = sample(1:2, 1, p=c(probSamp1, 1-probSamp1), replace=FALSE)
        
        if (OneOrTwo == 1) {
          UpdateBeta = UpdateBetaOne(tempZeta=tempZeta, f_jhi_nc=f_jhi[nc,,],
                                     betaC = betaCPost[nc,ni-1,], sigmaP=sigmaPost[nc,ni],
                                     tau=tauPost[nc,ni,], k=k, sigB=sigB, Xstar=Xstar,
                                     tempBeta = tempBeta, h=h, designC=designC, ns=ns)
        } else {
          UpdateBeta = UpdateBetaTwo(tempZeta=tempZeta, f_jhi_nc=f_jhi[nc,,],
                                     betaC = betaCPost[nc,ni-1,], sigmaP=sigmaPost[nc,ni],
                                     tau=tauPost[nc,ni,], k=k, sigB=sigB, Xstar=Xstar,
                                     tempBeta = tempBeta, h=h, designC=designC, ns=ns)
        }
        
        tempZeta = UpdateBeta$zeta
        tempBeta = UpdateBeta$beta
        f_jhi[nc,,] = UpdateBeta$f_jhi_nc
      }
      
      betaList[[ni]][[nc]] = tempBeta
      
      ## Remove redundant columns and combine coefficients
      for (h in 1 : k) {
        for (h2 in (1 : k)[-h]) {
          
          wh = which(tempZeta[,h] == 1)
          wh2 = which(tempZeta[,h2] == 1)
          if (all(wh2 %in% wh) == TRUE & length(wh2) > 0) {
            
            nsMat = expand.grid(rep(list(1 : (ns + 1)), length(wh)))
            nsMat2 = data.frame(matrix(1, (ns+1)^length(wh2), length(wh)))
            names(nsMat) = paste("X", wh, sep='')
            names(nsMat2) = paste("X", wh, sep='')
            
            nsMat2[,which(wh %in% wh2)] = expand.grid(rep(list(1 : (ns + 1)), length(wh2)))
            
            Names = do.call(paste, c(nsMat, sep="-"))
            Names2 = do.call(paste, c(nsMat2, sep="-"))
            
            wCombine = match(Names2, Names)
            
            betaList[[ni]][[nc]][[h]][wCombine] = betaList[[ni]][[nc]][[h]][wCombine] + 
              betaList[[ni]][[nc]][[h2]]
            
            
            if(length(wh) == 1) {
              tempXstar = Xstar[,wh,]
            } else {
              tempXstar = matrix(NA, n, (ns+1)^length(wh))
              nsMat = expand.grid(rep(list(1 : (ns + 1)), length(wh)))
              for (ii in 1 : nrow(nsMat)) {
                tempXstar[,ii] = Xstar[,wh[1], nsMat[ii,1]]
                for (jj in 2 : length(wh)) {
                  tempXstar[,ii] = tempXstar[,ii]*Xstar[,wh[jj],nsMat[ii,jj]]
                }
              }
            }
            
            f_jhi[nc,,h] = tempXstar %*%  as.vector(betaList[[ni]][[nc]][[h]])
            
            tempZeta[,h2] = rep(0,p)
            f_jhi[nc,,h2] = rep(0, n)
            
            betaList[[ni]][[nc]][[h2]] = c(0)
          }
        }
      }
      
      
      zetaPost[nc,ni,,] = tempZeta

      
      if (ni %% 50 == 0) {
        sumsB = 0
        sumsZ = 0
        for (ni2 in (ni - 49) : ni) {
          for (h in 1 : k) {
            sumsB = sumsB + sum(unlist(betaList[[ni2]][[nc]][[h]])[-1]^2)/sigmaPost[nc,ni2]
            sumsZ = sumsZ + length(which(unlist(betaList[[ni2]][[nc]][[h]])[-1] != 0))
          }
        }
        
        if (sumsZ == 0) {
          sigBPost[nc,ni] = sigBPost[nc,ni-1]
        } else {
          sigBPost[nc,ni] = max(sumsB / sumsZ, SigMin)
          sigB = sigBPost[nc,ni]
        }
      } else {
        sigBPost[nc,ni] = sigBPost[nc,ni-1]
        sigB = sigBPost[nc,ni]
      }
      
      ## Update betaC
      tempY = Y - apply(f_jhi[nc,,], 1, sum)
      
      muBetaC = solve((t(designC) %*% designC)/sigmaPost[nc,ni-1] + solve(SigmaC)) %*% 
        ((t(designC) %*% tempY)/sigmaPost[nc,ni-1] + solve(SigmaC) %*% muC)
      covBetaC = solve((t(designC) %*% designC)/sigmaPost[nc,ni-1] + solve(SigmaC))
      
      betaCPost[nc,ni,] = mvtnorm::rmvnorm(1, mean=muBetaC, sigma = covBetaC)
    }
  }
  
  keep = nBurn + (1:floor((nIter - nBurn) / thin))*thin
  
  return(mean(sigBPost[,nIter]))
}

MCMCmixtureMinSig = function(Y, X, C, Xstar, nPerms = 10, nIter = 500,
                             c = 0.001, d = 0.001, sigBstart = 1, muB, SigmaC, muC,
                             k = 10, ns = 2, threshold=0.25) {
  
  Ysave = Y
  sigSave = rep(NA, nPerms)
  
  for (ss in 1 : nPerms) {
    
    nChains = 1
    
    designC = cbind(rep(1,n), C)
    
    n = dim(X)[1]
    p = dim(X)[2]
    pc = dim(designC)[2] - 1
    
    ## Permute the data
    perms = sample(1:n, n, replace=FALSE)
    Y = Ysave[perms]
    
    zetaPost = array(NA, dim=c(nChains, nIter, p, k))
    tauPost = array(NA, dim=c(nChains, nIter, k))
    sigmaPost = array(NA, dim=c(nChains, nIter))
    betaPost = array(NA, dim=c(nChains, nIter, p, k, ns+1))
    betaCPost = array(NA, dim=c(nChains, nIter, pc+1))
    
    sigB = sigBstart
    
    ## set initial values
    
    for (nc in 1 : nChains) {
      zetaPost[nc,1,,] = 0
    }
    
    tauPost[,1,] = 0.5
    
    sigmaPost[,1] = rgamma(nChains, 1, 1)
    betaPost[,1,,,] = 0
    betaPost[,1,,,1] = 1
    betaCPost[,1,] = rnorm(nChains*(pc+1), sd=0.1)
    
    ## array that calculates value of functions for each subject for any cluster or covariate
    f_jhi = array(NA, dim=c(nChains,n,k))
    for (nc in 1 : nChains) {
      for (h in 1 : k) {
        f_jhi[nc,,h] = 0
      }
    }
    
    betaList = list()
    betaList[[1]] = list()
    for (nc in 1 : nChains) {
      betaList[[1]][[nc]] = list()
      for (h in 1 : k) {
        betaList[[1]][[nc]][[h]] = c(0)
      }
    }
    
    ## begin MCMC
    for (ni in 2 : nIter) {
      betaList[[ni]] = list()
      for (nc in 1 : nChains) {
        
        betaList[[ni]][[nc]] = list()
        for (h in 1 : k) {
          betaList[[ni]][[nc]][[h]] = c(0)
        }
        
        ## sample sigma squared
        tempSum0 = 0
        tempSum2 = 0
        for (h in 1 : k) {
          tempSum0 = tempSum0 + sum(unlist(betaList[[ni-1]][[nc]][[h]])[-1] != 0)
          tempSum2 = tempSum2 + sum(unlist(betaList[[ni-1]][[nc]][[h]])[-1]^2)
        }
        SampleDiffSq = (Y - apply(f_jhi[nc,,], 1, sum) - 
                          (designC %*% betaCPost[nc,ni-1,]))^2
        sigmaPost[nc,ni] = 1/rgamma(1, c + n/2 + tempSum0/2, d + sum(SampleDiffSq)/2 + tempSum2/(2*sigB))
        
        ## update tau_h
        for (h in 1 : k) {
          tauPost[nc,ni,h] = 0.5
        }
        
        tempZeta = zetaPost[nc,ni-1,,]
        
        if (ni < nIter) {
          tempZeta = matrix(0, p, k)
          zetaPost[nc,ni,,] = tempZeta
        } else {
          numZero = length(which(apply(tempZeta[,-h], 2, sum) == 0))
          h = 1
          groups = 1
          tempY = Y - apply(f_jhi[nc,,-h], 1, sum) -
            (designC %*% betaCPost[nc,ni-1,])
          
          tempZeta10 = tempZeta00 = tempZeta01 = tempZeta
          tempZeta10[groups,h] = 1
          tempZeta00[groups,h] = 0
          
          wv = which(tempZeta00[,h] == 1)
          
          gridSig = .75^(1:30)
          
          distinguish = TRUE
          counter = 1
          while (distinguish == TRUE) {
            
            sigB = gridSig[counter]
            
            p00 = p10 = -Inf
            
            ## First calculate probability for reduced model
            
            if (length(wv) == 0) {
              p00 = log(1 - tauPost[nc,ni,h])
            } else if (length(wv) == 1) {
              SigmaB = diag(sigmaPost[nc,ni]*sigB, ns+1)
              muB = rep(0, ns+1)
              muBeta = solve((t(Xstar[,wv,]) %*% Xstar[,wv,])/
                               sigmaPost[nc,ni] + solve(SigmaB)) %*%
                ((t(Xstar[,wv,]) %*% tempY)/
                   sigmaPost[nc,ni] + solve(SigmaB) %*% muB)
              covBeta = solve((t(Xstar[,wv,]) %*% Xstar[,wv,])/
                                sigmaPost[nc,ni] + solve(SigmaB))
              
              p00 = log(tauPost[nc,ni,h] * (1 - tauPost[nc,ni,h])) + 
                mvtnorm::dmvnorm(rep(0, ns), mean=muB[-1], sigma=as.matrix(SigmaB[-1,-1]), log=TRUE) - 
                mvtnorm::dmvnorm(rep(0, ns), mean=muBeta[-1], sigma=as.matrix(covBeta[-1,-1]), log=TRUE)
            } else {
              tempXstar = matrix(NA, n, (ns+1)^length(wv))
              nsMat = expand.grid(rep(list(1 : (ns + 1)), length(wv)))
              for (ii in 1 : nrow(nsMat)) {
                tempXstar[,ii] = Xstar[,wv[1], nsMat[ii,1]]
                for (jj in 2 : length(wv)) {
                  tempXstar[,ii] = tempXstar[,ii]*Xstar[,wv[jj],nsMat[ii,jj]]
                }
              }
              
              tempXstar = cbind(rep(1, n), scale(tempXstar[,-1]))
              
              SigmaB = diag(sigmaPost[nc,ni]*sigB, (ns+1)^length(wv))
              muB = rep(0, (ns+1)^length(wv))
              muBeta = solve((t(tempXstar) %*% tempXstar)/
                               sigmaPost[nc,ni] + solve(SigmaB)) %*%
                ((t(tempXstar) %*% tempY)/
                   sigmaPost[nc,ni] + solve(SigmaB) %*% muB)
              covBeta = solve((t(tempXstar) %*% tempXstar)/
                                sigmaPost[nc,ni] + solve(SigmaB))
              
              p00 = log((tauPost[nc,ni,h]^length(wv)) * (1 - tauPost[nc,ni,h])) + 
                mvtnorm::dmvnorm(rep(0, length(muB)-1), mean=muB[-1], sigma=as.matrix(SigmaB[-1,-1]), log=TRUE) - 
                mvtnorm::dmvnorm(rep(0, length(muB)-1), mean=muBeta[-1], sigma=as.matrix(covBeta[-1,-1]), log=TRUE)
            }
            
            ## Now calculate probability for full interactions
            if ((numZero > 0 & sum(tempZeta10[,h]) > 1) |
                sum(tempZeta10[,h]) == 1) {
              if (length(wv) == 0) {
                SigmaB = diag(sigmaPost[nc,ni]*sigB, ns+1)
                muB = rep(0, ns+1)
                muBeta = solve((t(Xstar[,groups,]) %*% Xstar[,groups,])/
                                 sigmaPost[nc,ni] + solve(SigmaB)) %*%
                  ((t(Xstar[,groups,]) %*% tempY)/
                     sigmaPost[nc,ni] + solve(SigmaB) %*% muB)
                covBeta = solve((t(Xstar[,groups,]) %*% Xstar[,groups,])/
                                  sigmaPost[nc,ni] + solve(SigmaB))
                
                p10 = log(tauPost[nc,ni,h]) + 
                  mvtnorm::dmvnorm(rep(0, ns), mean=muB[-1], sigma=as.matrix(SigmaB[-1,-1]), log=TRUE) - 
                  mvtnorm::dmvnorm(rep(0, ns), mean=muBeta[-1], sigma=as.matrix(covBeta[-1,-1]), log=TRUE)
              } else {
                wv2 = c(wv, groups)
                tempXstar = matrix(NA, n, (ns+1)^length(wv2))
                nsMat = expand.grid(rep(list(1 : (ns + 1)), length(wv2)))
                for (ii in 1 : nrow(nsMat)) {
                  tempXstar[,ii] = Xstar[,wv2[1], nsMat[ii,1]]
                  for (jj in 2 : length(wv2)) {
                    tempXstar[,ii] = tempXstar[,ii]*Xstar[,wv2[jj],nsMat[ii,jj]]
                  }
                }
                
                tempXstar = cbind(rep(1, n), scale(tempXstar[,-1]))
                
                SigmaB = diag(sigmaPost[nc,ni]*sigB, (ns+1)^length(wv2))
                muB = rep(0, (ns+1)^length(wv2))
                muBeta = solve((t(tempXstar) %*% tempXstar)/
                                 sigmaPost[nc,ni] + solve(SigmaB)) %*%
                  ((t(tempXstar) %*% tempY)/
                     sigmaPost[nc,ni] + solve(SigmaB) %*% muB)
                covBeta = solve((t(tempXstar) %*% tempXstar)/
                                  sigmaPost[nc,ni] + solve(SigmaB))
                
                p10 = log((tauPost[nc,ni,h]^length(wv2))) + 
                  mvtnorm::dmvnorm(rep(0, length(muB)-1), mean=muB[-1], sigma=as.matrix(SigmaB[-1,-1]), log=TRUE) - 
                  mvtnorm::dmvnorm(rep(0, length(muB)-1), mean=muBeta[-1], sigma=as.matrix(covBeta[-1,-1]), log=TRUE)
              }
            }
            
            maxlog = max(c(p00, p10))
            updated_p = exp(c(p00, p10) - maxlog)
            p0 = updated_p[1] / sum(updated_p)
            p1 = 1 - p0
            
            if (p1 > threshold) {
              distinguish = FALSE
            }
            counter = counter + 1
          }
        }
        
        ## Update betaC
        tempY = Y - apply(f_jhi[nc,,], 1, sum)
        
        muBetaC = solve((t(designC) %*% designC)/sigmaPost[nc,ni-1] + solve(SigmaC)) %*% 
          ((t(designC) %*% tempY)/sigmaPost[nc,ni-1] + solve(SigmaC) %*% muC)
        covBetaC = solve((t(designC) %*% designC)/sigmaPost[nc,ni-1] + solve(SigmaC))
        
        betaCPost[nc,ni,] = mvtnorm::rmvnorm(1, mean=muBetaC, sigma = covBetaC)
      }
    }
    
    sigSave[ss] = sigB
  }
  
  return(median(sigSave))
}

