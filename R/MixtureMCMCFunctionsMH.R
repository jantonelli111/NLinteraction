ActivePred_MH = function(zeta, k) {
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

numextract_MH <- function(string){
  unlist(regmatches(string,gregexpr("[[:digit:]]+\\.*[[:digit:]]*",string)))
}

ChooseModels2_MH = function(zeta, dim1, dim2, k=k, intMax) {
  tempZeta0 = tempZeta1 = tempZeta = tempZetaColOnly0 = tempZetaColOnly1 = zeta
  tempZeta0[dim1,dim2] = 0
  tempZeta1[dim1,dim2] = 1
  
  ## remove columns that are subsets of the one being considered
  wv = which(tempZeta1[,dim2] == 1)
  for (k1 in (1 : k)[-dim2]) {
    wv2 = which(tempZeta1[,k1] == 1)
    if (all(wv2 %in% wv)) {
      tempZetaColOnly0[,k1] = 0
      tempZetaColOnly1[,k1] = 0
    }
  }
  
  tempZetaColOnly0[dim1,dim2] = 0
  tempZetaColOnly1[dim1,dim2] = 1
  
  model0 = ActivePred_MH(tempZeta0, k=k)
  model1 = ActivePred_MH(tempZeta1, k=k)
  
  modelColOnly0 = ActivePred_MH(tempZetaColOnly0, k=k)
  modelColOnly1 = ActivePred_MH(tempZetaColOnly1, k=k)
  
  modelSet0 = as.character(unique(unlist(model0)))
  modelSet1 = as.character(unique(unlist(model1)))
  
  modelSetColOnly0 = as.character(unique(unlist(modelColOnly0)))
  modelSetColOnly1 = as.character(unique(unlist(modelColOnly1)))
  
  ## Now add back in the elements of modelSet0 that aren't in ColOnly1
  tempSet = modelSet0[which(!modelSet0 %in% modelSetColOnly1)]
  modelSetColOnly0 = c(modelSetColOnly0, tempSet)
  modelSetColOnly1 = c(modelSetColOnly1, tempSet)
  stlVecColOnly1 = unlist(lapply(lapply(modelSetColOnly1, numextract_MH), length))
  
  
  newSet = modelSetColOnly1[which(!modelSetColOnly1 %in% modelSetColOnly0)]
  newSetReduced = newSet[-which.max(stringr::str_length(newSet))]
  
  stlVec = unlist(lapply(lapply(newSet, numextract_MH), length))
  
  wKeep = which(stlVec <= intMax)
  
  stlVec = stlVec[wKeep]
  newSet = newSet[wKeep]

  lSet = length(newSet)
  
  considerModels = list()
  considerModels[[1]] = sort(modelSetColOnly0)
  
  ## Keep track of the column that would need to be updated to get back to original matrix
  whichColUpdate = c()
  
  modToZeta = modelToZeta_MH(considerModels[[1]], k=k)
  isValid = rep(NA, ncol(modToZeta))
  for (jj in 1 : ncol(modToZeta)) {
    isValid[jj] = 1*(all(modToZeta[,jj] == tempZeta0[,dim2]))
  }
  whichColUpdate = c(whichColUpdate, which(isValid == 1)[1])
  
  ## only add in the model for tempZeta1 if it is different and less than intMax
  if (length(modelSetColOnly0) != length(modelSetColOnly1) &
      all(stlVecColOnly1 <= intMax)) {
    considerModels[[2]] = sort(modelSetColOnly1)
    
    modToZeta = modelToZeta_MH(considerModels[[2]], k=k)
    isValid = rep(NA, ncol(modToZeta))
    for (jj in 1 : ncol(modToZeta)) {
      isValid[jj] = 1*(all(modToZeta[,jj] == tempZeta1[,dim2]))
    }
    whichColUpdate = c(whichColUpdate, which(isValid == 1)[1])
  }
  
  numModelsTemp = length(considerModels)
  
  if (lSet > 0) {

    subsetTerms = list()
    for (j in 1 : lSet) {
      tempElements = as.numeric(numextract_MH(newSet[j]))
      subsetTerms[[j]] = which(sapply(lapply(lapply(newSet, numextract_MH), as.numeric), 
                                      function(x) return(all(x %in% tempElements))) &
                                 (1 : lSet) != j)
    }
    
    
    counter = numModelsTemp + 1
    for (jj in lSet : 1) {
      jjActive = which(!1:lSet %in% as.numeric(subsetTerms[[jj]]) & 1:lSet < jj)
      
      if (length(jjActive) == 0) {
        proposedSet = unique(c(modelSetColOnly0, newSet[[jj]]))
        
        addZeta = matrix(0, p,length(proposedSet))
        for (j in 1 : length(proposedSet)) {
          w = as.numeric(numextract_MH(proposedSet[j]))
          addZeta[w,j] = 1
        }
        zetaPlusZeta = cbind(tempZetaColOnly0, addZeta)
        
        newModel = as.character(unique(unlist(ActivePred_MH(zetaPlusZeta, k=dim(zetaPlusZeta)[2]))))
        
        keep = TRUE
        for (j in 1 : length(considerModels)) {
          if (length(newModel) == length(considerModels[[j]]) & 
              all(newModel %in% considerModels[[j]])) {
            keep = FALSE
          }
        }
        if (keep == TRUE) {
          considerModels[[counter]] = sort(newModel)
          counter = counter + 1
        }
      } else {
        tempGrid = expand.grid(rep(list(c(0,1)), length(jjActive)))
        
        for (ii in 1 : nrow(tempGrid)) {
          proposedSet = unique(c(modelSetColOnly0, newSet[[jj]], 
                                 newSet[jjActive[which(tempGrid[ii,] == 1)]]))
          
          addZeta = matrix(0, p,length(proposedSet))
          for (j in 1 : length(proposedSet)) {
            w = as.numeric(numextract_MH(proposedSet[j]))
            addZeta[w,j] = 1
          }
          zetaPlusZeta = cbind(tempZetaColOnly0, addZeta)
          
          newModel = as.character(unique(unlist(ActivePred_MH(zetaPlusZeta, k=dim(zetaPlusZeta)[2]))))
          
          keep = TRUE
          for (j in 1 : length(considerModels)) {
            if (length(newModel) == length(considerModels[[j]]) & 
                all(newModel %in% considerModels[[j]])) {
              keep = FALSE
            }
          }
          if (keep == TRUE) {
            considerModels[[counter]] = sort(newModel)
            counter = counter + 1
          }
        }
        
      }
      
    }
  }
  
  if (length(considerModels) > numModelsTemp) {
    for (j in (numModelsTemp + 1) : length(considerModels)) {
      modToZeta = modelToZeta_MH(considerModels[[j]], k=k)
      isValid = rep(NA, ncol(modToZeta))
      for (jj in 1 : ncol(modToZeta)) {
        isValid[jj] = 1*(all(modToZeta[,jj] == tempZeta0[,dim2]))
      }
      whichColUpdate = c(whichColUpdate, which(isValid == 1)[1])
    }
  }
  
  l = list(considerModels = considerModels,
           whichColUpdate = whichColUpdate)
  
  
  return(l)
}

ChooseModelsTwo2_MH = function(zeta, dim1, dim2, k=k, intMax) {
  tempZeta0 = tempZeta1 = tempZeta = tempZetaColOnly0 = tempZetaColOnly1 = zeta
  tempZeta0[dim1,dim2] = 0
  tempZeta1[dim1,dim2] = 1
  
  ## remove columns that are subsets of the one being considered
  wv = which(tempZeta1[,dim2] == 1)
  for (k1 in (1 : k)[-dim2]) {
    wv2 = which(tempZeta1[,k1] == 1)
    if (all(wv2 %in% wv)) {
      tempZetaColOnly0[,k1] = 0
      tempZetaColOnly1[,k1] = 0
    }
  }
  
  tempZetaColOnly0[dim1,dim2] = 0
  tempZetaColOnly1[dim1,dim2] = 1
  
  model0 = ActivePred_MH(tempZeta0, k=k)
  model1 = ActivePred_MH(tempZeta1, k=k)
  
  modelColOnly0 = ActivePred_MH(tempZetaColOnly0, k=k)
  modelColOnly1 = ActivePred_MH(tempZetaColOnly1, k=k)
  
  modelSet0 = as.character(unique(unlist(model0)))
  modelSet1 = as.character(unique(unlist(model1)))
  
  modelSetColOnly0 = as.character(unique(unlist(modelColOnly0)))
  modelSetColOnly1 = as.character(unique(unlist(modelColOnly1)))
  
  ## Now add back in the elements of modelSet0 that aren't in ColOnly1
  tempSet = modelSet0[which(!modelSet0 %in% modelSetColOnly1)]
  modelSetColOnly0 = c(modelSetColOnly0, tempSet)
  modelSetColOnly1 = c(modelSetColOnly1, tempSet)
  stlVecColOnly1 = unlist(lapply(lapply(modelSetColOnly1, numextract_MH), length))
  
  newSet = modelSetColOnly1[which(!modelSetColOnly1 %in% modelSetColOnly0)]
  newSetReduced = newSet[-which.max(stringr::str_length(newSet))]
  
  stlVec = unlist(lapply(lapply(newSet, numextract_MH), length))
  
  wKeep = which(stlVec <= intMax)
  
  stlVec = stlVec[wKeep]
  newSet = newSet[wKeep]
  
  lSet = length(newSet)
  
  considerModels = list()
  considerModels[[1]] = sort(modelSetColOnly0)
  
  ## Keep track of the column that would need to be updated to get back to original matrix
  whichColUpdate = c()
  
  modToZeta = modelToZeta_MH(considerModels[[1]], k=k)
  isValid = rep(NA, ncol(modToZeta))
  for (jj in 1 : ncol(modToZeta)) {
    isValid[jj] = 1*(all(modToZeta[,jj] == tempZeta0[,dim2]))
  }
  whichColUpdate = c(whichColUpdate, which(isValid == 1)[1])
  
  ## only add in the model for tempZeta1 if it is different
  if (length(modelSetColOnly0) != length(modelSetColOnly1) &
      all(stlVecColOnly1 <= intMax)) {
    considerModels[[2]] = sort(modelSetColOnly1)
    
    modToZeta = modelToZeta_MH(considerModels[[2]], k=k)
    isValid = rep(NA, ncol(modToZeta))
    for (jj in 1 : ncol(modToZeta)) {
      isValid[jj] = 1*(all(modToZeta[,jj] == tempZeta1[,dim2]))
    }
    whichColUpdate = c(whichColUpdate, which(isValid == 1)[1])
  }
  
  numModelsTemp = length(considerModels)
  
  if (lSet > 0) {
    
    subsetTerms = list()
    for (j in 1 : lSet) {
      tempElements = as.numeric(numextract_MH(newSet[j]))
      subsetTerms[[j]] = which(sapply(lapply(lapply(newSet, numextract_MH), as.numeric), 
                                      function(x) return(all(x %in% tempElements))) &
                                 (1 : lSet) != j)
    }
    
    
    counter = numModelsTemp + 1
    for (jj in lSet : 1) {
      jjActive = which(!1:lSet %in% as.numeric(subsetTerms[[jj]]) & 1:lSet < jj)
      
      if (length(jjActive) == 0) {
        proposedSet = unique(c(modelSetColOnly0, newSet[[jj]]))
        
        addZeta = matrix(0, p,length(proposedSet))
        for (j in 1 : length(proposedSet)) {
          w = as.numeric(numextract_MH(proposedSet[j]))
          addZeta[w,j] = 1
        }
        zetaPlusZeta = cbind(tempZetaColOnly0, addZeta)
        
        newModel = as.character(unique(unlist(ActivePred_MH(zetaPlusZeta, k=dim(zetaPlusZeta)[2]))))
        
        keep = TRUE
        for (j in 1 : length(considerModels)) {
          if (length(newModel) == length(considerModels[[j]]) & 
              all(newModel %in% considerModels[[j]])) {
            keep = FALSE
          }
        }
        if (keep == TRUE) {
          considerModels[[counter]] = sort(newModel)
          counter = counter + 1
        }
      } else {
        
        tempGrid = expand.grid(rep(list(c(0,1)), length(jjActive)))
        
        for (ii in 1 : nrow(tempGrid)) {
          proposedSet = unique(c(modelSetColOnly0, newSet[[jj]], 
                                 newSet[jjActive[which(tempGrid[ii,] == 1)]]))
          
          addZeta = matrix(0, p,length(proposedSet))
          for (j in 1 : length(proposedSet)) {
            w = as.numeric(numextract_MH(proposedSet[j]))
            addZeta[w,j] = 1
          }
          zetaPlusZeta = cbind(tempZetaColOnly0, addZeta)
          
          newModel = as.character(unique(unlist(ActivePred_MH(zetaPlusZeta, k=dim(zetaPlusZeta)[2]))))
          
          keep = TRUE
          for (j in 1 : length(considerModels)) {
            if (length(newModel) == length(considerModels[[j]]) & 
                all(newModel %in% considerModels[[j]])) {
              keep = FALSE
            }
          }
          if (keep == TRUE) {
            considerModels[[counter]] = sort(newModel)
            counter = counter + 1
          }
        }
        
      }
      
    }
  }
  
  if (length(considerModels) > numModelsTemp) {
    for (j in (numModelsTemp + 1) : length(considerModels)) {
      modToZeta = modelToZeta_MH(considerModels[[j]], k=k)
      isValid = rep(NA, ncol(modToZeta))
      for (jj in 1 : ncol(modToZeta)) {
        isValid[jj] = 1*(all(modToZeta[-dim1,jj] == tempZeta0[-dim1,dim2]))
      }
      whichColUpdate = c(whichColUpdate, which(isValid == 1)[1])
    }
  }
  
  l = list(considerModels = considerModels,
           whichColUpdate = whichColUpdate)
  
  
  return(l)
}

modelToZeta_MH = function(model, k) {
  if (length(model) == 1) {
    finalZeta = matrix(0, p, k)
  } else {
    tempZeta = matrix(0, p, length(model)+2)
    if (length(model) > 1) {
      for (j in 2 : length(model)) {
        activeExposures = as.numeric(numextract_MH(model[j]))
        tempZeta[activeExposures,j-1] = 1
      }
    }
    
    ## Now loop through and remove unnecessary exposures
    for (j1 in 1 : (length(model)+2)) {
      for (j2 in 1 : (length(model)+2)) {
        w1 = which(tempZeta[,j1] == 1)
        w2 = which(tempZeta[,j2] == 1)
        
        if (all(w1 %in% w2) & length(w1) > 0 & j1 != j2) tempZeta[,j1] = 0
      }
    }
    
    ## Now move features to the left and make the matrix p by k
    for (j in 1 : (length(model) + 2)) {
      if (sum(tempZeta[,j]) > 0) {
        wZ = which(apply(tempZeta, 2, sum) == 0)[1]
        tempZeta[,wZ] = tempZeta[,j]
        tempZeta[,j] = 0
      }
    }
    
    wNon = which(apply(tempZeta, 2, sum) > 0)
    
    if (length(wNon) > k) {
      print("Increase k, the number of model components")
      
      finalZeta = matrix(0, p, k)
      finalZeta[, 1:k] = tempZeta[,wNon[1:k]]
    } else {
      
      finalZeta = matrix(0, p, k)
      finalZeta[, 1:length(wNon)] = tempZeta[,wNon]
    }
  }
  
  return(finalZeta)
}

zetaPosterior_MH = function(Y, zeta, model, f_jhi_nc, betaC, sigmaP, tau,
                         k, sigB, Xstar, designC, ns, intMax) {
  
  PZ = -Inf
  tempY = Y - (designC %*% betaC)
  
  tauMat = t(matrix(rep(tau, p), k,p))
  tauMat = tauMat*zeta + (1-tauMat)*(1-zeta)
  
  ## remove intercept from model
  model = model[-1]
  
  if (length(model) == 0) {
    PZ = sum(log(1-tau))*p
    betaDraw = c(0)
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
    
    if(any(unlist(lapply(modelNumeric2, length)) > intMax)) {
      PZ=-Inf
      betaDraw = c(0)
      f_jhi_nc = rep(0,n)
    } else {
      
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
      
      SigmaB = diag(sigmaP*sigB, dim(tempXstar2)[2])
      muB = rep(0, dim(tempXstar2)[2])
      muBeta = solve((t(tempXstar2) %*% tempXstar2)/
                       sigmaP + solve(SigmaB)) %*%
        ((t(tempXstar2) %*% tempY)/
           sigmaP + solve(SigmaB) %*% muB)
      covBeta = solve((t(tempXstar2) %*% tempXstar2)/
                        sigmaP + solve(SigmaB))
      PZ = mvtnorm::dmvnorm(rep(0, length(muB)), mean=muB, sigma=as.matrix(SigmaB), log=TRUE) - 
        mvtnorm::dmvnorm(rep(0, length(muB)), mean=muBeta, sigma=as.matrix(covBeta), log=TRUE) +
        sum(log(tauMat))
      
      ## Draw beta from its distribution, in case you use this model
      betaDraw = mvtnorm::rmvnorm(1, mean=muBeta, sigma=as.matrix(covBeta))
      
      f_jhi_nc = as.vector(tempXstar2 %*% as.vector(betaDraw))
    }
  }
  
  ll = list(PZ=PZ, betaDraw=betaDraw, f_jhi_nc=f_jhi_nc)
  
  return(ll)
  
}

updateBetaOne_MH = function(Y, zeta, f_jhi_nc, betaC, sigmaP, tau,
                         k, sigB, Xstar, designC, ns, intMax) {
  
  beta = c(0)
  
  p = dim(Xstar)[2]
  
  sxx = sample(1:p, min(p,4), replace=FALSE)
  
  for (sx in sxx) {
    ## sample column to look at
    model = sort(as.character(unique(unlist(ActivePred_MH(zeta=zeta, k=k)))))
    
    wNonZero = which(apply(zeta, 2, sum) > 0)
    wZero = which(apply(zeta, 2, sum) == 0)
    
    hx = sample(c(wNonZero, wZero[1]), 1)
    
    CM2 = ChooseModels2_MH(zeta=zeta, dim1=sx, dim2=hx, k=k, intMax=intMax)
    PossibleModels = CM2$considerModels
    whichColUpdate = CM2$whichColUpdate
    PossibleZetas = lapply(PossibleModels, modelToZeta_MH, k=k)
    
    if (length(PossibleModels) == 1) {
      tempResults = zetaPosterior_MH(Y=Y, zeta=PossibleZetas[[1]], model=PossibleModels[[1]], 
                                  f_jhi_nc = f_jhi_nc, betaC=betaC, 
                                  sigmaP = sigmaP, tau = tau,
                                  k=k, sigB=sigB, Xstar=Xstar, designC=designC, 
                                  ns=ns, intMax=intMax)
      
      zeta = PossibleZetas[[1]]
      beta = tempResults$betaDraw
      f_jhi_nc = tempResults$f_jhi_nc
    } else {
      wOrig = which(unlist(lapply(PossibleModels, 
                                  function(x) return(paste(x, collapse="--")))) == 
                      paste(model, collapse="--"))
      
      ## Find which model to jump to and probability of going there
      if (length(PossibleModels) == 2) {
        wNew = (1:2)[-wOrig]
        pOldToNew = 1
      } else {
        wNew = sample((1:length(PossibleModels))[-wOrig], 1)
        pOldToNew = 1 / (length(PossibleModels) -1)
      }
      
      hx2 = whichColUpdate[wNew]
      
      ## Find probability of coming back to current model from new model
      PossibleModelsNewToOld = ChooseModels2_MH(zeta=PossibleZetas[[wNew]], 
                                             dim1=sx, dim2=hx2, k=k,
                                             intMax=intMax)$considerModels
      PossibleZetasNewToOld = lapply(PossibleModelsNewToOld, modelToZeta_MH, k=k)
      wNewOld = which(unlist(lapply(PossibleModelsNewToOld, 
                                    function(x) return(paste(x, collapse="--")))) == 
                        paste(model, collapse="--"))
      wNewNew = which(unlist(lapply(PossibleModelsNewToOld, 
                                    function(x) return(paste(x, collapse="--")))) == 
                        paste(PossibleModels[[wNew]], collapse="--"))
      
      if (length(wNewOld) == 0) {
        pNewToOld = 0
      } else {
        pNewToOld = 1 / (length(PossibleModelsNewToOld)-1)
      }
      
      ## Now do MH to get probability of accepting a move
      
      tempResultsOld = zetaPosterior_MH(Y=Y, zeta=PossibleZetas[[wOrig]], 
                                     model=PossibleModels[[wOrig]], 
                                     f_jhi_nc = f_jhi_nc, betaC=betaC, 
                                     sigmaP = sigmaP, tau = tau,
                                     k=k, sigB=sigB, Xstar=Xstar, designC=designC, 
                                     ns=ns, intMax=intMax)
      tempResultsNew = zetaPosterior_MH(Y=Y, zeta=PossibleZetas[[wNew]], 
                                     model=PossibleModels[[wNew]], 
                                     f_jhi_nc = f_jhi_nc, betaC=betaC, 
                                     sigmaP = sigmaP, tau = tau,
                                     k=k, sigB=sigB, Xstar=Xstar, designC=designC, 
                                     ns=ns, intMax=intMax)
      
      maxlog = max(tempResultsOld$PZ, tempResultsNew$PZ)
      updated_p = exp(c(tempResultsOld$PZ, tempResultsNew$PZ) - maxlog)
      
      move = "reject"
      if (pNewToOld == 0) {
        move = "reject"
      } else if (updated_p[1] == 0) {
        move = "accept"
      } else {
        ratio = updated_p[2]*pNewToOld / (updated_p[1]*pOldToNew)
        unif = runif(1)
        if (ratio > unif) {
          move = "accept"
        }
      }
      
      ## accept or reject
      if (move == "accept") {
        zeta = PossibleZetas[[wNew]]
        beta = tempResultsNew$betaDraw
        f_jhi_nc = tempResultsNew$f_jhi_nc
      } else {
        zeta = PossibleZetas[[wOrig]]
        beta = tempResultsOld$betaDraw
        f_jhi_nc = tempResultsOld$f_jhi_nc
      }
    }
  }
  
  ll = list(zeta=zeta, beta=beta, f_jhi_nc=f_jhi_nc)
  
  return(ll)
}

updateBetaTwo_MH = function(Y, zeta, f_jhi_nc, betaC, sigmaP, tau,
                         k, sigB, Xstar, designC, ns, intMax) {
  
  beta = c(0)
  
  p = dim(Xstar)[2]
  
  nUpdates = min(floor(p/2),4)
  
  sxx = sample(1:p, nUpdates*2, replace=FALSE)
  
  for (nu in 1 : nUpdates) {
    ## sample column to look at
    model = sort(as.character(unique(unlist(ActivePred_MH(zeta=zeta, k=k)))))
    
    wNonZero = which(apply(zeta, 2, sum) > 0)
    wZero = which(apply(zeta, 2, sum) == 0)
    
    hx = sample(c(wNonZero, wZero[1]), 1)
    
    ## sample exposures to look at
    sx = sxx[((nu-1)*2 + 1) : (nu*2)]
    
    CM2 = ChooseModelsTwo2_MH(zeta=zeta, dim1=sx, dim2=hx, k=k, intMax=intMax)
    PossibleModels = CM2$considerModels
    whichColUpdate = CM2$whichColUpdate
    PossibleZetas = lapply(PossibleModels, modelToZeta_MH, k=k)
    
    if (length(PossibleModels) == 1) {
      tempResults = zetaPosterior_MH(Y=Y, zeta=PossibleZetas[[1]], model=PossibleModels[[1]], 
                                  f_jhi_nc = f_jhi_nc, betaC=betaC, 
                                  sigmaP = sigmaP, tau = tau,
                                  k=k, sigB=sigB, Xstar=Xstar, designC=designC, 
                                  ns=ns, intMax=intMax)
      
      zeta = PossibleZetas[[1]]
      beta = tempResults$betaDraw
      f_jhi_nc = tempResults$f_jhi_nc
    } else {
      wOrig = which(unlist(lapply(PossibleModels, 
                                  function(x) return(paste(x, collapse="--")))) == 
                      paste(model, collapse="--"))
      
      ## Find which model to jump to and probability of going there
      if (length(PossibleModels) == 2) {
        wNew = (1:2)[-wOrig]
        pOldToNew = 1
      } else {
        wNew = sample((1:length(PossibleModels))[-wOrig], 1)
        pOldToNew = 1 / (length(PossibleModels) -1)
      }
      
      hx2 = whichColUpdate[wNew]
      
      ## Find probability of coming back to current model from new model
      PossibleModelsNewToOld = ChooseModelsTwo2_MH(zeta=PossibleZetas[[wNew]], 
                                             dim1=sx, dim2=hx2, k=k,
                                             intMax=intMax)$considerModels
      PossibleZetasNewToOld = lapply(PossibleModelsNewToOld, modelToZeta_MH, k=k)
      wNewOld = which(unlist(lapply(PossibleModelsNewToOld, 
                                    function(x) return(paste(x, collapse="--")))) == 
                        paste(model, collapse="--"))
      wNewNew = which(unlist(lapply(PossibleModelsNewToOld, 
                                    function(x) return(paste(x, collapse="--")))) == 
                        paste(PossibleModels[[wNew]], collapse="--"))
      
      if (length(wNewOld) == 0) {
        pNewToOld = 0
      } else {
        pNewToOld = 1 / (length(PossibleModelsNewToOld)-1)
      }
      
      ## Now do MH to get probability of accepting a move
      
      tempResultsOld = zetaPosterior_MH(Y=Y, zeta=PossibleZetas[[wOrig]], 
                                     model=PossibleModels[[wOrig]], 
                                     f_jhi_nc = f_jhi_nc, betaC=betaC, 
                                     sigmaP = sigmaP, tau = tau,
                                     k=k, sigB=sigB, Xstar=Xstar, designC=designC, 
                                     ns=ns, intMax=intMax)
      tempResultsNew = zetaPosterior_MH(Y=Y, zeta=PossibleZetas[[wNew]], 
                                     model=PossibleModels[[wNew]], 
                                     f_jhi_nc = f_jhi_nc, betaC=betaC, 
                                     sigmaP = sigmaP, tau = tau,
                                     k=k, sigB=sigB, Xstar=Xstar, designC=designC, 
                                     ns=ns, intMax=intMax)
      
      maxlog = max(tempResultsOld$PZ, tempResultsNew$PZ)
      updated_p = exp(c(tempResultsOld$PZ, tempResultsNew$PZ) - maxlog)
      
      move = "reject"
      if (pNewToOld == 0) {
        move = "reject"
      } else if (updated_p[1] == 0) {
        move = "accept"
      } else {
        ratio = updated_p[2]*pNewToOld / (updated_p[1]*pOldToNew)
        unif = runif(1)
        if (ratio > unif) {
          move = "accept"
        }
      }
      
      ## accept or reject
      if (move == "accept") {
        zeta = PossibleZetas[[wNew]]
        beta = tempResultsNew$betaDraw
        f_jhi_nc = tempResultsNew$f_jhi_nc
      } else {
        zeta = PossibleZetas[[wOrig]]
        beta = tempResultsOld$betaDraw
        f_jhi_nc = tempResultsOld$f_jhi_nc
      }
    }
  }
  
  ll = list(zeta=zeta, beta=beta, f_jhi_nc=f_jhi_nc)
  
  return(ll)
}


MCMCmixture_MH = function(Y, X, C, Xstar, nChains = 2, nIter = 30000, nBurn = 10000, thin = 20,
                       c = 0.001, d = 0.001, sigB, muB, SigmaC, muC,
                       k = 10, ns = 3, alph, gamm, intMax=3, probSamp1=0.95) {
  n = dim(X)[1]
  p = dim(X)[2]
  
  designC = cbind(rep(1,n), C)
  
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
  f_jhi = array(NA, dim=c(nChains,n))
  for (nc in 1 : nChains) {
    f_jhi[nc,] = 0
  }
  
  betaList = list()
  betaList[[1]] = list()
  for (nc in 1 : nChains) {
    betaList[[1]][[nc]] = c(0)
  }
  
  ## begin MCMC
  for (ni in 2 : nIter) {
    betaList[[ni]] = list()
    for (nc in 1 : nChains) {
      
      betaList[[ni]][[nc]] = c()
      
      if(ni %% 100 == 0 & nc == 1)  {
        print(paste(ni, " MCMC scans have finished", sep=""))
      }
      
      alphaPost[nc,ni] = alph
      gammaPost[nc,ni] = gamm
      
      ## sample sigma squared
      tempSum0 = 0
      tempSum2 = 0
      tempSum0 = tempSum0 + sum(unlist(betaList[[ni-1]][[nc]]) != 0)
      tempSum2 = tempSum2 + sum(unlist(betaList[[ni-1]][[nc]])^2)
      SampleDiffSq = (Y - f_jhi[nc,] - 
                        (designC %*% betaCPost[nc,ni-1,]))^2
      sigmaPost[nc,ni] = 1/rgamma(1, c + n/2 + tempSum0/2, 
                                  d + sum(SampleDiffSq)/2 + tempSum2/(2*sigB))
      
      ## update tau_h
      for (h in 1 : k) {
        tauPost[nc,ni,h] = rbeta(1, (alphaPost[nc,ni]*gammaPost[nc,ni]/k) + 
                                   sum(zetaPost[nc,ni-1,,h]), 
                                 gammaPost[nc,ni] + sum(1 - zetaPost[nc,ni-1,,h]))
      }
      
      
      ## Update zeta and beta
      ## randomly choose whether to update 1 or 2 exposures at a time
      oneOrTwo = sample(1:2, 1, prob = c(probSamp1, 1-probSamp1))
      
      if (oneOrTwo == 1) {
        updateBetaZeta = updateBetaOne_MH(Y=Y, zeta = zetaPost[nc,ni-1,,], f_jhi_nc=f_jhi[nc,], 
                                       betaC = betaCPost[nc,ni-1,], sigmaP=sigmaPost[nc,ni], 
                                       tau=tauPost[nc,ni,], k=k, sigB=sigB, Xstar=Xstar, 
                                       designC=designC, ns=ns, intMax=intMax)
      } else {
        updateBetaZeta = updateBetaTwo_MH(Y=Y, zeta = zetaPost[nc,ni-1,,], f_jhi_nc=f_jhi[nc,], 
                                       betaC = betaCPost[nc,ni-1,], sigmaP=sigmaPost[nc,ni], 
                                       tau=tauPost[nc,ni,], k=k, sigB=sigB, Xstar=Xstar, 
                                       designC=designC, ns=ns, intMax=intMax)
      }
      
      betaList[[ni]][[nc]] = updateBetaZeta$beta
      zetaPost[nc,ni,,] = updateBetaZeta$zeta
      f_jhi[nc,] = updateBetaZeta$f_jhi_nc
      
      
      ## Update betaC
      tempY = Y - f_jhi[nc,]
      
      muBetaC = solve((t(designC) %*% designC)/sigmaPost[nc,ni-1] + solve(SigmaC)) %*% 
        ((t(designC) %*% tempY)/sigmaPost[nc,ni-1] + solve(SigmaC) %*% muC)
      covBetaC = solve((t(designC) %*% designC)/sigmaPost[nc,ni-1] + solve(SigmaC))
      
      betaCPost[nc,ni,] = mvtnorm::rmvnorm(1, mean=muBetaC, sigma = covBetaC)
    }
  }
  
  keep = nBurn + (1:floor((nIter - nBurn) / thin))*thin
  
  posterior = list(sigma = sigmaPost[,keep],
                   tau = tauPost[,keep,],
                   zeta = zetaPost[,keep,,],
                   beta = betaList[keep],
                   betaC = betaCPost[,keep,])
  
  return(posterior)
}

MCMCmixtureMinSig_MH = function(Y, X, C, Xstar, nPerms = 10, nIter = 500,
                             c = 0.001, d = 0.001, sigBstart = 1, muB, SigmaC, muC,
                             k = 10, ns = 2, threshold=0.25) {
  
  n = dim(X)[1]
  p = dim(X)[2]
  
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
          
          distinguish = TRUE
          counter = 1
          while (distinguish == TRUE) {
            
            sigB = .75^counter
            
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


MCMCmixtureEB_MH = function(Y, X, C, Xstar, nChains = 2, nIter = 30000, nBurn = 10000, thin = 20,
                         c = 0.001, d = 0.001, sigBstart, muB, SigmaC, muC,
                         k = 10, ns = 3, alph, gamm, intMax=3, SigMin = 0.1,
                         probSamp1 = 0.95) {
  n = dim(X)[1]
  p = dim(X)[2]
  
  designC = cbind(rep(1,n), C)
  
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
  f_jhi = array(NA, dim=c(nChains,n))
  for (nc in 1 : nChains) {
    f_jhi[nc,] = 0
  }
  
  betaList = list()
  betaList[[1]] = list()
  for (nc in 1 : nChains) {
    betaList[[1]][[nc]] = c(0)
  }
  
  ## begin MCMC
  for (ni in 2 : nIter) {
    betaList[[ni]] = list()
    for (nc in 1 : nChains) {
      
      betaList[[ni]][[nc]] = c()
      
      if(ni %% 100 == 0 & nc == 1)  {
        print(paste(ni, " MCMC scans have finished", sep=""))
      }
      
      alphaPost[nc,ni] = alph
      gammaPost[nc,ni] = gamm
      
      ## sample sigma squared
      tempSum0 = 0
      tempSum2 = 0
      tempSum0 = tempSum0 + sum(unlist(betaList[[ni-1]][[nc]]) != 0)
      tempSum2 = tempSum2 + sum(unlist(betaList[[ni-1]][[nc]])^2)
      SampleDiffSq = (Y - f_jhi[nc,] - 
                        (designC %*% betaCPost[nc,ni-1,]))^2
      sigmaPost[nc,ni] = 1/rgamma(1, c + n/2 + tempSum0/2, 
                                  d + sum(SampleDiffSq)/2 + tempSum2/(2*sigB))
      
      ## update tau_h
      for (h in 1 : k) {
        tauPost[nc,ni,h] = rbeta(1, (alphaPost[nc,ni]*gammaPost[nc,ni]/k) + 
                                   sum(zetaPost[nc,ni-1,,h]), 
                                 gammaPost[nc,ni] + sum(1 - zetaPost[nc,ni-1,,h]))
      }
      
      
      ## Update zeta and beta
      oneOrTwo = sample(1:2, 1, prob = c(probSamp1, 1-probSamp1))
      
      if (oneOrTwo == 1) {
        updateBetaZeta = updateBetaOne_MH(Y=Y, zeta = zetaPost[nc,ni-1,,], f_jhi_nc=f_jhi[nc,], 
                                       betaC = betaCPost[nc,ni-1,], sigmaP=sigmaPost[nc,ni], 
                                       tau=tauPost[nc,ni,], k=k, sigB=sigB, Xstar=Xstar, 
                                       designC=designC, ns=ns, intMax=intMax)
      } else {
        updateBetaZeta = updateBetaTwo_MH(Y=Y, zeta = zetaPost[nc,ni-1,,], f_jhi_nc=f_jhi[nc,], 
                                       betaC = betaCPost[nc,ni-1,], sigmaP=sigmaPost[nc,ni], 
                                       tau=tauPost[nc,ni,], k=k, sigB=sigB, Xstar=Xstar, 
                                       designC=designC, ns=ns, intMax=intMax)
      }
      
      betaList[[ni]][[nc]] = updateBetaZeta$beta
      zetaPost[nc,ni,,] = updateBetaZeta$zeta
      f_jhi[nc,] = updateBetaZeta$f_jhi_nc
      
      
      if (ni %% 50 == 0) {
        sumsB = 0
        sumsZ = 0
        for (ni2 in (ni - 49) : ni) {
          sumsB = sumsB + sum(unlist(betaList[[ni2]][[nc]])^2)/sigmaPost[nc,ni2]
          sumsZ = sumsZ + length(which(unlist(betaList[[ni2]][[nc]]) != 0))
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
      tempY = Y - f_jhi[nc,]
      
      muBetaC = solve((t(designC) %*% designC)/sigmaPost[nc,ni-1] + solve(SigmaC)) %*% 
        ((t(designC) %*% tempY)/sigmaPost[nc,ni-1] + solve(SigmaC) %*% muC)
      covBetaC = solve((t(designC) %*% designC)/sigmaPost[nc,ni-1] + solve(SigmaC))
      
      betaCPost[nc,ni,] = mvtnorm::rmvnorm(1, mean=muBetaC, sigma = covBetaC)
      
    }
  }
  
  keep = nBurn + (1:floor((nIter - nBurn) / thin))*thin
  
  return(mean(sigBPost[,nIter]))
}


