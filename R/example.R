n = 100
p = 10
pc = 1

X = matrix(rnorm(n*p), n, p)

C = matrix(rnorm(n*pc), nrow=n)

TrueH = function(X) {
  return(1.5*(X[,2]*X[,3]) - 1.6*(X[,4]^2 * X[,5]))
}

Y = 5 + C + TrueH(X) + rnorm(n)

NLmod = NLint(Y=Y, X=X, C=C, nIter=1000, nBurn=200, thin=2, nChains=2, probSamp1 = 1)

## Print posterior inclusion probabilities
NLmod$MainPIP

## Show the two way interaction probability matrix
NLmod$InteractionPIP
