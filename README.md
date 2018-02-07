# NLinteraction

This is an R package to implement the ideas in Antonelli et. al (2018), which can be found at the following link:

https://arxiv.org/pdf/1711.11239.pdf

Please don't hesitate to contact Joseph Antonelli with any questions at jantonelli111@gmail.com. Please report any bugs if you encounter any as well!

Below is a quick example of how to use the programs:

```{r, eval=FALSE, message=FALSE}
library(devtools)
install_github(repo = "jantonelli111/NLinteraction")
library(NLinteraction)

n = 100
p = 10
pc = 1

X = matrix(rnorm(n*p), n, p)

C = matrix(rnorm(n*pc), nrow=n)

TrueH = function(X) {
  return(1.5*(X[,2]*X[,3]) - 1.6*(X[,4]^2 * X[,5]))
}

Y = 5 + C + TrueH(X) + rnorm(n)

NLmod = NLint(Y=Y, X=X, C=C, nIter=20, nBurn=10, thin=2, nChains=2)

## Print posterior inclusion probabilities
NLmod$MainPIP

## Show the two way interaction probability matrix
NLmod$InteractionPIP
```
