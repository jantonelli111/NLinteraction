# NLinteraction

This is an R package to implement the ideas in Antonelli et. al (2018), which can be found at the following link:

https://arxiv.org/pdf/1711.11239.pdf

Please don't hesitate to contact Joseph Antonelli with any questions at jantonelli111@gmail.com. Please report any bugs if you encounter any as well!

# How to use NLinteraction

The first thing we need to do is install the library, which can be done using the following

```
library(devtools)
install_github(repo = "jantonelli111/NLinteraction")
library(NLinteraction)
```

So let's generate some data that we can use to illustrate the different functionalities in the package:

```
n = 100
p = 10
pc = 1

X = matrix(rnorm(n*p), n, p)

C = matrix(rnorm(n*pc), nrow=n)

TrueH = function(X) {
  return(1.5*(X[,2]*X[,3]) - 1.6*(X[,4]^2 * X[,5]))
}

Y = 5 + C + TrueH(X) + rnorm(n)
```

So we've generated p=10 exposures that are of interest for this software, one additional covariate that we control for (typically demographic variables, etc.), and an outcome conditional on both of these things. We see that the true relationship between X and Y is nonlinear as there is a squared term. We have intentionally chosen strong signals for this illustration, so that we can clearly see in all of the plots what is happening.

### Main MCMC function

Now we will build our model. The first thing to note is that the number of degrees of freedom of the splines in the model is very important and therefore the model should be built for a few different values of the degrees of freedom and then the WAIC can be used to evaluate which model to proceed with. Here we will choose 2 and 5 degrees of freedom, and we know that the model with 2 should be preferred. Degrees of freedom is controlled by the ns parameter in the model. We also must select the number of MCMC iterations we want (nIter), the number to throw out with the burn-in (nBurn), how many samples to thin (thin), and how many chains to run (nChains). Other parameters that control things such as the hyperparameters of the prior distributions in the model can be changed as well, but we recommend using the pre-set values unless their is strong prior information and the user fully understands the prior distributions.

```
NLmod2 = NLint(Y=Y, X=X, C=C, nIter=1000, nBurn=500, thin=2, nChains=2, ns=2)
NLmod5 = NLint(Y=Y, X=X, C=C, nIter=1000, nBurn=500, thin=2, nChains=2, ns=5)

```


**References**

Joseph Antonelli, Maitreyi Mazumdar, David Bellinger, David C. Christiani, Robert Wright, Brent A. Coull. **Bayesian variable selection for multi-dimensional semiparametric regression models**. 2018. arXiv:1711.11239
