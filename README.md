
<!-- README.md is generated from README.Rmd. Please edit that file -->

# tsbugs

<!-- badges: start -->
<!-- badges: end -->
<!-- The goal of tsbugs is to ... -->

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("guyabel/tsbugs")
```

The package is no longer on CRAN.

## Details

The functions in the tsbugs package are aimed to automate the writing of
time series models to run in
<a href="http://www.mrc-bsu.cam.ac.uk/bugs/welcome.shtml">WinBUGS</a> or
<a href="http://www.openbugs.net/w/FrontPage">OpenBUGS</a>. I created
these functions when working on [model averaging for time series
models](https://www.demographic-research.org/volumes/vol29/43/default.htm).
I found it a lot easier to build R functions to write the BUGS models
than the more error-inducing process of copy and pasting BUGS scripts,
and then making slight alterations to create new models. It also allowed
me to add arguments to specify different lag lengths, prior
distributions, variance assumptions and data lengths. Below are examples
for three types of time series models; autorgressive models with

-   [Constant
    variance](https://github.com/guyabel/tsbugs#autoregressive-models)
-   [Stochastic
    volatility](https://github.com/guyabel/tsbugs#stochastic-volatility-models)
-   [Random variance shift
    models](https://github.com/guyabel/tsbugs#random-variance-shift-models)

## Autoregressive Models

The `ar.bugs` command builds a BUGS script for autoregressive (AR)
models ready to use in
<a href="http://cran.r-project.org/web/packages/R2OpenBUGS">R2OpenBUGS</a>.
For example, consider the `LakeHuron` data.

``` r
LH <- LakeHuron
par(mfrow=c(2,1))
plot(LH, main="Level (ft)")
plot(diff(LH), main="Differenced Level")
```

<a href="http://gjabel.files.wordpress.com/2014/03/tsbugs1-cvdata.png"><img class="aligncenter size-full wp-image-1219" src="http://gjabel.files.wordpress.com/2014/03/tsbugs1-cvdata.png" alt="tsbugs1-cvdata" width="600" height="600" /></a><br />
We can construct a AR(1) model for this data (after differencing the
data to obtain a stationary mean) as such:

``` r
library(tsbugs)
#> Welcome to the tsbugs package
#> Plesae cite as:
#> Abel, G.J., Bijak, J., Forster, J.J., Raymer J., Smith P.W.F. and Wong J.S.T. (2013). Integrating uncertainty in time series population forecasts: An illustration using a simple projection model. Demographic Research. 43 (29) 187--1226
ar1 <- ar.bugs(y=diff(LH), ar.order=1)
print(ar1$bug)
#>  [1] "model{"                              "#likelihood"                        
#>  [3] "for(t in 2:97){"                     "\ty[t] ~ dnorm(y.mean[t], isigma2)"  
#>  [5] "}"                                   "for(t in 2:97){"                    
#>  [7] "\ty.mean[t] <- phi1*y[t-1]"           "}"                                  
#>  [9] "#priors"                             "phi1 ~ dnorm(0,1)"                  
#> [11] "isigma2 ~ dgamma(0.000001,0.000001)" "sigma <- pow(isigma2,-0.5)"         
#> [13] "}"
```

The `ar.bugs` function allows for alternative specifications for prior
distributions, forecasts and the inclusion of mean term:

``` r
ar2 <- ar.bugs(y=diff(LH), ar.order=2, ar.prior="dunif(-1,1)", var.prior="dgamma(0.001,0.001)",
               k = 10, mean.centre = TRUE)
print(ar2$bug)
#>  [1] "model{"                                                      
#>  [2] "#likelihood"                                                 
#>  [3] "for(t in 3:107){"                                            
#>  [4] "\ty[t] ~ dnorm(y.mean[t], isigma2)"                           
#>  [5] "}"                                                           
#>  [6] "for(t in 3:107){"                                            
#>  [7] "\ty.mean[t] <- phi0 + phi1*(y[t-1]-phi0) + phi2*(y[t-2]-phi0)"
#>  [8] "}"                                                           
#>  [9] "#priors"                                                     
#> [10] "phi0 ~ dunif(-1,1)"                                          
#> [11] "phi1 ~ dunif(-1,1)"                                          
#> [12] "phi2 ~ dunif(-1,1)"                                          
#> [13] "sigma2 ~ dgamma(0.001,0.001)"                                
#> [14] "isigma2 <- pow(sigma2,-1)"                                   
#> [15] "#forecast"                                                   
#> [16] "for(t in 98:107){"                                           
#> [17] "\ty.new[t] <- y[t]"                                           
#> [18] "}"                                                           
#> [19] "}"
```

The tsbugs objects can be used with R2OpenBUGS to easily run models from
R. This is made even easier using the `inits` and `nodes` functions
(also in the tsbugs package). For example:

``` r
writeLines(ar2$bug, "ar2.txt")
library("R2OpenBUGS")
ar2.bug <- bugs(data = ar2$data,
                inits = list(inits(ar2)),
                param = c(nodes(ar2, "prior")$name, "y.new"),
                model = "ar2.txt",
                n.iter = 11000, n.burnin = 1000, n.chains = 1)
```

Note, 1) the model is written to a `.txt` file (as required by
R2OpenBUGS), 2) the data used is part of the `tsbugs` object. The
`ar.bugs` command cleans the data and adds missing values at the end of
the series for foretasted values, 3) the initial values offered by the
inits function are very crude, and with more complicated data or models,
users might be better off specifying there own list of initial values.
The parameter traces and posterior distributions can be plotted using
the <a href="http://cran.r-project.org/web/packages/coda">coda</a>
package:

``` r
library(coda)
param.mcmc <- as.mcmc(ar2.bug$sims.matrix[,nodes(ar2, "prior")$name])
plot(param.mcmc[,1:4])
```

<a href="http://gjabel.files.wordpress.com/2014/03/tsbugs1-cvparam.png"><img class="aligncenter size-full wp-image-1221" src="http://gjabel.files.wordpress.com/2014/03/tsbugs1-cvparam.png" alt="tsbugs1-cvparam" width="600" height="600" /></a><br />

The fanplot package can be used to plot the entire series of posterior
predictive distributions. We may also plot (after deriving using the
`diffinv` function) the posterior predictive distributions of the lake
level:

``` r
# derive future level
ynew.mcmc <- ar2.bug$sims.list$y.new
lhnew.mcmc <- apply(ynew.mcmc, 1, diffinv, xi = tail(LH,1))
lhnew.mcmc <- t(lhnew.mcmc[-1,])

# plot differenced
par(mfrow=c(2,1))
plot(diff(LH), xlim = k0 + c(-50, 10), main="Differenced Level")

# add fan
library("fanplot")
k0 <- end(LH)[1]
fan(ynew.mcmc, start=k0+1, rcex=0.5)

# plot undifferenced
plot(LH, xlim=k0+c(-50,10), main="Level")
fan(lhnew.mcmc, start=k0+1, rcex=0.5)
```

<a href="http://gjabel.files.wordpress.com/2014/03/tsbugs1-cvforc.png"><img class="aligncenter size-full wp-image-1220" src="http://gjabel.files.wordpress.com/2014/03/tsbugs1-cvforc.png" alt="tsbugs1-cvforc" width="600" height="600" /></a>

## Stochastic Volatility Models

The `sv.bugs` command builds a BUGS script for stochastic volatility SV
models ready to use in R2OpenBUGS. For example, consider the `svpdx`
data.

``` r
# plot
plot(svpdx$pdx, type = "l",
     main = "Return of Pound-Dollar exchange rate data from 2nd October 1981 to 28th June 1985", 
     cex.main = 0.8)
```

<a href="http://gjabel.files.wordpress.com/2014/03/tsbugs1-svdata.png"><img class="aligncenter size-full wp-image-1216" src="http://gjabel.files.wordpress.com/2014/03/tsbugs1-svdata.png" alt="tsbugs1-svdata" width="600" height="600" /></a><br />
We can construct a AR(0)-SV model for this data, and also obtain
posterior simulations using the `sv.bugs` command:

``` r
y <- svpdx$pdx
sv0 <- sv.bugs(y, sim=TRUE)
print(sv0$bug)
#>  [1] "model{"                                     
#>  [2] "#likelihood"                                
#>  [3] "for(t in 1:945){"                           
#>  [4] "\ty[t] ~ dnorm(y.mean[t], isigma2[t])"       
#>  [5] "\tisigma2[t] <- exp(-h[t])"                  
#>  [6] "\th[t] ~ dnorm(h.mean[t], itau2)"            
#>  [7] "}"                                          
#>  [8] "for(t in 1:945){"                           
#>  [9] "\ty.mean[t] <- 0"                            
#> [10] "}"                                          
#> [11] "for(t in 1:1){"                             
#> [12] "\th.mean[t] <- psi0"                         
#> [13] "}"                                          
#> [14] "for(t in 2:945){"                           
#> [15] "\th.mean[t] <- psi0 + psi1*(h[t-1]-psi0)"    
#> [16] "}"                                          
#> [17] "#priors"                                    
#> [18] "psi0 ~ dnorm(0,0.001)"                      
#> [19] "psi1.star ~ dunif(0,1)"                     
#> [20] "psi1 <- 2*psi1.star-1"                      
#> [21] "itau2 ~ dgamma(0.01,0.01)"                  
#> [22] "tau <- pow(itau2,-0.5)"                     
#> [23] "#simulation"                                
#> [24] "for(t in 1:945){"                           
#> [25] "\ty.mean.c[t] <- cut(y.mean[t])"             
#> [26] "\tisigma2.c[t] <- cut(isigma2[t])"           
#> [27] "\ty.sim[t] ~ dnorm(y.mean.c[t],isigma2.c[t])"
#> [28] "}"                                          
#> [29] "}"
```

This model closely matches those presented in
<a href="http://onlinelibrary.wiley.com/doi/10.1111/1368-423X.00046/abstract">Meyer
and Yu (2002)</a>. There are further options in the tsbugs package to
incorporate different priors that do not involve transformations such as
those for `psi1` above. Using R2OpenBUGS we can fit the model,

``` r
# decent initial value for variance in first period
init <- inits(sv0, warn=FALSE)
init$psi0 <- log(var(y))
# write bug
writeLines(sv0$bug, "sv0.txt")
# might take a while to compile
sv0.bug <- bugs(data = sv0$data,
                inits = list(init),
                param = c(nodes(sv0, "prior")$name,"y.sim","h"),
                model = "sv0.txt",
                n.iter = 11000, n.burnin = 1000, n.chains = 1)
```

The volatility and estimates can be easily extracted,

``` r
h.mcmc <- sv0.bug$sims.list$h
```

Which allows us to directly view the estimated volatility process or the
time-dependent standard deviation using the
<a href="https://guyabel.github.io/fanplot/">fanplot</a> package,

``` r
# plot
plot(NULL, xlim = c(1, 945)+c(0,40), ylim = c(-4,2), main="Estimated Volatility from SV Model")
# fan
fan(h.mcmc, type = "interval")
```

<a href="http://gjabel.files.wordpress.com/2014/03/tsbugs1-svh.png"><img class="aligncenter size-full wp-image-1217" src="http://gjabel.files.wordpress.com/2014/03/tsbugs1-svh.png" alt="tsbugs1-svh" width="600" height="600" /></a><br />We
can also plot the posterior simulations from the model:

``` r
# derive percentiles
y.mcmc <- sv0.bug$sims.list$y.sim
# plot
plot(NULL, type = "l", xlim = c(1, 945)+c(0,20), ylim = range(y), 
     main = "Posterior Model Simulations and Data")
fan(y.mcmc)
lines(y)
```

<a href="http://gjabel.files.wordpress.com/2014/03/tsbugs1-svysim.png"><img class="aligncenter size-full wp-image-1218" src="http://gjabel.files.wordpress.com/2014/03/tsbugs1-svysim.png" alt="tsbugs1-svysim" width="600" height="600" /></a>

## Random Variance Shift Models

The `rv.bugs` command builds a BUGS script for random variance (RV)
shift models, similar to that of
<a href="http://www.tandfonline.com/doi/abs/10.1080/01621459.1993.10476364#.UygPt4U7bIc">McCulloch
and Tsay (1993)</a> ready to use in R2OpenBUGS. Consider the `ew` data.

``` r
r <- ts(ew[2:167]/ew[1:166]-1, start=1841)
y <- diff(r)
plot(y, main="Difference in England and Wales Population Growth Rate")
```

<a href="http://gjabel.files.wordpress.com/2014/03/tsbugs1-rvdata.png"><img class="aligncenter size-full wp-image-1223" src="http://gjabel.files.wordpress.com/2014/03/tsbugs1-rvdata.png" alt="tsbugs1-rvdata" width="600" height="600" /></a><br />We
can create a BUGS script to fit a RV model to this data, including
posterior simulations, using the `rv.bugs` command:

``` r
rv0 <- rv.bugs(y, sim=TRUE)
print(rv0)
```

and then run the script in R2OpenBUGS (this can take a couple of hours):

``` r
# decent inital value for variance in first period
init <- inits(rv0, warn=FALSE)
init$isig02<-sd(y)^-2
# write bug
writeLines(rv0$bug,"rv0.txt")
# might take a while to compile
rv0.bug <- bugs(data = rv0$data,
                inits = list(init),
                param = c(nodes(rv0, "prior")$name,"y.sim",
                          "h","delta","beta"),
                model = "rv0.txt",
                n.iter = 11000, n.burnin = 1000, n.chains = 1)
```

We can plot the posterior simulations from the model using the fanplot
package:

``` r
# derive percentiles
y0 <- tsp(y)[1]
y.mcmc <- rv0.bug$sims.list$y.sim

# plot
plot(NULL, xlim=tsp(y)[1:2]+c(-5,5), ylim = range(y),
     main="Posterior Simulations")
fan(y.mcmc, start = y0, rlab=c(10,50,90), llab=TRUE)
lines(y)
```

<a href="http://gjabel.files.wordpress.com/2014/03/tsbugs1-rvysim.png"><img class="aligncenter size-full wp-image-1225" src="http://gjabel.files.wordpress.com/2014/03/tsbugs1-rvysim.png" alt="tsbugs1-rvysim" width="600" height="600" /></a><br />Alongside
the posterior distributions of the standard deviations,

``` r
# derive sigma
h.mcmc <- rv0.bug$sims.list$h
sigma.mcmc <- sqrt(exp(h.mcmc))
# plots
plot(NULL, xlim =tsp(y)[1:2]+c(-5,5), ylim = c(0,0.008), main="Standard Deviation")
fan(sigma.mcmc, start = y0, rlab=c(5,50,95), llab = c(5,50,95))
```

<a href="http://gjabel.files.wordpress.com/2014/03/tsbugs1-rvhsigma.png"><img class="aligncenter size-full wp-image-1224" src="http://gjabel.files.wordpress.com/2014/03/tsbugs1-rvhsigma.png" alt="tsbugs1-rvhsigma" width="600" height="600" /></a><br />The
posterior distributions of the probability of a variance shift and
multiplier effect of the shift in variance (`delta[t]` and `beta[t]` in
the BUGS model) can also be plotted. Note, when there is no variance
shift, the posterior of the `beta[t]` is similar to the prior
distribution.

``` r
#extract data
delta.mcmc <- rv0.bug$sims.list$delta
beta.mcmc <- rv0.bug$sims.list$beta

# plots
par(mfrow=c(2,1))
plot(NULL, xlim = tsp(y)[1:2]+c(-5,5), ylim = c(0,1), main="Probability of Variance Change Point")
fan(delta.mcmc, start=y0, ln = NULL, rlab = NULL)

plot(NULL, xlim = tsp(y)[1:2]+c(-5,5), ylim = c(-2,2), main="Variance Multiplier")
fan(beta.mcmc, start=y0)
```

<a href="http://gjabel.files.wordpress.com/2014/03/tsbugs1-rvbeta.png"><img class="aligncenter size-full wp-image-1222" src="http://gjabel.files.wordpress.com/2014/03/tsbugs1-rvbeta.png" alt="tsbugs1-rvbeta" width="600" height="600" /></a>
