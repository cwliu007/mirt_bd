# Instructions: (Updating the code ...)
Instructions for running the models proposed in my paper (UNDER REVISION): *Multidimensional Item Response Theory Models for Testlet-Based Doubly Bounded Data*. 

# Step 1:
Install the following R packages: `install.packages(c("mvtnorm","randcorr","MASS","nimble","evaluate","Rcpp","coda"))`

Make sure you have installed RTools from https://cran.r-project.org/bin/windows/Rtools/.

# Step 2:
Download the two files: `multidimensional beta model correlations.R` and `sim_beta.cpp`.

Source the file: `source("multidimensional beta model correlations.R")`

----
Following are the arguments of the main function: `mcmc_main`

`response`:
* Item-response matrix (row: *sample size*, column: *number of items* by *number of testlets*). By default, `NULL`. 
* For instance, *number of items* = 2 and *number of testlets* = 3 means that the first two columns contain the item responses of the first testlet, and the second and third testlets are concatenated in column order. The resulting size of the `response` is *sample size* by *number of items x number of testlets*. 

`rand.seed`:
* Set the seed of R‘s random number generator. In effect, it is `set.seed(rand.seed)`. 

`testlet`:
* Number of testlets. (Used for simulating data)

`itemnum`:
* Number of items. (Used for simulating data)

`people`:
* Sample size of persons. (Used for simulating data)

`simulation`:
* `TRUE`: Carry out simulation study.
* `FALSE`: Parameter estimation of real data (it also returns `logProb_response` for calculating WAIC by `loo` package)

`return_sim = TRUE`:
* Return the simulated data and true values of the parameters

`burnin`:
* Number of initial, MCMC iterations to discard. Default value is 1.

`keep`:
* Number of iterations to run each MCMC chain after `burnin`. Default value is 2.

`nchains`:
* Number of MCMC chains to run. Default value is 1.

`model`:
 * model = 0 # logit-normal model (proposed model)
 * model = 1 # multidimensional beta model (without item-correlation parameters in the model)
 * model = 2 # beta copula model with item correlations (proposed model)
 * model = 3 # beta copula model with no item correlations (equivalent to `model = 1`)

`est_theta = TRUE`:
* Obtain the estimates of latent traits.

`sim_data = TRUE`:
* Simulate responses from posteriori of parameters

`use_real_data_pars`:
* For internal experiments, which can be omitted.

`use_real_data_pars_prior`:
* For internal experiments, which can be omitted.

`get_item_fit = TRUE`:
* Obtain outfit and infit for items. (Effective for model = 1 or 2)

`waic_nimble = TRUE`:
* Obtain the WAIC from NIMBLE.

# Step 3: Carry out an example
1. Simulate data:

```
source("multidimensional beta model correlations.R")
d = mcmc_main(rand.seed = 1
                ,testlet = 3
                ,itemnum = 2
                ,people = 1000
                ,simulation = TRUE
                ,return_sim=TRUE   # <--- return simulated data and true values of parameters
                ,model = 2)
```
`d` contains simulated data and true values of parameters.

```
out = mcmc_main(response = d$response  # <---- your data matrix
                ,rand.seed = 1
                ,testlet = 3
                ,itemnum = 2
                ,burnin = 1000         
                ,keep = 2000
                ,nchains = 3
                ,model = 2
                ,est_theta = TRUE
                ,sim_data = TRUE
                ,get_item_fit=TRUE
                ,waic_nimble=TRUE)
```
If you want to analyze your data, just replace `d$response` with your data. Make sure to change `testlet`, `itemnum`, and other arguments if needed. 

Check the non-convergence of MCMC chains by `gelman.diag` function from `coda` when `nchains` > 1:
```
gelman.diag(out$samples_items, confidence = 0.95, transform=FALSE, autoburnin=TRUE, multivariate=FALSE)
```

Model comparison criterion: `out$WAIC`

Item fit:
* outfit: `out$ppp_outfit_item`
* infit: `out$ppp_infit_item`

Compare the parameter estimates with true values if simulation study was carried out.
1. For BCM (the first column is the true values):

`alpha` (slope): `cbind(c(do.call(rbind,d$alpha)), out$est$alpha)`  

`delta` (location): `cbind(c(do.call(cbind,d$delta)), out$est$delta)`

`omega` (precision): `cbind(c(do.call(cbind,d$tau)), out$est$tau)`

`correlations between traits`: `cbind(c(d$sigma), out$est$corr)`

`correlations between item errors`: `cbind(do.call(c,d$sigma_item_cor), out$est$sigma_item_cor)`


2. For LNM (the first column is the true values):

`lambda` (slope): `cbind(do.call(c,d$alpha), out$est$alpha)`

`tau` (location): `cbind(do.call(c,d$delta), out$est$delta)`

`correlations between traits`: `cbind(c(d$sigma), out$est$corr)`

`variance-covariance matrix between item errors`: `cbind(do.call(c,d$sigma_item_cor), out$est$sigma_item_cor)`


# Notes:
If you have testlets that have different number of items, the NIMBLE code of `multidimensional beta model correlations.R` needs to be further modified accordingly.

