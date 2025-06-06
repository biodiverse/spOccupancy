# spOccupancy 0.8.1

+ Fixed a bug in `stMsPGOcc()` and `svcTMsPGOcc()` that prevented the models from running when there was only 1 within-season replicate. Thanks to Kevin Lafferty for bringing this to my attention. 
+ Fixed a bug in cross-validation of integrated models that led to the function crashing in certain situations, particularly on Windows machines. Thanks to Frederico Martins for bringing this to my attention.
+ Fixed a bug in the `predict()` function for `svcPGOcc()` (e.g., `predict.svcPGOcc()`) when using the `grid.index` functionality. When using `predict()`, a warning message would be displayed saying that there was a mismatch in the lengths of different vectors. This was a problem that resulted in inaccurate predictions when using `svcPGOcc()` and specifying the resolution of the SVCs as something higher than the point locations using the `grid.index` functionality. If you did this, I strongly encourage you to rerun the functions with v0.8.1, and I apologize for the inconvience this caused. Thanks very much to Ben Goldstein for pointing this out. Note that this bug did not affect results if using `svcPGOcc()` without specifying `grid.index` in the data list (which is likely the vast majority of instances for when `svcPGOcc()` is used). 

# spOccupancy 0.8.0

+ For queries on anything related to `spOccupancy` (and `spAbundance`), please use the new [spOccupancy/spAbundance mailing list](https://groups.google.com/g/spocc-spabund-users). 
+ New functionality for fitting multi-season, single-species integrated occupancy models. The function `tIntPGOcc()` fits a non-spatial multi-season integrated occupancy model, `stIntPGOcc()` fits a spatial multi-season integrated occupancy model, and `svcTIntPGOcc()` fits a spatially-varying coefficient multi-season occupancy model. Random intercepts are supported in both the occurrence and detection formulas for both model types. I am behind on adding vignettes for some of the newer functionality (sorry!), but adding a vignette for this new functionality is on my todo list. If interested in using these functions and you're having problems fitting them, please send your questions to the mailing list. 
+ Added in functionality for both occupancy and detection random intercepts in single-species single-season integrated models (`intPGOcc()` and `spIntPGOcc()`) using `lme4` syntax (e.g., `(1 | observer)` for a random effect of observer).
+ `simTIntPGOcc()` is a new function that allows simulation of single-species multi-season detection-nondetection data from multiple data sources.
+ Updated `simMsIntPGOcc()` to now include simulation of data sets with spatially-varying coefficients and unstructured random effects on both occurrence and detection.
+ Fixed a bug in the k-fold cross-validation for spatial integrated occupancy models (NNGP models only) that could lead to incorrect model deviance results under certain situations depending on how the spatial coordinates were ordered on the user-side relative to how they are re-ordered when fitting the model. If using `spIntPGOcc()` with `NNGP = TRUE` and using cross-validation, I suggest re-running the analysis. Apologies for the inconvenience.
+ Added in a `residuals()` function to extract occupancy and detection residuals following the approach of [Wilson et al. (2019)](https://esajournals.onlinelibrary.wiley.com/doi/abs/10.1002/ecy.2703) for single-season, single-species occupancy models (functions `PGOcc()`, `spPGOcc()`, and `svcPGOcc()`). I'm hoping to implement this for all model functions and improve GoF functionality. If anyone has any interest in helping out with this, then please let me know!
+ `waicOcc()` for integrated single-species models is now substantially faster.
+ `updateMCMC()` now works with lfJSDM.
+ Fixed a bug in `updateMCMC()` that prevented it from working with `spAbundance::msAbund()` when there were random effects in the model. Also added the `save.fitted` argument to `updateMCMC()` to allow it to work with `msAbund()` and not save the replicate/fitted data values in cases where the amount of RAM is an important consideration.
+ Added the `include.w` argument in the `predict()` function for `lfMsPGOcc()` models that enables predicting without the latent factors. This also allows prediction to occur without needing to supply the coordinates, which is useful when generating conditional probability plots.
+ Updated `lfJSDM()` to give an error more quickly when there are memory limitations.
+ Fixed a bug in all multi-season, multi-species models that caused the model to crash upon initialization of the MCMC algorithm when data were supplied in a way such that for a given data set, the maximum number of times a specific site was sampled was less than the total number of "replicate periods" (i.e., the fourth dimension of the data list). This may happen when the "replicates" are structured as specific time periods (i.e., weeks, years) instead of a specific "replicate". Thanks to José Ribeiro for bringing this to my attention. 
+ Fixed a bug in multi-species cross-validation that could cause an error when using a smaller number of threads for cross-validation compared to the number of folds used.

# spOccupancy 0.7.6

+ Fixed a memory problem in the saving of the tuning values for `svcTPGOcc` models that required updating v0.7.3 on CRAN to pass valgrind checks, as well as a memory leak in the calculation of the nearest neighbors, and a small problem with the DESCRIPTION file for including on CRAN. 

# spOccupancy 0.7.5

+ See v0.7.6

# spOccupancy 0.7.4

+ See v0.7.5

# spOccupancy 0.7.3

+ Fixed a problem that could arise when calculating Rhat in all models when running multiple chains (but usually only happened in multispecies models) when there was a high amount of correlation between parameter estimates. This would lead to the model running completely, but then failing after all chains have been run. This most often occurred when fitting a multispecies model with a lot of rare species. Thanks to Marc Kery for bringing this to my attention.
+ Added in a check at the top of all model fitting functions to return an error when the number of posterior samples saved based on the MCMC criteria (`n.batch`, `batch.length`, `n.samples`, `n.burn`, `n.thin`, `n.chains`) are specified in a way that leads to a non-integer value. In such situations, models would previously run and return without an error, but sometimes the last posterior sample in any given chain could have widely inaccurate values, or values that prevented subsequent functions from working. Thanks to Wendy Leuenberger and Colin Swider for bringing this to my attention. 
+ Added in functionality for fitting spatially-explicit models where the spatial random effects (or spatially varying coefficients) are not specified at the individual site, but rather are specified at a larger spatial resolution. This is accomplished using a new component of the `data` list supplied to model fitting functions called `grid.index`. This is useful for data sets where there is some sort of nested structuring among the data collection protocol, such that you may wish to specify the spatial random effects at a lower resolution than each individual location. Further, it can be particularly useful for SVC models where you only want to specify nonstationarity at a lower spatial resolution (e.g., across a set of grid cells). This is currently implemented for the following functions: `spPGOcc`, `sfMsPGOcc`, `stMsPGOcc`, `stPGOcc`, `svcPGOcc`, `svcTMsPGOcc`, `svcTPGOcc`. See the documentation for a given model function for how to specify this. I am hoping to eventually write up a small example that shows how to do this, but for now documentation is fairly limited to just the manual pages for each function. Feel free to contact me if you want to use this functionality and have any questions. 
+ Added in the `updateMCMC()` function. This function is in active development, but it will ultimately allow for all `spOccupancy` and `spAbundance` model objects to be updated with additional MCMC samples, instead of having to completely rerun an MCMC analysis if adequate burn-in/convergence was not reached. It currently works for the function `sfJSDM()` in `spOccupancy` and `msAbund()` in `spAbundance`. 
+ Added in the ability to specify independent priors for the species-level regression coefficients for two functions: `svcTMsPGOcc` and `sfJSDM`. This is done by setting the tags `independent.betas` and `independent.alphas` to TRUE. This will fix the values of the community-level mean and variance parameters to the initial values specified in `inits`. This is equivalent to setting an independent Gaussian prior for each of the species-specific regression coefficients, which may potentially be useful in certian situations where the assumption of normality in the distribution of the species-level effects is not well met. This functionality will eventually be incorporated for all multi-species models.
+ Fixed a bug in `intMsPGOcc()` that caused the model to crash upon initialization of the MCMC algorithm when data were supplied in a way such that for a given data set, the maximum number of times a specific site was sampled was less than the total number of "replicate periods" (i.e., the third dimension of the data list). This may happen when the "replicates" are structured as specific time periods (i.e., weeks, years) instead of a specific "replicate". This was previously fixed in all other model fitting functions.  
+ Wrote a new "vignette" (really more of a blog post) on some recommendations to help improve interpretability of inferences in SVC models.
+ Fixed a few typos in the MCMC sampler vignettes for factor models and SVC models.
+ Fixed a bug that prevented cross-validation from working properly in multi-species models when setting `k.fold.only = TRUE`. Thanks to Zack Steel for pointing this out.
+ Fixed a typo in the generation of initial values for latent unstructured random effects in all model functions. The typo had no major ramifications, if anything it would have just led to slower convergence, as it resulted in very large (or very small) initial values for the latent random effects that are not really viable on the logit scale.

# spOccupancy 0.7.2

+ Added in functionality for using the `plot()` function to generate simple traceplots using `spOccupancy` model objects. Details can be found in the help page (e.g., for `spPGOcc` models, type `?plot.spPGOcc` in the console).  
+ Not an update to the package, but a [new vignette](https://doserlab.com/files/spoccupancy-web/articles/identifiability) has been posted on testing model identifiability using `spOccupancy`. Thanks to Sara Stoudt for writing this!
+ Added in the ability to fit `lfJSDM()` without residual species correlations by setting `n.factors = 0`. This is a model analogous to `msPGOcc()`, but without the detection component. 
+ Added in the `shared.spatial` argument to `sfJSDM()`. If set to `TRUE`, this argument estimates a common spatial process for all species instead of using the default spatial factor modeling approach. 
+ Fixed a bug in `predict.svcTMsPGOcc()` when same variable was used for a fixed and random effect (e.g., if including a linear year trend and also an unstructured random intercept for year). Thanks to Liam Kendall for pointing this out.  

# spOccupancy 0.7.1

+ Small changes in C++ code to pass CRAN additional issues

# spOccupancy 0.7.0

spOccupancy v0.7.0 contains a variety of substantial updates, most notably functionality for fitting non-spatial and spatial multi-species multi-season occupancy models, as well as multi-species spatially-varying coefficient models. There are also a variety of smaller bug fixes/additional error handling that will help eliminate some common hard-to-interpret errors that users encountered.

+ New functionality for fitting multi-species, multi-season occupancy models. The function `tMsPGOcc()` fits non-spatial, multi-season, multi-species occupancy models, and the function `stMsPGOcc()` fits spatial, multi-season occupancy models. The spatially-explicit function also inherently accounts for species correlations with a spatial factor modeling approach (e.g., they are joint species distribution models with imperfect detection and species correlations). See [Doser et al. 2023](https://doi.org/10.1002/ecy.4137) for statistical details on the spatial factor modeling approach. A vignette will be posted that details fitting these models in depth in the coming months, but the syntax is essentially a combination of multi-species models (e.g. `msPGOcc()`, `sfMsPGOcc()`) and multi-season single-species models (i.e., `tPGOcc()` and `stPGOcc()`), so the recommendations provided in the vignettes for those models is applicable for these models as well. 
+ New functionality for multi-species spatially-varying coefficient occupancy models for single-season (`svcMsPGOcc()`) and multi-season (`svcTMsPGOcc()`) models. These approaches use a spatial factor modeling approach for each of the SVCs to make the models relatively computationally efficient. The functions inherently account for species correlations. The vignette on spatially-varying coefficients provides an example for `svcMsPGOcc()`, with an example for `svcTMsPGOcc()` coming soon.  
+ The function `simTMsOcc()` simulates multi-season, multi-species occupancy models.
+ Updated `getSVCSamples()` to now work with multi-species spatially-varying coefficient models.
+ Added in a new check in all spatially-explicit models to see if all the spatial coordinates in the `data$coords` object were unique, as this is a requirement for `spOccupancy` spatially-explicit models. In previous versions, this resulted in an error of `c++ error: dpotrf failed`, or something along those lines, which was a common source of confusion.
+ Updated all model fitting functions to avoid running for a long time, just to eventually crash. Now, if trying to run models and save an object that is too large for memory, R should crash at the beginning. This occurred in situations where the `n.burn` argument was greater than 0 and/or `n.thin` was greater than 1. Thanks to Alex Bacjz for bringing this to my attention.
+ Added in the `by.sp` argument to `waicOcc()` to allow users to calculate WAIC separately for individual species in all multi-species model types in `spOccupancy`.
+ Minor updates to multiple vignettes to reflect changes since their original versions.

# spOccupancy 0.6.1

+ Small changes in C++ code to pass CRAN additional issues

# spOccupancy 0.6.0

+ Incorporates new functionality to fit a non-spatial integrated multi-species occupancy model using the function `intMsPGOcc()`. This fits a single-season version of the "intgrated community occupancy model" from [Doser et al. 2022](https://doi.org/10.1111/2041-210X.13811). The function `intMsPGOcc()` should be considered [experimental](https://lifecycle.r-lib.org/articles/stages.html#experimental) and is still under development. We have done adequate testing of the function and users can be confident the resulting estimates are correct. Rather, we consider this "experimental" because it lacks all the functionality currently supported for other `spOccupancy` model types. In particular, `intMsPGOcc` model objects do not currently work with `ppcOcc()` (posterior predictive checks), `fitted()` (generated fitted values), or k-fold cross-validation, and there may be specific data set situations that cause the function to break. Please contact us if you use the function and have any feedback or run into any problems. We are in active development of the associated spatial versions of the function (both without spatial factors and with spatial factors), as well as the above mentioned limitations. `intMsPGOcc()` does not currently support random effects in the detection models, which we are actively working on. See `vignette("integratedMultispecies")` for more details.
+ New functionality to fit posthoc linear models to parameter estimates using the function `postHocLM()`. The function `postHocLM()` fits a basic linear (mixed) model to a response variable that is assumed to come from a previous model fit, and thus each value in the data response variable has a full set of posterior MCMC samples associated with it. While this function can be used for a variety of situations (including objects that don't come from `spOccupancy`), `postHocLM()` may be particularly useful for use with multi-species occupancy models to explore associations between species-specific covariate effect estimates from a multi-species occupancy model with species-level covariates, while fully accounting for uncertainty in the estimates. A vignette displaying how to do this will be posted in the coming months, but see the documentation for the function for basic instructions on how to use the function. 
+ Updated `sfMsPGOcc()` to allow a half-t prior on the community-level variance parameters. See documentation for more information on how to specify this. All multi-species occupancy model fitting functions will eventually be updated to allow for this prior, which can be a less informative prior when sample sizes (i.e., the number of species in this case) is low.
+ Updated `intPGOcc()` and `spIntPGOcc()` to remove an error that may occur if a data set only has site level detection covariates.
+ Updated `getSVCSamples()` to eliminate errors that prevented the function from working under certain circumstances depending on which covariates in the design matrix were modelled as spatially-varying coefficients.
+ Updated `tPGOcc()` and `stPGOcc()` to eliminate an error that occurred when trying to run these models with single-visit data sets.
+ Added in the `mis.spec.type` and `scale.param` arguments to the `simTOcc()` function to simulate multi-season detection-nondetection data under varying forms of model mis-specification. See `simTOcc()` documentation for detials. Thanks to Sara Stoudt for her help with this. 
+ Updated a typo in the MCMC sampler documentation for multi-species occupancy models. Specifically, the **T** in the mean component of Equations 23 and 24 from the [MCMC samplers vignette](https://doserlab.com/files/spoccupancy-web/articles/mcmcSamplers.pdf) was incorrect, and instead is now correctly **T**$^{-1}$. Similarly, Equations 9 and 10 were updated in the [MCMC samplers for factor models vignette](https://doserlab.com/files/spoccupancy-web/articles/mcmcFactorModels.pdf). Note that these were just typos in the vignettes, the underlying models are correct.
+ Fixed a bug that prevented `ppcOcc()` from working when there were only site-level random effects on detection. This also sometimes caused problems with cross-validation functionality as well. Thanks to Jose Luis Mena for bringing this to my attention. 

# spOccupancy 0.5.2

spOccupancy v0.5.2 contains an important bug fix in the cross-validation functionality for single-season occupancy models with unbalanced sampling across replicates in the data set. Specifically, the reported cross-validation deviance metrics may be inaccurate when one or more sites had a detection history where a missing value came before a non-missing value. For example, if one or more sites had a detection history of `c(NA, 1, 0, 0, 1)`, this would lead to the problem occurring, but this would not occur if all missing values were at the end of the detection history (e.g., `c(1, 0, 0, 1, NA)`). The affected functions include the following: `PGOcc()`, `spPGOcc()`, `msPGOcc()`, `spMsPGOcc()`, `lfMsPGOcc()`, `sfMsPGOcc()`, `intPGOcc()`, `spIntPGOcc()`. We strongly encourage users who have performed cross-validation with these models and unbalanced sampling across replicates in the manner described to rerun their analyses using v0.5.2. We apologize for any troubles this has caused.

# spOccupancy 0.5.1

+ Fixed issues with unicode text in the manual for passing CRAN checks on Windows
+ Fixed a bug in the k-fold cross-validation for models that include unstructured random intercepts on the occupancy portion of the model. This bug could have led to inacurrate cross-validation metrics when comparing a model with the unstructured random effect and without the unstructured random effect. We strongly encourage users who have performed cross-validation under such a scenario to rerun their analyses using v0.5.1. 

# spOccupancy 0.5.0

spOccupancy v0.5.0 contains numerous substantial updates that provide new functionality, improved run times for models with unstructured random effects, an important bug fix for cross-validation with unstructured random effects under certain scenarios, and some other minor bug fixes. The changes include: 

+ New functionality for fitting spatially-varying coefficient occupancy models. The function `svcPGOcc()` fits a single-season spatially-varying coefficient model, and `svcTPGOcc()` fits a multi-season spatially-varying coefficient model. We also include the functions `svcPGBinom()` and `svcTPGBinom()` for fitting spatially-varying coefficient generalized linear models when ignoring imperfect detection. We also include the helper function `getSVCSamples()` to more easily extract the SVC samples from the resulting model objects if they are desired.
+ Updated the underlying `C++` code to reduce run times for models that include unstructured random intercepts. 
+ Added the `k.fold.only` argument to all model-fitting functions, which allows users to only perform k-fold cross-validation instead of having to run the model first with the entire data set.
+ Adjusted how random intercepts in the detection model were being calculated, which resulted in unnecessary massive objects when fitting a model with a large number of random effect levels and spatial locations. See [GitHub issue 14](https://github.com/biodiverse/spOccupancy/issues/14). 
+ Fixed a bug that prevented prediction from working for multi-species models when `X.0` was supplied as a data frame and not a matrix. See [GitHub issue 13](https://github.com/biodiverse/spOccupancy/issues/13).
+ Fixed an error that occurred when the detection-nondetection data were specified in a specific way. See [GitHub issue 12](https://github.com/biodiverse/spOccupancy/issues/12).


# spOccupancy 0.4.0

+ Major new functionality for fitting multi-season (i.e., spatio-temporal) single-species occupancy models using the functions `tPGOcc()` and `stPGOcc()`. 
+ Fixed a bug in calculation of the detection probability values in `fitted()` functions for all spOccupancy model objects. See [this Github issue](https://github.com/biodiverse/spOccupancy/issues/10) for more details. 
+ Fixed an error that occurred when predicting for multi-species models and setting `ignore.RE = TRUE`.  
+ Fixed other small bugs that caused model fitting functions to break under specific circumstances.

# spOccupancy 0.3.2

+ Fixed a bug in `waicOcc()` for integrated models (`intPGOcc()` and `spIntPGOcc()`) that sometimes resulted in incorrect estimates of WAIC for data sets other than the first data set. We strongly encourage users who have used `waicOcc()` with an integrated model to rerun their analyses using v0.3.2. 
+ Fixed a bug introduced in v0.3.0 that sometimes resulted in incorrect predictions from a spatially-explicit model with non-spatial random effects in the occurrence portion of the model. We strongly encourage users who have used `predict()` on a spatially-explicit model with non-spatial random effects in the occurrence portion of the model to rerun their analyses using v0.3.2.
+ Users can now specify a uniform prior on the spatial variance parameter instead of an inverse-Gamma prior. We also allow users to fix the value of the spatial variance parameter at the initial value. See the reference pages of spatially-explicit functions for more details. 
+ Slight changes in the information printed when fitting spatially-explicit models. 
+ Removed dependency on spBayes to pass CRAN checks. 

# spOccupancy 0.3.1

+ Fixed two small problems with `intPGOcc()` and `spIntPGOcc()` that were accidentally introduced in v0.3.0. See [this Github issue](https://github.com/biodiverse/spOccupancy/issues/8) for more details.
+ Adapted C/C++ code to properly handle characters strings when calling Fortran BLAS/LAPACK routines following the new requirements for R 4.2.0.  

# spOccupancy 0.3.0

spOccupancy Version 0.3.0 contains numerous substantial updates that provide new functionality, improved computational performance for model fitting and subsequent model checking/comparison, and minor bug fixes. The changes include: 

+ Additional functionality for fitting spatial and non-spatial multi-species occupancy models with residual species correlations (i.e., joint species distribution models with imperfect detection). See documentation for `lfMsPGOcc()` and `sfMsPGOcc()`. We also included the functions `lfJSDM()` and `sfJSDM()` which are more typical joint species distribution models that fail to explicitly account for imperfect detection.
+ All single-species and multi-species models allow for unstructured random intercepts in both the occurrence and detection portions of the occupancy model. Prior to this version, random intercepts were not supported in the occurrence portion of spatially-explicit models. 
+ `predict()` functions for single-species and multi-species models now include the argument `type`, which allows for prediction of detection probability (`type = 'detection'`) at a set of covariate values as well as predictions of occurrence (`type = 'occupancy'`). 
+ All models are substantially faster than version 0.2.1. We improved performance by implementing a change in how we sample the latent Polya-Gamma variables in the detection component of the model. This results in substantial increases in speed for models where the number of replicates varies across sites. We additionally updated how non-spatial random effects were sampled, which also contributes to improved computational performance.
+ All model fitting functions now include the object `like.samples` in the resulting model object, which contains model likelihood values needed for calculation of WAIC. This leads to much shorter run times for `waicOcc()` compared to previous versions.
+ All `fitted.*()` functions now return both the fitted values and the estimated detection probability samples from a fitted `spOccupancy` model. 
+ Improved error handling for models with missing values and random effects.
+ Added the argument `ignore.RE` to all `predict()` functions. If non-spatial random intercepts are included when fitting the model, setting `ignore.RE = TRUE` will yield predictions that ignore the values of the random effects. If `ignore.RE = FALSE`, the model will predict new values using the random intercepts for both sampled and non-sampled levels of the effects.  
+ Fixed a bug in the cross-validation component of all `spOccupancy` model fitting functions that occurred when random effects were included in the occurrence and/or detection component of the model.
+ Fixed minor bug in `simOcc()` and `simMsOcc()` that prevented simulating data with multiple random intercepts on detection. 
+ Fixed minor bug in spatially-explicit models that resulted in an error when setting `NNGP = FALSE` and not specifying initial values for the spatial range parameter `phi`. 
+ Fixed a bug in the `predict()` functions for `spMsPGOcc` and `spPGOcc` objects that resulted in potentially inaccurate predictions when `n.omp.threads` > 1. 

# spOccupancy 0.2.1

+ Minor changes related to arguments in C++ code in header files to pass CRAN additional issues.

# spOccupancy 0.2.0

+ Added an `n.chains` argument to all model-fitting functions for running multiple chains in sequence.
+ Added posterior means, standard deviations, Gelman-Rubin diagnostic (Rhat) and Effective Sample Size (ESS) to `summary` displays for each model-fitting function.
+ Fixed spatially-explicit `predict` functions to return occurrence probabilities at sampled sites instead of NAs.

# spOccupancy 0.1.3

+ Minor bug fixes related to memory allocation in C++ code.

# spOccupancy 0.1.2

* This is the first release of `spOccupancy`.
