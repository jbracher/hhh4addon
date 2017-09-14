# hhh4addon

This repository contains a development version of `hhh4addon`, an add-on package to the `surveillance` package. `hhh4addon` adds functionality to `surveillance:hhh4` which implements endemic-epidemic modelling for infectious disease count time series. Specifically it extends the model class towards distributed lags and adds functions for the analytical calculation of probabilistic path forecasts and cyclostationary moments. Note that it is currently *only partly compatible* with the `hhh4contacts` package, another add-on package for `surveillance::hhh4` by Sebastian Meyer.

The exported functions are by now largely documented, including short examples. to get an overview over the functionality of the package please consider the vignette `vignette(hhh4addon)`.

The theory behind longterm predictions and cyclostationary moments can be found in Held, Meyer and Bracher (2017): [Probabilistic forecasting in infectious disease epidemiology: the 13th Armitage lecture](http://onlinelibrary.wiley.com/doi/10.1002/sim.7363/full#references)  and Bracher and Held (2017): [Periodically stationary multivariate non-Gaussian autoregressive models](https://arxiv.org/abs/1707.04635).

To install the package directly from github type:

`library(devtools)`

`install_github("jbracher/hhh4addon")`