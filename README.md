# hhh4addon

This repository contains a development version of `hhh4addon`, an add-on package to the `surveillance` package. `hhh4addon` adds additional functionality to `surveillance:hhh4` which implements endemic-epidemic modelling for infectious disease count time series. Specifically it extends the model class towards distributed lags and adds functions for the analytical calculation of probabilistic path forecasts. Note that it is currently *only partly compatible* with the `hhh4contacts` package, another add-on package for `surveillance::hhh4` by Sebastian Meyer.

The exported functions are by now largely documented, including short examples. to get an overview over the functionality of the package please consider the vignette `vignette(hhh4addon)`.
