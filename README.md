<!-- Top banner -->
<p align="center">
  <img src="man/figures/designbgm_banner.svg" width="100%" alt="designbgm">
</p>

<p align="center">
  <em>Plan the sample size of a prospective study for a Bayesian graphical
  model, given a prior elicited from a previous study.</em>
</p>

<!-- badges: start -->
[![github-repo-status](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![R-package-version](https://img.shields.io/github/r-package/v/Bayesian-Graphical-Modelling-Lab/designbgm)](https://github.com/Bayesian-Graphical-Modelling-Lab/designbgm)
[![R-CMD-check](https://github.com/Bayesian-Graphical-Modelling-Lab/designbgm/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Bayesian-Graphical-Modelling-Lab/designbgm/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/Bayesian-Graphical-Modelling-Lab/designbgm/graph/badge.svg)](https://app.codecov.io/gh/Bayesian-Graphical-Modelling-Lab/designbgm)
<!-- badges: end -->

An informative prior elicited from a previous study carries a quantifiable amount of information, the prior effective sample size. A prospective study should be large enough that its data outweigh that prior. `designbgm` estimates the prior effective sample size for Bayesian graphical models and plans the prospective sample size accordingly, under a data-to-prior information ratio or a Bayes factor design analysis.

## Installation

```r
# install.packages("remotes")
remotes::install_github("Bayesian-Graphical-Modelling-Lab/designbgm")
```

## Usage

Start from the parameters of a previous study: a precision matrix `K`, its
graph `G`, and the study size `nu`.

```r
library(designbgm)

params <- ggm_parameters(K, G, nu = 100)
ep     <- elicit_prior(params)
```

How much information does this elicited prior carry?

```r
prior_ess(ep)
```

Plan the size of the prospective study. `"DPIR"` targets a data-to-prior information ratio over the
model parameters; `"BFDA"` targets the Bayes factor power at a representative edge:

```r
plan <- design(ep, method = "DPIR")
plan
```

Validate the recommendation:

```r
validate(plan)
```

For simulation studies, `simulate_prior_study()` generates prior studies from the family likelihood:

```r
st  <- ggm_study(p = 10, nu = 100, structure = "smallworld")
pst <- simulate_prior_study(st, n_studies = 1)
```

## Background

Prior effective sample size (ESS) is a pre-data measure of how much information an elicited prior carries. 
It matters because it sets a baseline to the planning: a study should collect more than *ESS(θ)* observations, or the prior
determines the posterior. The ESS estimators and the planning methods implemented in this package are described in Arena et al. (2026).

## Scope

Currently, the package supports Gaussian graphical models with a Wishart (complete graph) or G-Wishart (sparse graph) prior. The package is structured by family model, and ordinal Markov random fields are planned.

## References

- Arena, G., et al. (2026). *What is your Prior Worth? Effective Sample Size and Sample Size Planning for Gaussian Graphical Models* arXiv. 
[arXiv:2606.22687](https://doi.org/10.48550/arXiv.2606.22687)

