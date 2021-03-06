% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rsaeGWG.R
\docType{package}
\name{rsaeGWR}
\alias{rsaeGWR}
\alias{rsaeGWR-package}
\title{rsaeGWR: A package for robust small area estimation under spatial non-stationarity}
\description{
The package provides two categories of important functions:
fit and predict.
}
\section{The fit functions}{

The fit functions \code{\link{gwlmm}} and \code{\link{rgwlmm}} fit a geographically weighted
linear mixed model (GWLMM) to data. The GWLMM is a nested error regression model that takes
into account spatial non-stationarity. The function \code{\link{gwlmm}} assumes normality
for the error term components. This assumption can be violated in the presence
of outliers. The function \code{\link{rgwlmm}} fits an outlier robust version of the GWLMM to the data.
}

\section{The predrict function}{

The predict functions \code{\link{predict.gwlmm}} and \code{\link{predict.rgwlmm}} estimate
small area means based on the model fits from the fit funtions when population data is provided.
As predcision measures MSE estiamtes are provided.
The predict functions can handle aggregated population information with centroid coordinates and
 unit-level population data with coordinates for each unit.
If no data is provided these functions provide the in-sample predictions
(random effects, residuals, etc.) based on the model fits.
}

