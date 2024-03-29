% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimateEffectContinuous.R
\name{plotCompareEffects}
\alias{plotCompareEffects}
\title{Draw plot comparing IPW and linear regression}
\usage{
plotCompareEffects(
  data,
  dep = "y",
  xvars = NULL,
  treatvar = NULL,
  weights = "IPW",
  seed = 1234,
  se = TRUE,
  print = TRUE
)
}
\arguments{
\item{data}{A data.frame as a result of addIPW()}

\item{dep}{Name of dependent variable}

\item{xvars}{Names of covariates}

\item{treatvar}{Name of treatment variable}

\item{weights}{Name of weight variable}

\item{seed}{A single integer}

\item{se}{logical}

\item{print}{logical}
}
\description{
Draw plot comparing IPW and linear regression
}
\examples{
set.seed(1234)
mydata=addIPW(treat~x1+x2,data=simData2)
plotCompareEffects(mydata)
}
