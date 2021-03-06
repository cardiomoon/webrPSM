% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimateEffectContinuous.R
\name{makePPTList_IPW}
\alias{makePPTList_IPW}
\title{make pptList with Matching with IPW}
\usage{
makePPTList_IPW(x, dep = "", seed = 1234)
}
\arguments{
\item{x}{Character}

\item{dep}{Name of dependent variable}

\item{seed}{Numeric}
}
\description{
make pptList with Matching with IPW
}
\examples{
result=makePPTList_IPW(x="addIPW(treat~x1+x2,data=simData2)",dep="y",seed=1234)
result=makePPTList_IPW(x="addIPW(educ~age+race+married,data=lalonde)",dep="re78",seed=1234)
}
