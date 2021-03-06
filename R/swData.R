#'Sample Data with Sampling Weight
#'
#'A toy dataset containing treatment variable(A), numeric variables(X1-X9), continuous outcome(Y_C), binary outcome(Y_B) and survival outcome(Y_S).
#'
#'@format A data frame with 2000 rows and 13 variables
#'\describe{
#'  \item{A}{treatment variable, 0 or 1}
#'  \item{X1}{covariate 1}
#'  \item{X2}{covariate 2}
#'  \item{X3}{covariate 3}
#'  \item{X4}{covariate 4}
#'  \item{X5}{covariate 5}
#'  \item{X6}{covariate 6}
#'  \item{X7}{covariate 7}
#'  \item{X8}{covariate 8}
#'  \item{X9}{covariate 9}
#'  \item{Y_C}{continuous outcome}
#'  \item{SW}{sampling weight}
#'}
#'@source \url{https://cran.r-project.org/web/packages/MatchIt/vignettes/sampling-weights.html}
"swData"
