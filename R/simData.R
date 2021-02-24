#' simulation data by Young Min Baek
#'
#' A toy dataset containing numeric variables(V1-V3), treatment variable(treat), Reverse code of treat(Rtreat) and continuous outcome(y).
#'
#'@format A data frame with 1000 rows and 6 variables
#'\describe{
#'  \item{V1}{covariate 1}
#'  \item{V2}{covariate 2}
#'  \item{V3}{covariate 3}
#'  \item{treat}{treatment variable, 0 or 1}
#'  \item{Rtreat}{reverse code of treat. simData$Rtreat=!simData$treat}
#'  \item{y}{continuous outcome}
#'}
#'@source \url{https://sites.google.com/site/ymbaek}
"simData"
