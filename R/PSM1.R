#'Describe names as natural language
#'@param x A charcter vector
#'@export
#'@examples
#'descNames("A")
#'descNames(c("A","B"))
#'descNames(LETTERS[1:5])
descNames=function(x=LETTERS[1:5]){
    no=length(x)
    if(no==0) {
        "None"
    } else if(no==1) {
        x
    } else{
        temp=paste0(x[-no],collapse=", ")
        temp=paste0(temp," and ",x[no])
        temp
    }
}

#' Report adequacy of balance
#' @param out An object of class matchit
#' @export
#' @examples
#'library(MatchIt)
#' formula=treat ~ age + educ + race+ married +nodegree + re74 + re75
#' out=matchit(formula, data =lalonde, method= "nearest")
#' reportSMD(out)
#' out1=matchit(formula, data =lalonde, method= "full",link="probit")
#' reportSMD(out1)
reportSMD=function(out){
    temp=""
    if(is.null(summary(out)$sum.matched)){
      matchSum=summary(out)$sum.across
      matchAllSum=summary(out,interactions=TRUE)$sum.across

    } else{
      matchSum=summary(out)$sum.matched
      matchAllSum=summary(out,interactions=TRUE)$sum.matched
    }

    SMD=matchSum[,3]

    highSMD=which(abs(SMD)>0.1)
    highSMDNames=names(SMD)[highSMD]
    covarno=nrow(matchSum)
    allcovarno=nrow(matchAllSum)
    SMD2=matchAllSum[(covarno+1):allcovarno,3]
    highSMD2=which(abs(SMD2)>0.15)
    highSMDNames2=names(SMD2)[highSMD2]
    if(length(highSMD)==0) {
        temp=paste0(temp,
                    " After matching, all standardized mean differences for the covariates were below 0.1")
        if(length(highSMD2)==0){
            temp=paste0(temp,
                        " and all standardized mean differences for squares and two-way interactions between",
                        " covariates were below 0.15, indicating adequate balance.")
        } else{
            temp=paste0(temp,
                        " but some standardized mean differences for squares and two-way interactions between",
                        " covariates ",descNames(highSMDNames2)," were above 0.15.")
        }
    } else{
        temp=paste0(temp,
                    "After matching, standardized mean differences for the covariates ",
                    descNames(highSMDNames), " were above 0.1 indicating poor balance.")
    }
    highSMD
    highSMD2
    list(
        #adequate=length(c(highSMD,highSMD2))==0,highSMD=highSMD,highSMD2=highSMD2,temp=temp
      adequate=length(c(highSMD))==0,highSMD=highSMD,highSMD2=highSMD2,temp=temp
    )
}

#' Report propensity sore matching
#' @param out An object of class matchit
#' @param depvar Character name of dependent variable
#' @param compare logical Whether or not compare to nearest model
#' @importFrom stats weights
#' @importFrom lmtest coeftest
#' @importFrom sandwich vcovCL
#' @importFrom stats lm as.formula
#' @export
#' @examples
#'library(MatchIt)
#' formula=treat ~ age + educ + race+ married +nodegree + re74 + re75
#' out=matchit(formula, data =lalonde, method= "nearest")
#' reportPSM(out)
#' out=matchit(formula, data = lalonde, method= "full",distance="glm",link="probit")
#' reportPSM(out,depvar="re78")
#' reportPSM(out,depvar="re78",compare=FALSE)
reportPSM=function(out,depvar="",compare=NULL){
        # depvar="re78"; compare=NULL

    xvars=attr(out$model$terms,"term.labels")
    yvar=names(out$model$model)[1]
    xvars
    length(xvars)
    yvar
    result=call2param(out$call)
    result
    temp=as.character(out$call)
    dfname=temp[3]
    dfname

    out$estimand
    result$method

    temp=paste0("We used propensity score matching to estimate the average marginal effect of the '",yvar,"'")
    if(depvar!="") temp=paste0(temp," on '",depvar,"'")
    if(out$estimand=="ATT"){
       temp=paste0(temp," on those who received it ")
    } else if(out$estimand=="ATE"){
      temp=paste0(temp," for all units in the target population ")
    } else if(out$estimand=="ATM"){
      temp=paste0(temp," in the remaining matched sample ")
    }
    temp=paste0(temp,"accounting for confounding by the included covariates. ")
    resultSMD=reportSMD(out)
    resultSMD
    if(is.null(compare)){
        if(result$method!="nearest") {
          compare=TRUE
        } else if(!is.null(out$info$replace)){
            if(out$info$replace) {
              compare=TRUE
            } else{
              compare=FALSE
            }
        } else if(!is.null(out$info$ratio)){
           compare=TRUE
        } else{
           compare=FALSE
        }
    }
    if(compare){
        temp=paste0(temp," We first attempted 1:1 nearest neighbor propensity score matching without replacement",
                    " with a propensity score estimated using logistic regression of the treatment on the covariates.")
        temp1=paste0("matchit(",yvar,"~",paste0(xvars,collapse="+"),",data=",dfname,")")
        out1=eval(parse(text=temp1))
        result1=reportSMD(out1)
        result1
        if(!result1$adequate) {
            temp=paste0(temp," But this matching yielded poor balance, so we tried ")
        }
    } else{
        temp=paste0(temp," We tried ")
    }
    if(!is.null(out$info$ratio)){
      temp=paste0(temp,out$info$ratio,":1 ")
    }
    temp=paste0(temp,result$method," matching on the propensity score ")
    if(!is.null(out$info$replace)){
        if(out$info$replace) {
          temp=paste0(temp,"with replacement")
        } else{
          temp=paste0(temp,"without replacement")
        }
    }
    temp=paste0(temp,", which yielded")

    if(resultSMD$adequate) {
                temp=paste0(temp," adequate balance, as indicated in Table and Figure. " )
            } else{
                temp=paste0(temp," inadequate balance, as indicated in Table and Figure. " )
            }



    temp=paste0(temp,resultSMD$temp)
    temp=paste0(temp," The propensity score was estimated using a ",
                ifelse(result$link=="logit","logistic",result$link),
                " regression of the '",yvar,"' on the covariates")
    if(result$link!="logit"){
        temp=paste0(temp,", which yielded better balance than did a logistic regression. ")
    } else{
        temp=paste0(temp,". ")
    }

    ## Matching method
    if(result$method=="full"){
        temp=paste0(temp," Full matching uses all treated and all control units, so no units were discarded by the matching. ")
    }

    temp
}



#'Make balance table
#'@param out An object of a class matchit
#'@param print logical
#'@importFrom cobalt bal.tab
#'@importFrom MatchIt matchit
#'@export
#'@examples
#'library(MatchIt)
#'formula=treat ~ age + educ + race+ married +nodegree + re74 + re75
#'out=matchit(formula, data =lalonde, method= "full",link="probit")
#'makeCompareBalTab(out)
makeCompareBalTab=function(out,print=TRUE){
  out1=MatchIt::matchit(out$formula,data=out$model$data)
  weights=list(full=out,nn=out1)
  names(weights)[1]=out$info$method
  res=cobalt::bal.tab(out$formula,data=out$model$data,un=TRUE,weights=weights)
  if(print){
  cat("Balance Meausures\n")
  cat("-----------------\n")
  print(res$Balance[c(1,2,4,6)])
  cat("\nEffective sample sizes\n")
  cat("----------------------\n")
  print(res$Observations)
  }
  invisible(res)
}

#'Drow love plot comparing to nearest match
#'@param out An object of a class matchit
#'@param stats character; which statistic(s) should be reported
#'@importFrom cobalt love.plot
#'@importFrom MatchIt matchit
#'@importFrom stringr str_to_title
#'@export
#'@examples
#'data(lalonde,package="MatchIt")
#'formula=treat ~ age + educ + race+ married +nodegree + re74 + re75
#'out=MatchIt::matchit(formula, data =lalonde, method= "full",link="probit")
#'compareLove.plot(out)
#'compareLove.plot(out,stats=c("m","ks"))
compareLove.plot=function(out,stats=c("m")){

  out1=MatchIt::matchit(out$formula,data=out$model$data)
  sample.names = c("Full Matching", "NN Matching", "Original")
  cobalt::love.plot(out, stats = stats, poly = 2, abs = TRUE,
            weights = list(nn = out1),
            drop.distance = TRUE, thresholds = c(m = .1),
            var.order = "unadjusted", binary = "std",
            shapes = c("triangle", "square", "circle"),
            colors = c("blue", "darkgreen", "red"),
            sample.names = c(paste0(str_to_title(out$info$method)," Matching"), "NN Matching", "Original"),
            position = "bottom")
}

# bal.plot(out,var.name="distance",which="both",type="histogram",mirror=TRUE)
