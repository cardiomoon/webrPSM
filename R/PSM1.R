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
    SMD=summary(out)$sum.matched[,3]
    highSMD=which(abs(SMD)>0.1)
    highSMDNames=names(SMD)[highSMD]
    covarno=nrow(summary(out)$sum.matched)
    allcovarno=nrow(summary(out,interactions=TRUE)$sum.matched)
    SMD2=summary(out,interactions=TRUE)$sum.matched[(covarno+1):allcovarno,3]
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
                        " covariates",descNames(highSMDNames2)," were above 0.15, indicating inadequate balance.")
        }
    } else{
        temp=paste0(temp,
                    "After matching, standardized mean differences for the covariates ",
                    descNames(highSMDNames), " were above 0.1 indicating poor balance.")
    }
    highSMD
    highSMD2
    list(
        adequate=length(c(highSMD,highSMD2))==0,highSMD=highSMD,highSMD2=highSMD2,temp=temp
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
#' out1=matchit(formula, data = lalonde, method= "full",distance="glm",link="probit")
#' reportPSM(out1,depvar="re78")
reportPSM=function(out,depvar="",compare=TRUE){
     # depvar="re78"; compare=TRUE

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

    result$method

    temp=paste0("We used propensity score matching to estimate the average marginal effect of the '",yvar,"'")
    if(depvar!="") temp=paste0(temp," on '",depvar,"'")
    temp=paste0(temp," on those who received it accounting for confounding by the included covariates. ")
    resultSMD=reportSMD(out)
    resultSMD
    if((result$method!="nearest")&(compare)){
        temp=paste0(temp," We first attempted 1:1 nearest neighbor propensity score matching without replacement",
                    " with a propensity score estimated using logistic regression of the treatment on the covariates.")
        temp1=paste0("matchit(",yvar,"~",paste0(xvars,collapse="+"),",data=",dfname,")")
        out1=eval(parse(text=temp1))
        result1=reportSMD(out1)
        result1
        if(!result1$adequate) {
            temp=paste0(temp," This matching yielded poor balance, so we instead tried ",
                        result$method," matching on the propensity score, which yielded")

            if(resultSMD$adequate) {
                temp=paste0(temp," adequate balance, as indicated in Table and Figure." )
            } else{
                temp=paste0(temp," inadequate balance again, as indicated in Table and Figure." )
            }
        }
    } else{

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

    if(depvar!=""){
      fit1=lm(as.formula(paste0(depvar,"~",yvar,"+",paste0(xvars,collapse='+'))),data=match.data(out),weights=weights)
      x=lmtest::coeftest(fit1,vcov.=vcovCL,cluster=~subclass)[yvar,,drop=FALSE]
      result=reportEffect(fit1,yvar=yvar,depvar=depvar)
      temp=paste0(temp,result)
    }
    temp
}



#' Report treatment Effect
#' @param out An object of class lm
#' @param yvar The name of group variable
#' @param depvar The name of dependent variable
#' @param digits integer indicating the number of decimal places
#' @param sedigits digits for standard error
#' @param pdigits digits for p value
#' @param se logical If true, report se. If false report confidence interval
#' @export
#' @examples
#'library(MatchIt)
#'library(sandwich)
#' formula= treat ~ age + educ + race+ married +nodegree + re74 + re75
#' out=matchit(formula, data =lalonde, method= "full",link="probit")
#' formula1= re78~treat +age + educ + race+ married +nodegree + re74 + re75
#' out1=lm(formula1,data=match.data(out),weights=weights)
#' reportEffect(out1,depvar="re78")
#' reportEffect(out1,depvar="re78",se=FALSE)
reportEffect=function(out,yvar="treat",depvar="re78",digits=2,sedigits=2,pdigits=4,se=TRUE){
    # out=fit1
    # yvar="treat";depvar="re78"
    x=lmtest::coeftest(out,vcov.=vcovCL,cluster=~subclass)[yvar,,drop=FALSE]
    ci=lmtest::coefci(out,vcov.=vcovCL,cluster=~subclass)[yvar,,drop=FALSE]
    temp=paste0("To estimate the '",yvar,"' effect and its standard error, we fit a linear regression model with '",
                 depvar,"' as the outcome and the '",yvar,"' and the covariates as additive predictors and ",
                "included the matching weights in the estimation. The coefficient on the '",yvar,
                "' was taken to be the estimate of the '",yvar,"' effect. The lm() function was used to estimate the effect",
                ", and a cluster-robust variance as implemented in the vcovCL() function in the sandwich package was used",
                " to estimate its standard error with matching stratum membership as the clustering variable. ")

        if(se){
          temp=paste0(temp,"The estimated effect was ",round(x[1],digits),"(SE = ",round(x[2],sedigits),
                      ", p = ",round(x[4],pdigits),"), ")
        } else {
          temp=paste0(temp,"The estimated effect was ",round(x[1],digits),"(95% CI: ",round(ci[1],digits),
                      " - ",round(ci[2],digits),", p = ",round(x[4],pdigits),"), ")
        }


        temp=paste0(temp,"indicating that the average effect of the '",
                    yvar,"' for those who received it is ")
        if(x[4]<0.05){
           temp=paste0(temp,"to ",ifelse(x[1]>=0,"increase","decrease")," '",depvar,"'.")
        } else{
           temp=paste0(temp,"insignificant.")
        }
        temp

}

#'Make balance table
#'@param out An object of a class matchit
#'@importFrom cobalt bal.tab
#'@importFrom MatchIt matchit
#'@export
#'@examples
#'library(MatchIt)
#'formula=treat ~ age + educ + race+ married +nodegree + re74 + re75
#'out=matchit(formula, data =lalonde, method= "full",link="probit")
#'makeCompareBalTab(out)
makeCompareBalTab=function(out){
  out1=MatchIt::matchit(out$formula,data=out$model$data)
  weights=list(full=out,nn=out1)
  names(weights)[1]=out$info$method
  cobalt::bal.tab(out$formula,data=out$model$data,un=TRUE,weights=weights)
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
