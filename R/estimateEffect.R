utils::globalVariables(c("id", "subclass"))

#' Estimate Effect after matching
#' @param out An object of class "matchit"
#' @param mode One of c("continuous","binary","survival")
#' @param multiple logical Whether or not perform multiple regression
#' @param dep Name of dependent variable
#' @param time Name of time variable
#' @param status Name of status variable
#' @param covarCentering logical
#' @param withinSubclass logical
#' @param digits integer indicating the number of decimal places
#' @param sedigits digits for standard error
#' @param pdigits digits for p value
#' @param se logical If true, report se. If false report confidence interval
#' @param print logical
#' @importFrom survival Surv
#' @importFrom margins margins
#' @importFrom stats binomial glm quasibinomial
#' @importFrom lmtest coefci
#' @importFrom sandwich vcovHC
#' @importFrom survival coxph
#' @importFrom MatchIt get_matches
#' @export
#' @examples
#' library(MatchIt)
#' library(survival)
#' out <- matchit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9, data = exData)
#' estimateEffect(out,dep=c("Y_C"),multiple=FALSE)
#' estimateEffect(out,dep=c("Y_C"))
#' estimateEffect(out,dep=c("Y_C","Y_B"))
#' estimateEffect(out,mode="binary",dep=c("Y_B"),multiple=FALSE)
#' estimateEffect(out,mode="binary",dep=c("Y_B"))
#' estimateEffect(out,mode="survival",dep=c("Y_S"),multiple=FALSE)
#' estimateEffect(out,mode="survival",dep=c("Y_S"))
#' \dontrun{
#' out=matchit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9, data = exData,
#'    link = 'linear.logit', caliper = .1, ratio = 3, replace = TRUE)
#' estimateEffect(out,dep=c("Y_C"),multiple=FALSE)
#' estimateEffect(out,dep=c("Y_C"))
#' estimateEffect(out,mode="binary",dep=c("Y_B"),multiple=FALSE)
#' estimateEffect(out,mode="binary",dep=c("Y_B"))
#' estimateEffect(out,mode="survival",dep=c("Y_S"),multiple=FALSE)
#' estimateEffect(out,mode="survival",dep=c("Y_S"))
#' out=matchit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9, data = exData,
#' method = "full", estimand = "ATE")
#' estimateEffect(out,dep=c("Y_C"),multiple=FALSE)
#' estimateEffect(out,dep=c("Y_C"),covarCentering=TRUE)
#' estimateEffect(out,mode="binary",dep=c("Y_B"),multiple=FALSE)
#' estimateEffect(out,mode="survival",dep=c("Y_S"),multiple=FALSE)
#' out=matchit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9, data = exData,
#' method = "subclass", estimand = "ATT",subclass = 8)
#' estimateEffect(out,dep=c("Y_C"),multiple=FALSE)
#' estimateEffect(out,dep=c("Y_C"),multiple=FALSE,withinSubclass=TRUE)
#' estimateEffect(out,mode="binary",dep=c("Y_B"),multiple=FALSE)
#' estimateEffect(out,mode="binary",dep=c("Y_B"))
#' estimateEffect(out,mode="survival",dep=c("Y_S"),multiple=FALSE)
#' out=matchit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9, data = exData,
#' method = "full", estimand = "ATT")
#' estimateEffect(out,mode="binary",dep=c("Y_B"))
#' estimateEffect(out,mode="continuous",dep=c("Y_B"))  # for conditional effect
#' out=matchit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9, data = swData,
#' method = "full", estimand = "ATE",s.weights="SW")
#' estimateEffect(out,dep=c("Y_C"))
#' x=matchit(treat ~ age + educ + race + married+nodegree + re74 + re75, data =lalonde,
#'    method="exact")
#' estimateEffect(x,dep=c("re78"))
#' }
estimateEffect=function(out,mode="continuous",multiple=TRUE,dep="",time="",status="",
                        covarCentering=FALSE,withinSubclass=FALSE,
                        digits=2,sedigits=2,pdigits=4,se=TRUE,print=TRUE){

  # mode="continuous";multiple=FALSE;dep="Y_C"
  # mode="continuous";multiple=TRUE;dep="re78"
  # covarCentering=FALSE;withinSubclass=FALSE;multiple=FALSE
  # digits=2;sedigits=2;pdigits=4;se=TRUE;print=TRUE
  # mode="survival";dep="";time="time";status="cens"

  temp1=formula2vars(out$formula)
  xvars=temp1$xvars
  xvars
  yvar=temp1$yvar
  yvar
  replace=FALSE
  deselect=c()


  if(is.null(out$info$replace)){
    md=match.data(out)
  } else if(out$info$replace){
    md=get_matches(out)
    replace=TRUE
  } else{
    md=match.data(out)
  }

  myCentering=function(x){
    if(is.numeric(x)){
      scale(x,scale=FALSE)
    } else{
      x
    }
  }
  if((dep[1]=="")&(time!="")&(status!="")) mode="survival"

  if(covarCentering){
    md[xvars]<-lapply(md[xvars],myCentering)
  }

  for(i in seq_along(dep)){

    if(mode=="survival"){
      if(dep[i]==""){
          temp=paste0("survival::Surv(",time,",",status,")~",yvar)
      } else if("Surv" %in%class(md[[dep[i]]])){
          temp=paste0(dep[i],"~",yvar)
      } else{
          temp=paste0("survival::Surv(",dep[i],")~",yvar)
      }
    } else{
      temp=paste0(dep[i],"~",yvar)
    }
    temp
    if(covarCentering){
      temp=paste0(temp,"*(",paste0(xvars,collapse="+"),")")
    } else if(multiple){

      for(j in seq_along(xvars)){
          if((length(unique(md[[xvars[j]]]))==1)&(!is.numeric(md[[xvars[j]]]))){
             deselect=c(deselect,j)
          }
      }
      if(length(deselect)>0) {
        xvars1=xvars[-deselect]
      } else{
        xvars1=xvars
      }
      temp=paste0(temp,"+",paste0(xvars1,collapse="+"))
    }
    if(withinSubclass & out$info$method=="subclass"){
      temp=paste0(dep[i],"~ subclass + subclass:",yvar,"-1")
    }
    form1=as.formula(temp)
    #str(form1)
    if(mode=="continuous"){
      if(withinSubclass & out$info$method=="subclass"){
        fit <- lm(form1,data = md)
        if(call2param(out$call)$estimand=="ATT"){
          result3=summary(margins::margins(fit,variables=yvar,
                                           data=md[md[[yvar]]==1,],
                                           vcov=vcovHC(fit)))
        } else{  #ATE
          result3=summary(margins::margins(fit,variables=yvar,
                                           vcov=vcovHC(fit)))
        }
      } else {
        fit <- lm(form1,data = md, weights = weights)
        if(out$info$method=="subclass"){
          result1=coeftest(fit,vcov.=vcovHC)[yvar,,drop=FALSE]
          result2=coefci(fit,vcov.=vcovHC)[yvar,,drop=FALSE]
        } else if(replace) {
          result1=coeftest(fit,vcov.=vcovCL,cluster=~subclass+id)[yvar,,drop=FALSE]
          result2=coefci(fit,vcov.=vcovCL,cluster=~subclass+id)[yvar,,drop=FALSE]
        } else{
          result1=coeftest(fit,vcov.=vcovCL,cluster=~subclass)[yvar,,drop=FALSE]
          result2=coefci(fit,vcov.=vcovCL,cluster=~subclass)[yvar,,drop=FALSE]
        }
        result3=cbind(result1,result2)
        result3=as.data.frame(result3)
      }

    } else if(mode=="binary"){
      if(replace | out$info$method %in% c("full","subclass")) {
        fit <- glm(form1,data = md, family=quasibinomial(link="logit"),weights = weights)
      } else{
        fit <- glm(form1,data = md, family=binomial(link="logit"),weights = weights)
      }
      if(out$info$method =="subclass"){
        result1=coeftest(fit,vcov.=vcovHC,cluster=~subclass)[yvar,,drop=FALSE]
        result2=coefci(fit,vcov.=vcovHC,cluster=~subclass)[yvar,,drop=FALSE]
      } else if(replace){
        result1=coeftest(fit,vcov.=vcovCL,cluster=~subclass+id)[yvar,,drop=FALSE]
        result2=coefci(fit,vcov.=vcovCL,cluster=~subclass+id)[yvar,,drop=FALSE]
      } else{
        result1=coeftest(fit,vcov.=vcovCL,cluster=~subclass)[yvar,,drop=FALSE]
        result2=coefci(fit,vcov.=vcovCL,cluster=~subclass)[yvar,,drop=FALSE]
      }
      result3=cbind(result1,result2)
      result3=as.data.frame(result3)
      result3$OR=exp(result3[[1]])
      result3$lower=exp(result3[[5]])
      result3$upper=exp(result3[[6]])
      result3=result3[c(1:4,7:9)]
      result3

    } else if(mode=="survival"){
      if(replace) {
        fit=austinCafri(form1,data=md)
      } else if(out$info$method=="subclass"){
        fit=coxph(form1,data=md,robust=TRUE,weights=weights)
      } else {
        fit=coxph(form1,data=md,robust=TRUE,cluster=subclass,weights=weights)
      }
      res=summary(fit)
      result1=cbind(res$coef,res$conf.int)[1,,drop=FALSE]
      result3=as.data.frame(result1)[c(1,3:7,9,10)]
      names(result3)[6]="HR"

    }
    if(i==1) {
      result=result3
    } else{
      result=rbind(result,result3)
    }

  }
  result
  if("summary.margins" %in% class(result)) {
    result=result[-1]
    class(result)="data.frame"
  }
  dep
  if(dep[1]!="") rownames(result)=dep
  report=""
  if(out$info$method=="subclass"){
    report="To estimate marginal effects after stratum matching, marginal mean weighting through stratification(MMWS) was used. "
  }
  report=paste0(report,"To estimate the '",yvar,"' effect and its standard error")
  if(!is.null(out$s.weights)){
    report=paste0(report," in the matched and sampling weighted sample")
  }
  report=paste0(report,", we fit ")
  if(mode=="continuous") {
    report=paste0(report,"a linear regression model ")
  } else if(mode=="binary"){
    report=paste0(report,"a weighted generalized linear regression model ")
  } else{  # "survival"
    report=paste0(report,"a Cox proprotional harzards model ")
  }
  if(dep[1]!=""){
  report=paste0(report,"with ",descNames(dep)," as the outcome and the '",yvar,"' ")
  } else{
  report=paste0(report,"with a survival object as the outcome and the '",yvar,"' ")
  }
  if(multiple) {
    report=paste0(report,"and the covariates ")
    if(length(deselect)>0) {
      report=paste0(report,descNames(xvars[-deselect])," ")
    } else{
      report=paste0(report,descNames(xvars)," ")
    }
    report=paste0(report,"as additive predictors ")
  } else{
    report=paste0(report,"as a predictor ")
  }
  report=paste0(report,"and included the matching weights in the estimation. ",
                "The coefficient on the '",yvar,"' was taken to be the estimate of the '",yvar,"' effect. ")
  if(mode=="continuous") {
    report=paste0(report,"The lm() function was used to estimate the effect, ")
  } else if(mode=="binary") {
    report=paste0(report,"The glm() function with a logit link function was used to estimate the marginal odds ratio, ")
  } else{
    report=paste0(report,"The coxph() function in the survival package was used to estimate the marginal harzard ratio, ")
  }
  if(out$info$method=="subclass"){
    if(mode=="survival"){
      report=paste0(report,
                    "and a regular robust standard error was used because of the MMWS weights. ")
    } else{
      report=paste0(report,
                    "and a regular robust standard error as implemented in the vcovHC() function in the sandwich package was used",
                    " to estimate its standard error. ")
    }
    if(withinSubclass){
      report=paste0(report,
                    "The within-subclass effects was estimated as the result of marginal effects procedure as "
                    ,"implemented in the margins() function in the marigins package. ")
    }

  } else{
    if(mode=="survival" & !is.null(out$info$replace)){
      report=paste0(report,"and the Austin and Cafri's(2020) SE estimator was used to estimate robust standard error. ")
    } else{
      report=paste0(report,
                    "and a cluster-robust variance as implemented in the vcovCL() function in the sandwich package was used",
                    " to estimate its standard error with matching stratum membership as the clustering variable. ")
    }
  }
  for(i in 1:nrow(result)){
    x=result[i,]
    if(mode=="continuous"){
      if(se){
        report=paste0(report,"For the outcome '",dep[i],"', the estimated effect was ",round(x[1],digits),"(SE = ",round(x[2],sedigits))
      } else {
        report=paste0(report,"For the outcome '",dep[i],"',The estimated effect was ",round(x[1],digits),"(95% CI: ",round(x[5],digits),
                      " - ",round(x[6],digits))
      }
    } else if(mode=="binary"){
      report=paste0(report,"The estimated odds ratio was ",round(x[5],digits),"(95% CI: ",round(x[6],digits),
                    " - ",round(x[7],digits))
    } else{
      report=paste0(report,"The estimated harzards ratio was ",round(x[6],digits),"(95% CI: ",round(x[7],digits),
                    " - ",round(x[8],digits))
    }
    pno=4
    if(mode=="survival") pno=5
    if(x[pno]<0.0001){
      report=paste0(report,", p < .0001), ")
    }  else{
      report=paste0(report,", p = ",sprintf(paste0("%0.",pdigits,"f"),x[pno]),"), ")
    }


    report=paste0(report,"indicating that the average effect of the '",
                  yvar,"' for those who received it is ")
    if(x[pno]<0.05){
      if(dep[1]!=""){
      report=paste0(report,"to ",ifelse(x[1]>=0,"increase","decrease")," '",dep[1],"'.")
      } else{
        report=paste0(report,"to ",ifelse(x[1]>=0,"increase","decrease")," event.")
      }
    } else{
      report=paste0(report,"insignificant.")
    }
  }
  attr(result,"report")=report
  if(print) {
    cat("Estimate Effect\n\n")
    print(result)
    cat("\nInterpretation\n\n",report)
    invisible(result)
  } else{
    result
  }

}

#' Austin and Cafri's SE estimator
#' @param formula A formula
#' @param data A data.frame which is the result of get_mactchs() function
#' @importFrom survival coxph
#' @export
#' @source Austin, Peter C., and Guy Cafri. 2020.
#' “Variance Estimation When Using Propensity-Score Matching with Replacement with Survival or Time-to-Event Outcomes.”
#' Statistics in Medicine 39 (11): 1623–40. \url{https://doi.org/10.1002/sim.8502}.
#' @examples
#' library(MatchIt)
#' library(survival)
#' out=matchit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9, data = exData,
#'     link = 'linear.logit', caliper = .1,ratio = 3, replace = TRUE)
#' austinCafri(Surv(Y_S) ~ A,data=get_matches(out))
austinCafri=function(formula, data){
  fs <- coxph(formula, data = data, robust = TRUE,
              weights = weights, cluster = subclass)
  Vs <- fs$var
  ks <- nlevels(data$subclass)

  fi <- coxph(formula, data = data, robust = TRUE,
              weights = weights, cluster = id)
  Vi <- fi$var
  ki <- length(unique(data$id))

  fc <- coxph(formula, data = data, robust = TRUE,
              weights = weights)
  Vc <- fc$var
  kc <- nrow(data)

  #Compute the variance
  V <- (ks/(ks-1))*Vs + (ki/(ki-1))*Vi - (kc/(kc-1))*Vc

  #Sneak it back into the fit object
  fc$var <- V

  fc
}
