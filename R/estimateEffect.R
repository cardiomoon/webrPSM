utils::globalVariables(c("id", "subclass"))

#' Estimate Effect after matching
#' @param out An object of class "matchit"
#' @param mode One of c("continuous","binary","survival")
#' @param multiple logical Whether or not perform multiple regression
#' @param dep Name of dependent variable
#' @param covarCentering logical
#' @param withinSubclass logical
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
#' estimateEffect(out,mode="binary",dep=c("Y_B"),multiple=FALSE)
#' estimateEffect(out,mode="binary",dep=c("Y_B"))
#' estimateEffect(out,mode="survival",dep=c("Y_S"),multiple=FALSE)
#' estimateEffect(out,mode="survival",dep=c("Y_S"))
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
estimateEffect=function(out,mode="continuous",multiple=TRUE,dep,covarCentering=FALSE,withinSubclass=FALSE){

      # mode="survival";multiple=FALSE;dep="Y_C"
    xvars=attr(out$model$terms,"term.labels")
    yvar=names(out$model$model)[1]
    replace=FALSE
    if(is.null(out$info$replace)){
      md=match.data(out)
    } else if(out$info$replace){
        md=get_matches(out)
        replace=TRUE
    } else{
        md=match.data(out)
    }
    if(covarCentering){
        md[xvars]<-scale(md[xvars],scale=FALSE)
    }

    for(i in seq_along(dep)){

        if(mode=="survival"){
            temp=paste0("survival::Surv(",dep[i],")~",yvar)
        } else{
            temp=paste0(dep[i],"~",yvar)
        }
        if(covarCentering){
           temp=paste0(temp,"*(",paste0(xvars,collapse="+"),")")
        } else if(multiple){
            temp=paste0(temp,"+",paste0(xvars,collapse="+"))
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
            result1=cbind(res$coef,res$conf.int)[yvar,,drop=FALSE]
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
    rownames(result)=dep
    result
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
