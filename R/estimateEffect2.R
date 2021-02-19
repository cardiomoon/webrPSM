#' Estimate Effect after matching using the bootstrap
#' @param out An object of class "matchit"
#' @param mode One of c("continuous","binary","survival")
#' @param multiple logical Whether or not perform multiple regression
#' @param dep Name of dependent variable
#' @param covarCentering logical
#' @param withinSubclass logical
#' @param seed numeric
#' @param digits integer indicating the number of decimal places
#' @param sedigits digits for standard error
#' @param pdigits digits for p value
#' @param se logical If true, report se. If false report confidence interval
#' @param print logical
#' @importFrom survival Surv
#' @importFrom boot boot boot.ci
#' @importFrom MatchIt get_matches
#' @importFrom stats coef predict
#' @export
#' @examples
#' library(MatchIt)
#' library(survival)
#' out <- matchit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9, data = exData)
#' estimateEffect2(out,dep="Y_C",multiple=FALSE)
#' \dontrun{
#' estimateEffect2(out,mode="binary",dep="Y_B",multiple=FALSE)
#' estimateEffect2(out,mode="binary",dep="Y_B")
#' estimateEffect2(out,mode="survival",dep="Y_S",multiple=FALSE)
#' out <- matchit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9, data = exData,
#'             link = 'linear.logit', caliper = .1, ratio = 3, replace = TRUE)
#' estimateEffect2(out,dep="Y_C",multiple=FALSE)
#' estimateEffect2(out,dep="Y_B",mode="binary",multiple=TRUE)
#' }
estimateEffect2=function(out,mode="continuous",multiple=TRUE,dep,
                         covarCentering=FALSE,withinSubclass=FALSE,seed=1,
                         digits=2,sedigits=2,pdigits=4,se=TRUE,print=TRUE){
# mode="continuous";multiple=FALSE;dep="Y_C"
# mode="binary";multiple=TRUE;dep="Y_B"
# mode="survival";multiple=FALSE;dep="Y_S"
# covarCentering=FALSE;withinSubclass=TRUE
# digits=2;sedigits=2;pdigits=4;se=TRUE

xvars=attr(out$model$terms,"term.labels")
yvar=names(out$model$model)[1]
yvar
replace=FALSE

set.seed(seed)
if(is.null(out$info$replace)){
    md=match.data(out)
} else if(out$info$replace){
    md=get_matches(out)
    replace=TRUE
} else{
    md=match.data(out)
}


changeDataCall=function(out,data){
    paramList=call2param(out$call)
    paramList$replace=as.logical(paramList$replace)
    paramList$data=data
    paramList$formula=as.formula(paramList$formula)
    do.call(matchit,paramList)
}

est_fun2=function(data,i,multiple=FALSE,mode="continuous"){
     md2=changeDataCall(out,data=data[i,])
     md_boot=match.data(md2)

     if(mode=="survival"){
         temp=paste0("Surv(",dep,")~",yvar)
     } else{
         temp=paste0(dep,"~",yvar)
     }

     if(multiple){
         if(length(xvars)>0)temp=paste0(temp,"+",paste0(xvars,collapse="+"))
     }

     if(mode=="continuous"){
         fit_boot<-lm(as.formula(temp),data=md_boot,weights=weights)
         return(coef(fit_boot)[yvar])
     } else if(mode=="binary"){
         fit_boot<-glm(as.formula(temp),data=md_boot,
                       family=quasibinomial(link="logit"),
                       weights=weights)
         #Estimate potential outcomes for each unit
         md_boot[[yvar]] <- 0
         P0 <- weighted.mean(predict(fit_boot, md_boot, type = "response"),
                             w = md_boot$weights)
         Odds0 <- P0 / (1 - P0)

         md_boot[[yvar]] <- 1
         P1 <- weighted.mean(predict(fit_boot, md_boot, type = "response"),
                             w = md_boot$weights)
         Odds1 <- P1 / (1 - P1)

         #Return marginal odds ratio
         return(Odds1 / Odds0)

     }


}

pair_ids=levels(md$subclass)

est_fun<-function(pairs,i,multiple=FALSE,mode="continuous"){
     numreps <- table(pairs[i])
     ids<-unlist(lapply(pair_ids[pair_ids %in% names(numreps)],
                        function(p) rep(which(md$subclass==p),numreps[p])))

     md_boot<-md[ids,]

     if(mode=="survival"){
         temp=paste0("Surv(",dep,")~",yvar)
     } else{
         temp=paste0(dep,"~",yvar)
     }

     if(multiple){
         if(length(xvars)>0)temp=paste0(temp,"+",paste0(xvars,collapse="+"))
     }

     if(mode=="continuous"){
     fit_boot<-lm(as.formula(temp),data=md_boot,weights=weights)
     return(coef(fit_boot)[yvar])
     } else if(mode=="binary"){
         fit_boot<-lm(as.formula(temp),data=md_boot,
                      family=binomial(link="logit"),
                      weights=weights)
         #Estimate potential outcomes for each unit
         #Under control
         md_boot[[yvar]] <- 0
         P0 <- weighted.mean(predict(fit_boot, md_boot, type = "response"),
                             w = md_boot$weights)
         Odds0 <- P0 / (1 - P0)

         #Under treatment
         md_boot[[yvar]] <- 1
         P1 <- weighted.mean(predict(fit_boot, md_boot, type = "response"),
                             w = md_boot$weights)
         Odds1 <- P1 / (1 - P1)

         #Return marginal odds ratio
         return(Odds1 / Odds0)
     } else if(mode=="survival"){
         cox_fit_boot <- coxph(as.formula(temp), data = md_boot)

         #Compute the marginal HR by exponentiating the coefficient
         #on treatment
         HR <- exp(coef(cox_fit_boot)[yvar])

         #Return the HR
         return(HR)
     }
}

replace
if(replace){
    boot_est<-boot(out$model$data,est_fun2,R=499,multiple=multiple,mode=mode)
} else{
   boot_est<-boot(pair_ids,est_fun,R=499,multiple=multiple,mode=mode)
}
boot_est
if(replace) {
    ci=boot.ci(boot_est,type="perc")
} else{
    ci=boot.ci(boot_est,type="bca")
}
ci
list(est=boot_est,ci=ci)
}
