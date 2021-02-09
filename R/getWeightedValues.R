#' Compute Weight Variance
#' @param x A numeric vector
#' @param weights a numeric vector of weights
#' @param normwt logical specify normwt=TRUE to make weights sum to length(x) after deletion of NAs. If weights are frequency weights, then normwt should be FALSE, and if weights are normalization (aka reliability) weights, then normwt should be TRUE. In the case of the former, no check is made that weights are valid frequencies.
#' @param na.rm logical set to FALSE to suppress checking for NAs
#' @param method determines the estimator type; if 'unbiased' (the default) then the usual unbiased estimate (using Bessel's correction) is returned, if 'ML' then it is the maximum likelihood estimate for a Gaussian distribution. In the case of the latter, the normwt argument has no effect. Uses stats:cov.wt for both methods.
#' @export
#' @examples
#' set.seed(1)
#' x <- runif(500)
#' wts <- sample(1:6, 500, TRUE)
#' weighted.var(x,wts)
weighted.var=function (x, weights = NULL, normwt = FALSE, na.rm = TRUE,
                       method = c("unbiased","ML"))
{
    method <- match.arg(method)
    if (!length(weights)) {
        if (na.rm)
            x <- x[!is.na(x)]
        return(var(x))
    }
    if (na.rm) {
        s <- !is.na(x + weights)
        x <- x[s]
        weights <- weights[s]
    }
    if (normwt)
        weights <- weights * length(x)/sum(weights)
    if (normwt || method == "ML")
        return(as.numeric(stats::cov.wt(cbind(x), weights, method = method)$cov))
    sw <- sum(weights)
    if (sw <= 1)
        warning("only one effective observation; variance estimate undefined")
    xbar <- sum(weights * x)/sw
    sum(weights * ((x - xbar)^2))/(sw - 1)
}

#' Get weighted value or p value
#' @param x An object of class matchit
#' @param xvar character colume name
#' @param treat numeric  If 0, control. If 1, treat. If 2 p value.
#' @param digits numeric
#' @importFrom survey svydesign svychisq
#' @importFrom stats lm weighted.mean var
getWeightedValues=function(x,xvar,treat=1,digits=1){
    yvar=names(x$model$model)[1]
    df1=match.data(x)
    summatch=as.data.frame(summary(x)$sum.matched[-1,])
    yvalues=sort(unique(df1[[yvar]]))
    data1=df1[df1[[yvar]]==yvalues[1],]
    data2=df1[df1[[yvar]]==yvalues[2],]
    ncontrol=nrow(data1)
    ntreat=nrow(data2)
    temp=df1[[xvar]]
    form1=paste0("%.",digits,"f")
    form2=paste0("%",digits+2,".",digits,"f")
    if(is.numeric(temp) & length(unique(temp))>2) {
        if(treat==0){
            result=paste0(sprintf(form1,weighted.mean(data1[[xvar]],data1$weights))," \u00b1 ",
                          sprintf(form1,sqrt(weighted.var(data1[[xvar]],data1$weights))))
        } else if(treat==1){
            result=paste0(sprintf(form1,weighted.mean(data2[[xvar]],data2$weights))," \u00b1 ",
                          sprintf(form1,sqrt(weighted.var(data2[[xvar]],data2$weights))))
        } else{
            formula=as.formula(paste0(xvar,"~",yvar))
            result=summary(lm(formula,data=df1,weights=weights))$coef[2,4]
            if(result<0.001) {
                result="< 0.001"
            }else{
                result=sprintf("%.3f",result)
            }
        }
    } else if(is.numeric(temp)){
        if(treat==2){
            formula=as.formula(paste0("~",xvar,"+",yvar))
            result=svychisq(formula,design=svydesign(ids=~1,weights=~weights,data=df1))$p.value
            if(result<0.001) {
                result="< 0.001"
            }else{
                result=sprintf("%.3f",result)
            }
        } else{
            ratio=summatch[xvar,ifelse(treat,1,2)]
            result=paste0(round(ifelse(treat,ntreat,ncontrol)*ratio)," (",sprintf(form2,ratio*100),"%)")
        }
    } else if(length(unique(temp))>2){
        xvalues=sort(unique(df1[[xvar]]))
        if(treat==2){
            formula=as.formula(paste0("~",xvar,"+",yvar))
            result=svychisq(formula,design=svydesign(ids=~1,weights=~weights,data=df1))$p.value
            if(result<0.001) {
                result="< 0.001"
            }else{
                result=sprintf("%.3f",result)
            }
            result=c(result,rep("",length(xvalues)))
        } else{
            result=""
            for(j in seq_along(xvalues)){
                ratio=summatch[paste0(xvar,xvalues[j]),ifelse(treat,1,2)]
                result=c(result,paste0(round(ifelse(treat,ntreat,ncontrol)*ratio)," (",sprintf(form2,ratio*100),"%)"))
            }
        }
    } else{
        if(treat==2){
            formula=as.formula(paste0("~",xvar,"+",yvar))
            result=svychisq(formula,design=svydesign(ids=~1,weights=~weights,data=df1))$p.value
            if(result<0.001) {
                result="< 0.001"
            }else{
                result=sprintf("%.3f",result)
            }
        } else{
            ratio=summatch[xvar,ifelse(treat,1,2)]
            result=paste0(round(ifelse(treat,ntreat,ncontrol)*ratio)," (",sprintf(form2,ratio*100),"%)")
        }
    }
    result
}


