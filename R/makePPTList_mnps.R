#' Make PPT List with an object of class mnps
#' @param x Character
#' @param dep Character Name of dependent variable
#' @param adjustCovar logical
#' @param covars string
#' @param method Method of PS estimate. One of c("GBM","lm","glm","chisq")
#' @export
#' @examples
#' library(twang)
#' data(AOD)
#' x <- "mnps(treat~illact+crimjust+subprob+subdep+white,data=AOD,verbose=FALSE,n.trees=3000)"
#' result=makePPTList_mnps(x,dep="suf12")
makePPTList_mnps=function(x,dep="",adjustCovar=FALSE,covars="",method="GBM"){
    # dep="suf12";adjustCovar=FALSE;covars="";method="GBM"
    #out<-eval(parse(text=x))
    title=c("Summary","Balance Check","Probability of Receiving Each Treatment",
            "Balance Check","t and chi-squared p-value",
            "Balance Table(1)","Balance Table(2)","Balance Table(3)")
    type=c("Rcode","plot","plot","plot","plot","Rcode","Rcode","Rcode")
    code=c(paste0("out<-",x,";summary(out)"),"p<-plot(out,plots=1)","p<-plot(out, plots = 2)",
           "p<-plot(out, plots = 3,pairwiseMax = FALSE)","print(plot(out, plots = 4))",
           "bal.table(out)","bal.table(out, collapse.to = 'covariate')","bal.table(out, collapse.to = 'stop.method')")
    if(dep!=""){
        title=c(title,"Estimation of Treatment Effect")
        type=c(type,"Rcode")
        code=c(code,paste0("estimateEffectTwang(out,dep='",dep,"',adjustCovar = ",adjustCovar,
                           ",covars=c('",paste0(covars,collapse="','"),"'),method='",method,"')"))
    }
    data.frame(title,type,code,stringsAsFactors = FALSE)
}


#' Estimate Treatment Effect
#' @param out An object of class ps
#' @param dep Name of dependent variable
#' @param stop.method A method or methods of measuring and summarizing balance across pretreatment variables.
#' @param adjustCovar logical
#' @param covars string
#' @param method Method of PS estimate. One of c("GBM","lm","glm","chisq").
#' @importFrom twang get.weights
#' @importFrom survey svydesign svyglm
#' @importFrom stats confint
#' @export
#' @examples
#'library(twang)
#'data(AOD)
#'out=mnps(treat ~ illact + crimjust + subprob + subdep + white,data = AOD,n.trees = 3000)
#'estimateEffectTwang(out,dep="suf12",adjustCovar=FALSE)
#'estimateEffectTwang(out,dep="suf12",adjustCovar=TRUE)
estimateEffectTwang=function(out,dep,stop.method="es.mean",adjustCovar=FALSE,covars="",method="GBM"){

    # out=ps.lalonde;stop.method = "es.mean";adjustCovar=TRUE;dep="re78";covars=""
    # out=mnps.AOD;dep="suf12";stop.method = "es.mean";adjustCovar=FALSE;covars=""
    # method="GBM"
    # method="lm"
    # method="glm"

    if(class(out)=="ps"){
        yvar=out$gbm.obj$response.name
        xvars=out$gbm.obj$var.names
    } else if(class(out)=="mnps"){
        yvar=colnames(out$data)[1]
        xvars=colnames(out$data)[-1]
        stop.method=out$stopMethods[1]
    }
    temp=paste0(dep,"~",yvar)
    if(covars!=""){
        temp=paste0(temp,"+",paste0(covars,collapse="+"))
    } else if(adjustCovar) {
        temp=paste0(temp,"+",paste0(xvars,collapse="+"))
    }
    temp
    form1=as.formula(temp)
    data1=out$data

    if(method=="GBM"){
        data1$w<-twang::get.weights(out,stop.method = stop.method)
        design.ps <- survey::svydesign(ids=~1, weights=~w, data=data1)
        model<-survey::svyglm(form1,design=design.ps)
    } else if(method=="lm"){
        model<-lm(form1,data=data1)
    } else if(method=="glm"){  # glm
        form2=paste0(yvar,"~",paste0(xvars,collapse="+"))
        out1=glm(form2,data=data1,family=binomial)
        data1$w.logit<-rep(1,nrow(data1))
        data1$w.logit[data1[[yvar]]==0] <-exp(predict(out1,subset(data1,data1[[yvar]]==0)))
        design.logit=svydesign(ids=~1,weights=~w.logit,data=data1)
        model<-svyglm(form1,design=design.logit)
    } else if(method=="chisq"){
        data1$w<-twang::get.weights(out,stop.method = stop.method)
        design.ps <- survey::svydesign(ids=~1, weights=~w, data=data1)
        # result=svychisq(as.formula(paste0("~",dep,"+",yvar)),design=design.ps)
        result=eval(parse(text=paste0("svychisq(~",dep,"+",yvar,",design=design.ps)")))

    }
    if(method!="chisq"){
        if(class(out)=="mnps"){
            result=as.data.frame(summary(model)$coef)
            result$OR=exp(result$Estimate)

            result$lower=exp(confint(model)[,1])
            result$upper=exp(confint(model)[,2])
        } else{
            result=as.data.frame(summary(model)$coef[yvar,,drop=FALSE])
            result$lower=confint(model)[yvar,1]
            result$upper=confint(model)[yvar,2]
        }
    }
    result
}
