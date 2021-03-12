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
    title=c("Matching with twang::mnps()","Summary","Balance Check","Probability of Receiving Each Treatment",
            "Balance Check","t and chi-squared p-value",
            "Balance Table(1)","Balance Table(2)","Balance Table(3)")
    type=c("Pre","Rcode","plot","plot","plot","plot","Rcode","Rcode","Rcode")
    code=c(paste0("out<-",x),"summary(out)","p<-plot(out,plots=1)","p<-plot(out, plots = 2)",
           "p<-plot(out, plots = 3,pairwiseMax = FALSE)","print(plot(out, plots = 4))",
           "bal.table(out)","bal.table(out, collapse.to = 'covariate')","bal.table(out, collapse.to = 'stop.method')")
    if(dep!=""){
        title=c(title,"Estimation of Treatment Effect")
        type=c(type,"Rcode")
        code=c(code,paste0("estimateEffectTwang(out,dep='",dep,"',adjustCovar = ",adjustCovar,
                           ",method='",method,"')"))

    }
    data.frame(title,type,code,stringsAsFactors = FALSE)
}


#' Estimate Treatment Effect
#' @param out An object of class ps
#' @param dep Name of dependent variable
#' @param stop.method A method or methods of measuring and summarizing balance across pretreatment variables.
#' @param adjustCovar logical
#' @param method Method of PS estimate. One of c("GBM","lm","glm","chisq").
#' @importFrom twang get.weights
#' @importFrom survey svydesign svyglm
#' @importFrom stats confint contr.sum
#' @export
#' @examples
#'library(twang)
#'data(AOD)
#'out=mnps(treat ~ illact + crimjust + subprob + subdep + white,data = AOD,
#'stop.method=c("es.mean","ks.mean"),n.trees = 3000)
#'estimateEffectTwang(out,dep="suf12")
#'out1=mnps(treat ~ illact + crimjust + subprob + subdep + white,data = AOD,
#'estimand="ATT",treatATT="community",stop.method=c("es.mean","ks.mean"),n.trees = 3000)
#'estimateEffectTwang(out1,dep="suf12")
estimateEffectTwang=function(out,dep,stop.method="es.mean",adjustCovar=FALSE,method="GBM"){

    # out=ps.lalonde;stop.method = "es.mean";adjustCovar=TRUE;dep="re78";covars=""
    # dep="suf12"
    # stop.method = "es.mean";adjustCovar=FALSE;covars="";method="GBM"
    # method="GBM"
    # method="lm"
    # method="glm"
    if(class(out)=="ps"){
        yvar=out$gbm.obj$response.name
        xvars=out$gbm.obj$var.names
    } else if(class(out)=="mnps"){

       yvar=out$treat.var
       xvars=attr(out[[1]][[1]]$gbm.obj$Terms,"term.labels")
       stop.method=out$stopMethods[1]
    }
    temp=paste0(dep,"~",yvar)
    if(adjustCovar) {
        temp=paste0(temp,"+",paste0(xvars,collapse="+"))
    }

    form1=as.formula(temp)
    data1=out$data

    data1$w<-twang::get.weights(out,stop.method = stop.method)
    design.ps <- survey::svydesign(ids=~1, weights=~w, data=data1)
    if(method=="GBM"){

    if(out$estimand=="ATE"){
        txnames=levels(data1[[yvar]])
        txnames
        no=length(txnames)
        mylist=list()

        for(i in 1:no){
            temp1=paste0("svyglm(",temp,",design=design.ps,contrast=list(",yvar,"=contr.treatment(n=",no,",base=",i,")))")
            glm=eval(parse(text=temp1))
            result=as.data.frame(summary(glm)$coeff)
            result$OR=exp(glm$coef)
            result$lower=confint(glm)[,1]
            result$upper=confint(glm)[,2]
            row.names(result)=c("(Intercept)",paste0(txnames[-i],"-",txnames[i]))
            mylist[[i]]=result
        }
        mylist
        names(mylist)=paste0("Ref : ", levels(data1[[yvar]]))
        temp1=paste0("svyglm(",temp,",design=design.ps,contrast=list(",yvar,"=contr.sum))")
        glm2=eval(parse(text=temp1))
        result2=as.data.frame(summary(glm2)$coeff)
        result3=c(-sum(coef(glm2)[-1]),sqrt(c(-1,-1) %*% summary(glm2)$cov.scaled[-1,-1] %*% c(-1,-1)),
                  NA,NA)
        result4=rbind(result2,result3)
        result4$OR=exp(result4$Estimate)
        result4$lower=exp(result4[[1]]-1.96*result4[[2]])
        result4$upper=exp(result4[[1]]+1.96*result4[[2]])

        rownames(result4)=c("(Intercept)",levels(out$data[[yvar]]))
        result4
        mylist[[no+1]]=result4
        names(mylist)[no+1]="causal effect of each tx relative to the average potential outcome of all tx."
        result= mylist

    } else {
        glm=svyglm(form1,design=design.ps)
        result=as.data.frame(summary(glm)$coeff)
        result$OR=exp(glm$coef)
        result$lower=confint(glm)[,1]
        result$upper=confint(glm)[,2]
        result
    }
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
    if(method %in% c("lm","glm")){
        # summary(model)$coef
        result=as.data.frame(summary(model)$coef)
        result$lower=confint(model)[,1]
        result$upper=confint(model)[,2]
    }
    result


}
