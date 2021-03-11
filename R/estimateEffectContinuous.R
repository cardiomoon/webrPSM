#' Add IPW(inverse probability weight) to a data.frame
#' @param formula A formula
#' @param data A data.frame
#' @importFrom stats dnorm
#' @export
#' @examples
#' formula=treat~x1+x2
#' mydata=addIPW(treat~x1+x2,data=simData2)
addIPW=function(formula, data){

    formula=treat~x1+x2
    temp=gsub(" ","",deparse(formula))
    temp=unlist(strsplit(temp,"~"))
    treatvar=temp[1]
    xvars=unlist(strsplit(temp[2],"+",fixed=TRUE))
    mydata=data
    lmGPS=lm(formula,data=mydata)

    mydata$GPS=dnorm(mydata[[treatvar]],mean=lmGPS$fitted,sd=summary(lmGPS)$sigma)
    mydata$numerator=dnorm(mydata[[treatvar]],mean=mean(mydata$treat),sd=sd(mydata$treat))
    mydata$IPW = mydata$numerator/mydata$GPS
    attr(mydata,"xvars")=xvars
    attr(mydata,"treatvar")=treatvar
    mydata
}


#'Standardize
#'@param x A vector
#'@param max.ylev Numeric maximum length of unique values to be treated as a dummy variable
#'@importFrom stats sd
#'@export
standardize=function(x,max.ylev=2){
    if(is.numeric(x) && (length(unique(x))>max.ylev)){
        result=(x-mean(x))/sd(x)
    } else{
        result=x
    }
    result
}



#' Covariates balance check
#' @param xvars Names of covariates
#' @param treatvar Name of treatment var
#' @param data A data.frame
#' @importFrom dplyr mutate_at vars all_of
#' @export
#' @examples
#' mydata=addIPW(treat~x1+x2,data=simData2)
#' balanceCheck(mydata)
balanceCheck=function(data,xvars=NULL,treatvar=NULL){
      # xvars=NULL;treatvar=NULL
    if(is.null(xvars)) xvars=attr(data,"xvars")
    if(is.null(treatvar)) treatvar=attr(data,"treatvar")
    all_of(c(xvars,treatvar))

    stdData=data %>% mutate_at(all_of(c(xvars,treatvar)),standardize)
    coef=c()

    for(i in seq_along(xvars)){
        form2=paste0(xvars[i],"~",treatvar)
        coef=c(coef,lm(as.formula(form2),data=stdData,weights=stdData[["IPW"]])$coef[2])
    }
    df=data.frame(covariate=xvars,coef=coef)
    df$balanced=df$coef<0.1
    df
}

#' estimate effect based on IPW(inverse probability weight)
#' @param data A data.frame as a result of addIPW()
#' @param dep Name of dependent variable
#' @param xvars Names of covariates
#' @param treatvar Name of treatment var
#' @param seed A single integer
#' @param probs numeric vector of probabilities with values in [0,1].
#' @param num an integer specifying the number of simulations to compute
#' @param weights Name of weight variable
#' @importFrom Zelig zelig
#' @export
#' @examples
#' mydata=addIPW(treat~x1+x2,data=simData2)
#' estimateEffectContinuous(mydata,dep="y",weights="IPW")
#' estimateEffectContinuous(mydata,dep="y")
estimateEffectContinuous=function(data,dep,xvars=NULL,treatvar=NULL,seed=1234,probs=0.1*(1:9),num=10000,weights=NULL){
    # data=addIPW(simData2,xvars=c("x1","x2"),treatvar="treat")
    # dep="y"
    # xvars=NULL;treatvar=NULL;seed=1234;probs=0.1*(1:9);num=10000;weights=NULL
    set.seed(seed)
    if(is.null(xvars)) xvars=attr(data,"xvars")
    if(is.null(treatvar)) treatvar=attr(data,"treatvar")
    form1=paste0(dep,"~",treatvar,"+",paste0(xvars,collapse="+"))
    #print(form1)
    if(is.null(weights)){
        out=zelig(as.formula(form1),data=data,model="ls",cite=FALSE)
    } else{
        out=zelig(as.formula(form1),data=data,model="ls",weights=weights,cite=FALSE)
    }
    as.data.frame(zelig2est(out))
}

#' Estimate dose-response function
#' @param out A object as a result of zelig()
#' @param probs numeric vector of probabilities with values in [0,1]
#' @param num an integer specifying the number of simulations to compute
#' @importFrom dplyr vars
#' @export
#' @examples
#' library(Zelig)
#' mydata=addIPW(treat~x1+x2,data=simData2)
#' out=zelig(y~treat+x1+x2,data=mydata,model="ls",weights="IPW",cite=FALSE)
#' zelig2est(out)
zelig2est=function(out,probs=0.1*(1:9),num=10000){
    result=data.frame()
    orgData=out$originaldata
    treatVar=formula2vars(out$formula)$xvars[1]
    range1=quantile(orgData[[treatVar]],probs=probs)
    #print(range1)
    for(i in seq_along(range1)){
        x=setx(out,treat=range1[i],data=orgData)
        s=sim(out,x,num=num)
        ev=data.frame(t(get_qi(s,"ev")))
        result=rbind(result,ev)
    }

    names(result)=paste0("sim",1:num)
    result[["treat"]]=range1

    result %>%
        pivot_longer(cols=starts_with("sim")) %>%
        group_by(.data$treat) %>%
        summarize(lower=quantile(.data$value,p=0.025),
                  est=quantile(.data$value,p=0.5), upper=quantile(.data$value,p=0.975))

}

#' Draw plot comparing IPW and linear regression
#' @param data A data.frame as a result of addIPW()
#' @param dep Name of dependent variable
#' @param seed A single integer
#' @param print logical
#' @importFrom dplyr bind_rows
#' @importFrom ggplot2 geom_ribbon
#' @export
#' @examples
#' set.seed(1234)
#' mydata=addIPW(treat~x1+x2,data=simData2)
#' plotCompareEffects(mydata)
plotCompareEffects=function(data,dep="y",seed=1234,print=TRUE){
    # if(is.null(xvars)) xvars=attr(data,"xvars")
    # if(is.null(treatvar)) treatvar=attr(data,"treatvar")
    # set.seed(seed)
    # OLS_estimate=zelig2est(zelig(y~treat+x1+x2,data=mydata,model="ls",weights="IPW",cite=FALSE))
    # set.seed(seed)
    # IPW_estimate=zelig2est(zelig(y~treat+x1+x2,data=mydata,model="ls",cite=FALSE))
    # data=addIPW(simData2,xvars=c("x1","x2"),treatvar="treat")
    # dep="y";seed=1234;print=TRUE
    OLS_estimate=estimateEffectContinuous(data,dep=dep,seed=seed)
    IPW_estimate=estimateEffectContinuous(data,dep=dep,weights="IPW",seed=seed)

    df=bind_rows(OLS_estimate %>% mutate(model="OLS"),
             IPW_estimate %>% mutate(model="IPW"))
    df
    p<-ggplot(df,aes_string(x="treat",y="est",fill="model"))+
    geom_line(aes_string(color="model"))+
    geom_ribbon(aes_string(ymin="lower",ymax="upper"),alpha=0.3)+
    theme_bw()+
    theme(legend.position="top")

    if(print) print(p)

    invisible(list(
        result=df,
        plot=p))

}

#' make pptList with Matching with IPW
#' @param x Character
#' @param dep Name of dependent variable
#' @param seed Numeric
#' @export
#' @examples
#' result=makePPTList_IPW(x="addIPW(treat~x1+x2,data=simData2)",dep="y",seed=1234)
makePPTList_IPW=function(x,dep,seed=1234){
       # x="addIPW(treat~x1+x2,data=simData2)";dep="y";seed=1234
    title<-type<-code<-c()

    title="Calculate IPW and add to data"
    type="Rcode"
    code=paste0("mydata<-",x,";head(mydata)")

    title=c(title,"Covariates Balance Check")
    type=c(type,"Rcode")
    code=c(code,"balanceCheck(mydata)")

    title=c(title,"Estimate Dose-Response Function")
    type=c(type,"Rcode")
    code=c(code,paste0("estimateEffectContinuous(mydata,dep='",dep,"',weights='IPW',seed=",seed,")"))

    title=c(title,"Comparing Effect with Linear Model")
    type=c(type,"Rcode")
    code=c(code,paste0("plotCompareEffects(mydata,dep='",dep,"',seed=",seed,",print=FALSE)$result"))

    title=c(title,"Comparing Effect with Linear Model")
    type=c(type,"ggplot")
    code=c(code,paste0("plotCompareEffects(mydata,dep='",dep,"',seed=",seed,",print=FALSE)$plot"))

    data.frame(title=title,type=type,code=code,stringsAsFactors = FALSE)
}
