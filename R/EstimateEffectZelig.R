#' Estimate Effect after matching with simulation
#' @param out An object of class "matchit"
#' @param mode One of c("continuous","binary","survival")
#' @param multiple logical Whether or not perform multiple regression
#' @param dep Name of dependent variable
#' @param seed numeric seed number
#' @param n_sim numeric number of simuations
#' @param reverse logical
#' @importFrom Zelig zelig setx sim get_qi
#' @importFrom stats quantile
#' @export
#' @examples
#' library(MatchIt)
#' out=matchit(formula=treat~V1+V2+V3,data=simData,link="linear.logit",
#' caliper=0.15,ratio=2)
#' estimateEffectZelig(out,dep="y")
#' out1=matchit(formula=treat~V1+V2+V3,data=simData,link="linear.logit",estimand="ATC",
#' caliper=0.15,ratio=2)
#' estimateEffectZelig(out1,dep="y")
#' out2=matchit(formula=Rtreat~V1+V2+V3,data=simData,link="linear.logit",
#' caliper=0.15,ratio=2,replace=TRUE)
#' estimateEffectZelig(out2,dep="y",reverse=TRUE)
#' fullATT=matchit(formula=treat~V1+V2+V3,data=simData,
#' link="linear.logit",method="full")
#' estimateEffectZelig(fullATT,dep="y")
#' fullATC=matchit(formula=treat~V1+V2+V3,data=simData,link="linear.logit",
#' method="full",estimand="ATC")
#' estimateEffectZelig(fullATC,dep="y")
#' fullATC2=matchit(formula=Rtreat~V1+V2+V3,data=simData,link="linear.logit",
#' method="full")
#' estimateEffectZelig(fullATC2,dep="y",reverse=TRUE)
#' fullATE=matchit(formula=treat~V1+V2+V3,data=simData,link="linear.logit",
#' method="full",estimand="ATE")
#' estimateEffectZelig(fullATE,dep="y")
estimateEffectZelig=function(out,mode="continuous",multiple=TRUE,dep,seed=1224,n_sim=10000,reverse=FALSE){
    # mode="continuous";multiple=TRUE;dep="y";seed=1224;n_sim=10000
    set.seed(seed)
    md=match.data(out)
    temp1=formula2vars(out$formula)
    xvars=temp1$xvars
    yvar=temp1$yvar
    estimand=call2param(out$call)$estimand

    for(i in seq_along(dep)){

        deselect=c()
        temp=paste0(dep[i],"~",yvar)
        if(multiple){
            form1=out$formula
        } else{
            form1=as.formula(temp)
        }
        z_model=zelig(form1,
                      data=md,
                      model='ls',weights="weights",cite=FALSE)


        #x0=setx(z_model,treat=0,data=md)
        #x1=setx(z_model,treat=1,data=md)
        temp1=paste0("setx(z_model,",yvar,"=0,data=md)")
        temp2=paste0("setx(z_model,",yvar,"=1,data=md)")

        x0=eval(parse(text=temp1))
        x1=eval(parse(text=temp2))
        if(reverse){
            x0=eval(parse(text=temp2))
            x1=eval(parse(text=temp1))
        }

        s0=sim(z_model,x0,num=n_sim)
        s1=sim(z_model,x1,num=n_sim)

        est=get_qi(s1,"ev")-get_qi(s0,"ev")
        res=quantile(est,p=c(0.025,0.5,0.975))
        res
        res=data.frame(est=res[2],lower=res[1],upper=res[3])
        if(i==1) {
            result=res
        } else{
            result=rbind(result,res)
        }
    }
    rownames(result)=dep
    result$estimand=estimand
    if(reverse) result$estimand="ATC"
    result$method=out$info$method
    attr(result,"est")=est
    result
}
