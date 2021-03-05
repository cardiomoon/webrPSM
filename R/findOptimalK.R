#' Find optimal number of groups for covariates for coarsened exact matching
#' @param out An object of class matchit with coarsened exact matching
#' @param plot logical Whether or not draw plot
#' @param k The number of groups for continuous covariates
#' @param mdthreshold numeric thershold for absolute standardized mean difference
#' @param nthreshold numeric Threshold of ratio of matched treated group and treated group
#' @importFrom MatchIt matchit
#' @importFrom ggplot2 aes_string geom_hline labs
#' @importFrom dplyr across
#' @importFrom egg ggarrange
#' @return A list
#' \itemize{
#' \item df - A data.frame of standardized mean differences
#' \item cutpoints - A list of cutpoints
#' \item k Recommended number of groups for continuous covariates
#' \item p1 The first ggplot
#' \item p2 The 2nd ggplot
#' }
#' @export
#' @examples
#' library(MatchIt)
#' out=matchit(formula=treat~V1+V2+V3,data=simData,method="cem")
#' estimateEffectZelig(out,dep="y")
#' result=findOptimalK(out)
#' out=matchit(formula=treat~V1+V2+V3,data=simData,method="cem",cutpoints=result$cutpoints)
#' estimateEffectZelig(out,dep="y")
#' out=matchit(treat ~ age + educ + race + married+nodegree + re74 + re75, data =lalonde,
#' method="cem",estimand="ATE")
#' findOptimalK(out)
findOptimalK=function(out,plot=FALSE,k=2:10,mdthreshold=0.1,nthreshold=0.97){
     # plot=TRUE;k=2:10;mdthreshold=0.1;nthreshold=0.97

    cat("\nInitial call: ")
    print(out$call)
    cat("\nTested optimal number of k :",k,"\n")
    cat("\nThreshold for standardized mean difference :", mdthreshold,"\n")
    cat("\nThreshold of ratio matched treated group and treated group :",nthreshold,"\n")
    temp1=formula2vars(out$formula)
    orgData=eval(parse(text=call2param(out)$data))
    xvars1=temp1$xvars
    xvars=c()
    for(i in seq_along(xvars1)){
        if(is.numeric(orgData[[xvars1[i]]])) xvars=c(xvars,xvars1[i])
    }
    xvars
    smd=matrix(NA,nrow=length(k),ncol=length(xvars))
    colnames(smd)=xvars
    rownames(smd)=paste0("k",k)
    cases=rep(NA,length(k))
    smd
    estimand=call2param(out$call)$estimand

    for(i in seq_along(k)){
        cutpoints=list()
        for(j in seq_along(xvars)){
            cutpoints[[xvars[j]]]=k[i]
        }
        out=matchit(formula=out$formula,data=orgData,method="cem",estimand=estimand,cutpoints=cutpoints)
        summary(out,standardize=TRUE)$sum.matched[,3]
        for(j in seq_along(xvars)){
        smd[k[i]-1,xvars[j]]=abs(summary(out,standardize=TRUE)$sum.matched[xvars[j],3])

        }
        cases[k[i]-1]=out$nn[4,2]
    }
    df1=data.frame(smd)
    df1
    df2=df1 %>% mutate(across(1:length(xvars),~ifelse(.x<mdthreshold,1,0)))
    df2$ok1=rowSums(df2)
    df2
    df1$k=k
    df3<-df1 %>% pivot_longer(cols=-k,names_to = "covariates")
    p1<-ggplot(df3,aes_string(x="k",y="value",color="covariates"))+geom_line()+geom_point()+
        geom_hline(yintercept=mdthreshold,linetype=2)+theme_bw()+
        labs(y="standardized mean differences")+
        theme(legend.position="top")
    p1
    df1$cases=cases
    p2<-ggplot(df1,aes_string(x="k",y="cases"))+geom_line()+geom_point()+theme_bw()+
        geom_hline(yintercept=max(df1$cases)*nthreshold,linetype=2)
    if(plot) egg::ggarrange(p1,p2)
    df1
    df1$ok1<-df2$ok1
    df1
    df1<-df1 %>% mutate(ok2=ifelse(cases>max(cases)*nthreshold,1,0),ok=.data$ok1+.data$ok2)
    df1

    k=min(df1$k[df1$ok==max(df1$ok)])
    cutpoints=list()
    for(j in seq_along(xvars)){
        cutpoints[[xvars[j]]]=k
    }
    result=list(df=df1,cutpoints=cutpoints,k=k,p1=p1,p2=p2)
    cat("\nThe optimal number of groups for continuous covariates : ",k,"\n\n")
    print(result$df)
    invisible(result)
}
