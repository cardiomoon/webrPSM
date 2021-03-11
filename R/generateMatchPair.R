#' Generate data.frame for sensitivity test
#' @param out An object of class matchit
#' @param dep Name of dependent variable
#' @importFrom dplyr select mutate_all
#' @importFrom tidyr drop_na
#' @importFrom tidyselect everything starts_with
#' @export
#' @examples
#' library(MatchIt)
#' out=matchit(formula=treat~V1+V2+V3,data=simData,link="linear.logit",caliper=0.15,ratio=2)
#' generateMatchPair(out,dep="y")
#' out=matchit(treat~age+educ+race+married+nodegree+re74+re75,data=lalonde)
#' generateMatchPair(out,dep="re78")
generateMatchPair=function(out,dep="y"){
    # dataName=call2param(out$call)$data
    # y=eval(parse(text=paste0(dataName,"$",dep)))
     # out=matchit(treat~age+educ+race+married+nodegree+re74+re75,data=lalonde)
     # dep="re78"
     # out=matchit(formula=treat~V1+V2+V3,data=simData,link="linear.logit",caliper=0.15,ratio=2)
     # dep="y"
    y=out$model$data[dep]

    myModify=function(x){
        x1=suppressWarnings(as.numeric(as.character(x)))
        if(sum(is.na(x1))==length(x1)){
            result=x
        } else{
            x1[x1<0]=NA
            result=x1
        }
        result
    }
    findX=function(x){
        y[x,dep]
    }
    temp=out$match.matrix

    df=data.frame(temp)
    df[["treat"]]=row.names(df)

    df %>% dplyr::select(starts_with("treat"),everything()) %>%
        mutate_all(myModify) %>%
        drop_na() %>%
        mutate_all(findX)
}

#' Find significant gamma range
#' @param out An object of class matchit
#' @param dep Character Name of dependent variable
#' @param start numeric start gamma value
#' @param threshold numeric p value threshold
#' @param method method to be passed to senmw()
#' @importFrom sensitivitymw senmw
#' @export
#' @examples
#' library(MatchIt)
#' out=matchit(formula=treat~V1+V2+V3,data=simData,link="linear.logit",caliper=0.15,ratio=2)
#' tail(gammaRangeSearch(out,dep="y"))
#' out=matchit(formula=treat~V1+V2+V3,data=simData,link="linear.logit",ratio=2,method="optimal")
#' tail(gammaRangeSearch(out,dep="y"))
#' out=matchit(treat~age+educ+race+married+nodegree+re74+re75,data=lalonde)
#' tail(gammaRangeSearch(out,dep="re78"))
gammaRangeSearch=function(out,dep="y",start=1,threshold=0.025,method="w"){
    # start=1;method="w"
    data=generateMatchPair(out,dep=dep)
    mygamma=start
    pvalue=0
    myresult=data.frame()
    while(pvalue<threshold){
        result=data.frame(senmw(data,gamma=mygamma,method=method))
        pvalue=result$pval
        temp=data.frame(gamma=mygamma,pval=result$pval)
        myresult=rbind(myresult,temp)
        mygamma=mygamma+0.1
    }
    myresult
}
