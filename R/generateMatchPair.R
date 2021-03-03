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
generateMatchPair=function(out,dep="y"){
    dataName=call2param(out$call)$data
    y=eval(parse(text=paste0(dataName,"$",dep)))

    myModify=function(x){
        x=as.numeric(as.character(x))
        x[x<0]=NA
        x
    }
    findX=function(x){
        y[x]
    }
    temp=out$match.matrix
    df=data.frame(temp)
    df[["treat"]]=row.names(df)
    df
    df %>% dplyr::select(starts_with("treat"),everything()) %>%
        mutate_all(myModify) %>%
        drop_na() %>%
        mutate_all(findX)
}

#' Find significant gamma range
#' @param data A data.frame as a result of generateMatchPair()
#' @param start numeric start gamma value
#' @param threshold numeric p value threshold
#' @param method method to be passed to senmw()
#' @importFrom sensitivitymw senmw
#' @export
#' @examples
#' library(MatchIt)
#' out=matchit(formula=treat~V1+V2+V3,data=simData,link="linear.logit",caliper=0.15,ratio=2)
#' df=generateMatchPair(out,dep="y")
#' tail(gammaRangeSearch(df))
#' out=matchit(formula=treat~V1+V2+V3,data=simData,link="linear.logit",ratio=2,method="optimal")
#' df=generateMatchPair(out,dep="y")
#' tail(gammaRangeSearch(df))
gammaRangeSearch=function(data,start=1,threshold=0.025,method="w"){
    # start=1;method="w"
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