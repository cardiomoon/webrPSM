#' Generate list for sensitivity test
#' @param out An object of class matchit
#' @param dep Name of dependent variable
#' @importFrom dplyr count group_by mutate ungroup row_number filter
#' @importFrom rlang .data
#' @export
#' @examples
#' library(MatchIt)
#' out=matchit(formula=treat~V1+V2+V3,data=simData,link="linear.logit",method="full")
#' generateMatchPairFull(out,dep="y")
generateMatchPairFull=function(out,dep){
    md=match.data(out)
    data=out$model$data
    data$full_subclass<-md$subclass
    temp_row_count<-data %>% count(.data$full_subclass)
    temp_col_count<-data %>% group_by(.data$full_subclass) %>%
        mutate(rid=row_number()) %>% ungroup()

    row_number=dim(temp_row_count)[1]
    col_number=max(temp_col_count$rid)
    mymatrix=matrix(NA,nrow=row_number,ncol=col_number)
    mymatrix
    myTCstatus=rep(NA,row_number)
    i=1
    for(i in 1:row_number){
        temp_subclass=data %>% dplyr::filter(.data$full_subclass==i)
        temp_subclass
        y_trt=temp_subclass[[dep]][temp_subclass$treat==1]
        y_ctrl=temp_subclass[[dep]][temp_subclass$treat==0]

        CtoT=c(y_ctrl,y_trt)
        TtoC=c(y_trt,y_ctrl)

        lengthC=rep(length(y_ctrl),length(CtoT))
        lengthT=rep(length(y_trt),length(TtoC))
        y_row=ifelse(lengthC>=lengthT,TtoC,CtoT)
        whichFirst=length(y_ctrl)>=length(y_trt)

        mymatrix[i,1:length(y_row)]=y_row
        myTCstatus[i]=whichFirst
    }
    list(mymatrix,myTCstatus)
}


#' Find significant gamma range for full matching
#' @param out An object of class matchit
#' @param dep Name of dependent variable
#' @param start numeric start gamma value
#' @param threshold numeric p value threshold
#' @importFrom sensitivityfull senfm
#' @export
#' @examples
#' library(MatchIt)
#' out=matchit(formula=treat~V1+V2+V3,data=simData,link="linear.logit",method="full")
#' tail(gammaRangeSearchFull(out,dep="y"))
#' \dontrun{
#' out=matchit(formula=treat~V1+V2+V3,data=simData,link="linear.logit",estimand="ATC",method="full")
#' tail(gammaRangeSearchFull(out,dep="y"))
#' out=matchit(treat ~ age + educ + race + married+nodegree + re74 + re75, data =lalonde,
#'    method="full")
#' tail(gammaRangeSearchFull(out,dep="re78"))
#' }
gammaRangeSearchFull=function(out,dep,start=1,threshold=0.025){
    mygamma=start
    pvalue=0
    res=generateMatchPairFull(out,dep=dep)
    myresult=data.frame()

    while(pvalue<threshold){
    result=data.frame(sensitivityfull::senfm(res[[1]],res[[2]],gamma=mygamma))
    pvalue=result$pval
    temp=data.frame(gamma=mygamma,pval=result$pval)
    myresult=rbind(myresult,temp)

    mygamma=mygamma+0.1
    }
    myresult
}

