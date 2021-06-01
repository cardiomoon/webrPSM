#' Draw balance plot
#' @param x An object of class ps.cont
#' @importFrom twangContinuous bal.table
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 aes_string element_blank geom_text
#' @examples
#' library(twangContinuous)
#' model=ps.cont( treat ~ x1 + x2,data=simData2,n.trees=5000)
#' balance.plot(model)
#' \dontrun{
#' data(dat)
#' model=ps.cont(tss_0~sfs8p_0+sati_0+sp_sm_0+recov_0+subsgrps_n+treat,data=dat,n.trees = 500)
#' balance.plot(model)
#' }
#' @export
balance.plot=function(x){
    res=twangContinuous::bal.table(x, digits = 3)
    colnames(res)=c("Unweighted","Weightd")
    res$id=row.names(res)
    res1<-res %>% pivot_longer(cols=1:2)
    res1$value=abs(res1$value)
    res1$hjust=0
    res1$hjust[res1$name=="Unweighted"]=1.1
    res1$hjust[res1$name=="Weightd"]=-0.1
    ggplot(res1,aes_string(x="name",y="value",group="id",color="id"))+
        geom_hline(yintercept=0.1,color="grey80",alpha=0.5)+
        geom_line()+
        geom_text(aes_string(label="id",hjust="hjust"))+
        guides(color=FALSE)+
        theme_bw()+
        labs(x=NULL,y="Absolute Correlation")+
        theme(panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),
        )

}

#' Make PPT list with continuous variable
#' @param x string
#' @param dep character Name of dependent variable
#' @param seed numeric random seed
#' @importFrom survey svydesign svyglm
#' @export
#' @examples
#' x="ps.cont( treat ~ x1 + x2,data=simData2,n.trees=5000)"
#' result=makePPTList_twangCont(x,dep="y")
#' x="ps.cont(tss_0~sfs8p_0+sati_0+sp_sm_0+recov_0+subsgrps_n+treat,data=dat,n.trees = 500)"
#' result=makePPTList_twangCont(x,dep="sfs8p_3",seed=1234)
makePPTList_twangCont=function(x,dep="",seed=1234){
    x=gsub(" ","",x)
    temp=unlist(strsplit(x,"(",fixed=TRUE))[2]
    temp2=unlist(strsplit(temp,"~"))
    yvar=temp2[1]
    temp3=unlist(strsplit(temp2[2],","))[2]
    dataname=unlist(strsplit(temp3,"="))[2]

    title=c("Set seed number","PS estimation","Assessing Diagnostics","Summary of model","Balance Assessment","Balance plot")

    type=c("Rcode","out","plot","Rcode","Rcode","ggplot")
    code=c(paste0("set.seed(",seed,")"),paste0("mod<-",x),"print(plot(mod,plots='optimize'))","summary(mod)",
           "twangContinuous::bal.table(mod, digits = 3)","balance.plot(mod)")

    if(dep!=""){
        title=c(title,"Attach weights","Make design object","Outcome Model","Estimating Causal Effect","Confidence Interval")
        type=c(type,"Rcode","Rcode","Rcode","Rcode","Rcode")
        code=c(code,paste0(dataname,"$wts<-twangContinuous::get.weights(mod)"),
               paste0("design.ps <- svydesign(ids=~1, weights=~wts, data=",dataname,")"),
               paste0("outcome.model <- svyglm(",dep,"~",yvar,",design = design.ps, family = gaussian())"),
               "summary(outcome.model)","confint(outcome.model)")
    }
    data.frame(title,type,code)
}
