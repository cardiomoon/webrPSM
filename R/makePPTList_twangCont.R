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

    type=c("Rcode","out","plot","Rcode","Rcode","plot")
    code=c(paste0("set.seed(",seed,")"),paste0("mod<-",x),"print(plot(mod,plots='optimize'))","summary(mod)",
           "twangContinuous::bal.table(mod, digits = 3)","print(plot(mod, plots='es'))")

    if(dep!=""){
        title=c(title,"Attach weights","Make design object","Outcome Model","Estimating causal Effect")
        type=c(type,"Rcode","Rcode","Rcode","Rcode")
        code=c(code,paste0(dataname,"$wts<-twangContinuous::get.weights(mod)"),
               paste0("design.ps <- svydesign(ids=~1, weights=~wts, data=",dataname,")"),
               paste0("outcome.model <- svyglm(",dep,"~",yvar,",design = design.ps, family = gaussian())"),
               "summary(outcome.model)")
    }
    data.frame(title,type,code)
}
