#' Draw balance Plot
#' @param out An object of class mnps
#' @param show.label logical Whether or not show label
#' @param show.legend logical Whether or not show legend
#' @param rowByGroup logical If true, draw row by group
#' @importFrom ggplot2 aes_string facet_grid geom_hline labs
#' @examples
#' library(twang)
#' data(AOD)
#' out=mnps(treat~illact+crimjust+subprob+subdep+white,data=AOD,n.trees=3000,verbose=FALSE)
#' balancePlotTwang(out)
#' out=mnps(treat~illact+crimjust+subprob+subdep+white,data=AOD,n.trees=3000,verbose=FALSE,
#' stop.method=c("es.mean","ks.max"))
#' balancePlotTwang(out)
#' @export
balancePlotTwang=function(out,show.label=FALSE,show.legend=FALSE,rowByGroup=TRUE){
    res=twang::bal.table(out)
    res$group=paste0(res$tmt1," vs ",res$tmt2)
    method=unique(res$stop.method)
    method=setdiff(method,"unw")
    count=length(method)
    res$method=res$stop.method
    res$method[res$stop.method=="unw"]=method[1]
    if(count>1){
        for(i in 1:(count-1)){
            res1=res[res$stop.method=="unw",]
            res1$method=method[i+1]
            res=rbind(res,res1)
        }
    }
    res$stop.method=factor(res$stop.method,levels=unique(res$stop.method))
    res$hjust=1.1
    res$hjust[res$stop.method!="unw"]=-0.1
    res$method1="Unweighted"
    res$method1[res$stop.method!="unw"]="Weighted"
    p<-ggplot(res,aes_string(x="method1",y="std.eff.sz",group="var",color="var"))+
        geom_point()+
        geom_line(alpha=0.5)+
        geom_hline(yintercept=0.1,color="grey80",alpha=0.5)+
        theme_bw()+
        theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
        labs(y="Absolute Standard Difference",x=NULL)
    if(show.label) p<-p+  geom_text(aes_string(label="var",hjust="hjust"))
    if(!show.legend) p<-p+guides(color=FALSE)
    if(count>1) {
        if(rowByGroup){
          p<-p+facet_grid(method~group,scales="free",space="free")
        } else{
          p<-p+facet_grid(group~method,scales="free",space="free")
        }
    } else{
        if(rowByGroup){
          p<-p+facet_grid(~group,scales="free",space="free")
        } else{
            p<-p+facet_grid(group~.,scales="free",space="free")
        }
    }

    p
}
