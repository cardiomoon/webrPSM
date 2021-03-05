#' Estimate Effect after subclass matching with simulation
#' @param out An object of class "matchit"
#' @param mode One of c("continuous","binary","survival")
#' @param multiple logical Whether or not perform multiple regression
#' @param dep Name of dependent variable
#' @param seed numeric seed number
#' @param n_sim numeric number of simuations
#' @param reverse logical
#' @importFrom dplyr summarize
#' @export
#' @examples
#' library(MatchIt)
#' out=matchit(formula=treat~V1+V2+V3,data=simData,estimand="ATE",link="linear.logit",
#' method="subclass",subclass=8)
#' estimateEffectZeligSubclass(out,dep="y")
estimateEffectZeligSubclass=function(out,mode="continuous",multiple=TRUE,dep,seed=1224,n_sim=10000,reverse=FALSE){
     # mode="continuous";multiple=TRUE;dep="y";seed=1224;n_sim=10000;reverse=FALSE
    set.seed(seed)
    md=match.data(out)
    temp1=formula2vars(out$formula,allowInteraction=TRUE)
    xvars=temp1$xvars
    yvar=temp1$yvar
    estimand=call2param(out$call)$estimand

    final=list()
    for(i in seq_along(dep)){

        deselect=c()
        temp=paste0(dep[i],"~",yvar)
        if(multiple){
            temp=paste0(temp,"+",paste0(xvars,collapse="+"))
        }
        form1=as.formula(temp)
        result=list()
        j=1
        for(j in 1:out$info$subclass){

            temp_md=md[md[["subclass"]]==j,]
            z_model=zelig(form1,
                          data=temp_md,
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
            est_portion=sum(temp_md$weights)/sum(md$weights)
            est_sub=data.frame(subclass=rep(j,n_sim),sim_ev=est*est_portion,sim_num=1:n_sim)
            est_sub
            if(j==1) {
                result=est_sub
            } else{
                result=rbind(result,est_sub)
            }
        }
        names(result)[2]="sim_ev"
        result_agg <- result %>% group_by(.data$sim_num) %>% summarize(est=sum(.data$sim_ev))
        result_agg
        res=quantile(result_agg$est,p=c(0.025,0.5,0.975))
        res=data.frame(est=res[2],lower=res[1],upper=res[3])
        if(i==1) {
            final=res
        } else{
            final=rbind(final,res)
        }
    }
    final
    row.names(final)=dep
    final$estimand=estimand
    if(reverse) final$estimand="ATC"
    final$method=out$info$method
    final
}
