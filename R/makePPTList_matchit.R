#' Draw summary plot of propensity score matching
#' @param x An object of class matchit
#' @param show.table logical whether or not show table
#' @param xpos numeric Min x-axis position of table
#' @param ypos numeric Min y-axis position of table
#' @importFrom ggplot2 annotation_custom guides geom_line
#' @importFrom ggrepel geom_text_repel
#' @importFrom gridExtra tableGrob ttheme_minimal
#' @export
#' @examples
#' require(MatchIt)
#' formula=treat ~ age + race+educ + nodegree + re74 + re75
#' x=matchit(formula, data =lalonde, method= "nearest",ratio=1,caliper=0.25)
#' ggPSMSummary(x,show.table=FALSE)
#' ggPSMSummary(x)
ggPSMSummary=function(x,show.table=TRUE,xpos=NULL,ypos=NULL){

     # show.table=TRUE;xpos=NULL;ypos=NULL
     # xvars=attr(x$model$terms,"term.labels")
     # yvar=names(x$model$model)[1]
     #
     # data1=x$model$data
     # data2=match.data(x)
     # summary(x)
     # summary(x)$sum.across
     # str(summary(x))

     res=summary(x)$sum.all[,3]
     if(is.null(summary(x)$sum.matched)){
         res2=summary(x)$sum.across[,3]
     } else{
         res2=summary(x)$sum.matched[,3]
     }

     res
     res2
     df1=data.frame(value=abs(res))
     df1$name="All Data"
     df2=data.frame(value=abs(res2))
     df2$name="Matched Data"
     df=rbind(df1,df2)
     rownames(df1)
     df$var=rep(rownames(df1),2)
     df
     df3=cbind(abs(res),abs(res2))
     colnames(df3)=c("All Data","Matched")
     df3
     df3=df3[order(df3[,1],decreasing=TRUE),]
     if(is.null(ypos)){
         maxy=max(res)
         miny=max(res2)
         ypos=miny+(maxy-miny)*1/2
     }
     if(is.null(xpos)) xpos=1.5
     table_grob=gridExtra::tableGrob(round(df3,3),theme=gridExtra::ttheme_minimal())
     p<-ggplot(data=df,aes_string(x="name",y="value",color="var",group="var",label="var"))+
         geom_point()+geom_line()+
         xlab("")+ylab("Absolute Standardized Diff in Means")+
         guides(colour=FALSE)+
         ggrepel::geom_text_repel(data=df[df$name=="All Data",],hjust=1.2)+
         ggrepel::geom_text_repel(data=df[df$name!="All Data",],hjust=-0.2)+
         theme_bw()
     if(show.table) p<-p+annotation_custom(grob=table_grob,xmin=xpos,ymin=ypos)
     p


}

#' Make p format
#' @param x string
#' @param digits numeric
#' @export
#'@examples
#'x=c("0.000","0.123","","0.34")
#'pformat(x)
pformat=function(x,digits=3){
     for(i in seq_along(x)){
         if(is.na(as.numeric(x[i]))) {
           next
         } else if(as.numeric(x[i])<10^(-digits)){
            x[i]=paste0("< 0.",paste0(rep("0",digits-1),collapse=""),"1")
         }
     }
     x
}

#'Make flextable summarizing propensity score matching
#' @param x An object of class matchit
#' @param digitsstd integer indicating the number of decimal places
#' @param grouplabel Optional group label
#' @export
#' @importFrom magrittr `%>%`
#' @importFrom officer fp_border
#' @importFrom moonBook mytable compress
#' @importFrom ztable ztable
#' @importFrom flextable flextable delete_part add_header_row hline_top hline align merge_at width autofit hline_bottom as_paragraph footnote
#' @examples
#' require(MatchIt)
#' formula=treat ~ age + educ + race+married+nodegree + re74 + re75
#' x=matchit(formula, data =lalonde, method= "full",link="probit")
#' PSMTable(x,grouplabel=c("Control","Treated"))
#' x=matchit(formula, data =lalonde, method= "nearest")
#' PSMTable(x,grouplabel=c("Control","Treated"))
PSMTable=function(x,digitsstd=3,grouplabel=NULL){

    # digitsstd=3;grouplabel=NULL
    temp1=formula2vars(x$formula)
    xvars=temp1$xvars
    yvar=temp1$yvar

    data1=x$model$data

    if(is.null(data1)) {
      result=call2param(x$call)
      data1=eval(parse(text=result$data))
    }
    data2=match.data(x)[1:ncol(data1)]
    data2

    data1$matchGroup="Before"
    data2$matchGroup="After"
    data3=rbind(data1,data2)


    temp=paste0("mytable(matchGroup+",yvar,"~",paste0(xvars,collapse="+"),",data=data3)")
    restable=eval(parse(text=temp))
    restable
    restable=compress(restable)
    res=list()
    form1=paste0("%.",digitsstd,"f")
    res[[1]]=sprintf(form1,summary(x)$sum.all[,3])

    if(is.null(summary(x)$sum.matched)){
      res[[2]]=sprintf(form1,summary(x)$sum.across[,3])
    } else{
      res[[2]]=sprintf(form1,summary(x)$sum.matched[,3])
    }
    if(rownames(summary(x)$sum.all)[1]=="distance"){
    res[[1]]=res[[1]][-1]
    res[[2]]=res[[2]][-1]
    }
    start=1
    restable[[1]]$res
    end=ncol(restable[[1]]$res)-7
    end

    px=restable[[1]]$res$ptest
    df=list()
    for( k in 1:2){
        df[[k]]=restable[[k]]$res[1:end]
        temp=c()
        j=1
        for(i in 1:length(px)){
            if(i==length(px)){
              temp=c(temp,res[[k]][j])
              j=j+1
            } else if((px[i]!="")&(px[i+1]=="")){
                temp=c(temp,"")
            } else{
                temp=c(temp,res[[k]][j])
                j=j+1
            }
        }
        df[[k]]$std=temp
    }

    df=cbind(df[[1]],df[[2]])
    df[[6]]<-""
    df[[4]]<-pformat(df[[4]])
    df[[9]]<-pformat(df[[9]])
    df
    weighted=FALSE
    if(length(unique(match.data(x)$weights))>1){

        result=getWeightedValues(x)

      df[[7]]=result$control
      df[[8]]=result$treat
      df[[9]]=result$p
      weighted=TRUE
    }

    subnames=colnames(df)

    colnames(df)=1:10

    header1=c("Covariates",paste0("N=",restable[[1]]$count),"","","",paste0("N=",restable[[2]]$count),"","")
    if(is.null(grouplabel)){
       name2=paste0(subnames[1],"=",subnames[2:3])
    } else{
       name2=grouplabel
    }
    header2=c("",name2,"p","standardized\ndifference","",name2,"p","standardized\ndifference")
    header3=c("Variables","Before Propensity Score Matching","","After Propensity Score Matching")

    big_border = fp_border(color="black", width = 2)
    small_border = fp_border(color="gray50", width = 1)
    ft<-flextable(df)
    ft
    ft<-ft %>% delete_part("header") %>%
        add_header_row(values=header1)%>%
        add_header_row(values=header2) %>%
        add_header_row(values=header3,colwidths=c(1,4,1,4)) %>%
        hline_top(part="header",border=big_border)%>%
        hline(i=2,j=2:3,border=small_border,part="header") %>%
        hline(i=2,j=7:8,border=small_border,part="header") %>%
        hline(i=1,j=2:5,border=small_border,part="header") %>%
        hline(i=1,j=7:10,border=small_border,part="header") %>%
        align(align="center",part="header") %>%
        align(align="right",part="body") %>%
        align(j=1,align="left",part="body") %>%
        merge_at(i=1:2,j=1,part="header") %>%
        merge_at(i=2:3,j=4,part="header") %>%
        merge_at(i=2:3,j=5,part="header") %>%
        merge_at(i=2:3,j=9,part="header") %>%
        merge_at(i=2:3,j=10,part="header") %>%
        hline_bottom(part="header",border=big_border)%>%
        hline_top(part="body",border=big_border)%>%
        hline(i=2,j=4:5,part="header",border=big_border)%>%
        hline(i=2,j=9:10,part="header",border=big_border)%>%
        autofit() %>%
        width(j=6,width=0.05) %>%
        width(j=1,width=0.5) %>%
        width(j=2:5,width=1.2) %>%
        width(j=7:10,width=1.2)
    if(weighted) ft <-footnote(ft,i=2,j=7:8,
                               ref_symbols=c("a"),
                               value=as_paragraph("Values are weighted mean \u00b1 weighted sd or weighted percentages"),
                               part="header",inline=T)
    ft
}

#' Make ggplot between covariates and propensity score
#' @param x An object of class matchit
#' @param xlabs labels for group variable
#' @param ncol Numeric Number of facet column
#' @importFrom tidyr pivot_longer
#' @importFrom MatchIt match.data matchit
#' @importFrom ggplot2 geom_smooth facet_wrap element_rect ggplot aes_string
#' @importFrom ggplot2 geom_point xlab ylab theme_bw theme
#' @export
#' @examples
#' require(MatchIt)
#' require(ggplot2)
#' formula=treat ~ age + race+educ + nodegree + re74 + re75
#' x=matchit(formula, data =lalonde, method= "nearest",ratio=1,caliper=0.25)
#' ggPS(x,xlabs=c("control","treated"))+theme(legend.position=c(0.8,0.15))
ggPS=function(x,xlabs=NULL,ncol=2){
      # xlabs=NULL;ncol=2
     temp1=formula2vars(x$formula)
     xvars=temp1$xvars
     yvar=temp1$yvar
     # xvars=attr(x$model$terms,"term.labels")
     # yvar=names(x$model$model)[1]
     # yvar
     dta_m=match.data(x)
     dta_m
     if(!is.factor(dta_m[[yvar]])) {
          dta_m[[yvar]]=factor(dta_m[[yvar]])
     }

     if(!is.null(xlabs)) levels(dta_m[[yvar]])=xlabs
    for(i in seq_along(xvars)){
        if(is.factor(dta_m[[xvars[i]]])) dta_m[[xvars[i]]]<-as.numeric(dta_m[[xvars[i]]])
    }

     longdata<-pivot_longer(dta_m,
                            cols=xvars,
                            names_to="variable")
     ggplot(longdata, aes_string(x = "distance", y = "value", colour = yvar)) +
          geom_point(alpha = 0.2, size = 1.3) +
          geom_smooth(method = "loess", se = F,formula="y~x") +
          facet_wrap("variable",ncol=ncol,scales="free")+
          xlab("Propensity score") +ylab("")+
          theme_bw() +
          theme(legend.background = element_rect(color="gray50"))
}

#'Extract parameter from call
#'@param call An object of class call
#'@importFrom stringr str_replace str_replace_all str_trim
#'@return a list of parameters
#'@export
call2param=function(call){
    # temp=deparse(call)[2]
    temp=paste(deparse(call),sep="",collapse="")
    if(is.na(temp)){
        result=list()
        result$method="nearest"
    } else{
        temp=str_replace(temp,"matchit\\(","")
        temp=str_replace_all(temp," ","")
        temp
        temp=unlist(strsplit(temp,","))
        temp=strsplit(temp,"=")
        result=lapply(temp,function(x){
            x=str_replace_all(x,'\\"|\\)',"")
            x=str_trim(x)
            x
        })
        res=list()
        for(i in 1:length(result)){
            res[[result[[i]][1]]]=result[[i]][2]

        }
        res2=lapply(res,function(x){
            ifelse(is.na(suppressWarnings(as.numeric(x))),x,as.numeric(x))
        })
        result=res2
        result
    }
    if(is.null(result$method)) result$method="nearest"
    if(is.null(result$ratio)) result$ratio=1
    if(is.null(result$caliper)) result$caliper=0.25
    if(is.null(result$link)) result$link="logit"
    if(is.null(result$distance)) result$distance="glm"
    if(is.null(result$estimand)) result$estimand="ATT"
    result
}


#' Extract variable names with formula
#' @param formula An object of class formula
#' @param allowInteraction logical
#' @export
#' @examples
#'formula = treat ~ age + educ + race + married+nodegree + re74 + re75
#'formula = treat ~ age+ race+ age:race
#'formula2vars(formula)
formula2vars=function(formula,allowInteraction=FALSE){
    temp=deparse(formula)
    temp=gsub(" ","",temp)
    temp=unlist(strsplit(temp,"~"))
    yvar=temp[1]
    yvar
    if(allowInteraction){
      xvars=unlist(strsplit(temp[2],"+",fixed=TRUE))
    } else{
       xvars=unique(unlist(strsplit(temp[2],"[+]|[*]|:")))
    }
    xvars
    list(yvar=yvar,xvars=xvars)
}

#' make pptList with an object of class matchit
#' @param x An object of class matchit
#' @param depvar Variable name serves as dependent variables
#' @param time  Name of time variable
#' @param status Name of status variable
#' @param seed Integer
#' @param m.threshold numeric The default value is 0.1
#' @param v.threshold numeric  The default value is 2
#' @param compare logical
#' @param report logical
#' @param multiple logical
#' @param depKind character One of c("continuous","binary","survival")
#' @param covarCentering logical
#' @param withinSubclass logical
#' @param analyzeSens logical
#' @importFrom lmtest coeftest
#' @importFrom sandwich vcovCL
#' @importFrom MatchIt matchit
#' @importFrom cobalt bal.plot bal.tab
#' @export
#' @examples
#' require(MatchIt)
#' x=matchit(treat ~ age + educ + race + married+nodegree + re74 + re75, data =lalonde,
#'    method="subclass",subclass=4)
#' x=matchit(treat ~ age + educ + race + married+nodegree + re74 + re75, data =lalonde,
#'    method="exact")
#' x=matchit(treat ~ age + educ + race + married+nodegree + re74 + re75, data =lalonde,
#'    method="full",link='probit')
#' x=matchit(treat ~ age + educ + race + married+nodegree + re74 + re75, data =lalonde,
#'    method="nearest",link='probit')
#' result=makePPTList_matchit(x)
#' result=makePPTList_matchit(x,depvar="re78",analyzeSens=TRUE)
#' data("GBSG2",package="TH.data")
#' x=matchit(horTh~age+menostat+tsize+tgrade+pnodes+progrec+estrec,data=GBSG2)
#' result=makePPTList_matchit(x,time="time",status="cens")
makePPTList_matchit=function(x,depvar=NULL,time="",status="",seed=1234,
                             m.threshold=0.1,v.threshold=2,
                             compare=TRUE,report=TRUE,
                             multiple=TRUE, depKind="continuous",
                             covarCentering=FALSE,withinSubclass=FALSE,
                             analyzeSens=FALSE){
     # depvar=NULL;time="time";status="cens";seed=1234;compare=TRUE;report=TRUE
     # multiple=TRUE; depKind="continuous"
     # covarCentering=FALSE;withinSubclass=FALSE
     # m.threshold=0.1;v.threshold=2

     `%!in%` <- Negate(`%in%`)

     if("character" %in% class(x)) {
         matched=eval(parse(text=x))
         matchedCall=x
     } else if("matchit" %in% class(x)){
         matched=x
         temp=paste0(deparse(x$call),collapse="")
         matchedCall=temp
    }
     result=call2param(matched$call)
     result
     matchMethod=result$method
     ratio=result$ratio
     caliper=result$caliper
     link=result$link
     distance=result$distance
     estimand=result$estimand
     temp=as.character(matched$call)
     # matchedCall=paste0(temp[1],"(",temp[2],",data=",temp[3],",method='",matchMethod,"'")
     # if(ratio!=1) matchedCall=paste0(matchedCall,",ratio=",ratio)
     # if(caliper!=0.25) matchedCall=paste0(matchedCall,",cailper=",caliper)
     # if(link!='logit') matchedCall=paste0(matchedCall,",link='",link,"'")
     # if(distance!='glm') matchedCall=paste0(matchedCall,",distance='",distance,"'")
     # if(estimand!='ATT') matchedCall=paste0(matchedCall,",estimand='",estimand,"'")
     # matchedCall=paste0(matchedCall,")")

     eq=temp[2]
     df=matched$model$data
     dfname=temp[3]

     temp1=formula2vars(matched$formula)
     xvars=temp1$xvars
     yvar=temp1$yvar

     count=length(xvars)

     title<-type<-code<-c()

     title="Set Seed Number"
     type="Rcode"
     code=paste0("set.seed(",seed,")")

     title=c(title,"Check Initial Imbalance")
     type=c(type,"Rcode")
     temp=paste0("out<-matchit(",paste0(deparse(matched$formula),collapse=""),",data=",dfname,",method=NULL,distance='glm')\nsummary(out)")
     code=c(code,temp)


     # title=c(title,"Chi-square test before matching")
     # type=c(type,"Rcode")
     # temp=paste0("myxBalance(",eq,",data=",dfname,",report='chisquare.test')")
     # code=c(code,temp)

     # if(!is.null(depvar)){
     # title=c(title,"Outcome Model using Regression Analysis")
     # type=c(type,"Rcode")
     # temp=paste0("fit1=lm(",depvar,"~",yvar,"+",paste0(xvars,collapse='+'),",data=",dfname,")\nsummary(fit1)")
     # code=c(code,temp)
     # }

     title=c(title,paste0("Matching using ",matchMethod))
     type=c(type,"Rcode")
     temp=paste0("matched <-",matchedCall,"\nmatched")
     code=c(code,temp)

     title=c(title,paste0("Summary of Matching using ",matchMethod))
     type=c(type,"Rcode")
     temp=paste0("summary(matched)")
     code=c(code,temp)

     if(matchMethod=="subclass"){
       title=c(title,paste0("Summary of Subclass"))
       type=c(type,"Rcode")
       temp=paste0("summary(matched,subclass=TRUE,un=FALSE)")
       code=c(code,temp)
     }

     if(matchMethod %!in% c("exact","cem")){
     title=c(title,"Covariates vs. Propensity Score")
     type=c(type,"ggplot")
     temp=paste0("ggPS(matched)")
     code=c(code,temp)
     }

     title=c(title,"Change of Absolute Standardised Differences")
     type=c(type,"ggplot")
     temp=paste0("ggPSMSummary(matched)")
     code=c(code,temp)

     if(matchMethod %!in% c("exact","cem")){
     title=c(title,"Distribution of Propensity Score")
     type=c(type,"plot")
     temp=paste0("plot(matched,type='jitter',col='blue',interactive=FALSE)")
     code=c(code,temp)

     title=c(title,"Histogram of Propensity Score")
     type=c(type,"plot")
     temp=paste0("plot(matched,type='hist',interactive=FALSE)")
     code=c(code,temp)
     }

     title=c(title,"Matched data")
     type=c(type,"Rcode")
     temp=paste0("match.data = match.data(matched)\nhead(match.data)")

     eval(parse(text=temp))

     code=c(code,temp)

     title=c(title,"Balance Table")
     type=c(type,"Rcode")
     code=c(code,paste0("options(crayon.enabled=FALSE);bal.tab(matched,thresholds=c(m=",m.threshold,",v=",v.threshold,"))"))

     if(matchMethod=="subclass"){
         title=c(title,"Love plot")
         type=c(type,"plot")
         temp=paste0("plot(summary(matched,subclass=TRUE),var.order='unmatched',abs=FALSE)")
         code=c(code,temp)
     } else{
        title=c(title,"Love plot","Love plot(variance)")
        type=c(type,"ggplot","ggplot")
          temp=paste0("love.plot(matched,threshold=",m.threshold,",sample.names=c('Unmatched','matched'),var.order='unadjusted',stars='raw')")
          code=c(code,temp)
          temp=paste0("love.plot(matched,stats='variance.ratios',threshold=",v.threshold,",sample.names=c('Unmatched','matched'),var.order='unadjusted')")
          code=c(code,temp)
     }


     title=c(title,"Summary of Propensity Score Matching")
     type=c(type,"table")
     temp=paste0("PSMTable(matched)")
     code=c(code,temp)

     if(matchMethod %!in% c("exact","cem")){
     title=c(title,"Summary of Propensity Score Matching")
     type=c(type,"ggplot")
     code=c(code,"cobalt::bal.plot(matched,var.name='distance',which='both',type='histogram',mirror=TRUE)")
     }

     if(compare & (matchMethod!="nearest")){
       title=c(title,"Compare Balance Table")
       type=c(type,"Rcode")
       code=c(code,"makeCompareBalTab(matched)")
       title=c(title,"compare Love Plot")
       type=c(type,"ggplot")
       code=c(code,"compareLove.plot(matched)")

     }

     if(report){
       title=c(title,"Report PS Matching")
       type=c(type,"html")
       if(is.null(depvar)){
         code=c(code,paste0("cat(reportPSM(matched,compare=",compare,"))"))
       } else{
         code=c(code,paste0("cat(reportPSM(matched,depvar='",depvar,"',compare=",compare,"))"))
       }

     }

     if(!is.null(depvar)){
       # str(depvar)
       title=c(title,"Estimating Treatment Effect")
       type=c(type,"Rcode")
       #      temp=paste0(
       # "fit1=lm(",depvar,"~",yvar,"+",paste0(xvars,collapse='+'),",data=match.data,weights=weights)\n",
       # "coeftest(fit1,vcov.=vcovCL,cluster=~subclass)")
       temp1=paste0("c('",paste0(depvar,collapse="','"),"')")
       temp=paste0("effect=estimateEffect(matched,mode='",depKind,"',multiple=",multiple,
                   ",dep=",temp1,",covarCentering=",covarCentering,",withinSubclass=",withinSubclass,",print=FALSE);effect")
       code=c(code,temp)

       if(report){
         title=c(title,"Report Treatment Effect")
         type=c(type,"html")
         code=c(code,"cat(attr(effect,'report'))")

       }


     if(analyzeSens){
     if(matchMethod %in% c("nearest","optimal","genetic")) {
       title=c(title,"Sensitivity Analysis")
       type=c(type,"Rcode")
       code=c(code,paste0("gammaRangeSearch(matched,dep='",depvar,"')"))
     } else if(matchMethod=="full"){
       title=c(title,"Sensitivity Analysis")
       type=c(type,"Rcode")
       code=c(code,paste0("gammaRangeSearchFull(matched,dep='",depvar,"')"))
     }
     }
     }
     if((time!="")&(status!="")){
       title=c(title,"Estimating Treatment Effect")
       type=c(type,"Rcode")
       #      temp=paste0(
       # "fit1=lm(",depvar,"~",yvar,"+",paste0(xvars,collapse='+'),",data=match.data,weights=weights)\n",
       # "coeftest(fit1,vcov.=vcovCL,cluster=~subclass)")
       temp=paste0("effect=estimateEffect(matched,mode='survival',multiple=",multiple,
                   ",time='",time,"',status='",status,"',covarCentering=",covarCentering,",withinSubclass=",withinSubclass,",print=FALSE);effect")
       code=c(code,temp)

       if(report){
         title=c(title,"Report Treatment Effect")
         type=c(type,"html")
         code=c(code,"cat(attr(effect,'report'))")

       }
     }

     data.frame(title,type,code,stringsAsFactors = FALSE)
}



