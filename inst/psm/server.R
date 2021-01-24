library(shiny)
library(RColorBrewer)
library(rmarkdown)
library(gcookbook)
library(moonBook)
library(ggthemes)
library(ggfortify)
#library(twang)
library(MatchIt)
library(optmatch)
library(tidyr)
library(lmtest)
library(sandwich)
library(cobalt)
library(webrSub)
library(webrPSM)

options(rgl.useNULL = TRUE)



cath=function(string="",size=3){
    cat(paste("<h",size,">",string,"</h",size,">",sep=""))
}


shinyServer(function(input,output,session){

    dataEx=c("lalonde","apm","exData","swData","acs", "radial","iris","dirty","band_members")

    dirty<<-rio::import("dirty_data.xlsx")

    apm<<-read.csv("data/apm.csv",stringsAsFactors = FALSE)

    dataname=reactiveValues(dataname="",preprocessing="")

    savedPPT=reactiveValues(type=c(),title=c(),code=c())

    # data(colon)
    # data(lung,package="survival")
    # data(lalonde,package="twang")
    data(lalonde,package="MatchIt")

    langchoice=function(en,kor){
        ifelse(input$language=="en",en,kor)
    }

    langchoice1=function(id){

        temp=dic[dic$id==id,]
        ifelse(input$language=="en",temp$en,temp$kor)
    }

    my_theme=function(p,theme=""){
        if(theme==""){
          my_theme=input$theme
        } else{
          my_theme=theme
        }
        if(my_theme!="gray") {
              p<-eval(parse(text=paste("p+theme_",my_theme,"()",sep="")))
        }
      p
    }


    output$title=renderUI({

        tagList(
            h1(langchoice("Web-based Analysis with R 4.0","웹에서 하는 R 통계분석 4.0")),
            hr(),
            # if(input$main=="DataSelect")
            myp(langchoice("With this app, you can perform analysis `without` R in your computer. You can `analyze` data, make `tables` and `plots` and download the report as a `pdf`, `docx` or `pptx` file. You can also download the high-quality plots with desired size and resolution.","자신의 컴퓨터에 R을 설치할 필요 없이 R을 이용한 통계분석을 할 수 있읍니다. 그룹변수와 행변수를 선택하여 쉽게 표를 만들 수 있으며 그래프를 통한 자료 탐색과 여러가지 통계분석이 가능합니다. `자신의 데이타`를 xlsx 또는 csv형식으로 업로드하여 분석을 할 수 있을 뿐 아니라 그 결과를 `pdf`,`docx`,`powerpoint` 파일로 다운로드할 수 있읍니다. 또한 Plot을 원하는 크기로 저장할 수 있읍니다. `표가 보일 때까지 잠시만` 기다려주세요.")),
            # if(input$main=="DataSelect")
            hr()
        )
    })
    result=callModule(dataSelect,"data",
                      dataEx=reactive(dataEx),
                      lang=reactive(input$language),
                      dataname=reactive(dataname$dataname),
                      preprocessing0=reactive(dataname$preprocessing))

    df=callModule(prep,"pre",dataname=reactive(result()$name),
                  preprocessing=reactive(result()$preprocessing),
                  lang=reactive(input$language))

    origin=reactive({df()})

    resPSM=callModule(PSMModule,"PSM",dataname=reactive(result()$name),
                       df=reactive(df()),
                       preprocessing=reactive(result()$preprocessing),
                       PPTdata=reactive(PPTdata()),
                       lang=reactive(input$language))


    pptdf=reactive({
        if(length(savedPPT$code)==0) {
            result=""
        } else{

            result<-data.frame(type=savedPPT$type,title=savedPPT$title,code=savedPPT$code,
                               stringsAsFactors = FALSE)
        }
        result
    })

    PPTdata=reactive({
        data.frame(type=savedPPT$type,title=savedPPT$title,code=savedPPT$code,
                   stringsAsFactors = FALSE)
    })

    pptdf2=callModule(pptxList,"List1",data=reactive(pptdf()),
                      preprocessing=reactive(dataname$preprocessing))

    output$text=renderPrint({
        # cat("str(df())\n")
        # str(df())
        cat("dataname$dataname=",dataname$dataname,"\n")
        if(!is.null(result()$name)){
            if(result()$name!=""){
                temp=result()$name
                dataname$dataname<-temp
                # cat("temp=",temp,"\n")
                assign(temp,df())
                temp1=paste0("str(",temp,")")
                cat(temp1,"\n")
                eval(parse(text=temp1))
            }}
        cat("\ndataname$dataname=",dataname$dataname,"\n")
    })

    output$table3=renderTable({head(df(),10)})

    output$table4=renderPrint({
        result<-NULL
        try(result<-pptdf2())
        if(!is.null(result)){
            # cat("str(pptdf2())\n")
            # str(result)
            # cat("is.null(pptdf2()$code)\n")
            # is.null(result$code)
            savedPPT$type=result$type
            savedPPT$title=result$title
            savedPPT$code=result$code
        }

    })


    observeEvent(resPSM(),{
      savedPPT$type=resPSM()$df$type
      savedPPT$title=resPSM()$df$title
      savedPPT$code=resPSM()$df$code
      if(!is.null(resPSM()$exData)){
        if(resPSM()$exData!="None") dataname$dataname<-resPSM()$exData
      }
    })

    output$PPTxListTable1=renderUI({
      df2flextable(PPTdata()) %>% autofit() %>% htmltools_value()
    })
    output$showList=renderUI({
      data=PPTdata()

      if(input$main!="PPTxList"){
      tagList(
        hr(),
        checkboxInput("showPPTxList","show PPTxList"),
        conditionalPanel(sprintf("input['%s']==true","showPPTxList"),
                         if(nrow(data)>0) uiOutput("PPTxListTable1"),
                         if(nrow(data)==0) p("There is no saved data")
        )
      )
      }
    })

})
