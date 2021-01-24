library(shiny)
library(shinythemes)
library(webrSub)
library(rrtable)
library(moonBook)
library(ggplot2)
library(dplyr)
library(shinyFiles)
library(shinybusy)
library(shinyjs)
library(twang)

options(shiny.sanitize.errors = FALSE)

dic<-webrSub:::dic

cath=function(string="",size=3){
  cat(paste("<h",size,">",string,"</h",size,">",sep=""))
}

ui=fluidPage(
    shinyjs::useShinyjs(),
    add_busy_gif(src = "https://jeroen.github.io/images/banana.gif", height = 70, width = 70),
    uiOutput("title"),
    singleton(
        tags$head(tags$script(src = "message-handler.js"))
    ),
    radioButtons(inputId = "language", label = "Select Language",
                 choices = list("English" = "en", "한국어(Korean)" = "kor"),
                 selected = "en",inline=TRUE),
    navbarPage( "Web-R.org",
                tabPanel("DataSelect",
                         dataSelectInput("data"),
                         tableOutput("table3")
                         ,verbatimTextOutput("text")

                ),
                tabPanel("dataWrangling",
                         prepInput("pre")),
                tabPanel("Propensity Score",
                         PSMModuleInput("PSM")),
                tabPanel("PPTxList",
                         pptxListInput("List1")),


                id='main',
                theme=shinytheme("cerulean")
    )
    ,verbatimTextOutput("table4")

)

server=function(input,output,session){

    dataEx=c("lalonde","apm","acs", "radial","iris","dirty","band_members")

     dirty<<-rio::import("dirty_data.xlsx")

     dataname=reactiveValues(dataname="",preprocessing="")

     savedPPT=reactiveValues(type=c(),title=c(),code=c())

     langchoice=function(en,kor){
          ifelse(input$language=="en",en,kor)
     }

     langchoice1=function(id){

       temp=dic[dic$id==id,]
       ifelse(input$language=="en",temp$en,temp$kor)
     }

     apm<<-read.csv("data/apm.csv",stringsAsFactors = FALSE)

     data(lalonde,package="twang")

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
                       myp(langchoice("With this app, you can perform propensity score matching `without` R in your computer. ",
                                      "자신의 컴퓨터에 R을 설치할 필요 없이 R을 이용한 성향점수맞추기를 할 수 있읍니다. `자신의 데이타`를 xlsx 또는 csv형식으로 업로드할 수 있습니다.")),
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

     observeEvent(resPSM(),{
       savedPPT$type=resPSM()$df$type
       savedPPT$title=resPSM()$df$title
       savedPPT$code=resPSM()$df$code
       if(!is.null(resPSM()$exData)){
         if(resPSM()$exData!="None") dataname$dataname<-resPSM()$exData
       }
     })


     output$text=renderPrint({
         # cat("str(df())\n")
         # str(df())
         if(!is.null(result()$name)){
         if(result()$name!=""){
           temp=result()$name
            # cat("temp=",temp,"\n")
            assign(temp,df())
            temp1=paste0("str(",temp,")")
            cat(temp1,"\n")
            eval(parse(text=temp1))
         }}
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


}

shinyApp(ui=ui,server=server)
