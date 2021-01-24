library(shinybusy)
library(shinyjs)

shinyUI(fluidPage(
    shinyjs::useShinyjs(),
    uiOutput("title"),
    add_busy_gif(src = "https://jeroen.github.io/images/banana.gif", height = 70, width = 70),
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
                         pptxListInput("List1"),icon=icon("shopping-cart")),



                id='main',
                theme=shinytheme("cerulean")
    ),
    uiOutput("showList"),
    verbatimTextOutput("table4")

)
)
