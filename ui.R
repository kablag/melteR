
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#
library(shiny)

shinyUI(fluidPage(

  # Application title
  titlePanel("melteR"),

  # Sidebar with a slider input for number of bins
  sidebarLayout(
    sidebarPanel(      
      uiOutput("ui.controls.source"),
      conditionalPanel(
        condition = "input['controls.source'] == 'standard.set'",
        uiOutput("ui.controls.source.standard")),
      conditionalPanel(
        condition = "input['controls.source'] == 'file'",
        uiOutput("ui.controls.source.file")),      
      fileInput("rdml.file", h4("Файл с неизвестными образцами"),
                accept=c("application/zip", ".rdml")),      
      sliderInput("sigma.threshold",
                  h4("Порог \u03C3"),
                  min = 0,
                  max = 1,
                  value = 0.1,
                  step = 0.01),
      h4("Обозначения образцов"),
      textInput("wt.name", 
                "Совпавшие с контролями", 
                "Чувствительный"),
      textInput("mut.name", 
                "Несовпавшие с контролями", 
                "Устойчивый"),      
      radioButtons("lang",
                   h4("Interface Language"),
                   c("Русский" = "ru",
                     "English" = "en")
      )
    ),

    # Show a plot of the generated distribution
    mainPanel(
      conditionalPanel(
        condition = "output.fileUploaded",
        tabsetPanel(
          tabPanel("Результат",                                      
                   dataTableOutput("short.result.table")),
          tabPanel("Результат подробно",
                   fluidRow(                     
                     column(6,                            
                            h5("Контрольные образцы"),                   
                            tableOutput("controls.table")),
                     column(6,
                            h5("Графики"),
#                             conditionalPanel(
#                               condition = "output.rowNotPlotted",
#                               helpText("Нажмите на реееед в таблице ниже длееее отображениеее")),
                            plotOutput("row.plot"))
                   ),
                   h5("Результат анализа"),
                   dataTableOutput("result.table")),
          tabPanel("Создать отчёт",
                   h5("Включить в отчёт"),
                   checkboxInput("rep.settings", "Параметры анализа", value = TRUE),
                   checkboxInput("rep.results.short", "Таблицу с результатами", value = TRUE),
                   checkboxInput("rep.controls", "Таблицу с контрол\u044Fми", value = TRUE),
                   checkboxInput("rep.results.long", "Таблицу с подробными результатами", value = TRUE),
                   checkboxInput("rep.plots", "Графики", value = TRUE),
                   conditionalPanel(
                     condition = "input['rep.plots']",
                     radioButtons("rep.plots.type",
                                  "Графики",
                                  c("Все" = "all",
                                    "Только дл\u044F образцов с диким типом" = "wt",
                                    "Только дл\u044F образцов с мутаци\u044Fми" = "mut"))
                   ),
                   radioButtons('format', 'Document format', c('PDF', 'HTML', 'Word'),
                                inline = TRUE),
                   downloadButton('downloadReport')
          )
        )
      )
    )
  )
))