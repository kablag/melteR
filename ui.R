
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#
library(shiny)

files <- list.files("./positive_controls")
names(files) <- files

shinyUI(fluidPage(

  # Application title
  titlePanel("melteR"),
  sidebarLayout(
    sidebarPanel(      
      radioButtons("controls.source",
                   h4("Источник контрольных образцов"),
                   c("Стандартный набор" = "standard.set",
                     "Из внешнего файла" = "file",
                     "Из файла с неизвестными образцами" = "inner")),
      conditionalPanel(
        condition = "input['controls.source'] == 'standard.set'",
        selectInput("controls.source.standard",
                    "Стандартный набор контролей:",
                    files)),
      conditionalPanel(
        condition = "input['controls.source'] == 'file'",
        fileInput("rdml.controls.file", "Файл с контрольными образцами",
                  accept=c("application/zip", ".rdml"))),      
      fileInput("rdml.file", h4("Файл с неизвестными образцами"),
                accept=c("application/zip", ".rdml")),      
      sliderInput("cor.r.threshold",
                  h4("Порог r"),
                  min = 0.8,
                  max = 1,
                  value = 0.98,
                  step = 0.01),
      h4("Обозначения образцов"),
      textInput("wt.name", 
                "Совпавшие с контролями", 
                "Чувствительный"),
      textInput("mut.name", 
                "Несовпавшие с контролями", 
                "Устойчивый")      
#       radioButtons("lang",
#                    h4("Interface Language"),
#                    c("Русский" = "ru",
#                      "English" = "en")
#       )
    ),
    
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
                            checkboxInput("show.cor", "Отображать корреляцию",
                                          value = TRUE),
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
                   checkboxInput("rep.controls", "Таблицу с контролями", value = TRUE),
                   checkboxInput("rep.results.long", "Таблицу с подробными результатами", value = TRUE),
                   checkboxInput("rep.plots", "Графики", value = TRUE),
                   conditionalPanel(
                     condition = "input['rep.plots']",
                     radioButtons("rep.plots.type",
                                  "Графики",
                                  c("Все" = "all",
                                    "Только для образцов с диким типом" = "wt",
                                    "Только для образцов с мутациями" = "mut"))
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