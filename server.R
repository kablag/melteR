
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)

if(!("RDML" %in% rownames(installed.packages()))) {
  library(devtools)
  install_github("kablag/RDML")  
}

library(RDML)
library(MBmca) 
library(plotrix)
library(ggplot2)

shinyServer(function(input, output) {
  vals <- reactiveValues()
  
  if("NormMeltData" %in% names(RDML$public_methods) == FALSE) {
    RDML$set("public", "NormMeltData", function() {  
       
      temps <- as.numeric(rownames(private$.melt.fdata))
      out <- c()
      n.fdata <- ncol(private$.melt.fdata)
      progress.step <- 1/n.fdata
      withProgress(message = "Обработка образцов",
                   value = 0, {
                     for(i in 1:n.fdata) {
#                        data <- cbind(temps, as.numeric(private$.melt.fdata[, i]))
#                        mcurve<-invisible(meltcurve(data,
#                                                    temps=1,
#                                                    fluo=2,
#                                                    norm=TRUE,
#                                                    plot=FALSE))    
#                        mcurve.norm <- qpcR:::rescale(mcurve[[1]]$df.dT, tomin=0, tomax=1)                       
                       mcurve.norm <- mcaSmoother(temps, private$.melt.fdata[, i], 
                                                  minmax = TRUE)
                       mcurve.diff <- suppressMessages(diffQ(mcurve.norm,
                                            verbose = TRUE)$xy[["d(F) / dT"]])
#                      
                       out <- cbind(out, mcurve.diff)                       
                       #                        rownames(out) <- mcurve[[1]]$Temp
                       incProgress(progress.step)
                     }
                   })
      rownames(out) <- temps[2:(length(temps) - 1)]
      private$.melt.fdata <- out  
    })
  }
  
  rdml.obj <- reactive({
    if(is.null(input$rdml.file))      
      return(NULL)    
#     cat("loading data\n")    
    isolate({
      withProgress(message = "Загрузка RDML данных",
                   value = 0, {
                     rdml.obj <- RDML$new(input$rdml.file$datapath,
                                          name.pattern = "%TUBE% %NAME% %TARGET%")                     
                   })     
      
    })    
    rdml.obj$NormMeltData()                 
    return(rdml.obj)
  })
  output$fileUploaded <- reactive({
    return(!is.null(rdml.obj()))
  })
  outputOptions(output, 'fileUploaded', suspendWhenHidden=FALSE)
  
  vals$controls <- reactive({
    if(is.null(rdml.obj()))
      return()
#     cat("process controls\n") 
    
    # get controls from standard set
    if(input$controls.source == "standard.set") {      
      isolate({
        withProgress(message = "Загрузка RDML данных для контролей",
                     value = 0, {
                       file <- paste0("./positive_controls/",
                                      input$controls.source.standard)
                       controls.obj <- RDML$new(file,
                                            name.pattern = "%TUBE% %NAME% %TARGET%")                     
                     })    
        })
      controls.obj$NormMeltData()
      return(
        controls.obj$GetFData(
          filter = list(method = "melt",
                        types = "pos"),
          as.table = FALSE)
        )
    }
    
    # get controls from other file
    if(input$controls.source == "file") {
      if(is.null(input$rdml.controls.file))      
        return(NULL)
      isolate({
        withProgress(message = "Загрузка RDML данных для контролей",
                     value = 0, {
                       file <- input$rdml.controls.file$datapath                       
                       controls.obj <- RDML$new(file,
                                                name.pattern = "%TUBE% %NAME% %TARGET%")                     
                     })    
      })
      controls.obj$NormMeltData()
      return(
        controls.obj$GetFData(
          filter = list(method = "melt",
                        types = "pos"),
          as.table = FALSE)
      )
    }
    
    # get controls from file with unknowns
    rdml.obj()$GetFData(
      filter = list(method = "melt",
                    types = "pos"),
      as.table = FALSE)
  })
  
  output$controls.table <- renderTable({
    if(is.null(vals$controls()))
      return()
#     cat("generating controls table\n")
    controls.table <- lapply(vals$controls(), function(control) {
      c(control$Tube,
        control$TubeName,
        control$Dye,
        control$Target)
    })
    controls.table <- matrix(unlist(controls.table),
                             ncol = 4,
                             byrow = TRUE,
                             dimnames = list(c(),
                                             c("Пробирка",
                                               "Название",
                                               "Канал",
                                               "Мишень")))
    vals$controls.table <- controls.table
    return(controls.table)
  })
  outputOptions(output, 'controls.table', suspendWhenHidden=FALSE)
  
  vals$results <- reactive({
    if(is.null(vals$controls()))
      return()
#     cat("calc results\n")
    unkns <- rdml.obj()$GetFData(
      filter = list(method = "melt",
                    types = "unkn"),
      as.table = FALSE)
    unkns.names <- c()
    unkns <- lapply(unkns,
                    function(unkn) {
                      unkns.names <<- c(unkns.names,
                                       unkn$FDataName)
                      unkn$result <- input$mut.name
                      unkn$result.i <- NA
                      unkn$cor.r <- NA
                      unkn$substracted <- NA
                      for(i in 1:length(vals$controls())) {                        
                        if(unkn$Target == vals$controls()[[i]]$Target){
                          substr <- unkn$Melt - vals$controls()[[i]]$Melt                          
                          substr.model <- lm(substr ~ as.numeric(names(unkn$Melt)))        
                         cor.r <- cor(unkn$Melt,
                                       vals$controls()[[i]]$Melt)
                         
                          if(is.na(unkn$cor.r)) {                            
                            if(cor.r > input$cor.r.threshold){                              
                              unkn$result <- vals$controls()[[i]]$FDataName 
                            }
                            unkn$result.i <- i
                            unkn$cor.r <- cor.r
                            unkn$substracted <- substr
                          }                                    
                          else {
#                             if(cor.r < unkn$cor.r) {
                            if(cor.r > unkn$cor.r) {
                              unkn$cor.r <- cor.r
                              unkn$substracted <- substr
                              if(cor.r > input$cor.r.threshold){
                                unkn$result.i <- i
                                unkn$result <- vals$controls()[[i]]$FDataName
                              }
                            }
                          }
                        }
                        
                      }                                                                  
                    return(unkn)                    
                    })
    names(unkns) <- unkns.names    
    return(unkns)
  })
  
  res.tbl <- reactive({
    if(is.null(vals$results()))
      return()
#     cat("generating result table\n")    
    unkns.table <- lapply(vals$results(), function(unkn) {
      c(unkn$Tube,
        unkn$TubeName,
        unkn$Dye,
        unkn$Target,
        unkn$result.i,
        round(unkn$cor.r, digits = 2),        
        unkn$result
      )
    })    
    unkns.table <- matrix(unlist(unkns.table),
                          ncol = 7,
                          byrow = TRUE,
                          dimnames = list(c(),
                                          c("Пробирка",
                                            "Название",
                                            "Канал",
                                            "Мишень",
                                            "Контроль",
                                            "r",
                                            "Результат")))
    vals$results.table <- unkns.table
    return(unkns.table)
  })
  
  output$result.table <- renderDataTable(
    res.tbl(),
    callback = "function(table) {
    table.on('click.dt', 'tr', function() {
    $(this).toggleClass('selected').siblings().removeClass('selected');
    Shiny.onInputChange('selected.row',
    table.rows('.selected').data().toArray());
    });
    }"
  )
  outputOptions(output, 'result.table', suspendWhenHidden=FALSE)
  
  output$row.plot <- renderPlot({    
    if(length(input$selected.row) != 7)
      return(NULL)    
    control.i <- as.integer(input$selected.row[5])
    sample.name <- paste(input$selected.row[1],
                         input$selected.row[2],
                         input$selected.row[4],
                         sep = " ")
#     to.min.max.fluor <- c(vals$controls()[[control.i]]$Melt,
#                           vals$results()[[sample.name]]$Melt,
#                           vals$results()[[sample.name]]$substracted)
#     min.f <- min(to.min.max.fluor)
#     max.f <- max(to.min.max.fluor)
    to.min.max.t <- as.numeric(names(vals$controls()[[control.i]]$Melt))
    min.t <- min(to.min.max.t)
    max.t <- max(to.min.max.t)
    scaled.result <- rescale(vals$results()[[sample.name]]$Melt,
                             c(min.t, max.t))
    
    tmp <- data.frame(temp = c(to.min.max.t,
                               as.numeric(names(vals$results()[[sample.name]]$Melt))),
                               #to.min.max.t),
#                                scaled.result),
                      fluor = c(vals$controls()[[control.i]]$Melt * -1,
                                vals$results()[[sample.name]]$Melt * -1),
                                #vals$results()[[sample.name]]$substracted * -1),
#                                 vals$controls()[[control.i]]$Melt * -1),
                      name = rep(c(sample.name,
                                   vals$controls()[[control.i]]$FDataName),
                                 each = length(to.min.max.t)))
if(input$show.cor) {
  tmp <- do.call("rbind",list(tmp,
               data.frame(temp = scaled.result,
                          fluor = vals$controls()[[control.i]]$Melt * -1,
                          name = rep("Кор.", length(to.min.max.t)))))
}
    
    tmplm <- data.frame(control = vals$controls()[[control.i]]$Melt * -1,
                        sample = scaled.result)
    model.y <- lm(control ~ sample, data = tmplm)
    coef.y <- coef(model.y)
    p <- qplot(temp, fluor, data = tmp, geom = "path", color = name,
          main = paste0(sample.name,
                       "; r = ",
                       round(vals$results()[[sample.name]]$cor.r,
                             digits = 2)
                       ),
          ylab = "RFU",
          xlab = "Температура, °C") +
  theme(legend.position="bottom") +
  theme(legend.title=element_blank())
if(input$show.cor){
p <- p +
  geom_abline(intercept=coef.y[1],
              slope=coef.y[2])
}
p

  })
#   output$rowNotPlotted <- reactive({
# #     print(output$row.plot())
#     return(is.null(output$row.plot))
#   })
#   outputOptions(output, 'rowNotPlotted', suspendWhenHidden=FALSE)
  
  output$short.result.table <- renderDataTable({
    if(is.null(vals$results()))
      return()
#     cat("generating short result table\n")
    unkn.unique.names <- unique(
      sapply(vals$results(), function(unkn) {
        unkn$TubeName
      }))
    targets <- rdml.obj()$targets
    result.table <- matrix(nrow = length(unkn.unique.names),
                           ncol = length(targets) + 1,
                           dimnames = list(unkn.unique.names,
                                           c("Название",
                                             targets)))
    for(unkn in vals$results()) {
      result.table[unkn$TubeName, "Название"] <- unkn$TubeName
      result <- ifelse(unkn$result == input$mut.name,
                       input$mut.name,
                       input$wt.name)
      # check if result for this sample already exists at atble
      # (two or more replications of sample)
      if(!is.na(result.table[unkn$TubeName, unkn$Target])) {
        if(result.table[unkn$TubeName, unkn$Target] != result) {
          result.table[unkn$TubeName, unkn$Target] <- 
            paste0(result.table[unkn$TubeName, unkn$Target], " !")
        }
      }
      else {
        result.table[unkn$TubeName, unkn$Target] <- result
      }
    }
    vals$short.result.table <- result.table
    return(result.table)
  })
  
output$downloadReport <- downloadHandler(
  filename = function() {
    paste('my-report', sep = '.', switch(
      input$format, PDF = 'pdf', HTML = 'html', Word = 'docx'
    ))
  },
  
  content = function(file) {
    report.template <- "report_en.Rnw"
    src <- normalizePath(report.template)
    
    # temporarily switch to the temp dir, in case you do not have write
    # permission to the current working directory
    owd <- setwd(tempdir())
    on.exit(setwd(owd))
    file.copy(src, report.template, overwrite = TRUE)
    
#     library(rmarkdown)
#     out <- render('test.Rmd', switch(
#       input$format,
#       PDF = pdf_document(), HTML = html_document(), Word = word_document(),
#       encoding = "UTF-8"
#     ))
    library(knitr)
#     kn <- Sweavy2knitr("report_ru.Rnw")
    out<- knit2pdf(report.template, texi2dvi="pdflatex")
    #out <- texi2dvi("report_ru.tex", pdf = TRUE)
    file.rename(out, file)
  }
)
})
