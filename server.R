
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



shinyServer(function(input, output) {
  vals <- reactiveValues()
  
  if("NormMeltData" %in% names(RDML$public_methods) == FALSE) {
    RDML$set("public", "NormMeltData", function() {  
      library(qpcR)  
      temps <- rownames(private$.melt.fdata)
      out <- c()
      n.fdata <- ncol(private$.melt.fdata)
      progress.step <- 1/n.fdata
      withProgress(message = "Обработка образцов",
                   value = 0, {
                     for(i in 1:n.fdata) {
                       data <- cbind(as.numeric(temps), as.numeric(private$.melt.fdata[, i]))
                       mcurve<-invisible(meltcurve(data,
                                                   temps=1,
                                                   fluo=2,
                                                   norm=TRUE,
                                                   plot=FALSE))    
                       mcurve.norm <- qpcR:::rescale(mcurve[[1]]$df.dT, tomin=0, tomax=1)
                       out <- cbind(out, mcurve.norm)
                       rownames(out) <- mcurve[[1]]$Temp
                       incProgress(progress.step)
                     }
                   })      
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
                      unkn$sigma <- NA
                      unkn$substracted <- NA
                      for(i in 1:length(vals$controls())) {                        
                        if(unkn$Target == vals$controls()[[i]]$Target){
                          substr <- unkn$Melt - vals$controls()[[i]]$Melt                          
                          substr.model <- lm(substr ~ as.numeric(names(unkn$Melt)))        
                          sigma <- summary(substr.model)$sigma
                          if(is.na(unkn$sigma)) {                            
                            if(sigma < input$sigma.threshold){                              
                              unkn$result <- vals$controls()[[i]]$FDataName 
                            }
                            unkn$result.i <- i
                            unkn$sigma <- sigma
                            unkn$substracted <- substr
                          }                                    
                          else {
                            if(sigma < unkn$sigma) {
                              unkn$sigma <- sigma
                              unkn$substracted <- substr
                              if(sigma < input$sigma.threshold){
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
        round(unkn$sigma, digits = 2),        
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
                                            "sigma",
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
    plot(as.numeric(names(vals$controls()[[control.i]]$Melt)),
         vals$controls()[[control.i]]$Melt,
         type = "l",
         ylim = c(1, -1), #invert y-axe
         xlim = c(40, 80),
         col = "green",
         lwd = 2,
         xlab = "T",
         ylab = "norm(d(RFU)/dT)")    
    lines(as.numeric(names(vals$results()[[sample.name]]$Melt)),
          vals$results()[[sample.name]]$Melt,                      
          col = "blue",
          lwd = 2)
    lines(as.numeric(names(vals$results()[[sample.name]]$Melt)),
          vals$results()[[sample.name]]$substracted,          
          col = "red",
          lwd = 2)
    title.els <- c(sample.name,
                   "=",
                         round(vals$results()[[sample.name]]$sigma,
                         digits = 2))
    title(bquote(.(title.els[1]) ~~~ hat(sigma) ~ .(title.els[2]) ~ .(title.els[3])))
    legend(70, -1, legend = c(vals$controls()[[control.i]]$TubeName,
                              vals$results()[[sample.name]]$TubeName,
                              paste(vals$results()[[sample.name]]$TubeName,
                                    "-",
                                    vals$controls()[[control.i]]$TubeName,
                                    sep = " ")),
           col = c("green", "blue", "red"),
           bty = "n",
           lty = 1,
           lwd = 2)
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
