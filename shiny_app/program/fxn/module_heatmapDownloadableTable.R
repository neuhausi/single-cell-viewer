heatmap_brks <- seq(0, 100, by = 1)
heatmap_clrs <- colorRampPalette(c("deepskyblue3", "white", "firebrick3"), space = "rgb")(length(heatmap_brks) + 1)


heatmap_downloadableTableUI <- function(id, downloadtypes = c("csv"), hovertext = NULL, contentHeight = "200px",
                                         singleSelect = FALSE) {
    ns <- shiny::NS(id)
    list(shiny::span(id = ns("dtableButtonDiv"), class = "periscope-downloadable-table-button",
                     style = "display:none", downloadFileButton(ns("dtableButtonID"),
                                                                downloadtypes, hovertext)), DT::dataTableOutput(ns("dtableOutputID")),
         shiny::tags$input(id = ns("dtableOutputHeight"), type = "text",
                           class = "shiny-input-container hidden", value = contentHeight),
         shiny::tags$input(id = ns("dtableSingleSelect"), type = "text",
                           class = "shiny-input-container hidden", value = singleSelect))
}


heatmap_downloadableTable <- function(input, output, session, logger, filenameroot,
                                       downloaddatafxns = list(), tabledata, rownames = TRUE, caption = NULL) {
    shiny::callModule(downloadFile, "dtableButtonID", logger,
                      filenameroot, downloaddatafxns)
    dtInfo <- shiny::reactiveValues(selected = NULL, tabledata = NULL,
                                    downloaddatafxns = NULL)
    shiny::observe({
        dtInfo$selected <- input$dtableOutputID_rows_selected
    })
    shiny::observe({
        dtInfo$tabledata <- tabledata()
    })
    shiny::observe({
        dtInfo$downloaddatafxns <- lapply(downloaddatafxns, do.call,
                                          list())
        rowct <- lapply(dtInfo$downloaddatafxns, nrow)
        session$sendCustomMessage("downloadbutton_toggle", message = list(btn = session$ns("dtableButtonDiv"),
                                                                          rows = sum(unlist(rowct))))
    })
    output$dtableOutputID <- DT::renderDataTable({
        sourcedata <- dtInfo$tabledata
        if (!is.null(sourcedata) && nrow(sourcedata) > 0) {
            DT_RowId <- paste0("rowid_", seq(1:nrow(sourcedata)))
            sourcedata <- DT::datatable(cbind(DT_RowId, sourcedata),
                                        options   = list(deferRender     = FALSE,
                                                         scrollX         = TRUE,
                                                         scrollY         = input$dtableOutputHeight,
                                                         paging          = FALSE,
                                                         dom             = "<\"periscope-downloadable-table-header\"f>tr",
                                                         processing      = TRUE,
                                                         rowId           = 1,
                                                         columnDefs      = list(list(targets = 0, visible = FALSE, searchable = FALSE)),
                                                         searchHighlight = TRUE),
                                        class     = paste("periscope-downloadable-table table-condensed", "table-striped table-responsive"),
                                        rownames  = rownames,
                                        selection = ifelse(input$dtableSingleSelect == "TRUE", "single", "multi"),
                                        caption   = caption, escape = FALSE,
                                        style     = "bootstrap") %>%
                DT::formatStyle(c(6, 7), backgroundColor = DT::styleInterval(heatmap_brks, heatmap_clrs))
        }
        sourcedata
    })

    selectedrows <- shiny::reactive({return(shiny::isolate(dtInfo$tabledata)[dtInfo$selected, ])})
    return(selectedrows)
}
