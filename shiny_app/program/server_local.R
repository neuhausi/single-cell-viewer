# ----------------------------------------
# --       PROGRAM server_local.R       --
# ----------------------------------------
# USE: Session-specific variables and
#      functions for the main reactive
#      shiny server functionality.  All
#      code in this file will be put into
#      the framework inside the call to
#      shinyServer(function(input, output, session)
#      in server.R
#
# NOTEs:
#   - All variables/functions here are
#     SESSION scoped and are ONLY
#     available to a single session and
#     not to the UI
#
#   - For globally scoped session items
#     put var/fxns in server_global.R
#
# FRAMEWORK VARIABLES
#     input, output, session - Shiny
#     ss_userAction.Log - Reactive Logger S4 object
# ----------------------------------------

# -- IMPORTS --


# -- VARIABLES --
userData <- reactiveValues()
    userData$object              <- NULL
    userData$objectId            <- 0
    userData$modalOpen           <- FALSE
    userData$filterList          <- NULL
    userData$filteredObject      <- NULL

userSel <- reactiveValues()
    userSel$genes                <- NULL
    userSel$scatterCluster       <- "All"
    userSel$scatterPanelColumns  <- NULL
    userSel$runScatterPlot       <- FALSE
    userSel$violinPanelColumns   <- NULL
    userSel$runViolinPlot        <- FALSE
    userSel$addDotGenes          <- NULL
    userSel$runDotPlot           <- FALSE
    userSel$addHeatmapGenes      <- NULL
    userSel$runHeatmap           <- FALSE
    userSel$diffCluster1         <- NULL
    userSel$diffCluster2         <- NULL
    userSel$diffClusterCells1    <- NULL
    userSel$diffClusterCells2    <- NULL
    userSel$runDiffCalc          <- FALSE


# -- FUNCTIONS --

disable_plots <- function() {
    userSel$runScatterPlot <- FALSE
    userSel$runViolinPlot  <- FALSE
    userSel$runDotPlot     <- FALSE
    userSel$runHeatmap     <- FALSE
    userSel$runDiffCalc    <- FALSE
}

filtered_data <- reactive({
    result  <- userData$object
    filters <- list()

    if (!is.null(result)) {

        # Check if filters available in input and that values have changed since last time
        input_names        <- names(input)
        input_filters      <- get_filter_input_Fields(input_names, isolate(userData$objectId))

        if (!is.null(result$metadata) && !identical(input_filters, character(0))) {

            changed_fields   <- NULL
            data_filtering   <- FALSE
            # detect if there has been a change since the last time
            for (input_filter in input_filters) {
                filter_values <- input[[input_filter]]
                filter_name   <- gsub(get_filter_prefix(isolate(userData$objectId)), "", input_filter)
                if (is.null(filter_values) || !identical(isolate(userData$filterList)[[filter_name]], filter_values)) {
                    data_filtering <- TRUE
                    changed_fields <- c(changed_fields, filter_name)
                }
                ({userData$filterList[[filter_name]] <- filter_values})
            }
            # Filter the data
            if (data_filtering) {
                metadata         <- result$metadata %>%
                    mutate(id = seq(1, nrow(result$metadata)))
                for (field in changed_fields) {
                    filter_values <- isolate(userData$filterList)[[field]]
                    if (is.null(filter_values)) {
                        metadata <- metadata[0,]
                        break
                    }
                    metadata <- metadata %>%
                        filter(as.vector(!!(rlang::sym(field))) %in% filter_values)
                }
                if (nrow(metadata) > 0) {
                    filtered_rowids <- metadata %>% pull(id)
                } else {
                    filtered_rowids <- c()
                }
                result$tsne       <- get_filtered_cell_data(result$tsne, row_ids = filtered_rowids)
                result$clusters   <- get_filtered_cell_data(result$clusters, row_ids = filtered_rowids)
                result$cells      <- result$cells[filtered_rowids]

                # filter expression and detection data
                seurat_object <- result$seurat
                if (!is.null(filtered_rowids)) {
                    assay_use     <- seurat_object@active.assay
                    
                    seurat_object@assays[[assay_use]]@data <- seurat_object@assays[[assay_use]]@data[, filtered_rowids]
                    seurat_object@active.ident              <- droplevels(seurat_object@active.ident[filtered_rowids])
                    result$expression                       <- future({avg.ex.scale(seurat_object) %>% as.data.frame()}, stdout = FALSE)
                    result$detection                        <- future({local_AverageDetectionRate(seurat_object) %>% as.data.frame()}, stdout = FALSE)
                    
                    #cleanup
                    rm(seurat_object)
                } else {
                    result$expression <- future({NULL})
                    result$detection  <- future({NULL})
                }

                # save filtered data and disable plots
                isolate(userData$filteredObject <- result)
                disable_plots()
            }
            # save filtered data
            if (!is.null(isolate(userData$filteredObject))) {
                result <- isolate(userData$filteredObject)
            }
        }
    }
    result
})

top10_genes <- reactive({
    get_top_genes_data(userData$object, "top10")
})

top30_genes <- reactive({
    get_top_genes_data(userData$object, "top30")
})

differentials_table_content <- reactive({
    c1 <- userSel$diffCluster1
    s1 <- userSel$diffClusterCells1
    c2 <- userSel$diffCluster2
    s2 <- userSel$diffClusterCells2

    result <- NULL
    if (userSel$runDiffCalc) {
        result <- calculate_differentials(filtered_data(), c1, s1, c2, s2)
        if (is.null(result)) {
            createAlert(session,
                        "bodyAlert",
                        "zeroDiffOutputAlertID",
                        style   = "warning",
                        content = paste("Differential Analysis resulted in no output for the current selections and logFC > ",
                                        g_differential_logfc_threshold),
                        append  = FALSE)
        }
    }
    result
})


# plot name reactives

base_filename <- reactive({
    current_time <- format(Sys.time(), "%Y.%m.%d_%H.%M")
    paste0(current_time, "_SCV")
})

base_plot_filename <- reactive({
    paste0(base_filename(), "_", userData$object$meta$object_name)
})

overview_plot_filename <- reactive({
    paste0(base_plot_filename(), "_Overview")
})

scatter_plot_filename <- reactive({
    paste0(base_plot_filename(), "_Scatter")
})

violin_plot_filename <- reactive({
    paste0(base_plot_filename(), "_Violin")
})

dot_plot_filename <- reactive({
    paste0(base_plot_filename(), "_Dot")
})

heatmap_filename <- reactive({
    paste0(base_plot_filename(), "_Heatmap")
})

differentials1_plot_filename <- reactive({
    paste0(base_plot_filename(), "_Differentials1")
})

differentials2_plot_filename <- reactive({
    paste0(base_plot_filename(), "_Differentials2")
})

# download filenames for tables

top10_DE_filename <- reactive({
    get_top_DE_download_filename(userData$object$meta$object_name, "_Top10DE_Genes")
})

top30_DE_filename <- reactive({
    get_top_DE_download_filename(userData$object$meta$object_name, "_Top30DE_Genes")
})

differentials_filename <- reactive({
    get_differentials_filename(userData$object$meta$object_name)
})

# -- MODULES --
callModule(downloadableTable, "top10DE", ss_userAction.Log,
           filenameroot     = top10_DE_filename,
           downloaddatafxns = list(csv = top10_genes,
                                   tsv = top10_genes),
           tabledata        = top10_genes,
           rownames         = FALSE)

callModule(downloadFile,     "top30DE", ss_userAction.Log,
           filenameroot     = top30_DE_filename,
           datafxns         = c(csv = top30_genes,
                                tsv = top30_genes))

callModule(heatmap_downloadableTable, "differentialsTable", ss_userAction.Log,
           filenameroot     = differentials_filename,
           downloaddatafxns = list(csv = differentials_table_content,
                                   tsv = differentials_table_content),
           tabledata        = differentials_table_content,
           rownames         = FALSE)


# ----------------------------------------
# --          SHINY SERVER CODE         --
# ----------------------------------------

observeEvent(userData$object, {
    updateSelectizeInput(session,
                         "genesSel",
                         choices  = userData$object$genes,
                         selected = character(0),
                         server = TRUE)

    updateSelectizeInput(session,
                         "genesOnSel",
                         choices  = userData$object$genes,
                         selected = character(0),
                         server = TRUE)

    updateSelectizeInput(session,
                         "genesOffSel",
                         choices  = userData$object$genes,
                         selected = character(0),
                         server = TRUE)
})

observeEvent(c(userData$object, userData$filteredObject), {
    if (is.null(userData$filteredObject)) {
        data_object <- userData$object
    } else {
        data_object <- userData$filteredObject
    }
    cluster_options <- c("All", as.character(sort(unique(data_object$clusters$Cluster))))

    updateSelectizeInput(session,
                         "scatterClusterSel",
                         choices  = cluster_options,
                         selected = "All",
                         server   = FALSE)

    for (item in c("differentialsCluster1Sel", "differentialsCluster2Sel")) {
        updateSelectizeInput(session,
                             item,
                             choices  = cluster_options,
                             selected = "",
                             server   = FALSE)
    }
})

output$summaryTitle <- renderText({
    paste("Summary:", userData$object$meta$title)
})

output$filterOptions <- renderUI({
    body        <- NULL
    filter_list <- userData$object$filter_list

    if (length(names(filter_list)) > 0) {
        header_text <- "Select/Deselect the below options to change the cells included in the visualization data.  All the dataset cells are included (checked) by default."
        body <- lapply(names(filter_list), FUN = function(name) {
                    checkboxGroupInput(inputId  = paste0(get_filter_prefix(userData$objectId), name),
                                       label    = name,
                                       choices  = filter_list[[name]],
                                       selected = filter_list[[name]],
                                       width    = "80%")
        })
    } else {
        header_text <- "No global filters were defined for this object by the object author"
    }
    list(tags$div(tags$br(),
                       tags$h4("Global Filtering"),
                       tags$p(style = "margin:10px;", header_text)),
         tags$div(id = "filtersDiv", body))
})
# Since the filterOptions is on the second (inactive) tab, it's not rendered automatically.
# When switching to this tab, the plot will be rendered again though there is no change. Line below forces it to render.
outputOptions(output, "filterOptions", suspendWhenHidden = FALSE)

output$datasetSummary <- renderUI({
    get_dataset_summary(userData$object)
})

output$differentialsText <- renderUI({
    diff_text <- NULL
    select_area_text <- "select cells on each chart."
    html_diff_string_start <- "<center><em>For cluster differential expression analysis select two different clusters.<br>
                               For sub-cluster analysis choose the same cluster and then"
    if (!is.null(input$differentialsCluster1Sel) && !is.null(input$differentialsCluster2Sel) &&
       (input$differentialsCluster1Sel == input$differentialsCluster2Sel)) {
        diff_text <- paste(html_diff_string_start, paste0("<b>", select_area_text, "</b>"))
    } else {
        diff_text <- paste(html_diff_string_start, select_area_text)
    }
    diff_text <- paste0(diff_text, "</em></center>")
    HTML(diff_text)
})

output$differentialsTableTitle <- renderUI({
    title <- "Differentials Between Selected Clusters"
    if (userSel$runDiffCalc) {
        if (userSel$diffCluster1 == userSel$diffCluster2) {
            title <- paste("Differentials Between Selected Cells in Cluster", userSel$diffCluster1)
        } else {
            title <- paste("Differentials Between Cluster", userSel$diffCluster1, "and Cluster", userSel$diffCluster2)
        }
    }
    tags$strong(title)
})

output$differentialsTableAlternativeText <- renderUI({
    result <- NULL
    diff_table_content <- differentials_table_content()
    if (is.null(diff_table_content) || nrow(diff_table_content) == 0) {
        result <- HTML("<center><em>Press the Calculate Differentials button to perform differential expression analysis on the selected clusters or sub-clusters.</em></center>")
    }
    result
})

output$cxOverviewPlot <- renderCanvasXpress({
    plot <- get_overview_plot(filtered_data(),
                              overview_plot_filename())
    if (is.null(plot)) {
        return(canvasXpress(destroy = TRUE))
    }
    else {
        return(plot)
    }
})

output$cxScatterPlot <- renderCanvasXpress({
    plot <- NULL
    if (userSel$runScatterPlot && !is.null(input$genesSel)) {
        plot <- get_scatter_panel_plot(filtered_data(),
                                       userSel$genes,
                                       userSel$scatterCluster,
                                       scatter_plot_filename(),
                                       userSel$scatterPanelColumns)
    }
    if (is.null(plot)) {
        return(canvasXpress(destroy = TRUE))
    }
    else {
        return(plot)
    }
})

output$cxViolinPlot <- renderCanvasXpress({
    plot <- NULL
    if (userSel$runViolinPlot && !is.null(input$genesSel)) {
        plot <- get_violin_panel_plot(filtered_data(),
                                      userSel$genes,
                                      violin_plot_filename(),
                                      userSel$violinPanelColumns)
    }
    if (is.null(plot)) {
        return(canvasXpress(destroy = TRUE))
    }
    else {
        return(plot)
    }
})

output$cxDotPlot <- renderUI({
    gene_count <- 0

    if (userSel$runDotPlot) {
        plot_result <- get_dot_plot(filtered_data(),
                                    input$genesSel,
                                    get_additional_genes(input$addDotGenes, top10_genes(), top30_genes()),
                                    input$addDotGenes,
                                    dot_plot_filename())
        output$dotplot1 <- renderCanvasXpress({plot_result[[1]]})
        gene_count <- plot_result[[2]]
    } else {
        output$dotplot1 <- renderCanvasXpress({canvasXpress(destroy = TRUE)})
    }

    tagList(
        canvasXpressOutput("dotplot1", height = get_dynamic_plot_height(gene_count))
    )
})

output$cxHeatmapPlot <- renderUI({
    gene_count <- 0

    if (userSel$runHeatmap) {
        plot_result <- get_heatmap_plot(filtered_data(),
                                        input$genesSel,
                                        get_additional_genes(input$addHeatmapGenes, top10_genes(), top30_genes()),
                                        input$addHeatmapGenes,
                                        heatmap_filename())
        output$heatmap1 <- renderCanvasXpress({plot_result[[1]]})
        gene_count <- plot_result[[2]]
    } else {
        output$heatmap1 <- renderCanvasXpress({canvasXpress(destroy = TRUE)})
    }

    tagList(
        canvasXpressOutput("heatmap1", height = get_dynamic_plot_height(gene_count))
    )
})

output$cxDifferentialsScatterPlot1 <- renderCanvasXpress({
    if (!is.null(input$differentialsCluster1Sel) && input$differentialsCluster1Sel != "") {
        get_differential_scatter_plot(isolate(filtered_data()),
                                      input$differentialsCluster1Sel,
                                      paste("Cluster", input$differentialsCluster1Sel),
                                      "cxDifferentialsSelected1",
                                      differentials1_plot_filename())
    } else {
        return(canvasXpress(destroy = TRUE))
    }
})

output$cxDifferentialsScatterPlot2 <- renderCanvasXpress({
    if (!is.null(input$differentialsCluster2Sel) && input$differentialsCluster2Sel != "") {
        get_differential_scatter_plot(isolate(filtered_data()),
                                      input$differentialsCluster2Sel,
                                      paste("Cluster", input$differentialsCluster2Sel),
                                      "cxDifferentialsSelected2",
                                      differentials2_plot_filename())
    } else {
        return(canvasXpress(destroy = TRUE))
    }
})

# observe inputs

observeEvent(input$top10Genes, {
    if (input$top10Genes) {
        updateCheckboxInput(session, "top30Genes", value = FALSE)
    }
})

observeEvent(input$top30Genes, {
    if (input$top30Genes) {
        updateCheckboxInput(session, "top10Genes", value = FALSE)
    }
})

observeEvent(c(input$differentialsCluster1Sel, input$differentialsCluster2Sel), {
    if (input$differentialsCluster1Sel != "" && input$differentialsCluster2Sel != "") {
        enable("diffCalculateBtn")
    } else {
        disable("diffCalculateBtn")
    }
})

# blank out plots when options have changed
observeEvent(input$genesSel, {
    userSel$runScatterPlot <- FALSE
    userSel$runViolinPlot  <- FALSE
    userSel$runDotPlot     <- FALSE
    userSel$runHeatmap     <- FALSE
})

observeEvent(input$scatterClusterSel, {
    if (input$scatterClusterSel != userSel$scatterCluster) {
        userSel$runScatterPlot <- FALSE
    }
})

observeEvent(input$scatterPanelPlotColumns, {
    if (!is.null(userSel$genes) && userSel$genes != "") {
        check_result <- check_panel_plot_columns(userSel$genes, input$scatterPanelPlotColumns)
        if (!check_result[[1]] || (check_result[[1]]) && check_result[[2]] != userSel$scatterPanelColumns) {
            userSel$runScatterPlot <- FALSE
        }
    }
})

observeEvent(input$violinPanelPlotColumns, {
    if (userSel$runViolinPlot) {
        check_result <- check_panel_plot_columns(userSel$genes, input$violinPanelPlotColumns)
        if (!check_result[[1]] || (check_result[[1]]) && check_result[[2]] != userSel$violinPanelColumns) {
            userSel$runViolinPlot <- FALSE
        }
    }
})

observeEvent(input$addDotGenes, {
    if (userSel$runDotPlot) {
        if (input$addDotGenes != userSel$addDotGenes) {
            userSel$runDotPlot <- FALSE
        }
    }
})

observeEvent(input$addHeatmapGenes, {
    if (userSel$runHeatmap) {
        if (input$addHeatmapGenes != userSel$addHeatmapGenes) {
            userSel$runHeatmap <- FALSE
        }
    }
})

observeEvent(c(input$differentialsCluster1Sel, input$differentialsCluster2Sel), {
    if (userSel$runDiffCalc) {
        userSel$runDiffCalc <- FALSE
    }
})


observeEvent(c(input$cxDifferentialsSelected1, input$cxDifferentialsSelected2), {
    if ((input$differentialsCluster1Sel == input$differentialsCluster2Sel) && userSel$runDiffCalc) {
        userSel$runDiffCalc <- FALSE
    }
})

# File to be uploaded has been chosen by user
observeEvent(input$fileChosen, {
    if ((input$fileChosen / (1024*1024)) > 10) {
        toggleModal(session, "loading_modal", toggle = "open")
        userData$modalOpen <- TRUE
    }
})

# File modal closed by user
observeEvent(input$fileModalClosed, {
    userData$modalOpen <- FALSE
})

# File upload has finished
observeEvent(input$fileInputDialog, {
    input_file <- input$fileInputDialog
    if (!is.null(input_file)) {
        load_data_result <- load_data(input_file)
        error_messages   <- load_data_result[["errors"]]
        if (is.null(error_messages)) {
            userData$objectId       <- userData$objectId + 1
            userData$object         <- load_data_result[["object"]]
            userData$filteredObject <- NULL
            disable_plots()
        }

        if (userData$modalOpen) {
            toggleModal(session, "loading_modal", toggle = "close")
            userData$modalOpen <- FALSE
        }

        if (!is.null(error_messages)) {
            output$loading_error_message <- renderText({get_file_dialog_error_message(input_file, error_messages)})
            toggleModal(session, "file_error_modal", toggle = "open")
        }
    }
})

# About app link clicked
observeEvent(input$about_link, {
    session$sendCustomMessage("openTitleInfoBox", runif(1))
})

# Scatter Plot Button Pushed
observeEvent(input$scatterPlotBtn, {
    #clear out old alerts and plots
    try({
        closeAlert(session, "invalidScatterInputAlertID")
        userSel$runScatterPlot <- FALSE
    })

    # check input
    invalidInputMessages <- list()
    if (is.null(input$genesSel) || input$genesSel == "") {
        userSel$genes <- NULL
    } else {
        userSel$genes <- input$genesSel
    }

    if (is.null(userSel$genes)) {
        invalidInputMessages <- append(invalidInputMessages, "No Genes selected in Chart Options.")
    }

    if (is.null(input$scatterClusterSel) || input$scatterClusterSel == "") {
        userSel$scatterCluster <- NULL
        invalidInputMessages <- append(invalidInputMessages, "No Cluster selected.")
    }
    else {
        userSel$scatterCluster <- input$scatterClusterSel
    }

    if (is.null(input$scatterPanelPlotColumns) || input$scatterPanelPlotColumns == "") {
        userSel$scatterPanelColumns <- NULL
        invalidInputMessages <- append(invalidInputMessages, "No Panel Plot Columns selected.")
    }
    else {
        check_result <- check_panel_plot_columns(userSel$genes, input$scatterPanelPlotColumns)
        if (check_result[[1]]) {
            userSel$scatterPanelColumns <- check_result[[2]]
        } else {
            userSel$scatterPanelColumns <- 0
            updateTextInput(session, "scatterPanelPlotColumns", value = "auto")
        }
    }

    if (length(invalidInputMessages) > 0) {
        createAlert(session,
                    "bodyAlert",
                    "invalidScatterInputAlertID",
                    style   = "warning",
                    content = paste(invalidInputMessages, collapse = "<br>"),
                    append  = FALSE)
    } else {
        # trigger plots
        userSel$runScatterPlot <- TRUE
    }
})

# Violin Plot Button Pushed
observeEvent(input$violinPlotBtn, {
    #clear out old alerts and plots
    try({
        closeAlert(session, "invalidViolinInputAlertID")
        userSel$runViolinPlot <- FALSE
    })

    # check input
    invalidInputMessages <- list()
    if (is.null(input$genesSel) || input$genesSel == "") {
        userSel$genes <- NULL
    } else {
        userSel$genes <- input$genesSel
    }

    if (is.null(userSel$genes)) {
        invalidInputMessages <- append(invalidInputMessages, "No Genes selected in Chart Options.")
    }

    if (is.null(input$violinPanelPlotColumns) || input$violinPanelPlotColumns == "") {
        userSel$violinPanelColumns <- NULL
        invalidInputMessages <- append(invalidInputMessages, "No Panel Plot Columns selected.")
    }
    else {
        check_result <- check_panel_plot_columns(userSel$genes, input$violinPanelPlotColumns)
        if (check_result[[1]]) {
            userSel$violinPanelColumns <- check_result[[2]]
        } else {
            userSel$violinPanelColumns <- 0
            updateTextInput(session, "violinPanelPlotColumns", value = "auto")
        }
    }

    if (length(invalidInputMessages) > 0) {
        createAlert(session,
                    "bodyAlert",
                    "invalidViolinInputAlertID",
                    style   = "warning",
                    content = paste(invalidInputMessages, collapse = "<br>"),
                    append  = FALSE)
    } else {
        # trigger plots
        userSel$runViolinPlot <- TRUE
    }
})

# Dot Plot Button Pushed
observeEvent(input$dotPlotBtn, {
    #clear out old alerts and plots
    try({
        closeAlert(session, "invalidDotInputAlertID")
        userSel$runDotPlot <- FALSE
    })

    # check input
    invalidInputMessages <- list()
    if (is.null(input$genesSel) && input$addDotGenes == "off") {
        userSel$addDotGenes <- NULL
        invalidInputMessages <- append(invalidInputMessages, "No Genes selected in Chart Options.")
    } else {
        userSel$addDotGenes <- input$addDotGenes
    }

    if (length(invalidInputMessages) > 0) {
        createAlert(session,
                    "bodyAlert",
                    "invalidDotInputAlertID",
                    style   = "warning",
                    content = paste(invalidInputMessages, collapse = "<br>"),
                    append  = FALSE)
    } else {
        # trigger plots
        userSel$runDotPlot <- TRUE
    }
})

# Heatmap Plot Button Pushed
observeEvent(input$heatmapPlotBtn, {
    #clear out old alerts and plots
    try({
        closeAlert(session, "invalidHeatmapInputAlertID")
        userSel$runHeatmap <- FALSE
    })

    # check input
    invalidInputMessages <- list()
    if (is.null(input$genesSel) && input$addHeatmapGenes == "off") {
        userSel$addHeatmapGenes <- NULL
        invalidInputMessages <- append(invalidInputMessages, "No Genes selected in Chart Options.")
    } else {
        userSel$addHeatmapGenes <- input$addHeatmapGenes
    }

    if (length(invalidInputMessages) > 0) {
        createAlert(session,
                    "bodyAlert",
                    "invalidHeatmapInputAlertID",
                    style   = "warning",
                    content = paste(invalidInputMessages, collapse = "<br>"),
                    append  = FALSE)
    } else {
        # trigger plots
        userSel$runHeatmap <- TRUE
    }
})

# Differential Button Pushed
observeEvent(input$diffCalculateBtn, {
    #clear out old alerts and plots
    try({
        closeAlert(session, "invalidDiffInputAlertID")
        closeAlert(session, "zeroDiffOutputAlertID")
    })

    userSel$runDiffCalc       <- FALSE
    userSel$diffCluster1      <- NULL
    userSel$diffCluster2      <- NULL
    userSel$diffClusterCells1 <- NULL
    userSel$diffClusterCells2 <- NULL

    # check input
    invalidInputMessages <- list()
    if (is.null(input$differentialsCluster1Sel) || input$differentialsCluster1Sel == "" ||
        is.null(input$differentialsCluster2Sel) || input$differentialsCluster2Sel == "") {
        invalidInputMessages <- append(invalidInputMessages, "Both Cluster 1 and Cluster 2 must be selected.")
    }
    else {
        userSel$diffCluster1 <- input$differentialsCluster1Sel
        userSel$diffCluster2 <- input$differentialsCluster2Sel

        # Combination of All and cluster X not possible
        if ((userSel$diffCluster1 == "All" && userSel$diffCluster2 != "All") ||
            (userSel$diffCluster1 != "All" && userSel$diffCluster2 == "All")) {
            cluster <- setdiff(c(userSel$diffCluster1, userSel$diffCluster2), "All")
            message <- paste("Differential analysis between All and Cluster", cluster, "is not available.", "
                             If you wish to perform sub-cluster analysis on the entire dataset choose All for both plots and select cells on each plot to compute the differentials.")
            invalidInputMessages <- append(invalidInputMessages, message)
        }
        else if (userSel$diffCluster1 == userSel$diffCluster2) {
            # Checks if the same cluster selected
            if (is.null(input$cxDifferentialsSelected1) || length(input$cxDifferentialsSelected1) < 1 ||
                is.null(input$cxDifferentialsSelected2) || length(input$cxDifferentialsSelected2) < 1) {
                message <- paste("Cells on both plots must be selected to perform differential sub-cluster analysis.")
                invalidInputMessages <- append(invalidInputMessages, message)
            }
            else if (all(input$cxDifferentialsSelected1 %in% input$cxDifferentialsSelected2) &&
                     all(input$cxDifferentialsSelected2 %in% input$cxDifferentialsSelected1)) {
                message <- paste("The same cells are selected on both plots.  To perform differential sub-cluster analysis there must be a difference in the cells selected.")
                invalidInputMessages <- append(invalidInputMessages, message)
            }
            else {
                t1 <- input$cxDifferentialsSelected1
                t2 <- input$cxDifferentialsSelected2

                userSel$diffClusterCells1 <- input$cxDifferentialsSelected1
                userSel$diffClusterCells2 <- input$cxDifferentialsSelected2
            }
        }
    }

    if (length(invalidInputMessages) > 0) {
        createAlert(session,
                    "bodyAlert",
                    "invalidDiffInputAlertID",
                    style   = "warning",
                    content = paste(invalidInputMessages, collapse = "<br>"),
                    append  = FALSE)
    } else {
        userSel$runDiffCalc <- TRUE
    }
})
