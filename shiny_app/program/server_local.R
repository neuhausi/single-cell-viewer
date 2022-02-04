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
#     ss_userid - userID, authenticated
#     ss_userAction.Log - Reactive Logger S4 object
# ----------------------------------------

# -- IMPORTS --


# -- VARIABLES --

userData <- reactiveValues()
    userData$object              <- NULL

userSel <- reactiveValues()
    userSel$clusterColumn        <- NULL
    userSel$genes                <- NULL
    userSel$geneSignatures       <- NULL
    userSel$scatterCluster       <- "All"
    userSel$scatterPanelColumns  <- NULL
    userSel$runScatterPlot       <- FALSE
    userSel$violinPanelColumns   <- NULL
    userSel$runViolinPlot        <- FALSE
    userSel$runHeatmap           <- FALSE
    userSel$diffCluster1         <- NULL
    userSel$diffCluster2         <- NULL
    userSel$diffClusterCells1    <- NULL
    userSel$diffClusterCells2    <- NULL
    userSel$runDiffCalc          <- FALSE
    userSel$runCoExpPlot         <- FALSE
    userSel$coExpCluster         <- "All"
    userSel$xGene                <- NULL
    userSel$yGene                <- NULL
    userSel$xGate                <- 0.5
    userSel$yGate                <- 0.5
    userSel$noise                <- 0.01

# cx expression calculated/input fields initialization
cx_exp <- reactiveValues()
    cx_exp$coExTableData      <- NULL
    cx_exp$coExTableFileName  <- NULL
    cx_exp$cells              <- NULL
    cx_exp$coGeneXData        <- NULL
    cx_exp$coGeneYData        <- NULL
    cx_exp$pearson            <- NULL
    cx_exp$plot               <- NULL

# differentials Table heatmap colors
heatmap_brks <- seq(0, 100, by = 1)
heatmap_clrs <- colorRampPalette(c("deepskyblue3", "white", "firebrick3"), space = "rgb")(length(heatmap_brks) + 1)

# -- FUNCTIONS --

disable_plots <- function() {
    if (g_debug) message('function: ', 'disable_plots')

    userSel$runScatterPlot <- FALSE
    userSel$runViolinPlot  <- FALSE
    userSel$runHeatmap     <- FALSE
    userSel$runDiffCalc    <- FALSE
    userSel$runCoExpPlot   <- FALSE
}


differentials_table_content <- reactive({
    if (g_debug) message('reactive: ', 'differentials_table_content')

    c1 <- isolate(userSel$diffCluster1)
    s1 <- isolate(userSel$diffClusterCells1)
    c2 <- isolate(userSel$diffCluster2)
    s2 <- isolate(userSel$diffClusterCells2)

    result <- NULL
    if (userSel$runDiffCalc) {
        result <- calculate_differentials(isolate(userData$object), c1, s1, c2, s2)
        if (is.null(result)) {
            createAlert(session,
                        "bodyAlert",
                        "chartAlert",
                        style   = "warning",
                        content = paste("Differential Analysis resulted in no output for the current selections and logFC > ",
                                        g_differential_logfc_threshold),
                        append  = FALSE)
        }
    }
    if (!is.null(result)) {
        result <- result %>%
            arrange(desc(abs(log_FC)), p_val)
    }
    result
})

coEx_table_content <- reactive({
    if (g_debug) message('reactive: ', 'coEx_table_content')
    cx_exp$coExTableData
})

coEx_table_file_name <- reactive({
    if (g_debug) message('reactive: ', 'coEx_table_file_name')
    isolate(cx_exp$coExTableFileName)
})

base_filename <- reactive({
    if (g_debug) message('reactive: ', 'base_filename')

    current_time <- format(Sys.time(), "%Y.%m.%d_%H.%M")
    paste0(current_time, "_SCV_", userData$object$meta$object_name)
})


# -- Downloadable Tables --

downloadableTable("differentialsTable",
                  logger           = ss_userAction.Log,
                  filenameroot     = base_filename,
                  downloaddatafxns = list(csv = differentials_table_content,
                                          tsv = differentials_table_content),
                  tabledata        = differentials_table_content,
                  formatStyle      = list(columns = c(5, 6),
                                          backgroundColor = styleInterval(heatmap_brks, heatmap_clrs)))

downloadableTable("coExpTable",
                  logger           = ss_userAction.Log,
                  filenameroot     = coEx_table_file_name,
                  downloaddatafxns = list(csv = coEx_table_content, tsv = coEx_table_content),
                  tabledata        = coEx_table_content,
                  rownames         = FALSE)

# ----------------------------------------
# --          SHINY SERVER CODE         --
# ----------------------------------------

# Disable Plots at Startup
disable_plots()

# Load Dataset
observeEvent(TRUE, {
    if (any(is.null(g_clusters_path),
            length(g_clusters_path) == 0)) {
        stopApp("Invalid 'g_clusters_path' has been specified.")
    } else {
        allowed_clusters <- load_clusters(g_clusters_path)
    }

    if (any(is.null(g_default_cluster),
            g_default_cluster == "",
            !is.character(g_default_cluster),
            !(g_default_cluster %in% names(allowed_clusters)))) {
        # Business requirement: must have a valid default cluster in the app configuration. The default cluster must
        # exist in the database, but doesn't have to exist in each individual dataset.
        stopApp("Invalid 'g_default_cluster' has been specified.")
    } else {
        object <- load_data(get_url_parameters(session))

        if (is.null(object)) {
            parms <- get_url_parameters(session)
            msg   <- case_when(
                !is.null(parms$isetID) && is.null(parms$projectID) && is.null(parms$sampleID)  ~ paste0('Unable to load data for isetID: ',    parms$isetID),
                is.null(parms$isetID) && !is.null(parms$projectID) && !is.null(parms$sampleID) ~ paste0('Unable to load data for projectID: ', parms$projectID, ' and sampleID: ', parms$sampleID),
                TRUE                                                                           ~ "URL parameters were specified incorrectly for use in this application")

            # alert that the data was not available
            createAlert(session, "bodyAlert", "fileAlert", style = "danger", content = msg)
        } else {
            # select the columns if present
            cells_columns   <- names(object$cells)
            cluster_columns <- names(allowed_clusters)
            clusters        <- intersect(cells_columns, cluster_columns)

            if (any(is.null(clusters),
                    length(clusters) == 0)) {
                msg <- "The requested dataset does not contain the cluster information needed for display in this viewer.  Contact the dataset owner for further information."
                createAlert(session, "bodyAlert", "fileAlert", style = "danger", content = msg)

            } else {
                if (g_default_cluster %in% clusters) {
                    default_cluster <- g_default_cluster
                } else {
                    default_cluster <- clusters[1]
                }

                # format the cluster choices, swap the names and the values to work for updateSelectizeInput
                named_clusters         <- allowed_clusters[clusters]
                cluster_choices        <- names(named_clusters)
                names(cluster_choices) <- named_clusters

                updateSelectizeInput(session,
                                     "clustersSel",
                                     choices  = cluster_choices,
                                     selected = default_cluster,
                                     server = TRUE)
                userSel$clusterColumn <- default_cluster
                userData$object       <- update_cluster_options(session, object, default_cluster)
                updateSelectizeInput(session,
                                     "genesSel",
                                     choices  = userData$object$genes$Symbol,
                                     selected = character(0),
                                     server = TRUE)
            }
        }
    }
}, once = TRUE)

observeEvent(userData$object, {
    if (g_debug) message('observeEvent: ', 'userData', '->', 'setup application after load of data')

    output$cxOverviewPlot <- renderCanvasXpress({
        plot <- get_overview_plot(userData$object,
                                  base_filename())
        if (is.null(plot)) {
            return(canvasXpress(destroy = TRUE))
        } else {
            return(plot)
        }
    })
})

observeEvent(input$clustersSel, {
    if (g_debug) message("observeEvent: ", "input$clustersSel", "->", "update cluster signature")

    if (all(!is.null(input$clustersSel),
            input$clustersSel != "",
            input$clustersSel != isolate(userSel$clusterColumn),
            input$clustersSel %in% names(userData$object$cells))) {
        open_cluster_change_modal(session)
    } else {
        updateSelectizeInput(session, "clustersSel", selected = isolate(userSel$clusterColumn))
    }

}, ignoreInit = TRUE)

observeEvent(userData, {
    if (g_debug) message('observeEvent: ', 'userData', '->', 'update signature selections')

    signatures         <- names(g_gene_signatures)
    sig_choices        <- unlist(lapply(signatures, FUN = function(x) paste(g_gene_signatures[[x]], collapse = ", ")))
    names(sig_choices) <- signatures
    updateSelectizeInput(session,
                         "geneSignaturesSel",
                         choices  = sig_choices,
                         selected = character(0),
                         server   = FALSE,
                         options  = list(
                             placeholder = "Type/Click then Select",
                             searchField = "value",
                             plugins = list('remove_button'),
                             render  = I( "{ option: function(item, escape) {
                                          return '<div><b>' + item.label + '</b><br>' + '<i><small>' + item.value + '</small></i></div>'; }}")
                         ))
})

output$summaryTitle <- renderText({
    paste(userData$object$meta$title)
})

output$datasetSummary <- renderUI({
    get_dataset_summary(userData$object)
})

output$differentialsTableTitle <- renderUI({
    title <- "Differential Expression Analysis Results"
    if (userSel$runDiffCalc) {
        if (userSel$diffCluster1 == userSel$diffCluster2) {
            title <- paste("Differentials Between Selected Cells in Cluster", userSel$diffCluster1)
        } else {
            title <- paste("Differentials Between Cluster", userSel$diffCluster1, "and Cluster", userSel$diffCluster2)
        }
    }
    tags$strong(title)
})

output$coExpTableTitle <- renderUI({
    title <- "Population Combinations"
    tags$strong(title)
})

output$selected_gene_signatures_text <- renderUI({
    result <- ""
    if (!is.null(input$geneSignaturesSel) && input$geneSignaturesSel != "") {
        selected_signatures <- input$geneSignaturesSel
        non_expressed_genes <- FALSE
        for (signature_genes in selected_signatures) {
            signature_genes <- unlist(strsplit(signature_genes,", "))
            signature <- unlist(lapply(names(g_gene_signatures), FUN = function(x) if (all(g_gene_signatures[[x]] %in% signature_genes)) {x}))
            result <- paste(result,"<br><b>", signature, "</b><br>")
            sig_genes_length <- length(signature_genes)
            for (i in 1:sig_genes_length) {
                gene <- signature_genes[i]
                if (gene %in% userData$object$genes$Symbol) {
                    result <- paste(result, gene)
                } else {
                    non_expressed_genes <- TRUE
                    result <- paste(result, "<i>", paste0("*", gene), "</i>")
                }
                result <- paste0(result, ifelse(i == sig_genes_length, "", ","))
            }
        }
        if (non_expressed_genes) {
            result <- paste(result, "<br><br><i>*Not Expressed in Dataset</i>")
        }
    }
    HTML(result)
})

output$cxScatterPlot <- renderCanvasXpress({
    plot <- NULL
    if (userSel$runScatterPlot && (!is.null(input$geneSignaturesSel) || !is.null(input$genesSel))) {
        plot <- get_scatter_panel_plot(userData$object,
                                       userSel$genes,
                                       userSel$geneSignatures,
                                       userSel$scatterCluster,
                                       base_filename(),
                                       userSel$scatterPanelColumns)
        if (is.null(plot)) {
            createAlert(session,
                        "bodyAlert",
                        "chartAlert",
                        style   = "warning",
                        content = "There are no expressed genes in the Signature(s) selected.",
                        append  = FALSE)
        }
    }
    if (is.null(plot)) {
        return(canvasXpress(destroy = TRUE))
    } else {
        return(plot)
    }
})

output$cxViolinPlot <- renderCanvasXpress({
    plot <- NULL
    if (userSel$runViolinPlot && (!is.null(input$geneSignaturesSel) || !is.null(input$genesSel))) {
        plot <- get_violin_panel_plot(userData$object,
                                      userSel$genes,
                                      userSel$geneSignatures,
                                      base_filename(),
                                      userSel$violinPanelColumns)
        if (is.null(plot)) {
            createAlert(session,
                        "bodyAlert",
                        "chartAlert",
                        style   = "warning",
                        content = "There are no expressed genes in the Signature(s) selected.",
                        append  = FALSE)
        }
    }
    if (is.null(plot)) {
        return(canvasXpress(destroy = TRUE))
    } else {
        return(plot)
    }
})

output$cxCoExpPlot <- renderCanvasXpress({
    if (userSel$runCoExpPlot &&
        !any(is.null(input$cxRefreshBtn),
             is.null(input$coExpPlotBtn))) {
        hideFeedback("xGate")
        hideFeedback("yGate")
        hideFeedback("noise")
        isolate({
            co_expressions_calc_output <- calculate_co_expressions_genes_data()
            if (co_expressions_calc_output$error != "") {
                msg <- glue("Could not calculate Pearson correlation between '{co_expressions_calc_output$result$xGene}' and '{co_expressions_calc_output$result$yGene}' due to: {co_expressions_calc_output$error}")
                logerror(msg, logger = ss_userAction.Log)
            }

            # the controls are not ready for display yet
            if (!any(is.null(input$xGate),
                     is.null(input$yGate),
                     is.null(input$noise))) {
                xGate <- input$xGate
                yGate <- input$yGate
                noise <- input$noise

                # check xGate
                if (any(is.na(xGate), is.null(xGate), xGate == "")) {
                    showFeedbackDanger("xGate", "X Gate value cannot be empty")
                    return(cx_exp$plot)
                }

                xGate <- as.numeric(xGate)

                if (is.na(xGate)) {
                    showFeedbackDanger("xGate", "X Gate must be numeric")
                    return(cx_exp$plot)
                }

                if (xGate <= 0) {
                    showFeedbackDanger("xGate", "X Gate must be greater than 0")
                    return(cx_exp$plot)
                }

                if (xGate > co_expressions_calc_output$result$xGateMax) {
                    showFeedbackWarning("xGate", "X Gate is out of range for the chart values")
                }

                # check yGate
                if (any(is.na(yGate), is.null(yGate), yGate == "")) {
                    showFeedbackDanger("yGate", "Y Gate value cannot be empty")
                    return(cx_exp$plot)
                }

                yGate <- as.numeric(yGate)

                if (is.na(yGate)) {
                    showFeedbackDanger("yGate", "Y Gate must be numeric")
                    return(cx_exp$plot)
                }

                if (yGate <= 0) {
                    showFeedbackDanger("yGate", "Y Gate must be greater than 0")
                    return(cx_exp$plot)
                }

                if (yGate > co_expressions_calc_output$result$yGateMax) {
                    showFeedbackWarning("yGate", "Y Gate is out of range for the chart values")
                }

                # check noise
                if (any(is.na(noise), is.null(noise), noise == "")) {
                    showFeedbackDanger("noise", "Noise value cannot be empty")
                    return(cx_exp$plot)
                }

                noise <- as.numeric(noise)
                if (is.na(noise)) {
                    showFeedbackDanger("noise", "Noise must be numeric")
                    return(cx_exp$plot)
                }

                if (noise < 0) {
                    showFeedbackDanger("noise", "Noise value must be greater than or equal 0")
                    return(cx_exp$plot)
                }
            }

        })

        cx_exp$cells             <- co_expressions_calc_output$result$cells
        cx_exp$coGeneXData       <- co_expressions_calc_output$result$coGeneXData
        cx_exp$coGeneYData       <- co_expressions_calc_output$result$coGeneYData
        cx_exp$pearson           <- co_expressions_calc_output$result$pearson
        plot_data <- get_co_exp_plot(userData$object,
                                     co_expressions_calc_output$result$cells,
                                     co_expressions_calc_output$result$coGeneXData,
                                     co_expressions_calc_output$result$coGeneYData,
                                     base_filename(),
                                     xGate,
                                     yGate,
                                     noise,
                                     co_expressions_calc_output$result$pearson)

        cx_exp$coExTableData     <- plot_data$percentages
        cx_exp$plot              <- plot_data$plot
        cx_exp$coExTableFileName <- glue("{base_filename()}_{plot_data$geneX}_{plot_data$geneY}_Population_Combinations")
        cx_exp$plot

    } else {
        cx_exp$coExTableData <-  data.frame()
        canvasXpress(destroy = TRUE)
    }
})

output$cxNoiseControls <- renderUI({
    if (userSel$runCoExpPlot) {
        tagList(
            tags$div(
                tags$br(),
                tags$br(),
                tags$br(),
                tags$br(),
                tags$br(),
                tags$br(),
                tags$br(),
                tags$br(),
                numericInput(inputId = "xGate",
                             label   = ui_tooltip("xGate_lbl",
                                                  "X Gate",
                                                  "Expression gate value for the x-axis gene (red line on the plot)"),
                             value   = userSel$xGate,
                             width   = "100%"),
                numericInput(inputId = "yGate",
                             label   = ui_tooltip("yGate_lbl",
                                                  "Y Gate",
                                                  "Expression gate value for the y-axis gene (blue line on the plot)"),
                             value   = userSel$yGate,
                             width   = "100%"),
                numericInput(inputId    = "noise",
                             label      = ui_tooltip("noise_lbl",
                                                     "Noise ",
                                                     "Amount of random noise added to expression values to enhance plot visibility"),
                             value      = userSel$noise,
                             width      = "100%",
                             min        = 0),
                bsButton(inputId  = "cxRefreshBtn",
                         label    = "Refresh",
                         style    = "primary",
                         width    = "100%",
                         icon     = icon("arrow-left"))
            )
        )
    } else {
        tagList()
    }
})

output$cxHeatmapPlot <- renderUI({
    plot       <- NULL
    gene_count <- 0

    if (userSel$runHeatmap) {
        plot_result <- get_heatmap_plot(userData$object,
                                        input$genesSel,
                                        input$geneSignaturesSel,
                                        NULL,
                                        base_filename())
        plot       <- plot_result[[1]]
        gene_count <- plot_result[[2]]
        if (is.null(plot)) {
            createAlert(session,
                        "bodyAlert",
                        "chartAlert",
                        style   = "warning",
                        content = "There are no expressed genes in the Signature(s) selected.",
                        append  = FALSE)
        } else {
            output$heatmap1 <- renderCanvasXpress({plot})
        }

    }

    if (is.null(plot)) {
        output$heatmap1 <- renderCanvasXpress({canvasXpress(destroy = TRUE)})
    }

    tagList(
        canvasXpressOutput("heatmap1", height = get_dynamic_plot_height(gene_count))
    )
})

output$cxDifferentialsScatterPlot1 <- renderCanvasXpress({
    if (!is.null(input$differentialsCluster1Sel) && input$differentialsCluster1Sel != "") {
        get_differential_scatter_plot(isolate(userData$object),
                                      input$differentialsCluster1Sel,
                                      paste("Cluster", input$differentialsCluster1Sel),
                                      "cxDifferentialsSelected1",
                                      paste0(base_filename(), '_Differentials1'))
    } else {
        return(canvasXpress(destroy = TRUE))
    }
})

output$cxDifferentialsScatterPlot2 <- renderCanvasXpress({
    if (!is.null(input$differentialsCluster2Sel) && input$differentialsCluster2Sel != "") {
        get_differential_scatter_plot(isolate(userData$object),
                                      input$differentialsCluster2Sel,
                                      paste("Cluster", input$differentialsCluster2Sel),
                                      "cxDifferentialsSelected2",
                                      paste0(base_filename(), "_Differentials2"))
    } else {
        return(canvasXpress(destroy = TRUE))
    }
})

# observe inputs

observeEvent(c(input$differentialsCluster1Sel, input$differentialsCluster2Sel), {
    if (g_debug) message('observeEvent: ', 'c(input$differentialsCluster1Sel, input$differentialsCluster2Sel)', '->', 'enable/disable calculate button')

    if (input$differentialsCluster1Sel != "" && input$differentialsCluster2Sel != "") {
        enable("diffCalculateBtn")
    } else {
        disable("diffCalculateBtn")
    }
})

# blank out plots when options have changed
observeEvent(c(input$genesSel, input$geneSignaturesSel), {
    if (g_debug) message('observeEvent: ', 'c(input$genesSel, input$geneSignaturesSel)', '->', 'blank out plots if genes/sigs change')

    userSel$runScatterPlot <- FALSE
    userSel$runViolinPlot  <- FALSE
    userSel$runHeatmap     <- FALSE
})

observeEvent(input$scatterClusterSel, {
    if (g_debug) message('observeEvent: ', 'input$scatterClusterSel', '->', 'run scatterplot')

    if (is.null(userSel$scatterCluster) || input$scatterClusterSel != userSel$scatterCluster) {
        userSel$runScatterPlot <- FALSE
    }
})

observeEvent(input$scatterPanelPlotColumns, {
    if (g_debug) message('observeEvent: ', 'input$scatterPanelPlotColumns', '->', 'change scatter plot column preference')

    if (!is.null(userSel$genes) && userSel$genes != "") {
        check_result <- check_panel_plot_columns(userSel$genes, input$scatterPanelPlotColumns)
        if (!check_result[[1]] || (check_result[[1]]) && check_result[[2]] != userSel$scatterPanelColumns) {
            userSel$runScatterPlot <- FALSE
        }
    }
})

observeEvent(input$violinPanelPlotColumns, {
    if (g_debug) message('observeEvent: ', 'input$violinPanelPlotColumns', '->', 'change violin plot column preference')

    if (userSel$runViolinPlot) {
        check_result <- check_panel_plot_columns(userSel$genes, input$violinPanelPlotColumns)
        if (!check_result[[1]] || (check_result[[1]]) && check_result[[2]] != userSel$violinPanelColumns) {
            userSel$runViolinPlot <- FALSE
        }
    }
})

# Scatter Plot Button Pushed
observeEvent(input$scatterPlotBtn, {
    if (g_debug) message('observeEvent: ', 'input$scatterPlotBtn', '->', 'scatterplot button')

    #clear out old alerts and plots
    try({
        closeAlert(session, "chartAlert")
        userSel$runScatterPlot <- FALSE
    })

    # check input
    invalidInputMessages <- list()
    if (is.null(input$genesSel) || input$genesSel == "") {
        userSel$genes <- NULL
    } else {
        userSel$genes <- input$genesSel
    }

    if (is.null(input$geneSignaturesSel) || input$geneSignaturesSel == "") {
        userSel$geneSignatures <- NULL
    } else {
        userSel$geneSignatures <- input$geneSignaturesSel
    }

    if (is.null(userSel$genes) && is.null(userSel$geneSignatures)) {
        invalidInputMessages <- append(invalidInputMessages, "There are no Genes or Signatures selected.")
    }

    if (is.null(input$scatterClusterSel) || input$scatterClusterSel == "") {
        userSel$scatterCluster <- NULL
        invalidInputMessages <- append(invalidInputMessages, "No Cluster selected.")
    } else {
        userSel$scatterCluster <- input$scatterClusterSel
    }

    if (is.null(input$scatterPanelPlotColumns) || input$scatterPanelPlotColumns == "") {
        userSel$scatterPanelColumns <- NULL
        invalidInputMessages <- append(invalidInputMessages, "No Panel Plot Columns selected.")
    } else {
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
                    "chartAlert",
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
    if (g_debug) message('observeEvent: ', 'input$violinPlotBtn', '->', 'violinplot button')

    #clear out old alerts and plots
    try({
        closeAlert(session, "chartAlert")
        userSel$runViolinPlot <- FALSE
    })

    # check input
    invalidInputMessages <- list()
    if (is.null(input$genesSel) || input$genesSel == "") {
        userSel$genes <- NULL
    } else {
        userSel$genes <- input$genesSel
    }

    if (is.null(input$geneSignaturesSel) || input$geneSignaturesSel == "") {
        userSel$geneSignatures <- NULL
    } else {
        userSel$geneSignatures <- input$geneSignaturesSel
    }

    if (is.null(userSel$genes) && is.null(userSel$geneSignatures)) {
        invalidInputMessages <- append(invalidInputMessages, "There are no Genes or Signatures selected.")
    }

    if (is.null(input$violinPanelPlotColumns) || input$violinPanelPlotColumns == "") {
        userSel$violinPanelColumns <- NULL
        invalidInputMessages <- append(invalidInputMessages, "No Panel Plot Columns selected.")
    } else {
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
                    "chartAlert",
                    style   = "warning",
                    content = paste(invalidInputMessages, collapse = "<br>"),
                    append  = FALSE)
    } else {
        # trigger plots
        userSel$runViolinPlot <- TRUE
    }
})

# Heatmap Plot Button Pushed
observeEvent(input$heatmapPlotBtn, {
    if (g_debug) message('observeEvent: ', 'input$heatmapPlotBtn', '->', 'heatmap button')

    #clear out old alerts and plots
    try({
        closeAlert(session, "chartAlert")
        userSel$runHeatmap <- FALSE
    })

    # check input
    invalidInputMessages <- list()
    if (is.null(input$genesSel) && is.null(input$geneSignaturesSel)) {
        invalidInputMessages <- append(invalidInputMessages, "There are no Genes or Signatures selected.")
    }

    if (length(invalidInputMessages) > 0) {
        createAlert(session,
                    "bodyAlert",
                    "chartAlert",
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
    if (g_debug) message('observeEvent: ', 'input$diffCalculateBtn', '->', 'DE calculate button')

    #clear out old alerts and plots
    try({
        closeAlert(session, "chartAlert")
    })

    userSel$diffCluster1      <- NULL
    userSel$diffCluster2      <- NULL
    userSel$diffClusterCells1 <- NULL
    userSel$diffClusterCells2 <- NULL
    userSel$runDiffCalc       <- FALSE

    # check input
    invalidInputMessages <- list()
    if (is.null(input$differentialsCluster1Sel) || input$differentialsCluster1Sel == "" ||
        is.null(input$differentialsCluster2Sel) || input$differentialsCluster2Sel == "") {
        invalidInputMessages <- append(invalidInputMessages, "Both Cluster 1 and Cluster 2 must be selected.")
    } else {
        userSel$diffCluster1 <- input$differentialsCluster1Sel
        userSel$diffCluster2 <- input$differentialsCluster2Sel

        # Combination of All and cluster X not possible
        if ((userSel$diffCluster1 == "All" && userSel$diffCluster2 != "All") ||
            (userSel$diffCluster1 != "All" && userSel$diffCluster2 == "All")) {
            cluster <- setdiff(c(userSel$diffCluster1, userSel$diffCluster2), "All")
            message <- paste("Differential analysis between All and Cluster", cluster, "is not available.", "
                             If you wish to perform sub-cluster analysis on the entire dataset choose All for both plots and select cells on each plot to compute the differentials.")
            invalidInputMessages <- append(invalidInputMessages, message)
        } else if (userSel$diffCluster1 == userSel$diffCluster2) {
            # Checks if the same cluster selected
            if (is.null(input$cxDifferentialsSelected1) || length(input$cxDifferentialsSelected1) < 1 ||
                is.null(input$cxDifferentialsSelected2) || length(input$cxDifferentialsSelected2) < 1) {
                message <- paste("Cells on both plots must be selected to perform differential sub-cluster analysis.")
                invalidInputMessages <- append(invalidInputMessages, message)
            } else if (all(input$cxDifferentialsSelected1 %in% input$cxDifferentialsSelected2) &&
                     all(input$cxDifferentialsSelected2 %in% input$cxDifferentialsSelected1)) {
                message <- paste("The same cells are selected on both plots.  To perform differential sub-cluster analysis there must be a difference in the cells selected.")
                invalidInputMessages <- append(invalidInputMessages, message)
            } else {
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
                    "chartAlert",
                    style   = "warning",
                    content = paste(invalidInputMessages, collapse = "<br>"),
                    append  = FALSE)
    } else {
        de_total <- length(input$cxDifferentialsSelected1) + length(input$cxDifferentialsSelected2)
        if (userSel$diffCluster1 != userSel$diffCluster2) {
            de_total <- isolate(userData$object$cells) %>%
                dplyr_filter(Cluster %in% c(userSel$diffCluster1, userSel$diffCluster2)) %>%
                NROW()
        }

        if ((de_total > g_de_message_low) && (de_total < g_de_message_high)) {
            output$dacText <- renderUI({
                HTML("You have requested an analysis of over",
                     g_de_message_low, "cells.<br/><br/>",
                     "This may take 30s or more to calculate...")
            })
            open_dac_modal()
        } else if (de_total > g_de_message_high) {
            output$dacText <- renderUI({
                HTML("You have requested an analysis of over",
                     g_de_message_high, "cells.<br/><br/>",
                     "This may take several minutes to calculate...")
            })
            open_dac_modal()
        } else {
            userSel$runDiffCalc <- TRUE
        }
    }
})

observeEvent(input$proceed, {
    disable("proceed")
    disable("cancel")
    toggleModal(session, "dac_modal", toggle = "close")
    userSel$runDiffCalc <- TRUE
})

observeEvent(input$cancel, {
    toggleModal(session, "dac_modal", toggle = "close")
})

# Either coExpPlotBtn or cxRefreshBtn buttons got pushed
observeEvent(any(input$cxRefreshBtn > 0,
                 input$coExpPlotBtn > 0), {
    if (g_debug) message('observeEvent: ', 'input$coExpPlotBtn', '->', 'coExpPlotBtn button')
    # clear out old alerts and plots
    try({
        closeAlert(session, "chartAlert")
    })

    invalidInputMessages <- validate_co_expressions_genes_input()

    if (length(invalidInputMessages) > 0) {
        userSel$runCoExpPlot <- FALSE
        createAlert(session,
                    "bodyAlert",
                    "chartAlert",
                    style   = "warning",
                    content = paste(invalidInputMessages, collapse = "<br>"),
                    append  = FALSE)
    } else {
        # trigger plots
        userSel$runCoExpPlot <- TRUE
    }
}, ignoreInit = TRUE)

# helper methods
open_dac_modal <- function() {
    enable("proceed")
    enable("cancel")
    toggleModal(session, "dac_modal", toggle = "open")
    runjs("$('#dac_modal').on('shown.bs.modal', function(e) {$('#proceed').focus();});")
}

validate_co_expressions_genes_input <- function() {
    invalidInputMessages <- list()
    if (is.null(input$gene1Sel) || input$gene1Sel == "") {
        invalidInputMessages <- append(invalidInputMessages, "Co-Expression Analysis: Gene 1 must be selected to create the chart.")
    } else {
        userSel$xGene <- input$gene1Sel
    }

    if (is.null(input$gene2Sel) || input$gene2Sel == "") {
        invalidInputMessages <- append(invalidInputMessages, "Co-Expression Analysis: Gene 2 must be selected to create the chart.")
    } else {
        userSel$yGene <- input$gene2Sel
    }

    if (is.null(input$coExpClusterSel) || input$coExpClusterSel == "") {
        invalidInputMessages <- append(invalidInputMessages, "Co-Expression Analysis: No Cluster selected. Cluster must be selected to create the chart.")
    } else {
        userSel$coExpCluster <- input$coExpClusterSel
    }

    invalidInputMessages
}

calculate_co_expressions_genes_data <- function() {
    error <- ""
    result <- list()
    cells    <- get_cluster_filtered_cells(userData$object, userSel$coExpCluster)
    gene_data <- get_gene_filtered_expression(userData$object, c(userSel$xGene, userSel$yGene)) %>%
        dplyr_filter(Cell %in% cells$Cell)

    geneX_data <- gene_data %>% select(Cell, all_of(userSel$xGene))
    geneY_data <- gene_data %>% select(Cell, all_of(userSel$yGene))

    output <- calculate_gene_x_gene_pearson_correlation(geneX_data[[userSel$xGene]],
                                                        geneY_data[[userSel$yGene]])
    if (output$error != "") {
        output$result <- NA
        error         <- output$error
    }

    result$coGeneXData  <- geneX_data
    result$coGeneYData  <- geneY_data
    result$pearson      <- output$result
    result$xGateMax     <- max(geneX_data[[userSel$xGene]])
    result$xGateMin     <- min(geneX_data[[userSel$xGene]])
    result$yGateMax     <- max(geneY_data[[userSel$yGene]])
    result$yGateMin     <- min(geneY_data[[userSel$yGene]])
    result$cells        <- cells
    result$xGene        <- userSel$xGene
    result$yGene        <- userSel$yGene
    list(result = result, error = error)
}


observeEvent(input$cluster_proceed, {
    disable("cluster_proceed")
    disable("cluster_cancel")
    toggleModal(session, "cluster_change_modal", toggle = "close")
    updateTabsetPanel(session, "tabBodyElement2", selected = "Overview")
    userSel$clusterColumn <- input$clustersSel
    userData$object <- update_cluster_options(session, isolate(userData$object), input$clustersSel)
    try({
        closeAlert(session, "chartAlert")
    })
    disable_plots()
})

observeEvent(input$cluster_cancel, {
    toggleModal(session, "cluster_change_modal", toggle = "close")
    updateSelectizeInput(session, "clustersSel", selected = isolate(userSel$clusterColumn))
})

observeEvent(input$cluster_hide, {
    if (input$cluster_hide) {
        updateSelectizeInput(session, "clustersSel", selected = isolate(userSel$clusterColumn))
    }
})
