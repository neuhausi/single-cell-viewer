# ----------------------
# Plot Related Functions
# ----------------------
suppressPackageStartupMessages(library(htmlwidgets))

get_overview_plot <- function(dataobj, base_file_name, report_modus = FALSE) {
    if (g_debug) message('function: ', 'get_overview_plot')

    plot <- NULL

    if (!is.null(dataobj)) {
        if (!is.null(dataobj$cells$red_axis1) && !is.null(dataobj$cells$red_axis2)) {
            ydata <- data.frame(PC1 = dataobj$cells$red_axis1,
                                PC2 = dataobj$cells$red_axis2)
            rownames(ydata) <- dataobj$cells$Cell

            zdata <- dataobj$cells %>%
                mutate(Cluster = as.character(Cluster)) %>%
                select(-Cell, -CellType.select, #must be present
                       -matches("red_axis"), -matches("barcode"), -matches("CLID"), -matches("ident")) #others
            rownames(zdata) <- dataobj$cells$Cell

            events <- JS("{ 'mousemove' : function(o, e, t) {
                          if (o != null) {
                            if (o.objectType == null) {
                             t.showInfoSpan(e, '<b>' + o.z.Cluster + '</b>');
                            } else {
                             t.showInfoSpan(e, o.display);
                          }; }}}")

            plot <- canvasXpress(
                data                     = ydata,
                varAnnot                 = zdata,
                graphType                = "Scatter2D",
                colorBy                  = "Cluster",
                colorKey                 = dataobj$colorKey,
                title                    = "Cluster Overview",
                titleScaleFontFactor     = 0.6,
                xAxisTitle               = "UMAP1",
                yAxisTitle               = "UMAP2",
                xAxisMinorTicks          = FALSE,
                yAxisMinorTicks          = FALSE,
                axisTitleScaleFontFactor = 1.5,
                legendOrder              = list("Cluster" = dataobj$meta$clusters),
                transparency             = 0.5,
                dataPointSize            = 9,
                scatterOutlineThreshold  = 0,
                broadcast                = FALSE,
                zoomDisable              = TRUE,
                disableWheel             = TRUE,
                noValidate               = TRUE,
                saveFilename             = paste0(base_file_name, "_Overview"),
                events                   = events,
                disableToolbar           = report_modus,
                disableTouchToolbar      = report_modus
            )
        } else {
            warning('REDUCTION data not present')
        }
    }
    plot
}

get_scatter_panel_plot <- function(dataobj, genelist, gene_signatures, cluster, base_file_name,
                                   panel_plot_columns = 0, report_modus = FALSE) {
    if (g_debug) message('function: ', 'get_scatter_panel_plot')

    plot  <- NULL

    genelist <- get_total_genelist(dataobj$genes, genelist, gene_signatures)
    gene_count <- length(genelist)
    if (gene_count >= 1) {
        cells    <- get_cluster_filtered_cells(dataobj, cluster)
        data     <- get_gene_filtered_expression(dataobj, genelist)

        if (!is.null(data) && NROW(data) > 0) {
            exp_data  <- data %>%
                right_join(cells %>% select(Cell, Cluster, matches("red_axis")), by = "Cell", copy = TRUE) %>%
                select(-Cell) %>%
                gather("Gene", "Expression", -matches("red_axis"), -Cluster) %>%
                arrange(Expression)

            ydata <- exp_data %>%
                select(matches("red_axis")) %>%
                as.data.frame()
            rownames(ydata) <- make.names(1:NROW(ydata))

            zdata <- exp_data %>%
                select(Gene, Expression, Cluster) %>%
                as.data.frame()
            rownames(zdata) <- rownames(ydata)

            layout <-  get_layout_topology(genelist, panel_plot_columns)
            events <- JS("{'mousemove' : function(o, e, t) {
                            if (o != null) {
                                if (o.objectType == null) {
                                     t.showInfoSpan(e, '<b>' + o.z.Gene + '</b><br/>' +
                                                       '<i>Cluster: </i>' + o.z.Cluster + '<br>' +
                                                       'Expression: ' + Math.round(o.z.Expression*10000)/10000);
                                } else {
                                    t.showInfoSpan(e, o.display);
                            }; }}}")

            plot <- canvasXpress(
                data                    = ydata,
                varAnnot                = zdata,
                graphType               = "Scatter2D",
                colorBy                 = list("Expression"),
                colorSpectrum           = list("#E0E0E0", "darkred"),
                title                   = get_scatter_plot_title(cluster),
                titleScaleFontFactor    = 0.6,
                xAxisTitle              = "UMAP1",
                yAxisTitle              = "UMAP2",
                xAxisMinorTicks         = FALSE,
                yAxisMinorTicks         = FALSE,
                dataPointSize           = 14,
                scatterOutlineThreshold = 0,
                segregateVariablesBy    = list("Gene"),
                layoutTopology          = layout,
                broadcast               = FALSE,
                zoomDisable             = TRUE,
                disableWheel            = TRUE,
                noValidate              = TRUE,
                saveFilename            = paste0(base_file_name, "_Scatter"),
                events                  = events,
                disableToolbar          = report_modus,
                disableTouchToolbar     = report_modus
            )
        }
    }
    plot
}

get_violin_panel_plot <- function(dataobj, genelist, gene_signatures, base_file_name, panel_plot_columns = 0, report_modus = FALSE) {
    if (g_debug) message('function: ', 'get_violin_panel_plot')

    plot  <- NULL

    genelist <- get_total_genelist(dataobj$genes, genelist, gene_signatures)
    gene_count <- length(genelist)

    if (gene_count >= 1) {
        data <- get_gene_filtered_expression(dataobj, genelist)

        if (!is.null(data) && NROW(data) > 0) {
            exp_data <- data %>%
                right_join(dataobj$cells %>% select(Cell, Cluster), by = "Cell", copy = TRUE) %>%
                select(-Cell) %>%
                arrange(match(Cluster, dataobj$meta$clusters))

            ydata <- t(exp_data %>% select(-Cluster)) %>% as.data.frame()

            xdata <- t(exp_data %>% select(Cluster)) %>% as.data.frame(stringsAsFactors = F)
            zdata <- data.frame("Gene" = colnames(exp_data %>% select(-Cluster)), stringsAsFactors = F)
            rownames(zdata) <- rownames(ydata)
            layout <-  get_layout_topology(genelist, panel_plot_columns)
            events <- JS("{ 'mousemove' : function(o, e, t) {
                         if (o != null) {
                            if (o.objectType == null) {
                                t.showInfoSpan(e, '<b>' + o.w.smps + ': ' + o.w.vars + '</b><br/>' +
                                               'Q1 = '   + Math.round(o.w.qtl1*10000)/10000 + '<br/>' +
                                               'Mean = ' + Math.round(o.w.mean*10000)/10000 + '<br/>' +
                                               'Q3 = '   + Math.round(o.w.qtl3*10000)/10000);
                            } else {
                                t.showInfoSpan(e, o.display);
                         }; }}}")

            plot <- canvasXpress(
                data                         = ydata,
                smpAnnot                     = xdata,
                varAnnot                     = zdata,
                graphType                    = "Boxplot",
                graphOrientation             = "vertical",
                colorBy                      = "Cluster",
                colorKey                     = dataobj$colorKey,
                showViolinBoxplot            = TRUE,
                violinTrim                   = "true",
                violinScale                  = "width",
                title                        = "Distribution by Gene",
                titleScaleFontFactor         = 0.6,
                showSampleNames              = FALSE,
                xAxisTitle                   = character(0),
                yAxisTitle                   = character(0),
                transparency                 = 0.5,
                showLegend                   = TRUE,
                groupingFactors              = list("Cluster"),
                segregateVariablesBy         = list("Gene"),
                layoutTopology               = layout,
                legendOrder                  = list("Cluster" = dataobj$meta$clusters),
                broadcast                    = FALSE,
                zoomDisable                  = TRUE,
                disableWheel                 = TRUE,
                noValidate                   = TRUE,
                afterRender                  = list(list("switchNumericToString", list("Cluster", TRUE))),
                showFunctionNamesAfterRender = FALSE,
                saveFilename                 = paste0(base_file_name, "_Violin"),
                events                       = events,
                disableToolbar               = report_modus,
                disableTouchToolbar          = report_modus
            )
        }
    }
    plot
}

get_heatmap_plot <- function(dataobj, genelist, gene_signatures, additional_genes, base_file_name,
                             report_modus = FALSE) {
    if (g_debug) message('function: ', 'get_heatmap_plot')

    plot  <- NULL

    genelist   <- get_total_genelist(dataobj$genes, genelist, gene_signatures, additional_genes)
    gene_count <- length(genelist)
    if (gene_count >= 1) {
        data <- get_gene_filtered_expression(dataobj, genelist)

        if (!is.null(data) && NROW(data) > 0) {
            ydata <- right_join(data, dataobj$cells %>% select(Cell, Cluster), by = "Cell", copy = TRUE) %>%
                mutate(total_exp = rowSums(select(., -Cell, -Cluster)),
                       Cluster = as.character(Cluster)) %>%
                arrange(match(Cluster, dataobj$meta$clusters), total_exp) %>%
                select(-total_exp)
            rownames <- ydata$Cell

            zdata <- ydata %>% select(Cluster) %>% as.data.frame()
            ydata <- ydata %>% select(-Cluster, -Cell) %>% as.data.frame()

            rownames(zdata) <- rownames(ydata) <- rownames

            overlay_scale_font_factor <- get_overlay_scale_font_factor(gene_count, as.character(unique(zdata$Cluster)))
            events <- JS("{ 'mousemove' : function(o, e, t) {
                         if (o != null) {
                              if (o.objectType == null) {
                                 t.showInfoSpan(e, '<b>' + o.z.Cluster + ': ' + o.y.smps[0] + '</b><br/>' +
                                 o.y.data[0]);
                              } else {
                                t.showInfoSpan(e, o.display);
                         }; } }}")

            convertClusters <- ifelse(!is.na(suppressWarnings(sum(as.numeric(dataobj$meta$clusters)))), list("Cluster"), list())

            plot <- canvasXpress(
                data                     = ydata,
                varAnnot                 = zdata,
                stringVariableFactors    = convertClusters,
                graphType                = "Heatmap",
                colorKey                 = dataobj$colorKey,
                colorSpectrum            = list("white", "darkred"),
                title                    = "Expression by Gene",
                titleScaleFontFactor     = 0.6,
                subtitleScaleFontFactor  = 0.45,
                xAxisTitle               = "Expression",
                varOverlays              = list("Cluster"),
                varOverlayProperties     = list(Cluster = list(position = "top", fontstyle = "bold", thickness = 25)),
                overlayScaleFontFactor   = overlay_scale_font_factor,
                showVariableNames        = FALSE,
                heatmapIndicatorPosition = "top",
                heatmapIndicatorWidth    = 400,
                broadcast                = FALSE,
                zoomDisable              = TRUE,
                disableWheel             = TRUE,
                noValidate               = TRUE,
                saveFilename             = paste0(base_file_name, "_Heatmap"),
                events                   = events,
                disableToolbar           = report_modus,
                disableTouchToolbar      = report_modus
            )
        }
    }
    return(list(plot, gene_count))
}

get_differential_scatter_plot <- function(dataobj, cluster, plottitle, selectedInputID, file_name) {
    if (g_debug) message('function: ', 'get_differential_scatter_plot')

    plot  <- NULL
    if (!is.null(dataobj)) {
        ydata <- dataobj$cells %>%
            select(Cell, Cluster, matches("red_axis")) %>%
            rename(PC1 = red_axis1,
                   PC2 = red_axis2)
        rownames(ydata) <- ydata$Cell

        zdata <- dataobj$cells %>%
            mutate(Cluster = as.character(Cluster)) %>%
            select(-CellType.select, #must be present
                   -matches("red_axis"), -matches("barcode"), -matches("CLID"), -matches("ident")) #others

        show_legend     <- TRUE

        if (cluster != "All") {
            ydata           <- ydata %>% dplyr_filter(Cluster == cluster)
            zdata           <- zdata %>% dplyr_filter(Cluster == cluster)
            show_legend     <- FALSE
        }

        ydata$Cluster <- NULL

        rownames(ydata) <- ydata$Cell; ydata$Cell <- NULL
        rownames(zdata) <- zdata$Cell; zdata$Cell <- NULL

        events <- JS(
                     paste0("{'mousemove' : function(o, e, t) {
                                  if (o != null) {
                                      if (o.objectType == null) {
                                          t.showInfoSpan(e, null);
                                      } else {
                                          t.showInfoSpan(e, o.display);
                                      };
                                  }
                              },
                              'select': function(o, e, t) {
                                   if (CanvasXpress.selector.selections > 0) {
                                       Shiny.onInputChange('", selectedInputID, "', Object.keys(CanvasXpress.selector.vars));
                                   } else {
                                       Shiny.onInputChange('", selectedInputID, "', null);
                                   };
                              } }"))

        plot <- canvasXpress(
            data                     = ydata,
            varAnnot                 = zdata,
            graphType                = "Scatter2D",
            colorBy                  = "Cluster",
            colorKey                 = dataobj$colorKey,
            legendOrder              = list("Cluster" = dataobj$meta$clusters),
            title                    = plottitle,
            titleScaleFontFactor     = 0.6,
            xAxisMinorTicks          = FALSE,
            yAxisMinorTicks          = FALSE,
            xAxisTitle               = "UMAP1",
            yAxisTitle               = "UMAP2",
            axisTitleScaleFontFactor = 1.5,
            transparency             = 0.5,
            dataPointSize            = 14,
            scatterOutlineThreshold  = 0,
            showLegend               = show_legend,
            broadcast                = FALSE,
            zoomDisable              = TRUE,
            disableWheel             = TRUE,
            noValidate               = TRUE,
            saveFilename             = file_name,
            events                   = events,
            disableTouchToolbar      = TRUE
        )
    }
    plot
}

get_scatter_plot_title <- function(cluster) {
    if (g_debug) message('function: ', 'get_scatter_plot_title')

    title <- ""
    if (cluster == "All") {
        title <- "Expression by Gene - All Cells"
    } else {
        title <- paste("Cluster", cluster, "Expression by Gene")
    }
    title
}

get_overlay_scale_font_factor <- function(gene_count, cluster_names) {
    if (g_debug) message('function: ', 'get_overlay_scale_font_factor')

    result <- 1
    if (gene_count > 10) {
        cluster_name_length <- round(mean(nchar(cluster_names)))
        result <- 1 - 0.0012 * gene_count * cluster_name_length
    }
    if (result < 0.2) {
        result <- 0.2
    }
    result
}

get_co_exp_plot <- function(dataobj,
                            cells,
                            geneX_data,
                            geneY_data,
                            file_base_name,
                            xGate,
                            yGate,
                            noise,
                            pearson) {
    if (g_debug) message('function: ', 'get_co_exp_plot')

    plot  <- NULL
    percentages <- data.frame()

    cells_count <- NROW(geneX_data)

    if (cells_count >= 1) {
        xGate <- as.numeric(xGate)
        yGate <- as.numeric(yGate)
        geneX_col_names <- colnames(geneX_data)
        geneY_col_names <- colnames(geneY_data)
        geneX <- geneX_col_names[geneX_col_names != "Cell"]
        geneY <- geneY_col_names[geneY_col_names != "Cell"]

        # if user selected the same gene for both input
        if (geneX == geneY) {
            geneX <- paste0(geneX, ".x")
            geneY <- paste0(geneY, ".y")
        }

        gene_data <- geneX_data %>%
            left_join(geneY_data, by = "Cell")

        gene_data$noisyX <- gene_data[, geneX, drop = T] + rnorm(cells_count, mean = 0, sd = noise)
        gene_data$noisyY <- gene_data[, geneY, drop = T] + rnorm(cells_count, mean = 0, sd = noise)

        exp_data  <- gene_data %>%
            right_join(cells, by = "Cell") %>%
            select(-Cell,
                   -CellType.select,
                   -matches("red_axis"),
                   -matches("barcode"),
                   -matches("CLID"),
                   -matches("ident"))

        var_data <- exp_data %>% select(-noisyX, -noisyY)
        cx_data  <- exp_data %>% select(noisyX, noisyY)

        x_above_y_above <- exp_data %>%
            dplyr_filter(.data[[geneX]] >= xGate, .data[[geneY]] > yGate) %>%
            NROW()
        x_above_y_below <- exp_data %>%
            dplyr_filter(.data[[geneX]] >= xGate, .data[[geneY]] <= yGate) %>%
            NROW()
        x_below_y_above <- exp_data %>%
            dplyr_filter(.data[[geneX]] <= xGate, .data[[geneY]] > yGate) %>%
            NROW()
        x_below_y_below <- exp_data %>%
            dplyr_filter(.data[[geneX]] <= xGate, .data[[geneY]] <= yGate) %>%
            NROW()

        top_right_quadrant    <- round((x_above_y_above / cells_count) * 100, 2) %>% format(nsmall = 2)
        bottom_right_quadrant <- round((x_above_y_below / cells_count) * 100, 2) %>% format(nsmall = 2)
        top_left_quadrant     <- round((x_below_y_above / cells_count) * 100, 2) %>% format(nsmall = 2)
        bottom_left_quadrant  <- round((x_below_y_below / cells_count) * 100, 2) %>% format(nsmall = 2)

        percentages <- data.frame(geneX = c("TRUE", "TRUE", "FALSE", "FALSE"),
                                  geneY = c("TRUE", "FALSE", "TRUE", "FALSE"),
                                  Frequency	= c(x_above_y_above,
                                                x_above_y_below,
                                                x_below_y_above,
                                                x_below_y_below),
                                  Percentage = c(top_right_quadrant,
                                                 bottom_right_quadrant,
                                                 top_left_quadrant,
                                                 bottom_left_quadrant) %>% as.numeric())
        colnames(percentages) <- c(geneX, geneY, "Frequency", "Percentage")

        xMax <- cx_data %>%
            select(noisyX) %>%
            summarise(max(.)) %>%
            pull()
        yMax <- cx_data %>%
            select(noisyY) %>%
            summarise(max(.)) %>%
            pull()

        ## calculate offset for placing quarters text
        noise_decimals <- decimal_places(noise)
        ## try to update the text to be just beside noise points and not too far
        if (noise_decimals > 1) {
            noise_offset <- ceiling_decimals(noise, noise_decimals - 1)
        } else {
            noise_offset <- ceiling_decimals(noise, noise_decimals + 1)
        }

        if (!is.null(exp_data) && NROW(exp_data) > 0) {
            events <- JS(
                glue("{'mousemove' : function(o, e, t) {
                            if (o != null) {
                                if (o.objectType == null) {
                                     gates = o.z;
                                     geneX = '{{geneX}}';
                                     geneY = '{{geneY}}';
                                     gene1_data = parseFloat(gates[geneX][0]).toFixed(4);
                                     gene2_data = parseFloat(gates[geneY][0]).toFixed(4);
                                     t.showInfoSpan(e,
                                         '<b>' + o.z.Cluster[0] + '</b><br/>' +
                                         '<b>' + geneX + '</b>: ' + gene1_data + '<br/>' +
                                         '<b>' + geneY + '</b>: ' + gene2_data + '<br/>');
                                } else {
                                    t.showInfoSpan(e, o.display);
                            }; }}}", .open = "{{", .close = "}}"))

            title     <- glue("{geneX} vs {geneY} Expression")
            file_name <- glue("{file_base_name}_{geneX}_{geneY}_CoExpression")
            subtitle  <- glue("Pearson Correlation: {round(pearson, 5)}")
            min_value <- min(-0.1, noise * -10)

            plot <- canvasXpress(
                data                    = cx_data,
                varAnnot                = var_data,
                xAxisTitle              = geneX,
                yAxisTitle              = geneY,
                graphType               = "Scatter2D",
                colorBy                 = "Cluster",
                colorKey                = dataobj$colorKey,
                legendOrder             = list("Cluster" = dataobj$meta$clusters),
                transparency            = 0.8,
                title                   = title,
                subtitle                = subtitle,
                titleScaleFontFactor    = 0.8,
                subtitleScaleFontFactor = 0.45,
                scatterOutlineThreshold = 0,
                broadcast               = FALSE,
                zoomDisable             = TRUE,
                noValidate              = TRUE,
                events                  = events,
                setMinX                 = min_value,
                setMinY                 = min_value,
                fixedAspectRatio        = 1,
                decorations             = list(
                    line = list(list(color = "red",
                                     width = 2,
                                     x     = xGate,
                                     type  = "DashedLine"),
                                list(color = "blue",
                                     width = 2,
                                     y     = yGate,
                                     type  = "DashedLine")),
                    text = list(list(label = glue("{top_right_quadrant}%"),
                                     x     = floor_decimals(xMax) - noise_offset,
                                     y     = floor_decimals(yMax) - noise_offset),
                                list(label = glue("{bottom_right_quadrant}%"),
                                     x     = floor_decimals(xMax) - noise_offset,
                                     y     = ceiling_decimals(min_value)),
                                list(label = glue("{top_left_quadrant}%"),
                                     x     = ceiling_decimals(min_value),
                                     y     = floor_decimals(yMax) - noise_offset),
                                list(label = glue("{bottom_left_quadrant}%"),
                                     x     = ceiling_decimals(min_value),
                                     y     = ceiling_decimals(min_value)))),
                saveFilename            = file_name
            )
        }
    }
    list(plot = plot, percentages = percentages, geneX = geneX, geneY = geneY)
}
