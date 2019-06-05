# ----------------------
# Plot Related Functions
# ----------------------
suppressPackageStartupMessages(library(htmlwidgets))

get_overview_plot <- function(dataobj, file_name) {
    plot <- NULL

    if (!is.null(dataobj) && nrow(dataobj$tsne) > 0) {
        ydata <- dataobj$tsne %>% select(-Cell)
        rownames(ydata) <- dataobj$tsne$Cell

        zdata <- dataobj$clusters %>% select(-Cell)
        zdata <- add_annotation_metadata(zdata, dataobj)
        rownames(zdata) <- dataobj$clusters$Cell

        events <- htmlwidgets::JS("{ 'mousemove' : function(o, e, t) {
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
            xAxisTitle               = "tSNE_1",
            yAxisTitle               = "tSNE_2",
            xAxisMinorTicks          = FALSE,
            yAxisMinorTicks          = FALSE,
            axisTitleScaleFontFactor = 1.5,
            legendOrder              = list("Cluster" = dataobj$clustOrder),
            transparency             = 0.5,
            dataPointSize            = 9,
            scatterOutlineThreshold  = 0,
            broadcast                = FALSE,
            zoomDisable              = TRUE,
            saveFilename             = file_name,
            events                   = events
        )
    }
    plot
}

get_scatter_panel_plot <- function(dataobj, genelist, cluster, file_name, panel_plot_columns = 0) {
    plot  <- NULL

    genelist <- get_total_genelist(dataobj$genes, genelist)
    gene_count <- length(genelist)
    if (gene_count >= 1) {
        cells    <- get_cluster_filtered_cells(dataobj, cluster)
        data     <- get_gene_filtered_expression(dataobj, genelist, cells)

        if (!is.null(data) && nrow(data) > 0) {
            exp_data  <- data %>% left_join(dataobj$tsne[cells, ], by = "Cell") %>%
                                  left_join(dataobj$clusters %>% filter(Cell %in% cells), by = "Cell") %>%
                                  select(-Cell) %>%
                                  gather("Gene", "Expression", -tSNE_1, -tSNE_2, -Cluster) %>%
                                  arrange(Expression)

            ydata <- exp_data %>% select(tSNE_1, tSNE_2)
            rownames(ydata) <- make.names(1:NROW(ydata))

            zdata <- exp_data %>% select(Gene, Expression, Cluster)
            rownames(zdata) <- rownames(ydata)
            layout <-  get_layout_topology(genelist, panel_plot_columns)
            events <-  htmlwidgets::JS("{'mousemove' : function(o, e, t) {
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
                xAxisTitle              = "tSNE_1",
                yAxisTitle              = "tSNE_2",
                xAxisMinorTicks         = FALSE,
                yAxisMinorTicks         = FALSE,
                dataPointSize           = 14,
                scatterOutlineThreshold = 0,
                segregateVariablesBy    = list("Gene"),
                layoutTopology          = layout,
                broadcast               = FALSE,
                zoomDisable             = TRUE,
                saveFilename            = file_name,
                events                  = events
            )
        }
    }
    plot
}


get_violin_panel_plot <- function(dataobj, genelist, file_name, panel_plot_columns = 0) {
    plot  <- NULL

    genelist <- get_total_genelist(dataobj$genes, genelist)
    gene_count <- length(genelist)

    if (gene_count >= 1) {
        data     <- get_gene_filtered_expression(dataobj, genelist, dataobj$cells)

        if (!is.null(data) && nrow(data) > 0) {
            exp_data  <- left_join(data, dataobj$clusters, by = "Cell") %>%
                select(-Cell) %>%
                arrange(Cluster, Cluster %in% c(dataobj$clustOrder))

            ydata <- t(exp_data %>% select(-Cluster)) %>% as.data.frame()

            xdata <- t(exp_data %>% select(Cluster)) %>% as.data.frame(stringsAsFactors = F)
            zdata <- data.frame("Gene" = colnames(exp_data %>% select(-Cluster)), stringsAsFactors = F)
            rownames(zdata) <- rownames(ydata)
            layout <-  get_layout_topology(genelist, panel_plot_columns)
            events <- htmlwidgets::JS("{ 'mousemove' : function(o, e, t) {
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
                legendOrder                  = list("Cluster" = dataobj$clustOrder),
                broadcast                    = FALSE,
                zoomDisable                  = TRUE,
                afterRender                  = list(list("switchNumericToString", list("Cluster", TRUE))),
                showFunctionNamesAfterRender = FALSE,
                saveFilename                 = file_name,
                events                       = events
            )
        }
    }
    plot
}


get_dot_plot <- function(dataobj, genelist, additional_genes, top_DEG_option, file_name) {
    plot  <- NULL

    genelist   <- get_total_genelist(dataobj$genes, genelist, additional_genes)
    gene_count <- length(genelist)

    if (gene_count >= 1) {
        data       <- get_dotplot_data(dataobj, genelist)

        if (!is.null(data)) {
            ydata.1 <- data %>%
                select(-pct_Exp) %>%
                spread(Cluster, mean_Exp) %>%
                select(c("Gene", dataobj$clustOrder)) %>%
                arrange(match(Gene, genelist))

            rownames(ydata.1) <- ydata.1$Gene
            ydata.1 <- ydata.1 %>% select(-Gene) %>% t()

            ydata.2 <- data %>%
                select(-mean_Exp) %>%
                spread(Cluster, pct_Exp) %>%
                select(c("Gene", dataobj$clustOrder)) %>%
                arrange(match(Gene, genelist))

            rownames(ydata.2) <- ydata.2$Gene
            ydata.2 <- ydata.2 %>% select(-Gene) %>% t()

            v.annot <- data.frame(Cluster = rownames(ydata.1), stringsAsFactors = F)
            rownames(v.annot) <- v.annot$Cluster

            overlay_scale_font_factor <- get_overlay_scale_font_factor(gene_count, v.annot$Cluster)
            titles <- get_dot_plot_titles(top_DEG_option)
            events <- htmlwidgets::JS("{ 'mousemove' : function(o, e, t) {
                          if (o != null) {
                              if (o.objectType == null) {
                                 t.showInfoSpan(e, '<b>' + o.z.Cluster + ': ' + o.y.smps[0] + '</b><br/>' +
                                 o.y.data2[0] + ' % Expressing' + '<br/>'+
                                 'Mean = ' + o.y.data[0]);
                               } else {
                                 t.showInfoSpan(e, o.display);
                           }; }}}")

            convertClusters <- ifelse(!is.na(suppressWarnings(sum(as.numeric(dataobj$clustOrder)))), list("Cluster"), list())

            plot <- canvasXpress(
                data                     = list(y = ydata.1, data2 = ydata.2),
                varAnnot                 = v.annot,
                stringVariableFactors    = convertClusters,
                graphType                = "Heatmap",
                sizeBy                   = "% Expressing",
                sizes                    = list(4, 8, 12, 16, 20, 24, 28, 32),
                sizeByData               = "data2",
                sizeByContinuous         = TRUE,
                colorSpectrum            = list("white", "darkred"),
                colorKey                 = dataobj$colorKey,
                title                    = titles[[1]],
                subtitle                 = titles[[2]],
                titleScaleFontFactor     = 0.6,
                subtitleScaleFontFactor  = 0.45,
                xAxisTitle               = "Mean Expression",
                varOverlays              = list("Cluster"),
                varOverlayProperties     = list(Cluster = list(position = "top", fontstyle = "bold", thickness = 25)),
                overlayScaleFontFactor   = overlay_scale_font_factor,
                showSmpDendrogram        = TRUE,
                showVariableNames        = FALSE,
                legendPosition           = "right",
                heatmapIndicatorPosition = "top",
                heatmapIndicatorWidth    = 400,
                broadcast                = FALSE,
                zoomDisable              = TRUE,
                saveFilename             = file_name,
                events                   = events
            )
        }
    }
    return(list(plot, gene_count))
}


get_heatmap_plot <- function(dataobj, genelist, additional_genes, top_DEG_option, file_name) {
    plot  <- NULL

    genelist   <- get_total_genelist(dataobj$genes, genelist, additional_genes)
    gene_count <- length(genelist)
    if (gene_count >= 1) {
        data       <- get_gene_filtered_expression(dataobj, genelist, dataobj$cells)

        if (!is.null(data) && nrow(data) > 0) {
            ydata <- left_join(data, dataobj$clusters, by = "Cell") %>%
                mutate(total_exp = rowSums(select(., -Cell, -Cluster))) %>%
                group_by(Cluster) %>%
                arrange(Cluster, desc(total_exp), .by_group = TRUE) %>%
                ungroup() %>%
                select(-total_exp)

            zdata <- ydata %>%
                select(Cluster)

            rownames(ydata) <- ydata$Cell
            rownames(zdata) <- ydata$Cell

            ydata <- ydata %>% select(-Cluster, -Cell)

            overlay_scale_font_factor <- get_overlay_scale_font_factor(gene_count, as.character(unique(zdata$Cluster)))
            titles <- get_heatmap_titles(top_DEG_option)
            events <- htmlwidgets::JS("{ 'mousemove' : function(o, e, t) {
                         if (o != null) {
                              if (o.objectType == null) {
                                 t.showInfoSpan(e, '<b>' + o.z.Cluster + ': ' + o.y.smps[0] + '</b><br/>' +
                                 o.y.data[0]);
                              } else {
                                t.showInfoSpan(e, o.display);
                         }; } }}")

            convertClusters <- ifelse(!is.na(suppressWarnings(sum(as.numeric(dataobj$clustOrder)))), list("Cluster"), list())

            plot <- canvasXpress(
                data                     = ydata,
                varAnnot                 = zdata,
                stringVariableFactors    = convertClusters,
                graphType                = "Heatmap",
                colorKey                 = dataobj$colorKey,
                colorSpectrum            = list("white", "darkred"),
                title                    = titles[[1]],
                subtitle                 = titles[[2]],
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
                saveFilename             = file_name,
                events                   = events
            )
        }
    }
    return(list(plot, gene_count))
}

get_differential_scatter_plot <- function(dataobj, cluster, plottitle, selectedInputID, file_name) {
    plot  <- NULL

    if (!is.null(dataobj)) {
        ydata         <- dataobj$tsne
        ydata$Cluster <- dataobj$clusters$Cluster
        zdata         <- dataobj$clusters

        if (cluster != "All") {
            ydata <- ydata %>%
                filter(Cluster == cluster) %>%
                column_to_rownames('Cell')
            zdata <- zdata %>%
                filter(Cluster == cluster) %>%
                column_to_rownames('Cell')
            show_legend <- FALSE
        } else {
            ydata           <- ydata %>% select(-Cell)
            zdata           <- zdata %>% select(-Cell)
            rownames(zdata) <- dataobj$clusters$Cell
            show_legend     <- TRUE
        }
        zdata  <- add_annotation_metadata(zdata, dataobj)
        events <- htmlwidgets::JS(paste0("{'mousemove' : function(o, e, t) {
                                    if (o != null) {
                                        if (o.objectType == null) {
                                            t.showInfoSpan(e, null);
                                        } else {
                                            t.showInfoSpan(e, o.display);
                                    }; }
                                },
                              'select': function(o, e, t) {
                                    if (typeof o === 'boolean' && o === false) {
                                        Shiny.onInputChange('", selectedInputID, "', null);
                                    } else {
                                        Shiny.onInputChange('", selectedInputID, "', t.selectDataPointObject.y.vars);
                                    };
                                }
                            }"))

        plot <- canvasXpress(
            data                     = ydata,
            varAnnot                 = zdata,
            graphType                = "Scatter2D",
            colorBy                  = "Cluster",
            colorKey                 = dataobj$colorKey,
            title                    = plottitle,
            titleScaleFontFactor     = 0.6,
            xAxisMinorTicks          = FALSE,
            yAxisMinorTicks          = FALSE,
            xAxisTitle               = "tSNE_1",
            yAxisTitle               = "tSNE_2",
            axisTitleScaleFontFactor = 1.5,
            transparency             = 0.5,
            dataPointSize            = 14,
            scatterOutlineThreshold  = 0,
            showLegend               = show_legend,
            broadcast                = FALSE,
            zoomDisable              = TRUE,
            saveFilename             = file_name,
            events                   = events,
            disableTouchToolbar      = TRUE
        )
    }
    plot
}

get_scatter_plot_title <- function(cluster) {
    title <- ""
    if (cluster == "All") {
        title <- "Expression by Gene - All Cells"
    } else {
        title <- paste("Cluster", cluster, "Expression by Gene")
    }
    title
}

get_dot_plot_titles <- function(top_DEG_option) {
    title <- "Mean and Percent Expression by Gene"
    subtitle <- get_top_DEG_subtitle(top_DEG_option)
    return(list(title, subtitle))
}

get_heatmap_titles <- function(top_DEG_option) {
    title <- "Expression by Gene"
    subtitle <- get_top_DEG_subtitle(top_DEG_option)
    return(list(title, subtitle))
}

get_top_DEG_subtitle <- function(top_DEG_option) {
    result <- ifelse(top_DEG_option == "off", "", paste("(Including Top", gsub("top", "", top_DEG_option), "Differentially Expressed Genes)"))
    result
}

get_overlay_scale_font_factor <- function(gene_count, cluster_names) {
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

add_annotation_metadata <- function(zdata, dataobj) {
    metadata <- dataobj$metadata
    if (!is.null(metadata) && nrow(metadata) > 0) {
        # zdata might be filtered, if so also filter metadata
        if (nrow(zdata) != nrow(metadata)) {
            metadata <- metadata[rownames(metadata) %in% rownames(zdata),, drop = FALSE]
        }
        zdata <- zdata %>% cbind(metadata)
    }
    zdata
}
