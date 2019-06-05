# ----------------------
# Data Related Functions
# ----------------------
library(future)
plan(multiprocess)
options(future.globals.maxSize = +Inf)


check_required_field <- function(error_messages, field, field_name, check_type = "null") {
    if (check_type == "null" && is.null(field) ||
        check_type == "empty" && field == "") {
        if (is.null(error_messages[[g_missing_fields]])) {
            error_messages[[g_missing_fields]] <- field_name
        } else{
            error_messages[[g_missing_fields]] <- c(error_messages[[g_missing_fields]], field_name)
        }
    }
    error_messages
}

# Load the data from given file and check the content
load_data <- function(file) {
    result         <- NULL
    object         <- NULL
    error_messages <- list()

    tryCatch({
        object         <- readRDS(file$datapath)
    }, error = function(e) {
        if (e$message == "unknown input format") {
            error_messages[[g_error_field]] <<- get_file_upload_error_message(file, "is not in the correct format.")
        } else {
            error_messages[[g_error_field]] <<- get_file_upload_error_message(file, e$message)
        }
    })

    if (is.null(error_messages[[g_error_field]])) {
        if (is.null(object)) {
            error_messages[[g_error_field]] <- get_file_upload_error_message(file, "is empty.")
        }
        else if (tolower(class(object)) != "seurat") {
            error_messages[[g_error_field]] <- get_file_upload_error_message(file, "should be a Seurat object.")
        }
        else if (substr(object@version, 1, 1) != "3") {
            error_messages[[g_error_field]] <- get_file_upload_error_message(file, "should be a Seurat version 3 object.")
        }
        else {
            meta_info         <- object@misc$meta.info
            meta_fields_top   <- c("Title", "Authors", "Publication", "Summary")
            meta_fields_table <- setdiff(names(meta_info), meta_fields_top)
            meta_info_top     <- meta_info[meta_fields_top]

            # Check required fields
            meta_fields_required <- meta_fields_top[1]
            for (field in meta_fields_required) {
                error_messages <- check_required_field(error_messages, object@misc$meta.info[field], paste(field, paste0("(@misc$meta.info$", field, ")")))
            }
            for (field in c("top10", "top30")) {
                error_messages    <- check_required_field(error_messages, object@misc$DE[field], paste(field, paste0("(@misc$DE$", field, ")")))
            }
            error_messages    <- check_required_field(error_messages, colnames(object), "Cellnames (colnames(object))", check_type = "empty")
            error_messages    <- check_required_field(error_messages, Embeddings(object, reduction = "tsne"), "Dimensionality Reduction (Embeddings(object, reduction = 'tsne'))")
            error_messages    <- check_required_field(error_messages, object@active.ident, "Ident (@ident)")
            error_messages    <- check_required_field(error_messages, object@active.assay, "Assay (@assay)")
            error_messages    <- check_required_field(error_messages, GetAssayData(object), "Data (GetAssayData(object))")

            # Data segregation fields should match with meta_data
            data_segregation <- object@misc$DataSegregation
            if (!is.null(data_segregation)) {
                fields <- names(data_segregation)
                meta_fields <- colnames(object@meta.data)
                if (!identical(setdiff(fields, meta_fields), character(0))) {
                    error_messages[[g_error_field]] <- get_file_upload_error_message(file, "Data segregation fields doesn't match with meta data fields.")
                } else {
                    result$filter_list <- data_segregation
                    result$metadata    <- object@meta.data[, fields, drop = FALSE]
                }
            }

            if (length(error_messages) == 0) {
                result$seurat     <- object
                result$tsne       <- as.data.frame(Embeddings(object, reduction = "tsne"))
                result$tsne$Cell  <- rownames(result$tsne)

                result$cells      <- colnames(object)
                result$genes      <- sort(rownames(object))
                result$clusters   <- data.frame("Cluster" = as.factor(object@active.ident),
                                                "Cell"    = result$cells,
                                                stringsAsFactors = F)

                result$clustOrder <- unique(result$clusters$Cluster)
                if (is.factor(result$clustOrder)) {
                    result$clustOrder <- levels(result$clustOrder)
                }

                result$expression <- future({avg.ex.scale(object) %>% as.data.frame()}, stdout = FALSE)
                result$detection  <- future({local_AverageDetectionRate(object) %>% as.data.frame()}, stdout = FALSE)

                result$meta       <- list(object_name = gsub("\\..*", "", file$name),
                                          title       = meta_info_top$Title,
                                          author      = meta_info_top$Authors,
                                          publication = meta_info_top$Publication,
                                          summary     = meta_info_top$Summary,
                                          cells       = length(result$cells),
                                          genes       = length(result$genes),
                                          clusters    = levels(result$clusters$Cluster),
                                          top10       = object@misc$DE$top10,
                                          top30       = object@misc$DE$top30)

                result$meta$table <- meta_info[meta_fields_table]
                result$meta$table <- result$meta$table[!is.na(result$meta$table) & result$meta$table != "NA" & result$meta$table != ""]

                # setup colors for the clusters globally
                nClusters <- nlevels(result$clusters$Cluster)
                set1_max_colors <- 9

                # if nClusters is larger than the max colors supported by the palette, interpolate colors
                if (nClusters > set1_max_colors) {
                    get_palette  <- colorRampPalette(brewer.pal(set1_max_colors, "Dark2"))
                    cluster.cols <- get_palette(nClusters)
                } else {
                    cluster.cols <- brewer.pal(nClusters, "Dark2")
                }
                result$colorKey <- list("Cluster" = list())
                for (l in levels(result$clusters$Cluster)) {
                    result$colorKey[["Cluster"]][l] = cluster.cols[which(levels(result$clusters$Cluster) == l)]
                }
            }
        }
    }
    if (length(error_messages) == 0) {
        result <- list(errors = NULL, object = result)
    } else {
        result <- list(errors = error_messages, object = NULL)
    }
    result
}


get_cluster_filtered_cells <- function(dataobj, cluster) {
    cells <- dataobj$cells
    if (cluster %in% levels(dataobj$clusters$Cluster)) {
        cells <- dataobj$clusters %>% filter(Cluster == cluster)
        cells <- cells$Cell
    }
    cells
}

get_filtered_cell_data <- function(data, row_ids = NULL, column_names = NULL) {
    if (is.null(column_names)) {
        column_names <- colnames(data)
        data$id <- seq(1, nrow(data))
        result <- data[data$id %in% row_ids, column_names, drop = FALSE]
    }
    else if (is.null(row_ids)) {
        result <- data[, column_names, drop = FALSE]
    } else {
        result <- data
    }
    result
}


get_gene_filtered_expression <- function(dataobj, genelist, cells) {
    object_data   <- as.matrix(GetAssayData(dataobj$seurat))
    genes.checked <- intersect(genelist, dataobj$genes)
    gene.data <- NULL

    if (length(genes.checked) > 0) {
        gene.data <- t(object_data[genes.checked, cells, drop = FALSE])
        gene.data <- as.data.frame(as.matrix(gene.data))
        gene.data$Cell <- rownames(gene.data)
    }

    gene.data
}


get_dotplot_data <- function(dataobj, genelist) {
    genes.checked <- intersect(genelist, dataobj$genes)
    result <- NULL

    expression <- future::value(dataobj$expression)
    detection  <- future::value(dataobj$detection)
    if (length(genes.checked) > 0 && !is.null(expression) && !is.null(detection)) {
        genes.exp <- expression[genes.checked, , drop = FALSE] %>%
            rownames_to_column("Gene") %>%
            gather(key = "Cluster", value = "mean_Exp", -Gene)

        genes.pct <- detection[genes.checked, , drop = FALSE] %>%
            rownames_to_column("Gene") %>%
            gather(key = "Cluster", value = "pct_Exp", -Gene) %>%
            mutate(pct_Exp = pct_Exp * 100)

        result <- left_join(genes.exp, genes.pct, by = c("Gene", "Cluster"))
    }

   result
}


calculate_differentials <- function(dataobj, cluster1, cells1, cluster2, cells2) {
    differentials <- NULL

    # both clusters must be chosen
    # All cannot be compared to another cluster
    # Cells must be selected if clusters are the same
    if (is.null(cluster1) || is.null(cluster2) ||
        ((cluster1 != cluster2) && ((cluster1 == "All") || cluster2 == "All")) ||
        ((cluster1 == cluster2) && (is.null(cells1) || is.null(cells2)))) {
            return(differentials)
        }

    if (cluster1 != cluster2) {
        # cells, even if selected, will be ignored
        # this is a full-cluster comparison
        cells1 <- get_cluster_filtered_cells(dataobj, cluster1)
        cells2 <- get_cluster_filtered_cells(dataobj, cluster2)
    }
    else {
        # same cluster, cells determined by selections
        if (is.null(cells1) || length(cells1) < 1 ||
            is.null(cells2) || length(cells2) < 1 ||
            (all(cells1 %in% cells2) && all(cells2 %in% cells1))) {
            return(differentials)
        }
        else {
            # remove any cells that are selected in both clusters
            cells1.a <- cells1[!(cells1 %in% cells2)]
            cells2.a <- cells2[!(cells2 %in% cells1)]
            cells1 <- cells1.a
            cells2 <- cells2.a
            rm(cells1.a, cells2.a)
        }
    }
    # check that we have the minimum required cells
    if (length(cells1) < g_differential_min_no_cells || length(cells2) < g_differential_min_no_cells) {
        return(differentials)
    }

    # check that we have at least 1 comparable set
    if (all(cells1 %in% cells2) || all(cells2 %in% cells1)) {
        return(differentials)
    }

    differentials <- Fast.WC.DE(dataobj$seurat, cells.1 = cells1, cells.2 = cells2,
                                logfc.threshold = g_differential_logfc_threshold,
                                num.gene2test   = g_differential_gene_threshold,
                                min.pct.diff    = g_differential_pct_threshold)

    if (!is.null(differentials)) {
        if (cluster1 == cluster2) {
            colnames(differentials)[which(colnames(differentials) == "pct.1")] <- paste0("pct_expr_", cluster1, "_sel1")
            colnames(differentials)[which(colnames(differentials) == "pct.2")] <- paste0("pct_expr_", cluster2, "_sel2")
        }
        else {
            colnames(differentials)[which(colnames(differentials) == "pct.1")] <- paste0("pct_expr_", cluster1)
            colnames(differentials)[which(colnames(differentials) == "pct.2")] <- paste0("pct_expr_", cluster2)
        }

        differentials <- differentials %>%
            mutate(p_val     =  signif(p_val, 4),
                   p_val_adj = signif(p_val_adj, 4))
    }
    differentials
}

# scale average expression data
avg.ex.scale <- function(seurat_obj) {
    avg.exp <- AverageExpression(seurat_obj, assays = seurat_obj@active.assay, verbose = FALSE)
    exp.scale <- t(x = scale(x = t(x = avg.exp[[seurat_obj@active.assay]])))
    exp.scale <- MinMax(data = exp.scale, max = 2.5, min = (-1) * 2.5)
    return(exp.scale)
}

get_file_upload_error_message <- function(file, message) {
    paste0("<h3>The uploaded file (", file$name, paste(")", message), "</h3><br>")
}

get_file_dialog_error_message <- function(file, error_messages) {
    if (length(error_messages) > 0) {
        error          <- error_messages[[g_error_field]]
        missing_fields <- error_messages[[g_missing_fields]]
        if (!is.null(error)) {
            error
        }
        else if (length(missing_fields) > 0) {
            header  <- get_file_upload_error_message(file, "is missing required fields:")
            content <- paste0("<ul id='missing_list'>", paste(paste0("<li>", missing_fields, "</li>"), collapse = ""), "</ul>")
            paste0(header, content)
        }
    }
}

get_filter_prefix <- function(objectId) {
    paste0(objectId, "filter")
}

get_filter_input_Fields <- function(input_fields, objectId) {
    input_fields[grepl(paste0("^", get_filter_prefix(objectId)), input_fields)]
}

local_AverageDetectionRate <- function(seurat_obj, thresh.min = 0) {
    ident.use <- seurat_obj@active.ident
    data.all <- data.frame(row.names = rownames(x = GetAssay(seurat_obj, seurat_obj@active.assay)))
    for (i in sort(x = unique(x = ident.use))) {
        temp.cells <- WhichCells(object = seurat_obj, ident = i)
        data.temp <- apply(
            X = GetAssay(seurat_obj, seurat_obj@active.assay)[, temp.cells],
            MARGIN = 1,
            FUN = function(x) {
                return(sum(x > thresh.min)/length(x = x))
            }
        )
        data.all <- cbind(data.all, data.temp)
        colnames(x = data.all)[ncol(x = data.all)] <- i
    }
    colnames(x = data.all) <- sort(x = unique(x = ident.use))
    return(data.all)
}
