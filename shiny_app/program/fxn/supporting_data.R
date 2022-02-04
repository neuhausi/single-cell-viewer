suppressPackageStartupMessages({
    library(tictoc)
    library(httr) #temporary due to an issue with revealsc v1.0.4
})

# ----------------------
# Data Related Functions
# ----------------------
options(dplyr.summarise.inform = FALSE)

load_data <- function(parms) {
    if (g_debug) message('function: ', 'load_data')

    result         <- NULL
    error_messages <- list()
    src_data.all   <- NULL
    object_name    <- NULL

    tryCatch({
        # NOTE: revealsc is large and messy so we are going to only call
        # it with :: throughout this app to avoid many of the conflicts
        # introduced by sourcing via library() or require()

        if (!requireNamespace('revealsc')) {
            # package is not installed
            stopApp('the revealsc package must be installed/available')
        } else {
            # revealsc has a bad call to as.scidb that isn't double-colon
            # declared so we need to actually attach scidb at this point
            suppressPackageStartupMessages({
                require(scidb)
            })
            message("scidb    = v", packageVersion('scidb'))
            message("revealsc = v", packageVersion('revealsc'))
        }

        try({
            env_user <- Sys.getenv("SHINY_SCV_SCIDB_USER")
            env_pass <- Sys.getenv("SHINY_SCV_SCIDB_PASS")

            if (is.null(env_pass) || is.null(env_user) ||
                any("" %in% c(env_pass, env_user))) {
                stop('Cannot retrieve credentials to connect to SCIDB, username or token are missing.')
            }
        })

        server <- Sys.getenv("SHINY_SCV_SCIDB_HOST")
        message('server   = ', server)
        if (is.null(server) || is.na(server) || length(server) < 1) {
            stop('API server hostname is not defined for SCIDB')
        }

        g_api_con <<- revealsc::revealsc_connect(host              = server,
                                                 username          = env_user,
                                                 token             = env_pass,
                                                 result_size_limit = 16*1048,
                                                 local             = FALSE)

        if (!('isetID' %in% names(parms))) {
            # project/sample being used

            if (!('projectID' %in% names(parms)) || !('sampleID' %in% names(parms))) {
                error_messages <- append(error_messages,
                                         'projectID and sampleID must both be given as URL parameters')
            } else {
                g_api_pipelines <<- revealsc::get_matrix_counts(con = g_api_con, by = 'Sample', return_db_indices = TRUE)
                src_data.all    <- get_api_data(parms)
                object_name     <- paste0("Project: ", parms$projectID, " and Sample: ", parms$sampleID)
            }
        } else {
            # iset being used

            # ensure someone hasn't sent all/conflicting parameters
            if (('projectID' %in% names(parms)) || ('sampleID' %in% names(parms))) {
                error_messages <- append(error_messages,
                                         'isetID conflicts with specifying (projectID + sampleID) as URL parameters')
            } else {
                g_api_pipelines <<- revealsc::get_matrix_counts(con = g_api_con, by = 'iSet', return_db_indices = TRUE)
                src_data.all    <- get_api_data(parms)
                object_name     <- paste0("iSet: ", parms$isetID)
            }
        }

        error_messages <- append(error_messages,
                                 check_required_fields(src_data.all))

        if (length(error_messages) == 0) {
            result       <- list()

            #ensure all genes have a symbol and gene Symbols are unique
            result$genes <- data.frame(gene_id = src_data.all$genes$feature_id,
                                       Symbol = ifelse(is.na(src_data.all$genes$gene_symbol),
                                                       src_data.all$genes$name,
                                                       src_data.all$genes$gene_symbol),
                                       stringsAsFactors = F) %>% arrange(Symbol)

            result$genes$Symbol <- make.names(result$genes$Symbol, unique = T)

            result$api_selector <- src_data.all$selector

            result$cells <- src_data.all$cells %>%
                left_join(src_data.all$red, by = "cellid") %>%
                mutate(Cell    = cellid) %>%
                select(Cell, everything(), -cellid)

            result$meta <- list(
                object_name = object_name,
                title       = object_name,
                cells       = NROW(result$cells),
                genes       = NROW(result$genes))

            if (!('iset' %in% names(src_data.all))) {
                if (NROW(src_data.all$sample$tags_df) < 1) {
                    result$meta$table <- data.frame()
                } else {
                    meta_tbl <- data.frame("field" = colnames(src_data.all$sample$tags_df),
                                           "value" = t(src_data.all$sample$tags_df)[, 1],
                                           stringsAsFactors = F)

                    meta_tbl <- lookup_meta_values(meta_tbl)

                    result$meta$table <- meta_tbl
                }
            } else {
                result$meta$table <- data.frame("field" = "Description",
                                                "value" = src_data.all$iset$description,
                                                stringsAsFactors = F) %>% as.data.frame()
            }
        } else {
            # log out error messages to console for potential debugging
            message('Errors loading data for (', paste(parms, collapse = ", "), '): ', paste(error_messages, collapse = ';'))
        }
    },
    error = function(e) {
        warning('Error in load_data: ', e)
    })
    result
}

calc_cluster_data <- function(cells) {
    result       <- list()

    cluster_order <- cells %>%
        select(Cell, Cluster) %>%
        group_by(Cluster) %>%
        summarise(cells = n()) %>%
        arrange(desc(cells)) %>%
        pull(Cluster)

    # setup colors for the clusters globally
    suppressWarnings({
        nClusters <- length(cluster_order)

        if (nClusters > 9) {
            # we must interpolate the colors
            get_palette  <- colorRampPalette(brewer.pal(9, "Set1"))
            cluster.cols <- get_palette(nClusters)
        } else {
            cluster.cols <- brewer.pal(nClusters, "Set1")
            cluster.cols <- cluster.cols[1:nClusters]
        }

        cluster.cols <- setNames(cluster.cols, cluster_order)
    })

    result$colorKey <- list("Cluster" = as.list(cluster.cols))
    result$clusters <- cluster_order
    result
}

get_api_data <- function(parms) {
    result.obj <- list()
    selector   <- list(projectid     = parms$projectID,
                       sampleid      = parms$sampleID,
                       isetname      = parms$isetID,
                       matrix_counts = "normalized_counts")

    tryCatch({
        tic(paste0('API retrieval'))

        if (is.null(selector$isetname)) {
            # using project/sample ids

            addl_selectors <- g_api_pipelines %>%
                dplyr_filter(projectid    == selector$projectid,
                             sampleid     == selector$sampleid,
                             matrix_count == selector$matrix_counts) %>%
                select(sample_id, measurementset_id, reference)

            if (NROW(addl_selectors) != 1) {
                errorCondition(paste0("Project/Sample <", selector$projectid, "/",
                                      selector$sampleid, " does not exist in pipelines table"))
            } else {
                selector$sample_id         <- as.numeric(addl_selectors$sample_id[1])
                selector$measurementset_id <- as.numeric(addl_selectors$measurementset_id[1])
                selector$reference         <- addl_selectors$reference[1]
            }

            result.obj          <- revealsc::get_project(con       = g_api_con,
                                                         projectid = selector$projectid)
            result.obj$samples  <- NULL
            result.obj$selector <- selector

            result.obj$sample   <- revealsc::get_sample(con       = g_api_con,
                                                        projectid = selector$projectid,
                                                        sampleid  = selector$sampleid)

            result.obj$cells   <- revealsc::get_cells(con       = g_api_con,
                                                      projectid = selector$projectid,
                                                      sampleid  = selector$sampleid,
                                                      silent    = T)

            result.obj$genes   <- revealsc:::.get_feature_ids_at_sample_matrix_count(con          = g_api_con,
                                                                                     projectid    = selector$projectid,
                                                                                     sampleid     = selector$sampleid,
                                                                                     matrix_count = selector$matrix_counts,
                                                                                     returnType   = "feature_metadata")

            result.obj$red     <- revealsc::get_reduction(con        = g_api_con,
                                                          projectid  = selector$projectid,
                                                          sampleid   = selector$sampleid,
                                                          reduction  = "umap",
                                                          returnType = 'data.frame', )
            colnames(result.obj$red) <- c("cellid", "red_axis1", "red_axis2")
        } else {
            # using iset ids

            addl_selectors <- g_api_pipelines %>%
                dplyr_filter(isetid       == selector$isetname,
                             matrix_count == selector$matrix_counts) %>%
                select(iset_id, measurementset_id, reference)

            if (NROW(addl_selectors) != 1) {
                errorCondition(paste0("iSet <", selector$isetname, " does not exist in pipelines table"))
            } else {
                selector$iset_id           <- as.numeric(addl_selectors$iset_id)
                selector$measurementset_id <- as.numeric(addl_selectors$measurementset_id[1])
                selector$reference         <- addl_selectors$reference[1]
            }

            result.obj$iset <- revealsc::get_iset(con    = g_api_con,
                                                      isetid = selector$isetname)

            # temporarily add the description which has to be retrieved differently
            desc <- revealsc::get_isets(con = g_api_con) %>%
                dplyr_filter(isetid == selector$isetname) %>%
                pull(description)

            result.obj$iset$description <- ifelse(length(desc) == 0, "", desc)

            result.obj$selector <- as.list(selector)

            result.obj$projects <- NULL
            for (p in result.obj$iset$projectid) {
                if (!is.null(p) && !is.na(p)) {
                    pinfo <- revealsc::get_project(con       = g_api_con,
                                                   projectid = p)
                    result.obj$projects <- bind_rows(result.obj$projects, pinfo$project)
                }
            }

            result.obj$cells <- revealsc::get_cells(con       = g_api_con,
                                                    isetid    = selector$isetname,
                                                    silent    = T)

            result.obj$genes <- revealsc:::.get_feature_ids_at_iset_matrix_count(con          = g_api_con,
                                                                                 isetid       = selector$isetname,
                                                                                 matrix_count = selector$matrix_counts,
                                                                                 returnType   = "feature_metadata")

            result.obj$red <- revealsc::get_reduction(con        = g_api_con,
                                                      isetid     = selector$isetname,
                                                      reduction  = "umap",
                                                      returnType = 'data.frame')

            colnames(result.obj$red) <- c("cellid", "red_axis1", "red_axis2")
        }

        toc(log = T)
        result.obj$api_load <- gsub(' elapsed', '', unlist(tic.log()))
        tic.clearlog()
    },
    warning = function(w) {
        message('Warning in get_api_data: ', w)
    },
    error = function(e) {
        warning(paste0('Error in get_api_data: ', e))
    })

    result.obj
}

check_required_fields <- function(datalist) {
    msgs <- list()
    toplevel.req <- c("cells", "genes", "red")

    fields.req <- list(
        "project" = c("projectid"),
        "sample"  = c("sampleid", "cell_count", "tags_df"),
        "iset"    = c("iset_id"),
        "cells"   = c("cellid", "CellType.select"),
        "genes"   = c("feature_id", "gene_symbol", "name"),
        "red"     = c("cellid", "red_axis1", "red_axis2"))

    exp.req <- c("cell_id", "feature_id", "value")

    tryCatch({
        if (is.null(datalist)) {
            msgs <- append(msgs,
                           'Data was not able to be retrieved, is NULL')
        }

        if (any(!(toplevel.req %in% names(datalist)))) {
            msgs <- append(msgs,
                           paste0('Missing top level item(s): ', paste(toplevel.req[!toplevel.req %in% names(datalist)], collapse = ", ")))
        }

        for (fieldtype in names(fields.req)) {
            if (!is.null(datalist[[fieldtype]]) &&
                any(!(fields.req[[fieldtype]] %in% names(datalist[[fieldtype]])))) {
                msgs <- append(msgs,
                               paste0('Missing ', fieldtype, ' field(s): ',
                                      paste(fields.req[[fieldtype]][!(fields.req[[fieldtype]] %in% names(datalist[[fieldtype]]))], collapse = ", ")))
            }
        }

        if (is.null(datalist$selector)) {
            msgs <- append(msgs, 'Missing api selector field')
        }
    },
    error = function(e) {
        warning('Error in check_required_fields: ', e)
    })
    msgs
}

#' load_clusters
#'     reads and checks the contents of the yaml file for valid clusters that would be used in the application
#'
#' @param filepath - path to the yaml file
#'
#' @return named list
load_clusters <- function(filepath) {

    default_clusters <- list("CellType.select" = "SingleR",
                             "CellType.author" = "Author",
                             "seurat_clusters" = "Seurat")

    allowed_clusters <- tryCatch({
        suppressWarnings(yaml::read_yaml(filepath))
    },
    error = function(e) {
        message(glue::glue("{e}\n"))
        default_clusters})

    #check if the clusters are valid
    clusters_are_character  <- unlist(lapply(allowed_clusters, is.character))
    clusters_no_space       <- !grepl("\\s", names(allowed_clusters))

    if (all(clusters_are_character, clusters_no_space)) {
        return(allowed_clusters)
    }

    message("The clusters specified are not allowed. Please check if the text is formatted properly. Using the default configuration.")
    default_clusters

}

get_gene_signatures <- function() {
    if (g_debug) message('function: ', 'get_gene_signatures')

    file_name     <- "signatures.csv"
    file_location <- paste0("program/data/", file_name)

    result <- list()
    tryCatch({
        lines <- readLines(file_location)
        for (line in lines) {
            elements <- trimws(gsub("\"", "", unlist(strsplit(line, ","))))
            elements_count <- length(elements)
            signature  <- elements[1]
            sig_genes <- list(elements[2:elements_count])
            names(sig_genes) <- signature
            result <- c(result, sig_genes)
        }
    },
    error = function(e) {
        warning('Error in get_gene_signatures: ', e)
    })
    result
}


get_cluster_filtered_cells <- function(dataobj, cluster, id_only = FALSE) {
    if (g_debug) message('function: ', 'get_cluster_filtered_cells')

    cells <- dataobj$cells
    if (cluster %in% dataobj$meta$clusters) {
        cells <- dataobj$cells %>%
            dplyr_filter(Cluster == cluster)
    }

    if (id_only) {
        cells <- cells %>%
            pull(Cell)
    }
    cells
}


get_gene_filtered_expression <- function(dataobj, genelist) {
    if (g_debug) message('function: ', 'get_gene_filtered_expression')

    genes.checked <- c()
    gene.data     <- NULL

    if (length(genelist) > 0) {
        genes.checked <- dataobj$genes %>%
            dplyr_filter(Symbol %in% genelist) %>%
            pull(gene_id) %>%
            as.numeric()
    }

    if (length(genes.checked) > 0) {
        expdata <- NULL

        if (!is.null(g_api_con) && !is.null(dataobj$api_selector)) {
            # load from API
            tic(paste0('Expression Retrieval for ', length(genes.checked), ' genes'))
            if (is.null(dataobj$api_selector$isetname)) {
                apidata <- revealsc:::.search_expression_by_sample_or_iset(
                    con               = g_api_con,
                    entity            = "Sample",
                    measurementset_id = dataobj$api_selector$measurementset_id,
                    sis_id            = dataobj$api_selector$sample_id,
                    feature_id        = genes.checked)
            } else {
                apidata <- revealsc:::.search_expression_by_sample_or_iset(
                    con               = g_api_con,
                    entity            = "iSet",
                    measurementset_id = dataobj$api_selector$measurementset_id,
                    sis_id            = dataobj$api_selector$iset_id,
                    feature_id        = genes.checked)
            }
            toc()

            expdata <- apidata %>%
                left_join(dataobj$genes, by = c("feature_id" = "gene_id")) %>%
                select(cell_id, Symbol, value)
        }

        other.cells <- dataobj$cells %>%
            dplyr_filter(!(Cell %in% expdata$cell_id)) %>%
            select(cell_id = Cell)

        gene.data <- expdata %>%
            mutate(Symbol = as.character(Symbol)) %>%
            pivot_wider(id_cols = c(cell_id), names_from = Symbol, values_from = value) %>%
            bind_rows(other.cells) %>%
            mutate_all(replace_na, replace = 0) %>%
            mutate(Cell = cell_id) %>%
            select(-cell_id)
    }

    gene.data
}

get_cell_filtered_expression <- function(dataobj, cells) {
    if (g_debug) message('function: ', 'get_cell_filtered_expression')

    cells.checked <- list()
    cells.data    <- NULL

    if (!is.null(cells) && length(cells) > 0) {
        cells.checked <- dataobj$cells %>%
            dplyr_filter(Cell %in% cells) %>%
            mutate(Cell = as.numeric(Cell)) %>%
            pull(Cell)
    }

    if (length(cells.checked) > 0) {
        expdata <- NULL

        if (!is.null(g_api_con) && !is.null(dataobj$api_selector)) {
            # load from API
            #TODO - remove tictoc
            tic(paste0('Expression Retrieval for ', length(cells.checked), ' cells'))

            apidata <- revealsc::search_expression(con          = g_api_con,
                                                   matrix_count = dataobj$api_selector$matrix_counts,
                                                   cellid       = cells.checked,
                                                   returnType   = 'data.frame_raw')

            toc()

            expdata <- apidata %>%
                left_join(dataobj$genes, by = c("feature_id" = "gene_id")) %>%
                select(cell_id, Symbol, value)

            cells.data <- expdata %>%
                mutate(Symbol = as.character(Symbol)) %>%
                pivot_wider(id_cols = c(cell_id), names_from = Symbol, values_from = value) %>%
                mutate_all(replace_na, replace = 0) %>%
                mutate(Cell = cell_id) %>%
                select(-cell_id)
        }
    }

    cells.data
}

calculate_differentials <- function(dataobj, cluster1, cells1, cluster2, cells2) {
    if (g_debug) message('function: ', 'calculate_differentials')

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
        cells1 <- get_cluster_filtered_cells(dataobj, cluster1, id_only = TRUE)
        cells2 <- get_cluster_filtered_cells(dataobj, cluster2, id_only = TRUE)
    } else {
        # same cluster, cells determined by selections
        if (is.null(cells1) || length(cells1) < 1 ||
            is.null(cells2) || length(cells2) < 1 ||
            (all(cells1 %in% cells2) && all(cells2 %in% cells1))) {
            return(differentials)
        } else {
            # remove any cells that are selected in both clusters
            cells1 <- cells1[!(cells1 %in% cells2)]
            cells2 <- cells2[!(cells2 %in% cells1)]
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
    #TODO - remove tictoc
    tic('Fast.DE2 Calculation: ')
    differentials <- Fast.DE2(dataobj,
                              cells.1         = cells1,
                              cells.2         = cells2,
                              logfc.threshold = g_differential_logfc_threshold,
                              num.gene2test   = g_differential_gene_threshold,
                              min.pct.diff    = g_differential_pct_threshold)
    toc()
    if (!is.null(differentials)) {
        if (cluster1 == cluster2) {
            colnames(differentials)[which(colnames(differentials) == "pct.1")] <- paste0("pct_expr_", cluster1, "_sel1")
            colnames(differentials)[which(colnames(differentials) == "pct.2")] <- paste0("pct_expr_", cluster2, "_sel2")
        } else {
            colnames(differentials)[which(colnames(differentials) == "pct.1")] <- paste0("pct_expr_", cluster1)
            colnames(differentials)[which(colnames(differentials) == "pct.2")] <- paste0("pct_expr_", cluster2)
        }

        differentials <- differentials %>%
            mutate(p_val     = signif(p_val, 4),
                   p_val_adj = signif(p_val_adj, 4))
    }
    differentials
}

lookup_meta_values <- function(metadatatbl) {
    result <- metadatatbl

    if (!is.null(result)) {
        if ("DOID" %in% result$field) {
            doid <- revealsc::search_ontology(ontology = 'DOID',
                                              id = result[result$field == "DOID", "value"])
            if (NROW(doid) > 0) {
                result[result$field == "DOID", "value"] <- doid$term
            }
        }

        if ("UBERONID" %in% result$field) {
            uberon <- revealsc::search_ontology(ontology = 'UBERON',
                                                id = result[result$field == "UBERONID", "value"])
            if (NROW(uberon) > 0) {
                result[result$field == "UBERONID", "value"] <- uberon$term
            }
        }
    }
    result
}

calculate_gene_x_gene_pearson_correlation <- function(geneX, geneY) {
    pearson <- NA
    tryCatch({
        if (is.null(geneX) || length(geneX) < 1) {
            return(list(result = pearson, error = "gene 1 data is null or empty."))
        }

        if (is.null(geneY) || length(geneY) < 1) {
            return(list(result = pearson, error = "gene 2 data is null or empty."))
        }

        if (length(geneX) != length(geneY)) {
            return(list(result = pearson, error = "gene 1 data length does not equal gene 2 data length."))
        }

        cof <- suppressWarnings(cor.test(x = geneX, y = geneY, method = "pearson"))
        list(result = cof$estimate[["cor"]], error = "")
    },

    error = function(e) {
        return(list(result = pearson, error = e$message))
    })
}

update_cluster_options <- function(session, object, cluster_column) {
    object$cells <- object$cells %>%
        mutate(Cluster = .data[[cluster_column]])

    cluster_data <- calc_cluster_data(object$cells)
    object$colorKey <- cluster_data$colorKey
    object$meta$clusters <- cluster_data$clusters


    cluster_options <- object$cells %>%
        select(Cluster) %>%
        unique() %>%
        pull()

    cluster_options <- c("All", cluster_options)
    update_clusters_inputs(session, object, cluster_options)
    object
}

update_clusters_inputs <- function(session, object, cluster_options) {
    updateSelectizeInput(session,
                         "scatterClusterSel",
                         choices  = cluster_options,
                         selected = "All",
                         server   = FALSE)

    updateSelectizeInput(session,
                         "coExpClusterSel",
                         choices  = cluster_options,
                         selected = "All",
                         server   = FALSE)

    updateSelectizeInput(session,
                         "gene1Sel",
                         choices  = object$genes$Symbol,
                         selected = character(0),
                         server   = TRUE)

    updateSelectizeInput(session,
                         "gene2Sel",
                         choices  = object$genes$Symbol,
                         selected = character(0),
                         server   = TRUE)

    for (item in c("differentialsCluster1Sel", "differentialsCluster2Sel")) {
        updateSelectizeInput(session,
                             item,
                             choices  = cluster_options,
                             selected = "",
                             server   = FALSE)
    }
}

open_cluster_change_modal <- function(session) {
    enable("cluster_proceed")
    enable("cluster_cancel")
    toggleModal(session, "cluster_change_modal", toggle = "open")
    runjs("$('#cluster_change_modal').one('shown.bs.modal', function(e) {
                                                                         Shiny.onInputChange('cluster_hide', false);
                                                                         $('#cluster_proceed').focus();});")
    runjs("$('#cluster_change_modal').one('hidden.bs.modal', function(e) {Shiny.onInputChange('cluster_hide', true);});")
}

