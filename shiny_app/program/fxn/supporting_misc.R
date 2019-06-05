# ----------------------
# Misc Functions
# ---------------------

check_panel_plot_columns <- function(genelist, panel_plot_columns) {
    is_valid <- TRUE
    value    <- 0

    # check if number
    if (suppressWarnings(!is.na(as.numeric(panel_plot_columns)))) {
        gene_count <- as.numeric(panel_plot_columns)
        if (gene_count <= 0 || round(gene_count) != gene_count) {
            is_valid <- FALSE
        } else {
            value <- gene_count
        }
    } else {
        if (panel_plot_columns != "auto") {
            is_valid <- FALSE
        }
    }
    return(list(is_valid, value))
}

get_layout_topology <- function(genelist, panel_plot_columns) {
    layout <-  FALSE
    if (panel_plot_columns > 0) {
        if (panel_plot_columns > length(genelist)) {
            panel_plot_columns <- length(genelist)
        }
        layout <- paste0(ceiling(length(genelist)/panel_plot_columns), "X", panel_plot_columns)
    }
    layout
}

get_total_genelist <- function(all_genes, selected_genes, additional_genes = NULL) {
    intersect(unique(c(selected_genes, additional_genes)), all_genes)
}

get_additional_genes <- function(add_genes_option, top10, top30) {
    additional_genes <- NULL
    if (add_genes_option == "top10") {
        additional_genes <- as.character(top10$Gene)
    } else if (add_genes_option == "top30") {
        additional_genes <- as.character(top30$Gene)
    }
    additional_genes
}

get_dynamic_plot_height <- function(gene_count) {
    600 + 20 * gene_count
}

get_top_DE_download_filename <- function(project, top_DE_name) {
    current_time <- format(Sys.time(), "%Y.%m.%d_%H.%M")
    paste0(current_time, "_", project, top_DE_name)
}

get_differentials_filename <- function(project) {
    current_time <- format(Sys.time(), "%Y.%m.%d_%H.%M")
    paste0(current_time, "_", project, "_DE_Analysis")
}

get_top_genes_data <- function(userData, top) {
    result <- data.frame(
        Cluster    = character(0),
        Gene       = character(0),
        avg_logFC  = character(0),
        p_val      = character(0),
        stringsAsFactors = F)

    if (!is.null(userData)) {
        top_genes <- userData$meta[[top]]
        if (!is.null(top_genes)) {
            #visual cleanup
            colnames(top_genes) <- tolower(colnames(top_genes))
            if ('p_val' %in% colnames(top_genes)) {
                top_genes$p_val <- signif(top_genes$p_val, 4)
            }
            if ('p_val_adj' %in% colnames(top_genes)) {
                top_genes$p_val_adj <- signif(top_genes$p_val_adj, 4)
            }
            if ('avg_logfc' %in% colnames(top_genes)) {
                top_genes$avg_logfc <- round(top_genes$avg_logfc, 2)
                top_genes <- rename(top_genes, avg_logFC = avg_logfc)
            }

            result <- top_genes %>%
                select("cluster", "gene", names(top_genes)[!(names(top_genes) %in% c('gene', 'cluster'))]) %>%
                rename(Gene = gene, Cluster = cluster) %>%
                mutate(Cluster = as.character(Cluster))
        }
    }

    result
}

get_metadata_field <- function(field, italics) {
    ifelse(is.null(field) || is.na(field) || field == "NA",
           paste0(italics[1], '&lt;missing&gt;', italics[2]),
           iconv(field[[1]], "latin1", "ASCII", sub = ""))
}

get_dataset_summary <- function(ud) {
    bold_start    <- "<b>"
    bold_end      <- "</b>"
    italic_start  <- "<i>"
    italic_end    <- "</i>"
    italics       <- c(italic_start, italic_end)
    colon         <- ":"
    break_c       <- "<br>"

    # start section
    start_section <- ""
    start_section <- paste(start_section, paste0(bold_start, "Contributors: ", bold_end, get_metadata_field(ud$meta$author, italics), break_c, break_c))
    start_section <- paste(start_section, paste0(bold_start, "Publication: ", bold_end, get_metadata_field(ud$meta$publication, italics), break_c, break_c))
    start_section <- paste(start_section, paste0(bold_start, "Summary: ", bold_end, get_metadata_field(ud$meta$summary, italics), break_c, break_c))
    start_section <- paste(start_section, paste0(bold_start, "About the data", colon, bold_end, break_c, break_c))

    # table
    table <- ""
    meta_data_table <- ud$meta$table
    if (!is.null(meta_data_table) && length(meta_data_table) > 0) {
        labels   <- names(meta_data_table)
        table   <- paste(table, paste0("<table style='width:80%'>"))
        for (i in 1:length(labels)) {
            label    <- labels[i]
            if (i %% 2 == 1) {
                table   <- paste(table, paste0("<tr><td>", bold_start, label, bold_end, "</td><td>", get_metadata_field(meta_data_table[label], italics), "</td>"))
            } else {
                table   <- paste(table, paste0("<td>", bold_start, label, bold_end, "</td><td>", get_metadata_field(meta_data_table[label], italics), "</td></tr>"))
            }
        }
        table   <- paste(table, "</table><br>")
    }
    return(HTML(paste0(start_section, table)))
}
