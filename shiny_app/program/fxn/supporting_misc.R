# ----------------------
# Misc Functions
# ---------------------

check_panel_plot_columns <- function(genelist, panel_plot_columns) {
    if (g_debug) message('function: ', 'check_panel_plot_columns')

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
    if (g_debug) message('function: ', 'get_layout_topology')

    layout <-  FALSE
    if (panel_plot_columns > 0) {
        if (panel_plot_columns > length(genelist)) {
            panel_plot_columns <- length(genelist)
        }
        layout <- paste0(ceiling(length(genelist)/panel_plot_columns), "X", panel_plot_columns)
    }
    layout
}

get_total_genelist <- function(all_genes, selected_genes, signature_genes, additional_genes = NULL) {
    if (g_debug) message('function: ', 'get_total_genelist')

    gene_signatures_genes <- c()
    if (!is.null(signature_genes) && signature_genes != "") {
        gene_signatures_genes <- trimws(unlist(strsplit(signature_genes, ",")))
        gene_signatures_genes <- intersect(gene_signatures_genes, all_genes$Symbol)
    }
    total_genes <- unique(c(selected_genes, gene_signatures_genes, additional_genes))
    total_genes
}

get_dynamic_plot_height <- function(gene_count) {
    if (g_debug) message('function: ', 'get_dynamic_plot_height')

    600 + 20 * gene_count
}

get_metadata_field <- function(field, italics) {
    if (g_debug) message('function: ', 'get_metadata_field')

    ifelse(is.null(field) || is.na(field) || field == "NA",
           paste0(italics[1], '&lt;missing&gt;', italics[2]),
           iconv(field[[1]], "latin1", "ASCII", sub = ""))
}

get_dataset_summary <- function(ud, report_modus = FALSE) {
    if (g_debug) message('function: ', 'get_dataset_summary')

    bold          <- c(ifelse(report_modus, " **", "<b>"),   ifelse(report_modus, "** ", "</b>"))
    italics       <- c(ifelse(report_modus, " *", "<i>"),    ifelse(report_modus, "* ", "</i>"))
    h4            <- c(ifelse(report_modus, "####", "<h4>"), ifelse(report_modus, "", "</h4>"))
    colon         <- ifelse(report_modus, "\\:", ":")
    break_c       <- ifelse(report_modus, "\n", "<br>")

    # start section
    start_section <- ifelse(!(report_modus), "",
                            paste(h4[1], paste0(bold[1], get_metadata_field(ud$meta$title, italics), bold[2]),
                                  h4[2], break_c))
    if (report_modus) {
        header <- list(start_section)
        meta_data_table <- ud$meta$table
        labels_length <- 0
        df_table <- data.frame(Field = character(labels_length), Value = character(labels_length), stringsAsFactors = FALSE)
        if (!is.null(meta_data_table) && NROW(meta_data_table) > 0) {
            labels  <- meta_data_table$field
            labels_length <- length(labels)
            for (i in 1:length(labels)) {
                label <- labels[i]
                df_table[i, "Field"] <- paste0(italics[1], label, italics[2])
                df_table[i, "Value"] <- get_metadata_field(meta_data_table$value[[i]], italics)
            }
        }

        df_table[labels_length + 1, "Field"] <- paste0(italics[1], "Cell Count", italics[2])
        df_table[labels_length + 1, "Value"] <- get_metadata_field(ud$meta$cells, italics)
        df_table[labels_length + 2, "Field"] <- paste0(italics[1], "Gene Count", italics[2])
        df_table[labels_length + 2, "Value"] <- get_metadata_field(ud$meta$genes, italics)
        df_table[labels_length + 3, "Field"] <- paste0(italics[1], "Clusters", italics[2])
        df_table[labels_length + 3, "Value"] <- get_metadata_field(paste(ud$meta$clusters, collapse = '; '), italics)
        colnames(df_table) <- c("", "")
        header <- list(start_section, df_table)
        return(header)

    } else {
        table <- ""
        meta_data_table <- ud$meta$table

        table   <- paste(table, paste0("<table style='width:80%;margin-left:auto;margin-right:auto;'>"))

        if (!is.null(meta_data_table) && length(meta_data_table) > 0) {
            labels   <- meta_data_table$field
            for (i in 1:length(labels)) {
                label <- labels[i]
                if (i %% 2 == 1) {
                    table   <- paste(table, paste0("<tr><td>", bold[1], label, bold[2], "</td><td>", get_metadata_field(meta_data_table[i, 2], italics), "</td>"))
                } else {
                    table   <- paste(table, paste0("<td>", bold[1], label, bold[2], "</td><td>", get_metadata_field(meta_data_table[i, 2], italics), "</td></tr>"))
                }
            }
            #extra row added on for spacing between statistics and metadata
            table <- paste(table, paste0("<tr><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td></tr>"))
        }

        table <- paste(table, paste0("<tr><td>", bold[1], "Cell Count:", bold[2], "</td><td>", get_metadata_field(ud$meta$cells), "</td><td>&nbsp;</td><td>&nbsp;</td></tr>"))
        table <- paste(table, paste0("<td>", bold[1], "Gene Count:", bold[2], "</td><td>", get_metadata_field(ud$meta$genes), "</td><td>&nbsp;</td><td>&nbsp;</td></tr>"))
        table <- paste(table, "</table><br>")

        # add standard information statistics
        table <- paste(table,
                       "<table style='width:80%;margin-left:auto;margin-right:auto;'>",
                           "<tr style='vertical-align:top;'><td style='padding-right:10px;'>",
                               bold[1], "Clusters: ", bold[2], "</td><td>", get_metadata_field(paste(ud$meta$clusters, collapse = '; ')),
                           "</td></tr>",
                       "</table><br>")

        end_section <- ""

        return(HTML(paste0(start_section, table, end_section)))
    }
}

create_bs_button <- function(id, label, width = "100%", disabled = FALSE) {
    bsButton(inputId  = id,
             label    = label,
             style    = "primary",
             width    = width,
             icon     = icon("arrow-right"),
             disabled = disabled)
}

decimal_places <- function(x) {
    #length zero input
    if (length(x) == 0) {
        return(numeric())
    }
    # get total number char count
    x_nchr <-  x %>% abs() %>% format(scientific = FALSE) %>% nchar() %>% as.numeric()
    # get integer part count
    x_int <-  floor(x) %>% abs() %>% nchar()
    # get final result result, 1 is for the decimal point
    x_nchr <-  x_nchr - 1 - x_int
    # return 1 at least
    x_nchr[x_nchr <= 0] <-  1
    x_nchr
}

floor_decimals <- function(x, level = 1) {
    round(x - 5 * 10^(-level - 1), level)
}

ceiling_decimals <- function(x, level = 1){
    round(x + 5 * 10^(-level - 1), level)
}
