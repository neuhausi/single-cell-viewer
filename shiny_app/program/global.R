# ----------------------------------------
# --          PROGRAM global.R          --
# ----------------------------------------
# USE: Global variables and functions
#
# NOTEs:
#   - All variables/functions here are
#     globally scoped and will be available
#     to server, UI and session scopes
# ----------------------------------------


# -- Setup your Application --
g_app_version <- "0.1.0"

set_app_parameters(title = "Single Cell Viewer",
                   titleinfo = HTML(sprintf(readLines("program/data/about.html"), g_app_version)),
                   loglevel = "DEBUG",
                   showlog = FALSE,
                   app_version = g_app_version)

# -- PROGRAM --
suppressPackageStartupMessages(library(canvasXpress))

source('program/fxn/module_heatmapDownloadableTable.R')
source('program/fxn/supporting_data.R')
source('program/fxn/supporting_plots.R')
source('program/fxn/supporting_misc.R')
source('program/fxn/diff_expression.R')

# Variables

g_differential_logfc_threshold <- 0.5
g_differential_gene_threshold  <- 1000
g_differential_pct_threshold   <- 1
g_differential_min_no_cells    <- 3

# ui
g_add_top_genes_options  <- c('Off' = 'off', 'Top 10' = 'top10', 'Top 30' = 'top30')

# File upload
options(shiny.maxRequestSize = 500*1024^2) #max 500 MB
g_error_field    <- "error"
g_missing_fields <- "missing_fields"
