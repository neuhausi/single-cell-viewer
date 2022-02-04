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
set_app_parameters(title = "Single Cell Viewer",
                   titleinfo   = HTML(sprintf(readLines("program/data/about.html"), "3.2.0")),
                   loglevel    = "DEBUG",
                   showlog     = FALSE,
                   app_version = "3.2.0")

# -- PROGRAM --
suppressPackageStartupMessages({
    library(canvasXpress)
    library(glue)
    library(periscope)
    library(DT)
})

source('program/fxn/supporting_data.R')
source('program/fxn/supporting_plots.R')
source('program/fxn/supporting_misc.R')
source('program/fxn/diff_expression.R')

# Variables
g_debug         <- FALSE
g_api_con       <- NULL  #setup on first connection
g_api_pipelines <- NULL  #setup on first connection

# RSConnect deployment
g_clusters_path <- "program/data/allowed_clusters.yaml"

g_default_cluster              <- "CellType.select"
g_gene_signatures              <- get_gene_signatures()
g_differential_logfc_threshold <- 0.5
g_differential_gene_threshold  <- 1000
g_differential_pct_threshold   <- 1
g_differential_min_no_cells    <- 3
g_de_message_low               <- 500
g_de_message_high              <- 5000
