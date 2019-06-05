# ----------------------------------------
# --      PROGRAM server_global.R       --
# ----------------------------------------
# USE: Server-specific variables and
#      functions for the main reactive
#      shiny server functionality.  All
#      code in this file will be put into
#      the framework outside the call to
#      shinyServer(function(input, output, session)
#      in server.R
#
# NOTEs:
#   - All variables/functions here are
#     SERVER scoped and are available
#     across all user sessions, but not to
#     the UI.
#
#   - For user session-scoped items
#     put var/fxns in server_local.R
#
# FRAMEWORK VARIABLES
#     none
# ----------------------------------------

# -- IMPORTS --
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(Seurat))

if (packageVersion("Seurat") < "3.0.0") { stop("Need to upgrade Seurat to >= version 3.0.0") }

# -- VARIABLES --


# -- FUNCTIONS --


