# ----------------------------------------
# --       PROGRAM ui_sidebar.R         --
# ----------------------------------------
# USE: Create UI elements for the
#      application sidebar (left side on
#      the desktop; contains options) and
#      ATTACH them to the UI by calling
#      add_ui_sidebar_basic() or
#      add_ui_sidebar_advanced()
#
# NOTEs:
#   - All variables/functions here are
#     not available to the UI or Server
#     scopes - this is isolated
# ----------------------------------------

# -- IMPORTS --
suppressPackageStartupMessages(library(shinyjs))


# ----------------------------------------
# --     SIDEBAR ELEMENT CREATION       --
# ----------------------------------------

# -- Create Basic Elements

sidebar_header <- tags$div(align = "center",
                           tags$h3('Chart Options'), br(),
                           useShinyjs())

clusters_select <- selectizeInput("clustersSel",
                               label    = ui_tooltip("clustersTooltip", "Cluster",
                                                     "Choose from clusters available in the dataset."),
                               choices  = NULL,
                               multiple = FALSE,
                               options  = list(placeholder = "Type/Click then Select",
                                               searchField = "value",
                                               plugins     = list('remove_button')))

genes_select <- selectizeInput("genesSel",
                   label    = ui_tooltip('genesTooltip', "Genes",
                                         'Choose from genes expressed in the dataset. <i>Unexpressed genes are not available to be chosen.</i>'),
                   choices  = NULL,
                   multiple = TRUE,
                   options  = list(placeholder = "Type/Click then Select",
                                   searchField = "value",
                                   plugins     = list('remove_button')))

gene_signatures_select <- selectizeInput("geneSignaturesSel",
                                         label    = ui_tooltip('geneSignaturesTooltip', "Gene Signatures",
                                                               'Choose from custom gene signatures. <i>Only genes expressed in the dataset will be available for plotting and analysis.</i>'),
                                         choices  = NULL,
                                         multiple = TRUE,
                                         options  = list(placeholder = "Type/Click then Select",
                                                         searchField = "value",
                                                         plugins     = list('remove_button')))

selected_gene_signatures_text <- uiOutput("selected_gene_signatures_text")


# -- Register Basic Elements in the ORDER SHOWN in the UI
add_ui_sidebar_basic(list(sidebar_header,
                          clusters_select,
                          genes_select,
                          gene_signatures_select,
                          selected_gene_signatures_text))

# -- Create Advanced Elements

# add_ui_sidebar_advanced(list())

