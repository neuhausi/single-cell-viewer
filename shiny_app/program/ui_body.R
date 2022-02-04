# ----------------------------------------
# --          PROGRAM ui_body.R         --
# ----------------------------------------
# USE: Create UI elements for the
#      application body (right side on the
#      desktop; contains output) and
#      ATTACH them to the UI by calling
#      add_ui_body()
#
# NOTEs:
#   - All variables/functions here are
#     not available to the UI or Server
#     scopes - this is isolated
# ----------------------------------------

# -- IMPORTS --
library(shinyFeedback)

# ----------------------------------------
# --      BODY ELEMENT CREATION         --
# ----------------------------------------

# -- Create Elements

# Add css class to radiobuttons div to position them on same row as the text.
radiobuttons_js <- '$(document).on("shiny:connected", function(e) {
                        var elements = document.getElementsByClassName("shiny-options-group");
                        for (var i = 0; i < elements.length; i++)
                        {
                            elements[i].classList.add("radio-inline");
                        }
                    });'

loading_modal <- bsModal("loading_modal",
                         title = "Single Cell Viewer Loading...",
                         trigger = NULL,
                         size = "small",
                         tags$head(tags$script("$(document).ready(function(){$('#loading_modal').modal();});"),
                                   tags$style("#loading_modal .modal-title {color:#3380A9; font-weight:bold;}"),
                                   tags$style("#loading_modal .modal-footer {display:none;}")),
                         h4(align = 'center', 'The requested data is being downloaded and will load momentarily'),
                         br(),
                         p(align = 'center', em('Click anywhere outside this dialog to continue')))
dac_modal <- bsModal("dac_modal",
                     title = "Differential Analysis Calculation",
                     trigger = NULL,
                     size = "small",
                     tags$head(tags$style("#dac_modal .modal-title {color:#3380A9; font-weight:bold;}"),
                               tags$style("#dac_modal .modal-footer {display:none;}")),
                     uiOutput("dacText"),
                     br(),
                     tags$div(actionButton("cancel", "Cancel"),
                              actionButton("proceed", "Proceed", style = "float:right")),
                     br())
cluster_change_modal <- bsModal("cluster_change_modal",
                                title = "Changing Clusters",
                                trigger = NULL,
                                size = "small",
                                tags$head(tags$style("#cluster_change_modal .modal-title {color:#3380A9; font-weight:bold;}"),
                                          tags$style("#cluster_change_modal .modal-footer {display:none;}")),
                                p("Changing clusters will clear all existing charts and analysis results..."),
                                br(),
                                tags$div(actionButton("cluster_cancel", "Cancel"),
                                         actionButton("cluster_proceed", "Proceed", style = "float:right")),
                                br())

body1 <- box(id       = "bodyElement1",
             title    = textOutput("summaryTitle"),
             width    = 12,
             collapsible = TRUE,
             collapsed   = TRUE,
             uiOutput("datasetSummary")
            )


build_gene_select <- function(id, label) {
    selectizeInput(id,
                   label    = ui_tooltip(glue("{id}_Tooltip"), label,
                                         'Choose from genes expressed in the dataset. <i>Unexpressed genes are not available to be chosen.</i>'),
                   choices  = NULL,
                   multiple = FALSE,
                   options  = list(placeholder = "Type/Click then Select",
                                   searchField = "value",
                                   plugins     = list('remove_button')))
}

gene1_select <- build_gene_select("gene1Sel", "Gene 1")
gene2_select <- build_gene_select("gene2Sel", "Gene 2")

body2 <- box(id          = "bodyElement2",
             title       = "Gene Expression Plots",
             width       = 12,
             collapsible = FALSE,
             tags$head(tags$script(HTML(radiobuttons_js))),
             useShinyFeedback(),
             tabBox(id     = "tabBodyElement2",
                    width  = 12,
                    selected = "Overview",
                    tabPanel("Overview",
                             canvasXpressOutput("cxOverviewPlot", width = "100%", height = "800px"),),
                    tabPanel("Scatterplot",
                             fluidRow(
                                 column(width = 3, selectizeInput("scatterClusterSel",
                                                                  label    = "Cluster Plotted",
                                                                  choices  = NULL,
                                                                  multiple = FALSE)),
                                 column(width = 3, textInput("scatterPanelPlotColumns", "Panel Plot Columns", "auto")),
                                 column(width = 3, tags$label(HTML("&zwnj;")), tags$br(),
                                                    create_bs_button(id    = "scatterPlotBtn",
                                                                     label = " Plot/Refresh"))),
                             tags$hr(),
                             canvasXpressOutput("cxScatterPlot", width = "100%", height = "800px")),
                    tabPanel("Violin",
                             fluidRow(
                                 column(width = 3, textInput("violinPanelPlotColumns", "Panel Plot Columns", "auto")),
                                 column(width = 3, tags$label(HTML("&zwnj;")), tags$br(),
                                        create_bs_button(id    = "violinPlotBtn",
                                                         label = " Plot/Refresh"))),
                             tags$hr(),
                             canvasXpressOutput("cxViolinPlot", width = "100%", height = "800px")),
                    tabPanel("Heatmap",
                             fluidRow(
                                 column(width = 3, create_bs_button(id    = "heatmapPlotBtn",
                                                                    label = " Plot/Refresh"))),
                             tags$hr(),
                             uiOutput("cxHeatmapPlot")),
                    tabPanel("Co-Expression Analysis",
                             fluidRow(
                                 column(width = 3, selectizeInput("coExpClusterSel",
                                                                  label    = "Cluster Plotted",
                                                                  choices  = NULL,
                                                                  multiple = FALSE)),
                                 column(width = 3, gene1_select),
                                 column(width = 3, gene2_select),
                                 column(width = 3, tags$label(HTML("&zwnj;")), tags$br(),
                                        create_bs_button(id    = "coExpPlotBtn",
                                                         label = " Plot/Refresh"))),
                             tags$hr(),
                             fluidRow(column(width = 1),
                                      column(width = 8,
                                             tags$div(align = "right",
                                                      style = "padding-right:5%",
                                                canvasXpressOutput("cxCoExpPlot", width = "100%", height = "600px"))),
                                      column(width = 2,
                                             uiOutput("cxNoiseControls")),
                                      column(width = 1)),
                             tags$hr(),
                             tags$h4(align = "center",
                                     uiOutput("coExpTableTitle")),
                             tags$div(style = "padding-left:15%; padding-right:15%;",
                                      downloadableTableUI("coExpTable", list("csv", "tsv"),
                                                          "Download population combinations data"))),
                    tabPanel("Differential Analysis",
                             dac_modal,
                             fluidRow(
                                 column(width = 3, selectizeInput("differentialsCluster1Sel",
                                                                  label    = "Cluster 1",
                                                                  choices  = NULL,
                                                                  multiple = FALSE)),
                                 column(width = 3, selectizeInput("differentialsCluster2Sel",
                                                                  label    = "Cluster 2",
                                                                  choices  = NULL,
                                                                  multiple = FALSE)),
                                 column(width = 6,
                                        tags$label(HTML("&zwnj;")), tags$br(),
                                        create_bs_button(id    = "diffCalculateBtn",
                                                         label = "Calculate Differentials",
                                                         width = "60%",
                                                         disabled = TRUE))),
                             tags$hr(),
                             fluidRow(
                                 column(width = 6, div(canvasXpressOutput("cxDifferentialsScatterPlot1", width = "100%", height = "400px"),
                                                       style = "padding-left:5%;padding-right:2.5%;")),
                                 column(width = 6, div(canvasXpressOutput("cxDifferentialsScatterPlot2", width = "100%", height = "400px"),
                                                       style = "padding-left:2.5%;padding-right:5%;"))),
                             tags$hr(),
                             tags$h4(align = "center", uiOutput("differentialsTableTitle")),
                             uiOutput("differentialsTableAlternativeText"),
                             div(downloadableTableUI(id            = "differentialsTable",
                                                     downloadtypes = list("csv", "tsv"),
                                                     hovertext     = "Download table data",
                                                     contentHeight = "200px"),
                                 style = "padding-left:10%;padding-right:10%;"))
                )
             )

# -- Register Elements in the ORDER SHOWN in the UI
add_ui_body(list(loading_modal, cluster_change_modal, body1, body2))
