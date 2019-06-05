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



# ----------------------------------------
# --      BODY ELEMENT CREATION         --
# ----------------------------------------

# -- Create Elements

# Add css class to radiobuttons div to position them on same row as the text & handle behaviour related to the file upload.
application_js <- '$(document).on("shiny:connected", function(e) {
                    var elements = document.getElementsByClassName("shiny-options-group");
                    for (var i = 0; i < elements.length; i++)
                    {
                        elements[i].classList.add("radio-inline");
                    }
                   });
                   $(document).on("change", "#fileInputDialog", function(e) {
                        Shiny.onInputChange("fileChosen", this.files[0].size);
                   });
                   $(document).on("focusout", "#loading_modal", function(e) {
                        Shiny.onInputChange("fileModalClosed", Math.random());
                   });
                   Shiny.addCustomMessageHandler("openTitleInfoBox",
                        function(message) {
                           $("#titleinfobox_trigger").trigger("click");
                        }
                    );'

loading_modal <- bsModal("loading_modal",
                           title   = "Single Cell Viewer Loading...",
                           trigger = NULL,
                           size    = "small",
                           tags$head(tags$style("#loading_modal .modal-title {color:#3380A9; font-weight:bold;}"),
                                     tags$style("#loading_modal .modal-footer {display:none;}")),
                           h4(align = 'center', 'The file you have selected is larger than 10 MB. It is now being uploaded and will load momentarily.'),
                           br(),
                           p(align = 'center', em('Click anywhere outside this dialog to continue')))

file_error_modal <- bsModal("file_error_modal",
                            title   = "Single Cell Viewer File Upload Error",
                            trigger = NULL,
                            size    = "large",
                            tags$head(tags$style("#file_error_modal .modal-title {color:#3380A9; font-weight:bold; text-align:center;}"),
                                      tags$style("#file_error_modal .modal-footer {display:none;}"),
                                      tags$style("#missing_list {width:60%; margin-left:25%}"),
                                      tags$style("#missing_list li {text-align:left;}")),
                            h4(align = 'center', htmlOutput("loading_error_message")),
                            br(),
                            p(align = 'center', em('Click anywhere outside this dialog to continue')))

create_bs_button <- function(id, label, width = "100%", disabled = FALSE) {
    bsButton(inputId  = id,
             label    = label,
             type     = "primary",
             value    = FALSE,
             style    = "primary",
             size     = "default",
             width    = width,
             icon     = icon("arrow-right"),
             block    = FALSE,
             disabled = disabled)
}

body1 <- box(id       = "bodyElement1",
             title    = textOutput("summaryTitle"),
             width    = 12,
             collapsible = TRUE,
             collapsed   = TRUE,
             uiOutput("datasetSummary")
            )

body2 <- box(id          = "bodyElement2",
             title       = "Gene Expression Plots",
             width       = 12,
             collapsible = FALSE,
             tags$head(tags$script(HTML(application_js))),
             tags$link(rel = "stylesheet", type = "text/css", href = "style.css"),
             tabBox(id     = "tabBodyElement2",
                    width  = 12,
                    selected = "Overview",
                    tabPanel("Overview",
                             canvasXpressOutput("cxOverviewPlot", width = "100%", height = "600px"),
                             tags$br(),
                             tags$h4(align = "center", tags$strong('Top 10 Differentially Expressed Genes By Cluster')),
                             downloadFileButton("top30DE", c("csv", "tsv"), "Download top 30"),
                             downloadableTableUI("top10DE", list("csv", "tsv"),
                                                 "Download top 10",
                                                 contentHeight = "200px")),
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
                             canvasXpressOutput("cxScatterPlot", width = "100%", height = "600px")),
                    tabPanel("Violin",
                             fluidRow(
                                 column(width = 3, textInput("violinPanelPlotColumns", "Panel Plot Columns", "auto")),
                                 column(width = 3, tags$label(HTML("&zwnj;")), tags$br(),
                                        create_bs_button(id    = "violinPlotBtn",
                                                         label = " Plot/Refresh"))),
                             tags$hr(),
                             canvasXpressOutput("cxViolinPlot", width = "100%", height = "600px")),
                    tabPanel("DotPlot",
                             fluidRow(
                                 column(width = 6, radioButtons("addDotGenes",
                                                                label    = "Add the top differentially expressed genes:",
                                                                choices  = g_add_top_genes_options,
                                                                selected = "off",
                                                                inline   = TRUE)),
                                 column(width = 3,
                                        create_bs_button(id    = "dotPlotBtn",
                                                         label = " Plot/Refresh"))),
                             tags$hr(),
                             uiOutput("cxDotPlot")),
                    tabPanel("Heatmap",
                             fluidRow(
                                 column(width = 6, radioButtons("addHeatmapGenes",
                                                                label    = "Add the top differentially expressed genes:",
                                                                choices  = g_add_top_genes_options,
                                                                selected = "off",
                                                                inline   = TRUE)),
                                 column(width = 3, create_bs_button(id    = "heatmapPlotBtn",
                                                                    label = " Plot/Refresh"))),
                             tags$hr(),
                             uiOutput("cxHeatmapPlot")),
                    tabPanel("Differential Analysis",
                             uiOutput("differentialsText"),
                             tags$p(),
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
                                                       style = "padding-left:2.5%;padding-right:5%;")))
                             ,
                             tags$hr(),
                             tags$h4(align = "center", uiOutput("differentialsTableTitle")),
                             uiOutput("differentialsTableAlternativeText"),
                             div(heatmap_downloadableTableUI("differentialsTable", list("csv", "tsv"),
                                                             "Download table data",
                                                             contentHeight = "200px"),
                                 style = "padding-left:10%;padding-right:10%;"))
                )
             )

# -- Register Elements in the ORDER SHOWN in the UI
add_ui_body(list(loading_modal, file_error_modal, body1, body2))
