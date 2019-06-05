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

sidebar_header <- tags$div(useShinyjs())

file_input <- tags$div(
    style = "margin-top:20px;margin-bottom:25px;",
    align = "center",
    fileInput(inputId = "fileInputDialog",
              label   = NULL,
              buttonLabel = "Load New Data Object",
              width   = "80%"))

genes_select <- tags$div(
    h4("Chart Options"),
    selectizeInput("genesSel",
                   label    = ui_tooltip('genesTooltip', "Genes",
                                         'Choose from genes expressed in the dataset. <i>Unexpressed genes are not available to be chosen.</i>'),
                   choices  = NULL,
                   multiple = TRUE,
                   options  = list(placeholder = "Type/Click then Select",
                                   searchField = "value",
                                   plugins     = list('remove_button'))))

help_text <- tags$div(
    tags$br(),
    tags$h4(tags$a(href = 'canvasxpress.org', "CanvasXpress"), "Tips"),
    tags$p(style = "margin:10px;",
           tags$strong('Toolbar'), "- the main toolbar is available if you hover over the top title of the chart. Additional functionality can also be accessed by right-clicking anywhere on the chart",
           tags$br(), tags$br(),
           tags$strong('Zoom'), "- select an area or use the mouse to zoom by scrolling.  Reset the canvas by hitting Esc",
           tags$br(), tags$br(),
           tags$strong('Focus'), "- select a legend item to toggle a fade on that item",
           tags$br(), tags$br(),
           tags$strong('Select'), "- select points on a plot using shift-drag to select an area.  Deselect by clicking a blank area or hitting Esc",
           tags$br(), tags$br(),
           tags$strong('Download'), "- hover over the main title and select the camera icon from the toolbar"))

about_text <- tags$div(
                tags$br(),
                tags$h4(style = "margin:-10px;",
                        actionLink("about_link", "About This App")))


add_ui_sidebar_basic(list(sidebar_header,
                          file_input,
                          genes_select,
                          help_text,
                          about_text),
                     tabname = "Application")

# -- Create Advanced Elements

add_ui_sidebar_advanced(list(uiOutput("filterOptions"),
                             hr()),
                        tabname = "Filtering")

