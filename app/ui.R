
ui <- dashboardPage(
  
  dashboardHeader(title = "CEPHALOPOD results"),
  
  dashboardSidebar(
    width = 350,
    radioButtons("data_type", "1. Select Data Type", choices = c("Diversity", "Ensemble")),
    fileInput("ncfile", "2. Choose the appropriate NetCDF File", accept = ".nc"),
    hr(),
    uiOutput("feature_select"),  # Dynamic UI for feature selection
    selectInput("variable", "4. Select Variable", choices = c("mean_values", "sd_values")),
    div(
      style = "display: flex; align-items: center; margin-bottom: 5px;",
      sliderInput("time", "5. Select Time", min = 1, max = 13, value = 1),
      div(
        textOutput("time_label"),
        style = "margin-left: 5px; margin-bottom: 10px; font-size: 13px; font-weight: bold;"
      )
    ),
    checkboxInput("contour_checkbox", "Show Contours", value = FALSE),
    checkboxInput("highlight_high_sd", "Highlight High SD Values", value = FALSE),
    uiOutput("highlight_type_ui"),
    actionButton("update", "Update Map")
  ),
  
  dashboardBody(
    tags$style(HTML("
      #map {
        width: 100% !important;
        height: 800px;
      }
      .centered-box {
        display: flex;
        justify-content: center;
        align-items: center;
      }
      .center-content {
        margin: auto;
        text-align: center;
      }
      .checkbox-inline {
        margin-bottom: 5px;
      }
      #contour_checkbox,
      #highlight_high_sd {
        margin-right: 10px;
      }
    ")),
    mainPanel(
      tabsetPanel(
        id = "tabselected",
        tabPanel(
          "Explanations",
          fluidRow(
            box(
              title = "Explanation",
              width = 14,
              "This page provides explanations and details about the NetCDF file and its variables."
            ),
            uiOutput("ensemble_box")
          )
        ),
        tabPanel(
          "Map",
          fluidRow(
            box(
              title = "Map Display",
              width = 14,
              leafletOutput("map", height = "400px")
            )
          ),
          fluidRow(
            column(
              width = 6,
              plotOutput("legend_plot", height = 100)
            ),
            column(
              width = 6,
              box(
                title = "Explanations",
                width = NULL,
                textOutput("explanation_text")
              )
            )
          )
        )
      )
    )
  ),
  skin = "black"
)
