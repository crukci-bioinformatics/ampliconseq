fluidPage(

  titlePanel(
    title = NULL,
    # title = "Amplicon Sequencing",
    windowTitle = "ampliconseq"
  ),

  tabsetPanel(
    id = "mainTabsetPanel",

    tabPanel(
      "Read counts",
      tags$div(style="line-height:50%;", br()),
      fileInput(
        'readCountFile',
        'Upload read count file',
        accept = c('text/tab-separated-values', ".txt")
      ),
      tags$p("There may be a short delay after uploading a large read count file before the table below is updated."),
      htmlOutput("readCountSummary"),
      tags$br(),
      DT::dataTableOutput("readCountTable")
    ),

    tabPanel(
      "SNVs",
      tags$div(style="line-height:50%;", br()),
      fileInput(
        'snvFile',
        'Upload SNV file',
        accept = c('text/tab-separated-values', ".txt")
      ),
      tags$p("Read counts and allele fractions are obtained/calculated from the read count file."),
      htmlOutput("snvSummary"),
      tags$br(),
      DT::dataTableOutput("snvTable")
    ),

    tabPanel(
      "Locations",
      tags$div(style="line-height:50%;", br()),

      sidebarLayout(

        sidebarPanel(
          htmlOutput("locationSummary"),
          tags$br(),
          DT::dataTableOutput("locationTable"),
          tags$br(),
          selectInput(
            "locationAmplicon",
            "Amplicon",
            choices = NULL
          ),
          radioButtons(
            "locationAlternateAllele",
            "Alternate allele",
            choices = bases,
            inline = TRUE
          ),
          numericInput(
            "locationMinimumDepth",
            "Minimum depth",
            100, min = 0, step = 500
          ),
          tags$div(style="line-height:50%;", br()),
          radioButtons(
            "locationDistribution",
            "Theoretical distribution",
            c("Normal", "Log-normal", "Beta"),
            selected = "Beta",
            inline = TRUE
          ),
          radioButtons(
            "locationThresholdProbability",
            "Threshold probability",
            thresholdProbabilities,
            selected = initialThresholdProbability,
            inline = TRUE
          ),
          numericInput(
            "locationExcludeHighestProportionForFitting",
            "Proportion with highest AF to exclude",
            0.1, min = 0.0, max = 0.1, step = 0.005
          ),
          numericInput(
            "locationMaximumAlleleFractionForFitting",
            "Maximum allele fraction for fitting",
            0.03, min = 0.0001, max = 1.0, step = 0.05
          ),
          checkboxInput("includeZeroAlleleFractionValues", "Include zero allele fraction values", value = FALSE),
          tags$div(style="line-height:50%;", br()),
          checkboxInput("setLocationMaximumAlleleFraction", "Set maximum displayed allele fraction", value = FALSE),
          conditionalPanel(
            condition = "input.setLocationMaximumAlleleFraction",
            numericInput(
              "locationMaximumAlleleFraction",
              "Maximum allele fraction",
              1.0, min = 0.0001, max = 1.0, step = 0.05
            )
          )
        ),

        mainPanel(
          htmlOutput("locationSubstitutionSummary"),
          tags$div(style="line-height:50%;", br()),
          tabsetPanel(
            tabPanel(
              "Scatter/box plot",
              tags$div(style="line-height:100%;", br()),
              highchartOutput("locationScatterBoxPlot", height = "500px"),
              tags$p("Drag within plot to zoom. Select a sample by clicking on a point."),
              actionButton("locationCheckSelectedButton", "Check selection"),
              actionButton("locationClearSelectedButton", "Clear selection"),
              tags$div(style="line-height:50%;", br()),
              DT::dataTableOutput("locationSelectedSampleTable"),
              tags$div(style="line-height:100%;", br())
            ),
            tabPanel(
              "Density plot",
              tags$div(style="line-height:100%;", br()),
              highchartOutput("locationDensityPlot", height = "600px"),
              tags$p("Click and drag within plot to zoom.")
            ),
            tabPanel(
              "Cullen & Frey graph",
              plotOutput("locationCullenFreyGraph", height = "650px")
            )
          )
        )
      )
    ),

    tabPanel(
      "Libraries",
      tags$div(style="line-height:50%;", br()),

      sidebarLayout(

        sidebarPanel(
          htmlOutput("librarySummary"),
          tags$br(),
          DT::dataTableOutput("libraryTable"),
          tags$br(),
          selectInput(
            "librarySubstitution",
            "Sustitution",
            choices = substitutions,
            selected = substitutions[1]
          ),
          numericInput(
            "libraryMinimumDepth",
            "Minimum depth",
            100, min = 0, step = 500
          ),
          tags$div(style="line-height:50%;", br()),
          radioButtons(
            "libraryDistribution",
            "Theoretical distribution",
            c("Normal", "Log-normal", "Beta"),
            selected = "Beta",
            inline = TRUE
          ),
          radioButtons(
            "libraryThresholdProbability",
            "Threshold probability",
            thresholdProbabilities,
            selected = initialThresholdProbability,
            inline = TRUE
          ),
          numericInput(
            "libraryExcludeHighestProportionForFitting",
            "Proportion with highest AF to exclude",
            0.1, min = 0.0, max = 0.1, step = 0.005
          ),
          numericInput(
            "libraryMaximumAlleleFractionForFitting",
            "Maximum allele fraction for fitting",
            0.03, min = 0.0001, max = 1.0, step = 0.05
          ),
          tags$div(style="line-height:50%;", br()),
          checkboxInput("setLibraryMaximumAlleleFraction", "Set maximum displayed allele fraction", value = FALSE),
          conditionalPanel(
            condition = "input.setLibraryMaximumAlleleFraction",
            numericInput(
              "libraryMaximumAlleleFraction",
              "Maximum allele fraction",
              1.0, min = 0.0001, max = 1.0, step = 0.05
            )
          )
        ),

        mainPanel(
          htmlOutput("librarySubstitutionSummary"),
          tags$div(style="line-height:50%;", br()),
          tabsetPanel(
            tabPanel(
              "Scatter/box plot",
              tags$div(style="line-height:100%;", br()),
              highchartOutput("libraryScatterBoxPlot", height = "500px"),
              tags$p("Drag within plot to zoom. Select a location by clicking on a point."),
              actionButton("libraryCheckSelectedButton", "Check selection"),
              actionButton("libraryClearSelectedButton", "Clear selection"),
              tags$div(style="line-height:50%;", br()),
              DT::dataTableOutput("librarySelectedLocationTable"),
              tags$div(style="line-height:100%;", br())
            ),
            tabPanel(
              "Density plot",
              tags$div(style="line-height:100%;", br()),
              highchartOutput("libraryDensityPlot", height = "600px"),
              tags$p("Click and drag within plot to zoom.")
            ),
            tabPanel(
              "Cullen & Frey graph",
              plotOutput("libraryCullenFreyGraph", height = "650px")
            )
          )
        )
      )
    )

  )
)

