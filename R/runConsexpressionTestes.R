# library(shiny)
# library(shinyFiles)
# library(fs)

runConsexpressionTestes <- function(){
  ui <- shiny::fluidPage(
    h1("consexpression", span("R", style = "font-weight: 300"),
       style = "color: #fff; text-align: center;
        background-color:#27296d;
        padding: 5%;
        margin-bottom: 2%;"),
    shiny::fluidRow(
      shiny::h2("Configure Experiment", style = "text-align: center;"),
      shiny::fluidRow(
        shiny::column(width = 1,
        ),
        shiny::column(width = 10,
                      shiny::wellPanel(shiny::h3("Experiment Design", style="color: #5e63b6;"),
                                       shiny::fluidRow(
                                         shiny::column(width = 6,
                                                       shiny::textInput(inputId = "experimentNameInp",
                                                                        label="Experiment Name",
                                                                        value="genericExperiment",
                                                                        placeholder = "genericExperiment"),
                                         ),
                                         shiny::column(width = 6,
                                                       shiny::numericInput(inputId = "numberReplicsInp",
                                                                           label="Number of Replics",
                                                                           value= 1,
                                                                           min = 1,),
                                         ),
                                         shiny::p(shiny::helpText(
                                           "Note: This tool expect the same number of replics in each group of treatment."
                                         )),
                                       ),



                                       shiny::textInput(inputId ="groupNameInp",
                                                        label="Treatment Names",
                                                        placeholder ="Treat1, Treat2"),
                                       shiny::helpText("Note: Comma separeted list by sample treatment names."),
                                       shiny::h3("Table count file", style="color: #5e63b6;"),
                                       # Input: Selector for choosing dataset ----
                                       shiny::fluidRow(
                                         shiny::column(width = 6,
                                                       shiny::selectInput(inputId = "sepCharcterInp",
                                                                          label = "Choose a separator:",
                                                                          choices = c("Comma-separated"=",","TAB" = "\t")),
                                         ),

                                         shiny::column(width = 6,
                                                       shiny::fileInput(inputId = "tableCountInp",
                                                                        label = "Select a table count file",
                                                                        multiple = FALSE,
                                                                        accept = c("text/csv/tsv",
                                                                                   "text/comma-separated-values,text/plain",".csv")),
                                                       helpText("Note: while the data view will show only the specified", "number of replics, the summary will still be based","on the full dataset."),
                                         )

                                       ),
                                       shiny::div( style= "text-align: center;",
                                                   shiny::actionButton(inputId = "go",
                                                                       label = "Load count dataset",
                                                                       width = '51%',
                                                                       class= "btn-primary" ),
                                       )

                      ),
        ),
        column(width = 1)
      ),
      shiny::fluidRow(
        shiny::column(width = 12,
                      shiny::h2("Details of data set load", style = "text-align: center;"),
                      shiny::wellPanel(
                        DT::dataTableOutput("sample")
                      )
        ),
      ),
      shiny::fluidRow(
        shiny::column(width = 5,
                      shiny::h2("Expression analysis parameters", style="text-aligner:center;"),
                      shiny::p("The consensus methodology of this tool was tested with the parameters shown below."),
                      shiny::wellPanel(
                        shiny::h3("DESEq2", style="color: #5e63b6;"),
                        shiny::selectInput(inputId = "fitTypeInp",
                                           label = "fitType",
                                           choices = c("parametric"="parametric", "local"="local", "mean"="mean", "glmGamPoi"="glmGamPoi"),
                                           selected = "local",
                                           width = NULL
                        ), #selectInput
                        shiny::fileInput("kallistoDirRuns", "Select a .tsv file in kallisto output folder",
                                         multiple = FALSE,
                                         accept = c(".tsv")),
                        helpText("Note: select the folder in default output by kallisto"),
                        shiny::verbatimTextOutput("selectedFolder")
                        # shiny::fileInput("kallistoDirRuns", "Directory that contais kallisto runs.",)
                        # pathDirRuns = ".",
                        # pathReportFile = "report_txi.txt",
                        # subDirRuns = "dir_runs",
                        # fileKallisto = "abundance.tsv"
                        # shiny::helpText("Use the default (selected) settings of the expression analysis methods, or configure manually")
                      ) #wellPanel
        ) #column(width = 5,
      ), # fluidRow
    ) # fluidRow
  ) #fluidPage



  # Define server logic required to draw a histogram
  server <- function(input, output) {
    options(shiny.maxRequestSize=30*1024^2)

    datasetCount <- eventReactive(input$go, {
      inFile <- input$tableCountInp
      readCountFile(inFile$datapath, input$sepCharcterInp)
    })

    output$sample <- DT::renderDataTable({
      DT::datatable(datasetCount())
    })

# tentar .zip indicdo pelo bing
    output$selectedFolder <- renderText({
      inFile <- input$kallistoDirRuns
      if(is.null(inFile)){
        return(NULL)
      }else{
        dirname(inFile$datapath)
      }
    })
  }

  # Run the application
  shiny::shinyApp(ui = ui, server = server)
}
