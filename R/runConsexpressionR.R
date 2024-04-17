#' Execute app shinny version
#'
#' @return void
#' @export
#'
#' @examples
runConsexpressionR <- function(){
  # Namespaces in Imports field not imported from:
  # ‘bslib’ ‘nlme’ ‘readr’ ‘shinythemes’ ‘tximportData’

  ui <- shinydashboard::dashboardPage( skin = "black",
    shinydashboard::dashboardHeader(title = "consexpressionR"),
    shinydashboard::dashboardSidebar(
      shinydashboard::sidebarMenu(
        shinydashboard::menuItem("Configure and View", tabName = "configAndView", icon = icon("upload")),
        shinydashboard::menuItem("Details of Dataset", tabName = "datasetView", icon = icon("gears")),
        shinydashboard::menuItem("Results By Methods", tabName = "resultsByMethods", icon = icon("magnifying-glass-chart"))
      )
    ),
    shinydashboard::dashboardBody(
      shinydashboard::tabItems(
        shinydashboard::tabItem(tabName="configAndView",
                                shiny::fluidRow(
                                  shinydashboard::box(
                                    title = "Experiment Design",
                                    status = "warning",
                                    solidHeader = TRUE,
                                    width = 12,
                                    style = "display: flex; flex-direction: row; flex-wrap: wrap;",
                                    # outDirPath="consexpression2_results/",
                                    # methodNorm = "TMM",
                                    # methodAdjPvalue = "BH",
                                    # numberTopTable = 1000000,
                                    # printResults=FALSE,
                                    # kallistoReport = "report.txt",
                                    # kallistoDir = "kallisto_quant",
                                    # kallistoSubDir = "expermient_kallisto",
                                    # kallistoOut = "abundance.tsv"
                                    shiny::column(
                                      width = 12,
                                      shiny::textInput(inputId = "experimentNameInp",
                                                       label="Experiment Name",
                                                       value="genericExperiment",
                                                       placeholder = "genericExperiment")
                                    ),

                                    shiny::column(
                                      width = 6,
                                      shiny::numericInput(inputId = "numberReplicsInp",
                                                        label="Number of Replics",
                                                        value= 1,
                                                        min = 1,),
                                      shiny::p(
                                        shiny::helpText("Note: This tool expect the same number of replics in each group of treatment.")
                                      )
                                      ),
                                    shiny::column(
                                      width = 6,
                                      shiny::numericInput(inputId ="qtdGroupInp",
                                                          label="Number of Treatment Groups",
                                                          value= 2,
                                                          min = 2)
                                    ),
                                    shiny::column(
                                      width = 12,
                                      shiny::textInput(inputId ="groupNameInp",
                                                       label="Treatment Names",
                                                       placeholder ="Treat1, Treat2"),
                                      shiny::helpText("Note: Comma separeted list by sample treatment names.")
                                    ),
                                    shiny::column(
                                      width = 6,
                                      # shiny::tags$p("Group and Sample configiguration in file need be:"),
                                      # Input: Selector for choosing dataset ----
                                      shiny::selectInput(inputId = "sepCharcterInp",
                                                         label = "Choose a separator:",
                                                         choices = c("Comma-separated"=",","TAB" = "\t")
                                                         )
                                    ),
                                    shiny::column(
                                      width = 6,
                                      # Input: Select a file ----
                                      shiny::fileInput(inputId = "tableCountInp",
                                                       label = "Select a table count file (extension .CSV)",
                                                       multiple = FALSE,
                                                       accept = c("text/csv",
                                                                  "text/comma-separated-values,text/plain",
                                                                  ".csv")),
                                      helpText("Note: while the data view will show only the specified", "number of replics, the summary will still be based","on the full dataset."),
                                    )
                                  ),
                                ),
                          shiny::fluidRow(
                                shinydashboard::box(
                                  title = "Intersection to DE indications",
                                  solidHeader = TRUE,
                                  width = 12,
                                  shiny::textOutput("designExperiment"),
                                  shiny::plotOutput(outputId = "upset")
                                ),
                          ),

        ),
        shinydashboard::tabItem(tabName = "datasetView",

                                h2("details of the data set read"),
                                shiny::fluidRow(
                                  title = "Read Data", status = "primary", solidHeader = TRUE,
                                  DT::dataTableOutput("sample")
                                ),
        ),
        shinydashboard::tabItem(tabName = "resultsByMethods",
                                h2("Expression Analysis"))
      )

    )
  )

  # Define server logic required to draw a histogram
  server <- function(input, output) {
    options(shiny.maxRequestSize=30*1024^2)

    df_upload <- shiny::reactive({
      inFile <- input$tableCountInp
      if (is.null(inFile))
        return(NULL)
      df <- readCountFile(inFile$datapath, input$sepCharcterInp)
      sink("log.txt")
      cons_result <- consexpressionR(numberReplics = input$numberReplicsInp,
                                      rDataFrameCount = df,
                                      groupName = c(unlist(strsplit(input$groupNameInp, ",")))
      )
      sink()
      expDef_result <- expressionDefinition(resultTool = cons_result)
      deByTool <- listDeByTool(cons_result, expDef_result)
      plotUp <- UpSetR::upset(deByTool,
                    sets = colnames(deByTool),
                    sets.bar.color = "#56B4E9",
                    order.by = "freq",
                    empty.intersections = "on")
      output$upset <- shiny::renderPlot(plotUp, res = 130)
      return(df)
    })

    designExperiment <- shiny::reactive({
      shiny::validate(
        shiny::need(input$groupNameInp != "", "Please type a list of treatment names separeted by comma.")
      )
      gpName <- c(unlist(strsplit(input$groupNameInp, ",")))
      nr <- input$numberReplicsInp
      d <- rep(gpName, each=nr)
      return(d)
    })

    output$designExperiment <- shiny::renderPrint({
      ds <- designExperiment()
      df <- df_upload()
      if(!is.null(df)){
        if(length(ds) != dim(df)[2]){
          ds <- "Error! File need has number of columns equal replics * groups length."
        }else{
          output$sample<- DT::renderDataTable({
            DT::datatable(df)
          })
        }
      }
      print(ds)

    })
  }

  # Run the application
  shiny::shinyApp(ui = ui, server = server)
}
