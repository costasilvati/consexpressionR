#' Execute app shinny version
#'
#' @return void
#' @export
#'
#' @examples
runConsexpressionR <- function(){
    # Namespaces in Imports field not imported from:
  # ‘bslib’ ‘nlme’ ‘readr’ ‘shinythemes’ ‘tximportData’

    ui <- shinydashboard::dashboardPage(
      shinydashboard::dashboardHeader(title = "consexpression2"),
      shinydashboard::dashboardSidebar(
        shinydashboard::sidebarMenu(
          shinydashboard::menuItem("Load Dataset", tabName = "load", icon = icon("upload")),
          shinydashboard::menuItem("Parameters", tabName = "parameters", icon = icon("gears")),
          shinydashboard::menuItem("Gene Expression", tabName = "geneExpression", icon = icon("magnifying-glass-chart"))
          )
        ),
      shinydashboard::dashboardBody(
        shinydashboard::tabItems(
          shinydashboard::tabItem(tabName="load",
                 shinydashboard::fluidRow(
                   shinydashboard::box(
                      title = "Experiment Design", status = "warning", solidHeader = TRUE,
                      # outDirPath="consexpression2_results/",
                      # methodNorm = "TMM",
                      # methodAdjPvalue = "BH",
                      # numberTopTable = 1000000,
                      # printResults=FALSE,
                      # kallistoReport = "report.txt",
                      # kallistoDir = "kallisto_quant",
                      # kallistoSubDir = "expermient_kallisto",
                      # kallistoOut = "abundance.tsv"
                      shiny::textInput(inputId = "experimentNameInp", label="Experiment Name", value="genericExperiment", placeholder = "genericExperiment"),
                      shiny::numericInput(inputId = "numberReplicsInp",label="Number of Replics", value= 1, min = 1),
                      shiny::helpText("Note: This tool expect the same number of replics in each group of treatment."),
                      shiny::numericInput(inputId ="qtdGroupInp",label="Number of Treatment Groups",value= 2, min = 2),
                      shiny::textInput("groupName", label="Treatment Names", placeholder ="Treat1, Treat2, Treat3"),
                      shiny::helpText("Note: Comma separeted list by sample treatment names."),
                      tags$p("Group and Sample configiguration in file need be:"),
                      shiny::textOutput("designExperiment"),
                      # Input: Selector for choosing dataset ----
                      ),
                   shinydashboard::box(
                        title = "File Reader Configuration",
                        status = "primary",
                        solidHeader = TRUE,
                        shiny::selectInput(inputId = "sepCharcterInp",
                                    label = "Choose a separator:",
                                    choices = c("TAB" = "\t", "Comma-separated"=","),),
                        # Input: Select a file ----
                        shiny::fileInput(inputId = "tableCountInp",
                                  label = "Select a table count file (extension .CSV)",
                                  multiple = FALSE,
                                  accept = c("text/csv",
                                             "text/comma-separated-values,text/plain",
                                             ".csv")),
                        # helpText("Note: while the data view will show only the specified", "number of replics, the summary will still be based","on the full dataset."),

                        shinyjs::disabled(actionButton("confirmModel", "Go to expression analysis", class="btn-primary", hidden=TRUE))

                        #actionButton("action", "Button"),
                        #actionButton("action2", "Button2", class = "btn-primary")
                        )
                    ),
                 shinydashboard::fluidRow(
                      title = "Read Data", status = "primary", solidHeader = TRUE,
                      DT::dataTableOutput("sample")

                    ),
                #useShinyjs()
            ),
          shinydashboard::tabItem(tabName = "parameters",
                    h2("Set parameters of gene expression tools")
            ),
          shinydashboard::tabItem(tabName = "geneExpression",
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
        observe({
          shinyjs::enable("confirmModel", condition)
        })
        return(df)
      })

      designExperiment <- shiny::reactive({
        shiny::validate(
          shiny::need(input$groupName != "", "Please type a list of treatment names separeted by comma.")
        )
        groupName <- c(unlist(strsplit(input$groupName, ",")))
        nr <- input$numberReplicsInp
        d <- rep(groupName, each=nr)
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
              df <- df_upload()
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




