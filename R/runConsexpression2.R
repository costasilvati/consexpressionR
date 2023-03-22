library(shiny)

#' Execute app shinny version
#'
#' @return void
#' @export
#'
#' @examples runConsexpression2()
runConsexpression2 <- function(){
    # Define UI for application that draws a histogram
    ui = fluidPage(
        #shinythemes::themeSelector(),  # <--- Add this somewhere in the UI
        theme = shinythemes::shinytheme("journal"),
        # App title ----
        titlePanel("Consexpression2"),
        sidebarPanel(
            # numberReplics,
            # groupName,
            # tableCountPath,
            # sepCharacter="\t",
            # experimentName="genericExperiment",
            # outDirPath="consexpression2_results/",
            # methodNorm = "TMM",
            # methodAdjPvalue = "BH",
            # numberTopTable = 1000000,
            # printResults=FALSE,
            # kallistoReport = "report.txt",
            # kallistoDir = "kallisto_quant",
            # kallistoSubDir = "expermient_kallisto",
            # kallistoOut = "abundance.tsv"
            textInput(inputId = "experimentNameInp", label="Experiment Name", value="genericExperiment", placeholder = "genericExperiment"),
            numericInput(inputId = "numberReplicsInp",label="Number of Replics", value= 1, min = 1),
            numericInput(inputId ="qtdGroupInp",label="Number of Groups",value= 2, min = 2),
            # Input: Selector for choosing dataset ----
            selectInput(inputId = "sepCharcterInp",
                        label = "Choose a separator:",
                        choices = c("TAB" = "\t", "Comma-separated"=","),),
            # Input: Select a file ----
            fileInput(inputId = "tableCountInp", label = "Select a table count file (extension .CSV)",
                      multiple = FALSE,
                      accept = c("text/csv",
                                 "text/comma-separated-values,text/plain",
                                 ".csv")),
            helpText("Note: while the data view will show only the specified",
                     "number of replics, the summary will still be based",
                     "on the full dataset."),

            actionButton("update", "Update View"),

            #actionButton("action", "Button"),
            #actionButton("action2", "Button2", class = "btn-primary")
        ),
        mainPanel(
            tabsetPanel(
                id = 'dataset',
                tabPanel("Input Data", DT::dataTableOutput("tableCount")),
                tabPanel("Tab 2")
            )
        )
    )

    # Define server logic required to draw a histogram
    server <- function(input, output) {

        #options(shiny.maxRequestSize=30*1024^2)
        output$tableCount <- DT::renderDataTable(
            consexpression2::readCountFile(input$tableCountInp$datapath,
                                           input$sepCharcterInp),
            options = list(pageLength = 50,
                           initComplete = I("function(settings, json) {alert('Done.');}"))
            )
    }

    # Run the application
    shinyApp(ui = ui, server = server)
}



