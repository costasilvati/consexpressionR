# library(shiny)
# library(shinyFiles)
# library(fs)

runConsexpressionTestes <- function(){
  ui <- shiny::fluidPage(
    shiny::tags$style(
      ".well-panel-tools{
          min-height:600px;
      }
      .well-panel-tools2{
        min-height:630px;
      }
      .center{
        margin-left:8%;
        margin-top:2%;
      }
      "

    ),
    shiny::fluidRow(
      h1("consexpression", span("R", style = "font-weight: 300"),
         style = "color: #fff; text-align: center;
        background-color:#27296d;
        padding: 5%;
        margin-bottom: 2%;"),
    ),
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
        shiny::column(width = 1)
      ),
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
        shiny::column(width = 1),
        shiny::column(width = 10,
                      shiny::h2("Expression analysis parameters", style="text-aligner:center;"),
                      shiny::p("The consensus methodology of this tool was tested with the parameters shown below."),
                      shiny::fluidRow(
                        #----- LIMMA-------------
                        shiny::column( width = 3,
                                       shiny::wellPanel( class="well-panel-tools",
                                         shiny::h3("limma", style="color: #5e63b6;"),
                                         shiny::selectInput(inputId = "methNormLimmaInp",
                                                            label = "Method to normalize library sizes",
                                                            choices = c("TMM", "TMMwsp", "RLE", "upperquartile", "none"),
                                                            selected = "TMM",
                                                            width = NULL),
                                         shiny::selectInput(inputId = "adjPvalueLimma",
                                                            label = "Method to adjust P-Values",
                                                            choices = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none"),
                                                            selected = "BH"),
                                         shiny::numericInput(inputId = "topTableLimmaInp",
                                                             label = "Number lines in Top Table",
                                                             value = 1000000),
                                         shiny::h4("Differential Expression Metrics"),
                                         shiny::numericInput(inputId = "lfcMinLimmaInp",
                                                             label = "Log Fold Change less or equal to",
                                                             value = -2.0,
                                                             step = 2),
                                         shiny::numericInput(inputId = "lfcMaxLimmaInp",
                                                             label = "Log Fold Change greater or equal to",
                                                             value = 2.0,
                                                             step = 2),
                                         shiny::numericInput(inputId = "pValueLimmaInp",
                                                             value = 0.05,
                                                             label = "Maximum P-value",
                                                             step = 2)
                                       )#wellPanel
                        ), #column 3
                        #----- SAMSEQ-------------
                        shiny::column( width = 3,
                                       shiny::wellPanel( class="well-panel-tools",
                                         shiny::h3("SAMSeq", style="color: #5e63b6;"),
                                         shiny::selectInput(inputId = "respTypeSamseqInp",
                                                            choices = c("Quantitative", "Two class unpaired","Survival", "Multiclass", "Two class paired"),
                                                            selected = "Two class unpaired",
                                                            label = "Problem type"),
                                         shiny::numericInput(inputId = "npermSamseq",
                                                             value = 100,
                                                             label = "Number of Permutations"),
                                         shiny::h4(""),
                                         shiny::h4("Differential Expression Metrics"),
                                         shiny::numericInput(inputId = "lfcMinSamseqInp",
                                                             label = "Log Fold Change less or equal to",
                                                             value = -2.0,
                                                             step = 2),
                                         shiny::numericInput(inputId = "lfcMaxSamseqInp",
                                                             label = "Log Fold Change greater or equal to",
                                                             value = 2.0,
                                                             step = 2),
                                         shiny::numericInput(inputId = "pValueSamseqInp",
                                                             value = 0.8,
                                                             max = 1.0,
                                                             min = 0.0,
                                                             label = "Minimum q-value(%)",
                                                             step = 2),
                                       )
                        ),
                        #----- DESEQ2 -------------
                        shiny::column( width = 3,
                                       shiny::wellPanel( class="well-panel-tools",
                                         shiny::h3("DESeq2", style="color: #5e63b6;"),
                                         shiny::selectInput(inputId = "fitTypeDeseq2Inp",
                                                            label = "fitType",
                                                            choices = c("parametric"="parametric", "local"="local", "mean"="mean", "glmGamPoi"="glmGamPoi"),
                                                            selected = "local",
                                                            width = NULL
                                         ), #selectInput

                                         shiny::h4("Differential Expression Metrics", class="spacey"),
                                         shiny::numericInput(inputId = "lfcMinDeseq2Inp",
                                                             label = "Log Fold Change less or equal to",
                                                             value = -2.0,
                                                             step = 2),
                                         shiny::numericInput(inputId = "lfcMaxDeseq2Inp",
                                                             label = "Log Fold Change greater or equal to",
                                                             value = 2.0,
                                                             step = 2),
                                         shiny::numericInput(inputId = "pValueDeseq2Inp",
                                                             value = 0.05,
                                                             label = "Maximum P-value ",
                                                             max = 1.0,
                                                             min = 0.0,
                                                             step = 2),

                                         #shiny::fileInput("kallistoDirRuns", "Select a .tsv file in kallisto output folder",
                                         #                 multiple = FALSE,
                                         #                 accept = c(".tsv")),
                                         #helpText("Note: select the folder in default output by kallisto"),
                                         #shiny::verbatimTextOutput("selectedFolder")
                                         # shiny::fileInput("kallistoDirRuns", "Directory that contais kallisto runs.",)
                                         # pathDirRuns = ".",
                                         # pathReportFile = "report_txi.txt",
                                         # subDirRuns = "dir_runs",
                                         # fileKallisto = "abundance.tsv"
                                         # shiny::helpText("Use the default (selected) settings of the expression analysis methods, or configure manually")
                                       ), #wellPanel DESeq2
                        ),
                        #----- EDGER -------------
                        shiny::column(width = 3,
                                      shiny::wellPanel( class="well-panel-tools",
                                        shiny::h3("edgeR", style="color: #5e63b6;"),
                                        shiny::selectInput(inputId = "methNormEdgerInp",
                                                           label = "Method to normalize library sizes",
                                                           choices = c("TMM", "TMMwsp", "RLE", "upperquartile", "none"),
                                                           selected = "TMM",
                                                           width = NULL
                                        ),
                                        shiny::h4("Differential Expression Metrics", class="spacey"),
                                        shiny::numericInput(inputId = "lfcMinEdgerInp",
                                                            label = "Log Fold Change less or equal to",
                                                            value = -2,
                                                            step = 2),
                                        shiny::numericInput(inputId = "lfcMaxEdgerInp",
                                                            label = "Log Fold Change greater or equal to",
                                                            value = 2,
                                                            step = 2),
                                        shiny::numericInput(inputId = "pValueEdgerInp",
                                                            value = 0.05,
                                                            label = "Maximum P-value ",
                                                            max = 1.0,
                                                            min = 0.0,
                                                            step = 2),
                                      ),#wellPanel edgeR
                        ),

        ), #column(width = 10
        ),
        shiny::column(width = 1)
      ), # fluidRow
      shiny::fluidRow(
        shiny::column(width = 1),
        shiny::column(width = 10,
                shiny::fluidRow(
                  #----- NOISEQ -------------
                  shiny::column(width = 3,
                                shiny::wellPanel( class="well-panel-tools2",
                                  shiny::h3("NOISeq", style="color: #5e63b6;"),
                                  shiny::selectInput(inputId = "methNormNoiseqInp",
                                                     label = "Method to normalize library sizes",
                                                     choices = c("rpkm","uqua","tmm","n"),
                                                     selected = "rpkm",
                                                     width = NULL ),
                                  shiny::selectInput(inputId = "replicatesNoiseqInp",
                                                     label = "Type of replicates",
                                                     choices = c("technical", "biological","no"),
                                                     selected = "technical",
                                                     width = NULL
                                  ),
                                  shiny::numericInput(inputId = "lcNoiseqInp",
                                                      label = "lc: Length correction",
                                                      value = 0),
                                  shiny::textInput(inputId = "factorNoiseqInp",
                                                   label = "factor: Factor name",
                                                   value = "Tissue"),
                                  shiny::numericInput(inputId = "kNoiseqInp",
                                                      label = "k: Values with 0 are replaced by",
                                                      value = 0.5,
                                                      step= 2),
                                  shiny::h4("Differential Expression Metrics"),
                                  shiny::numericInput(inputId = "probNoiseqInp",
                                                      value = 0.8,
                                                      max = 1.0,
                                                      min = 0.0,
                                                      label = "Minimum Probability",
                                                      step = 2)
                                ),#wellPanel NOISeq
                  ),
                  #----- KNOWSEQ-------------
                  shiny::column( width = 3,
                                 shiny::wellPanel( class="well-panel-tools2",
                                   shiny::h3("KnowSeq", style="color: #5e63b6;"),
                                   shiny::selectInput(inputId = "filterKnowseqInp",
                                                      label = "Filter",
                                                      choices = c("ensembl_gene_id", "external_gene_name", "percentage_gene_gc_content","entrezgene_id"),
                                                      selected = "ensembl_gene_id",
                                                      width = NULL),
                                   shiny::radioButtons(
                                     inputId = "notHumanKnowseq",
                                     label = "Not Human",
                                     choices = c("TRUE" =TRUE, "FALSE" =FALSE),
                                     selected = FALSE),
                                   shiny::h4("Differential Expression Metrics", class="spacey-m"),
                                   shiny::numericInput(inputId = "lfcMinKnowseqInp",
                                                       label = "Log Fold Change less or equal to",
                                                       value = -2,
                                                       step = 2),
                                   shiny::numericInput(inputId = "lfcMaxKnowseqInp",
                                                       label = "Log Fold Change greater or equal to",
                                                       value = 2,
                                                       step = 2),
                                   shiny::numericInput(inputId = "pValueKnowseqInp",
                                                       value = 0.05,
                                                       label = "Maximum P-value ",
                                                       max = 1.0,
                                                       min = 0.0,
                                                       step = 2),
                                 ),

                  ),
                  #----- EBSEQ-------------
                  shiny::column( width = 3,
                                 shiny::wellPanel( class="well-panel-tools2",
                                   shiny::h3("EBSeq", style="color: #5e63b6;"),
                                   shiny::numericInput(inputId = "pValueEbseqInp",
                                                       value = 0.05,
                                                       label = "FDR",
                                                       max = 1.0,
                                                       min = 0.0,
                                                       step = 2),
                                   shiny::numericInput(inputId = "maxRoundEbseqInp",
                                                       value = 50,
                                                       label = "Number of iterations (maxround)"),
                                  shiny::selectInput("methodEbseqInp",
                                                     label = "Method to GetDEResults",
                                                     choices = c("robust","classic"),
                                                     selected = "robust"),
                                  shiny::h4("Differential Expression Metrics", class="spacey"),
                                  shiny::selectInput( "classDeEbseqInp",
                                                     label = "DE class",
                                                     choices = c("DE","EE"),
                                                     selected = "DE"),
                                 )
                  ),
                ),
                shiny::div( style= "text-align: center;",
                            shiny::actionButton(inputId = "goDeg",
                                                label = "Execute Diffrencial expression analysis",
                                                width = '51%',
                                                class= "btn-primary" ),
                )
        ),
        shiny::column(width = 1)
      ),
    shiny::fluidRow(
      shiny::column(width = 10, class="center",
                    shiny::verbatimTextOutput("execResults"))
    ),
    shiny::fluidRow(
      shiny::column(width = 10, class="center",

                    shiny::div( style= "text-align: center;",
                                shiny::actionButton(inputId = "goUpsetPlot",
                                                    label = "View consensus plot",
                                                    width = '51%',
                                                    class= "btn-primary" ),
                    ),
                    shiny::plotOutput("upsetPlot"),
                    ),
    )
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

    # tentar .zip ou root
    # output$selectedFolder <- renderText({
    #   inFile <- input$kallistoDirRuns
    #   if(is.null(inFile)){
    #     return(NULL)
    #   }else{
    #     dirname(inFile$datapath)
    #   }
    # })

    output$sample <- DT::renderDataTable({
      DT::datatable(datasetCount())
    })

    degList <- eventReactive(input$goDeg, {
      consexpressionR(numberReplics = input$numberReplicsInp,
                                   rDataFrameCount = datasetCount(),
                                   groupName = c(unlist(strsplit(input$groupNameInp, ","))),
                                   experimentName=input$experimentNameInp,
                                   #outDirPath="consexpression2_results/",
                                   #printResults=FALSE,
                                   methodNormLimma = input$methNormLimmaInp,
                                   methodAdjPvalueLimma = input$adjPvalueLimma,
                                   numberTopTableLimma = input$topTableLimmaInp,
                                   respTypeSamseq = input$respTypeSamseqInp,
                                   npermSamseq = input$npermSamseq,
                                   fitTypeDeseq2 = input$fitTypeDeseq2Inp,
                                   methodNormEdgeR = input$methNormEdgerInp,
                                   normNoiseq = input$methNormNoiseqInp,
                                   kNoiseq = input$kNoiseqInp,
                                   factorNoiseq=input$factorNoiseqInp,
                                   lcNoiseq = input$lcNoiseqInp,
                                   replicatesNoiseq = input$replicatesNoiseqInp,
                                   filterIdKnowseq=input$filterKnowseqInp,
                                   notSapiensKnowseq = as.logical(input$notHumanKnowseq),
                                   fdrEbseq=input$pValueEbseqInp,
                                   maxRoundEbseq = input$maxRoundEbseqInp,
                                   methodDeResultsEbseq = input$methodEbseqInp)
    })

  output$execResults <- renderPrint({
    degList()
  })

  }

  # Run the application
  shiny::shinyApp(ui = ui, server = server)
}
