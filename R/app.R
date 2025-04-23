# # app.R
#
# # Carregar pacotes necessários
# library(shiny)
# library(devtools)  # ou library(remotes)
#
# # Instalar pacote do GitHub
# devtools::install_github("costasilvati/consexpressionR")
#
# # Carregar o pacote
# library(consexpressionR)

#' @export
consexpressionR <- function(){
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
      shiny::h1("consexpressionR", shiny::span("R", style = "font-weight: 200"),
                style = "color: #fff; text-align: center;
        background-color:#27296d;
        padding: 2%;
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
                                                                           min = 1),
                                         ),
                                         shiny::p(shiny::helpText(
                                           "Note: This tool expect the same number of replics in each group of treatment."
                                         )),
                                       ),
                                      shiny::fluidRow(
                                         shiny::column(width = 6,
                                                       shiny::textInput(inputId ="groupNameInp",
                                                                        label="Treatment Names",
                                                                        placeholder ="Treat1, Treat2",
                                                                        value="Control,Treat"),
                                                       shiny::helpText("Note: Comma separeted list by sample treatment names. First group name was consider treatment (reference) group. Order of groups need be the same in column file"),
                                         ),
                                         shiny::column(width = 6,
                                                       shiny::radioButtons(
                                                         inputId = "deNovoAanalysisInp",
                                                         label = "De novo assembly RNA-seq data?",
                                                         choices = c("Yes" =TRUE, "No" =FALSE),
                                                         selected = FALSE,
                                                         inline = TRUE),
                                                       shiny::conditionalPanel(
                                                         condition = "input.groupNameInp.split(',').filter(function(e){ return e.trim().length > 0 }).length > 2",
                                                         shiny::radioButtons(inputId = "condition1Inp", "Select first condition to be compared by the differential expression algorithm:",
                                                                             choices = c(""),
                                                                             inline = TRUE),
                                                         shiny::radioButtons(inputId = "condition2Inp", "Select second condition to be compared by the differential expression algorithm:",
                                                                             choices = c(""),
                                                                             inline = TRUE)
                                                       ),
                                         )
                                      ),

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
                                                       shiny::helpText("Note: while the data view will show only the specified", "number of replics, the summary will still be based","on the full dataset."),
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
                                                       shiny::conditionalPanel(
                                                         condition = "input.groupNameInp.split(',').filter(function(e){ return e.trim().length > 0 }).length > 2",
                                                         shiny::numericInput(inputId = "FMinLimma",
                                                                             label = "Minimum F-statistic value",
                                                                             value = 1.0,
                                                                             step = 2,
                                                                             max = 1.0,
                                                                             min = 0.0),
                                                       ),
                                                       shiny::conditionalPanel(
                                                         condition = "input.groupNameInp.split(',').filter(function(e){ return e.trim().length > 0 }).length <= 2",
                                                         shiny::numericInput(inputId = "lfcMinLimmaInp",
                                                                             label = "Log Fold Change less or equal to",
                                                                             value = -2.0,
                                                                             step = 2),
                                                         shiny::numericInput(inputId = "lfcMaxLimmaInp",
                                                                             label = "Log Fold Change greater or equal to",
                                                                             value = 2.0,
                                                                             step = 2),
                                                       ),
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
                                                       shiny::conditionalPanel(
                                                         condition = "input.groupNameInp.split(',').filter(function(e){ return e.trim().length > 0 }).length > 2",
                                                         shiny::numericInput(inputId = "scoreDSamseqInp",
                                                                             label = "Score(d)",
                                                                             value = 0.8,
                                                                             step = 2),
                                                       ),
                                                       shiny::conditionalPanel(
                                                         condition = "input.groupNameInp.split(',').filter(function(e){ return e.trim().length > 0 }).length <= 2",
                                                         shiny::numericInput(inputId = "lfcMinLimmaInp",
                                                                             label = "Log Fold Change less or equal to",
                                                                             value = -2.0,
                                                                             step = 2),
                                                         shiny::numericInput(inputId = "lfcMaxLimmaInp",
                                                                             label = "Log Fold Change greater or equal to",
                                                                             value = 2.0,
                                                                             step = 2),
                                                       ),
                                                       shiny::numericInput(inputId = "lfcMinSamseqInp",
                                                                           label = "Log Fold Change less or equal to",
                                                                           value = -2.0,
                                                                           step = 2),
                                                       shiny::numericInput(inputId = "lfcMaxSamseqInp",
                                                                           label = "Log Fold Change greater or equal to",
                                                                           value = 2.0,
                                                                           step = 2),
                                                       shiny::numericInput(inputId = "qValueSamseqInp",
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

                                                       shiny::selectInput(inputId = "controlDeseq2Inp",
                                                                          choices = c(""),
                                                                          label = "Reference level (control sample)"),

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
                                                                         selected = "biological",
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
                                                         label = "Human",
                                                         choices = c("FALSE","TRUE")),
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
                    shinybusy::use_busy_spinner(spin = "fading-circle"),
                    shiny::textOutput("execResults"))
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
    ),
    shiny::fluidRow(
      shiny::column(width = 10 , class="center",
                    shiny::wellPanel(
                      shiny::h3("Consensus Option", style="color: #5e63b6;"),
                      shiny::numericInput(inputId = "thresoldDe",
                                          label = "Minimum number of methods to indicate a gene as DE?",
                                          value = 5,
                                          min = 1,
                                          max = 7),
                      shiny::div( style= "text-align: center;",
                                  shiny::actionButton(inputId = "goConsList",
                                                      label = "View DE genes by consensus",
                                                      width = '51%',
                                                      class= "btn-success" ),
                      ),
                    ),
      )
    ),
    shiny::fluidRow(
      shiny::column(width = 12,
                    shiny::h2("Details of genes inidcated as DE ", style = "text-align: center;"),
                    shiny::helpText("Lines are genes, columns show values of each method executed in dataset"),
                    shiny::downloadButton("downloadData", "Download consensus result"),
                    shiny::wellPanel(
                      DT::dataTableOutput("tableConsensus")
                    )
      ),
    ),
    shiny::fluidRow(
      shiny::column(width = 10, class = "center",
                    shiny::h2("Heat Map of Differential Expressed Genes", style = "text-align: center;"),
                    shiny::helpText("Based in lines of table count, normalized by log2."),
                    shiny::div( style= "text-align: center;",
                                shiny::actionButton(inputId = "goHeatMap",
                                                    label = "View Heat Map to DE genes by consensus",
                                                    width = '51%',
                                                    class= "btn-success" ),
                    ),
                      shiny::wellPanel(
                        plotly::plotlyOutput("heatMapPlot")
                      )
                    )
    )
  ) #fluidPage



  # Define server logic required to draw a histogram
  server <- function(input, output, session) {
    options(shiny.maxRequestSize=30*1024^2)
    options('future.debug' = TRUE)
    consResult <- shiny::reactiveValues()
    deByTool <- NULL
    varCondNoiseq <- reactiveValues(cond = c(""))
    consListFinal <- shiny::reactiveValues()
    expDef_result <- shiny::reactiveValues()
    countData <- shiny::reactiveValues()

    groupUpdate <- shiny::eventReactive(input$groupNameInp, {
      groups = c(unlist(strsplit(input$groupNameInp, ",")))
      shiny::updateSelectInput(inputId = "controlDeseq2Inp",
                               choices = groups)
      if (length(groups) > 2) {
        shiny::updateRadioButtons(inputId = "condition1Inp",
                                  choices = groups,
                                  inline = TRUE)
        shiny::updateRadioButtons(inputId = "condition2Inp",
                                  choices = groups,
                                  inline = TRUE)
        varCondNoiseq$cond <- c(input$condition1Inp, input$condition2Inp)
        print(paste("---- condition1Noiseq - ", varCondNoiseq))
      } else {
        shiny::updateRadioButtons(inputId = "condition1Inp",
                                  choices = c(""))
        shiny::updateRadioButtons(inputId = "condition2Inp",
                                  choices = c(""))
        varCondNoiseq$cond <- c("")
      }
    })

    observeEvent(input$condition1Inp, {
      if (input$condition1Inp != "" && input$condition2Inp != "") {
        varCondNoiseq$cond <- c(input$condition1Inp, input$condition2Inp)
        print(paste("---- condition1Noiseq - ", varCondNoiseq$cond))
      }
    })

    observeEvent(input$condition2Inp, {
      if (input$condition1Inp != "" && input$condition2Inp != "") {
        varCondNoiseq$cond <- c(input$condition1Inp, input$condition2Inp)
        print(paste("---- condition1Noiseq - ", varCondNoiseq$cond))
      }
    })

    observeEvent(input$groupNameInp, {
      groupUpdate()
    })

    datasetCount <- shiny::eventReactive(input$go, {
      inFile <- input$tableCountInp
      countData$readData <- readCountFile(inFile$datapath, input$sepCharcterInp)
      qtdGroups <- length(c(unlist(strsplit(input$groupNameInp, ","))))
      informedColumns <- input$numberReplicsInp * qtdGroups
      colsReadData <- length(colnames(countData$readData))
      shiny::validate(
        need(informedColumns == colsReadData,
             label = paste("ERROR: infromed ",qtdGroups," groups of tretment and ",
                                             input$numberReplicsInp, "replics. File need be ",
                                             (input$numberReplicsInp * qtdGroups),
                                             "columns. But only ",colsReadData, "was founded."))
      )
      countData$readData
    })

    output$sample <- DT::renderDataTable({
      DT::datatable(datasetCount())

    })

    cons_res <- shiny::eventReactive(input$goDeg, {
      progress <- shiny::Progress$new()
      progress$set(message = "Computing expression", value = 0)
      on.exit(progress$close())

      updateProgress <- function(value = NULL, detail = NULL) {
        if (is.null(value)) {
          value <- progress$getValue()
          value <- value + (progress$getMax() - value) / 5
        }
        progress$set(value = value, detail = detail)
      }
      groupsName = c(unlist(strsplit(input$groupNameInp, ",")));
      print(paste("---- notHumanKnowseq - ",input$notHumanKnowseq ))
      consResult$exp <- runExpression(numberReplics = input$numberReplicsInp,
                                  rDataFrameCount = datasetCount(),
                                  groupName = groupsName,
                                  experimentName=input$experimentNameInp,
                                  printResults = FALSE,
                                  methodNormLimma = input$methNormLimmaInp,
                                  methodAdjPvalueLimma = input$adjPvalueLimma,
                                  numberTopTableLimma = input$topTableLimmaInp,
                                  respTypeSamseq = input$respTypeSamseqInp,
                                  npermSamseq = input$npermSamseq,
                                  fitTypeDeseq2 = input$fitTypeDeseq2Inp,
                                  controlDeseq2 = input$controlDeseq2Inp,
                                  methodNormEdgeR = input$methNormEdgerInp,
                                  normNoiseq = input$methNormNoiseqInp,
                                  kNoiseq = input$kNoiseqInp,
                                  factorNoiseq=input$factorNoiseqInp,
                                  lcNoiseq = input$lcNoiseqInp,
                                  replicatesNoiseq = input$replicatesNoiseqInp,
                                  condExpNoiseq = varCondNoiseq$cond,
                                  filterIdKnowseq=input$filterKnowseqInp,
                                  notSapiensKnowseq = as.logical(input$notHumanKnowseq),
                                  deNovoAanalysis = input$deNovoAanalysisInp,
                                  fdrEbseq=input$pValueEbseqInp,
                                  maxRoundEbseq = input$maxRoundEbseqInp,
                                  methodDeResultsEbseq = input$methodEbseqInp,
                                  progressShiny=updateProgress)

      return(consResult)
    })

    output$execResults <- shiny::renderText({
      consResult <- cons_res()
      print("DIFERENTIAL EXPRESSION ANALYSIS WAS COMPLETE!")
    })

    consensusPlot <- shiny::eventReactive(input$goUpsetPlot, {
      expDef_result$df <- expressionDefinition(resultTool = consResult$exp,
                                               groups = c(unlist(strsplit(input$groupNameInp, ","))),
                                            lfcMinLimma = input$lfcMinLimmaInp,
                                            lfcMaxLimma = input$lfcMaxLimmaInp,
                                            pValueLimma = input$pValueLimmaInp,
                                            lfcMinSamseq = input$lfcMinSamseqInp,
                                            lfcMaxSamseq = input$lfcMaxSamseqInp,
                                            qValueSamseq = input$qValueSamseqInp,
                                            scoreDSamseq = input$scoreDSamseqInp,
                                            lfcMinDeseq2 = input$lfcMinDeseq2Inp,
                                            lfcMaxDeseq2 = input$lfcMaxDeseq2Inp,
                                            pValueDeseq2 = input$pValueDeseq2Inp,
                                            lfcMinEdger = input$lfcMinEdgerInp,
                                            lfcMaxEdger = input$lfcMaxEdgerInp,
                                            pValueEdger = input$pValueEdgerInp,
                                            probNoiseq = input$probNoiseqInp,
                                            lfcMinKnowseq = input$lfcMinKnowseqInp,
                                            lfcMaxKnowseq = input$lfcMaxKnowseqInp,
                                            pValueKnowseq = input$pValueKnowseqInp,
                                            deClassEbseq = input$classDeEbseqInp)
      m <- as.matrix(consResult$exp)
      deByTool <- listDeByTool(consexpressionList = consResult$exp,
                               geneNames = row.names(m),
                               deList = expDef_result$df)

      if(length(deByTool > 0)){
        deByTool_filtered <- deByTool[, apply(deByTool, 2, function(col) sum(col) != 0)]
        UpSetR::upset(deByTool_filtered,
                      sets = colnames(deByTool_filtered), #(deByTool_filtered),
                      sets.bar.color = "#56B4E9",
                      order.by = "freq",
                      empty.intersections = "off")
      }else{
        print("ERROR: Problems to execute listDeByTool function.")
      }
    })

    output$upsetPlot <- shiny::renderPlot(consensusPlot(), res = 130)

    deConsList <- shiny::eventReactive(input$goConsList, {
      deByTool <- listDeByTool(consResult$exp, row.names(m),expDef_result$df)
      remove(m)
      consListFinal$dfList <- consensusList(consexpressionList = consResult$exp,
                                     deTool = deByTool,
                                     threshold = input$thresoldDe)
      if(length(consListFinal$dfList) <= 0){
        consListFinal$dfList <- paste("No genes identified as differentially expressed by ",
                                      input$thresoldDe ,
                                      " methods were found.")
      }
      return(consListFinal$dfList)
    })

    output$tableConsensus <- DT::renderDataTable({
      DT::datatable(as.data.frame(deConsList()))
    })

    output$downloadData <- shiny::downloadHandler(
      filename = function() {
        paste(input$tableConsensus, ".csv", sep = "")
      },
      content = function(file) {
        write.csv(consListFinal$dfList, file, row.names = FALSE)
      }
    )

    heatPlot <- shiny::eventReactive(input$goHeatMap, {
      req(input$tableCountInp)
      deList <- consListFinal$dfList
      countDf <- countData$readData
      index <- rownames(countDf) %in% rownames(deList$limma)
      subsetCountDF <- countDf[index, ]
      subsetDf <- as.data.frame(subsetCountDF)
      return(subsetDf)
    })

    output$heatMapPlot <- plotly::renderPlotly({
      df <- heatPlot()
      log_df <- log2(df)
      p <- plotly::plot_ly(x=colnames(df),
                           y=rownames(df),
                           z = as.matrix(df),
                           type = "heatmap",
                           # colors = "Blues"
                           )
      return(p)
    })

  }
  shiny::shinyApp(ui = ui, server = server)
}

#consexpressionR::consexpressionR()
