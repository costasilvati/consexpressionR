#' consexpressionR: Differential Expression Analysis using consensus of multiple methods
#'
#' Comprehensive workflow for differential expression (DE) analysis using multiple
#' methods (e.g., DESeq2, edgeR, limma, EBSeq, NOISeq, SAMSeq, KnowSeq) and
#' consensus strategies to derive robust DE gene lists.
#'
#' @docType package
#' @name consexpressionR
#'
#' @examples
#' obj <- createExpressionResultSet(
#'   results = list(edgeR = data.frame(gene = "g1", logFC = 1)),
#'   methodNames = "edgeR"
#' )
#' obj
"_PACKAGE"

## usethis namespace: start
#' @importFrom attempt is_try_error
#' @importFrom cqn cqn
#' @importFrom DT dataTableOutput
#' @importFrom mclust mclustBIC
#' @importFrom plotly renderPlotly
#' @importFrom shiny fileInput
#' @importFrom shiny fluidRow
#' @importFrom shiny icon
#' @importFrom shiny need
#' @importFrom shiny numericInput
#' @importFrom shiny plotOutput
#' @importFrom shiny reactive
#' @importFrom shiny renderPrint
#' @importFrom shiny selectInput
#' @importFrom shiny shinyApp
#' @importFrom shiny tags
#' @importFrom shiny textInput
#' @importFrom shiny textOutput
#' @importFrom shiny validate
#' @importFrom shinydashboard box
#' @importFrom shinydashboard dashboardBody
#' @importFrom shinydashboard dashboardHeader
#' @importFrom shinydashboard dashboardPage
#' @importFrom shinydashboard dashboardSidebar
#' @importFrom shinydashboard menuItem
#' @importFrom shinydashboard sidebarMenu
#' @importFrom shinydashboard tabItem
#' @importFrom shinydashboard tabItems
#' @importFrom shinyjs disable
#' @importFrom shinyjs enable
#' @importFrom stats relevel
#' @importFrom stats var
#' @importFrom UpSetR upset
#' @importFrom utils write.csv
