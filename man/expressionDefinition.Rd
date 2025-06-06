% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/expressionDefinition.R
\name{expressionDefinition}
\alias{expressionDefinition}
\title{Define who genes are deferentially expressed by method}
\usage{
expressionDefinition(
  resultTool,
  groups = c(""),
  lfcMinLimma = -2,
  lfcMaxLimma = 2,
  pValueLimma = 0.05,
  FLimma = 0.8,
  lfcMinSamseq = -2,
  lfcMaxSamseq = 2,
  qValueSamseq = 0.05,
  lfcMinDeseq2 = -2,
  scoreDSamseq = 0.8,
  lfcMaxDeseq2 = 2,
  pValueDeseq2 = 0.05,
  lfcMinEdger = -2,
  lfcMaxEdger = 2,
  pValueEdger = 0.05,
  probNoiseq = 0.8,
  lfcMinKnowseq = -2,
  lfcMaxKnowseq = 2,
  pValueKnowseq = 0.05,
  deClassEbseq = "DE",
  ppThresholdEbseq = 0.8,
  printResults = FALSE,
  pathOutput = "."
)
}
\arguments{
\item{resultTool}{table results by consexpressionR function}

\item{groups}{text, name of samples or treatment}

\item{lfcMinLimma}{minimum value to consider of log Fold Change (default: -2).}

\item{lfcMaxLimma}{maximum value to consider of log Fold Change (default: 2)}

\item{pValueLimma}{minimum P-Value to consider (default: 0.05)}

\item{FLimma}{The F statistics is used to test the null hypotheses of all groups have the same median of expression}

\item{lfcMinSamseq}{minimum value to consider of log Fold Change (default: -2).}

\item{lfcMaxSamseq}{maximum value to consider of log Fold Change (default: 2)}

\item{qValueSamseq}{q-value is a measure of the statistical significance of the difference in expression between the compared groups, taking the problem of multiple comparisons into account (default: 0.8)}

\item{lfcMinDeseq2}{minimum value to consider of log Fold Change (default: -2).}

\item{scoreDSamseq}{statistic used to evaluate the significance of each gene in multi class analysis}

\item{lfcMaxDeseq2}{maximum value to consider of log Fold Change (default: 2)}

\item{pValueDeseq2}{minimum P-Value to consider (default: 0.05)}

\item{lfcMinEdger}{minimum value to consider of log Fold Change (default: -2).}

\item{lfcMaxEdger}{maximum value to consider of log Fold Change (default: 2)}

\item{pValueEdger}{minimum P-Value to consider (default: 0.05)}

\item{probNoiseq}{floating point, minimum probability that a read count comes from a real gene expression peak, rather than background noiseq (default: 0.95)}

\item{lfcMinKnowseq}{minimum value to consider of log Fold Change (default: -2).}

\item{lfcMaxKnowseq}{maximum value to consider of log Fold Change (default: 2)}

\item{pValueKnowseq}{minimum P-Value to consider (default: 0.05)}

\item{deClassEbseq}{name of class to consider by EBSeq (default: "DE")}

\item{ppThresholdEbseq}{description}

\item{printResults}{logical variable: TRUE print report by each tool, FALSE print only consensus result}

\item{pathOutput}{path to write output, need be a directory (default: ".")
#'}
}
\description{
Define who genes are deferentially expressed by method
}
\examples{
data(gse95077)
treats = c("BM", "JJ")
cons_result <- runExpression(numberReplics = 3,
                              groupName = treats,
                              rDataFrameCount = gse95077,
                              sepCharacter = ",",
                              experimentName = "test_cons",
                              controlDeseq2 = "BM",
                              contrastDeseq2 = "JJ",
                              outDirPath = "." )
expDef_result <- expressionDefinition(resultTool = cons_result, groups = treats)
}
