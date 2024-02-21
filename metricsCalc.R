#' GSE95077_metrics2 <- metricsCalc(deToolList = GSE95077_deList, goldList = GSE95077_qRTPCR, pathOut = "/Volumes/SD128/consexpression2_testesOutput/data/GSE95077/", experimentName = "GSE95077", goldColGeneName = 2, goldColValue = 1, goldPositiveValue = c("Selected"), goldNegativeValue = "NO")
#' SRA010153_metrics2 <- metricsCalc(deToolList = SRA010153_deList, goldList = SRA010153_qRTPCR, pathOut = "/Volumes/SD128/consexpression2_testesOutput/data/SRA010153/", experimentName = "SRA010153")
metricsCalc <- function(deToolList,
                    goldList,
                    goldColGeneName = 1,
                    goldColValue = 2,
                    goldPositiveValue = c(1,-1),
                    goldNegativeValue = 0,
                    experimentName = "experiment",
                    pathOut = ".",
                    writeData = FALSE
                    ){
  cols = c("TP","FP", "TN", "FN", "TPR","Specificity","PPV", "ACC", "F1-Score")
  toolNames <- names(deToolList)
  metrics <- data.frame(matrix(0, ncol = length(cols),
                                        nrow = length(toolNames)))
  colnames(metrics) <- cols
  row.names(metrics) <- toolNames

  for (i in 1:length(toolNames)) {
    cat(i,"\n")
    condicao_TP <- (goldList[,goldColValue] %in% goldPositiveValue) & (goldList[,goldColGeneName] %in% row.names(deToolList[[i]]))
    tp <- sum(condicao_TP)
    condicao_FP <- (goldList[,goldColValue] == goldNegativeValue) & (goldList[,goldColGeneName] %in% row.names(deToolList[[i]]))
    fp <- sum(condicao_FP)
    condicao_TN <- (goldList[,goldColValue] == goldNegativeValue) & (!(goldList[,goldColGeneName] %in% row.names(deToolList[[i]])))
    tn <- sum(condicao_TN)
    condicao_FN <- (goldList[,goldColValue] %in% goldPositiveValue) & (!(goldList[,goldColGeneName] %in% row.names(deToolList[[i]])))
    fn <- sum(condicao_FN)
    #TPR (True Positive Rate), SPC (Specificity), PPV (Positive Predict Value), ACC (Accuracy) and F1 measure
    tpr <- tp/(tp+fn)
    spc <- tn/(tn+fp)
    ppv <- tp/(tp+fp)
    acc <- (tp + tn)/(tp+fp+tn+fn)
    f1score <- 2*(ppv * tpr)/(ppv+tpr)
    metrics[i,] <- c(tp, fp, tn, fn, tpr, spc, ppv, acc, f1score)
  }
  if(writeData){
    write.csv(metrics, file = paste0(pathOut,experimentName,"_metrics.csv"))
  }
  return(metrics)
}


