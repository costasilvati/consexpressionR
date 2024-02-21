#--- GSE95077 pipeline

GSE95077_consexpression2 <- consexpression2(numberReplics = 3, groupName = c("BM", "JJ"),tableCountPath = "data/GSE95077_filtred.csv",sepCharacter = ",",experimentName = "GSE95077",outDirPath = "/Volumes/SD128/consexpression2_testesOutput/GSE95077", printResults = TRUE)
GSE95077_deList <- expressionDefinition(GSE95077_consexpression2)

toolNames <- names(GSE95077_deList)
geneNames <- GSE95077_qRTPCR[[2]]
t <- 1
for (toolDe in GSE95077_deList) {
  genesDe <- rownames(toolDe)
  for (i in 1:length(geneNames)) {
    if(geneNames[i] %in% genesDe){
      linhaPCR <- which(GSE95077_qRTPCR$EnsemblID == geneNames[i])
      if(GSE95077_qRTPCR[linhaPCR,1] == "Selected"){
        GSE95077_DE_Tool[geneNames[i], t] <- 1
      }
    } else {
      print("Termo não encontrado.")
    }
  }
  t <- t + 1
}


#--- SRA010153 pipeline

SRA010153_consexpression2 <- consexpression2(numberReplics = 7, groupName = c("UHR", "Brain"),tableCountPath = "data/SRA010153_filtred.csv",sepCharacter = ",",experimentName = "SRA010153",outDirPath = "/Volumes/SD128/consexpression2_testesOutput/SRA010153/", printResults = TRUE)

SRA010153_deList <- expressionDefinition(SRA010153_consexpression2, lfcMax = 2, lfcMin = -2, pValue = 0.05)

SRA010153_DE_Tool <- data.frame(matrix(0, ncol = length(SRA010153_deList),
                                       nrow = nrow(SRA010153_qRTPCR)))

colnames(SRA010153_DE_Tool) <- names(SRA010153_deList)
row.names(SRA010153_DE_Tool) <- SRA010153_qRTPCR$Gene
toolNames <- names(SRA010153_deList)
geneNames <- SRA010153_qRTPCR[[1]]
t <- 1
for (toolDe in SRA010153_deList) {
  genesDe <- rownames(toolDe)
  for (i in 1:length(geneNames)) {
    if(geneNames[i] %in% genesDe){
      linhaPCR <- which(SRA010153_qRTPCR$Gene == geneNames[i])
      if(SRA010153_qRTPCR[linhaPCR,2] == 1 | SRA010153_qRTPCR[linhaPCR,2] == -1){
        SRA010153_DE_Tool[geneNames[i], t] <- 1
      }
    } else {
      cat(".")
    }
  }
  t <- t + 1
}
remove(toolNames, geneNames, toolDe, genesDe, i, linhaPCR, t)

# upSetPlotTools()

SRA010153_consensus <- list()
namesCons <- c()
for (i in 1: length(colnames(SRA010153_DE_Tool))) {
  namesCons[i] <- paste("Number os Tools:" ,i)
  SRA010153_consensus[[i]] <- consensusList(SRA010153_DE_Tool, i)
}
names(SRA010153_consensus) <- namesCons
