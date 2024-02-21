# SRA010153_DE_Tool <- listDeByTool(consexpressionList = SRA010153_consexpression2, deList = SRA010153_deList)
# upSetPlotTools(SRA010153_DE_Tool, "SRA010153", pathOut = "/Volumes/SD128/consexpression2_testesOutput/SRA010153/", TRUE)
upSetPlotTools <- function(df,
                          condition = "condition_default",
                          pathOut = ".",
                          writeData = TRUE){
  df[is.na(df)] <- 0
  if(writeData){
    meu_grafico <- UpSetR::upset(df,
                                 sets = colnames(df),
                                 sets.bar.color = "#56B4E9",
                                 order.by = "freq",
                                 empty.intersections = "on")
    pdf(paste0(pathOut,
               "upset_plot_",
               condition ,
               ".pdf"),
        width = 10,
        height = 7.5,
        onefile = FALSE)
    print(meu_grafico)
    dev.off()
    df$`number of methods` <- rowSums(df)
    write.csv(df,
              file = paste0(pathOut,
                           "upset_plot_",
                           condition ,
                           ".csv"))
  }else{
    df$`number of methods` <- rowSums(df)
    UpSetR::upset(df,
                  sets = colnames(df),
                  sets.bar.color = "#56B4E9",
                  order.by = "freq",
                  empty.intersections = "off",
                  main = paste("upset_plot_",
                               condition))
  }
  return(df)
}

