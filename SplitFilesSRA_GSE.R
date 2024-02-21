colunas_a_manter <- grep("JJ", names(GSE95077_TG003), invert = TRUE)
GSE95077_TG003_KM <- GSE95077_TG003[, colunas_a_manter]

colunas_a_manter <- grep("BM", names(GSE95077_TG003), invert = TRUE)
GSE95077_TG003_JJ <- GSE95077_TG003[, colunas_a_manter]

names(GSE95077__JJN3_and_CTRL)[1] <- "GeneID"
names(GSE95077_KMS12BM_and_CTRL)[1] <- "GeneID"
names(GSE95077_TG003)[1] <- "GeneID"
names(GSE95077_TG003_JJ)[1] <- "GeneID"
names(GSE95077_TG003_KM)[1] <- "GeneID"

GSE95077__JJN3_and_CTRL <- merge(GSE95077__JJN3_and_CTRL,
                                 GSE95077_TG003_JJ,
                                 by = "GeneID",
                                 all.x = TRUE)


GSE95077_KMS12BM_and_CTRL <- merge(GSE95077_KMS12BM_and_CTRL,
                                 GSE95077_TG003_KM,
                                 by = "GeneID",
                                 all.x = TRUE)

#---- Read counts RAW
# Specify the path to your folder
folder_path <- "~/Downloads/GSE95077_RAW/"

# Get a list of all files in the folder
file_list <- list.files(folder_path, pattern = "\\.gz$", full.names = TRUE)
pattern <- "(_)[A-Z]{2}(_)([A-Z]{4}|[A-Z]{2})(_)[0-9]{6}"
treats <- substring(regmatches(file_list, regexpr(pattern, file_list)), 2)
treats
library(data.table)
# Read all .gz files into a list
GSE95077_raw_count_data <- lapply(file_list, function(file_path) {
  con <- gzfile(file_path, "rb")
  data <- read_tsv(con, col_names = FALSE)
  close(con)
  return(data)
})
names(GSE95077_raw_count_data) <- treats
# Renomeando colunas
for (i in 1:length(treats)) {
  colnames(GSE95077_raw_count_data[[i]]) <- c("GeneID", treats[i])
}
GSE95077_count_df <- data.frame()
GSE95077_count_df <- Reduce(function(x, y) merge(x, y,
                                                 by = "GeneID",
                                                 all = TRUE), GSE95077_raw_count_data)
save(GSE95077_raw_count_data, file = "data/GSE95077_raw_count_data.RData")
rownames(GSE95077_count_df) <- GSE95077_count_df$GeneID
GSE95077_count_df <- GSE95077_count_df[,-1]

# separando BM
pattern <- "^BM(_)"
matching_columns <- grep(pattern, colnames(GSE95077_count_df))
GSE95077_count_BM <- GSE95077_count_df[, matching_columns]

# separando JJ
pattern <- "^JJ(_)"
matching_columns <- grep(pattern, colnames(GSE95077_count_df))
GSE95077_count_JJ <- GSE95077_count_df[, matching_columns]


# --- extraindo apenas os genes avaliados GSE95077
library(tidyverse)
qRT_PCR_selected_genes <- read.csv("~/Library/CloudStorage/OneDrive-Personal/Doutorado/consexpression_2/datasets/GSE95077_RAW/qRT_PCR_selected_genes.csv")
save(qRT_PCR_selected_genes, file = "data/qRT_PCR_selected_genes.RData")

colnames(qRT_PCR_selected_genes)[2] <- "EnsemblID"
colnames(qRT_PCR_selected_genes)[1] <- "select"
qrtpcr_filtered <- qRT_PCR_selected_genes %>%
  filter(select == "Selected")

genes_qrtpcr <- qrtpcr_filtered$EnsemblID
GSE95077_filtered <- GSE95077_count_df[genes_qrtpcr,]
save(GSE95077_count_df, file = "data/GSE95077_count.RData")
save(GSE95077_filtered, file = "data/GSE95077_selected.RData")

GSE95077_filtered_BMJJ <- subset(GSE95077_filtered, select = c(4,5,6,10,11,12))
save(GSE95077_filtered_BMJJ, file = "/Volumes/SD128/consexpression2_testesOutput/data/GSE95077_filtred_BMJJ.RData")

remove(colunas,
       colunas_a_manter,
       linhas,
       matching_columns,
       GSE95077_count_BM,
       GSE95077_count_df,
       GSE95077_count_JJ,
       GSE95077_TG003_KM,
       GSE95077_raw_count_data,
       GSE95077__JJN3_and_CTRL,
       GSE95077_KMS12BM_and_CTRL,
       qRT_PCR_selected_genes,
       count_data,
       pattern,
       folder_path,
       qrtpcr_filtered,
       genes_qrtpcr)
# --- extraindo apenas os genes avaliados SRA010153
SRA010153_count <- as.data.frame(read.csv("data/SRA010153_count.txt"))
SRA010153_selected <- read.csv("data/SRA010153_selected.txt")
# ---- Remover .[0-9] no final do gene ID
gene_names <- row.names(SRA010153_count)
row.names(SRA010153_count) <- gsub("\\..*", "", gene_names)
SRA010153_filtred <- SRA010153_count[SRA010153_selected$Gene,]
# Ao escrever o .csv a linha 499 fica toda com NA (remover a linha)

save(SRA010153_count, file = "data/SRA010153_count.RData")
save(SRA010153_selected, file = "data/SRA010153_selected.RData")
remove(genes_sra,
       gene_names,
       SRA010153_count,
       SRA010153_selected)
