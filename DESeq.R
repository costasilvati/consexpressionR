# Title     : run DESeq
# Objective : execute DESeq by informations user
# Created by: julianacostasilva
# Created on: 10/03/22

# Fonte: https://github.com/labbcb/rnaseq/wiki/Performing-RNA-Seq-Data-Analysis-with-Bioconductor

# pkgs = c('XML', 'DESeq', 'DESeq2', 'pasilla', 'gplots', 'RColorBrewer')

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(pkgs)
# BiocManager::install()
install.packages("DESeq", repos = NULL, source("DESeq.tgz"))
library("DESeq")
table_path = '/Users/julianacostasilva/OneDrive/ProjetoDoutorado/2021-2/bioconvergencia/reads_RNApa/kallisto_quant_RNApa_asca_0B_pb_tpm_tabular.csv'
separator = ','

countData <- read.table(table_path, sep = separator, row.names = 1, header = TRUE, stringsAsFactors=FALSE)
# ,gon0b1,mus0b1,gon1b4,mus1b4,gon0b5,mus0b5,gon1b6,mus1b6,gon1b8,mus1b8,gon0b9,mus0b9
apaDesign <- data.frame(row.names = colnames(countData),
                        condition = c( "gon0b", "mus0b", "gon1b", "mus1b", "gon0b", "mus0b", "gon1b", "mus1b", "gon1b", "mus1b","gon0b", "mus0b"),
                        libType = c( "paired-end", "paired-end", "paired-end", "paired-end", "paired-end", "paired-end","paired-end", "paired-end", "paired-end", "paired-end", "paired-end", "paired-end") )

pairedSamples <- apaDesign$libType == "paired-end"
countTable <- countData[ , pairedSamples ]
conds <- apaDesign$condition[ pairedSamples ]
head(countTable)
cds <- newCountDataSet( countTable, conds )

