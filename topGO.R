library(topGO)
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
library(data.table)
library(tidyverse)

GO_terms <- readMappings("~/Desktop/fly_to_GO.delim")
#gtf <- rtracklayer::import("~/Desktop/dmel-all-r6.23.gtf")
genes_of_interest <- read.table("~/Desktop/SNPs_of_interest_names.txt")
genes_of_interest <- genes_of_interest[,2]


remove <- grep('-', genes_of_interest)
genes_of_interest <- genes_of_interest[-remove]

allgenes <- data.frame(GenomicFeatures::genes(TxDb.Dmelanogaster.UCSC.dm6.ensGene))
fbgn <- allgenes$gene_id
present <- as.integer(fbgn %in% genes_of_interest)

gene_list <- setNames(as.numeric(present), fbgn)

gene_filter <- function(allScore){
  return(allScore == 1)
}


data <- new("topGOdata",
            ontology = "BP",
            allGenes = gene_list,
            geneSel = gene_filter,
            annotationFun = annFUN.gene2GO,
            nodeSize = 10,
            #geneSelectionFun = gene_filter,
            gene2GO = GO_terms)


resultFisher <- runTest(data, algorithm = "weight01", statistic = "fisher")

allRes <- GenTable(data, classicFisher = resultFisher, topNodes = 20)

showSigOfNodes(data, score(resultFisher), firstSigNodes = 5, useInfo = 'all')


