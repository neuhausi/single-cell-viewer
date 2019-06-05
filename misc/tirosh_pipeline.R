#!/bin/env Rscript

# Shuoguo Wang
# 2019/05/02
# Copyright: Bristol Myers Squibb


# Source Data Files:
#   - melanoma_expression.txt
#   - melanoma_cluster_assignment_portal.txt
#   - melanoma_coordinates_portal.txt
#
# Can all be downloaded from Broad Institute:
# https://portals.broadinstitute.org/single_cell/study/SCP11/melanoma-intra-tumor-heterogeneity


### -------- functions --------- ###

replace.minus <- function(input_list) {
    input_list.collapse <- paste(unlist(input_list), collapse = "|")
    input_list.collapse.corr <- gsub("-", replacement = "_", x = input_list.collapse)
    input_list.corr <- strsplit(x = input_list.collapse.corr,split = "[|]")
    unlist(input_list.corr)
}



### ---------- setup ----------- ###

library(Seurat)
library(plyr)
library(dplyr)
data_loc   <- "./data/"
output_loc <- "./"

if(packageVersion("Seurat") < "3.0.0") { stop("Need to upgrade Seurat to >= version 3.0.0") }


### --------- pipeline --------- ###


## Load normalized data ##
df.data.temp  <- read.csv(paste0(data_loc, "melanoma_expression.txt"),
                          row.names = 1,
                          header = T,
                          check.names = F,
                          sep = "\t")
## Load classification ##
df.ident.temp <- read.csv(paste0(data_loc, "melanoma_cluster_assignment_portal.txt"),
                          row.names = 1,
                          skip = 1,
                          col.names = c("Name","Cluster","Sub_Cluster"),
                          check.names = F,
                          header = T,
                          sep = "\t")
## Load coordinates ##
df.tsne       <- read.csv(paste0(data_loc, "melanoma_coordinates_portal.txt"),
                          row.names = 1,
                          skip = 1,
                          col.names = c("Name","tSNE_1","tSNE_2"),
                          check.names = F,
                          header = T,
                          sep = "\t")

## Subset immune cells using TSNE cell names
cells <- row.names(df.tsne)
idx <- (row.names(df.ident.temp) %in% cells)
df.ident <- df.ident.temp[idx,]
idx <- names(df.data.temp) %in% cells
df.data <- df.data.temp[,idx]

## Fix cell_names which has "-"
row.names(df.ident) <- replace.minus(row.names(df.ident))
names(df.data)  <- replace.minus(names(df.data))
row.names(df.tsne)    <- replace.minus(row.names(df.tsne))

# Create seurat 3 object
myproject <- CreateSeuratObject(counts = df.data)
myproject@active.assay <- "RNA"

## Add meta.info
meta.info.tmp <- read.csv(paste0(data_loc, "Tirosh-MM-2016-CD45P.csv"),
                          check.names = F,
                          row.names = 1,
                          na.strings = F,
                          stringsAsFactors = F)
meta.info <- as.list(subset(meta.info.tmp, select = "Description", drop = T))
names(meta.info) <- row.names(meta.info.tmp)
myproject@misc$meta.info <- meta.info

## set cell ident
cell_ident <- droplevels(df.ident$Cluster)
names(cell_ident) <- rownames(df.ident)
myproject[["Tirosh.ident"]] <- cell_ident
Idents(myproject) <- myproject[["Tirosh.ident"]]

## set dr
dr_tsne <- CreateDimReducObject(embeddings = as.matrix(df.tsne),assay = "RNA",key = "tSNE_")
myproject[['tsne']] <- dr_tsne

## compute DEG
myproject@misc$all_genes <- rownames(myproject)
myproject@misc$DE$all <- FindAllMarkers(object = myproject,
                                        only.pos = FALSE,
                                        min.pct = 0.1)
myproject@misc$DE$all %>% group_by(cluster) %>% top_n(10, avg_logFC) %>% as.data.frame() -> myproject@misc$DE$top10
myproject@misc$DE$all %>% group_by(cluster) %>% top_n(30, avg_logFC) %>% as.data.frame() -> myproject@misc$DE$top30

## add meta
## The following meta.data was manually curirated from Tirosh supplementary table S1 and S2:
sample_id <- table(tolower(unlist(lapply(colnames(myproject),
                                         function(x) {substr(x,1,4)[[1]][[1]]}))))
sample_list <- tolower(unlist(lapply(colnames(myproject),
                                     function(x) {substr(x,1,4)[[1]][[1]]})))
PatientID <- mapvalues(sample_list,
                       from = names(sample_id),
                       to   = c("Mel53","Mel58","Mel60","Mel72","Mel74","Mel79","Mel80","Mel81","Mel84","Mel88","Mel89","Mel94"))
Mutation <- mapvalues(sample_list,
                      from = names(sample_id),
                      to   = c("Wild-type","Wild-type","BRAF","NRAS","na","Wild-type","NRAS","BRAF","Wild-type","NRAS","na","Wild-type"))
Sex <- mapvalues(sample_list,
                 from = names(sample_id),
                 to   = c("F","F","M","F","M","M","F","F","M","M","M","F"))
PreTreatment <- mapvalues(sample_list,
                          from = names(sample_id),
                          to = c("N","Y","Y","Y","Y","N","N","N","N","Y","N","Y"))
names(PatientID) <- colnames(myproject)
names(Mutation) <- colnames(myproject)
names(Sex) <- colnames(myproject)
names(PreTreatment) <- colnames(myproject)

myproject@meta.data$Mutation <- Mutation
myproject@meta.data$Sex <- Sex
myproject@meta.data$PreTreatment <- PreTreatment
myproject@meta.data$PatientID <- PatientID

myproject@misc$DataSegregation <- list("Mutation" = names(table(myproject@meta.data$Mutation)),
                                       "Sex" = names(table(myproject@meta.data$Sex)),
                                       "PatientID" = names(table(myproject@meta.data$PatientID)),
                                       "PreTreatment" = names(table(myproject@meta.data$PreTreatment)))

saveRDS(myproject, paste0(output_loc, "tirosh_seurat3.RDS"))
message('Finished, file output to: ', paste0(output_loc, "tirosh_seurat3.RDS"))
