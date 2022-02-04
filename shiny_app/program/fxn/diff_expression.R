# Identifies differentially expressed genes between two groups of cells
# Returns a p-value ranked matrix of putative differentially expressed features

suppressPackageStartupMessages(library(stats))
suppressPackageStartupMessages(library(pbapply))
suppressPackageStartupMessages(library(limma))


Fast.DE2 <- function(obj.use, cells.1 = NULL, cells.2 = NULL,
                     logfc.threshold  = 0.5,
                     pseudocount.use  = 1,
                     num.gene2test    = 1000,
                     min.pct          = 0.1,
                     min.pct.diff     = 1.25) {

    diff_genes <- NULL
    min.pct    <- min.pct*100

    cells.exp <- get_cell_filtered_expression(obj.use, c(cells.1, cells.2))

    cells.exp <- cells.exp %>%
        mutate(Cell_Cluster = ifelse(Cell %in% cells.1, 1, 2))

    genes.use <- setdiff(colnames(cells.exp), c('Cell', 'Cell_Cluster'))

    data.1 <- apply(cells.exp[cells.exp$Cell_Cluster == 1, ],
                    MARGIN = 2,
                    FUN = function(x) log(mean(expm1(x)) + pseudocount.use))
    data.2 <- apply(cells.exp[cells.exp$Cell_Cluster == 2, ],
                    MARGIN = 2,
                    FUN = function(x) log(mean(expm1(x)) + pseudocount.use))

    # Genes not expressed in either group
    genes.zeros <- union(names(which(data.1 == 0)), names(which(data.2 == 0)))

    # Fold change greater than threshold
    total.diff <- (data.1 - data.2)
    genes.diff <- names(x = which(x = abs(x = total.diff) > logfc.threshold))
    genes.diff <- union(genes.diff, genes.zeros)
    genes.use <- intersect(x = genes.use, y = genes.diff)

    # Remove ribisomal and mitochondria genes: "^RP[SL]" and "^MT-"
    idx_mt_rp <- grep(pattern = "^MT-|^RP[SL]", x = genes.use, ignore.case = T)
    if (any(idx_mt_rp)) {
        genes.use <- genes.use[-idx_mt_rp]
    }

    # Remove genes if cells% < min.cells in BOTH clusters
    pct.1 <- apply(cells.exp[cells.exp$Cell_Cluster == 1, ],
                   MARGIN = 2,
                   FUN = function(x) round(sum(x > 0)/length(cells.1)*100,2))
    pct.2 <- apply(cells.exp[cells.exp$Cell_Cluster == 2, ],
                   MARGIN = 2,
                   FUN = function(x) round(sum(x > 0)/length(cells.2)*100,2))

    genes.use <- genes.use[!((pct.1 < min.pct & pct.2 < min.pct) |
                                 (pmax(pct.1 / pct.2, pct.2 / pct.1) < min.pct.diff))]

    # Select maximum top N genes for DE gene test
    if (length(x = genes.use) > num.gene2test) {
        genes.use <- names(head(sort(abs(total.diff[genes.use]),decreasing = T), num.gene2test))
    }

    if (length(genes.use) > 0) {
        data.use <- t(cells.exp[, genes.use])
        colnames(data.use) <- cells.exp$Cell

        if (NROW(data.use) > 0) {
            p_val <- sapply(1:NROW(x = data.use),
                            function(x) {
                                min(2 * min(rankSumTestWithCorrelation(index = 1:length(cells.1), statistics = data.use[x, ])), 1)
                            }
            )

            # Add genes as the first column
            diff_genes <- data.frame(gene = rownames(x = data.use), p_val = p_val, row.names = rownames(x = data.use))

            # Append log FC and Bonferroni correction
            diff_genes$log_FC    <- round(total.diff[genes.use],2)
            diff_genes$p_val_adj <- p.adjust(diff_genes$p_val,
                                             method = "bonferroni",
                                             n = NROW(data.use))

            # Compute positive rate
            diff_genes$pct.1     <- pct.1[genes.use]
            diff_genes$pct.2     <- pct.2[genes.use]

            diff_genes$eff[diff_genes$log_FC > 0] <- "up.reg"
            diff_genes$eff[diff_genes$log_FC < 0] <- "down.reg"

            # Sort by log_FC
            diff_genes <- diff_genes[order(diff_genes$log_FC,decreasing = T),]
        }
    }

    diff_genes
}
