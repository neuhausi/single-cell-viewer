# Differential expression using Wilcoxon Rank Sum
#
# Modified from Seurat v2.3.4.
# Applicable to Seurat V3 objects
# Shuoguo Wang @ 05/11/2019
#
# Identifies differentially expressed genes between two groups of cells using
# a Wilcoxon Rank Sum test
#
# @param data.use Data matrix to test
# @param cells.1 Group 1 cells
# @param cells.2 Group 2 cells
#
# @return Returns a p-value ranked matrix of putative differentially expressed
# features
#
#' @importFrom pbapply pbsapply
#' @importFrom stats wilcox.test

suppressPackageStartupMessages(library(stats))
suppressPackageStartupMessages(library(pbapply))

Fast.WC.DE <- function(obj.use, cells.1 = NULL, cells.2 = NULL, out_csv = NULL,
                       logfc.threshold  = 0.5,
                       pseudocount.use  = 1,
                       num.gene2test    = 1000,
                       min.pct          = 0.1,
                       min.pct.diff     = 1.25) {

    diff_genes <- NULL

    # Get data, and all genes as candidate
    data.use  <- GetAssayData(obj.use, assay = obj.use@active.assay)
    genes.use <- rownames(data.use)
    min.pct   <- min.pct*100

    # Log mean -- adapted from seurat::differential_expression.R, lines 175-181
    data.1 <- apply(X = data.use[genes.use, cells.1, drop = F],
                    MARGIN = 1,
                    FUN = function(x) log(x = mean(x = expm1(x = x)) + pseudocount.use))
    data.2 <- apply(X = data.use[genes.use, cells.2, drop = F],
                    MARGIN = 1,
                    FUN = function(x) log(x = mean(x = expm1(x = x)) + pseudocount.use))

    # Genes not expressed in one or the other groups
    genes.zeros <- union(names(which(data.1 == 0)), names(which(data.2 == 0)))

    # Fold change greater than threshold [0.5]
    total.diff <- (data.1 - data.2)
    genes.diff <- names(x = which(x = abs(x = total.diff) > logfc.threshold))
    genes.diff <- union(genes.diff, genes.zeros)
    genes.use <- intersect(x = genes.use, y = genes.diff)

    # Remove ribisomal and mitochondria genes: "^RP[SL]" and "^MT-"
    idx_mt_rp <- grep(pattern = "^MT-|^RP[SL]", x = genes.use)
    if (any(idx_mt_rp)) {
        genes.use <- genes.use[-idx_mt_rp]
    }else{
        genes.use <- genes.use
    }

    # Remove genes if cells% < min.cells in BOTH clusters
    pct.1 <- apply(data.use[genes.use,cells.1],1, function(x) round(sum(x > 0)/length(cells.1)*100,2))
    pct.2 <- apply(data.use[genes.use,cells.2],1, function(x) round(sum(x > 0)/length(cells.2)*100,2))
    genes.use <- genes.use[!((pct.1 < min.pct & pct.2 < min.pct) |
                                 (pmax(pct.1 / pct.2, pct.2 / pct.1) < min.pct.diff) )]

    # Select maximum top N genes for DE gene test
    if (length(x = genes.use) > num.gene2test  ) {
        genes.use <- names(head(sort(abs(total.diff[genes.use]),decreasing = T), num.gene2test))
    }

    if (!is.null(cells.1) && !is.null(cells.2) && (length(cells.1) > 0) && (length(cells.2) > 0)) {

        # Compute p_val using Wilcoxon Rank Sum test
        group.info <- data.frame(row.names = c(cells.1, cells.2))
        group.info[cells.1, "group"] <- "Group1"
        group.info[cells.2, "group"] <- "Group2"
        group.info[, "group"] <- factor(x = group.info[, "group"])
        data.use <- data.use[genes.use, rownames(x = group.info), drop = FALSE]
        p_val <- pbsapply(
            X = 1:nrow(x = data.use),
            FUN = function(x) {
                return(wilcox.test(data.use[x, ] ~ group.info[, "group"], data = data.use)$p.value)
            }
        )

        # Add genes as the first column
        # Connie - Modified 5/13
        diff_genes <- data.frame(gene = rownames(x = data.use), p_val = p_val, row.names = rownames(x = data.use))

        # Append log FC and Bonferroni correction
        diff_genes$log_FC    <- round(total.diff[genes.use],2)
        diff_genes$p_val_adj <- p.adjust(diff_genes$p_val,
                                         method = "bonferroni",
                                         n = nrow(data.use))

        # Compute positive rate and odds
        diff_genes$pct.1     <- pct.1[genes.use]
        diff_genes$pct.2     <- pct.2[genes.use]
        # Removing the below two fields as they are not in the original calculation/app
        # Connie - Modified 5/13
        # diff_genes$odds.up   <- round(diff_genes$pct.1/diff_genes$pct.2,2)
        # diff_genes$odds.down <- round(diff_genes$pct.2/diff_genes$pct.1,2)
        diff_genes$eff[diff_genes$log_FC > 0] <- "up.reg"
        diff_genes$eff[diff_genes$log_FC < 0] <- "down.reg"

        # Sort by log_FC
        diff_genes <- diff_genes[order(diff_genes$log_FC,decreasing = T),]
    }

    # output
    if (!is.null(out_csv)) {
        write.csv(diff_genes, file = out_csv)
    }else{
        diff_genes
    }
}
