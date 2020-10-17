
#' Normalize matrix using knn-matrix
#' 
#' @param mat matrix to normalize
#' @param knn_mat knn_matrix to normalize with
#' 
#' @return normalized matrix
#' 
#' 
#' @export
normalize_mat_knn <- function(mat, knn_mat) {
    mat <- mat[, colnames(knn_mat)]
    mat_filt_na <- mat
    mat_filt_na[is.na(mat)] <- 0
    mat_exp <- as.matrix(mat_filt_na) %*% as.matrix(t(knn_mat))
    not_na_mat <- !is.na(mat)
    mat_exp_n <- as.matrix(not_na_mat) %*% as.matrix(t(knn_mat))
    mat_exp_norm <- mat_exp / mat_exp_n
    mat_oe <- mat - mat_exp_norm
    return(mat_oe)
}


#' Define immune and CAF signatures 
#' 
#' @param em_clust output of \code{cluster_em_cross_cor}
#' @param scale_feats scale features (recommended)
#' @param caf_gene name of gene that is inside the CAF cluster (default: "CAV1")
#' @param immune_gene name of gene that is inside the immune cluster (default: "CD3D")
#' 
#' @return matrix with 'CAF' and 'immune' score for each sample
#' 
#' @export
immune_cafs_signature_from_clust <- function(em_clust, scale_feats = TRUE, caf_gene = "CAV1", immune_gene = "CD3D") {
    expr_clust <- em_clust$expr_clust

    caf_mod <- expr_clust %>%
        filter(grepl(caf_gene, name)) %>%
        slice(1) %>%
        pull(clust)
    immune_mod <- expr_clust %>%
        filter(grepl(immune_gene, name)) %>%
        slice(1) %>%
        pull(clust)

    stopifnot(nrow(expr_clust %>% filter(clust == immune_mod, grepl(immune_gene, name))) > 0)
    stopifnot(nrow(expr_clust %>% filter(clust == caf_mod, grepl(caf_gene, name))) > 0)

    message(glue("immune module: {immune_mod}"))
    message(glue("cafs module: {caf_mod}"))

    feats <- em_clust$expr_mods[, c(caf_mod, immune_mod)]
    colnames(feats) <- c("caf", "immune")
    
    if (scale_feats) {
        feats <- scale(feats, center = TRUE, scale = TRUE)
    }

    return(feats)
} 


#' Calculate knn-matrix based on feature matrix
#' 
#' @param feat feature matrix (e.g. output of \code{define_immune_cafs_signature})
#' @param k k parameter for knn
#' 
#' @return knn-matrix: number of samples X number of features with 1 for being knn in the dimesion and 0 otherwise.
#' 
#' 
#' @export
feats_to_knn_matrix <- function(feats, k) {
    dist_mat <- tgs_dist(feats)    
    knn_df <- tgs_knn(100 - as.matrix(dist_mat), k)
    knn_df <- knn_df %>% mutate(col1 = factor(col1), col2 = factor(col2, levels = levels(col1)))

    knn_mat <- Matrix::sparseMatrix(as.numeric(knn_df$col1), as.numeric(knn_df$col2), x = 1)
    rownames(knn_mat) <- levels(knn_df$col1)
    colnames(knn_mat) <- levels(knn_df$col2)
    return(knn_mat)
}

calc_norm_immune_cafs_meth <- function(meth_mat, k = NULL, feats = NULL, knn_mat = NULL) {
    if (is.null(knn_mat)) {
        if (is.null(feats) || is.null(k)) {
            stop("Please provide either feats and k or knn_mat")
        }
        knn_mat <- feats_to_knn_matrix(feats, k = k)
    }
        
    meth_oe <- normalize_mat_knn(meth_mat, as.matrix(knn_mat))
    

    invisible(meth_oe)
}

#' Normalize methylation matrix for immune and CAF signatures
#' 
#' @inheritParams immune_cafs_signature_from_clust
#' @inheritParams feats_to_knn_matrix
#' 
#' @return normalized methylation matrix
#' 
#' @export
normalize_immune_cafs_meth <- function(em_clust, k, meth_mat = em_clust$meth_mat, scale_feats = TRUE, caf_gene = "CAV1", immune_gene = "CD3D") {
    feats <- immune_cafs_signature_from_clust(em_clust, scale_feats = scale_feats, caf_gene = caf_gene, immune_gene = immune_gene)
    knn_mat <- feats_to_knn_matrix(feats, k = k)    
    norm_meth <- calc_norm_immune_cafs_meth(meth_mat, k = k, feats = feats, knn_mat = knn_mat)
    
    return(norm_meth)
}


#' Plot scatter of raw expression-methylation correaltion vs normalized for specific genes
#' 
#' @param em_list_raw output of \code{em_cross_cor} on raw data
#' @param em_list_norm output of \code{em_cross_cor} on normalized data
#' @param gene name of the gene to plot
#' @param color points color
#' 
#' @return ggplot object
#' 
#' 
#' @export
plot_norm_vs_raw_gene_scatter <- function(em_list_raw, em_list_norm, gene, color="black"){
    df <- em_list_raw$em_cross[, gene] %>% enframe("locus", "raw_cor") %>% 
        left_join(em_list_norm$em_cross[, gene] %>% enframe("locus", "norm_cor"), by = "locus")
    limits <- c(min(c(df$raw_cor, df$norm_cor), na.rm=TRUE), max(c(df$raw_cor, df$norm_cor), na.rm=TRUE))
    p <- df %>% 
        ggplot(aes(x=raw_cor, y=norm_cor)) + 
            geom_point(size=0.01, color=color) + 
            geom_abline(linetype = "dashed", color="black") + 
            coord_cartesian(xlim = limits, ylim = limits) + 
            xlab("Raw correlation") + 
            ylab("Normzalied correlation") +
            ggtitle(gene) + 
            theme(aspect.ratio = 1)
    return(p)    
}