
#' Deconvolute TME (tumor microenvironment) effects from methylation data
#' 
#' @param meth_mat Matrix with methylation values to use for TME inference (usually - promoter methylation). Each row is a locus and each column is a sample. This can be a different matrix than \code{raw_meth_mat}.
#' @param expr_mat Matrix with expression values. Each row is a gene and each column is a sample.
#' @param raw_meth_mat Matrix with methylation values to deconvolute. Each row is a locus and each column is a sample. This can be a different matrix than \code{meth_mat}.
#' @param k_meth number of methylation clusters when clustering expression-methylation data
#' @param k_expr number of expression clusters when clustering expression-methylation data
#' @param caf_gene name of gene that is inside the CAF cluster (default: "CAV1")
#' @param immune_gene name of gene that is inside the immune cluster (default: "CD3D")
#' @param scores_cor_thresh minimal (anti) correlation of loci methylation to TME expression scores in order to create TME methylatoon scores
#' @param ... additional parameters to \code{em_cross_cor}
#' 
#' @return a list with:
#' \itemize{
#'  \item{norm_meth}{Matrix with normalized methylation values of \code{raw_meth_mat}}
#'  \item{tme_features}{Data frame with Immune and CAF expression and methylation signatures per sample}
#'  \item(tme_loci){list with intervals set of loci that were included in the immune and CAF methylation scores}
#'  \item{em_cross}{List as returned from \code{em_cross_cor}. Mainly includes a matrix with expression-methylation correlation values that were used for the deconvolution algorithm, where rows are loci and columns are genes}
#'  \item{em_cross_clust}{Clustering of expression-methylation correlation data, as returned from \code{cluster_em_cross_cor}}
#' }
#' 
#' @inheritParams em_cross_cor
#' @inheritParams normalize_immune_cafs_meth
#' @export
deconv_TME <- function(meth_mat, expr_mat, raw_meth_mat, min_meth = 0.1, max_meth = 0.9, min_expr = NULL, meth_cor_thresh = 0.25, expr_cor_thresh = 0.25, min_sd = NULL, k = NULL, k_meth = 30, k_expr = k_meth, caf_gene = "CAV1", immune_gene = "CD3D", scores_cor_thresh = 0.3, ...){
    
    if (!is.null(min_meth) && !is.null(max_meth)) {
        meth_mat <- filter_meth_mat_by_avg(meth_mat, min_meth, max_meth)
    }
    
    if (is.null(min_sd)){
        meth_locus_sd <- matrixStats::rowSds(as.matrix(meth_mat), na.rm = TRUE)
        min_sd <- quantile(meth_locus_sd, 0.1, na.rm=TRUE)
        message(glue("min_sd: {min_sd}"))
    }
    
    if (is.null(min_expr)){
        gene_maxs <- matrixStats::rowMaxs(expr_mat, na.rm=TRUE)
        min_expr <- quantile(gene_maxs[is.finite(gene_maxs)], 0.05, na.rm=TRUE)
        expr_mat <- expr_mat[is.finite(gene_maxs) & !is.na(gene_maxs) & gene_maxs >= min_expr, ]
        message(glue("min_expr: {min_expr}"))
    }    
   
    message("calculating em-cross")
    em_cross <- em_cross_cor(meth_mat, expr_mat, min_meth = min_meth, max_meth = max_meth, max_na = 0, meth_cor_thresh = meth_cor_thresh, expr_cor_thresh = expr_cor_thresh, ...)    
    
    message("clustering em-cross")
    em_cross_clust <- cluster_em_cross_cor(em_cross, k_meth = k_meth, k_expr = k_expr)
    
    if (is.null(k)){
        k <- round(min(15, 0.03 * ncol(meth_mat)))
        message(glue("k: {k}"))
    }
    
    message("normalizing methylation")
    feats <- immune_cafs_signature_from_clust(em_cross_clust, scale_feats = TRUE, caf_gene = caf_gene, immune_gene = immune_gene)
    knn_mat <- feats_to_knn_matrix(feats, k = k)    
    norm_meth <- calc_norm_immune_cafs_meth(raw_meth_mat, k = k, feats = feats, knn_mat = knn_mat)
    
    
    message("calculating TME methylation scores")
    tme_meth_scores <- get_TME_meth_scores(em_cross$meth_mat, feats, scores_cor_thresh)
    
    feats <- feats %>% 
        as.data.frame() %>% 
        rownames_to_column("samp") %>% 
        left_join(tme_meth_scores$immune_meth_score, by = "samp") %>% 
        left_join(tme_meth_scores$caf_meth_score, by = "samp")
    
    res <- list(
        norm_meth = norm_meth,
        tme_features = feats,
        tme_loci = list(immune_loci = tme_meth_scores$immune_loci, caf_loci = tme_meth_scores$caf_loci),
        em_cross = em_cross,
        em_cross_clust = em_cross_clust        
    )
    
    return(res) 
}


get_TME_meth_scores <- function(meth_mat, feats, scores_cor_thresh){
    loci_cors <- tgs_cor(t(meth_mat), feats, pairwise.complete.obs=TRUE)
    immune_loci <- loci_cors %>% mat_to_intervs() %>% filter(immune <= -scores_cor_thresh) %>% select(chrom:end)
    caf_loci <- loci_cors %>% mat_to_intervs() %>% filter(caf <= -scores_cor_thresh) %>% select(chrom:end)
    immune_meth_score <- meth_mat[immune_loci %>% intervs_to_mat() %>% rownames(), ] %>% colMeans(na.rm=TRUE) %>% enframe("samp", "immune.meth")
    caf_meth_score <- meth_mat[caf_loci %>% intervs_to_mat() %>% rownames(), ] %>% colMeans(na.rm=TRUE) %>% enframe("samp", "caf.meth")
    
    return(list(
        immune_loci = immune_loci,
        caf_loci = caf_loci,
        immune_meth_score = immune_meth_score,
        caf_meth_score = caf_meth_score        
    ))    
}




#' calculate expression of loci to specific genes before and after TME deconvolution
#' 
#' @param deconv_list output of \code{deconv_TME}.
#' @param raw_meth intervals set with raw methylation
#' @param genes vector with gene names
#' 
#' @return data frame with intervals set that are in raw_meth and normalized methylation and their correlations to the gene expression in the raw and normalized methylation datasets. 
#' 
#' 
#' @export
calc_gene_cor_before_after_deconv <- function(deconv_list, raw_meth, genes){
    expr_mat <- deconv_list$em_cross$expr_mat[genes, ]
    raw_meth_mat <- raw_meth %>% intervs_to_mat()
    samples <- intersect(colnames(deconv_list$norm_meth), colnames(raw_meth_mat)) %>% intersect(colnames(expr_mat))

    cm_raw <- tgs_cor(t(raw_meth_mat[, samples]), t(expr_mat[, samples]), pairwise.complete.obs=TRUE)
    cm_norm <- tgs_cor(t(deconv_list$norm_meth[, samples]), t(expr_mat[, samples]), pairwise.complete.obs=TRUE)
    colnames(cm_norm) <- paste0(colnames(cm_norm), ".norm")

    res <- cm_raw %>% mat_to_intervs() %>% inner_join(cm_norm %>% mat_to_intervs(), by = c("chrom", "start", "end"))

    return(res)
}


