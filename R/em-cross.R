#' Calculate Expression-methylation cross-correlation
#' 
#' @param meth_mat Matrix with methylation values. Each row is a locus and each column is a sample.
#' @param expr_mat Matrix with expression values. Each row is a gene and each column is a sample.
#' @param min_meth minimal locus average methylation
#' @param max_meth maximal locus average methylation
#' @param min_expr minimal expression level. Only genes with at least one sample with expression above min_expr would be included.
#' @param min_sd minimal standard deviation per locus. 
#' @param samples names of samples to include. If NULL - all samples that have both expression and methylation data would be included.
#' @param max_na maximal number of NAs allowed per locus.
#' @param meth_cor_thresh minimal correlation level of a locus. Only loci with at least one correlation above the threshold will be included. 
#' @param expr_cor_thresh minimal correlation level of a gene. Only genes with at least one correlation above the threshold will be included. 
#' 
#' @return a list with:
#' \itemize{
#'  \item{em_cross}{Matrix with expression-methylation correlation values. Rows are loci and columns are genes}
#'  \item{meth_mat}{Matrix with methylation values (filtered)}
#'  \item{expr_mat}{Matrix with expression values (filtered)}
#' }
#' 
#' @examples
#' 
#' \dontrun{
#'  init_tcga_brca_example()
#'  em_list <- get_em_cross(brca_meth_mat, brca_expr_mat, min_meth = 0.1, max_meth = 0.9, meth_cor_thresh = 0.25, expr_cor_thresh = 0.25, min_expr = 2, min_sd = 0.05, max_na = 0)
#'   
#' 
#' }
#' 
#' 
#' @export
em_cross_cor <- function(meth_mat, expr_mat, min_meth = NULL, max_meth = NULL, min_expr = NULL, min_sd = NULL, samples = NULL, max_na = NULL, meth_cor_thresh = NULL, expr_cor_thresh = NULL){
    
    samples <- samples %||% intersect(colnames(meth_mat), colnames(expr_mat))    

    message("# of samples: ", length(samples))
      
    if (!is.null(min_meth) && !is.null(max_meth)) {
        meth_mat <- filter_meth_mat_by_avg(meth_mat, min_meth, max_meth)
    }    
    
    if (!is.null(min_sd)) {
        meth_locus_sd <- matrixStats::rowSds(as.matrix(meth_mat), na.rm = TRUE)
        f_sd_meth <- meth_locus_sd > min_sd
        meth_mat <- meth_mat[f_sd_meth, ]
    }
    
    if (!is.null(min_expr)){
        gene_maxs <- rowMaxs(expr_mat, na.rm=TRUE)
        expr_mat <- expr_mat[gene_maxs >= min_expr, ]
    }
    
    message("expression (columns): ", nrow(expr_mat))
    message("methylation (rows): ", nrow(meth_mat))
    
    em_cross <- tgs_cor(t(meth_mat[, samples]), t(expr_mat[, samples]), spearman=TRUE, pairwise.complete.obs=TRUE)
    
    if (!is.null(meth_cor_thresh)){
        meth_maxcor_f <- apply(abs(em_cross), 1, max, na.rm = TRUE) > meth_cor_thresh    
        message(glue("{rows} rows had at least one cor > {meth_cor_thresh}", rows = sum(meth_maxcor_f)))
    } else {
        meth_maxcor_f <- rep(TRUE, nrow(em_cross))
    }
    
    if (!is.null(expr_cor_thresh)){
        expr_maxcor_f <- apply(abs(em_cross), 2, max, na.rm = TRUE) > expr_cor_thresh    
        message(glue("{columns} columns had at least one cor > {expr_cor_thresh}", columns = sum(expr_maxcor_f)))
    } else {
        expr_maxcor_f <- rep(TRUE, ncol(em_cross))
    }    

    if (!is.null(max_na)){
        meth_maxcor_f <- meth_maxcor_f & rowSums(is.na(em_cross)) <= max_na
        message(glue("{rows} rows did not have more than {max_na} Na's ", rows = sum(meth_maxcor_f)))
    }
    
    res <- list(
        em_cross = em_cross[meth_maxcor_f, expr_maxcor_f],
        meth_mat = meth_mat[meth_maxcor_f, samples],        
        expr_mat = expr_mat[expr_maxcor_f, samples]
    )
    
    return(res)        
}



#' Cluster Expression-methylation cross-correlation
#' 
#' @param em_list output of \code{em_cross_cor}
#' @param k_meth number of methylation clusters
#' @param k_expr number of expression clusters
#' @param hc_exp \code{hclust} object with existing clustering of the matrix columns
#' @param normal_meth matrix with methylation in normal samples
#' 
#' @return a list with:
#' \itemize{
#'  \item{"em_cross"}{Matrix with expression-methylation correlation values. Rows are loci and columns are genes}
#'  \item{meth_mat}{Matrix with methylation values (filtered)}
#'  \item{expr_mat}{Matrix with expression values (filtered)}
#'  \item{hc_meth}{hclust object with clustering of the loci}
#'  \item{hc_exp}{hclust object with clustering of the genes}
#'  \item{meth_clust}{tibble with cluster number for each row}
#'  \item{expr_clust}{tibble with cluster number for each column}
#'  \item{expr_mods}{matrix with mean expression in every cluster for each sample}
#'  \item{meth_mods}{matrix with mean methylation in every cluster for each sample}
#'  \item{expr_mods_locus}{matrix with mean expression in every cluster for each locus}
#'  \item{meth_mods_gene}{matrix with mean methylation in every cluster for each gene}
#' }
#' 
#' @export
cluster_em_cross_cor <- function(em_list, k_meth, k_expr=k_meth, hc_exp = NULL, normal_meth = NULL) {
    list2env(em_list, envir = environment())

    hc_meth <- hclust(tgs_dist(em_cross), "ward.D2")
    if (is.null(hc_exp)) {
        hc_exp <- hclust(tgs_dist(t(em_cross)), "ward.D2")
    }

    km_meth <- cutree_order(hc_meth, k = k_meth)
    meth_clust <- enframe(km_meth, name = "locus", value = "clust")

    km_expr <- cutree_order(hc_exp, k = k_expr)
    expr_clust <- enframe(km_expr, name = "name", value = "clust")

    if (!is.null(normal_meth)) {
        sort_mat <- normal_meth        
    } else {
        sort_mat <- meth_mat
    }
    
    expr_mods <- t(tgs_matrix_tapply(t(expr_mat[names(km_expr), ]), km_expr, mean, na.rm=TRUE))
    rownames(expr_mods) <- colnames(expr_mat)
    meth_mods <- t(tgs_matrix_tapply(t(meth_mat[names(km_meth), ]), km_meth, mean, na.rm = TRUE))
    rownames(meth_mods) <- colnames(meth_mat)

    meth_mods_gene <- t(tgs_matrix_tapply(t(em_cross), meth_clust$clust, mean, na.rm=TRUE))
    rownames(meth_mods_gene) <- colnames(em_cross)

    expr_mods_locus <- t(tgs_matrix_tapply(em_cross, expr_clust$clust, mean, na.rm=TRUE))
    rownames(expr_mods_locus) <- rownames(em_cross)

    res <- list(
        em_cross = em_cross,
        meth_mat = meth_mat,
        expr_mat = expr_mat,
        hc_meth = hc_meth,
        hc_exp = hc_exp,
        meth_clust = meth_clust,
        expr_clust = expr_clust,
        expr_mods = expr_mods,
        meth_mods = meth_mods,
        expr_mods_locus = expr_mods_locus,
        meth_mods_gene = meth_mods_gene
    )

    return(res)
}

#' Export Expression-methylation cross-correlation clustering to excel sheet
#' 
#' @param em_clust output of \code{cluster_em_cross_cor}
#' @param filename name of excel file name
#' @param promoters intervals set with additional column 'name' with promoter names
#' 
#' @export
export_em_cross_mods <- function(em_clust, filename, promoters=NULL){
    expr_sheet <- em_clust$meth_mods_gene %>% as.data.frame() %>% rownames_to_column("gene") %>% as.data.frame()
    expr_sheet <- expr_sheet %>% left_join(em_clust$expr_clust %>% rename(gene = name), by = "gene") %>% select(gene, clust, everything()) %>% arrange(clust, gene)
    colnames(expr_sheet)[-1:-2] <- paste0("METH_", colnames(expr_sheet)[-1:-2])

    meth_sheet <- em_clust$expr_mods_locus %>% mat_to_intervs()    
    meth_sheet <- meth_sheet %>% left_join(em_clust$meth_clust %>% column_to_rownames("locus") %>% mat_to_intervs(), by = c("chrom", "start", "end")) %>% select(chrom:end, clust, everything()) %>% arrange(clust, chrom, start, end)
    colnames(meth_sheet)[-1:-4] <- paste0("EXPR_", colnames(meth_sheet)[-1:-4])

    if (!is.null(promoters)){
        meth_sheet <- meth_sheet %>% left_join(promoters %>% select(chrom, start, end, name), by=c("chrom", "start", "end")) %>% select(chrom, start, end, name, clust, everything())
    }    

    writexl::write_xlsx(x = list("Loci (Expression)" = meth_sheet, "Genes (Methylation)" = expr_sheet), filename)
}



