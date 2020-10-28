
#' Get candidates for promoters that are regulated by methylation in cis
#' 
#' @param meth_mat Matrix with methylation values. Each row is a promoter and each column is a sample. Rownames should contain "chrom", "start" and "end" separated by "_"
#' @param expr_mat Matrix with expression values. Each row is a gene and each column is a sample.
#' @param promoter_intervs intervals set with additional "name" column. In case of duplicate names or coordinates the first one in the matrix will be used. Only intervals that are within promoter_intervs would be used. 
#' @param k_locus_rank rank of the correlation to extract for each promoter. For example, if \code{k_locus_rank = 2} the second best correlation would be extracted for every promoter in the field \code{kth}
#' @param min_samples minimal number of samples per gene (both expression and methylation). Default is 100. 
#' @param spearman use spearman correlation (if FALSE - use pearson)
#' 
#' @return dataframe with the following columns:
#' \itemize{
#'  \item{name}{name of the gene}
#'  \item{rank}{rank of the correlation of the gene with it's promoter}
#'  \item{cor}{correlation of the gene with it's promoter}
#'  \item{kth}{\code{k_locus_rank}th best correlation of the gene with a promoter}
#'  \item{fdr}{false discovery rate for detecting the cis-raget with rank value <= \code{rank}. Assuming independence between expression and methylation, the false discovery rate (FDR) for detecting the cis-target (diagonal value) with rank value <= k was defined as k / m, where m is the number of observed genes having cis-target (diagonal value) with rank value <= k.}
#' }
#' 
#' 
#' @export
cis_em_promoters <- function(meth_mat, expr_mat, promoter_intervs, k_locus_rank = 2, min_samples=100, spearman = FALSE) {

    meth_mat <- coord_to_promoter_names(meth_mat, promoter_intervs)
    
    samples <- intersect(colnames(meth_mat), colnames(expr_mat))    
    genes <- intersect(rownames(meth_mat), rownames(expr_mat))

    f <- rowSums(!is.na(meth_mat[genes, ])) >= min_samples & rowSums(!is.na(expr_mat[genes, ])) >= min_samples
    genes <- genes[f]
    
    message("# of samples: ", length(samples))
    message("# of genes: ", length(genes))
    message("calculating expression-methylation correlation...")
    m <- tgs_cor(t(meth_mat[genes, samples]), t(expr_mat[genes, samples]), spearman=spearman, pairwise.complete.obs=TRUE)

    genes <- intersect(colnames(m), rownames(m))
    m <- m[genes, genes]

    message("ranking...")
    m_rank <- matrixStats::colRanks(m, preserveShape = TRUE)
    colnames(m_rank) <- colnames(m)
    rownames(m_rank) <- rownames(m)

    cands <- enframe(diag(m), value = "cor") %>%
        left_join(enframe(diag(m_rank), value = "r"), by = "name") 

    cands <- cands %>%
        mutate(kth = apply(m, 2, function(x) x[order(x)][k_locus_rank])) %>%
        mutate(best = apply(m, 2, function(x) x[order(x)][1])) %>% 
        filter(!is.na(r))


    message("calculaing FDR...")
    min_ranks <- cands %>%
        group_by(name) %>%
        summarise(r = min(r)) %>%
        pull(r)
        
    fdr_ranks <- tibble(r = min_ranks) %>% 
        count(r, name="n_r") %>% 
        arrange(r) %>% 
        mutate(m = cumsum(n_r), fdr = pmin(1, r / m))
        
    cands <- cands %>% left_join(fdr_ranks %>% select(r, fdr, n_fdr=m), by = "r")

    return(cands)
}
