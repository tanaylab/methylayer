#' Extract TME genes from methylation expression clustering
#' 
#' Extracts names of genes that were clustered together with canonical immune and CAF genes in the expression-methylation correlation matrix. This is used mostly for filtering those genes from additional analysis such as cis-regulation. 
#' 
#' @param em_list output of \code{em_cross_cor}
#' @param immune_gene name of canonical immune gene (default: CD3D)
#' @param caf_gene name of canonical immune gene (default: CAV1)
#' 
#' @return vector with names of TME genes
#' 
#' @export
get_TME_genes <- function(em_list, immune_gene = "CD3D", caf_gene = "CAV1"){
	em_list$expr_clust %>%
        filter(grepl(glue("({caf_gene})|({immune_gene})"), name)) %>%
        select(clust) %>%
        inner_join(em_list$expr_clust, by = "clust") %>%
        pull(name)
}