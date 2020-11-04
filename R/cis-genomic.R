
#' Get candidates non-promoter cis regulation
#' 
#' @param meth_mat Matrix with methylation values. Each row is a locus and each column is a sample. Rownames should contain "chrom", "start" and "end" separated by "_"
#' @param expr_mat Matrix with expression values. Each row is a gene and each column is a sample.
#' @param gene_tss intervals set with additional "name" column with transcription start site (TSS) for each gene. Used for allocating a genomic coordinate for each gene. In case of duplicate names or coordinates the first one in the matrix will be used.
#' @param min_samples minimal number of samples per gene (both expression and methylation). Default is 100. 
#' @param max_k maximal number of genes to report per locus. Default is 50. Note that for efficiency (memory and computation) it is best to keep this number small. 
#' @param max_dist maximal distance from TSS to consider "cis" (used for FDR computations)
#' @param min_dist minimal distance from TSS (in order to exclude promoter regulation). Used for FDR computation. 
#' @param spearman use spearman correlation (if FALSE - use pearson)
#' @param parallel use mclapply for computation of ranks. Note that for large matrices forking might fail due to lack of memory and therefore parallel should be set to FALSE. 
#' 
#' 
#' @export
cis_em_genomic <- function(meth_mat, expr_mat, gene_tss, min_samples = 100, max_k = 50, max_dist = 5e5, min_dist = 200, spearman = FALSE, parallel = FALSE){
    samples <- intersect(colnames(meth_mat), colnames(expr_mat))   
    
    # calculate the number of samples with both expression and methylation per locus and gene
    message("calculating sample coverage...")
    cov_mat <- (1*!is.na(meth_mat[, samples])) %*% t(1*!is.na(expr_mat[, samples]))
    cov_mat_f <- cov_mat >= min_samples
    
    f_loci <- matrixStats::rowAnys(cov_mat_f)
    f_genes <- matrixStats::colAnys(cov_mat_f)
    message("# of samples: ", length(samples))
    message("# of genes: ", sum(f_genes))
    message("# of loci: ", sum(f_loci))
    
    message("calculating expression-methylation correlation...")
    m <- tgs_cor_large_matrices(t(meth_mat[f_loci, samples]), t(expr_mat[f_genes, samples]), spearman = spearman, pairwise.complete.obs=TRUE)
    
    message("ranking correlations of each locus (row) with all the genes (columns)...")
    top_k <- rank_genomic_cis_matrix(m, max_k, parallel = parallel)
    
    message("shuffling...")
    m_shuff <- shuffle_each_column(m)
    
    message("ranking correlations of each locus (row) with all the genes (columns)...")
    top_k_shuff <- rank_genomic_cis_matrix(m_shuff, max_k, parallel = parallel)
   
    cands <- bind_rows(
        top_k %>% mutate(type = "obs"),
        top_k_shuff %>% mutate(type = "shuff")
    )
    
    message("adding metadata")
    
    cands_df <- cands %>% separate(coord, c("chrom", "start", "end")) %>% mutate(start = as.numeric(start), end = as.numeric(end))  
    
    
    cands_df_genes <- cands_df %>% select(chrom:end, type, ends_with("name")) %>% gather("rank", "gene", -(chrom:end), -type) %>% mutate(rank = as.numeric(gsub("_name", "", rank)))
    
    cands_df_cors <- cands_df %>% select(chrom:end, type, ends_with("cor")) %>% gather("rank", "cor", -(chrom:end), -type) %>% mutate(rank = as.numeric(gsub("_cor", "", rank)))
    
    cands_df <- cands_df_genes %>% left_join(cands_df_cors, by = c("chrom", "start", "end", "type", "rank"))
    
    cands_df <- cands_df %>% left_join(gene_tss %>% distinct(name, .keep_all=TRUE) %>% select(gene = name, chrom_expr = chrom, start_expr = start, end_expr = end, strand_expr = strand), by = "gene")
    
    cands_df <- cands_df %>% mutate(dist = gintervals.distance(chrom1 = chrom, start1 = start, end1 = end, strand1 = 1, chrom2 = chrom_expr, start2 = start_expr, end2 = end_expr, strand2 = strand_expr))
    
    fdr_ranks <- compute_cis_fdr(cands_df, max_k, max_dist, min_dist)
    cands_df <- cands_df %>% left_join(fdr_ranks, by = "rank")
    
    return(cands_df)
}

rank_genomic_cis_matrix <- function(m, max_k, parallel=TRUE){
    if (parallel){
        buckets_df <- tibble(coord = rownames(m), bucket = ntile(1:nrow(m), tgutil:::num_physical_cores()))
        run_func <-   tgutil::run_wide      
    } else {
        buckets_df <- tibble(coord = rownames(m), bucket = 1)
        run_func <-   lapply      
    }

    l <- tgutil::run_wide(unique(buckets_df$bucket), function(x) {
        m_slice <- m[buckets_df$bucket == x, ]
        ranks <- matrixStats::rowRanks(m_slice, preserveShape = TRUE)
        rownames(ranks) <- rownames(m_slice)
        colnames(ranks) <- colnames(m_slice)
        
        gene_names <- apply(ranks, 1, function(x) names(sort(x)[1:max_k])) %>% t()
        colnames(gene_names) <- paste0(1:max_k, "_name")
       
        gene_cors <- lapply(1:nrow(ranks), function(i) m_slice[i, gene_names[i, ]] %>% rlang::set_names(paste0(1:max_k, "_cor"))) %>% do.call(rbind, .)
        rownames(gene_cors) <- rownames(m_slice)
        
        res <- cbind(as.data.frame(gene_names), as.data.frame(gene_cors))
        
        return(res)
    })
    
    cands_df <- do.call(rbind, l) %>% rownames_to_column("coord") %>% as_tibble()
    
    return(cands_df)
}

compute_cis_fdr <- function(cands_df, max_k, max_dist, min_dist){
    top_shuff <- cands_df %>% filter(type == "shuff", !is.na(dist), abs(dist) <= max_dist, abs(dist) > min_dist)
    top_obs <- cands_df %>% filter(type == "obs", !is.na(dist), abs(dist) <= max_dist, abs(dist) > min_dist)

    min_ranks <- top_obs %>%
        group_by(chrom, start, end) %>%
        summarise(r = min(rank)) %>%
        pull(r)    
    
    min_ranks_shuff <- top_shuff %>%
        group_by(chrom, start, end) %>%
        summarise(r = min(rank)) %>%
        pull(r)    

    fdr_ranks <- map_dfr(1:max_k, ~ tibble(rank = .x, n_obs = sum(min_ranks <= .x), n_shuff = sum(min_ranks_shuff <= .x))) %>% mutate(fdr = n_shuff / n_obs)
   
    return(fdr_ranks)
}

