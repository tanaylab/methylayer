#' Change rownames of matrix from coordinates to promoter names
#' 
#' @param mat matrix with coordinate rownames (chrom_start_end).
#' @param promoter_intervs intervals set with additional "name" column. In case of duplicate names or coordinates the first one in the matrix will be used. 
#' 
#' @return \code{mat} with rownames changed to promoter name taken from \code{promtoer_intervs}
#' 
#' @export
coord_to_promoter_names <- function(mat, promoter_intervs){
    promoter_intervs <- promoter_intervs %>%                 
        distinct(chrom, start, end, name)

    mat <- mat %>% 
        mat_to_intervs() %>% 
        inner_join(promoter_intervs, by = c("chrom", "start", "end")) %>%         
        distinct(name, .keep_all=TRUE) %>% 
        distinct(chrom, start, end, .keep_all=TRUE) %>% 
        select(-(chrom:end)) %>% 
        column_to_rownames("name") %>% 
        as.matrix()
    return(mat)    
}

promoter_names_to_coord <- function(mat, promoter_intervs){
    promoter_intervs <- promoter_intervs %>%                 
        distinct(chrom, start, end, name)
    
    mat <- mat %>% 
        rownames_to_column("name") %>% 
        as.data.frame() %>%        
        left_join(promoter_intervs, by = c("chrom", "start", "end")) %>%         
        distinct(name, .keep_all=TRUE) %>% 
        distinct(chrom, start, end, .keep_all=TRUE) %>% 
        select(-name) %>% 
        select(chrom, start, end, everything()) %>% 
        intervs_to_mat()
        
    return(mat)    
}

filter_meth_mat_by_avg <- function(meth_mat, min_meth, max_meth){
    meth_locus_avgs <- rowMeans(meth_mat, na.rm = TRUE)
    f_mid_meth <- !is.na(meth_locus_avgs) & meth_locus_avgs > min_meth & meth_locus_avgs < max_meth
    meth_mat <- meth_mat[f_mid_meth, ]
    return(meth_mat)
}

gintervals.distance <- function(chrom1, start1, end1, strand1, chrom2, start2, end2, strand2) {
    left_dist <- ifelse(strand2 == 1, 
        -1 * (start2 - end1), 
        start2 - end1
    )
    
    right_dist <- ifelse(strand2 == 1, 
        -1 * (end2 - start1), 
        end2 - start1
    )
    
    d <- ifelse(
        abs(left_dist) <= abs(right_dist), 
        left_dist, 
        right_dist
    )
    
    d <- ifelse(
        pmax(start1, start2) < pmin(end1, end2),
        0,
        d        
    )
    
    d <- ifelse(chrom1 == chrom2, 
        d,         
        NA)
    
    return(d)
}

########################################################################
# Split the matrix x to k-by-k blocks and apply func to each block
downscale_matrix <- function(x, k, func = mean) {
    row_new <- ceiling(nrow(x) / k)
    col_new <- ceiling(ncol(x) / k)
    grid <- get_downscale_grid(x, k)

    calc_block <- function(i) {
        block <- x[grid$rstart[i]:grid$rend[i], grid$cstart[i]:grid$cend[i]]
        return(func(block))
    }
    result <- run_wide(1:nrow(grid), calc_block)
    result <- unlist(result)
    dim(result) <- c(row_new, col_new)

    return(result)
}

get_downscale_grid <- function(x, k) {
    row_new <- ceiling(nrow(x) / k)
    col_new <- ceiling(ncol(x) / k)
    grid <- tidyr::crossing(cstart = (1:col_new - 1) * k + 1, rstart = (1:row_new - 1) * k + 1) %>%
        arrange(cstart, rstart)
    grid <- grid %>%
        mutate(cend = pmin(cstart + k - 1, ncol(x)), rend = pmin(rstart + k - 1, nrow(x)))
    return(grid)
}

get_downscale_block <- function(x, k, type = "row") {
    calc_row_block <- function(i) {
        data.frame(row = rownames(x)[grid$rstart[i]:grid$rend[i]], ds_row = i)
    }

    calc_column_block <- function(i) {
        data.frame(col = colnames(x)[grid$cstart[i]:grid$cend[i]], ds_col = i)
    }

    if (type == "row") {
        if (is.null(rownames(x))) {
            rownames(x) <- 1:nrow(x)
        }
        grid <- get_downscale_grid(x, k) %>% distinct(rstart, rend)
        result <- run_wide(1:nrow(grid), calc_row_block)
        result <- map_dfr(result, ~.x)
    }

    if (type == "column") {
        if (is.null(colnames(x))) {
            colnames(x) <- 1:ncol(x)
        }
        grid <- get_downscale_grid(x, k) %>% distinct(cstart, cend)
        result <- run_wide(1:nrow(grid), calc_column_block)
        result <- map_dfr(result, ~.x)
    }

    return(result)
}


########################################################################
tgs_cor_large_matrices <- function(x, y, pairwise.complete.obs = FALSE, spearman = FALSE, nbins = 10){
    bucket <- ntile(rownames(x), nbins)
    res <- map(unique(bucket), ~ {        
        tgs_cor(x[, bucket == .x], y, pairwise.complete.obs = pairwise.complete.obs, spearman = spearman)
    } )

    res <- do.call(rbind, res)
    return(res)    
}