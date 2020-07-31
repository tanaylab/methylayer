#' Transform an intervals set to matrix
#' 
#' @param df intervals set
#' 
#' @return matrix with coordinate rownames (chrom_start_end).
#' 
#' @export
intervs_to_mat <- function(df){    
    mat <- df %>% 
        unite("coord", chrom:end) %>% 
        as.data.frame() %>% 
        column_to_rownames("coord") %>% 
        as.matrix()

    return(mat)
}


#' Transform a matrix with coordinate rownames to intervals set 
#' 
#' @param mat matrix with coordinate rownames (chrom_start_end).
#' 
#' @return intervals set
#' 
#' @export
mat_to_intervs <- function(mat){
    df <- mat %>% 
        as.data.frame() %>% 
        rownames_to_column("coord") %>% 
        separate(coord, c("chrom", "start", "end")) %>% 
        mutate(start = as.numeric(start), end = as.numeric(end))

    return(df)
}

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
        left_join(promoter_intervs, by = c("chrom", "start", "end")) %>%         
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