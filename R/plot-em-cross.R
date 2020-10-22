#' Plot clustered Expression-methylation cross-correlation matrix
#' 
#' @param em_clust output of \code{cluster_em_cross_cor}
#' @param fig_ofn filename of the output figure
#' @param width width of the saved figure
#' @param height height of the saved figure
#' @param k_meth k for \code{cutree} of matrix rows
#' @param k_expr k for \code{cutree} of matrix columns
#' 
#' @export
plot_em_cross_cor <- function(em_clust, fig_ofn=NULL, width = 1000, height = 1400, k_meth = max(meth_clust$clust), k_exp = max(expr_clust$clust), km_meth = NULL, km_exp = NULL, downscale = FALSE, zlim=c(-0.6, 0.6), colors=c("black", "darkred", "white", "darkblue", "cyan")) {
    list2env(em_clust, envir = environment())

    if (is.null(km_meth)) {
        km_meth <- cutree_order(hc_meth, k = k_meth)
    }
    if (is.null(km_exp)) {
        km_exp <- cutree_order(hc_exp, k = k_exp)
    }

    km_x <- tapply(seq(0, 1, l = length(hc_meth$order)), km_meth[hc_meth$order], mean, na.rm = TRUE)
    km_y <- tapply(seq(0, 1, l = length(hc_exp$order)), km_exp[hc_exp$order], mean, na.rm = TRUE)
    km_x_border <- tapply(seq(0, 1, l = length(hc_meth$order)), km_meth[hc_meth$order], max, na.rm = TRUE)
    km_y_border <- tapply(seq(0, 1, l = length(hc_exp$order)), km_exp[hc_exp$order], max, na.rm = TRUE)

    shades <- colorRampPalette(colors)(1000)        
    if (downscale){
        message("downscaling matrix")
        k_downscale <- min(round(nrow(em_cross) / height), round(ncol(em_cross) / width)) / 2        
        ds_mat <- downscale_matrix(em_cross[hc_meth$order, hc_exp$order], k = k_downscale, func = mean)
    } else {
        ds_mat <- em_cross[hc_meth$order, hc_exp$order]
    }
    
    if (!is.null(fig_ofn)){
        png(fig_ofn, width = width, height = height)    
    }

    message("plotting em cross")    
    image(ds_mat, zlim = zlim, col = shades, xaxt = "n", yaxt = "n")
    mtext(1:k_meth, at = km_x, side = 1, las = 2)
    mtext(1:k_exp, at = km_y, side = 2, las = 2)
    abline(v = km_x_border)
    abline(h = km_y_border)
    if (!is.null(fig_ofn)){
        dev.off()
    }
}

#' @export
plot_meth_mat_cm <- function(cm, k, fig_ofn=NULL, width = 1000, height = 1000, hc_meth=NULL, downscale = FALSE, zlim = c(-0.3, 0.3), colors = c("black", "darkred", "white", "darkblue", "cyan")) {    

    if (is.null(hc_meth)){
        hc_meth <- as.dist(1-cm) %>% hclust(method = "ward.D2")        
    }
        
    km_meth <- cutree_order(hc_meth, k = k)    

    km_x <- tapply(seq(0, 1, l = length(hc_meth$order)), km_meth[hc_meth$order], mean, na.rm = TRUE)
    km_x_border <- tapply(seq(0, 1, l = length(hc_meth$order)), km_meth[hc_meth$order], max, na.rm = TRUE)

    if (downscale){
        message("downscaling matrix")
        k_downscale <- round(min(round(nrow(cm) / height), round(ncol(cm) / width)) / 2)
        message("downscale k: ", k_downscale)
        ds_mat <- downscale_matrix(cm[hc_meth$order, hc_meth$order], k = k_downscale, func = mean)
    } else {
        ds_mat <- cm[hc_meth$order, hc_meth$order]
    }

    shades <- colorRampPalette(colors)(1000)

    if (!is.null(fig_ofn)){
        png(fig_ofn, width = width, height = height)    
    }
    
    message("plotting")
    image(ds_mat, zlim = zlim, col = shades, xaxt = "n", yaxt = "n")
    mtext(1:k, at = km_x, side = 1, las = 2)
    mtext(1:k, at = km_x, side = 2, las = 2)
    abline(v = km_x_border)
    abline(h = km_x_border)

    if (!is.null(fig_ofn)){
        dev.off()
    }
}