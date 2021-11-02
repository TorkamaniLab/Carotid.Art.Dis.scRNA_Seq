## ---- Utility Code - Alternative Moncole 3 plot_genes_by_group ----
## Change plot_genes_by_group to plot_genes_by_group2
## Change code to allow for as_ordered for pass parameter ordering_type ----
## Convert/Save Notebook at standard r script file

## Slightly udpated plot_genes_by_group from Monocle 3
monocle_theme_opts <- function()
{
  theme(strip.background = element_rect(colour = 'white', fill = 'white')) +
    theme(panel.border = element_blank()) +
    theme(axis.line.x = element_line(size=0.25, color="black")) +
    theme(axis.line.y = element_line(size=0.25, color="black")) +
    theme(panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank()) +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_blank()) +
    theme(panel.background = element_rect(fill='white')) +
    theme(legend.key=element_blank())
}


#' Create a dot plot to visualize the mean gene expression and percentage of
#' expressed cells in each group of cells
#'
#' @param cds A cell_data_set for plotting.
#' @param markers A list of gene ids (or short names) to show in the plot
#' @param group_cells_by How to group cells when labeling them. Must be either
#'   the name of a column of colData(cds), or one of "clusters" or "partitions".
#'   If a column in colData(cds), must be a categorical variable.
#' @param reduction_method The dimensionality reduction method used for clusters
#'   and partitions.
#' @param norm_method Determines how to transform expression values prior to
#'   plotting. Options are "log" and "size_only". Default is "log".
#' @param lower_threshold The lowest gene expressed treated as expressed. By
#'   default, zero.
#' @param max.size The maximum size of the dot. By default, it is 10.
#' @param ordering_type How to order the genes / groups on the dot plot. Only
#'   accepts 'cluster_row_col' (use biclustering to cluster the rows and
#'   columns), 'maximal_on_diag' (position each column so that the maximal color
#'   shown on each column on the diagonal, if the current maximal is used in
#'   earlier columns, the next largest one is position), and 'none' (preserve
#'   the ordering from the input gene or alphabetical ordering of groups).
#'   Default is 'cluster_row_col'.
#' @param axis_order Whether to put groups on x-axis, genes on y-axis (option
#'   'group_marker') or the reverse order (option 'marker_group'). Default is
#'   "group_marker".
#' @param flip_percentage_mean Logical indicating whether to use color of the
#'   dot to represent the percentage (by setting flip_percentage_mean = FALSE,
#'   default) and size of the dot the mean expression, or the opposite (by
#'   setting flip_percentage_mean = TRUE).
#' @param pseudocount A pseudo-count added to the average gene expression.
#' @param scale_max The maximum value (in standard deviations) to show in the
#'   heatmap. Values larger than this are set to the max.
#' @param scale_min The minimum value (in standard deviations) to show in the
#'   heatmap. Values smaller than this are set to the min.
#'
#' @return a ggplot2 plot object
#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom reshape2 dcast
#' @importFrom viridis scale_color_viridis
#' @export
plot_genes_by_group2 <- function(cds,
                                markers,
                                group_cells_by="cluster",
                                reduction_method = "UMAP",
                                norm_method = c("log", "size_only"),
                                lower_threshold = 0,
                                max.size = 10,
                                # maybe be also do the maximum color on the
                                # diagonal; the axis change be switched too
                                ordering_type = c('cluster_row_col',
                                                  'maximal_on_diag',
                                                  'none'),
                                axis_order = c('group_marker', 'marker_group'),
                                flip_percentage_mean = FALSE,
                                pseudocount = 1,
                                scale_max = 3,
                                scale_min = -3) {

  assertthat::assert_that(methods::is(cds, "cell_data_set"))

  if(!is.null(group_cells_by)) {
    assertthat::assert_that(group_cells_by %in% c("cluster", "partition") |
                              group_cells_by %in% names(colData(cds)),
                            msg = paste("group_cells_by must be a column in",
                                        "the colData table."))
  }

  norm_method = match.arg(norm_method)

  gene_ids = as.data.frame(fData(cds)) %>%
    tibble::rownames_to_column() %>%
    dplyr::filter(rowname %in% markers | gene_short_name %in% markers) %>%
    dplyr::pull(rowname)
  if(length(gene_ids) < 1)
    stop(paste('Please make sure markers are included in the gene_short_name",
               "column of the fData!'))

  if(flip_percentage_mean == FALSE){
    major_axis <- 1
    minor_axis <- 2
  } else if (flip_percentage_mean == TRUE){
    major_axis <- 2
    minor_axis <- 1
  }

  exprs_mat <- t(as.matrix(exprs(cds)[gene_ids, ]))
  exprs_mat <- reshape2::melt(exprs_mat)
  colnames(exprs_mat) <- c('Cell', 'Gene', 'Expression')
  exprs_mat$Gene <- as.character(exprs_mat$Gene)


  if (group_cells_by == "cluster"){
    cell_group <- tryCatch({clusters(cds,
                                     reduction_method = reduction_method)},
                           error = function(e) {NULL})
  } else if (group_cells_by == "partition") {
    cell_group <- tryCatch({partitions(cds,
                                       reduction_method = reduction_method)},
                           error = function(e) {NULL})
  } else{
    cell_group <- colData(cds)[,group_cells_by]
  }

  if (length(unique(cell_group)) < 2) {
    stop(paste("Only one type in group_cells_by. To use plot_genes_by_group,",
               "please specify a group with more than one type. "))
  }

  names(cell_group) = colnames(cds)

  exprs_mat$Group <- cell_group[exprs_mat$Cell]
  exprs_mat = exprs_mat %>% dplyr::filter(is.na(Group) == FALSE)
  ExpVal <- exprs_mat %>% dplyr::group_by(Group, Gene) %>%
    dplyr::summarize(mean = mean(log(Expression + pseudocount)),
                     percentage = sum(Expression > lower_threshold) /
                       length(Expression))
  ExpVal$mean <- ifelse(ExpVal$mean < scale_min, scale_min, ExpVal$mean)
  ExpVal$mean <- ifelse(ExpVal$mean > scale_max, scale_max, ExpVal$mean)

  ExpVal$Gene <- fData(cds)[ExpVal$Gene, 'gene_short_name']

  res <- reshape2::dcast(ExpVal[, 1:4], Group ~ Gene,
                         value.var = colnames(ExpVal)[2 + major_axis])
  group_id <- res[, 1]
  res <- res[, -1]
  row.names(res) <- group_id

  if(ordering_type == 'cluster_row_col') {
    row_dist <- stats::as.dist((1 - stats::cor(t(res)))/2)
    row_dist[is.na(row_dist)] <- 1

    col_dist <- stats::as.dist((1 - stats::cor(res))/2)
    col_dist[is.na(col_dist)] <- 1

    ph <- pheatmap::pheatmap(res,
                             useRaster = T,
                             cluster_cols=TRUE,
                             cluster_rows=TRUE,
                             show_rownames=F,
                             show_colnames=F,
                             clustering_distance_cols=col_dist,
                             clustering_distance_rows=row_dist,
                             clustering_method = 'ward.D2',
                             silent=TRUE,
                             filename=NA)

    ExpVal$Gene <- factor(ExpVal$Gene,
                          levels = colnames(res)[ph$tree_col$order])
    ExpVal$Group <- factor(ExpVal$Group,
                           levels = row.names(res)[ph$tree_row$order])

  } else if(ordering_type == 'maximal_on_diag'){

    order_mat <- t(apply(res, major_axis, order))
    max_ind_vec <- c()
    for(i in 1:nrow(order_mat)) {
      tmp <- max(which(!(order_mat[i, ] %in% max_ind_vec)))
      max_ind_vec <- c(max_ind_vec, order_mat[i, tmp])
    }
    max_ind_vec <- max_ind_vec[!is.na(max_ind_vec)]

    if(major_axis == 1){
      max_ind_vec <- c(max_ind_vec, setdiff(1:length(markers), max_ind_vec))
      ExpVal$Gene <- factor(ExpVal$Gene ,
                            levels = dimnames(res)[[2]][max_ind_vec])
    }
    else{
      max_ind_vec <- c(max_ind_vec, setdiff(1:length(unique(exprs_mat$Group)),
                                            max_ind_vec))
      ExpVal$Group <- factor(ExpVal$Group,
                             levels = dimnames(res)[[1]][max_ind_vec])
    }
######### This is a block added by D. Evans 8.27.2019 (see comments at the end of this block)
 } else if(ordering_type == 'as_ordered'){

    order_mat <- t(apply(res, major_axis, order))
    max_ind_vec <- c()
    for(i in 1:nrow(order_mat)) {
      tmp <- max(which(!(order_mat[i, ] %in% max_ind_vec)))
      max_ind_vec <- c(max_ind_vec, order_mat[i, tmp])
    }
    max_ind_vec <- max_ind_vec[!is.na(max_ind_vec)]

    if(major_axis == 1){
      max_ind_vec <- c(max_ind_vec, setdiff(1:length(markers), max_ind_vec))
      ExpVal$Gene <- factor(ExpVal$Gene ,
                            levels = dimnames(res)[[2]][max_ind_vec])
      GlobExpVal <<- ExpVal
    }
    else{
      max_ind_vec <- c(max_ind_vec, setdiff(1:length(unique(exprs_mat$Group)),
                                            max_ind_vec))
      ExpVal$Group <- factor(ExpVal$Group,
                             levels = dimnames(res)[[1]][max_ind_vec])
      GlobExpVal <<- ExpVal
    }
  
    ## This is the code that was added to ordering_type == 'maximal_on_diag'
    ## to make it work with ordering_type == 'as_ordered' pass parameter 
    ## This re-orders the factors according to the marker HUGO names (names of markers vector)   
    ExpVal$Gene <- factor(ExpVal$Gene, levels = names(markers))
      
  } else if(ordering_type == 'none'){
     ExpVal$Group <- factor(ExpVal$Gene, levels = markers)
    
  }    
    
  if(flip_percentage_mean){
    g <- ggplot(ExpVal, aes(y = Gene,  x = Group)) +
      geom_point(aes(colour = percentage,  size = mean)) +
      viridis::scale_color_viridis(name = 'percentage') +
      scale_size(name = 'log(mean + 0.1)', range = c(0, max.size))
  } else {
      g <- ggplot(ExpVal, aes(y = Gene,  x = Group)) +
      geom_point(aes(colour = mean,  size = percentage)) +
      viridis::scale_color_viridis(name = 'log(mean + 0.1)') +
      scale_size(name = 'percentage', range = c(0, max.size))
  }
      
      
  if (group_cells_by == "cluster"){
    g <- g + xlab("Cluster")
  } else if (group_cells_by == "partition") {
    g <- g + xlab("Partition")
  } else{
    g <- g + xlab(group_cells_by)
  }

  g <- g + ylab("Gene") + monocle_theme_opts() +
   theme(axis.text.x = element_text(angle = 30, hjust = 1))
  if(axis_order == 'marker_group') {
    g <- g + coord_flip()
  }

  g
}
