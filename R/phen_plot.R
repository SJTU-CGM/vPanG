


#' phen_heatmap
#'
#' Plot the result of phenotypic association in a heatmap.
#'
#' @param pav_obj A PAV object.
#' @param phen_stat_res The result from \code{\link[vPan]{phen_stat}}.
#' @param add_gene_info A character vector of `gene_info` names.
#'
#' @param p_threshold The threshold of p_value/p_adjusted.
#' @param adjust_p A logical value indicating whether adjust p_value.
#' @param only_show_significant A logical value indicating whether only show p_value/p_adjusted that satisfies the condition.
#' @param flip A logical value indicating whether flip cartesian coordinates.
#'
#' @param p_colors A vector of colors or a color mapping function for p_value/p_adjusted, pass to \code{\link[ComplexHeatmap]{Heatmap}}.
#' @param na_col A string of color for NA values.
#' @param cell_border_color `NA` or a string of color for the border of cells.
#' @param gene_info_color_list A list contains named vector of colors for `gene_info` annotation.
#' e.g. list(source = c("reference" = "red", "novel" = "blue"), length = c("orange", "red"))
#'
#' @param cluster_rows A logical value indicating whether perform clustering on rows.
#' @param clustering_distance_rows Method of measuring distance when clustring on rows, pass to \code{\link[stats]{dist}}.
#' @param clustering_method_rows Method to perform hierarchical clustering on rows, pass to \code{\link[stats]{hclust}}.
#' @param row_dend_side The position of the row dendrogram ("left", "right").
#' @param row_dend_width A \code{\link[grid]{unit}} object for the width of the row dendrogram.
#' @param row_sorted A vector of sorted row names. It only works when `cluster_rows = F`.
#'
#' @param show_row_names  A logical value indicating whether show row names.
#' @param row_names_side The position of row names ("left", "right").
#' @param row_names_size The size of row names.
#' @param row_names_rot The rotation of row names.
#'
#' @param cluster_columns A logical value indicating whether perform clustering on columns.
#' @param clustering_distance_columns Method of measuring distance when clustring on columns, pass to \code{\link[stats]{dist}}.
#' @param clustering_method_columns Method to perform hierarchical clustering on columns, pass to \code{\link[stats]{hclust}}.
#' @param column_dend_side The position of the column dendrogram ("top", "bottom").
#' @param column_dend_height A \code{\link[grid]{unit}} object for the height of the column dendrogram.
#' @param column_sorted A vector of sorted column names. It only works when `cluster_columns = F` and `split_block = F`.
#'
#' @param show_column_names A logical value indicating whether show column names.
#' @param column_names_side The position of column names ("top", "column").
#' @param column_names_size The size of column names.
#' @param column_names_rot The rotation of column names.
#'
#' @param anno_param_gene A list contains parameters for the phenotype annotation. These can be any of the following:
#' "show", "width", "border", "name_size", "name_rot" and "name_side".
#'
#' @param legend_side The position of legend ("top", "bottom", "right", "left").
#' @param legend_title The text for the legend title.
#' @param legend_title_size The size of legend title.
#' @param legend_text_size The size of legend item labels.
#' @param legend_grid_size A \code{\link[grid]{unit}} object for the size of legend grid.
#'
#'
#' @export

phen_heatmap <- function(pav_obj,
                         phen_stat_res,
                         add_gene_info = NULL,

                         p_threshold = 0.01,
                         adjust_p = T,
                         only_show_significant = T,
                         flip = F,

                         p_colors = c("#1874CD", "#B0E2FF"),
                         na_col = "gray",
                         cell_border_color = NA,
                         gene_info_color_list = NULL,

                         cluster_rows = F,
                         clustering_distance_rows = "euclidean",
                         clustering_method_rows = "complete",
                         row_dend_side = "left",
                         row_dend_width = grid::unit(5, "mm"),
                         row_sorted = c(),

                         show_row_names = F,
                         row_names_side = "right",
                         row_names_size = 10,
                         row_names_rot = 0,

                         cluster_columns = F,
                         clustering_distance_columns = "euclidean",
                         clustering_method_columns = "complete",
                         column_dend_side = "top",
                         column_dend_height = grid::unit(5, "mm"),
                         column_sorted = c(),

                         show_column_names = F,
                         column_names_side = "bottom",
                         column_names_size = 10,
                         column_names_rot = 90,

                         anno_param_gene = list(show = T, width = 5, border = FALSE,
                                                name_size = NULL, name_rot = 90, name_side = "bottom"),

                         legend_side = "right",
                         legend_title = list(pav = "PAV", type = "gene"),
                         legend_title_size = NULL,
                         legend_text_size = NULL,
                         legend_grid_size = grid::unit(4, "mm")){


  check_class(pav_obj, "PAV")

  check_phen_stat_res(phen_stat_res)

  if(adjust_p){
    phen_stat_res$p <- phen_stat_res$p_adjusted
  } else {
    phen_stat_res$p <- phen_stat_res$p_value
  }

  phen_stat_res <-  phen_stat_res[!is.na(phen_stat_res$p), ,drop = F]

  p_data <- as.data.frame(data.table::dcast(data.table::data.table(phen_stat_res),
                                            gene ~ phen, value.var = c("p")))
  rownames_p_data <- p_data[, 1]
  p_data <- p_data[, -1, drop = F]
  p_data <- apply(p_data, 2, as.numeric)
  rownames(p_data) <- rownames_p_data
  p_data <- p_data[apply(p_data, 1, function(x){min(x, na.rm = T)}) < p_threshold, ,drop = F]
  if(nrow(p_data) == 0){
    stop("no result.")
  }
  p_data <- p_data[order(apply(p_data, 1, function(x){min(x, na.rm = T)})), , drop = F]

  if(only_show_significant){
    p_data <- apply(p_data, 2, function(x){x[x > p_threshold] = NA; x})
  }

  p_data <- p_data[, colSums(p_data, na.rm = T) != 0, drop = F]

  ##

  if(length(pav_obj@gene$info) > 0 && length(add_gene_info) > 0){

    add_gene_info <- match.arg(add_gene_info, names(pav_obj@gene$info), several.ok = T)

    gene_info_data <- as.data.frame(pav_obj@gene$info)
    rownames(gene_info_data) <- pav_obj@gene$name
    gene_info_data <- gene_info_data[rownames(p_data), add_gene_info, drop = F]

    color_info <- get_anno_palette(gene_info_color_list, as.list(gene_info_data))

    anno_param_gene_def_args = list(show = T, width = 5, border = FALSE,
                                    name_size = NULL, name_rot = 90, name_side = "bottom")
    anno_param_gene <- merge_args(anno_param_gene_def_args, anno_param_gene)

    if(anno_param_gene$show){
      if(flip){
        anno_param_gene$name_side <- ifelse(anno_param_gene$name_side == "top", "left", "right")
        anno_param_gene$height <- anno_param_gene$width
        anno_bottom <- get_anno_column(gene_info_data, color_info, anno_param_gene)
        anno_right <- NULL
      } else {
        anno_right <- get_anno_row(gene_info_data, color_info, anno_param_gene)
        anno_bottom <- NULL
      }
    } else {
      anno_right <- NULL
      anno_bottom <- NULL
      color_info <- NULL
    }
  } else {
    anno_right <- NULL
    anno_bottom <- NULL
    color_info <- NULL
  }

  if(flip) p_data <- t(p_data)
  ht <- ComplexHeatmap::Heatmap(p_data,
                                col = p_colors,
                                rect_gp = grid::gpar(col = cell_border_color),
                                name = ifelse(adjust_p, "p_adjusted", "p_value"),
                                right_annotation = anno_right,
                                bottom_annotation = anno_bottom,
                                na_col = na_col,
                                cluster_rows = cluster_rows,
                                clustering_distance_rows = clustering_distance_rows,
                                clustering_method_rows = clustering_method_rows,
                                cluster_columns = cluster_columns,
                                clustering_distance_columns = clustering_distance_columns,
                                clustering_method_columns = clustering_method_columns,
                                row_dend_side = row_dend_side,
                                row_dend_width = row_dend_width,
                                column_dend_side = column_dend_side,
                                column_dend_height = column_dend_height,
                                column_names_rot = column_names_rot,
                                column_names_gp = grid::gpar(fontsize = column_names_size),
                                column_names_side = column_names_side,
                                row_names_rot = row_names_rot,
                                row_names_side = row_names_side,
                                row_names_gp = grid::gpar(fontsize = row_names_size))

  lg_info <- get_legend(color_info, gene_info_data, legend_title_size, legend_text_size, legend_grid_size)

  ComplexHeatmap::draw(ht,
                       auto_adjust = FALSE,
                       heatmap_legend_list = lg_info,
                       merge_legend = T,
                       heatmap_legend_side = legend_side)
}


#' phen_block
#'
#' Plot the result of phenotypic association of a certain phenotype in a block chart.
#'
#' @param pav_obj A PAV object.
#' @param phen_stat_res The result from \code{\link[vPan]{phen_stat}}.
#' @param phen_name The name of phenotype used for grouping.
#' @param add_gene_info A character vector of `gene_info` names. The p_value/p_adjusted will also be added to `gene_info`.
#'
#' @param p_threshold The threshold of p_value/p_adjusted.
#' @param adjust_p A logical value indicating whether adjust p_value.
#' @param only_show_significant A logical value indicating whether only show p_value/p_adjusted that satisfies the condition.
#' @param flip A logical value indicating whether flip cartesian coordinates.
#'
#' @param per_colors A vector of colors or a color mapping function for absence percentage, pass to \code{\link[ComplexHeatmap]{Heatmap}}.
#' @param na_col A string of color for NA values.
#' @param cell_border_color `NA` or a string of color for the border of cells.
#' @param gene_info_color_list A list contains named vector of colors for `gene_info` annotation.
#' e.g. list(source = c("reference" = "red", "novel" = "blue"), length = c("orange", "red"))
#'
#' @param cluster_rows A logical value indicating whether perform clustering on rows.
#' @param clustering_distance_rows Method of measuring distance when clustring on rows, pass to \code{\link[stats]{dist}}.
#' @param clustering_method_rows Method to perform hierarchical clustering on rows, pass to \code{\link[stats]{hclust}}.
#' @param row_dend_side The position of the row dendrogram ("left", "right").
#' @param row_dend_width A \code{\link[grid]{unit}} object for the width of the row dendrogram.
#' @param row_sorted A vector of sorted row names. It only works when `cluster_rows = F`.
#'
#' @param show_row_names  A logical value indicating whether show row names.
#' @param row_names_side The position of row names ("left", "right").
#' @param row_names_size The size of row names.
#' @param row_names_rot The rotation of row names.
#'
#' @param cluster_columns A logical value indicating whether perform clustering on columns.
#' @param clustering_distance_columns Method of measuring distance when clustring on columns, pass to \code{\link[stats]{dist}}.
#' @param clustering_method_columns Method to perform hierarchical clustering on columns, pass to \code{\link[stats]{hclust}}.
#' @param column_dend_side The position of the column dendrogram ("top", "bottom").
#' @param column_dend_height A \code{\link[grid]{unit}} object for the height of the column dendrogram.
#' @param column_sorted A vector of sorted column names. It only works when `cluster_columns = F` and `split_block = F`.
#'
#' @param show_column_names A logical value indicating whether show column names.
#' @param column_names_side The position of column names ("top", "column").
#' @param column_names_size The size of column names.
#' @param column_names_rot The rotation of column names.
#'
#' @param anno_param_gene A list contains parameters for the phenotype annotation. These can be any of the following:
#' "show", "width", "border", "name_size", "name_rot" and "name_side".
#'
#' @param legend_side The position of legend ("top", "bottom", "right", "left").
#' @param legend_title The text for the legend title.
#' @param legend_title_size The size of legend title.
#' @param legend_text_size The size of legend item labels.
#' @param legend_grid_size A \code{\link[grid]{unit}} object for the size of legend grid.
#'
#' @export


phen_block <- function(pav_obj,
                       phen_stat_res,
                       phen_name,
                       add_gene_info = "p",

                       p_threshold = 0.01,
                       adjust_p = T,
                       only_show_significant = T,
                       flip = F,

                       per_colors = c("#1E90FF", "#EE6363"),
                       na_col = "gray",
                       cell_border_color = NA,
                       gene_info_color_list = NULL,

                       cluster_rows = T,
                       clustering_distance_rows = "euclidean",
                       clustering_method_rows = "complete",
                       row_dend_side = "left",
                       row_dend_width = grid::unit(5, "mm"),
                       row_sorted = c(),

                       show_row_names = F,
                       row_names_side = "right",
                       row_names_size = 10,
                       row_names_rot = 0,

                       cluster_columns = T,
                       clustering_distance_columns = "euclidean",
                       clustering_method_columns = "complete",
                       column_dend_side = "top",
                       column_dend_height = grid::unit(5, "mm"),
                       column_sorted = c(),

                       show_column_names = F,
                       column_names_side = "bottom",
                       column_names_size = 10,
                       column_names_rot = 90,

                       anno_param_gene = list(show = T, width = 5, border = FALSE,
                                              name_size = NULL, name_rot = 90, name_side = "bottom"),

                       legend_side = "right",
                       legend_title = list(pav = "PAV", type = "gene"),
                       legend_title_size = NULL,
                       legend_text_size = NULL,
                       legend_grid_size = grid::unit(4, "mm")){

  check_class(pav_obj, "PAV")

  check_phen_stat_res(phen_stat_res)

  if(adjust_p){
    phen_stat_res$p <- phen_stat_res$p_adjusted
  } else {
    phen_stat_res$p <- phen_stat_res$p_value
  }

  pav_data <- pav_obj@pav_data
  check_obj_phen(pav_obj)
  phen_data <- as.data.frame(pav_obj@sample$phen)
  rownames(phen_data) <- pav_obj@sample$name
  phen_name <- match.arg(phen_name, unique(phen_stat_res$phen))

  p_data <- subset(phen_stat_res, phen == phen_name & p < p_threshold)
  p_data <- p_data[order(p_data$p, decreasing = F), ]
  if(nrow(p_data) == 0){
    warning("no result.")
    return(NULL)
  }
  p_pav_data <- pav_data[p_data$gene, , drop = F]
  groups <- unique(phen_data[[phen_name]])

  groups_data <- as.data.frame(lapply(groups, function(x){
    curr_samples <- rownames(phen_data)[phen_data[[phen_name]] == x]
    rowSums(p_pav_data[, curr_samples, drop = F]) / length(curr_samples)
  }))

  colnames(groups_data) <- paste0(names(table(phen_data[[phen_name]])),
                                  "(n=", table(phen_data[[phen_name]]), ")")

  ##

  gene_data <- p_data[, c(ifelse(adjust_p, "p_adjusted", "p_value"), "gene"), drop = F]
  if(length(pav_obj@gene$info) > 0){
    gene_data <- merge(gene_data,
                       cbind(gene = pav_obj@gene$name, as.data.frame(pav_obj@gene$info)),
                       by = "gene")
  }
  rownames(gene_data) <- gene_data$gene
  gene_data <- gene_data[, -1]

  if(length(add_gene_info) > 0){

    add_gene_info <- match.arg(add_gene_info, colnames(gene_data), several.ok = T)
    gene_info_data <- gene_data[rownames(groups_data), add_gene_info, drop = F]

    color_info <- get_anno_palette(gene_info_color_list, as.list(gene_info_data))

    anno_param_gene_def_args = list(show = T, width = 5, border = FALSE,
                                    name_size = NULL, name_rot = 90, name_side = "bottom")
    anno_param_gene <- merge_args(anno_param_gene_def_args, anno_param_gene)

    if(anno_param_gene$show){
      if(flip){
        anno_param_gene$name_side <- ifelse(anno_param_gene$name_side == "top", "left", "right")
        anno_param_gene$height <- anno_param_gene$width
        anno_bottom <- get_anno_column(gene_info_data, color_info, anno_param_gene)
        anno_right <- NULL
      } else {
        anno_right <- get_anno_row(gene_info_data, color_info, anno_param_gene)
        anno_bottom <- NULL
      }
    } else {
      anno_right <- NULL
      anno_bottom <- NULL
      color_info <- NULL
    }
  } else {
    anno_right <- NULL
    anno_bottom <- NULL
    color_info <- NULL
  }

  if(flip)  groups_data <- t(groups_data)

  ht <- ComplexHeatmap::Heatmap(as.matrix(groups_data),
                                name = "Presence(%)",
                                col = per_colors,
                                rect_gp = grid::gpar(col = cell_border_color),
                                right_annotation = anno_right,
                                bottom_annotation = anno_bottom,
                                na_col = na_col,
                                cluster_rows = cluster_rows,
                                clustering_distance_rows = clustering_distance_rows,
                                clustering_method_rows = clustering_method_rows,
                                cluster_columns = cluster_columns,
                                clustering_distance_columns = clustering_distance_columns,
                                clustering_method_columns = clustering_method_columns,
                                row_dend_side = row_dend_side,
                                row_dend_width = row_dend_width,
                                column_dend_side = column_dend_side,
                                column_dend_height = column_dend_height,
                                column_names_rot = column_names_rot,
                                column_names_gp = grid::gpar(fontsize = column_names_size),
                                column_names_side = column_names_side,
                                row_names_rot = row_names_rot,
                                row_names_side = row_names_side,
                                row_names_gp = grid::gpar(fontsize = row_names_size)
  )

  lg_info <- get_legend(color_info, gene_info_data, legend_title_size, legend_text_size, legend_grid_size)

  ComplexHeatmap::draw(ht,
                       auto_adjust = FALSE,
                       heatmap_legend_list = lg_info,
                       merge_legend = T,
                       heatmap_legend_side = legend_side)
}



#' phen_manhattan
#'
#' Plot the result of phenotypic association of a certain phenotype in a Manhattan plot.
#'
#' @param pav_obj A PAV object.
#' @param phen_stat_res The result from \code{\link[vPan]{phen_stat}}.
#' @param phen_name The name of phenotype used for displaying.
#' @param chr The name in `gene_info` denoting chromosomes.
#' @param bp The name in `gene_info` denoting chromosomal position.
#'
#' @param adjust_p A logical value indicating whether adjust p_value.
#' @param highlight_top_n The top `n` points will be highlighted.
#' @param highlight_text_size The size of labels on highlight points.
#'
#' @param point_size The size of points.
#' @param x_text_size The size of tick labels on x-axis.
#' @param x_title_size The size of x-axis title.
#' @param y_text_size The size of tick labels on y-axis.
#' @param y_title_size The size of y-axis title.
#'
#' @importFrom ggrepel geom_text_repel
#'
#' @export

phen_manhattan <- function(pav_obj,
                           phen_stat_res,
                           phen_name,
                           chr,
                           bp,
                           adjust_p = T,
                           # line_value = -log10(1e-05),

                           highlight_top_n = 5,
                           highlight_text_size = 4,
                           point_size = 1.5,
                           x_text_size = NULL,
                           x_title_size = NULL,
                           y_text_size = NULL,
                           y_title_size = NULL
){

  check_class(pav_obj, "PAV")

  check_phen_stat_res(phen_stat_res)
  check_obj_phen(pav_obj)
  phen_name <- match.arg(phen_name, names(pav_obj@sample$phen))

  check_obj_gene(pav_obj)
  gene_list <- pav_obj@gene
  chr <- match.arg(chr, names(gene_list$info))
  if(!is.numeric(gene_list$info[[chr]])){
    stop("the column of `chr` should be numerical")
  }
  bp <- match.arg(bp, names(gene_list$info))
  if(!is.numeric(gene_list$info[[bp]])){
    stop("the column of `bp` should be numerical")
  }
  d <- data.frame(GENE = gene_list$name, CHR = gene_list$info[[chr]], BP = gene_list$info[[bp]])
  d <- merge(d, subset(phen_stat_res, phen == phen_name), by.x = "GENE", by.y = "gene")
  if(adjust_p){
    d$P <- d$p_adjusted
  } else {
    d$P <- d$p_value
  }
  d$CHR <- as.numeric(d$CHR)
  d$BP <- as.numeric(d$BP)
  d <- d[order(d$CHR, d$BP),]
  d[is.na(d$CHR), "CHR"] <- "NA"
  d$pos <- 1:nrow(d)
  lengths <- unlist(lapply(unique(d$CHR), function(x){nrow(subset(d, CHR == x))}))
  ticks <- cumsum(lengths) - lengths/2

  n_CHR <- length(unique(d$CHR))
  CHR_colors <- rep(c("black","gray"),ceiling(length(unique(d$CHR))/2))[1:n_CHR]

  d2 <- subset(d, P <= sort(d$P)[highlight_top_n])

  ggplot(subset(d, !is.na(P)), aes(x = pos, y = -log10(P))) +
    geom_point(aes(color = factor(CHR, levels = unique(CHR))), size = point_size) +
    geom_point(data = d2, aes(x = pos, y = -log10(P)), color = "red", size = point_size) +
    ggrepel::geom_text_repel(data = d2, aes(x = pos, y = -log10(P), label = GENE), size = highlight_text_size) +
    # geom_hline(yintercept = line_value, color = "red") +
    labs(x = "Chromosome", y = paste0("-log10(", ifelse(adjust_p, "p_adjusted", "p_value"), ")")) +
    theme_classic() +
    scale_color_manual(values = CHR_colors) +
    scale_x_continuous(expand = expansion(mult = 0), breaks = ticks, labels = unique(d$CHR)) +
    scale_y_continuous(expand = expansion(mult = c(0,.05))) +
    theme(legend.position = "none",
          strip.background = element_blank(),
          axis.ticks.x = element_blank(),
          axis.line.x = element_blank(),
          axis.text.x = element_text(size = x_text_size),
          axis.title.x = element_text(size = x_title_size),
          axis.text.y = element_text(size = y_text_size),
          axis.title.y = element_text(size = y_title_size)
    )
}



#' phen_bar
#'
#' Draw the presence/absence of a specified gene in a specified phenotype group in a barplot.
#'
#' @param pav_obj A PAV object.
#' @param phen_name The name of phenotype.
#' @param gene_name The name of gene.
#' .
#' @param pav_colors A vector of colors for presence and absence.
#' @param bar_width A numeric vector giving the relative width of bars, ranging from 0 to 1.
#'
#' @param x_text_size The size of tick labels on x-axis.
#' @param x_title_size The size of x-axis title.
#' @param y_text_size The size of tick labels on y-axis.
#' @param y_title_size The size of y-axis title.
#' @param legend_side The position of legend.
#' @param legend_title_size The size of legend title.
#' @param legend_text_size The size of legend item labels.
#'
#' @importFrom ggplot2 ggplot aes geom_bar labs scale_fill_manual scale_y_continuous theme_classic theme
#'
#' @export

phen_bar <- function(pav_obj,
                     phen_name,
                     gene_name,
                     pav_colors = c("#F8766D", "#00BFC4"),
                     bar_width = .8,
                     x_text_size = NULL,
                     x_title_size = NULL,
                     y_text_size = NULL,
                     y_title_size = NULL,
                     legend_side = "top",
                     legend_title_size = NULL,
                     legend_text_size = NULL
){

  check_class(pav_obj, "PAV")

  check_obj_phen(pav_obj)
  phen_name <- match.arg(phen_name, names(pav_obj@sample$phen))
  if(!is.character(pav_obj@sample$phen[[phen_name]])){
    stop("please select a character phen")
  }
  gene_name <- match.arg(gene_name, pav_obj@gene$name)
  pav_data <- pav_obj@pav_data

  p_data <- data.table::melt(
    data.table::data.table(table(data.frame(pav = pav_data[gene_name,],
                                            phen = pav_obj@sample$phen[[phen_name]]))),
    id.vars = c("pav", "phen"))

  ggplot(p_data, aes(x = phen, y = value, fill = as.character(pav))) +
    geom_bar(stat = "identity", position = "stack", width = bar_width) +
    labs(fill  = gene_name, x = phen_name, y = "Sample Number") +
    scale_fill_manual(values = pav_colors, breaks = c(0, 1), labels = c("Absence", "Presence")) +
    scale_y_continuous(expand = expansion(mult= c(0, .05))) +
    theme_classic() +
    theme(axis.text.x = element_text(size = x_text_size),
          axis.title.x = element_text(size = x_title_size),
          axis.text.y = element_text(size = y_text_size),
          axis.title.y = element_text(size = y_title_size),
          legend.position = legend_side,
          legend.title = element_text(size = legend_title_size),
          legend.text = element_text(size = legend_text_size))
}


#' phen_violin
#'
#' Compare a phenotype in absence/presence genes.
#'
#' @param pav_obj A PAV object.
#' @param phen_name The name of phenotype.
#' @param gene_name The name of gene.
#' .
#' @param pav_colors A vector of colors for presence and absence.
#'
#' @param x_text_size The size of tick labels on x-axis.
#' @param x_title_size The size of x-axis title.
#' @param y_text_size The size of tick labels on y-axis.
#' @param y_title_size The size of y-axis title.
#'
#' @param legend_side The position of legend.
#' @param legend_title_size The size of legend title.
#' @param legend_text_size The size of legend item labels.
#'
#' @importFrom ggplot2 ggplot aes geom_violin geom_jitter labs scale_color_manual theme_classic theme scale_x_discrete
#' @importFrom ggsignif geom_signif
#'
#' @export


phen_violin <- function(pav_obj,
                        phen_name,
                        gene_name,
                        pav_colors = c("#F8766D", "#00BFC4"),
                        x_text_size = NULL,
                        x_title_size = NULL,
                        y_text_size = NULL,
                        y_title_size = NULL,
                        legend_side = "top",
                        legend_title_size = NULL,
                        legend_text_size = NULL){

  check_class(pav_obj, "PAV")

  check_obj_phen(pav_obj)
  phen_name <- match.arg(phen_name, names(pav_obj@sample$phen))
  if(!is.numeric(pav_obj@sample$phen[[phen_name]])){
    stop("please select a numerical phen")
  }
  gene_name <- match.arg(gene_name, pav_obj@gene$name)
  pav_data <- pav_obj@pav_data

  p_data <- data.frame(pav = as.character(pav_data[gene_name,]),
                       phen = pav_obj@sample$phen[[phen_name]])
  p_data <- p_data[!is.na(p_data$phen),]

  ggplot(p_data, aes(x = pav, y = phen, color = pav)) +
    geom_violin() + geom_jitter() +
    labs(color = gene_name, y = phen_name, x = "PAV") +
    ggsignif::geom_signif(comparisons = list(c("0","1"))) +
    scale_x_discrete(labels = c("0" = "Absence", "1" = "Presence")) +
    scale_color_manual(values = pav_colors, breaks = c("0", "1"), labels = c("Absence", "Presence")) +
    theme_classic() +
    theme(axis.text.x = element_text(size = x_text_size),
          axis.title.x = element_text(size = x_title_size),
          axis.text.y = element_text(size = y_text_size),
          axis.title.y = element_text(size = y_title_size),
          legend.position = legend_side,
          legend.title = element_text(size = legend_title_size),
          legend.text = element_text(size = legend_text_size))
}









