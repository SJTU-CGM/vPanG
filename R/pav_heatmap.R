

#' pav_heatmap
#'
#' Plot a heatmap for a object of PAV class.
#'
#' @param pav_obj A PAV object.
#' @param gene_type A vector of gene types. These can be any of the following: "core", "softcore", "distributed" and "private".
#' @param add_phen_info A character vector of `phen_info` names.
#' @param add_gene_info A character vector of `gene_info` names.
#'
#' @param pav_colors  A named vector of colors for presence and absence. e.g. c(presence = "steelblue", absence = "gray70")
#' @param type_colors A named vector of colors for types. e.g. c("distributed" = "red")
#' @param gene_info_color_list A list contains named vector of colors for `gene_info` annotation.
#' e.g. list(source = c("reference" = "red", "novel" = "blue"), length = c("orange", "red"))
#' @param phen_info_color_list A list contains named vector of colors for `phen_info` annotation.
#' e.g. list(gender = c("Male" = "green", "Female" = "red"), age = c("yellow", "red"))
#'
#' @param border A logical value or a string of color indicating whether draw border.
#'
#' @param split_block A logical value indicating whether split columns based on gene types.
#' @param block_name_size The size of block name.
#' @param block_name_rot The rotation of block name.
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
#' @param anno_param_row_phen A list contains parameters for the phenotype annotation. These can be any of the following:
#' "show", "width", "border", "name_size", "name_rot" and "name_side".
#' @param anno_param_column_gene A list contains parameters for the gene annotation. These can be any of the following:
#' "show", "height", "border", "name_size", "name_rot" and "name_side".
#' @param anno_param_row_stat A list contains parameters for the stat annotation of rows. These can be any of the following:
#' "show", "width", "border", "title", "title_size", "title_rot", "title_side", "axis_side", "axis_at", "axis_labels" and "axis_labels_size".
#' @param anno_param_column_stat A list contains parameters for the stat annotation of columns. These can be any of the following:
#' "show", "height", "border", "title", "title_size", "title_rot", "title_side", "axis_side", "axis_at", "axis_labels" and "axis_labels_size".
#'
#' @param legend_side The position of legend ("top", "bottom", "right", "left").
#' @param legend_title The text for the legend title.
#' @param legend_title_size The size of legend title.
#' @param legend_text_size The size of legend item labels.
#' @param legend_grid_size A \code{\link[grid]{unit}} object for the size of legend grid.
#'
#' @param use_raster Whether render the heatmap body as a raster image. pass to \code{\link[ComplexHeatmap]{Heatmap}}
#'
#' @importFrom randomcoloR distinctColorPalette
#'
#' @export


pav_heatmap <- function(
  pav_obj,
  gene_type,
  add_phen_info = NULL,
  add_gene_info = NULL,

  border = T,
  split_block = T,  # when sort is null.
  block_name_size = NULL,
  block_name_rot = 0,

  cluster_rows = F,
  clustering_distance_rows = "euclidean",
  clustering_method_rows = "complete",
  row_dend_side = "left",
  row_dend_width = grid::unit(5, "mm"),
  row_sorted = c(),  #only work when `cluster_rows` = F

  show_row_names = F,
  row_names_side = "right",
  row_names_size = 10,
  row_names_rot = 0,

  cluster_columns = F,
  clustering_distance_columns = "euclidean",
  clustering_method_columns = "complete",
  column_dend_side = "top",
  column_dend_height = grid::unit(5, "mm"),
  column_sorted = c(), #only work when `columncluster` = F

  show_column_names = F,
  column_names_side = "bottom",
  column_names_size = 10,
  column_names_rot = 90,

  anno_param_row_phen = list(show = T,  width = 5, border = F,
                       name_size = NULL, name_rot = 90, name_side = "top"),
  anno_param_column_gene = list(show = T, height = 5, border = FALSE,
                          name_size = NULL, name_rot = 0, name_side = "right"),
  anno_param_row_stat = list(show = T, width = 10, bar_width = 1, border = FALSE,
                     title = "Presence\nNumber", title_size = 10, title_side = "bottom", title_rot = 0,
                     axis_side = "bottom",axis_at = NULL, axis_labels = NULL, axis_labels_size = 8),
  anno_param_column_stat = list(show = T, height = 10, bar_width = 1, border = FALSE,
                        title = "Presence\nNumber", title_size = 10, title_side = "left", title_rot = 0,
                        axis_side = "left",axis_at = NULL, axis_labels = NULL, axis_labels_size = 8),

  pav_colors = c(presence = "steelblue", absence = "gray70"),
  type_colors = NULL,
  gene_info_color_list = NULL,
  phen_info_color_list = NULL,

  legend_side = "right",
  legend_title = list(pav = "PAV", type = "gene"),
  legend_title_size = NULL,
  legend_text_size = NULL,
  legend_grid_size = grid::unit(4, "mm"),

  use_raster = NULL
){

  check_class(pav_obj, "PAV")

  pav_phen <- pav_obj@sample$phen
  if(length(pav_phen) > 0 && !is.null(add_phen_info)){
    add_phen_info <- match.arg(add_phen_info, names(pav_phen), several.ok = T)
    pav_phen <- pav_phen[names(pav_phen) %in% add_phen_info]
  } else {
    pav_phen <- NULL
  }

  pav_gene_info <- pav_obj@gene$info
  if(length(pav_gene_info) > 0 && !is.null(add_gene_info)){
    add_gene_info <- match.arg(add_gene_info, names(pav_gene_info), several.ok = T)
    pav_gene_info <- pav_gene_info[names(pav_gene_info) %in% add_gene_info]
  } else {
    pav_gene_info <- NULL
  }

  gene_type <- match.arg(gene_type, c("core", "softcore", "distributed", "private"), several.ok = TRUE)
  split_block <- match_logi("split_block", split_block)

  ## gene_data
  genes_data <- data.frame(pav_obj@gene[1:3], stringsAsFactors = F)
  if(length(pav_gene_info) != 0){
    genes_data <- cbind(genes_data, data.frame(pav_gene_info, stringsAsFactors = F))
  }
  rownames(genes_data) <- genes_data$name

  ## data_main
  data_main <- t(pav_obj@pav_data)

  if(!cluster_columns){
    if(length(column_sorted) == nrow(genes_data) && all(column_sorted %in% rownames(genes_data))){
      data_main <- data_main[, column_sorted]
      genes_data <- genes_data[column_sorted, ]
      split_block <- F
    } else {
      data_main <- data_main[, names(sort(colSums(data_main), decreasing = TRUE))]
      genes_data <- genes_data[colnames(data_main), ]
    }
  }

  genes_data <- subset(genes_data, type %in% gene_type)
  data_main <- data_main[, rownames(genes_data)]

  ##sample_data
  samples_data <- data.frame(pav_obj@sample[1], stringsAsFactors = F)
  if(length(pav_phen) != 0){
    samples_data <- cbind(samples_data, data.frame(pav_phen, stringsAsFactors = F))
  }
  rownames(samples_data) <- samples_data$name

  if(!cluster_rows){
    if(length(row_sorted) == nrow(data_main) && all(row_sorted %in% rownames(data_main))){
      data_main <- data_main[row_sorted,]
      samples_data <- samples_data[row_sorted,]
    } else {
      data_main <- data_main[names(sort(rowSums(data_main), decreasing = TRUE)), ]
      samples_data <- samples_data[rownames(data_main), ]
    }
  }

  ## colors
  color_pav_def <- c(presence = "steelblue", absence = "gray70")
  if(!is.vector(pav_colors) || is.null(names(pav_colors)) || length(pav_colors) != 2 ||
     !all(c("presence", "absence") %in% names(pav_colors))){
    warning("`pav_colors` shoud be a named vector.")
    color_pav <- color_pav_def
  } else {
    color_pav <- merge_args(color_pav_def, pav_colors)
  }

  color_type <- get_type_palette(type_colors)[gene_type]

  phen_info_color_list_info <- get_anno_palette(c(phen_info_color_list, gene_info_color_list), c(pav_phen, pav_gene_info))
  color_phen <- phen_info_color_list_info[names(pav_phen)]
  color_info <- phen_info_color_list_info[names(pav_gene_info)]

  ## anno param

  anno_param_row_phen_def_args <- list(show = T,  width = 5, border = F,
                                 name_size = NULL, name_rot = 90, name_side = "top")
  anno_param_row_phen <- merge_args(anno_param_row_phen_def_args, anno_param_row_phen)
  anno_param_column_gene_def_args = list(show = T, height = 5, border = FALSE,
                                   name_size = NULL, name_rot = 0, name_side = "right")
  anno_param_column_gene <- merge_args(anno_param_column_gene_def_args, anno_param_column_gene)
  anno_param_row_stat_def_args = list(show = T, width = 10, bar_width = 1, border = FALSE,
                              title = "Presence\nNumber", title_size = 10, title_side = "bottom", title_rot = 0,
                              axis_side = "bottom",axis_at = NULL, axis_labels = NULL, axis_labels_size = 8)
  anno_param_row_stat <- merge_args(anno_param_row_stat_def_args, anno_param_row_stat)
  anno_param_column_stat_def_args = list(show = T, height = 10, bar_width = 1, border = FALSE,
                                 title = "Presence\nNumber", title_size = 10, title_side = "left", title_rot = 0,
                                 axis_side = "left",axis_at = NULL, axis_labels = NULL, axis_labels_size = 8)
  anno_param_column_stat <- merge_args(anno_param_column_stat_def_args, anno_param_column_stat)

  ## anno_left

  data_samplePN <- lapply(gene_type, function(x){
    rowSums(data_main[, subset(genes_data, type == x)$name])})

  if(anno_param_row_stat$show){
    anno_param_row_stat_color <- color_type
    anno_left <- ComplexHeatmap::rowAnnotation(
      PN = ComplexHeatmap::anno_barplot(
        data_samplePN, border = anno_param_row_stat$border,
        width = grid::unit(anno_param_row_stat$width, 'mm'),
        bar_width = anno_param_row_stat$bar_width,
        gp = grid::gpar(fill = anno_param_row_stat_color, col = NA),
        axis_param = list(side = anno_param_row_stat$axis_side,
                          at = anno_param_row_stat$axis_at,
                          labels = anno_param_row_stat$axis_labels,
                          direction = "reverse",
                          gp = grid::gpar(fontsize = anno_param_row_stat$axis_labels_size))
      ), annotation_label = anno_param_row_stat$title,
      annotation_name_side = anno_param_row_stat$title_side,
      annotation_name_rot = anno_param_row_stat$title_rot,
      annotation_name_gp = grid::gpar(fontsize = anno_param_row_stat$title_size)
    )
  }else{
    anno_left <- NULL
  }

  ## anno_right

  if(length(pav_phen) > 0){
    data_phen <- samples_data[, names(pav_phen), drop = F]

    if(anno_param_row_phen$show){
      anno_right <- get_anno_row(data_phen, color_phen, anno_param_row_phen)
        } else {
          anno_right <- NULL
    }
  } else {
    anno_right <- NULL
  }


  ## anno_top

  data_presenceN <- colSums(data_main)

  if(anno_param_column_stat$show){
    anno_param_column_stat_color <- color_type[match(genes_data$type, gene_type)]
    anno_top <- ComplexHeatmap::HeatmapAnnotation(
      PN = ComplexHeatmap::anno_barplot(
        data_presenceN, bar_width = anno_param_column_stat$bar_width, border = anno_param_column_stat$border,
        height = grid::unit(anno_param_column_stat$height, 'mm'),
        gp = grid::gpar(fill = anno_param_column_stat_color, col = NA),
        axis_param = list(side = anno_param_column_stat$axis_side,
                          at = anno_param_column_stat$axis_at,
                          labels = anno_param_column_stat$axis_labels,
                          gp = grid::gpar(fontsize = anno_param_column_stat$axis_labels_size))
      ), annotation_label = anno_param_column_stat$title,
      annotation_name_side = anno_param_column_stat$title_side,
      annotation_name_rot = anno_param_column_stat$title_rot,
      annotation_name_gp = grid::gpar(fontsize = anno_param_column_stat$title_size)
    )
  } else {
    anno_top <- NULL
  }

  #anno_bottom

  if(length(pav_gene_info) > 0){
    data_info <- genes_data[, names(pav_gene_info), drop = F]

    if(anno_param_column_gene$show){
      anno_bottom <- get_anno_column(data_info, color_info, anno_param_column_gene)
    } else {
      anno_bottom <- NULL
    }
  } else {
    anno_bottom <- NULL
  }


  lg <- list(
    ComplexHeatmap::Legend(
      labels = names(color_pav),
      legend_gp = grid::gpar(fill = color_pav),
      title = ifelse(is.null(legend_title$pav), "PAV", legend_title$pav),
      title_gp = grid::gpar(fontsize = legend_title_size, fontface = "bold"),
      labels_gp = grid::gpar(fontsize = legend_text_size),
      grid_height = legend_grid_size,
      grid_width = legend_grid_size
    ), ComplexHeatmap::Legend(
      labels = gene_type,
      legend_gp = grid::gpar(fill = color_type),
      title = ifelse(is.null(legend_title$type), "Gene", legend_title$type),
      title_gp = grid::gpar(fontsize = legend_title_size, fontface = "bold"),
      labels_gp = grid::gpar(fontsize = legend_text_size),
      grid_height = legend_grid_size,
      grid_width = legend_grid_size
    ))

  lg_phen <- get_legend(color_phen, pav_phen, legend_title_size, legend_text_size, legend_grid_size)
  lg_info <- get_legend(color_info, pav_gene_info, legend_title_size, legend_text_size, legend_grid_size)

  if(split_block){
    column_split <- factor(genes_data$type, levels = gene_type)
  }else {
    column_split <- NULL
  }

  ht_main <- ComplexHeatmap::Heatmap(
    data_main,
    name = "main",
    use_raster = use_raster,
    col = as.vector(color_pav[c("absence", "presence")]),
    show_heatmap_legend = F,
    column_split =column_split,
    column_gap = grid::unit(0, "mm"),
    border = border,
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
    column_title_rot = block_name_rot,
    column_title_side = "top",
    column_title_gp = grid::gpar(fontsize = block_name_size, fontface = "bold"),
    show_column_names = show_column_names,
    cluster_column_slices = F,
    show_row_names = ifelse(show_row_names & row_names_side == "left", T, F),
    row_names_rot = row_names_rot,
    row_names_side = row_names_side,
    row_names_gp = grid::gpar(fontsize = row_names_size),
    left_annotation = anno_left,
    bottom_annotation = anno_bottom,
    top_annotation = anno_top
  )

  ht_right <- ComplexHeatmap::Heatmap(
    matrix(NA,ncol = 0, nrow=nrow(data_main),
           dimnames = list(rownames(data_main))),
    show_heatmap_legend = F,
    rect_gp =  grid::gpar(type = "none"),
    show_row_names = ifelse(show_row_names & row_names_side == "right", T, F),
    row_names_rot = row_names_rot,
    row_names_side = row_names_side,
    row_names_gp = grid::gpar(fontsize = row_names_size),
    show_column_names = F,
    right_annotation = anno_right
  )

  ComplexHeatmap::draw(ht_main + ht_right,
                       main_heatmap = "main",
                       auto_adjust = FALSE,
                       heatmap_legend_list = c(lg,lg_phen, lg_info),
                       merge_legend = T,
                       heatmap_legend_side = legend_side)

}













