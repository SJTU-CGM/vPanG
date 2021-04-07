

#' pav_stackbar
#'
#' The composition of genes in all samples can be viewed in `pav_stackbar()`.
#' The chart consists of a hierarchical cluster tree and a bar plot.
#'
#' @param pav_obj A PAV object.
#' @param show_relative A logical value indicating whether show relative value.
#' @param add_phen_info A character string of `phen_info` names.
#'
#' @param type_colors A named vector of colors for types. e.g. c("distributed" = "red")
#' @param phen_info_color_list A list contains named vector of colors for `phen_info` annotation.
#' e.g. list(gender = c("Male" = "green", "Female" = "red"))
#'
#' @param cluster_distance Method to measure distance, pass to \code{\link[stats]{dist}}.
#' @param cluster_method Method to perform hierarchical clustering, pass to \code{\link[stats]{hclust}}.
#'
#' @param bar_width A numeric vector giving the relative width of bars, ranging from 0 to 1.
#' @param sample_name_size The size of sample names.
#'
#' @param legend_side The position of legend ("top", "bottom", "right", "left").
#' @param legend_title The text for the legend title.
#' @param legend_title_size The size of legend title.
#' @param legend_text_size The size of legend item labels.
#'
#' @param dend_width The relative width of dendrogram, ranging from 0 to 1.
#' @param name_width The relative width of sample names, ranging from 0 to 1.
#'
#' @importFrom ggdendro dendro_data
#' @importFrom ggrepel geom_text_repel
#' @importFrom ggplot2 ggplot geom_bar scale_fill_manual geom_segment coord_flip labs
#' theme aes geom_text geom_point scale_color_manual element_text element_blank
#'
#' @export

pav_stackbar <- function(pav_obj,
                         show_relative = F,
                         add_phen_info = NULL,

                         type_colors = NULL,
                         phen_info_color_list = NULL,

                         cluster_distance = "euclidean",
                         cluster_method = "complete",

                         bar_width = 1,
                         sample_name_size = 4,

                         legend_side = "right",
                         legend_title = "Gene",
                         legend_title_size = NULL,
                         legend_text_size = NULL,

                         dend_width = .05,
                         name_width = .1){

  check_class(pav_obj, "PAV")

  show_relative <- match_logi("show_relative", show_relative)
  dend_width <- match_num("dend_width", dend_width, 0, 1)
  name_width <- match_num("name_width", name_width, 0, 1)
  if(dend_width + name_width > 1){
    stop("the sum of `dend_width` and `name_width` should be less than 1. ")
  }

  genes_data <- data.frame(pav_obj@gene[1:3], stringsAsFactors = F)
  sample <- pav_obj@sample
  sample_n <- length(sample$name)

  types <- unique(pav_obj@gene$type)

  p_data <- do.call(cbind, lapply(types, function(x){
    freq <- colSums(pav_obj@pav_data[subset(genes_data, type == x)$name, ])
    res <- data.frame(freq = freq)
    colnames(res) <- x
    res
  }))

  if(show_relative){
    p_data <- data.frame(t(apply(p_data,1,function(x){sum <- sum(x); x/sum*100})))
  }

  bar_data <- do.call(rbind,
                      lapply(1:length(types), function(x){
                        data.frame(freq = p_data[[x]],
                                   sample = rownames(p_data),
                                   type = colnames(p_data)[x])}))

  bar_len <- max(rowSums(p_data))
  total_len <-  bar_len / (1 - dend_width -name_width)

  types_col <- get_type_palette(type_colors)

  dend_data <- ggdendro::dendro_data(as.dendrogram(
    hclust(dist(p_data, method = cluster_distance), method = cluster_method)),
    type = "rectangle")

  segment_data <- dend_data$segments
  dend_max <- max(segment_data$y)
  dend_len <- total_len * dend_width
  segment_data$y <-  segment_data$y / dend_max * dend_len
  segment_data$yend <-  segment_data$yend / dend_max * dend_len
  label_data <- dend_data$labels
  name_len  <- total_len * name_width

  bar_data <- merge(bar_data, label_data, by.x = "sample", by.y = "label")

  y_breaks = round(unlist(lapply(1:length(types), function(x){mean(rowSums(p_data[, 1:x, drop = F]))})), 1)

  p <- ggplot() +
    geom_bar(data = bar_data, aes(x = x, y = freq, fill = factor(type, levels  = rev(types))), stat = "identity", width = bar_width) +
    scale_fill_manual(values = types_col) +
    geom_segment(data = segment_data, aes(x = x, y = -y - name_len, xend = xend, yend = -yend -name_len)) +
    coord_flip() +
    geom_segment(aes(x = 0.5, xend = sample_n + 0.5, y = y_breaks, yend = y_breaks), linetype = "dashed") +
    labs(fill = legend_title) +
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          panel.background = element_blank(),
          plot.background = element_blank(),
          legend.key = element_blank(),
          legend.title = element_text(size = legend_title_size),
          legend.text = element_text(size = legend_text_size),
          legend.position = legend_side,
          panel.grid = element_blank())

  if(show_relative){
    p <- p +
      ggrepel::geom_text_repel(aes(x = 0.5, y = y_breaks, label = paste0(y_breaks, "%")),
                               xlim=c(NA, -.5), direction = "x", hjust = .5, vjust = .5, size = sample_name_size)
  } else {
    p <- p +
      ggrepel::geom_text_repel(aes(x = 0.5, y = y_breaks, label = y_breaks), xlim=c(NA, -.5),
                               direction = "x", hjust = .5, vjust = .5, size = sample_name_size)
  }

  if(length(sample$phen) == 0 ||is.null(add_phen_info)){
    p <- p +
      geom_text(data = label_data, aes(x = x, y = -name_len/2, label = label), size = sample_name_size)
    phen_col <- NULL
  } else {
    add_phen_info <- match.arg(add_phen_info, names(sample$phen))
    phen_col <- get_anno_palette(phen_info_color_list, sample$phen[add_phen_info])
    label_data$phen <- sample$phen[[add_phen_info]][match(sample$name, label_data$label)]
    p <- p +
      geom_text(data = label_data, aes(x = x, y = -name_len/2, label = label, color = phen),  size = sample_name_size) +
      labs(color = add_phen_info) +
      geom_point(data = label_data, aes(x = x, y = -name_len, color = phen))
  }

  if(!is.null(phen_col)){
    p <- p +
      scale_color_manual(values = phen_col[[add_phen_info]])
  }


  print(p)
}
















