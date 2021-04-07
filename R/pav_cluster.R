

#' pav_cluster
#'
#' Cluster samples based on PAV table.
#'
#' @param pav_obj A PAV object.
#' @param cluster_distance Method to measure distance, pass to \code{\link[stats]{dist}}.
#' @param cluster_method Method to perform hierarchical clustering, pass to \code{\link[stats]{hclust}}.
#'
#' @param add_phen_info A character string of `phen_info` names.
#' @param phen_info_color_list A list contains named vector of colors for `phen_info` annotation.
#' e.g. list(gender = c("Male" = "green", "Female" = "red"))
#'
#' @param sample_name_size The size of labels.
#' @param mult A numebr of multiplicative range expansion factors.
#'
#' @param legend_side The position of legend ("top", "bottom", "right", "left").
#' @param legend_title_size The size of legend title.
#' @param legend_text_size The size of legend item labels.
#'
#' @importFrom ggdendro dendro_data
#' @importFrom ggplot2 ggplot geom_segment geom_text aes labs coord_polar scale_x_continuous theme scale_color_manual
#'
#' @export



pav_cluster <- function(pav_obj,
                        cluster_distance = "euclidean",
                        cluster_method = "complete",

                        add_phen_info = NULL,
                        phen_info_color_list = NULL,

                        sample_name_size = 4,
                        mult = .1,

                        legend_side = "right",
                        legend_title_size = NULL,
                        legend_text_size =NULL
                        ){

  check_class(pav_obj, "PAV")

  pav <- t(pav_obj@pav_data)
  sample <- pav_obj@sample

  dist <- dist(pav, method = cluster_distance)

  cluster <- hclust(dist, method = cluster_method)
  dend_data <- ggdendro::dendro_data(as.dendrogram(cluster), type = "rectangle")

  segment_data <- dend_data$segments
  ymax <- max(segment_data$y)
  label_data <- dend_data$labels
  xmax <- max(dend_data$labels$x)

  if(is.null(add_phen_info) || length(sample$phen) == 0){
    p <- ggplot() +
      geom_segment(data = segment_data, aes(x = x, y = ymax - y, xend = xend, yend = ymax - yend)) +
      geom_text(data = label_data, aes(x = x, y = ymax, label = label, angle = 360*x/xmax - 90), size = sample_name_size, hjust = 1)
    phen_col <- NULL
  } else {
    add_phen_info <- match.arg(add_phen_info, names(sample$phen))
    phen_col <- get_anno_palette(phen_info_color_list, sample$phen[add_phen_info])
    label_data$phen <- sample$phen[[add_phen_info]][match(label_data$label, sample$name)]
    segment_data_up <- subset(segment_data, yend != 0)
    segment_data_down <- subset(segment_data, yend == 0 & x == xend & x %in% label_data$x)
    segment_data_down$phen <- label_data$phen[match(label_data$x, segment_data_down$x)]

    p <- ggplot() +
      geom_segment(data = segment_data_up, aes(x = x, y = ymax - y, xend = xend, yend = ymax - yend), color = "gray50") +
      geom_segment(data = segment_data_down, aes(x = x, y = ymax - y, xend = xend, yend = ymax - yend, color = as.character(phen))) +
      geom_text(data = label_data, aes(x = x, y = ymax, label = label, color = as.character(phen),
                                       angle = 360*x/xmax - 90), hjust = 1, size = sample_name_size)
  }

  p <- p + labs(color = add_phen_info) +
    coord_polar("x", direction = -1) +
    scale_x_continuous(limits = c(0,xmax)) +
    scale_y_continuous(expand = expansion(mult = c(0, mult))) +
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          panel.background = element_blank(),
          plot.background = element_blank(),
          legend.background = element_blank(),
          legend.position = legend_side,
          legend.key = element_blank(),
          legend.title = element_text(size = legend_title_size),
          legend.text = element_text(size = legend_text_size),
          panel.grid = element_blank())

  if(!is.null(phen_col)){
    p <- p +
      scale_color_manual(values = phen_col[[add_phen_info]])
  }

  print(p)

}







