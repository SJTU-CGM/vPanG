

#' pav_hist
#'
#' The genes can be divided into multiple types based on how many samples are shared.
#'  `pav_hist()` integrate a ring chart and a histogram to showing the number of gene types.
#'
#' @param pav_obj A PAV object.
#' @param show_ring A logical value indicating whether draw ring chart.
#' @param ring_pos_x A numeric vector that specifies the x-location of the ring chart, ranging from 0 to 1.
#' @param ring_pos_y A numeric vector that specifies the y-location of the ring chart, ranging from 0 to 1.
#' @param ring_r The radius of the ring chart, ranging from 0 to 0.5.
#' @param ring_labels_size The size of labels on the ring chart.
#' @param type_colors A named vector of colors for types. e.g. c("distributed" = "red")
#' @param x_title The text for the x-axis title.
#' @param x_title_size The size of x-axis title.
#' @param y_title The text for the y-axis title.
#' @param y_title_size The size of y-axis title.
#' @param x_breaks A numeric vector of break values on the x-axis.
#' @param x_text_size The size of tick labels on the x-axis.
#' @param y_text_size The size of tick labels on the y-axis.
#'
#'
#' @importFrom grid viewport grid.newpage
#' @importFrom ggrepel geom_text_repel
#' @importFrom ggplot2 ggplot geom_tile coord_polar scale_y_continuous scale_fill_manual
#' theme labs theme_classic geom_bar
#'
#' @export


pav_hist <- function(pav_obj,
                     show_ring = T,
                     ring_pos_x = .5,
                     ring_pos_y = .6,
                     ring_r = .3,
                     ring_label_size = NA,
                     type_colors = NULL,
                     x_title = "Sample Number",
                     x_title_size = NULL,
                     y_title = "Count",
                     y_title_size = NULL,
                     x_breaks = NULL,
                     x_text_size = NULL,
                     y_text_size = NULL){

  check_class(pav_obj, "PAV")

  ring_pos_x <- match_num("ring_pos_x", ring_pos_x, 0, 1)
  ring_pos_y <- match_num("ring_pos_y", ring_pos_y, 0, 1)
  ring_r <- match_num("ring_r", ring_r, 0, .5)

  genes_data <- data.frame(pav_obj@gene[1:3], stringsAsFactors = F)
  types <-  c("core", "softcore", "distributed", "private")
  types_col <- get_type_palette(type_colors)

  p_data <- merge(data.frame(table(genes_data$sample_n)),
                  unique(genes_data[,c("sample_n", "type")]),
                  by.x = "Var1", by.y = "sample_n")
  p_data$Var1 <- as.integer(p_data$Var1)

  p1 <- ggplot() +
    geom_bar(data = p_data, aes(x = Var1, y = Freq, fill = factor(type, levels = types)),stat = "identity", width = 0.8) +
    theme_classic() +
    scale_fill_manual(values = types_col) +
    labs(x = "Sample Number", y = "Count", fill = "Gene") +
    scale_y_continuous(expand = expansion(mult=c(0,.05))) +
    theme(legend.position = "none",
          axis.title.x = element_text(size = x_title_size),
          axis.title.y = element_text(size = y_title_size),
          axis.text.x = element_text(size = x_text_size),
          axis.text.y = element_text(size  = y_text_size))

  if(is.numeric(x_breaks)){
    p1 <- p1 + scale_x_continuous(expand = expansion(mult = c(0, .01)), breaks = x_breaks)
  } else {
    p1 <- p1 + scale_x_continuous(expand = expansion(mult = c(0, .01)))
  }

  if(show_ring){
    grid::grid.newpage()
    print(p1)
    vp <- grid::viewport(x = ring_pos_x, y = ring_pos_y, width = 2 * ring_r, height = 2 * ring_r)
    ring_data <- data.frame(table(genes_data$type)[types])
    ring_data <- ring_data[!is.na(ring_data$Freq),]
    ring_data$x <- ring_data$Freq/2 + head(c(0, cumsum(ring_data$Freq)), -1)
    ring_data$per <- ring_data$Freq / sum(ring_data$Freq) * 100

    p2 <- ggplot(ring_data, aes(x = x, y = 1, width = Freq, height = 2)) +
      geom_tile(aes(fill = Var1)) +
      ggrepel::geom_text_repel(aes( y = 2.05, label = paste0(Var1,"\n(",round(per, 1)," %)")),
                               ylim = c(3,NA), direction = "x", nudge_y = 1, size = ring_label_size) +
      coord_polar("x") +
      scale_y_continuous(limits = c(-2, 4)) +
      scale_fill_manual(values = types_col) +
      theme(axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_blank(),
            panel.background = element_blank(),
            plot.background = element_blank(),
            panel.grid = element_blank(),
            legend.position = "none")
    print(p2, vp = vp)
  } else {
    print(p1)
  }

}



