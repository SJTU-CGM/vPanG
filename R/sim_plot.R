

#' sim_plot
#'
#' Plot the result of genome simulation.
#'
#' @param sim_stat_res The result from \code{\link[vPan]{sim_stat}}.
#' @param chart_type A character string of chart type. It should be one of "box", "errorbar", "jitter" and "ribbon".
#'
#' @param box_width A numeric vector giving the relative widths of the box and max/min line, ranging from 0 to 1.
#' @param box_size The size of lines in boxplot.
#' @param box_color A string of color or a named vector of colors for lines in boxplot.
#' @param box_fill A string of color or a named vector of colors for boxes in boxplot.
#' @param box_alpha The opacity of boxes, ranging from 0 to 1.
#' @param box_outlier_size The size of the outliers in boxplot.
#'
#' @param errorbar_width A numeric vector giving the relative width of errorbar, ranging from 0 to 1.
#' @param errorbar_size The size of errorbar.
#' @param errorbar_color A string of color or a named vector of colors for errorbar.
#' @param errorbar_alpha The opacity of errorbar, ranging from 0 to 1.
#' @param errorbar_point_size The size of point representing the mean value.
#' @param errorbar_point_color A string of color or a named vector of colors for point representing the mean value.
#'
#' @param jitter_width A numeric vector giving the relative width of jittered points, ranging from 0 to 1.
#' @param jitter_size The size of jittered points.
#' @param jitter_color A string of color or a named vector of colors for jittered points.
#' @param jitter_alpha The opacity of jittered points, ranging from 0 to 1.
#' @param jitter_point_size The size of point representing the mean value.
#' @param jitter_point_color A string of color or a named vector of colors for point representing the mean value.
#'
#' @param path_size The size of path.
#' @param path_color A string of color or a named vector of colors for path.
#' @param ribbon_fill A string of color or a named vector of colors for ribbon.
#' @param ribbon_alpha The opacity of ribbon, ranging from 0 to 1.
#'
#' @param x_title The text for the x-axis title.
#' @param y_title The text for the y-axis title.
#' @param x_title_size The size of x-axis title.
#' @param y_title_size The size of y-axis title.
#' @param x_breaks A numeric vector of break values on the x-axis.
#' @param y_breaks A numeric vector of break values on the y-axis.
#' @param x_text_size The size of tick labels on x-axis.
#' @param y_text_size The size of tick labels on y-axis.
#'
#' @param legend_show Whether show the legend.
#' @param legend_side The position of legend.
#' @param legend_title The text for the legend title.
#' @param legend_title_size The size of legend title.
#' @param legend_text_size The size of legend item labels.
#'
#' @importFrom ggplot2 stat_boxplot geom_boxplot element_rect geom_errorbar geom_jitter geom_ribbon geom_path
#'
#' @export


sim_plot <- function(
  sim_stat_res,
  chart_type = "box",

  box_width = c(.8, .8),
  box_size = 1,
  box_color = NULL,
  box_fill = NULL,
  box_alpha = .5,
  box_outlier_size = NULL,

  errorbar_width = .8,
  errorbar_size = 1,
  errorbar_color = NULL,
  errorbar_alpha = .8,
  errorbar_point_size = 2,
  errorbar_point_color = NULL,

  jitter_width = .8,
  jitter_size = 1,
  jitter_color = NULL,
  jitter_alpha = .1,
  jitter_point_size = 2,
  jitter_point_color = NULL,

  path_size = 1,
  path_color = NULL,
  ribbon_fill = NULL,
  ribbon_alpha = 0.5,

  x_title = "Times",
  y_title = "Count",
  x_title_size = NULL,
  y_title_size = NULL,
  x_breaks = NULL,
  y_breaks = NULL,
  x_text_size = NULL,
  y_text_size = NULL,

  legend_show = T,
  legend_side = "top",
  legend_title = "",
  legend_title_size = NULL,
  legend_text_size = NULL

){

  match.arg(colnames(sim_stat_res), c("Times", "pan", "core", "private"), several.ok = T)
  sim_data <- sim_stat_res

  sim_data_l <- data.table::melt(data.table::data.table(sim_data), id.vars = "Times")
  sim_data_stat <- sim_data_l[,.(Length=mean(value),SD=sd(value)),by=c("variable","Times")]
  gene_groups <- colnames(sim_data)[-1]

  chart_type_i <- match.arg(chart_type, choices = c("box", "errorbar", "jitter", "ribbon"))

  if(chart_type_i == "box"){

    if(length(box_width) == 1){
      widths <- rep(box_width, 2)
    } else if(length(box_width) == 2){
      widths <- box_width
    } else {
      stop("box_width should be one or two length vector.")
    }

    colors <- get_color(box_color, gene_groups)
    fills <- get_color(box_fill, gene_groups)

    p <- ggplot()
    for(i in gene_groups){
      p_data <- subset(sim_data_l, variable == i)
      p <- p +
        stat_boxplot(data = p_data,
                     aes(x = Times, y = value, group = Times),
                     geom = "errorbar",
                     width = widths[2],
                     size = box_size,
                     color = as.vector(colors[i])) +
        geom_boxplot(data = p_data,
                     aes(x = Times, y = value, group = Times),
                     width = widths[1],
                     size = box_size,
                     color = "white",
                     fill = "white",
                     outlier.shape = NA)
    }
    for(i in gene_groups){
      p_data <- subset(sim_data_l, variable == i)
      p <- p +
        geom_boxplot(data = p_data,
                     aes(x = Times, y = value, group = Times,
                         color = variable, fill = variable),
                     width = widths[1],
                     size = box_size,
                     alpha = box_alpha,
                     outlier.size = box_size,
                     outlier.color = as.vector(colors[i]),
                     outlier.shape = box_outlier_size)
    }
    p <- p + scale_fill_manual(values = fills) +
      scale_color_manual(values = colors)

  } else if(chart_type_i == "errorbar"){
    colors <- get_color(errorbar_color, gene_groups)
    point_colors <- get_color(errorbar_point_color, gene_groups)
    p <- ggplot(data = sim_data_stat,
                aes(x = Times, y = Length, color = variable)) +
      geom_errorbar(aes(ymin = Length - SD, ymax = Length + SD),
                    width = errorbar_width,
                    size = errorbar_size,
                    alpha =  errorbar_alpha)+
      geom_point(aes(fill = variable), shape = 21,
                 size = errorbar_point_size)
    p <- p + scale_fill_manual(values = point_colors) +
      scale_color_manual(values = colors)
  } else if(chart_type_i == "jitter"){
    colors <- get_color(jitter_color, gene_groups)
    point_colors <- get_color(jitter_point_color, gene_groups)
    p <- ggplot() +
      geom_jitter(data = sim_data_l,
                  aes(x = Times, y = value, color = variable),
                  width = jitter_width,
                  height = 0,
                  size = jitter_size,
                  alpha = jitter_alpha)+
      geom_point(data = sim_data_stat,
                 aes(x = Times, y = Length, fill = variable, color = variable),
                 shape = 21,
                 size = jitter_point_size)
    p <- p + scale_fill_manual(values = point_colors) +
      scale_color_manual(values = colors)
  } else if (chart_type_i == "ribbon"){
    colors <- get_color(path_color, gene_groups)
    fills <- get_color(ribbon_fill, gene_groups)
    p <- ggplot()+
      geom_ribbon(data = sim_data_stat,
                  aes(x = Times, y = Length,
                      ymin = Length - SD, ymax = Length + SD,
                      fill = variable),
                  alpha = ribbon_alpha)+
      geom_path(data = sim_data_stat,
                aes(x = Times, y = Length, color = variable),
                size = path_size)
    p <- p + scale_fill_manual(values = fills) +
      scale_color_manual(values = colors)
  }

  if(!is.null(x_breaks)){
    p <- p + scale_x_continuous(breaks = x_breaks)
  }

  if(!is.null(y_breaks)){
    p <- p + scale_y_continuous(breaks = y_breaks)
  }

  p <- p +
    labs(
      x = x_title,
      y = y_title,
      color = legend_title,
      fill = legend_title
    ) + theme_classic()

  p <- p +
    theme(
    axis.title.x = element_text(size = x_title_size),
    axis.title.y = element_text(size = y_title_size),
    axis.text.x = element_text(size = x_text_size),
    axis.text.y = element_text(size = y_text_size),
    panel.grid = element_blank(),
    panel.background = element_blank()
  )

  p <- p + theme(
    legend.position = ifelse(legend_show, legend_side, "none"),
    legend.direction = ifelse(legend_side %in% c("top", "bottom"),
                              "horizontal", "vertical"),
    legend.title = element_text(size = legend_title_size),
    legend.text = element_text(size = legend_text_size),
    legend.key = element_rect(fill = NA, color = NA)
  )

  print(p)

}



