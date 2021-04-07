



#' pav_halfviolin
#'
#' Plot a half-violin plot for a object of PAV class.
#'
#' @param pav_obj A PAV object.
#' @param violin_color A string of color for half-violin plot.
#' @param add_phen_info A character string of `phen_info` names.
#' @param phen_info_color_list A list contains named vector of colors for `phen_info` annotation.
#' e.g. list(gender = c("Male" = "green", "Female" = "red"))
#' @param x_text_size The size of tick labels on x-axis.
#' @param y_text_size The size of tick labels on y-axis.
#' @param x_title_size The size of x-axis title.
#' @param y_title_size The size of y-axis title.
#'
#' @importFrom ggplot2 geom_polygon geom_jitter scale_x_continuous theme_classic theme
#'
#' @export


pav_halfviolin <- function(pav_obj,
                           violin_color = "#F8766D",
                       add_phen_info = NULL,
                       phen_info_color_list = NULL,
                       x_text_size  = NULL,
                       y_text_size = NULL,
                       x_title_size = NULL,
                       y_title_size = NULL){

  pav_data <- pav_obj@pav_data

  sample_data <-
    data.frame(
      sample = colnames(pav_data),
      number = colSums(pav_data)
    )

  sample <- pav_obj@sample
  if(!is.null(add_phen_info) && length(pav_obj@sample$phen) > 0){
    add_phen_info <- match.arg(add_phen_info, names(sample$phen))
  } else {
    add_phen_info <- NULL
  }


  if(is.null(add_phen_info)){
    p_data <- data.frame(loc = density(sample_data$number)$x, den = density(sample_data$number)$y)
    p_data$den <- p_data$den / max(p_data$den) /2
    p_data <- subset(p_data, loc >= min(sample_data$number) & loc <= max(sample_data$number))
    p_data <- rbind(p_data,
                    c(max(sample_data$number), 0),
                    c(min(sample_data$number), 0))
    p <- ggplot() +
      geom_polygon(data = p_data, aes(x = -den , y = loc), fill = violin_color) +
      geom_jitter(data = sample_data, aes(x = .25, y = number), width = .125) +
      scale_x_continuous(limit = c(-.5, .5)) +
      theme_classic() +
        theme(axis.text.x = element_blank(),
              axis.title.x = element_blank())

  } else {
    sample_data$phen <- sample$phen[[add_phen_info]][match(sample$name, sample_data$sample)]
    phen_list <- sort(unique(sample_data$phen))
    p_data <- do.call(rbind, lapply(1:length(phen_list), function(x){
      curr_phen <- phen_list[x]
      curr_data <- subset(sample_data, phen == curr_phen)
      if(nrow(curr_data) == 1) return(NULL)
      den_data <- data.frame(loc = density(curr_data$number)$x, den = density(curr_data$number)$y)
      den_data <- subset(den_data, loc >= min(curr_data$number) & loc <= max(curr_data$number))
      data.frame(loc = c(den_data$loc, rev(range(curr_data$number))),
                 den = c(den_data$den / max(den_data$den) / 2, 0, 0),
                 phen = phen_list[x],
                 x = x)
    }))
    sample_data$x <- match(sample_data$phen, phen_list)
    p <- ggplot() +
      geom_polygon(data = p_data, aes(x = -den + x, y = loc, fill = phen)) +
      geom_jitter(data = sample_data, aes(x = x + .25, y = number, fill = phen), width = .125) +
      scale_x_continuous(breaks = 1:length(phen_list), labels = phen_list) +
      labs(x = add_phen_info) +
      theme_classic() + theme(legend.position = "none") +
      theme(axis.text.x = element_text(size = x_text_size),
            axis.title.x = element_text(size = x_title_size))

    color_phen <- get_anno_palette(phen_info_color_list, sample$phen[add_phen_info])
      p <- p + scale_fill_manual(values = color_phen[[add_phen_info]])

  }

  p + labs(y = "Gene Number") +
    theme(axis.ticks.x = element_blank(),
          axis.line.x = element_blank(),
          axis.text.y = element_text(size = y_text_size),
          axis.title.y = element_text(size = y_title_size))

}













