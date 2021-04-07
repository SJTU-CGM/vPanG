


#' sim_multi_stat
#'
#' Perform genome simulation calculations for groups.
#'
#' @param pav_obj A PAV object.
#' @param phen_name The name of phenotype used for grouping.
#' @param n_simulation The number of simulations.
#' @param parallel A logical value indicating whether to use parallel computing.
#' @param parallel_n The number of CPU cores used for parallel computing.
#'
#' @export

sim_multi_stat <- function(pav_obj,
                           phen_name,
                           n_simulation = 10,
                           parallel = FALSE,
                           parallel_n = parallel::detectCores() - 1){

  check_class(pav_obj, "PAV")

  sample_list <- pav_obj@sample
  pav_data <- pav_obj@pav_data
  args <- pav_obj@args
  show_phen <- match.arg(phen_name, names(sample_list$phen))

  do.call(rbind,
          lapply(unique(sample_list$phen[[show_phen]]), function(group){
              curr_samples <- sample_list$name[sample_list$phen[[show_phen]] == group]
              if(length(curr_samples) == 1){
                warning(paste0(group, "just have one sample, so ignore."))
                return(NULL)
              }
              curr_pav_data <- pav_data[, curr_samples, drop = F]
              curr_pav_data <- curr_pav_data[rowSums(curr_pav_data) != 0, ]
              sim_data <- sim_stat(get_pav_obj(curr_pav_data, add_softcore = args$add_softcore, add_private = args$add_private,
                                       softcore_p_value = args$softcore_p_value, softcore_loss_rate = args$softcore_loss_rate),
                                   n_simulation = n_simulation, parallel = parallel, parallel_n = parallel_n)
              sim_data_l <- data.table::melt(data.table::data.table(sim_data), id.vars = "Times")
              sim_data_stat <- sim_data_l[, .(Length = mean(value), SD = sd(value)), by = c("variable", "Times")]
              res_data <- cbind(sim_data_stat, group)
          })
  )
}



#' sim_multi_plot
#'
#' Plot the result of genome simulation.
#'
#' @param sim_multi_data The result from \code{\link[vPan]{sim_multi_stat}}.
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
#' @param legend_side The position of legend.
#' @param legend_title_size The size of legend title.
#' @param legend_text_size The size of legend item labels.
#'
#' @importFrom ggplot2 expansion
#'
#' @export

sim_multi_plot <- function(sim_multi_data,
                           path_size = 1,
                           path_color = NULL,
                           ribbon_fill = NULL,
                           ribbon_alpha = 0.3,

                           x_title = "Times",
                           y_title = "Count",
                           x_title_size = NULL,
                           y_title_size = NULL,
                           x_breaks = NULL,
                           y_breaks = NULL,
                           x_text_size = NULL,
                           y_text_size = NULL,

                           legend_side = "right",
                           legend_title_size = NULL,
                           legend_text_size = NULL){

  if(!all(colnames(sim_multi_data) == c("variable", "Times", "Length", "SD", "group"))){
    stop("please input result from function `sim_multi_stat`.")
  }

  p_data <- sim_multi_data
  gene_groups <- unique(sim_multi_data$group)
  colors <- get_color(path_color, gene_groups)
  fills <- get_color(ribbon_fill, gene_groups)

  p <- ggplot(subset(p_data), aes(x = Times, y = Length)) +
    geom_ribbon(aes(ymin = Length - SD, ymax = Length + SD, fill = group, group = paste0(group, variable)), alpha = ribbon_alpha) +
    geom_path(aes(color = group, group = paste0(group, variable), linetype = variable), size = path_size) +
    scale_x_continuous(expand = expansion(mult = c(0, .05))) +
    scale_fill_manual(values = fills) +
    scale_color_manual(values = colors) +
    labs(x = x_title, y = y_title, linetype = "Genome") +
    theme_classic() +
    theme(
      axis.title.x = element_text(size = x_title_size),
      axis.title.y = element_text(size = y_title_size),
      axis.text.x = element_text(size = x_text_size),
      axis.text.y = element_text(size = y_text_size),
      legend.position = legend_side,
      legend.title = element_text(size = legend_title_size),
      legend.text = element_text(size = legend_text_size)
    )
  if(!is.null(x_breaks)){
    p <- p + scale_x_continuous(breaks = x_breaks)
  }

  if(!is.null(y_breaks)){
    p <- p + scale_y_continuous(breaks = y_breaks)
  }
  print(p)
}


