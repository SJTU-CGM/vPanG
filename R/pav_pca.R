

#' pav_pca
#'
#' `pav_pca()` will perform PCA analysis of PAV data using `prcomp()`. The `center`, `scale`, and `rank` will pass to `prcomp()`.
#'
#' @param pav_obj A PAV object.
#' @param center A logical value indicating whether the variables should be shifted to be
#' zero centered, pass to \code{\link[stats]{prcomp}}.
#' @param scale A logical value indicating whether the variables should be scaled to have
#' unit variance before the analysis takes place, pass to \code{\link[stats]{prcomp}}.
#' @param rank A number specifying the maximal rank, pass to \code{\link[stats]{prcomp}}.
#'
#' @param add_phen_info A character string of `phen_info` names.
#' @param phen_info_color_list A list contains named vector of colors for `phen_info` annotation.
#' e.g. list(gender = c("Male" = "green", "Female" = "red"))
#'
#' @param axis_text_size The size of tick labels on axis.
#' @param axis_title_size The size of axis title.
#'
#' @param legend_side The position of legend ("top", "bottom", "right", "left").
#' @param legend_text_size The size of legend item labels.
#' @param legend_title_size The size of legend title.
#'
#' @importFrom ggplot2 ggplot geom_point labs aes xlab ylab theme_bw theme scale_color_manual
#'
#' @export

pav_pca <- function(pav_obj,
                    center = T,
                    scale = F,
                    rank = NULL,

                    add_phen_info = NULL,
                    phen_info_color_list = NULL,

                    axis_text_size = NULL,
                    axis_title_size = NULL,

                    legend_side = "right",
                    legend_text_size = NULL,
                    legend_title_size = NULL){

  check_class(pav_obj, "PAV")

  pav_data <- pav_obj@pav_data
  sample <- pav_obj@sample

  pca_res <- prcomp(t(pav_data), center = center, scale = scale, rank = rank)

  pca.var.per <- round(pca_res$sdev^2/sum(pca_res$sdev^2)*100, 2)

  p_data <- data.frame(sample = rownames(pca_res$x), PC1 = pca_res$x[,1], PC2 = pca_res$x[,2])

  if(length(sample$phen) == 0 || is.null(add_phen_info)){
    p <- ggplot(p_data, aes(x = PC1, y = PC2)) + geom_point()
    phen_col <- NULL
  } else {
    add_phen_info <- match.arg(add_phen_info, names(sample$phen))
    phen_col <- get_anno_palette(phen_info_color_list, sample$phen[add_phen_info])
    phen_data <- data.frame(sample = sample$name,
                            phen = sample$phen[[add_phen_info]])
    p_data <- merge(p_data, phen_data, by = "sample")

    p <- ggplot(p_data, aes(x = PC1, y = PC2, color = phen)) +
      geom_point() + labs(color = add_phen_info)
  }

  p <- p +
    xlab(paste("PC1(",pca.var.per[1],"%",")",sep=""))+
    ylab(paste("PC2(",pca.var.per[2],"%",")",sep=""))+
    theme_bw() +
    theme(axis.text = element_text(size = axis_text_size),
          axis.title = element_text(size = axis_title_size),
          legend.position = legend_side,
          legend.text = element_text(size = legend_text_size),
          legend.title = element_text(size = legend_title_size))

  if(!is.null(phen_col)){
    p <- p +
      scale_color_manual(values = phen_col[[add_phen_info]])
  }

  print(p)
}








