
setClass("PAV",
         slots = list(
           pav_data = "matrix",
           args = "list",
           gene = "list",
           sample = "list"
         )
)


#' get_pav_obj
#'
#' Get an object of the PAV class.
#'
#' @param pav_data  A numeric `matrix` or `data.frame` of PAV data. `0` represent "absence" and `1` represent "presence".
#' The row names are gene names. The column names are sample names.
#' @param phen_info The phenotype data. The row names are sample name.
#' @param gene_info The gene data. The row names are gene name.
#' @param add_softcore A logical value indicating whether to consider `soft-core` when determining the gene types.
#' @param add_private A logical value indicating whether to consider `private` when determining the gene types.
#' @param softcore_loss_rate A numeric vector of loss rate.
#' @param softcore_p_value A numeric vector of p-value (binomial test).
#'
#' `add_softcore` and `add_private` are logical values indicating whether to consider "softcore" and "private" when determining the gene types.
#'
#' If `add_softcore` is `TRUE`, the genes with loss rates not significantly larger than `softcore_loss_rate` will be considered as "softcore" gene. Binomial tests (with a null hypothesis of loss rate < `softcore_loss_rate`) are carried out for each gene. A p-vlaue below `softcore_p_value` means that this gene is lost in a significant proportion and is a distributed gene(loss rate > `softcore_loss_rate`).
#'
#' If `add_private` is `TRUE`, the genes present in only one sample will be condidered as "private" gene.
#'
#' @export

get_pav_obj <- function(pav_data,
                        phen_info = NULL,
                        gene_info = NULL,
                        add_softcore = T,
                        add_private = T,
                        softcore_loss_rate = .1,
                        softcore_p_value = .05){

  pav_data <- check_pav_cov_data("pav_data", pav_data)
  pav_data <- pav_data[order(rowSums(pav_data), decreasing = T), order(colSums(pav_data), decreasing = T)]

  if(!is.null(gene_info)) gene_info <- check_info(gene_info, rownames(pav_data), "gene", "gene")
  if(!is.null(phen_info)) phen_info <- check_info(phen_info, colnames(pav_data), "phen", "sample")

  add_softcore <- match_logi("add_softcore", add_softcore)
  add_private <- match_logi("add_private", add_private)
  softcore_loss_rate <- match_num("softcore_loss_rate", softcore_loss_rate, 0, 1)
  softcore_p_value <- match_num("softcore_p_value", softcore_p_value, 0, 1)


  type_info_res <- get_type_info(ncol(pav_data), add_private, add_softcore, softcore_loss_rate, softcore_p_value)
  type_info <- type_info_res$type_info
  dist_sample_n <- type_info_res$dist_sample_n

  type <- get_gene_type(pav_data, type_info)
  n <- as.integer(as.vector(rowSums(pav_data)))

  obj <- new("PAV")

  obj@pav_data <- as.matrix(pav_data, dimnames = NULL)

  obj@args <- list(
    add_softcore = add_softcore,
    add_private = add_private,
    softcore_loss_rate = softcore_loss_rate,
    softcore_p_value = softcore_p_value,
    dist_n = dist_sample_n
  )

  obj@gene <- list(
    name = rownames(pav_data),
    sample_n = n,
    type = type,
    info = as.list(gene_info)
  )

  obj@sample <- list(
    name = colnames(pav_data),
    phen = as.list(phen_info)
  )

  return(obj)

}


