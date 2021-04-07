


setClass("COV",
         slots = list(
           cov_data = "matrix",
           gene = "list",
           sample = "list"
         )
)



#' get_cov_obj
#'
#' Get an object of the COV class.
#'
#' @param cov_data A numeric `matrix` or `data.frame` of the coverage data.
#' The row names are gene names. The column names are sample names.
#' @param phen_info A `data.frame` of phenotype data. The row names are sample name.
#' @param gene_info A `data.frame` of gene data. The row names are gene name.
#'
#' The `cov_data` shouldn't contain gene that is 0 in all sample, and the same for sample.
#'
#'
#' @export

get_cov_obj <- function(cov_data,
                        phen_info = NULL,
                        gene_info = NULL){

  cov_data <- check_pav_cov_data("cov_data", cov_data)

  if(!is.null(gene_info)) gene_info <- check_info(gene_info, rownames(cov_data), "gene", "gene")
  if(!is.null(phen_info)) phen_info <- check_info(phen_info, colnames(cov_data), "phen", "sample")

  obj <- new("COV")

  obj@cov_data <- as.matrix(cov_data)
  obj@gene <- as.list(gene_info)
  obj@sample <- as.list(phen_info)

  obj
}



