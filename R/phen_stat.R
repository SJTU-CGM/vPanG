

#' phen_stat
#'
#' Perform phenotypic association calculations
#'
#' @param pav_obj A PAV object.
#' @param phen_name The name of phenotype used for calculation.
#' @param p_adjust_method The adjustment methods, pass to  \code{\link[stats]{p.adjust}}
#' @param pallel A logical value indicating whether to use parallel computing.
#' @param parallel_n The number of CPU cores used for parallel computing.
#'
#'
#' @export

phen_stat <- function(pav_obj,
                      phen_name,
                      p_adjust_method = "fdr",
                      pallel = F,
                      parallel_n = parallel::detectCores() - 1){

  check_class(pav_obj, "PAV")

  if(length(pav_obj@sample$phen) == 0) {
    stop(pav_obj, " doesn't have `phen_info`.")
  }
  phen_data <- data.frame(pav_obj@sample$phen)
  phen_name <- match.arg(phen_name, colnames(phen_data), several.ok = T)
  phen_data <- phen_data[, phen_name, drop = F]

  gene_data <- data.frame(gene = pav_obj@gene$name, type = pav_obj@gene$type)
  pav_data <- pav_obj@pav_data[subset(gene_data, type != "core")$gene,]

  genes <- rownames(pav_data)
  phens <- colnames(phen_data)
  phen_res <- do.call(
    rbind,
    lapply(phens, function(phen){
      if(pallel){
        snowfall::sfInit(parallel = TRUE, cpus = parallel_n )
        snowfall::sfExport("phen_stat_p")
        snowfall::sfExport("genes")
        snowfall::sfExport("phen_data")
        snowfall::sfExport("phen")
        snowfall::sfExport("pav_data")
        res <- snowfall::sfLapply(genes, function(gene){
          p_value <- phen_stat_p(phen_data[, phen], pav_data[gene, ])
          c(phen, gene, p_value)
        })
        snowfall::sfStop()
        res
      } else {
        res <- lapply(genes, function(gene){
          p_value <- phen_stat_p(phen_data[, phen], pav_data[gene, ])
          c(phen, gene, p_value)
        })
      }
      data <- data.frame(t(as.data.frame(res)))
      rownames(data) <- NULL
      colnames(data) <- c("phen", "gene", "p_value")
      data$p_value <- as.numeric(data$p_value)
      data$p_adjusted <- p.adjust(data$p_value, method = p_adjust_method)
      data
    })
  )
}






phen_stat_p <- function(cur_phen_data, cur_gene_data){
  if(is.numeric(cur_phen_data)){
    p_value <- tryCatch({
      presence_data <- cur_phen_data[cur_gene_data == 1]
      absence_data <- cur_phen_data[cur_gene_data == 0]
      wilcox_res <- wilcox.test(presence_data, absence_data)
      wilcox_res$p.value
    }, error = function(e) {
      NA
    })
  } else {
    p_value <- tryCatch({
      fisher_res <- fisher.test(cur_phen_data, cur_gene_data, simulate.p.value = F)
      fisher_res$p.value
    }, error = function(e) {
      NA
    })
  }
  p_value
}




