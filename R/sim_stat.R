
#' sim_stat
#'
#' Perform genome simulation calculations.
#'
#' @param pav_obj A PAV object.
#' @param genome_type A vector of genome types. These can be any of the following: "pan", "core" and "private".
#' @param n_simulation The number of simulations.
#' @param parallel A logical value indicating whether to use parallel computing.
#' @param parallel_n The number of CPU cores used for parallel computing.
#'
#' @import snowfall
#' @importFrom parallel detectCores
#' @importFrom data.table data.table
#'
#' @export


sim_stat <- function(pav_obj,
                     genome_type = c("pan", "core"),
                       n_simulation = 10,
                       parallel = FALSE,
                       parallel_n = parallel::detectCores() - 1
){

  check_class(pav_obj, "PAV")

  pav_data <- t(pav_obj@pav_data)
  if(nrow(pav_data) > 0){
    if(all(apply(pav_data, 2, function(x){ all(x %in% c(0,1)) })) ) {
      if(sum(pav_data) == nrow(pav_data) * ncol(pav_data)){
        stop("the input is all 1.")
      }else{
        all_data <- data.table::data.table(pav_data, keep.rownames = T)
      }
    }else{
      stop("the input should just include 0 and 1.")
    }
  }else{
    stop("`pav_data` don't have rows.")
  }

  if("private" %in% unique(pav_obj@gene$type)) {
    groups <- match.arg(genome_type, c("pan", "core", "private"), several.ok = T)
  } else {
    groups <- match.arg(genome_type, c("pan", "core"), several.ok = T)
  }

  if(!(is.numeric(n_simulation) & n_simulation > 1)){
    stop("'n_simulation' should be upper than 1.")
  }

  if(is.logical(parallel)){
    if(parallel){
      if(is.numeric(parallel_n) & (parallel_n <= parallel::detectCores())){
        res_data <- sim_n_pallel(all_data, groups, n_simulation, parallel_n)
      }else{
        stop("`parallel_n` should be less than CPU cores.")
      }
    }else{
      res_data <- sim_n(all_data, groups, n_simulation)
    }

  }else{
    stop("'parallel' should be a logical value.")
  }

  return(res_data)
}





cal_num <- function(data, groups, core_n){
  col_sum <- colSums(data[,-1,])
  stat_res <- sapply(groups, function(x){
    if(x == "pan"){
      sum(col_sum > 0)
    } else if( x == "core"){
      sum(col_sum == nrow(data))
    } else if( x == "private"){
      sum(col_sum == 1) - core_n
    }
  })
  return(stat_res)
}

sim_n <- function(data,groups,n){
  n_core_all <- colnames(data[,-1,])[colSums(data[,-1,])==nrow(data[,-1,])]
  data <- data[,!..n_core_all,]
  res <- data.frame(matrix(unlist(
    lapply(1:nrow(data),function(x){
      unlist(lapply(1:n,function(y){
        c(x,cal_num(data[sample(nrow(data),x),,], groups, length(n_core_all)) + length(n_core_all))
      }))
    })
  ),ncol=length(groups)+1,byrow = T))
  colnames(res) <- c("Times",groups)
  return(res)
}


sim_n_pallel <- function(data, groups, n, ncore){
  require(snowfall)
  snowfall::sfInit(parallel = TRUE, cpus = ncore )
  snowfall::sfExport("cal_num")
  snowfall::sfExport("data")
  snowfall::sfExport("n")
  n_core_all <- colnames(data[,-1,])[colSums(data[,-1,])==nrow(data[,-1,])]
  data <- data[,!..n_core_all,]
  res <- do.call(rbind,
                 snowfall::sfLapply(1:nrow(data),function(x){
                   t(replicate(n,c(x,cal_num(data[sort(sample(nrow(data),x)),,],groups, length(n_core_all)) + length(n_core_all))))
                 }))
  snowfall::sfStop()
  colnames(res) <- c("Times", groups)
  data.table::data.table(res)
}


