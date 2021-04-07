
check_class <- function(obj, class_name){
  if(class(obj) != class_name){
    stop("Please input a `", class_name, "` object")
  }
}



check_pav_cov_data <- function(data_name, data){
  if(!is.data.frame(data) && !is.matrix(data)){
    stop("The `", data_name, "` should be a data frame or a matrix")
  } else if(!all(apply(data, 2, is.numeric))){
    stop("The `", data_name, "` should only contain numeric.")
  } else if(any(rowSums(data) == 0)){
    stop("Existing gene in `", data_name, "` that is 0 in all samples. Please remove.")
  } else if(any(colSums(data) == 0)){
    stop("Existing sample in `", data_name, "` that is 0 in all genes. Please remove.")
  } else {
    data
  }
}


check_info <- function(data, names, info_name, unit_name){
  if(!is.data.frame(data)){
    stop("The `", info_name, "_info` should be a data frame.")
  }
  if(!all(names %in% rownames(data))){
    err_gene <- setdiff(names, rownames(data))
    if(length(err_gene) <= 5){
      err_gene_message <- paste0(err_gene, collapse = ",")
    } else {
      err_gene_message <- paste0(paste0(err_gene[1:5], collapse = ","), "...")
    }
    stop("Missing ", unit_name, ": '", err_gene_message, "' in ", info_name ,"_info.")
  } else {
    data[names, , drop = F]
  }
}


match_logi <- function(arg_name, arg){
  if(is.logical(arg)){
    arg
  } else {
    stop(paste0(arg_name, " should be a logical type."))
  }
}

match_num <- function(arg_name, arg, min = -Inf, max = Inf){
  if(is.numeric(arg)){
    if(arg > min && arg < max){
      arg
    } else {
      stop(paste0(arg_name, "should less than ", max, " and greater than ", min, "."))
    }
  } else {
    stop(paste0(arg_name, " should be a numerical type."))
  }
}



dist_n <- function(n, prob, p_value){
  for(i in 1:n) {
    binom_res <- binom.test(x = i, n = n, p = prob, alternative = 'greater')
    if(binom_res$p.value <= p_value){
      res <- i
      break
    }
  }
  return(res-1)
}



get_type_info <- function(sample_n, private, softcore, prob, p_value){
  if(softcore){
    sample_dist_n <- dist_n(sample_n, prob, p_value)
    type_info <- data.frame(
      sample_n = 1:sample_n,
      gene_type = c(rep("distributed",sample_n - 1 - sample_dist_n),
                    rep("softcore",sample_dist_n),"core")
    )
    dist_sample_n <- c(1, sample_n - 1 - sample_dist_n)
  } else {
    type_info <- data.frame(
      sample_n = 1:sample_n,
      gene_type = c(rep("distributed", sample_n - 1), "core")
    )
    dist_sample_n <- c(1, sample_n - 1)
  }
  if(private){
    type_info[1, 2] <- "private"
    dist_sample_n[1] <- 2
  }

  return(list(dist_sample_n = dist_sample_n, type_info = type_info))
}




get_gene_type <- function(pav, type_info){
  type_data <- merge(
    data.frame(
      gene = rownames(pav),
      col_sum = rowSums(pav)),
    type_info,
    by.x = "col_sum", all.x = TRUE,
    by.y = colnames(type_info)[1],
    sort = TRUE )
  type_data[match(rownames(pav), type_data$gene), 3]
}




merge_args <- function(def_args, args){
  c(args[names(args) %in% names(def_args)], def_args[!names(def_args) %in% names(args)])
}

get_type_palette <- function(colors_type){
  type_col_def <- structure(c("#7CAE00", "#00BFC4", "#F8766D", "#C77CFF"),
                            names = c("core", "softcore", "distributed", "private"))
  if(is.vector(colors_type) && !is.na(colors_type)) {
    if(is.null(names(colors_type))){
      warning("`colors_type` should be a named vector.")
      type_col_def
    } else {
      merge_args(type_col_def, colors_type)
    }
  } else {
    type_col_def
  }
}

get_palette <- function(data_list){
  if(length(data_list) == 0) return(NULL)

  ## continue
  data_list_continue <- data_list[unlist(lapply(data_list, is.numeric))]
  color_n <- length(data_list_continue)
  if(color_n > 0){
    colors <- scales::hue_pal()(color_n)
      res_continue <- lapply(1:color_n, function(x){
        if(length(unique(data_list_continue[[x]])) == 1){
          structure(colors[x], names = unique(data_list_continue[[x]]))
        } else {
          circlize::colorRamp2(breaks = c(max(data_list_continue[[x]], na.rm = T),
                                          min(data_list_continue[[x]], na.rm = T)),
                               colors = c(colors[x], "white"))
        }
      })
    names(res_continue) <- names(data_list_continue)
  } else {
    res_continue <- NULL
  }

  ## discrete
  data_list_discrete <- data_list[unlist(lapply(data_list, is.character))]
  if(length(data_list_discrete) > 0){
    named_colorN <- unlist(lapply(data_list_discrete, function(x){length(unique(x))}))
    colors <- randomcoloR::distinctColorPalette(sum(named_colorN))
    Ns <- c(0, cumsum(as.vector(named_colorN)))
    res_discrete <- lapply(1:length(named_colorN), function(x){
      cur_col <- colors[(Ns[x]+1):Ns[x+1]]
      names(cur_col) <- unique(data_list_discrete[[x]])
      cur_col
    })
    names(res_discrete) <- names(named_colorN)
  } else {
    res_discrete <- NULL
  }

  return(c(res_continue, res_discrete))
}

get_anno_palette <- function(input_colors, data_list){
  def_colors <- get_palette(data_list)
  res <- lapply(names(def_colors), function(x){
    x_def <- def_colors[[x]]
    x_input <- input_colors[[x]]
    x_data <- data_list[[x]]
    if(is.null(x_input)){
      x_def
    } else {
      if(is.numeric(x_data)){
        if(length(x_input) < 2){
          stop(x," should have at least two colors.")
        }
        circlize::colorRamp2(breaks = seq(min(x_data, na.rm = T),
                                          max(x_data, na.rm = T),
                                          length.out = length(x_input)),
                             colors = x_input)
      } else {
        if(!is.null(names(x_input))){
          merge_args(x_def, x_input)
        } else {
          x_def
        }
      }
    }
  })
  names(res) <- names(def_colors)
  res
}



is.color <- function(colors){
  all(unlist(lapply(colors, function(x){
    x %in% colors() || setequal(grep("#[0-9A-F]{6}", x), 1)
  })))
}


get_color <- function(input_color, gene_groups){
  def_colors <- structure(scales::hue_pal()(length(gene_groups)), names = gene_groups)
  if(is.null(input_color)){
    colors <- def_colors
  } else if(length(input_color) == 1 && is.color(input_color) && is.null(names(input_color))){
    colors <- structure(
      rep(input_color, length(gene_groups)),
      names = gene_groups
    )
  } else if(is.vector(input_color) && !is.null(names(input_color)) && is.color(input_color) ) {
    colors <- merge_args(def_colors, input_color)
  } else{
    stop("color should be one color or named vector")
  }
  return(colors)
}


get_legend <- function(color_info, info_data, legend_title_size, legend_text_size, legend_grid_size){

  if(length(color_info) > 0){
    lg_info <- lapply(1:length(color_info), function(x){
      if(is.numeric(info_data[[names(color_info)[x]]]) && length(unique(info_data[[names(color_info)[x]]])) > 1){
        ComplexHeatmap::Legend(col_fun = color_info[[x]],
                               title = names(color_info)[x],
                               title_gp = grid::gpar(fontsize = legend_title_size,
                                                     fontface = "bold"),
                               labels_gp = grid::gpar(fontsize = legend_text_size),
                               grid_height = legend_grid_size,
                               grid_width = legend_grid_size)
      } else {
        ComplexHeatmap::Legend(labels = names(color_info[[x]]),
                               legend_gp = grid::gpar(fill = color_info[[x]]),
                               title = names(color_info)[[x]],
                               title_gp = grid::gpar(fontsize = legend_title_size,
                                                     fontface = "bold"),
                               labels_gp = grid::gpar(fontsize = legend_text_size),
                               grid_height = legend_grid_size,
                               grid_width = legend_grid_size)
      }

    })
  } else {
    lg_info <- NULL
  }
}

get_anno_row <- function(data_phen, color_phen, row_anno_phen){
  eval(parse(text = paste0(
    "anno <- ComplexHeatmap::rowAnnotation(",
    paste0(lapply(names(data_phen), function(x){
      paste0(x, " = ComplexHeatmap::anno_simple(data_phen[['", x, "']], border = row_anno_phen$border,
             width = grid::unit(row_anno_phen$width, 'mm'),
             col = structure(color_phen[['", x, "']], names = names(color_phen[['", x, "']])))")
    }), collapse = ","),
    ", annotation_name_side = row_anno_phen$name_side, annotation_name_rot = row_anno_phen$name_rot,
    annotation_name_gp = grid::gpar(fontsize = row_anno_phen$name_size)    )"
    )))
  anno
}


get_anno_column <- function(data_info, color_info, column_anno_gene){
  eval(parse(text = paste0(
    "anno <- ComplexHeatmap::HeatmapAnnotation(",
    paste0(lapply(names(data_info), function(x){
      paste0(x, " = ComplexHeatmap::anno_simple(data_info[['", x, "']], border = column_anno_gene$border,
               height = grid::unit(column_anno_gene$height,'mm'),
               col = structure(color_info[['",x,"']], names = names(color_info[['",x,"']])))")
    }), collapse = ","),
    ", annotation_name_side = column_anno_gene$name_side, annotation_name_rot = column_anno_gene$name_rot,
      annotation_name_gp = grid::gpar(fontsize = column_anno_gene$name_size) )"
  )))
  anno
}


check_phen_stat_res <- function(data){
  if(!all(colnames(data) == c("phen", "gene", "p_value", "p_adjusted"))){
    stop("please input result from function `phen_stat`.")
  }
}


check_obj_phen <- function(pav_obj){
  if(length(pav_obj@sample$phen) == 0)
    stop("can't find `phen_info` in pav_obj.")
}

check_obj_gene <- function(pav_obj){
  if(length(pav_obj@gene$info) == 0)
    stop("can't find `gene_info` in pav_obj.")
}




