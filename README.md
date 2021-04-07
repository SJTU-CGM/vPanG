
# vPanG

## 1 Introduction
 
### 1.1 Pan-genome and PAV analysis

&emsp;&emsp;Pan-genome is the total genes found across members of a given population, revealing the diversity and functional potential within that population. PAV(Presence/absence variation) analysis is an essential step in pangenome analysis. The core genome contains genes shared by all individuals and the distributed genome consist of genes not shared by all. The distributed genome can be further divided into genes shared in most members (soft-core genes), genes shared between some members (distributed/accessory genes), and genes present in only one member (unique/private genes).

&emsp;&emsp;“Map-to-pan” is a common strategy for eukaryotic pan-genome analyses. A pan-genome is first constructed by integrating the reference genome and non-reference sequences. Then, reads are aligned to the pan-genome and the percentage of gene region and coding region covered by read alignment are examined. Finally, gene presence/absence variations (PAVs) are determined and consequent analysis is carried out based on this resulted PAV table.


### 1.2 The core processing workflow of vPan

&emsp;&emsp;The tool vPanG is efficient to explore and visualize the complex results in PAV analysis. It provides four modules to 1) display gene coverage distribution, 2) analyze and visualize PAV table, 3) estimate pan/core genome sizes by simulation, and 4) find phenotype-associated genes.

&emsp;&emsp;The workflow starts with an object of COV(module 1) or PAV(module 2-4) class. It can be produced by function `get_cov_obj()` and `get_pav_obj()`. It contains coverage/PAV matrix, arguments, phenotype information, and gene information, and will be the main input for the subsequent steps. 

* Module 1 focuses on showing the gene coverage and the distribution of gene coverage among individuals. The functions in module 1 are `cov_heatmap()` and `cov_density()`.

* Module 2 focuses on PAV analysis, including overview, classifying genes, clustering, PCA analysis. The functions for PAV analysis are `pav_heatmap()`, `pav_hist()`, `pav_stackbar()`, `pav_cluster()`, `pav_halfviolin()` and `pav_pca()`.

* Module 3 focuses on pan/core/private genome size estimation by simulation and drawing growth curve. The functions in module 3 are `sim_stat()`, `sim_plot()`, `sim_multi_stat()` and `sim_multi_plot`.

* Module 4 focuses on phenotype association and visualization. The functions in module 4 include `phen_heatmap()`, `phen_manhattan()`, `phen_block()`, `phen_bar()` and `phen_violin()`.


## 2 Installation

### 2.1 Installing R/RStudio
&emsp;&emsp;In order to run vPanG, we need the R software environment, which can be freely downloaded as follows:

* Install [R](https://www.r-project.org/)
* Install [RStudio](https://www.rstudio.com/)


### 2.2 Check or install packages

```
packages <- c("snowfall", "data.table", "ggdendro", "ggplot2", "ggrepel", "ggsignif", "randomcoloR")
lapply(packages, function(x) {
	if(!require(x, character.only = TRUE)) {
		install.packages(x, dependencies = TRUE)
	}})
if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}
BiocManager::install("ComplexHeatmap")

```

### 2.3 Install metaFunc from github.

```{r eval=FALSE}
if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")
library(devtools)
install_github("SJTU-CGM/vPanG", build_vignettes = TRUE)
```

## 3 User manual

&emsp;&emsp;You can use `vignette("vPan")` in R to view the vignette of vPan. It can help you start using vPan, including the format and examples of input data, and the use of functions. You can use `?function` to view the help documentation of the function, such as `?pav_heatmap`. 



