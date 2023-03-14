#' Single cell RNA-sequencing (scRNA-seq) training subset of the 10X Genomics MALT counts.
#'
#' A random subset of joint scRNA-seq/CITE-seq 10X Genomics human extranodal marginal zone B-cell tumor/mucosa-associated lymphoid tissue (MALT) training data. See the \code{dataProcessing.R} file from \code{data-raw} folder for code to recreate data subset.
#'
#' @format A scRNA-seq counts matrix of \code{dgCMatrix-class} from the \code{Matrix} package with 1000 cells and 33538 genes.
#' @source <https://www.10xgenomics.com/resources/datasets/10-k-cells-from-a-malt-tumor-gene-expression-and-cell-surface-protein-3-standard-3-0-0>
"train.malt.rna.mat"
#'
#'Cellular Indexing of Transcriptomes and Epitopes by Sequencing (CITE-seq) training subset of the 10X Genomics MALT counts.
#'
#' A random subset of joint scRNA-seq/CITE-seq 10X Genomics human extranodal marginal zone B-cell tumor/mucosa-associated lymphoid tissue (MALT) training data. See the \code{dataProcessing.R} file from \code{data-raw} folder for code to recreate data subset.
#'
#' @format A CITE-seq counts matrix of \code{dgeMatrix-class} from the \code{Matrix} package with 1000 cells and 17 genes.
#' @source <https://www.10xgenomics.com/resources/datasets/10-k-cells-from-a-malt-tumor-gene-expression-and-cell-surface-protein-3-standard-3-0-0>
"train.malt.adt.mat"
#'
#' Single cell RNA-sequencing (scRNA-seq) target subset of the 10X Genomics MALT counts.
#'
#' A random subset of joint scRNA-seq/CITE-seq 10X Genomics human extranodal marginal zone B-cell tumor/mucosa-associated lymphoid tissue (MALT) target data. See the \code{dataProcessing.R} file from \code{data-raw} folder for code to recreate data subset.
#'
#' @format A scRNA-seq counts matrix of \code{dgCMatrix-class} from the \code{Matrix} package with 4000 cells and 33538 genes.
#' @source <https://www.10xgenomics.com/resources/datasets/10-k-cells-from-a-malt-tumor-gene-expression-and-cell-surface-protein-3-standard-3-0-0>
"target.malt.rna.mat"
