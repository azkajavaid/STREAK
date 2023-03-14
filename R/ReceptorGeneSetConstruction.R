#' Gene sets weights membership matrix construction for receptor abundance estimation.
#'
#' Computes \eqn{n x h} gene sets weights membership matrix using associations learned between log-normalized and reduced rank reconstructed (RRR)
#' \eqn{m x n} scRNA-seq training data and \eqn{m x h} CITE-seq ADT training counts normalized using the centered log ratio (CLR) transformation. scRNA-seq counts are normalized and
#' RRR using the \code{\link[SPECK:randomizedRRR]{SPECK::randomizedRRR()}} function while CITE-seq counts are normalized using the
#' \code{\link[Seurat:NormalizeData]{Seurat::NormalizeData()}} function with the \code{normalization.method} parameter set to \code{CLR}. Spearman rank correlations are computed between the
#' normalized CITE-seq data and the normalized and RRR scRNA-seq data.
#'
#' @param train.rnaseq \eqn{m x n} scRNA-seq counts matrix for \eqn{m} cells
#' and \eqn{n} genes.
#' @param train.citeseq \eqn{m x h} CITE-seq ADT counts matrix for \eqn{m} cells (same cells as the \code{train.rnaseq matrix}) and \eqn{h} cell-surface proteins.
#' @param rank.range.end See documentation for the \code{randomizedRRR} function from the \code{SPECK} package.
#' @param min.consec.diff See documentation for the \code{randomizedRRR} function from the \code{SPECK} package.
#' @param rep.consec.diff See documentation for the \code{randomizedRRR} function from the \code{SPECK} package.
#' @param manual.rank See documentation for the \code{randomizedRRR} function from the \code{SPECK} package.
#' @param seed.rsvd See documentation for the \code{randomizedRRR} function from the \code{SPECK} package.
#' @return
#' \itemize{
#'   \item \code{receptor.geneset.matrix} - A \eqn{n x h} gene sets weights membership matrix where a column \eqn{i} from \eqn{h} corresponds to the weights for \eqn{n} genes from the scRNA-seq matrix trained against the corresponding CITE-seq ADT transcript \eqn{h}.
#' }
#'
#' @examples
#' data("train.malt.rna.mat")
#' data("train.malt.adt.mat")
#' receptor.geneset.matrix.out <- receptorGeneSetConstruction(train.rnaseq =
#'                                          train.malt.rna.mat[1:100,1:80],
#'                                          train.citeseq =
#'                                          train.malt.adt.mat[1:100,1:2],
#'                                          rank.range.end = 70,
#'                                          min.consec.diff = 0.01,
#'                                          rep.consec.diff = 2,
#'                                          manual.rank = NULL, seed.rsvd = 1)
#' dim(receptor.geneset.matrix.out)
#' head(receptor.geneset.matrix.out)
#' @export
receptorGeneSetConstruction <- function(train.rnaseq, train.citeseq, rank.range.end = 100, min.consec.diff = 0.01, rep.consec.diff = 2, manual.rank = NULL, seed.rsvd = 1) {
  if (missing(train.rnaseq)) {
    stop("scRNA-seq training expression matrix is missing.")
  }
  if (missing(train.citeseq)) {
    stop("CITE-seq training expression matrix is missing.")
  }
  message("Performing normalization and reduced rank reconstruction on training matrix.\n")
  counts.matrix.rrr.train <- randomizedRRR(counts.matrix = train.rnaseq, rank.range.end = rank.range.end, min.consec.diff = min.consec.diff, rep.consec.diff = rep.consec.diff, manual.rank = manual.rank, seed.rsvd = seed.rsvd)
  train.rnaseq.norm.rrr <- counts.matrix.rrr.train$rrr.mat
  train.citeseq.object <- CreateSeuratObject(counts = t(as.matrix(train.citeseq)), assay = "ADT")
  train.citeseq.norm <- NormalizeData(train.citeseq.object, normalization.method = "CLR", margin = 2)
  train.citeseq.norm.mat <- t(as.matrix(train.citeseq.norm@assays$ADT@data))
  message("Performing coexpression analysis.\n")
  receptor.geneset.matrix <- matrix(nrow = ncol(train.rnaseq.norm.rrr), ncol = ncol(train.citeseq.norm.mat))
  for (i in 1:ncol(train.citeseq.norm.mat)) {
    for(j in 1:ncol(train.rnaseq.norm.rrr)) {
      receptor.geneset.matrix[j,i] <- cor(as.numeric(train.citeseq.norm.mat[,i]),
                                      as.numeric(train.rnaseq.norm.rrr[,j]),
                                      method="spearman")
    }
    message("\rComputing co-expression for receptor ", i, " of ", ncol(train.citeseq.norm.mat))
  }
  colnames(receptor.geneset.matrix) <- colnames(train.citeseq)
  rownames(receptor.geneset.matrix) <- colnames(train.rnaseq)
  return(receptor.geneset.matrix)
}
