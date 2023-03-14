#' Receptor abundance estimation for single cell RNA-sequencing (scRNA-seq) data using gene set scoring and thresholding.
#'
#' Performs receptor abundance estimation for \eqn{m x n} scRNA-seq target data using gene set scoring and thresholding. scRNA-seq target counts are normalized and
#' reduced rank reconstructed (RRR) using the \code{\link[SPECK:randomizedRRR]{SPECK::randomizedRRR()}} function. Gene set scoring is next performed leveraging expression from
#' the top most weighted genes based on the gene sets weights membership matrix with the \code{\link[VAM:vam]{VAM::vam()}} function. The resulting cell-specific gene set scores
#' are then thresholded utilizing the \code{\link[Ckmeans.1d.dp:Ckmeans.1d.dp]{Ckmeans.1d.dp::Ckmeans.1d.dp()}} function. Note that this function only performs normalization
#' and does not perform any quality control (QC) checks on the inputted target scRNA-seq counts matrix. Any QC needed can be performed on the target matrix before passing it
#' as an input to the function.
#'
#' @param target.rnaseq \eqn{m x n} scRNA-seq counts matrix for \eqn{m} cells and \eqn{n} genes.
#' @param receptor.geneset.matrix \eqn{n x h} Gene sets weights membership matrix.
#' @param num.genes Number of top most weighted genes for subsequent gene set scoring and thresholding.
#' @param rank.range.end See documentation for the \code{randomizedRRR} function from the \code{SPECK} package.
#' @param min.consec.diff See documentation for the \code{randomizedRRR} function from the \code{SPECK} package.
#' @param rep.consec.diff See documentation for the \code{randomizedRRR} function from the \code{SPECK} package.
#' @param manual.rank See documentation for the \code{randomizedRRR} function from the \code{SPECK} package.
#' @param seed.rsvd See documentation for the \code{randomizedRRR} function from the \code{SPECK} package.
#' @param max.num.clusters See documentation for the \code{ckmeansThreshold} function from the \code{SPECK} package.
#' @param seed.ckmeans See documentation for the \code{ckmeansThreshold} function from the \code{SPECK} package.
#' @return
#' \itemize{
#'   \item \code{receptor.abundance.estimates} - A \eqn{m x h} matrix consisting of abundance estimates for \eqn{m} cells in \eqn{h} receptors.
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
#'                                          manual.rank = NULL,
#'                                          seed.rsvd = 1)
#' dim(receptor.geneset.matrix.out)
#' head(receptor.geneset.matrix.out)
#' data("target.malt.rna.mat")
#' receptor.abundance.estimates.out <- receptorAbundanceEstimation(target.rnaseq =
#'                                     target.malt.rna.mat[1:200,1:80],
#'                                     receptor.geneset.matrix =
#'                                     receptor.geneset.matrix.out,
#'                                     num.genes = 10, rank.range.end = 70,
#'                                     min.consec.diff = 0.01,
#'                                     rep.consec.diff = 2,
#'                                     manual.rank = NULL, seed.rsvd = 1,
#'                                     max.num.clusters = 4, seed.ckmeans = 2)
#' dim(receptor.abundance.estimates.out)
#' head(receptor.abundance.estimates.out)
#' @export
receptorAbundanceEstimation <- function(target.rnaseq, receptor.geneset.matrix, num.genes = 10, rank.range.end = 100, min.consec.diff = 0.01, rep.consec.diff = 2, manual.rank = NULL, seed.rsvd = 1, max.num.clusters = 4, seed.ckmeans = 2) {
  if (missing(target.rnaseq)) {
    stop("scRNA-seq target expression matrix is missing.")
  }
  if (missing(receptor.geneset.matrix)) {
    stop("Gene sets weights membership matrix is missing.")
  }
  message("Performing normalization and reduced rank reconstruction on target matrix.\n")
  counts.matrix.rrr.target <- randomizedRRR(counts.matrix = target.rnaseq, rank.range.end = rank.range.end, min.consec.diff = min.consec.diff, rep.consec.diff = rep.consec.diff, manual.rank = manual.rank, seed.rsvd = seed.rsvd)
  target.rnaseq.norm.rrr <- counts.matrix.rrr.target$rrr.mat
  message("Performing gene set scoring and thresholding.\n")
  target.rnaseq.norm.rrr <- target.rnaseq.norm.rrr[,which(colMeans(target.rnaseq.norm.rrr) != 0)]
  receptor.geneset.matrix <- receptor.geneset.matrix[which(rownames(receptor.geneset.matrix) %in% colnames(target.rnaseq.norm.rrr)),]
  receptor.abundance.estimates <- matrix(nrow = nrow(target.rnaseq.norm.rrr), ncol = ncol(receptor.geneset.matrix))
  for (i in 1:ncol(receptor.geneset.matrix)) {
    top.genes <- names(sort(receptor.geneset.matrix[,i], decreasing = TRUE))[1:num.genes]
    top.weights <- as.numeric(sort(receptor.geneset.matrix[,i], decreasing = TRUE))[1:num.genes]
    if (length(top.genes[!is.na(top.genes)]) == 0) {
      message(paste("No genes present in the gene sets weights membership matrix for ", colnames(receptor.geneset.matrix)[i], ".", "\n", sep=""))
      receptor.abundance.estimates[,i] <- NA
    } else {
      target.rnaseq.subset <- target.rnaseq.norm.rrr[,top.genes]
      vam.output <- vam(gene.expr = target.rnaseq.subset, gene.weights = top.weights, gamma = T, center = F)
      vam.cdf.values <- vam.output$cdf.value
      vam.sq.distances <- vam.output$distance.sq
      set.seed(seed.ckmeans)
      ckmeans.out <- suppressWarnings(Ckmeans.1d.dp(x = vam.sq.distances, k = c(1:max.num.clusters)))
      num.cluster <-  ckmeans.out$cluster
      val.centers <- ckmeans.out$centers
      if (length(val.centers) > 1) {
        min.val <- which(val.centers == min(val.centers))
        min.num <- which(num.cluster == min.val)
        vam.cdf.values[min.num] <- 0
      }
      receptor.abundance.estimates[,i] <- vam.cdf.values
    }
  }
  colnames(receptor.abundance.estimates) <- colnames(receptor.geneset.matrix)
  rownames(receptor.abundance.estimates) <- rownames(target.rnaseq)
  return(receptor.abundance.estimates)
}
