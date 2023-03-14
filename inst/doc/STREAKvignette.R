## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----message=FALSE, warning=FALSE---------------------------------------------
library(STREAK)

## ----message=FALSE, warning=FALSE---------------------------------------------
data("train.malt.rna.mat")
data("train.malt.adt.mat")
receptor.geneset.matrix.out <- receptorGeneSetConstruction(train.rnaseq = 
                                                  train.malt.rna.mat, 
                                                train.citeseq = 
                                                  train.malt.adt.mat[,1:5], 
                                                rank.range.end = 100, 
                                                min.consec.diff = 0.01, 
                                                rep.consec.diff = 2,
                                                manual.rank = NULL, 
                                                seed.rsvd = 1)
dim(receptor.geneset.matrix.out)
head(receptor.geneset.matrix.out)

## ----message=FALSE, warning=FALSE---------------------------------------------
data("target.malt.rna.mat")
receptor.abundance.estimates.out <- 
  receptorAbundanceEstimation(target.rnaseq = target.malt.rna.mat,
                              receptor.geneset.matrix = 
                                receptor.geneset.matrix.out,
                              num.genes = 10, rank.range.end = 100, 
                              min.consec.diff = 0.01, rep.consec.diff = 2,
                              manual.rank = NULL, seed.rsvd = 1, 
                              max.num.clusters = 4, seed.ckmeans = 2)
dim(receptor.abundance.estimates.out)
head(receptor.abundance.estimates.out)

