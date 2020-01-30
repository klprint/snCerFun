#' GO enrichment of selected genes
#'
#' @param interesting.genes Genes to test enrichment for
#' @param gene.universe All genes to test against
#' @param gene2go.list List of gene-name to GO-term
#' @param n.out Number of rows to return of the tibble
#'
#' @return Tibble with enrichment analysis, based on Fisher exact test.
#'
#' @export
run.topGO <- function(interesting.genes, gene.universe, gene2go.list, n.out = 25){
  library(topGO)
  gene.list <- as.factor(as.numeric(gene.universe %in% interesting.genes))
  names(gene.list) <- gene.universe

  sampleGO <- topGO::new("topGOdata",
                  ontology = "BP",
                  allGenes = gene.list,
                  # geneSel = names(treecut.pt)[treecut.pt == 2],
                  annot = annFUN.gene2GO,
                  gene2GO = gene2go.list)
  resultFisher <- runTest(sampleGO, algorithm = "classic", statistic = "fisher")

  return(
    GenTable(sampleGO,
             classicFisher = resultFisher,
             orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = n.out)
  )

}


#' Merge sparse matrices
#'
#' @param ... List of sparse matrices
#'
#' @return Merged sparse matrix
merge.sparse <- function(...) {

  cnnew <- character()
  rnnew <- character()
  x <- vector()
  i <- numeric()
  j <- numeric()

  icount <- 1
  cat("Merged matrices: ")
  for (M in list(...)) {
    cat(icount, " ")
    cnold <- colnames(M)
    rnold <- rownames(M)

    cnnew <- union(cnnew,cnold)
    rnnew <- union(rnnew,rnold)

    cindnew <- match(cnold,cnnew)
    rindnew <- match(rnold,rnnew)
    ind <- unname(Matrix::which(M != 0,arr.ind=T))
    i <- c(i,rindnew[ind[,1]])
    j <- c(j,cindnew[ind[,2]])
    x <- c(x,M@x)

    icount <- icount +1

  }
  cat("\n")

  sparseMatrix(i=i,j=j,x=x,dims=c(length(rnnew),length(cnnew)),dimnames=list(rnnew,cnnew))
}