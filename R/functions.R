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

#' Creates a pseudobulk using a grouped tibble
#'
#' @param umi UMI count matrix (genes x cells)
#' @param df.grouped Grouped tibble. The grouping should only be based on a single column
#' @param cell.col  Column in the tibble, where the cell names are stored.
#' @param min.cells Minimal number of cells which should be merged to a pseudobulk
#' @param n.cores Number of cores to run in parallel.
#'
#' @return Matrix of pseudobulks
get.pb.umi <- function(umi, df.grouped, cell.col = "cell_id",  min.cells = 50, n.cores = 1){
  df.list <- df.grouped %>% filter(n() >= min.cells) %>% group_split()
  names(df.list) <- df.grouped %>% filter(n() >= min.cells) %>% group_keys() %>% pull(1)

  df.list <- df.list

  pbmcapply::pbmclapply(df.list, function(x){
    cells <- x %>% pull(cell.col)
    rowSums(umi[ ,cells ])
  }, mc.cores = n.cores) %>%
    do.call(cbind, .)
}