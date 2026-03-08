### make df into a counts matrix and a metadata matrix per cell

#' Transform transcript molecule df into standard seurat object with counts matrix
#'
#' @param tx.dat dataframe of transcripts with at least gene, cell, and celltype columns
#'
#' @return Seurat object with the data in counts matrix form
#' @export
get_counts_meta_seurat <- function(tx.dat, normalize=FALSE) {
  gene.factors <- factor(tx.dat$gene)
  cell.factors <- factor(tx.dat$cell)

  # Aggregate duplicate (gene, cell) pairs — required for Matrix >= 1.5
  idx_agg <- data.frame(
    i = as.integer(gene.factors),
    j = as.integer(cell.factors)
  ) %>%
    dplyr::group_by(i, j) %>%
    dplyr::summarise(x = dplyr::n(), .groups = "drop")

  cm <- Matrix::sparseMatrix(
    i = idx_agg$i,
    j = idx_agg$j,
    x = idx_agg$x,
    dims = c(length(levels(gene.factors)), length(levels(cell.factors))),
    dimnames = list(
      levels(gene.factors),
      levels(cell.factors)
    )
  )

  meta <- tx.dat[,c("cell", "celltype")] %>%
    data.table::as.data.table() %>%
    unique() %>%
    as.data.frame()

  colnames(meta)[2] <- 'cell_type'
  rownames(meta) <- meta$cell
  meta <- meta[colnames(cm),]

  s.obj <- Seurat::CreateSeuratObject(cm, meta.data=meta)

  if (normalize) {
    s.obj <- Seurat::NormalizeData(s.obj)
  }
  return(s.obj)
}


