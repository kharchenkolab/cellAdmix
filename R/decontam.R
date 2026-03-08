#' @import dplyr
#' @importFrom magrittr %<>% %$%
#' @importFrom igraph V V<-
NULL

#' Create and adjacency matrix
#'
#' @param df Data frame with molecule coordinates
#' @param k Neighborhood size
#' @param weighted Logical whether to create a weighted matrix
#'
#' @return sparse matrix of molecule adjacencies
#' @export
knn_adjacency_matrix <- function(df, k, weighted=FALSE) {
  nns <- df[,c("x", "y", "z")] %>% FNN::get.knn(k=k)
  xn <- nns$nn.index

  i <- rep(1:nrow(xn), each=ncol(xn))
  j <- t(xn) %>% as.vector()
  x <- if(weighted) {t(nns$nn.dist) %>% as.vector()} else rep(1, length(i))
  adj_mat <- Matrix::sparseMatrix(i=i, j=j, x=x, dims=c(nrow(xn),nrow(xn)))
  return(adj_mat)
}

#' Get KNN counts for each molecule
#' @param df Data frame with molecule coordinates
#' @param k Neighborhood size
#' @param all.genes all genes to include as columns (defualt is all genes present in df)
#' @param include_i whether to include self-molecule gene as a nearest neighbor
#'
#' @return a molecule by gene matrix of knn counts
#' @export
knn_count_matrix <- function(df, k, all.genes=NULL, include_i=FALSE) {
  if (include_i) {
    adj_mat <- knn_adjacency_matrix(df, k)
    diag(adj_mat) <- 1
  } else {
    adj_mat <- knn_adjacency_matrix(df, k)
  }

  j <- if (is.null(all.genes)) factor(df$gene) else factor(df$gene, levels=all.genes)
  ncm <- adj_mat %*% Matrix::sparseMatrix(i=1:length(j), j=as.integer(j), x=rep(1, length(j)), dims=c(nrow(df),length(levels(j))))
  colnames(ncm) <- levels(j)
  rownames(ncm) <- rownames(df)

  return(ncm)
}

#' Get KNN count matrices for molecules in all cells separately
#' @param df Data frame with molecule coordinates
#' @param h Neighborhood size
#' @param n.cores Number of cores to use
#'
#' @return Data frame with KNN counts for each molecule
#' @export
get_knn_counts_all <- function(df, h, include_i=FALSE, n.cores=1) {
  h <- as.integer(h)

  cells <- df %>% count(cell) %>% filter(n > h) %>% pull(cell)

  N <- df %>%
    filter(cell %in% cells) %>%
    split(.$cell) %>%
    parallel::mclapply(\(df_cell) {
      cell <- unique(df_cell$cell)

      counts <- knn_count_matrix(df_cell, include_i=include_i, k=h)
      counts <- df_cell %>%
        dplyr::select(cell, id = mol_id, gene_i = gene) %>%
        cbind.data.frame(as.data.frame(as.matrix(counts)))
      counts
    }, mc.cores = n.cores) %>%
    bind_rows()

  N[is.na(N)] <- 0

  return(N)
}

##### Factorization #####

#' Downsample cells to balance their representation in the dataset
#'
#' @param df Data Frame of molecule coordinates, gene, cell, and celltype columns
#' @param cell_annot Data Frame of metadata
#' @param num_cells_samp Approximately the total number of cells to sample in the dataset
#'
#' @return the df downsampled to the selected cells
#' @export
samp_ct_equal <- function(df, cell.annot, num.cells.samp) {
  ## selecting cells to have balanced numbers of molecules coming from each cell type
  ct.counts <- table(cell.annot$celltype)
  total.lib.sizes <- table(df$celltype)
  av.lib.sizes <- total.lib.sizes / ct.counts
  # now do multinomial sampling with probabilities 1/av_lib_sizes, normalized
  probs <- 1 / av.lib.sizes
  probs <- probs / sum(probs)

  all.ct.sel.counts <- round(num.cells.samp*probs)

  cells.keep <- c()
  for (ct in names(all.ct.sel.counts)) {
    cell.annot.sub2 <- cell.annot[cell.annot$celltype == ct,]
    if (nrow(cell.annot.sub2) < all.ct.sel.counts[ct]) {
      cells.select <- cell.annot.sub2$cell
    } else {
      cells.select <- sample(cell.annot.sub2$cell, all.ct.sel.counts[ct])
    }
    cells.keep <- c(cells.keep, cells.select)
  }

  df.sub.samp <- df[df$cell %in% cells.keep,]
  return(df.sub.samp)
}

get_high_exp_genes <- function(tx.dat, cell.annot, frac.thresh=.1) {
  gene.factors <- factor(tx.dat$gene)
  cell.factors <- factor(tx.dat$cell)

  # Create a sparse matrix with counts
  mycounts.sparse <- Matrix::sparseMatrix(
    i = as.integer(cell.factors),  # Row indices (cells)
    j = as.integer(gene.factors),  # Column indices (genes)
    x = rep(1, nrow(tx.dat)),      # Each row in `tx_dat` contributes 1 count
    dims = c(length(levels(cell.factors)), length(levels(gene.factors))),
    dimnames = list(
      levels(cell.factors),        # Row names (cell names)
      levels(gene.factors)         # Column names (gene names)
    )
  )
  mycounts.sparse <- Matrix::t(mycounts.sparse)

  all.ct <- unique(cell.annot$celltype)
  all.g.keep <- c()
  for (ct in all.ct) {
    cells.keep <- rownames(cell.annot)[cell.annot$celltype == ct]
    cells.keep <- intersect(cells.keep, colnames(mycounts.sparse))
    counts.sub <- mycounts.sparse[, cells.keep]
    frac.exp <- Matrix::rowSums(counts.sub > 0) / ncol(counts.sub)
    g.keep <- names(frac.exp)[frac.exp >= frac.thresh]
    all.g.keep <- c(all.g.keep, g.keep)
  }
  all.g.keep <- unique(all.g.keep)
  return(all.g.keep)
}

get_high_var_genes <- function(tx.dat, n.hvgs=1000) {
  so.spat <- get_counts_meta_seurat(tx.dat)
  so.spat <- Seurat::NormalizeData(so.spat)
  n.hvgs <- pmin(nrow(so.spat), n.hvgs) # in case n_hvgs is bigger than total n genes
  so.spat <- Seurat::FindVariableFeatures(so.spat, selection.method="vst", nfeatures=n.hvgs)
  all.g.keep <- Seurat::VariableFeatures(so.spat)
  return(all.g.keep)
}

#' Downsample the data to a subset of genes before doing NMF
#'
#' @param df Data Frame of molecule coordinates, gene, cell, and celltype columns
#' @param cell_annot Data Frame of metadata
#' @param top_g_thresh Numeric between 0-1. Only select genes with mean fractional
#' expression above this threshold in any cell type (default=NULL)
#' @param n_hvgs Numeric number of high variance genes to select (default=NULL)
#'
#' @return the df subset to the selected genes.
#' @export
subset_genes <- function(df, cell.annot, top.g.thresh=NULL, n.hvgs=NULL) {
  if (!is.null(top.g.thresh)) {
    top.g <- get_high_exp_genes(df, cell.annot, frac.thresh=top.g.thresh)
    df <- df[df$gene %in% top.g,]
  }

  if (!is.null(n.hvgs)) {
    top.g <- get_high_var_genes(df, n.hvgs=n.hvgs)
    df <- df[df$gene %in% top.g,]
  }
  df$gene_sub <- TRUE
  return(df)
}

#' Run knn for each cell and NMF using molecules from many cells
#'
#' @param df Data frame with molecule coordinates
#' @param k Number of factors (default=5)
#' @param h Neighborhood size for KNN (default=20)
#' @param mols_test Specific molecules to test (default=NULL)
#' @param nmol_dsamp Number of molecules to downsample to (default=NULL)
#' @param n.cores Number of cores to run parallel with (defualt=1)
#'
#' @return the NMF results object
#' @export
run_knn_nmf <- function(df, k=5, h=20, mols.test=NULL, nmol.dsamp=NULL, n.cores=1) {
  rownames(df) <- df$mol_id

  if (!is.null(mols.test)) {
    # restrict df to just the cells with the molecules to test
    cells.test <- unique(df[mols.test, 'cell'])
    df <- df[df$cell %in% cells.test,]
  }

  cells <- df %>% count(cell) %>% filter(n > h) %>% pull(cell)
  df <- df[df$cell %in% cells,]
  genes <- unique(df$gene) %>% sort()
  N <- get_knn_counts_all(df, h=h, include_i=TRUE, n.cores=n.cores)
  rownames(N) <- N$id
  X <- as.matrix(N[, genes])

  # limit to only mols_test
  if (!is.null(mols.test)) {
    mols.test <- mols.test[mols.test %in% rownames(X)]
    X <- X[mols.test,]
  }

  # further downsampling X
  if (!is.null(nmol.dsamp)) {
    if (nrow(X) > nmol.dsamp) {
      rows.select <- sample(1:nrow(X), nmol.dsamp, replace=FALSE)
      X <- X[rows.select,]
    }
  }

  X <- X[rowSums(X) > 0, colSums(X) > 0]
  Z <- matrix(rep(1 / colSums(X), nrow(X)), nrow=nrow(X), byrow=TRUE)
  gc()

  n.cores <- paste0('vmp', n.cores)
  res <- NMF::nmf(X, rank=k, method='ls-nmf', weight=Z, nrun=30, .options=n.cores, seed=0)
  res$gene_sub <- ('gene_sub' %in% colnames(df))

  return(res)
}

#### MRF segmentation #####

#' Get graph from cell
#' @param df_cell Data frame with molecule coordinates
#' @param h Neighborhood size
#' @param distinct_edges Whether to keep only distinct edges
#' @return Graph object
get_graph <- function(df.cell, h=10, distinct.edges=FALSE) {
  h <- as.integer(h)

  if (h >= nrow(df.cell)) {
    h <- min(h, nrow(df.cell) - 1)
    message('h is too large, setting to ', h, '...')
  }

  if (h < 2) {
    stop('Need to remove cells with only two or fewer molecules')
  }

  xn <- df.cell %>% select(x, y, z) %>% FNN::get.knn(k=h) %>% .$nn.index

  gene.dict <- df.cell %>% mutate(id=1:n()) %>% {setNames(.$gene, .$mol_id)}
  id.dict <- df.cell %>% mutate(id=1:n()) %>% {setNames(.$mol_id, .$id)}

  xn <- apply(xn, 2, function(i) {id.dict[i]})
  rownames(xn) <- df.cell$mol_id

  edges <- xn %>% as.data.frame() %>%
    tibble::rownames_to_column('x') %>%
    reshape2::melt(id.var='x', value.name='y') %>%
    select(x, y)

  G <- igraph::graph_from_data_frame(edges, directed=FALSE)

  # only keep distinct edges
  if (distinct.edges) {
    G <- igraph::simplify(G)
  }

  V(G)$gene <- gene.dict[names(V(G))]

  V(G)$id <- names(V(G))

  # V(G)$centrality = igraph::estimate_betweenness(G, directed = F, cutoff = -1)

  return(G)
}

#' Run CRF for a cell
#' @param df_cell Data frame with coordinates and gene labels for one cell
#' @param lds Pre-normalized NMF loadings
#' @param k Number of factors
#' @param h Neighborhood size
#' @param t Transition probability
#'
#' @return Data frame with MRF segmentation
run_crf <- function(df.cell, lds, k, h=10, t=0.05, max.iter=10000) {
  G <- df.cell %>% get_graph(h=h)

  vert.mat <- igraph::as_data_frame(G, 'vertices') %>%
    left_join(lds, by=join_by(gene))

  adj <- igraph::as_adjacency_matrix(G)

  crf <- CRF::make.crf(adj, k)

  crf$node.pot <- vert.mat %>%
    select(as.character(1:k)) %>%
    as.matrix

  A <- matrix(rep(t, k^2), nrow=k)
  diag(A) <- 1-(k-1)*t

  for (i in 1:length(crf$edge.pot)) {
    crf$edge.pot[[i]] <- A
  }

  out <- CRF::decode.lbp(crf, max.iter=max.iter, verbose=TRUE)

  df.out <- data.frame(id=vert.mat$name, factor=out) %>%
    mutate(cell=unique(df.cell$cell))

  return(df.out)
}

# project nmf to all genes that are in the cells that were included
# note: this may still exclude some genes if they weren't present in any included cells
# using nnls to project more gene loadings
# use an h val which is quite bigger than the NMF one to reduce loadings sparsity
project_nmf <- function(df, res, nmf.h, n.cores=1) {
  project_gene <- function(new.gene, H) {
    nnls::nnls(H, new.gene)$x
  }

  rownames(df) <- df$mol_id
  mol.res <- res@fit@W # factor scores
  mols.tested <- rownames(mol.res)

  # get cells tested from the mols tested
  cells.tested <- unique(df[mols.tested, 'cell'])
  df <- df[df$cell %in% cells.tested,]

  cells <- df %>% count(cell) %>% filter(n > nmf.h) %>% pull(cell)
  df <- df[df$cell %in% cells,]
  mols.keep <- intersect(mols.tested, rownames(df))
  mol.res <- mol.res[mols.keep,]

  df$gene <- as.factor(df$gene)
  all.genes <- unique(df$gene)

  knn.res <- df %>%
    split(.$cell, drop=FALSE) %>%
    sccore::plapply(
      function(df.c) {
        print(df.c$cell[1])
        if (nmf.h >= nrow(df.c)) {
          nmf.h <- nrow(df.c) - 1
        }
        if (nmf.h < 2) {
          stop('Need to remove cells with only two or fewer molecules')
        }

        kmat <- knn_count_matrix(df.c, k=nmf.h, all.genes=all.genes, include_i=TRUE)
        kmat <- kmat[rownames(kmat) %in% mols.keep, , drop=FALSE]
        return(kmat)
      }, mc.preschedule=TRUE, n.cores=n.cores, progress=TRUE
    )

  knn.res <- do.call(rbind, knn.res)
  knn.res <- knn.res[mols.keep,]

  knn.res <- knn.res[, sparseMatrixStats::colSums2(knn.res) > 0]
  knn.res <- t(as.matrix(knn.res))

  lds.proj <- apply(knn.res, 1, function(gene.exp) project_gene(gene.exp, mol.res))

  # lds_proj <- apply(lds_proj, 1, function(x) {
  #   x[x==0] <- min(x[x!=0])
  #   return(x)
  # })
  # lds_proj <- t(lds_proj)

  # lds_proj <- apply(lds_proj, 2, function(x) {
  #   x[x==0] <- min(x[x!=0])
  #   return(x)
  # })

  lds.proj[lds.proj == 0] <- min(lds.proj[lds.proj != 0])

  return(lds.proj)
}

#' Run crf for all cells
#'
#' @param df Data Frame with molecule coordinates, gene, cell, and mol_id columns
#' @param res NMF results object
#' @param num_nn Neighborhood size to use in knn for CRF (default=10)
#' @param same.label.ratio Transition ratio P_no_change/P_to_FactorX (default=5)
#' @param normalize.by Method to use for normalizing NMF loadings (default='gene')
#' @param n.cores Number of cores to run parallel with (default=1)
#'
#' @return A Data Frame with a column for molecule factor assignment, ordered same as df
#' @export
run_crf_all <- function(
    df, res, num.nn=10, same.label.ratio=5,
    normalize.by=c('gene', 'factor', 'gene.factor'), proj.h=NULL, n.cores=1
  ) {
  normalize.by <- match.arg(normalize.by)
  rownames(df) <- df$mol_id

  ### if only hvgs are used, then need to project all genes onto NMF to get their loadings
  if (res$gene_sub) {
    if (is.null(proj.h)) {
      proj.h <- num.nn * 2
    }
    lds <- t(project_nmf(df, res, proj.h, n.cores))
  } else {
    lds <- t(res@fit@H)
  }

  k <- ncol(lds)

  # if there are genes missing from loadings, give them the min loading
  all.g <- unique(df$gene)
  g.miss <- all.g[!(all.g %in% rownames(lds))]
  if (length(g.miss) > 0) {
    lds.add <- matrix(min(lds), nrow=length(g.miss), ncol=ncol(lds))
    rownames(lds.add) <- g.miss
    lds <- rbind(lds, lds.add)
  }

  if (normalize.by == 'gene') { # normalize by sum over all genes for a factor
    lds %<>% {t(.) / colSums(.)} %>% t()
  } else if (normalize.by == 'factor') { # normalize by sum over all factors for a gene
    lds %<>% {. / rowSums(.)}
  } else if (normalize.by == 'gene.factor') {
    lds %<>% {. / rowSums(.)} %>% {t(.) / colSums(.)} %>% t()
  }
  lds <- cbind.data.frame(rownames(lds), lds)
  colnames(lds)[1] <- 'gene'

  # set the transition probabilities based on a ratio for (no transition) / (transition to f_i)
  tval <- 1 / (same.label.ratio + k - 1)

  max.iter <- 200
  crf.res <- df %>%
    split(.$cell, drop=TRUE) %>%
    sccore::plapply(
      function(df.c) {
        run_crf(df.c, lds, k=k, h=num.nn, max.iter=max.iter, t=tval)
      }, mc.preschedule=TRUE, n.cores=n.cores, progress=TRUE, fail.on.error=TRUE
    ) %>%
    bind_rows()

  # reordering to match cell order of df
  rownames(crf.res) <- crf.res$id
  crf.res <- crf.res[rownames(df),]
  crf.res <- crf.res[, 'factor', drop=FALSE]

  return(crf.res)
}













