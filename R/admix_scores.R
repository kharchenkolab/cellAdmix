### ME Admixture Scoring

#' @export
estimate_cell_adjacency <- function(df.spatial, random.shift=1e-3, edge.max.mad=4, n.cores=1) {
  required.cols <- c('x', 'y', 'celltype', 'cell')
  if (!all(required.cols %in% colnames(df.spatial))) {
    stop(paste0('Required columns: ', paste(required.cols, collapse=', ')))
  }

  df.spatial %<>% filter(!is.na(celltype))

  df.spatial$x %<>% {. + runif(length(.), -random.shift, random.shift)}
  df.spatial$y %<>% {. + runif(length(.), -random.shift, random.shift)}

  n.grid.cells <- if (n.cores > 1) ceiling(sqrt(n.cores) * 2) else 1
  df.spatial$grid_id <- as.integer(cut(df.spatial$x, n.grid.cells)) * n.grid.cells +
    as.integer(cut(df.spatial$y, n.grid.cells))

  grid.dfs <- split(df.spatial, ~grid_id) %>% .[sapply(., nrow) > 3]

  adj.df <- sccore::plapply(grid.dfs, \(cdf) {
    triangles <- cdf[,c('x', 'y')] %>% as.matrix() %>% geometry::delaunayn()
    edges <- triangles %>% {rbind(.[,c(1, 2)], .[,c(1, 3)], .[,c(2, 3)])} %>%
      {cbind(pmin(.[,1], .[,2]), pmax(.[,1], .[,2]))} %>% unique()

    # Molecule-based edges
    adj.df <- data.frame(is=edges[,1], ie=edges[,2]) %>% mutate(
      xs=cdf$x[is], xe=cdf$x[ie], ys=cdf$y[is], ye=cdf$y[ie],
      cts=cdf$celltype[is], cte=cdf$celltype[ie],
      cell_s=cdf$cell[is], cell_e=cdf$cell[ie]
    ) %>% select(-is, -ie)

    adj.df$dist <- (adj.df[,c('xs', 'ys')] - adj.df[,c('xe', 'ye')])^2 %>% rowSums() %>% sqrt()
    adj.df
  }, n.cores=n.cores, mc.preschedule=TRUE, fail.on.error=TRUE, progress=TRUE) %>% bind_rows()

  adj.df %<>% filter(dist < (median(dist) + edge.max.mad * mad(dist)))

  # Cell-cell edges
  adj.df %<>%
    filter(cell_s != cell_e) %>%
    select(cell_s, cell_e, cts, cte) %>%
    unique()

  # Add reverse edges
  adj.df %<>% select(cell_e, cell_s, cte, cts) %>% set_colnames(colnames(adj.df)) %>% rbind(adj.df) %>% unique()

  return(adj.df)
}

#' @export
estimate_cell_type_adjacency <- function(cell.adj.df) {
    adj.mat <- cell.adj.df %>%
        group_by(cell_s, cts, cte) %>% dplyr::count() %>% group_by(cts, cte) %>% dplyr::summarise(n=mean(n)) %>%
        tidyr::pivot_wider(names_from=cts, values_from=n) %>%
        tibble::column_to_rownames('cte') %>% as.matrix() %>%
        .[colnames(.), colnames(.)]

    adj.mat[is.na(adj.mat)] <- 0

    return(adj.mat)
}

#' @export
estimate_gene_prob_per_type <- function(cm, cell.annot, expr.thres=0, use.counts=FALSE) {
  cids.per.type <- cell.annot %>% {split(names(.), .)}
  if (use.counts) {
    probs <- cids.per.type %>% sapply(\(cids) Matrix::rowSums(cm[,cids,drop=FALSE])) %>% {. / colSums(.)}
  } else {
    probs <- cids.per.type %>% sapply(\(cids) Matrix::rowMeans(cm[,cids,drop=FALSE] > expr.thres))
  }

  return(t(probs))
}

#' Estimate gene scores
#'
#' @description Estimates terms of the Bayesian formula for  \eqn{p(m_i \in C | g_i, ...)}
#'
#' @param prob.gene.per.type Gene probabilities per cell type \eqn{p(g_i | t_x, M)}
#' @param cell.type.adj.mat Cell type adjacency matrix \eqn{K}
#' @param cell.type Current cell type \eqn{u}
#' @return Data.frame with columns `cont` = \eqn{\sum_{v \neq u} \left(p(g_i | t_v, M) \cdot \frac{k_{v,u}}{\sum_{x \neq u} k_{x, u}} \right)} and `ref`=\eqn{p(g_i | t_u, M)}
estimate_gene_scores <- function(prob.gene.per.type, cell.type.adj.mat, cell.type) {
  cont.fracs <- cell.type.adj.mat[,cell.type] %>% .[names(.) != cell.type] %>% {. / sum(.)}
  prob.gene.per.type.adj <- cont.fracs %>% {prob.gene.per.type[names(.),] * .}
  scores.per.gene.cont <- prob.gene.per.type.adj %>% .[rownames(.) != cell.type,] %>% colSums()
  scores.per.gene.ref <- prob.gene.per.type[cell.type, names(scores.per.gene.cont)]

  return(data.frame(cont=scores.per.gene.cont, ref=scores.per.gene.ref))
}

#' Estimate gene contamination probabilities
#'
#' @description Estimates \eqn{p(m_i \in C | g_i, ...)}
#' @param gene.scores Data.frame with columns `cont` = \eqn{\sum_{v \neq u} \left(p(g_i | t_v, M) \cdot \frac{k_{v,u}}{\sum_{x \neq u} k_{x, u}} \right)} and `ref`=\eqn{p(g_i | t_u, M)}
#' @param p.c Proportion of cells in the contamination cluster \eqn{p_c}. If null, forces \eqn{p_{c,u} = 0.5}
#' @param prior.cont.prob Prior probability of contamination \eqn{\frac{\sum_{x \neq u}k_{x,u}}{\sum_{x}k_{x,u}}}
estimate_gene_contamination_probabilities <- function(gene.scores, p.c=NULL, prior.cont.prob=NULL) {
  if (is.null(p.c)) {
    p.cu <- 0.5
  } else {
    p.cu <- p.c * prior.cont.prob
  }

  gene.cont.probs <- gene.scores %$%
      setNames(p.cu * cont / (p.cu * cont + (1 - p.cu) * ref), rownames(.))

  gene.cont.probs[is.nan(gene.cont.probs)] <- 0

  return(gene.cont.probs)
}

#' Estimate cell contamination fractions
#'
#' @description Estimates \eqn{E(f_c | ...)}
estimate_cell_admixture_fractions <- function(gene.cont.probs, cell.type.cm, mask=NULL) {
  if (!is.null(mask)) {
    gene.cont.probs[!mask[names(gene.cont.probs)]] <- 0
  }

  fracs <- cell.type.cm %>% {Matrix::colSums(. * gene.cont.probs[rownames(.)]) / Matrix::colSums(.)}

  return(fracs)
}

mask_out_genes <- function(cont.probs, genes, all.name='all') {
  if (is.null(genes)) {
    return(cont.probs)
  }

  if (!is.list(genes)) {
    cont.probs[intersect(rownames(cont.probs), genes),] <- 0
    return(cont.probs)
  }

  genes %<>% lapply(intersect, rownames(cont.probs)) %>% .[sapply(., length) > 0]

  for (n in intersect(names(genes), colnames(cont.probs))) {
    if (n == all.name) {
      cont.probs[genes[[n]],] <- 0
    } else {
      cont.probs[genes[[n]], n] <- 0
    }
  }

  return(cont.probs)
}

#' Title
#'
#' @param cm.rna counts matrix for scRNA-seq data
#' @param cm.spatial counts matrix for spatial data
#' @param annot.rna cell-type annotation for sc
#' @param annot.spatial cell-type annotation for spatial
#' @param cell.type.adj.mat ctype adjacency matrix
#' @param p.c prior probability of contamination
#' @param signal.thres threshold for selecting admixture genes
#' @param min.expr.frac min expression fraction
#' @param max.prob.high max prob
#' @param exclude.cell.types cell types to exclude
#' @param exclude.genes genes to exclude
#' @param use.counts whether to use counts (default=FALSE)
#' @param adjust whether to adjust (default=TRUE)
#'
#' @return contamination scores
#' @export
estimate_contamination_scores <- function(
    cm.rna, cm.spatial, annot.rna, annot.spatial, cell.type.adj.mat,
    p.c=0.25, signal.thres=0.25, min.expr.frac=0.05, max.prob.high=0.1,
    exclude.cell.types=NULL, exclude.genes=NULL, use.counts=FALSE, adjust=TRUE
  ) {
  common.genes <- intersect(rownames(cm.rna), rownames(cm.spatial))
  if (length(common.genes) != nrow(cm.rna) || length(common.genes) != nrow(cm.spatial)) {
    warning('Some genes not found in both cm.rna and cm.spatial, ', length(common.genes), ' genes left')
  }

  cm.rna <- cm.rna[common.genes,,drop=FALSE]
  cm.spatial <- cm.spatial[common.genes,,drop=FALSE]

  prob.gene.per.type <- estimate_gene_prob_per_type(cm.rna, annot.rna, use.counts=use.counts)
  prior.cont.probs <- cell.type.adj.mat %>% {1 - diag(.) / colSums(.)}
  expr.frac.mask <- (matrixStats::colMaxs(prob.gene.per.type) > min.expr.frac)

  prob.gene.per.type.high <- estimate_gene_prob_per_type(cm.rna, annot.rna, expr.thres=1)

  cont.probs.per.type <- colnames(cell.type.adj.mat) %>% setNames(., .) %>% setdiff(exclude.cell.types) %>%
    sapply(\(ct) {
      cont.probs <- prob.gene.per.type %>%
        estimate_gene_scores(cell.type.adj.mat, cell.type=ct) %>%
        estimate_gene_contamination_probabilities(p.c=p.c, prior.cont.prob=prior.cont.probs[ct])

      cont.probs[!((cont.probs > signal.thres) & expr.frac.mask[names(cont.probs)])] <- 0

      return(cont.probs)
    })

  cont.probs.per.type %<>% mask_out_genes(exclude.genes)

  cont.fracs.null <- colnames(cont.probs.per.type) %>% setNames(., .) %>% sapply(\(ct) {
    prob.gene.per.type[ct,] %>% {sum(cont.probs.per.type[names(.), ct] * .) / sum(.)}
  })

  doublet.scores <- colnames(cont.probs.per.type) %>% lapply(\(ct) {
    cm <- cm.spatial[,annot.spatial == ct, drop=FALSE]

    cont.probs <- cont.probs.per.type[,ct]
    cont.genes <- names(cont.probs)[cont.probs > signal.thres]
    if (any(prob.gene.per.type.high[ct, cont.genes] > max.prob.high)) {
      warning(
        'Genes with high expression in ', ct, ' detected in contamination: ',
        paste0("'", cont.genes[prob.gene.per.type.high[ct, cont.genes] > max.prob.high], "'", collapse=', '),
        '. Perhaps you should improve the annotation or mask these genes out (`exclude.genes`).'
      )
    }

    scores <- estimate_cell_admixture_fractions(cont.probs, cm)
    if (adjust) {
      scores <- pmax(scores - cont.fracs.null[ct], 0)
    }

    scores
  }) %>% unname() %>% unlist()

  return(list(
    cont.probs.per.type=cont.probs.per.type, doublet.scores=doublet.scores
    # prob.gene.per.type=prob.gene.per.type, prior.cont.probs=prior.cont.probs, expr.frac.mask=expr.frac.mask
  ))
}

#' Estimate contamination scores from Seurat objects
#'
#' @inheritDotParams estimate_contamination_scores
#' @export
estimate_contamination_scores_seurat <- function(so.rna, so.spatial, cell.type.adj.mat, idents='celltype', ...) {
  if (!is.null(idents)) {
    Seurat::Idents(so.rna) <- idents
    Seurat::Idents(so.spatial) <- idents
  }

  return(estimate_contamination_scores(
    cm.rna=so.rna[['RNA']]$counts, cm.spatial=so.spatial[['RNA']]$counts,
    annot.rna=Idents(so.rna), annot.spatial=Idents(so.spatial),
    cell.type.adj.mat=cell.type.adj.mat, ...
  ))
}

#' @export
sparse_cor <- function(x){
  # https://stackoverflow.com/questions/5888287/running-cor-or-any-variant-over-a-sparse-matrix-in-r
  n <- nrow(x)

  cmean.x <- Matrix::colMeans(x)
  covmat <- (Matrix::crossprod(x)-2*(cmean.x %o% Matrix::colSums(x)) + n*cmean.x %o% cmean.x)/(n-1)

  sdvec <- sqrt(Matrix::diag(covmat)) # standard deviations of columns
  cor.mat <- covmat / crossprod(t(sdvec)) # correlation matrix
  return(cor.mat)
}

#' @export
estimate_correlation_preservation <- function(cm1, cm2) {
  common.genes <- intersect(rownames(cm1), rownames(cm2)) %>%
    .[(Matrix::rowSums(cm1[.,]) > 0) & (Matrix::rowSums(cm2[.,]) > 0)]
  cors1 = cm1[common.genes,] %>% Matrix::t() %>% sparse_cor() %>% as.matrix()
  cors2 = cm2[common.genes,] %>% Matrix::t() %>% sparse_cor() %>% as.matrix()

  diag(cors1) = NA
  diag(cors2) = NA

  cor.cors = sapply(1:ncol(cors1), \(i) {
    cor(cors1[,i], cors2[,i], use="complete.obs")
  }) %>% setNames(common.genes)

  return(cor.cors)
}
