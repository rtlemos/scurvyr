#' Compute a Minimum Spanning Tree given a sparse similarity matrix
#'
#' @param similarity_matrix Sparse similarity matrix (= 1/distance_matrix)
#' @param verbose Print status to console
#' @param verbose_every Console notification of progress
#'
#' @return Sparse MST with 1 for connected nodes and 0 elsewhere
#' @export
#'
#' @examples
#' minimumSpanningTree(sparseMatrix(i = c(1,4,1,2,3,2,3,4,3,4),
#'                                  j = c(1,1,4,2,2,3,3,3,4,4),
#'                                  x = c(Inf,1,1,Inf,3,3,Inf,2,2,Inf)), 
#'                                  verbose = FALSE)
#'
minimumSpanningTree = function(similarity_matrix, verbose, verbose_every = 100){
  x = enforceM(similarity_matrix)
  small_number = 0
  nr = nrow(x)
  if (verbose) cat('Building Minimum Spanning Tree for', nr, 'points... ')
  tree = 1
  alt_tree = 1
  diag(x) = small_number
  maxmat = getColMax(x, 1:nr)
  for (i in 1:(nr - 1)) {
    m = maxmat[tree, 1]
    index.j = tree[which.max(m)[1]]
    index.i = maxmat[index.j,2]
    crit = x[index.i, tree] != 0
    pos = tree[crit]
    for (k in pos) {
      x[index.i, k] = small_number
      x[k, index.i] = small_number
    }
    tree = c(tree, index.i)
    alt_tree = c(alt_tree, index.j)
    maxmat[tree,] = getColMax(x, tree)
    if (verbose && i %% verbose_every == 0) cat(paste0(round(100 * i / nr, 1), '% '))
  }
  tree = tree[2:nr]
  alt_tree = alt_tree[2:nr]
  w = sparseMatrix(i=c(tree, alt_tree),
                   j=c(alt_tree, tree),
                   x=rep(1, 2 * (nr - 1)), dims = c(nr, nr))
  dimnames(w) = dimnames(x)
  if (verbose) cat('done.\n')
  w
}

#' Force a user-provided matrix to be of type sparseMatrix
#'
#' @param x User matrix
#'
#' @return sparseMatrix
#' @export
#'
#' @examples
#' enforceM(matrix(runif(6),2,3))
#'
enforceM = function(x) {
  tp = attr(class(x), "package")
  if (is.null(tp) || tp != "Matrix") {
    m = as.matrix(x)
    nr = nrow(m)
    nc = ncol(m)
    xx = unlist(lapply(1:nc, FUN = function(i) m[, i]))
    i = rep(1:nr,nc)
    j = unlist(lapply(1:nc, FUN = function(i) rep(i, nr)))
    crit = (xx != 0)
    i = i[crit]
    j = j[crit]
    xx = xx[crit]
    x = sparseMatrix(i=i, j=j, x=xx)
    x
  }
  return(x)
}

#' Fetch one column of a sparse matrix
#'
#' @param x sparse matrix
#' @param col_id column index
#'
#' @return List with non-empty rows and non-zero values
#' @export
#'
#' @examples
#' getColumn(sparseMatrix(i=1,j=2,x=9), 2)
#'
getColumn = function(x, col_id) {
  st = x@p[col_id] + 1
  en = x@p[col_id + 1]
  out = if (en < st) {
    list(i=NULL, x=NULL)
  } else {
    ii = x@i[st:en] + 1
    xx = x@x[st:en]
    list(i = ii, x = xx)
  }
  out
}


#' Compute max and which.max for columns of a sparse matrix
#'
#' @param x sparse matrix
#' @param cols column indices
#'
#' @return dense matrix with two columns: max and which.max
#' @export
#'
#' @examples
#' getColMax(sparseMatrix(i=1,j=2,x=9), 1:2)
#'
getColMax = function(x, cols) {
  mat = t(mapply(cols, FUN = function(j) {
    colinfo = getColumn(x, j)
    if (length(colinfo$x) > 0) {
      c(max(colinfo$x), colinfo$i[which.max(colinfo$x)])
    } else {
      c(0, j)
    }
  }))
  mat
}
