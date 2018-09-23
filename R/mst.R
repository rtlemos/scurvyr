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
    if (i %% verbose_every == 0) cat(paste0(round(100 * i / nr, 1), '% '))
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


#' Plot a path
#'
#' @param dataset List with data, path, lat, lon, etc.
#' @param plot_data Plot the data?
#' @param lat_bounds Latitude bounds (for zooming in)
#' @param lon_bounds Longitude bounds (for zooming in)
#' @param colored_line Color the path using a rainbow palette?
#' @param fixed_aspect_ratio Use a 1:1 aspect ratio in the plot?
#' @param background Background color
#'
#' @return ggplot of path
#' @export
#'
plotPath = function(dataset, plot_data = FALSE,
                    lat_bounds = NULL, lon_bounds = NULL,
                    colored_line = TRUE, fixed_aspect_ratio = TRUE, background = 'black') {
  path = dataset$path
  lat = dataset$lat
  lon = dataset$lon
  data = dataset$data
  restricted_bounds = !is.null(lat_bounds) & !is.null(lon_bounds)
  n = length(path)
  nr = length(lat)
  nc = length(lon)
  path_col = floor((path - 1) / nr) + 1
  path_row = (path - 1) %% nr + 1
  ri = path_row[1:(n - 1)]
  rj = path_row[2:n]
  ci = path_col[1:(n - 1)]
  cj = path_col[2:n]

  if (plot_data) {
    rst = cbind(expand.grid(lat=lat,lon=lon), z=as.numeric(data))
    if (restricted_bounds) {
      crit = rst$lat >= lat_bounds[1] & rst$lat <= lat_bounds[2] &
        rst$lon >= lon_bounds[1] & rst$lon <= lon_bounds[2]
      rst = rst[crit, ]
    }
    p = ggplot() + geom_raster(rst, mapping = aes(x=lon, y=lat, fill=z)) +
      scale_fill_gradient(low='dark gray', high='white')
  } else {
    p = ggplot()
  }

  df = data.frame(lat = lat[c(rbind(ri, rj))], lon = lon[c(rbind(ci, cj))], id = c(rbind(1:(n - 1), 1:(n - 1))))
  if (restricted_bounds) {
    crit = df$lat >= lat_bounds[1] & df$lat <= lat_bounds[2] &
      df$lon >= lon_bounds[1] & df$lon <= lon_bounds[2]
    df = df[crit,]
  }
  if (fixed_aspect_ratio) {
    p = p + coord_equal()
  }
  if (colored_line) {
    p = p + geom_line(data = df, mapping=aes(x=lon, y=lat, group=id, color = id)) +
      scale_color_gradientn(colours=rainbow(100))
  } else {
    p = p + geom_line(data=df, mapping=aes(x=lon, y=lat, group=id))
  }
  if (restricted_bounds) {
    p = p + scale_x_continuous(limits = lon_bounds, expand = c(0, 0)) +
      scale_y_continuous(limits = lat_bounds, expand = c(0, 0))
  }
  p = p + guides(color=FALSE, fill=FALSE) +
    theme(panel.background = element_rect(fill = background),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  p

}
