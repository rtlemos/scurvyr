#' Compute a square similarity matrix based on a raw rectangular matrix
#' with gridded information and possibly NAs
#'
#' @param data NxM matrix with gridded data
#' @param neighbors neighborhood structure of dual graph
#'
#' @return (NxM) x (NxM) sparse similarity matrix
#' @export
#'
#' @examples
#' n = 6
#' elevation = matrix(1, n, n)
#' elevation[1:2,1:2] = NA
#' elev_coarse = getCoarseData(elevation, n:1, 1:n, FALSE, verbose=TRUE)
#' elev_ngb = getNeighbors(elev_coarse$data, verbose=TRUE)
#' getSimilarityMatrix(elevation, elev_ngb, TRUE)
#'
getSimilarityMatrix = function(data, neighbors, verbose) {
  if (nrow(data) %% 2 == 1 || ncol(data) %% 2 == 1) {
    print('Data must have even number of rows and columns')
    stop()
  }
  if (verbose) cat('Building similarity matrix on dual graph... ')

  i = unlist(mapply(neighbors$oldid, neighbors$ngb, SIMPLIFY = FALSE, FUN = function(id, ngb) {
    rep(id, length(ngb))
  }))
  j = unlist(neighbors$ngb)

  nr = nrow(data) / 2
  ci = floor((i - 1) / nr) + 1
  ri = ((i - 1) %% nr) + 1
  cj = floor((j - 1) / nr) + 1
  rj = (j - 1) %% nr + 1

  distance = mapply(ri, ci, rj, cj, FUN = function(rri, cci, rrj, ccj) {
    d = if (cci != ccj) {
      c = min(cci, ccj)
      mean(c(abs(data[(rri - 1) * 2 + 1, c * 2] - data[(rri - 1) * 2 + 1, c * 2 + 1]),
             abs(data[rri * 2, c * 2] - data[rri * 2, c * 2 + 1])), na.rm = TRUE) -
        mean(c(abs(data[(rri - 1) * 2 + 1, c * 2] - data[rri * 2, c * 2]),
               abs(data[(rri - 1) * 2 + 1, c * 2 + 1] - data[rri * 2, c * 2 + 1])),
             na.rm=TRUE)
    } else if (rri != rrj) {
      r = min(rri, rrj)
      mean(c(abs(data[r * 2, (cci - 1) * 2 + 1] - data[r * 2 + 1, (cci - 1) * 2 + 1]),
             abs(data[r * 2, (ccj - 1) * 2 + 1] - data[r * 2 + 1, (ccj - 1) * 2 + 1])), na.rm = TRUE) -
        mean(c(abs(data[r * 2, (cci - 1) * 2 + 1] - data[r * 2, (ccj - 1) * 2 + 1]),
               abs(data[r * 2 + 1, (cci - 1) * 2 + 1] - data[r * 2 + 1, (ccj - 1) * 2 + 1])), na.rm = TRUE)
    } else 0
    max(0, d)
  })
  similarity = sparseMatrix(i = neighbors$newid[i], j = neighbors$newid[j], x = 1 / distance)
  if (verbose) cat('done.\n')
  similarity
}

#' Compute valid nearest neighbors of each valid point on a grid
#'
#' @param mat Matrix of gridpoints
#' @param mask_fun Function that identifies valid observations
#' For sparseMatrix where 0 denotes NA, use
#' mask_fun = function(x) as.logical(x != 0)
#'
#' @return List with IDs of valid points (old and new) and list (ngb) of their neighbors
#' @export
#'
#' @examples
#' mat = matrix(runif(6), 2, 3)
#' mat[2,2] = NA
#' getNeighbors(mat, verbose=F)
#' mat[2,2] = 0
#' getNeighbors(enforceM(mat), function(x) as.logical(x != 0), verbose=F)
#'
getNeighbors = function(mat, mask_fun = function(x) !is.na(x), verbose) {
  if (verbose) cat('Building neighborood structure for dual graph... ')
  nr = nrow(mat)
  nc = ncol(mat)
  n = nr * nc
  mask = mask_fun(mat)
  csum = cumsum(as.numeric(mask))
  newid = c(csum[1], csum[2:n] * (csum[2:n] - csum[1:(n-1)]))
  v = which(mask)
  ngb = lapply(v, FUN = function(k) {
    c = floor((k - 1) / nr) + 1
    r = (k - 1) %% nr + 1
    points = if (r == 1) c(0, 1) else if(r == nr) c(-1, 0) else c(-1, 0, 1)
    points = if (c == 1) c(points, nr) else if(c == nc) c(-nr, points) else c(-nr, points, nr)
    unlist(lapply(points + k, FUN = function(q) if (mask[q]) q else NULL))
  })
  neighbors = list(oldid = v, newid = newid, ngb = ngb)
  if (verbose) cat('done.\n')
  neighbors
}
