#' Fit a context-dependent surface filling curve to a dataset
#'
#' @param data Dataset in matrix format
#' @param lat Vector of latitudes
#' @param lon Vector of longitudes
#' @param verbose Print status to console?
#'
#' @return Surface filling curve object (list)
#' @export
#'
#' @examples
#' sfc(data = preprocessBathymetry(etopo30, neritic = -1000),
#'     lat = seq(89.75, -89.75, by=-0.5),
#'     lon = seq(-179.75, 179.75, by=0.5),
#'     verbose = TRUE)
#'
sfc = function(data = preprocessBathymetry(etopo30, neritic = -1000),
               lat = seq(89.75, -89.75, by=-0.5),
               lon = seq(-179.75, 179.75, by=0.5),
               verbose = TRUE) {

  if (nrow(data) %% 2 == 1 || ncol(data) %% 2 == 1) {
    print("Data matrix must have an even number of rows and columns")
    return()
  }
  if (nrow(data) < 6 || ncol(data) < 6) {
    print("Data matrix is too small, must be at least 6x6")
    return()
  }
  if (nrow(data) != length(lat)) {
    print('nrow(data) should match length(lat)')
    return()
  }
  if (ncol(data) != length(lon)) {
    print('ncol(data) should match length(lon)')
    return()
  }
  dual = getDual(data = data, lat = lat, lon = lon, verbose = verbose)
  connection_matrix = getConnectionMatrix(dual = dual, verbose = verbose)
  path = getPath(connection_matrix = connection_matrix, verbose = verbose)
  group = getGroups(data = data, path = path, lat = lat, lon = lon)
  dataset = list(data = data, lat = lat, lon = lon, connection_matrix = connection_matrix,
                 path = path, group = group)
  dataset

}

#' Preprocessing function to add NAs in regions where values exceed a threshold
#'
#' @param data Raw data with no NAs
#' @param thresh Upper bound for valid values
#' @param verbose Print status to console?
#'
#' @return Data with NAs in appropriate regions
#' @export
#'
preprocess = function(data, thresh, verbose) {
  if (verbose) cat('Preprocessing data... ')

  dt = data
  nr = nrow(dt)
  nc = ncol(dt)
  categ = mapply(1:nc, FUN = function(j) {
    h = 2 * (-0.5 + (j %% 2 == 1))
    mapply(1:nr, FUN = function(i) {
      k = 2 * (-0.5 + (i %% 2 == 1))
      chunk = dt[i:(i + k), j:(j + h)]
      if (all(chunk > thresh)) {
        0 # NA region
      } else {
        1 # OK region
      }
    })
  })
  dt[categ == 0] = NA
  if (verbose) cat('done.\n')
  dt
}


#' Special preprocessing algorithm for bathymetric data
#'
#' @param data Raw data with no NAs
#' @param neritic Depth (negative value) defining end of coastal region
#'
#' @return Scaled data (0-2) with NAs in appropriate regions, values
#'         between 0 and 1 in coastal regions, and between 1 and 2 in oceanic regions
#' @export
#'
preprocessBathymetry = function(data, neritic, verbose) {
  if (verbose) cat('Preprocessing bathymetric data... ')
  dt = data
  low = min(dt)
  dt[dt > 0] = 0
  dt = dt / low # (0,1) domain
  ubound = neritic / low

  nr = nrow(dt)
  nc = ncol(dt)
  categ = mapply(1:nc, FUN = function(j) {
    h = 2 * (-0.5 + (j %% 2 == 1))
    mapply(1:nr, FUN = function(i) {
      k = 2 * (-0.5 + (i %% 2 == 1))
      chunk = dt[i:(i + k), j:(j + h)]
      if (all(chunk > ubound)) {
        3 # oceanic region
      } else if (all(chunk == 0)) {
        1 # land
      } else
        2 # neritic region
    })
  })

  correction = mapply(3:(nc-2), FUN = function(j) mapply(3:(nr-2), FUN = function(i) {
    if (categ[i,j] == 3 && (categ[i+2,j] == 1 || categ[i,j + 2] == 1 ||
        categ[i-2,j] == 1 || categ[i,j-2] == 1 ||
        categ[i+2,j+2] == 1 || categ[i+2,j - 2] == 1 ||
        categ[i-2,j+2] == 1 || categ[i-2,j - 2] == 1)) -1 else 0
  }))
  categ[3:(nr-2), 3:(nc-2)] = categ[3:(nr-2), 3:(nc-2)] + correction

  dt[categ == 3] = dt[categ == 3] + 1
  dt[categ == 1] = NA
  if (verbose) cat('done.\n')
  dt
}


