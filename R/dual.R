#' Compute the dual graph object
#'
#' @param dataset
#'
#' @return List describing the dual graph object
#' @export
#'
getDual = function(data, lat, lon, verbose) {
  coarse = getCoarseData(data = data, lat = lat, lon = lon, verbose = verbose)
  neighbors = getNeighbors(mat = coarse$data, verbose = verbose)
  similarity_matrix = getSimilarityMatrix(data = data, neighbors = neighbors, verbose = verbose)
  mstree = minimumSpanningTree(similarity_matrix, verbose = verbose)
  dual = list(data = coarse$data, lat = coarse$lat, lon = coarse$lon,
              neighbors = neighbors, similarity_matrix = similarity_matrix, mstree = mstree)
  return(dual)
}

#' Build a coarser dataset
#'
#' @param data Original dataset
#' @param lat Vector of latitudes
#' @param lon Vector of longitudes
#' @param verbose Print statuts to console?
#'
#' @return List with coarser data, lat and lon
#' @export
#'
getCoarseData = function(data, lat, lon, verbose) {
  if (verbose) cat('Building coarser gridded dataset... ')
  nr = length(lat)
  nc = length(lon)
  res = abs(lat[2] - lat[1])
  clat = seq(max(lat) - res / 2, min(lat) + res / 2, by = -2 * res)
  clon = seq(min(lon) + res / 2, max(lon) + res / 2, by = 2 * res)
  vals = mapply(1:floor(nc/2), FUN = function(j){
    stj = (j - 1) * 2 + 1
    enj = j * 2
    mapply(1:floor(nr/2), FUN = function(i) {
      sti = (i - 1) * 2 + 1
      eni = i * 2
      dt = data[sti:eni,stj:enj]
      if (all(is.na(dt))) NA else mean(dt, na.rm = TRUE)
    })
  })
  coarse = list(data = vals, lat = clat, lon = clon)
  if (verbose) cat('done.\n')
  coarse
}

#' Build a connection matrix based on a minimum spanning tree
#'
#' @param mstree Minimum Spanning Tree
#' @param dual Dual object
#' @param Verbose Print status to console
#'
#' @return (4 x nr x nc) x (4 x nr x nc) sparse matrix of connections in dual tree
#' @export
#'
getConnectionMatrix = function(dual, verbose) {
  if (verbose) cat('Building connection matrix... ')

  neighbors = dual$neighbors
  mstree = dual$mstree
  nr = length(dual$lat)
  nc = length(dual$lon)
  nr2 = 2 * nr
  nc2 = 2 * nc

  boundary = t(matrix(nrow = 2, c(
    unlist(lapply(1:(nr * nc), FUN = function(i) {
      ri = (i - 1) %% nr + 1
      ci = floor((i - 1) / nr) + 1
      rri = (ri - 1) * 2 + 1
      cci = (ci - 1) * 2 + 1
      id = (cci - 1) * nr2 + rri
      lst = NULL
      if (neighbors$newid[i] != 0) {
        if (ri == 1 || neighbors$newid[i - 1] == 0) { #no valid neighbor above
          lst = cbind(lst, c(id, id + nr2))
        }
        if (ri == nr || neighbors$newid[i + 1] == 0) { #no valid neighbor below
          lst = cbind(lst, c(id + 1, id + nr2 + 1))
        }
        if (ci == 1 || neighbors$newid[i - nr] == 0) { #no valid neighbor to the left
          lst = cbind(lst, c(id, id + 1))
        }
        if (ci == nc || neighbors$newid[i + nr] == 0) { #no valid neighbor to the right
          lst = cbind(lst, c(id + nr2, id + nr2 + 1))
        }
      }
      lst
    }))
  )))
  boundary = rbind(boundary, cbind(boundary[,2], boundary[, 1]))
  ij = t(matrix(nrow = 2, unlist(mapply(neighbors$ngb, neighbors$oldid, FUN = function(ngb, i) {
    ii = neighbors$newid[i]
    ci = floor((i - 1) / nr) + 1
    ri = (i - 1) %% nr + 1
    unlist(lapply(ngb, FUN = function(j) {
      if (j <= i) return(NULL)
      jj = neighbors$newid[j]
      cj = floor((j - 1) / nr) + 1
      rj = (j - 1) %% nr + 1

      if (ri == rj) { #same row, different columns
        if (abs(ci - cj) != 1) stop(paste('ci and cj should differ by only 1 unit', ci, cj))
        r = (ri - 1) * 2 + 1
        c = min(ci, cj) * 2
        if (mstree[ii, jj] == 1) {
          matrix(nrow = 2, c((c - 1) * nr2 + r, c * nr2 + r, (c - 1) * nr2 + r + 1, c * nr2 + r + 1))
        } else {
          matrix(nrow = 2, c((c - 1) * nr2 + r, (c - 1) * nr2 + r + 1, c * nr2 + r, c * nr2 + r + 1))
        }
      } else { # same column, different rows
        if (abs(ri - rj) != 1) stop(paste('ci and cj should differ by only 1 unit', ri, rj))
        r = min(ri, rj) * 2
        c = (ci - 1) * 2 + 1
        if (mstree[ii, jj] == 1) {
          matrix(nrow = 2, c((c - 1) * nr2 + r, (c - 1) * nr2 + r + 1, c * nr2 + r, c * nr2 + r + 1))
        } else {
          matrix(nrow = 2, c((c - 1) * nr2 + r, c * nr2 + r, (c - 1) * nr2 + r + 1, c * nr2 + r + 1))
        }
      }
    }))
  }))))
  ij = rbind(ij, cbind(ij[,2], ij[,1]))

  links = rbind(ij, boundary)
  connection_matrix = sparseMatrix(i = links[,1], j = links[,2], x = rep(1, nrow(links)), dims=rep(nr2 * nc2, 2))
  if (verbose) cat('done.\n')
  connection_matrix
}

#' Compute a path that connects all the points in a matrix
#'
#' @param connection_matrix Connection matrix for dual graph
#'
#' @return Path that connects all elements in connection matrix
#' @export
#'
getPath = function(connection_matrix, verbose) {
  if (verbose) cat('Building path for dual tree... ')
  nr = nrow(connection_matrix)
  i = connection_matrix@i + 1
  j = unlist(lapply(1:nr, FUN = function(jj) {
    st = connection_matrix@p[jj]
    en = connection_matrix@p[jj + 1]
    if (en > st) rep(jj, en - st) else NULL
  }))

  recursivePath = trampoline(function(mypath, free) {
    neighb = j[i == mypath[1]]
    crit = free[neighb]
    if (sum(crit) == 0) {
      return(mypath)
    } else {
      point = neighb[crit][1]
      free[point] = FALSE
      mypath = c(point, mypath)
      recur(mypath, free)
    }
  })

  cmax = getColMax(connection_matrix, 1:nr)
  st = which(cmax[,1] > 0)[1]
  mypath = recursivePath(st, rep(TRUE, nr))
  if (verbose) cat('done.\n')
  mypath
}
