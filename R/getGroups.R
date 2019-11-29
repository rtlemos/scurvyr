#' Break path into groups, according to mean and variance
#'
#' @param data Dataset
#' @param path Path
#' @param lat Path latitudes
#' @param lon Path longitudes
#'
#' @return data.frame with 2 columns: group and mean value of group
#' @export
#'
getGroups <- function(data, path, lat, lon) {
  nr <- length(lat)
  path_col <- floor((path - 1) / nr) + 1
  path_row <- (path - 1) %% nr + 1
  values <- mapply(path_row, path_col, FUN=function(rr, cc) data[rr, cc])
  mycpt <- cpt.meanvar(data = values, method = "PELT")
  n <- length(mycpt@cpts)
  if (n > 1) {
    st <- c(1, mycpt@cpts[-n])
    cnt <- c(mycpt@cpts[1], mycpt@cpts[2:n] - mycpt@cpts[1:(n - 1)])
  } else {
    st <- 1
    cnt <- mycpt@cpts[1]
  }
  group_id <- unlist(mapply(cnt, 1:n, SIMPLIFY = FALSE, FUN = function(cc, ii) rep(ii, cc)))
  group_value <- unlist(mapply(cnt, mycpt@param.est$mean, SIMPLIFY = FALSE, FUN = function(cc, m) rep(m, cc)))
  group <- data.frame(id = group_id, value = group_value)
  return(group)
}
