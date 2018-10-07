#' ETOPO30 elevation data
#'
#' A dataset containing global elevation on a 0.5 degree grid, based on ETOPO5 (see getETOPO30).
#'
#' @format A data frame with 360 rows and 720 columns; etopo[1,1] corresponds to 89.75N, -179.75E,
#'         etopo30[1,2] to 89.75N, -179.25E, and so on.
#' @source \url{https://data.europa.eu/euodp/data/dataset/data_world-digital-elevation-model-etopo5}
"etopo30"

#' ETOPO30 bathymetry path
#'
#' An array connecting oceanic points, using the Surface Filling Curve algorithm in this package
#'
#' @format A data frame with 175469 rows and 3 columns: latitude, longitude, and path ID number
#' @source The path and the lat/lon can be created afresh by typing
#'         > data = preprocessBathymetry(etopo30, neritic = -2000, verbose = TRUE)
#'         > lat = seq(89.75, -89.75, by=-0.5)
#'         > lon = seq(-179.75, 179.75, by=0.5)
#'         > sEtopo = sfc(data = data, lat = lat, lon = lon, verbose = TRUE)
#'         > etopoPath = data.frame(lat = lat[(sEtopo$path - 1) %% 360 + 1],
#'         >                        lon = lon[floor((sEtopo$path - 1) / 360) + 1],
#'         >                        path = sEtopo$path)
"etopoPath"