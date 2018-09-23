#' Compute ETOPO30 from ETOPO5. Download and extract zipped geotiff from
#' https://data.europa.eu/euodp/data/dataset/data_world-digital-elevation-model-etopo5
#'
#' @param local_tiff Path to local ETOPO5 geotiff file ('alwdgg.tiff')
#'
#' @return ETOPO30
#' @export
#'
getETOPO30 = function(local_tiff) {
  require(raster)
  a = raster(local_tiff)
  dt = t(matrix(ncol=2160, getValues(a)))
  image(Matrix(dt[seq(1,2160,by=10),seq(1,4320,by=10)]))
  etopo30=mapply(1:720, FUN=function(lo) {
    st_j = 6 * (lo - 1) + 1
    en_j = 6 * lo
    mapply(1:360, FUN = function(la) {
      st_i = 6 * (la - 1) + 1
      en_i = 6 * la
      mean(dt[st_i:en_i, st_j:en_j])
    })
  })
  etopo30
}

#' Obtain a coarser version of ETOPO30
#'
#' @param coarsening_factor Integer >= 1; resolution = 0.5 * coarsening_factor
#'
#' @return Coarser version of ETOPO30
#' @export
#'
#' @examples
#' getCoarseETOPO(4)
#'
getCoarseETOPO = function(coarsening_factor) {
  data = etopo30[seq(1, 360, by=coarsening_factor), seq(1, 720, by=coarsening_factor)]
  lat = seq(89.75, -89.75, by=-0.5 * coarsening_factor)
  lon = seq(-179.75, 179.75, by=0.5 * coarsening_factor)
  out = list(data = data, lat = lat, lon = lon)
  out
}
