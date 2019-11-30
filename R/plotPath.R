#' Plot a path
#'
#' @param dataset List with data, path, lat, lon, etc.
#' @param plot_data Plot the data?
#' @param lat_bounds Latitude bounds (for zooming in)
#' @param lon_bounds Longitude bounds (for zooming in)
#' @param colored_line Color the path: NA=>no, "group"=> by group ID, "value"=> by group value
#' @param fixed_aspect_ratio Use a 1:1 aspect ratio in the plot?
#' @param background Background color
#' @param show_line_scale Display colorbar for line?
#' @param show_fill_scale Display colorbar for fill?
#'
#' @return ggplot of path
#' @export
#'
plotPath = function(dataset,
                    plot_data = FALSE,
                    lat_bounds = NULL,
                    lon_bounds = NULL,
                    colored_line = 'value', 
                    fixed_aspect_ratio = TRUE,
                    background = NA,
                    show_line_scale = FALSE,
                    show_fill_scale = FALSE) {
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
  gi = dataset$group$id[1:(n - 1)]
  gj = dataset$group$id[2:n]
  vi = dataset$group$value[1:(n - 1)]
  vj = dataset$group$value[2:n]
  
  if (plot_data) {
    rst = cbind(expand.grid(lat=lat,lon=lon), z=as.numeric(data))
    if (restricted_bounds) {
      crit = rst$lat >= lat_bounds[1] & rst$lat <= lat_bounds[2] &
        rst$lon >= lon_bounds[1] & rst$lon <= lon_bounds[2]
      rst = rst[crit, ]
    }
    p = ggplot() + geom_raster(rst, mapping = aes(x=lon, y=lat, fill=z)) +
      scale_fill_gradient(low=gray(0.01), high=gray(0.99))
  } else {
    p = ggplot()
  }
  
  df = data.frame(lat = lat[c(rbind(ri, rj))], 
                  lon = lon[c(rbind(ci, cj))],
                  id = c(rbind(gi, gj)),
                  val = c(rbind(vi, vj)),
                  idx = c(rbind(1:(n - 1), 1:(n - 1))))
  if (restricted_bounds) {
    crit = df$lat >= lat_bounds[1] & df$lat <= lat_bounds[2] &
      df$lon >= lon_bounds[1] & df$lon <= lon_bounds[2]
    df = df[crit,]
  }
  if (fixed_aspect_ratio) {
    p = p + coord_equal()
  }
  if (colored_line %in% c("group", "group_id", "id")) {
    p = p + geom_line(data = df, mapping=aes(x=lon, y=lat, group=idx, color = id)) +
      scale_color_gradientn(colours=rainbow(100))
  } else if (colored_line %in% c("val", "value", "v")) {
    p = p + geom_line(data = df, mapping=aes(x=lon, y=lat, group=idx, color = val)) +
      scale_color_viridis_c(option='D')
  } else {
    p = p + geom_line(data=df, mapping=aes(x=lon, y=lat, group=idx))
  }
  if (restricted_bounds) {
    p = p + scale_x_continuous(limits = lon_bounds, expand = c(0, 0)) +
      scale_y_continuous(limits = lat_bounds, expand = c(0, 0))
  }
  clr <- if (show_line_scale) 'colorbar' else 'none'
  fll <- if (show_fill_scale) 'colorbar' else 'none'
  p = p + guides(color=clr, fill=fll) +
    theme(panel.background = element_rect(fill = background),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  p
  
}