#' Great-circle distance between two points
#'
#' Haversine distance, in metres, on a sphere of radius 6,371,000 m: the radius the
#' whacky-position filters ([rb_whack_clean()] and friends) use, so a distance computed here
#' agrees with theirs. Vectorised, and the four arguments recycle to a common length, so it
#' works inside `mutate()` with `dplyr::lag()` and `dplyr::lead()`.
#'
#' @param lon1,lat1 Longitude and latitude of the first point, decimal degrees.
#' @param lon2,lat2 Longitude and latitude of the second point, decimal degrees.
#'
#' @return Numeric vector of distances in metres; `NA` where any coordinate is `NA`.
#'
#' @examples
#' rb_distance(0, 64, 0, 65)                 # one degree of latitude, 111,195 m
#' rb_distance(-21.93, 64.15, -23.93, 65.65) # Reykjavik to Talknafjordur, about 190 km
#' \dontrun{
#' pings |>
#'   group_by(vid) |>
#'   mutate(m = rb_distance(dplyr::lag(lon), dplyr::lag(lat), lon, lat)) |>
#'   ungroup()
#' }
#'
#' @export
rb_distance <- function(lon1, lat1, lon2, lat2) {
  k <- pi / 180
  a <- sin((lat2 - lat1) * k / 2)^2 + cos(lat1 * k) * cos(lat2 * k) * sin((lon2 - lon1) * k / 2)^2
  2 * 6371000 * asin(sqrt(pmin(1, a)))
}
