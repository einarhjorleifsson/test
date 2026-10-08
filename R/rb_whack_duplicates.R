#' Near-duplicate pings across AIS feeds
#'
#' Flags the second copy of a position report that two feeds both received: the
#' same vessel, a few seconds apart, metres apart. When feeds are merged (for
#' example a national receiver network and a commercial terrestrial or satellite
#' feed), a vessel heard by both is reported twice. The two copies do not share a
#' timestamp or a position exactly, so they survive a match on `(vid, time)` or on
#' whole rows; they are near-duplicates. In the Icelandic record 2013-2026 they
#' were 1.2 % of all pings, with a median gap of 1-4 s and a median distance of
#' 0.1-0.4 m by feed pair.
#'
#' Run it before [rb_whack_clean()]. With a short time-step floor a near-duplicate
#' pair that differs by a few metres looks like an impossible speed and gets
#' flagged as whacky; with `min_dt_s = 10` it does not, so it has to be removed
#' as what it is.
#'
#' @section Rule:
#' A ping is a near-duplicate when a ping of the same `vid` from a feed ranked
#' higher in `priority` lies within `max_dt_s` seconds and within the distance
#' `kn_max` allows in `max_dt_s` (128.6 m at the defaults, so the pair is one the
#' speed filter would itself call consistent). Only the `n` nearest pings on each
#' side, in time order, are compared. The top-ranked feed is never flagged. Feeds
#' not named in `priority` rank last and equal, so two of them never flag each
#' other; nor do two pings of the same feed: this is not a filter for exact or
#' same-feed duplicates.
#'
#' @param vid Vessel id.
#' @param time Time of the ping (`POSIXct`, or seconds).
#' @param lon,lat Position, decimal degrees.
#' @param feed The feed (provider) each ping came from.
#' @param priority Feeds in order of preference, the one to keep first.
#' @param max_dt_s Largest time gap, in seconds, for two pings to be copies (default 10,
#'   the time-step floor of [rb_whack_clean()]).
#' @param kn_max Speed, in knots, that sets the largest distance: `kn_max` over
#'   `max_dt_s` (default 25, as in [rb_whack_clean()]).
#' @param n Neighbours compared on each side, in time order (default 3).
#'
#' @return A logical vector, `TRUE` for the copy to drop, in the order of the input.
#'
#' @examples
#' t0 <- as.POSIXct("2020-06-01 12:00:00", tz = "UTC")
#' rb_whack_duplicates(vid  = c(1, 1, 1),
#'                     time = t0 + c(0, 2, 60),
#'                     lon  = c(-22, -22.000004, -22.01),
#'                     lat  = c(64, 64, 64),
#'                     feed = c("stk", "astd", "astd"))
#' # FALSE TRUE FALSE: the astd ping 2 s and 0.2 m from the stk one is its copy
#' \dontrun{
#' pings |>
#'   mutate(dup = rb_whack_duplicates(vid, time, lon, lat, provider)) |>
#'   filter(!dup) |>
#'   rb_whack_clean()
#' }
#'
#' @export
rb_whack_duplicates <- function(vid, time, lon, lat, feed,
                                priority = c("stk", "astd", "astdB", "emodnet"),
                                max_dt_s = 10, kn_max = 25, n = 3) {
  max_m <- kn_max * 1852 / 3600 * max_dt_s
  pr <- match(feed, priority); pr[is.na(pr)] <- length(priority) + 1L
  t  <- as.numeric(time)
  o  <- order(vid, t, pr, lon, lat)
  v <- vid[o]; t <- t[o]; x <- lon[o]; y <- lat[o]; p <- pr[o]
  m <- length(o); dup <- logical(m)
  if (m < 2L) return(dup)
  for (k in seq_len(min(n, m - 1L))) {
    i <- seq_len(m - k); j <- i + k
    near <- v[i] == v[j] & abs(t[j] - t[i]) <= max_dt_s & rb_distance(x[i], y[i], x[j], y[j]) <= max_m
    near[is.na(near)] <- FALSE
    dup[j[near & p[i] < p[j]]] <- TRUE
    dup[i[near & p[j] < p[i]]] <- TRUE
  }
  dup[order(o)]
}
