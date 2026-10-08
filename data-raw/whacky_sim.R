# whacky_sim.R — seven simulated vessel tracks, one scenario each, for the examples and the README.
#
# Output: data/whacky_sim.rda
#
# Tracks are drawn in kilometres from a start point and turned into lon/lat at 64 N, so every
# distance below is what the filter sees. Pings every 5 minutes, GPS noise of 10 m. `truth` is TRUE
# for a ping that was moved off its track on purpose: the error the filter should find.

set.seed(2026)
t0 <- as.POSIXct("2000-01-01 06:00:00", tz = "UTC")

# n pings every 5 min at kn knots; the heading turns by `turn` radians per ping
leg <- function(n, kn, start = t0, x0 = 0, y0 = 0, hdg0 = 0.4, turn = 0.04) {
  step <- kn * 1.852 / 12                       # km per 5 minutes
  h <- hdg0 + turn * (seq_len(n) - 1)
  data.frame(time = start + 300 * (seq_len(n) - 1),
             x = x0 + cumsum(c(0, step * sin(h[-1]))) + rnorm(n, 0, 0.01),
             y = y0 + cumsum(c(0, step * cos(h[-1]))) + rnorm(n, 0, 0.01))
}
# the next leg, gap_h hours after the last ping of `prev`, moved by (dx, dy) km
then <- function(prev, gap_h, n, kn, dx = 0, dy = 0, ...) {
  last <- prev[nrow(prev), ]
  leg(n, kn, start = last$time + gap_h * 3600, x0 = last$x + dx, y0 = last$y + dy, ...)
}
# move pings i by (dx, dy) km and mark them as the truth
push <- function(d, i, dx, dy) {
  d$x[i] <- d$x[i] + dx; d$y[i] <- d$y[i] + dy; d$truth[i] <- TRUE
  d
}
track <- function(...) { d <- rbind(...); d$truth <- FALSE; d }

a <- leg(30, 10)
b <- leg(20, 4, turn = 0.3)
f <- then(b, 5 / 60, 12, 30, hdg0 = 1, turn = 0)
sims <- list(
  "1 spike"                 = track(leg(60, 10)) |> push(30, 12, -6),          # one ping 13 km off
  "2 burst aside"           = track(leg(60, 10)) |> push(25:28, 0, 25),        # four pings 25 km to the side
  "3 burst ahead"           = track(leg(60, 10)) |> push(25:28, 25, 0),        # four pings 25 km along the heading
  "4 burst near"            = track(leg(60, 10)) |> push(25:28, 0, 9),         # four pings 9 km to the side
  "5 bad ping after a gap"  = track(a, then(a, 5, 30, 10, hdg0 = 1.2)) |> push(31, 15, 10),  # first ping after 5 h
  "6 real move after a gap" = track(a, then(a, 6, 30, 10, dx = 40, dy = 25, hdg0 = 2)),     # 47 km in 6 h
  "7 fast vessel"           = track(b, f, then(f, 5 / 60, 20, 4, turn = 0.3)))   # fishing 4 kn, steaming 30 kn

whacky_sim <- do.call(rbind, lapply(names(sims), function(s) {
  d <- sims[[s]]
  data.frame(vid = as.integer(substr(s, 1, 1)), scenario = s, time = d$time,
             lon = -22 + d$x / (111.195 * cos(64 * pi / 180)), lat = 64 + d$y / 111.195, truth = d$truth)
}))
save(whacky_sim, file = "data/whacky_sim.rda", compress = "xz")
