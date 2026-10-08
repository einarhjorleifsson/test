test_that("rb_whack_forward matches the labelled fixture whacks1", {
  d <- whacks1[order(whacks1$vid, whacks1$time), ]
  out <- rb_whack_forward(d, kn_max = 25)
  expect_equal(out$whack2, d$whacks)
})

test_that("rb_whack_forward's compiled core matches the original R loop", {
  # Pre-2026-09-23 reference implementation of .fwd() (time-step floor min_dt since
  # 2026-10-08, as in the port), kept here only to
  # verify the Rcpp port (src/whack_forward.cpp) is byte-identical, including
  # its NA-propagation quirk (see R/rb_whack_clean.R and the port's header comment).
  ref_fwd <- function(lon, lat, time, kn_max = 25, max_gap_h = 4, min_dt = 10) {
    ms_max      <- kn_max * 0.514444
    max_gap_sec <- max_gap_h * 3600
    r           <- 6371000
    n    <- length(lon)
    flag <- rep(FALSE, n)
    if (n < 2L) return(flag)
    prev <- 1L
    for (i in 2:n) {
      if (is.na(lon[i]) || is.na(lat[i]) || is.na(time[i])) {
        prev <- i
        next
      }
      dt <- as.numeric(time[i] - time[prev], units = "secs")
      if (is.na(dt) || dt > max_gap_sec) {
        prev <- i
        next
      }
      phi1 <- lat[prev] * pi / 180;  phi2 <- lat[i] * pi / 180
      dphi <- (lat[i] - lat[prev]) * pi / 180
      dlam <- (lon[i] - lon[prev]) * pi / 180
      d    <- 2 * r * asin(pmin(1, sqrt(
        sin(dphi/2)^2 + cos(phi1) * cos(phi2) * sin(dlam/2)^2)))
      if (is.na(d) || d / max(dt, min_dt) > ms_max) {
        flag[i] <- TRUE
      } else {
        prev <- i
      }
    }
    flag
  }

  set.seed(42)
  n <- 5000
  lon  <- -25 + cumsum(rnorm(n, 0, 0.0002))
  lat  <-  62 + cumsum(rnorm(n, 0, 0.0002))
  time <- as.POSIXct("2020-01-01", tz = "UTC") + sort(sample(1:200000, n))
  lon[100]         <- lon[100] + 5           # isolated spike
  lat[500:510]     <- lat[500:510] + 3       # an 11-point run
  time[2000:n]     <- time[2000:n] + 5 * 3600 # a >4h gap (tests max_gap_h reset)

  expect_identical(
    ref_fwd(lon, lat, time),
    test:::whack_forward_cpp(lon, lat, as.numeric(time), 25 * 0.514444, 4 * 3600, 10)
  )

  # A single NA lon (valid lat/time) — exercises the NA-propagation cascade
  # the port intentionally reproduces rather than "fixing" here.
  lon_na <- lon
  lon_na[50] <- NA_real_
  expect_identical(
    ref_fwd(lon_na, lat, time),
    test:::whack_forward_cpp(lon_na, lat, as.numeric(time), 25 * 0.514444, 4 * 3600, 10)
  )
})

test_that("rb_whack_forward handles degenerate inputs", {
  empty <- whacks1[0, ]
  expect_equal(nrow(rb_whack_forward(empty)), 0)

  one <- whacks1[1, ]
  out <- rb_whack_forward(one)
  expect_equal(out$whack2, FALSE)
})

test_that("a second report in the same second is judged by its distance", {
  t0 <- as.POSIXct("2022-10-15 00:00:00", tz = "UTC")
  d <- data.frame(vid = 1, time = t0 + 60 * c(0, 1, 1, 2, 3),
                  lon = -23 + 0.0003 * c(0, 1, 1, 2, 3), lat = 64 + 0.0002 * c(0, 1, 1, 2, 3))
  d$lon[3] <- d$lon[3] + 0.00001                      # duplicate, about 0.5 m away
  expect_false(any(rb_whack_forward(d)$whack2))
  d$lon[3] <- d$lon[3] + 0.1                          # same second, 5 km away
  expect_equal(which(rb_whack_forward(d)$whack2), 3L)
})
