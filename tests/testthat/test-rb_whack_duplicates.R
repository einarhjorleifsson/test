t0 <- as.POSIXct("2020-06-01 12:00:00", tz = "UTC")
dlon <- function(m) m / (111194.9 * cos(64 * pi / 180))   # degrees of longitude for m metres at 64 N

test_that("the lower-ranked copy of a report is flagged, the higher-ranked one kept", {
  d <- rb_whack_duplicates(c(1, 1), t0 + c(0, 1), c(-22, -22 + dlon(0.5)), c(64, 64), c("stk", "astd"))
  expect_equal(d, c(FALSE, TRUE))
  d <- rb_whack_duplicates(c(1, 1), t0 + c(0, 1), c(-22, -22 + dlon(0.5)), c(64, 64), c("astd", "stk"))
  expect_equal(d, c(TRUE, FALSE))
})

test_that("a pair outside the time or distance window, or of two vessels, is not a copy", {
  expect_equal(rb_whack_duplicates(c(1, 1), t0 + c(0, 11), c(-22, -22), c(64, 64), c("stk", "astd")), c(FALSE, FALSE))
  expect_equal(rb_whack_duplicates(c(1, 1), t0 + c(0, 5), c(-22, -22 + dlon(200)), c(64, 64), c("stk", "astd")), c(FALSE, FALSE))
  expect_equal(rb_whack_duplicates(c(1, 2), t0 + c(0, 1), c(-22, -22), c(64, 64), c("stk", "astd")), c(FALSE, FALSE))
  # the window is kn_max over max_dt_s: 128.6 m at the defaults
  expect_equal(rb_whack_duplicates(c(1, 1), t0 + c(0, 1), c(-22, -22 + dlon(120)), c(64, 64), c("stk", "astd")), c(FALSE, TRUE))
  expect_equal(rb_whack_duplicates(c(1, 1), t0 + c(0, 1), c(-22, -22 + dlon(140)), c(64, 64), c("stk", "astd")), c(FALSE, FALSE))
})

test_that("same-feed and unranked pairs are left alone; an unranked feed loses to a ranked one", {
  expect_equal(rb_whack_duplicates(c(1, 1), t0 + c(0, 0), c(-22, -22), c(64, 64), c("astd", "astd")), c(FALSE, FALSE))
  expect_equal(rb_whack_duplicates(c(1, 1), t0 + c(0, 1), c(-22, -22), c(64, 64), c("x", "y")), c(FALSE, FALSE))
  expect_equal(rb_whack_duplicates(c(1, 1), t0 + c(0, 1), c(-22, -22), c(64, 64), c("x", "emodnet")), c(TRUE, FALSE))
})

test_that("three feeds on one fix keep only the top one, whatever the input order", {
  vid  <- c(1, 1, 1, 1)
  time <- t0 + c(0, 2, 3, 300)
  lon  <- c(-22, -22 + dlon(1), -22 + dlon(2), -22.1)
  lat  <- rep(64, 4)
  feed <- c("stk", "astd", "astdB", "astd")
  expect_equal(rb_whack_duplicates(vid, time, lon, lat, feed), c(FALSE, TRUE, TRUE, FALSE))
  o <- c(4, 2, 1, 3)
  expect_equal(rb_whack_duplicates(vid[o], time[o], lon[o], lat[o], feed[o]), c(FALSE, TRUE, TRUE, FALSE)[o])
})

test_that("short and NA input do not fail", {
  expect_equal(rb_whack_duplicates(1, t0, -22, 64, "stk"), FALSE)
  expect_equal(rb_whack_duplicates(numeric(0), t0[0], numeric(0), numeric(0), character(0)), logical(0))
  expect_equal(rb_whack_duplicates(c(1, 1), t0 + c(0, 1), c(-22, NA), c(64, 64), c("stk", "astd")), c(FALSE, FALSE))
})
