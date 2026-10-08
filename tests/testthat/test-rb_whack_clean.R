mk <- function(n, start, lon0, vid = 1) {
  data.frame(vid = vid, time = start + 60 * (0:(n - 1)),
             lon = lon0 + 0.0003 * (0:(n - 1)), lat = 64 + 0.0002 * (0:(n - 1)))
}
t0 <- as.POSIXct("2022-10-15 00:00:00", tz = "UTC")

test_that("rb_whack_clean labels, never drops, and sorts by vid, time", {
  d <- whacks1[sample(nrow(whacks1)), ]
  out <- rb_whack_clean(d)
  expect_equal(nrow(out), nrow(d))
  expect_false(is.unsorted(order(out$vid, out$time)))
  expect_true(all(is.na(out$whack_stage) == !out$whack))
})

test_that("rb_whack_clean flags every labelled whack in the fixture", {
  d <- whacks1[order(whacks1$vid, whacks1$time), ]
  out <- rb_whack_clean(d)
  expect_true(all(out$whack[d$whacks]))
})

test_that("an isolated spike is flagged by the first stage (vmask or spike step)", {
  d <- mk(40, t0, -23); d$lon[20] <- d$lon[20] + 2
  out <- rb_whack_clean(d)
  expect_equal(which(out$whack), 20L)
  expect_true(out$whack_stage[20] %in% c("vmask", "spike"))
})

test_that("a bad first ping after a long gap is flagged, not adopted as the anchor", {
  d <- rbind(mk(15, t0, -23), mk(15, t0 + 4.5 * 3600, -22.99))
  d$lon[16] <- d$lon[16] + 6
  out <- rb_whack_clean(d)
  expect_equal(which(out$whack), 16L)
  # the forward scan alone gets this wrong: it keeps the bad ping and flags the good ones
  expect_equal(sum(rb_whack_forward(d)$whack2), 14L)
})

test_that("vessels are independent", {
  a <- mk(40, t0, -23, vid = 1); a$lon[20] <- a$lon[20] + 2
  b <- mk(40, t0, -20, vid = 2)
  both <- rb_whack_clean(rbind(a, b))
  expect_equal(rb_whack_clean(a)$whack, both$whack[both$vid == 1])
  expect_equal(rb_whack_clean(b)$whack, both$whack[both$vid == 2])
})

# argosfilter adds 1 s to every time step; rb_whack_sda floors the step at min_dt_s instead.
# With steps of 1-15 min the two speeds differ by under 2 %, and the flags should not.
test_that("rb_whack_sda flags as argosfilter::sdafilter at its own settings", {
  skip_if_not_installed("argosfilter")
  set.seed(1)
  n <- 400
  d <- data.frame(vid = 1, time = t0 + cumsum(sample(60:900, n, TRUE)),
                  lon = -23 + cumsum(rnorm(n, 0, 0.01)), lat = 64 + cumsum(rnorm(n, 0, 0.01)))
  k <- sample(10:390, 8); d$lon[k] <- d$lon[k] + rnorm(8, 0, 0.5)
  r <- argosfilter::sdafilter(d$lat, d$lon, d$time, rep("1", n), vmax = 12.86)
  ours <- rb_whack_sda(d, kn_max = 12.86 / 0.514444, ang = c(15, 25), distlim = c(2500, 5000),
                    speedlim_kn = c(0, 0), vmask_min_dist = 5000)
  ref <- r == "removed"
  # the R function labels the first/last two pings "end_location" and never removes them
  expect_equal(!is.na(ours$whack_sda)[-c(1:2, n - 1, n)], ref[-c(1:2, n - 1, n)])
})

test_that("rb_whack_clean on a DuckDB table equals the data-frame result", {
  skip_if_not_installed("duckdb"); skip_if_not_installed("dbplyr")
  a <- mk(60, t0, -23, vid = 1); a$lon[30] <- a$lon[30] + 2
  b <- rbind(mk(15, t0, -20, vid = 2), mk(15, t0 + 4.5 * 3600, -19.99, vid = 2)); b$lon[16] <- b$lon[16] + 6
  d <- rbind(a, b)
  con <- DBI::dbConnect(duckdb::duckdb()); on.exit(DBI::dbDisconnect(con, shutdown = TRUE))
  lz <- dplyr::copy_to(con, d, "pings")
  expect_s3_class(rb_whack_clean(lz, batch_pings = 70), "tbl_lazy")   # two batches
  got <- dplyr::collect(rb_whack_clean(lz, batch_pings = 70)) |> dplyr::arrange(vid, time)
  ref <- rb_whack_clean(d)
  expect_equal(got$whack, ref$whack)
  expect_equal(got$whack_stage, ref$whack_stage)
  expect_equal(nrow(got), nrow(d))
})

test_that("the result does not depend on the row order of the input (ties in time)", {
  d <- mk(40, t0, -23); d <- rbind(d, transform(d[20, ], lon = lon + 1, lat = lat + 1))   # two pings, same time
  set.seed(2)
  a <- rb_whack_clean(d); b <- rb_whack_clean(d[sample(nrow(d)), ])
  expect_equal(a$whack, b$whack)
  expect_equal(a$lon, b$lon)
})

test_that("the lazy result can be counted and sorted in DuckDB", {
  skip_if_not_installed("duckdb"); skip_if_not_installed("dbplyr")
  a <- mk(60, t0, -23, vid = 1); a$lon[30] <- a$lon[30] + 2
  con <- DBI::dbConnect(duckdb::duckdb()); on.exit(DBI::dbDisconnect(con, shutdown = TRUE))
  lz <- dplyr::copy_to(con, a, "pings")
  got <- rb_whack_clean(lz) |> dplyr::count(whack_stage) |> dplyr::arrange(whack_stage) |> dplyr::collect()
  expect_equal(sum(got$n), 60)
  expect_equal(got$n[!is.na(got$whack_stage)], 1)
})

test_that("a step shorter than min_dt_s counts as min_dt_s, in all three filters", {
  d <- data.frame(vid = 1, time = t0 + 5 * (0:29), lon = -23, lat = 64 + 0.000045 * (0:29))
  d$lon[15] <- d$lon[15] + 0.00205        # 100 m off the track, 5 s from either neighbour
  for (m in c(1, 10)) {
    hit <- m == 1                          # 100 m over 5 s is 20 m/s; over 10 s, 10 m/s (< 12.9)
    expect_equal(any(rb_whack_clean(d, min_dt_s = m)$whack), hit)
    expect_equal(any(!is.na(rb_whack_sda(d, min_dt_s = m)$whack_sda)), hit)
    expect_equal(any(rb_whack_forward(d, min_dt_s = m)$whack2), hit)
  }
  skip_if_not_installed("duckdb"); skip_if_not_installed("dbplyr")
  con <- DBI::dbConnect(duckdb::duckdb()); on.exit(DBI::dbDisconnect(con, shutdown = TRUE))
  lz <- dplyr::copy_to(con, d, "pings")
  for (m in c(1, 10)) {
    expect_equal(any(dplyr::pull(rb_whack_clean(lz, min_dt_s = m), whack)), m == 1)
  }
})
