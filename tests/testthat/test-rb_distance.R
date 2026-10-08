test_that("rb_distance() gives haversine metres on the filters' sphere", {
  expect_equal(rb_distance(0, 64, 0, 65), 2 * pi * 6371000 / 360)
  expect_equal(rb_distance(-21.9, 64.1, -21.9, 64.1), 0)
  expect_equal(rb_distance(10, 0, 11, 0), 2 * pi * 6371000 / 360)       # one degree of longitude at the equator
  expect_equal(rb_distance(-22, 64, -24, 65), rb_distance(-24, 65, -22, 64))
})

test_that("rb_distance() is vectorised, recycles and passes NA through", {
  d <- rb_distance(c(0, 0, NA), 64, 0, c(65, 64, 64))
  expect_length(d, 3)
  expect_equal(d[1:2], c(2 * pi * 6371000 / 360, 0))
  expect_true(is.na(d[3]))
  expect_equal(rb_distance(0, 0, 180, 0), pi * 6371000)                  # antipodes: asin() argument capped at 1
})
