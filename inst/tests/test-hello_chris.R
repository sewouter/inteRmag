library(testthat)

test_that("hello_chris outputs correctly", {
  expect_that(hello_chris("a"), prints_text("Are you a Chris ?"))
  expect_that(hello_chris("l"), prints_text("Hello Chris !"))
  expect_that(hello_chris("z"), prints_text("Hello Christian !"))
  expect_that(hello_chris("r"), prints_text("Hi Mr Rolf, nice to meet you !"))
})
