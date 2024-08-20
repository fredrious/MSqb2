library(testthat)

test_that("expandList imports objects from a named list to the parent environment", {
  # Define variables and create a named list
  a <- 1
  b <- 2
  my_list <- list(a = a, b = b)

  # Create a new environment to test the function
  env <- new.env()

  # Use the new environment as the parent environment
  with(env, {
    c <- 3
    # Expand the list into the parent environment
    expandList(my_list)
  })

  # Check if variables `a` and `b` are now in the parent environment of `env`
  expect_true(exists("a", envir = parent.env(env)))
  expect_true(exists("b", envir = parent.env(env)))

  # Check the values of the imported variables
  expect_equal(parent.env(env)$a, 1)
  expect_equal(parent.env(env)$b, 2)
})
