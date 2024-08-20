library(testthat)

test_that("expTopEnv exports objects to parent environment", {
  # Define variables in the current environment
  a <- 1
  b <- 2

  # Create a new environment to test the function
  env <- new.env()

  # Use the new environment as the parent environment
  with(env, {
    c <- 3
    # Export objects from the new environment to its parent
    expTopEnv()
  })

  # Check if variables `a` and `b` are now in the parent environment of `env`
  expect_true(exists("a", envir = parent.env(env)))
  expect_true(exists("b", envir = parent.env(env)))

  # Check the values of the exported variables
  expect_equal(parent.env(env)$a, 1)
  expect_equal(parent.env(env)$b, 2)
})
