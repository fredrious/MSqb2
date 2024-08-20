library(testthat)

test_that("`%<-1%` assigns value to parent environment", {
  # Define a variable in the global environment
  a <- NULL

  # Define a function to use as a test environment
  test_func <- function() {
    # Use the custom operator to assign a value to `a` in the global environment
    "a" %<-1% 1

    # Check if `a` has been assigned in the parent environment
    expect_equal(a, 1)
  }

  # Call the test function
  test_func()
})
