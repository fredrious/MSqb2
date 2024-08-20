library(testthat)

test_that("`%<-g%` assigns value to global environment", {
  # Define a variable in the global environment
  global_var <- NULL

  # Define a function to use as a test environment
  test_func <- function() {
    # Use the custom operator to assign a value to `global_var` in the global environment
    "global_var" %<-g% 100

    # Check if `global_var` has been assigned in the global environment
    expect_equal(global_var, 100)
  }

  # Call the test function
  test_func()
})
