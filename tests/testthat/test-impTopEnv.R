test_that("impTopEnv imports objects from parent environment", {

f1 <- function() {
  testvar <- 2
  f2()
}
f2 <- function() {
  MSqb2:::impTopEnv()
}
out <- f1()

# Check if variable `a` is now in the test environment
expect_true(exists("testvar", envir = out))

# Check the values of the imported variables
expect_equal(out$testvar, 2)

})

