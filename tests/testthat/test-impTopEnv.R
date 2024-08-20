test_that("impTopEnv imports objects from parent environment", {

f1 <- function() {
  a <- 2
  f2()
}
f2 <- function() {
  impTopEnv()
}
out <- f1()

# Check if variable `a` is now in the test environment
expect_true(exists("a", envir = out))

# Check the values of the imported variables
expect_equal(out$a, 2)

})

