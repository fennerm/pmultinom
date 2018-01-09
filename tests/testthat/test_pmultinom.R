context("pmultinom")

reps <- 100000
k <- 10
prob <- rep(1 / k, k)
n <- round(runif(1, 5, 20))
q <- c(1:n)
multinomials <- rmultinom(reps, n, prob)
maxs <- matrixStats::colMaxs(multinomials)

test_that("dens_func matches monte carlo expectations (lower.tail=TRUE))", {
            lapply(q, function(qi) {
                     print(paste("For q==", as.character(qi), sep = ""))
                     have_x_or_less <- which(maxs <= qi)
                     p_monte_carlo <- (length(have_x_or_less) + 1) / (reps + 1)
                     p <- pmultinom(qi, n, prob)
                     expect_true(abs(p_monte_carlo - p) < 0.02)
})
})

test_that("dens_func matches monte carlo expectations (lower.tail=FALSE))", {
            lapply(q, function(qi) {
                     print(paste("For q==", as.character(qi), sep = ""))
                     have_over_x <- which(maxs > qi)
                     p_monte_carlo <- (length(have_over_x) + 1) / (reps + 1)
                     p <- pmultinom(qi, n, prob, lower.tail = FALSE)
                     expect_true(abs(p_monte_carlo - p) < 0.02)
})
})

test_that("Multiple quantiles handled correctly", {
            k <- 2
            q <- c(0, 4)
            prob <- rep(0.5, 2)
            n <- 4
            multinomials <- rmultinom(reps, n, prob)
            have_x_or_less <- which(apply(multinomials <= q, 2, all))
            p_monte_carlo <- (length(have_x_or_less) + 1) / (reps + 1)
            p <- pmultinom(q, n, prob)
            expect_true(abs(p_monte_carlo - p) < 0.02)
})

test_that("Different probabilities handled correctly", {
            prob <- runif(k)
            prob <- prob / sum(prob)
            multinomials <- rmultinom(reps, n, prob)
            maxs <- matrixStats::colMaxs(multinomials)
            lapply(q, function(qi) {
                     print(paste("For q==", as.character(qi), sep = ""))
                     have_x_or_less <- which(maxs <= qi)
                     p_monte_carlo <- (length(have_x_or_less) + 1) / (reps + 1)
                     p <- pmultinom(qi, n, prob)
                     expect_true(abs(p_monte_carlo - p) < 0.02)
})
})

test_that("q == 0 handled correctly", {
            expect_equal(pmultinom(0, n, prob), 0)
})

test_that("q > n handled correctly", {
            expect_equal(pmultinom(n + 3, n, prob), 1)
})

test_that("lower.tail == FALSE and q == n handled correctly", {
            expect_equal(pmultinom(c(1, 1), 1, c(0.5, 0.5), lower.tail = FALSE),
                         0)
})

test_that("Underflow handled correctly", {
            expect_equal(pmultinom(1e6, 1e6 + 1, prob = c(0.5, 0.5),
                                   lower.tail = FALSE),
                         "<2.2e-16")
})
