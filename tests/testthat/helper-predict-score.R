# Shared fixtures + assertions for predict_score() tests.
# testthat sources helper-*.R before running the test files.

# Attach the modelling packages when available; individual tests still guard with
# skip_if_not_installed() so the suite degrades gracefully if one is missing.
if (requireNamespace("gamlss",  quietly = TRUE)) suppressMessages(library(gamlss))
if (requireNamespace("gamlss2", quietly = TRUE)) suppressMessages(library(gamlss2))
# attach gamlssTools (brings its dplyr etc. onto the search path) when available
if (requireNamespace("gamlssTools", quietly = TRUE)) suppressMessages(library(gamlssTools))

## ---- data generators -------------------------------------------------------

# Normal response: each site shifts mu (the "batch effect"). `extra_factor` adds
# a second, NON-batch factor `z` (used to test unseen-level handling).
sim_site_normal <- function(n_per = 400, seed = 1,
                            shift = c(site1 = 0, site2 = 0.8, site3 = -0.6, siteNEW = 1.5),
                            extra_factor = FALSE) {
  set.seed(seed)
  d <- do.call(rbind, lapply(names(shift), function(s) {
    x  <- runif(n_per, 0, 10)
    df <- data.frame(x = x, y = rnorm(n_per, 2 + 0.5 * x + shift[[s]], 1), site = s)
    if (extra_factor) df$z <- sample(c("A", "B", "C"), n_per, replace = TRUE)
    df
  }))
  d$site <- factor(d$site)
  if (extra_factor) d$z <- factor(d$z)
  d
}

# Positive response from a GG (generalised gamma); site shifts log-mu.
sim_site_gg <- function(n_per = 400, seed = 1,
                        shift = c(site1 = 0, site2 = 0.3, site3 = -0.2, siteNEW = 0.5)) {
  set.seed(seed)
  d <- do.call(rbind, lapply(names(shift), function(s) {
    x <- runif(n_per, 0, 10)
    data.frame(x = x,
               y = gamlss.dist::rGG(n_per, mu = exp(0.5 + 0.1 * x + shift[[s]]),
                                    sigma = 0.3, nu = 2),
               site = s)
  }))
  d$site <- factor(d$site)
  d
}

# Split into the subjects to score (siteNEW), the "known" sites, and the full data.
split_site <- function(d, new = "siteNEW") {
  list(new   = droplevels(subset(d, site == new)),
       train = droplevels(subset(d, site != new)),
       full  = d)
}

## ---- scoring + assertion helpers -------------------------------------------

# Score the unseen-site fit both ways (centile and quantile-residual / z-score)
# with the shared adjust=TRUE / rm.term="site" settings the agreement tests use.
score_both <- function(fitB, newdata, newformula, which.params, refdata = newdata) {
  common <- list(fitB, newdata = newdata, refdata = refdata, adjust = TRUE,
                 rm.term = "site", newformula = newformula, which.params = which.params)
  list(cent  = do.call(predict_score, c(common, type = "cent")),
       quant = do.call(predict_score, c(common, type = "quantile")))
}

# Gold-standard scores A (model INCLUDES siteNEW) and out-of-sample scores B
# (unseen-site adjustment) should agree up to sampling noise in shared coefs.
# Pass zA/zB to also check z-score (quantile-residual) agreement; the z-scale is
# wider than [0,1], so its tolerances are looser than the centile ones.
expect_agree <- function(a, b, cor_min = 0.99, max_diff = 0.05, mean_diff = 0.015,
                         zA = NULL, zB = NULL,
                         z_cor_min = 0.99, z_max_diff = 0.5, z_mean_diff = 0.05,
                         info = NULL) {
  testthat::expect_equal(length(a), length(b))
  testthat::expect_true(all(is.finite(a)) && all(is.finite(b)),
                        info = paste("non-finite scores", info))
  testthat::expect_gt(stats::cor(a, b), cor_min)
  testthat::expect_lt(max(abs(a - b)), max_diff)
  testthat::expect_lt(mean(a - b), mean_diff)
  if (!is.null(zA)) {
    testthat::expect_equal(length(zA), length(zB))
    testthat::expect_true(all(is.finite(zA)) && all(is.finite(zB)),
                          info = paste("non-finite z-scores", info))
    testthat::expect_gt(stats::cor(zA, zB), z_cor_min)
    testthat::expect_lt(max(abs(zA - zB)), z_max_diff)
    testthat::expect_lt(abs(mean(zA - zB)), z_mean_diff)
  }
}
