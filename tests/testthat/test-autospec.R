# Auto-derivation of newformula / which.params from the moments that contain
# rm.term. When both args are left NULL, predict_score should adjust exactly the
# moments carrying rm.term (offset(param)) and freeze the rest (offset(param) - 1).
# Supplying either arg disables auto-derivation (manual mode).

skip_if_not_installed("gamlss2")

# gold standard: a fit that INCLUDES the batch level, scored in-sample.
gold_g2 <- function(fitA, score_dat)
  fitA$family$cdf(q = score_dat$y,
                  par = predict(fitA, newdata = score_dat, type = "parameter"))

# resolved (auto) spec for a gamlss2 fit
auto_spec_g2 <- function(fit, rm.term = "site")
  .resolve_adjust_spec(NULL, NULL, .moment_formulas_gamlss2(fit), rm.term)

## ---- gamlss2: correct moment selection -------------------------------------

test_that("gamlss2 auto: site in mu only -> adjust mu, freeze sigma", {
  d  <- split_site(sim_site_normal())
  fB <- gamlss2(y ~ x + site | 1, family = NO, data = d$train)
  s  <- auto_spec_g2(fB)
  expect_equal(unname(s$which.params), "mu")
  expect_equal(deparse1(s$newformula), "y ~ offset(mu) | offset(sigma) - 1")
})

test_that("gamlss2 auto: site in mu AND sigma -> adjust both", {
  d  <- split_site(sim_site_normal())
  fB <- gamlss2(y ~ x + site | x + site, family = NO, data = d$train)
  s  <- auto_spec_g2(fB)
  expect_equal(unname(s$which.params), c("mu", "sigma"))
  expect_equal(deparse1(s$newformula), "y ~ offset(mu) | offset(sigma)")
})

test_that("gamlss2 auto (GG): non-rm.term moments frozen with offset()-1", {
  d  <- split_site(sim_site_gg())
  fB <- gamlss2(y ~ x + site | 1 | 1, family = GG, data = d$train)
  s  <- auto_spec_g2(fB)
  expect_equal(unname(s$which.params), "mu")
  expect_equal(deparse1(s$newformula),
               "y ~ offset(mu) | offset(sigma) - 1 | offset(nu) - 1")
})

## ---- gamlss2: end-to-end auto call agrees with the gold standard -----------

test_that("gamlss2 auto call (no formula args) agrees with in-sample gold", {
  d     <- split_site(sim_site_gg())
  fitA  <- gamlss2(y ~ x + site,         family = GG, data = d$full)
  fitB  <- gamlss2(y ~ x + site | 1 | 1, family = GG, data = d$train)
  centA <- gold_g2(fitA, d$new)
  centB <- predict_score(fitB, d$new, refdata = d$new, adjust = TRUE,
                         rm.term = "site", type = "cent")           # auto
  quantB <- predict_score(fitB, d$new, refdata = d$new, adjust = TRUE,
                          rm.term = "site", type = "quantile")      # auto
  # GG carries a small systematic mean offset -> looser mean tolerances
  expect_agree(centA, centB, mean_diff = 0.03,
               zA = qnorm(centA), zB = quantB, z_mean_diff = 0.12)
})

## ---- manual mode and fallback ----------------------------------------------

test_that("supplying either arg disables auto-derivation", {
  # which.params only -> newformula falls back to the legacy mu|sigma default
  s1 <- .resolve_adjust_spec(NULL, "mu",
                             list(mu = ~ site, sigma = ~ 1, nu = ~ 1), "site")
  expect_equal(s1$which.params, "mu")
  expect_equal(deparse1(s1$newformula), "y ~ offset(mu) | offset(sigma)")
  # newformula only -> which.params falls back to legacy c("mu","sigma")
  s2 <- .resolve_adjust_spec(y ~ offset(mu), NULL,
                             list(mu = ~ site, sigma = ~ 1), "site")
  expect_equal(s2$which.params, c("mu", "sigma"))
})

test_that("rm.term absent from every moment errors", {
  expect_error(
    .resolve_adjust_spec(NULL, NULL,
                         list(mu = ~ x, sigma = ~ 1), "site"),
    "not found in any moment")
})

## ---- gamlss (v1) method -----------------------------------------------------

test_that("gamlss v1 auto: site in mu only -> adjust mu, freeze sigma", {
  skip_if_not_installed("gamlss")
  d  <- split_site(sim_site_normal())
  fB <- gamlss::gamlss(y ~ x + site, family = NO, data = d$train, trace = FALSE)
  s  <- .resolve_adjust_spec(NULL, NULL,
                             .moment_formulas_gamlss(fB), "site")
  expect_equal(unname(s$which.params), "mu")
  expect_equal(deparse1(s$newformula), "y ~ offset(mu) | offset(sigma) - 1")
})

test_that("gamlss v1 auto: site in mu AND sigma -> adjust both", {
  skip_if_not_installed("gamlss")
  d  <- split_site(sim_site_normal())
  fB <- gamlss::gamlss(y ~ x + site, sigma.formula = ~ x + site,
                       family = NO, data = d$train, trace = FALSE)
  s  <- .resolve_adjust_spec(NULL, NULL,
                             .moment_formulas_gamlss(fB), "site")
  expect_equal(unname(s$which.params), c("mu", "sigma"))
  expect_equal(deparse1(s$newformula), "y ~ offset(mu) | offset(sigma)")
})

test_that("gamlss v1 auto (GG): non-rm.term moments frozen with offset()-1", {
  skip_if_not_installed("gamlss")
  d  <- split_site(sim_site_gg())
  fB <- gamlss::gamlss(y ~ x + site, sigma.formula = ~ 1,
                       family = GG, data = d$train, trace = FALSE)
  s  <- .resolve_adjust_spec(NULL, NULL,
                             .moment_formulas_gamlss(fB), "site")
  expect_equal(unname(s$which.params), "mu")
  expect_equal(deparse1(s$newformula),
               "y ~ offset(mu) | offset(sigma) - 1 | offset(nu) - 1")
})

test_that("gamlss v1 auto call (no formula args) agrees with in-sample gold", {
  skip_if_not_installed("gamlss")
  d     <- split_site(sim_site_gg())
  fitA  <- gamlss::gamlss(y ~ x + site, family = GG, data = d$full, trace = FALSE)
  pA    <- predictAll(fitA, newdata = d$new, data = d$full)
  centA <- gamlss.dist::pGG(d$new$y, mu = pA$mu, sigma = pA$sigma, nu = pA$nu)

  fitB   <- gamlss::gamlss(y ~ x + site, sigma.formula = ~ 1,
                           family = GG, data = d$train, trace = FALSE)
  centB  <- predict_score(fitB, d$new, refdata = d$new, adjust = TRUE,
                          rm.term = "site", type = "cent")           # auto
  quantB <- predict_score(fitB, d$new, refdata = d$new, adjust = TRUE,
                          rm.term = "site", type = "quantile")       # auto
  # GG carries a small systematic mean offset -> looser mean tolerances
  expect_agree(centA, centB, mean_diff = 0.03,
               zA = qnorm(centA), zB = quantB, z_mean_diff = 0.12)
})
