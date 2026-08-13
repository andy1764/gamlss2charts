# Same A-vs-B agreement checks for the gamlss method, on BOTH centiles and
# z-scores (quantile residuals). Gold-standard centiles come directly from the
# fitted distribution (no gamlssTools dependency); gold-standard z = qnorm(cent),
# unseen-site z = predict_score(type = "quantile").

skip_if_not_installed("gamlss")

test_that("Basic: mu-only site effect", {
  d <- split_site(sim_site_normal())
  fitA <- gamlss(y ~ x + site, sigma.formula = ~1, family = NO, data = d$full, trace = FALSE)
  pA   <- predictAll(fitA, newdata = d$new, data = d$full)
  centA <- gamlss.dist::pNO(d$new$y, mu = pA$mu, sigma = pA$sigma)

  fitB <- gamlss(y ~ x + site, sigma.formula = ~1, family = NO, data = d$train, trace = FALSE)
  b <- score_both(fitB, d$new, y ~ offset(mu), "mu")
  expect_agree(centA, b$cent, zA = qnorm(centA), zB = b$quant)
})

test_that("Batch effect in mu AND sigma (adjust both)", {
  d <- split_site(sim_site_normal())
  fitA <- gamlss(y ~ x + site, sigma.formula = ~ x + site, family = NO, data = d$full, trace = FALSE)
  pA   <- predictAll(fitA, newdata = d$new, data = d$full)
  centA <- gamlss.dist::pNO(d$new$y, mu = pA$mu, sigma = pA$sigma)

  fitB <- gamlss(y ~ x + site, sigma.formula = ~ x + site, family = NO, data = d$train, trace = FALSE)
  b <- score_both(fitB, d$new, y ~ offset(mu) | offset(sigma), c("mu", "sigma"))
  expect_agree(centA, b$cent, zA = qnorm(centA), zB = b$quant)
})

test_that("Penalised spline pb(x)", {
  d <- split_site(sim_site_normal())
  fitA <- gamlss(y ~ pb(x) + site, sigma.formula = ~1, family = NO, data = d$full, trace = FALSE)
  pA   <- predictAll(fitA, newdata = d$new, data = d$full)
  centA <- gamlss.dist::pNO(d$new$y, mu = pA$mu, sigma = pA$sigma)

  fitB <- gamlss(y ~ pb(x) + site, sigma.formula = ~1, family = NO, data = d$train, trace = FALSE)
  b <- score_both(fitB, d$new, y ~ offset(mu), "mu")
  expect_agree(centA, b$cent, zA = qnorm(centA), zB = b$quant)
})

test_that("Non-normal family (GG)", {
  d <- split_site(sim_site_gg())
  fitA <- gamlss(y ~ x + site, family = GG, data = d$full, trace = FALSE)
  pA   <- predictAll(fitA, newdata = d$new, data = d$full)
  centA <- gamlss.dist::pGG(d$new$y, mu = pA$mu, sigma = pA$sigma, nu = pA$nu)

  fitB <- gamlss(y ~ x + site, sigma.formula = ~1, family = GG, data = d$train, trace = FALSE)
  b <- score_both(fitB, d$new, y ~ offset(mu), "mu")
  # GG carries a small systematic mean offset -> looser mean tolerances
  expect_agree(centA, b$cent, mean_diff = 0.03,
               zA = qnorm(centA), zB = b$quant, z_mean_diff = 0.12)
})
