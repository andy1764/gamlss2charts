# Random-effect batch terms, and z-score / quantile-residual behaviour.

test_that("gamlss random(site) is supported and matches the gold standard", {
  skip_if_not_installed("gamlss")
  d <- split_site(sim_site_normal())
  # gold: fit includes siteNEW, so its random effect is estimated directly
  fitA <- gamlss(y ~ x + random(site), sigma.formula = ~1, family = NO, data = d$full, trace = FALSE)
  pA   <- predictAll(fitA, newdata = d$new, data = d$full)
  centA <- gamlss.dist::pNO(d$new$y, mu = pA$mu, sigma = pA$sigma)

  fitB <- gamlss(y ~ x + random(site), sigma.formula = ~1, family = NO, data = d$train, trace = FALSE)
  centB <- predict_score(fitB, d$new, d$new, type = "cent", adjust = TRUE, rm.term = "site",
                         newformula = y ~ offset(mu), which.params = "mu")
  expect_agree(centA, centB)
})
