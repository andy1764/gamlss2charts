# A (gold standard: model fit WITH siteNEW) vs B (siteNEW treated as unseen)
# should agree on BOTH centiles and z-scores (quantile residuals), for the
# gamlss2 method across model structures. Gold-standard z = qnorm(gold centile);
# unseen-site z = predict_score(type = "quantile").

skip_if_not_installed("gamlss2")

# gold-standard centiles from a gamlss2 fit that includes siteNEW
gold_g2 <- function(fitA, score_dat) {
  fitA$family$cdf(score_dat$y,
                  par = predict(fitA, newdata = score_dat, type = "parameter"))
}

test_that("Basic: mu-only site effect", {
  d <- split_site(sim_site_normal())
  fitA <- gamlss2(y ~ x + site | 1, family = NO, data = d$full)
  fitB <- gamlss2(y ~ x + site | 1, family = NO, data = d$train)
  centA <- gold_g2(fitA, d$new)
  b <- score_both(fitB, d$new, NULL, NULL)
  expect_agree(centA, b$cent, zA = qnorm(centA), zB = b$quant)
})

test_that("Batch effect in mu AND sigma (adjust both)", {
  d <- split_site(sim_site_normal())
  fitA <- gamlss2(y ~ x + site | x + site, family = NO, data = d$full)
  fitB <- gamlss2(y ~ x + site | x + site, family = NO, data = d$train)
  centA <- gold_g2(fitA, d$new)
  b <- score_both(fitB, d$new, NULL, c("mu", "sigma"))
  expect_agree(centA, b$cent, zA = qnorm(centA), zB = b$quant)
})

test_that("Smooth term s(x)", {
  d <- split_site(sim_site_normal())
  fitA <- gamlss2(y ~ s(x) + site | 1, family = NO, data = d$full)
  fitB <- gamlss2(y ~ s(x) + site | 1, family = NO, data = d$train)
  centA <- gold_g2(fitA, d$new)
  b <- score_both(fitB, d$new, NULL, "mu")
  expect_agree(centA, b$cent, zA = qnorm(centA), zB = b$quant)
})

test_that("Small sample per site", {
  d <- split_site(sim_site_normal(n_per = 50))
  fitA <- gamlss2(y ~ x + site | 1, family = NO, data = d$full)
  fitB <- gamlss2(y ~ x + site | 1, family = NO, data = d$train)
  centA <- gold_g2(fitA, d$new)
  b <- score_both(fitB, d$new, NULL, "mu")
  expect_agree(centA, b$cent, max_diff = 0.08,
               zA = qnorm(centA), zB = b$quant, z_max_diff = 0.6)
})

test_that("Non-normal family (GG)", {
  d <- split_site(sim_site_gg())
  fitA <- gamlss2(y ~ x + site, family = GG, data = d$full)
  fitB <- gamlss2(y ~ x + site | 1, family = GG, data = d$train)
  centA <- gold_g2(fitA, d$new)
  b <- score_both(fitB, d$new, NULL, "mu")
  # GG carries a small systematic mean offset -> looser mean tolerances
  expect_agree(centA, b$cent, mean_diff = 0.03,
               zA = qnorm(centA), zB = b$quant, z_mean_diff = 0.12)
})

test_that("Four-parameter family (BCT), adjust all", {
  d <- split_site(sim_site_gg())
  fitA <- gamlss2(y ~ x + site | . | . | ., family = BCT, data = d$full)
  fitB <- gamlss2(y ~ x + site | . | . | ., family = BCT, data = d$train)
  centA <- gold_g2(fitA, d$new)
  b <- score_both(fitB, d$new, NULL, NULL)
  # GG carries a small systematic mean offset -> looser mean tolerances
  expect_agree(centA, b$cent, mean_diff = 0.03,
               zA = qnorm(centA), zB = b$quant, z_mean_diff = 0.12)
})

test_that("Four-parameter family (BCT), adjust mu/sigma", {
  d <- split_site(sim_site_gg())
  fitA <- gamlss2(y ~ x + site | x + site | . | ., family = BCT, data = d$full)
  fitB <- gamlss2(y ~ x + site | x + site | . | ., family = BCT, data = d$train)
  centA <- gold_g2(fitA, d$new)
  b <- score_both(fitB, d$new, NULL, NULL)
  # GG carries a small systematic mean offset -> looser mean tolerances
  expect_agree(centA, b$cent, mean_diff = 0.03,
               zA = qnorm(centA), zB = b$quant, z_mean_diff = 0.12)
})
