# The gamlss method reconstructs offsets from the fitted object, so it should
# score new data with the original training frame removed from scope.

skip_if_not_installed("gamlss")

# gold standard on the full (incl-siteNEW) fit; `og` = the fit's training frame,
# passed to predictAll() so it isn't re-evaluated from the (local) call.
gold_no <- function(fitA, score_dat, og) {
  pA <- predictAll(fitA, newdata = score_dat, data = og)
  gamlss.dist::pNO(score_dat$y, mu = pA$mu, sigma = pA$sigma)
}

# fit B on a throwaway copy, remove it, then score -- if the data-free path
# failed and fell back to predict(), this would error (nothing to re-evaluate).
score_without_data <- function(formula, d) {
  train_only <- d$train
  fit <- gamlss(formula, sigma.formula = ~1, family = NO, data = train_only, trace = FALSE)
  rm(train_only)
  predict_score(fit, newdata = d$new, refdata = d$new, type = "cent", adjust = TRUE,
                rm.term = "site", newformula = y ~ offset(mu), which.params = "mu")
}

test_that("pb() model scores with no training data in scope", {
  d <- split_site(sim_site_normal())
  fitA <- gamlss(y ~ pb(x) + site, sigma.formula = ~1, family = NO, data = d$full, trace = FALSE)
  expect_agree(gold_no(fitA, d$new, d$full), score_without_data(y ~ pb(x) + site, d))
})

test_that("purely parametric model scores with no training data in scope", {
  d <- split_site(sim_site_normal())
  fitA <- gamlss(y ~ x + site, sigma.formula = ~1, family = NO, data = d$full, trace = FALSE)
  expect_agree(gold_no(fitA, d$new, d$full), score_without_data(y ~ x + site, d))
})

test_that("a non-reconstructable smooth (cs) requires traindata", {
  d <- split_site(sim_site_normal())
  fit <- gamlss(y ~ cs(x) + site, sigma.formula = ~1, family = NO, data = d$train, trace = FALSE)
  expect_error(
    predict_score(fit, d$new, d$new, type = "cent", adjust = TRUE, rm.term = "site",
                  newformula = y ~ offset(mu), which.params = "mu"),
    regexp = "cs/ps/ga/s|traindata"
  )
  expect_no_error(
    predict_score(fit, d$new, d$new, type = "cent", adjust = TRUE, rm.term = "site",
                  newformula = y ~ offset(mu), which.params = "mu", traindata = d$train)
  )
})
