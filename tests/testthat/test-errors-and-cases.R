# Informative errors/warnings, and the cases-vs-controls (refdata != newdata) use.

test_that("unseen level of a NON-rm.term factor errors informatively (gamlss)", {
  skip_if_not_installed("gamlss")
  d <- split_site(sim_site_normal(extra_factor = TRUE))
  fitB <- gamlss(y ~ x + z + site, sigma.formula = ~1, family = NO, data = d$train, trace = FALSE)
  nd <- d$new
  nd$z <- "D"   # level not seen in the fit, and NOT the batch term
  # message should name the offending factor and mention its unseen level
  expect_error(
    predict_score(fitB, newdata = nd, refdata = nd, type = "cent", adjust = TRUE,
                  rm.term = "site", newformula = y ~ offset(mu), which.params = "mu"),
    regexp = "level"
  )
})

test_that("rm.term with adjust=FALSE warns (removed but not re-estimated)", {
  skip_if_not_installed("gamlss2")
  d <- split_site(sim_site_normal())
  fitB <- gamlss2(y ~ x + site | 1, family = NO, data = d$train)
  expect_warning(
    predict_score(fitB, newdata = d$new, refdata = d$new, type = "cent", adjust = FALSE,
                  rm.term = "site", which.params = "mu"),
    regexp = "adjust = FALSE"
  )
})

test_that("Cases & controls: newdata (patients) differs from refdata (controls)", {
  skip_if_not_installed("gamlss2")
  d <- split_site(sim_site_normal())
  set.seed(2)
  idx <- sample(nrow(d$new), 150, replace = TRUE)   # different size than the controls
  patients <- d$new[idx, ]
  patients$y <- patients$y + rnorm(150, 0, 0.75)     # + patient noise
  rownames(patients) <- NULL

  fitA  <- gamlss2(y ~ x + site | 1, family = NO, data = d$full)
  centA <- fitA$family$cdf(patients$y, par = predict(fitA, newdata = patients, type = "parameter"))
  fitB  <- gamlss2(y ~ x + site | 1, family = NO, data = d$train)
  centB <- predict_score(fitB, newdata = patients, refdata = d$new, type = "cent",
                         adjust = TRUE, rm.term = "site",
                         newformula = y ~ offset(mu), which.params = "mu")
  expect_equal(length(centB), nrow(patients))   # scores follow newdata, not refdata
  expect_agree(centA, centB)
})
