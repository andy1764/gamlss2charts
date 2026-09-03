test_that("log-likelihood is accurate within training set", {
  d <- split_site(sim_site_gg())
  fitB <- gamlss2(y ~ x + site | . | ., family = GG, data = d$train)

  expect_equal(predict_score(fitB, d$train, adjust = FALSE, type = "logLik"), fitB$logLik)
})

test_that("log-likelihood of adjustment model is the same as logLik from predict_score for intercept only", {
  d <- split_site(sim_site_gg())
  fitB <- gamlss2(y ~ x + site | . | ., family = GG, data = d$train)

  # treat training set as newdata for simplicity
  object <- fitB
  newdata <- d$train
  pred <- predict(object, newdata = newdata, type = "link")
  newdata <- cbind(newdata, pred)
  fit2 <- gamlss2(y ~ offset(mu) | offset(sigma) | offset(nu), data = newdata, family = object$family)

  expect_equal(predict_score(fitB, d$train, adjust = TRUE, which.params = c("mu", "sigma", "nu"), type = "logLik"), fit2$logLik)
})

test_that("log-likelihood of adjustment model not the same as logLik from predict_score with covariates", {
  d <- split_site(sim_site_gg())
  fitB <- gamlss2(y ~ x + site | . | ., family = GG, data = d$train)

  d$train$x2 <- rbinom(nrow(d$train), 1, 0.5)
  refdata <- d$train
  params <- predict(fitB, newdata = d$train)
  refdata <- cbind(refdata, params)
  fit2 <- gamlss2(y ~ offset(mu) + x2 | offset(sigma) | offset(nu), family = GG, data = refdata)

  expect_false(isTRUE(all.equal(
    predict_score(fitB, d$train, newformula = y ~ offset(mu) + x2 | offset(sigma) | offset(nu),
                  type = "logLik"), fit2$logLik)))
})

test_that("log-likelihood is comparable versus intercept only model", {
  set.seed(8888)
  d <- split_site(sim_site_gg())
  fitB <- gamlss2(y ~ x + site | . | ., family = GG, data = d$train)

  d$train$x2 <- rbinom(nrow(d$train), 1, 0.5)

  l_null <- predict_score(fitB, d$train, newformula = y ~ offset(mu) | offset(sigma) | offset(nu), type = "logLik")
  l_full <- predict_score(fitB, d$train, newformula = y ~ offset(mu) + x2 | offset(sigma) | offset(nu), type = "logLik")

  expect_equal(l_null, l_full, tolerance = 2)

  expect_true(pchisq(-2*(l_null - l_full), 1) > 0.05)
})
